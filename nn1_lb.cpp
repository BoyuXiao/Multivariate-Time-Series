#include "nn1_lb.h"
#include "funcs.h"
#include "def.h"

#include <algorithm>
#include <cmath>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <optional>
#include <sstream>
#include <unordered_set>
#include <vector>

namespace {

enum class Nn1LbMethod {
    MethodMv,
    MethodPMv,
    MethodWebMv,
    MethodMc1,
    MethodMc2,
    MethodTi,
    MethodPc,
    MethodMc3,
    MethodMc3e
};

// 离线转置结果，不计入任何计时
struct OfflineTranspose {
    std::vector<std::vector<std::vector<double>>> test_t;
    std::vector<std::vector<std::vector<double>>> train_t;
};

struct QueryMethodCache {
    std::vector<std::vector<double>> q_l;
    std::vector<std::vector<double>> q_u;
    std::vector<std::vector<double>> l_uq;
    std::vector<std::vector<double>> u_lq;

    std::vector<double> mu_q;
    std::vector<double> sigma_q;
    std::vector<double> mu_q_l;
    std::vector<double> mu_q_u;
    std::vector<double> sigma_q_l;
    std::vector<double> sigma_q_u;

    std::vector<double> u_mu_q;
    std::vector<double> l_mu_q;
    std::vector<double> u_sigma_q;
    std::vector<double> l_sigma_q;
    std::vector<std::vector<std::vector<Cell>>> t_cells;
    std::vector<double> interval_mu;
    std::vector<double> interval_sigma;
    std::vector<int> segment_num;
    std::vector<double> dbl;
    std::vector<double> dbr;
    std::vector<double> dtl;
    std::vector<double> dtr;

    std::vector<std::vector<std::vector<std::vector<double>>>> pc_bboxes;
};

struct TrainMethodCache {
    std::vector<std::vector<double>> mu_a;
    std::vector<std::vector<double>> sigma_a;
};

double elapsed_sec(clock_t start, clock_t end) {
    return static_cast<double>(end - start) / CLOCKS_PER_SEC;
}

void offline_transpose_all(const std::vector<std::vector<std::vector<double>>>& test_set,
                           const std::vector<std::vector<std::vector<double>>>& train_set,
                           OfflineTranspose& offline) {
    offline.test_t.resize(test_set.size());
    for (size_t i = 0; i < test_set.size(); ++i) {
        offline.test_t[i] = transpose(test_set[i]);
    }
    offline.train_t.resize(train_set.size());
    for (size_t j = 0; j < train_set.size(); ++j) {
        offline.train_t[j] = transpose(train_set[j]);
    }
}

void build_mv_envelopes(const std::vector<std::vector<double>>& series_t,
                        int m,
                        int k_dims,
                        int r,
                        std::vector<std::vector<double>>& lower,
                        std::vector<std::vector<double>>& upper) {
    lower.assign(k_dims, std::vector<double>(m));
    upper.assign(k_dims, std::vector<double>(m));
    for (int k = 0; k < k_dims; ++k) {
        lower_upper_lemire(const_cast<std::vector<double>&>(series_t[k]), m, r, lower[k], upper[k]);
    }
}

void build_web_envelopes(const std::vector<std::vector<double>>& series_t,
                         int m,
                         int k_dims,
                         int r,
                         std::vector<std::vector<double>>& l_uq,
                         std::vector<std::vector<double>>& u_lq) {
    l_uq.assign(k_dims, std::vector<double>(m));
    u_lq.assign(k_dims, std::vector<double>(m));
    for (int k = 0; k < k_dims; ++k) {
        lower_lemire(const_cast<std::vector<double>&>(series_t[k]), m, r, l_uq[k]);
        upper_lemire(const_cast<std::vector<double>&>(series_t[k]), m, r, u_lq[k]);
    }
}

void prepare_mc3_query(const std::vector<std::vector<double>>& query,
                       int m,
                       int k_dims,
                       int r,
                       QueryMethodCache& cache) {
    cache.mu_q.assign(m, 0.0);
    cache.sigma_q.assign(m, 0.0);
    mu_sigma(const_cast<std::vector<std::vector<double>>&>(query), cache.mu_q, cache.sigma_q);

    cache.t_cells.assign(m, std::vector<std::vector<Cell>>());
    cache.interval_mu.assign(m, 0.0);
    cache.interval_sigma.assign(m, 0.0);
    cache.segment_num.assign(m, 0);
    cache.u_mu_q.assign(m, 0.0);
    cache.l_mu_q.assign(m, 0.0);
    cache.u_sigma_q.assign(m, 0.0);
    cache.l_sigma_q.assign(m, 0.0);

    auto dnn = Preprocessing_MC3(const_cast<std::vector<std::vector<double>>&>(query),
                               cache.mu_q, cache.sigma_q,
                               cache.u_mu_q, cache.l_mu_q, cache.u_sigma_q, cache.l_sigma_q,
                               cache.t_cells, cache.interval_mu, cache.interval_sigma, cache.segment_num, r);
    cache.dbl = dnn[0];
    cache.dbr = dnn[1];
    cache.dtl = dnn[2];
    cache.dtr = dnn[3];
}

void prepare_query_for_method(Nn1LbMethod kind,
                              const std::vector<std::vector<double>>& query,
                              const std::vector<std::vector<double>>& q_t,
                              int m,
                              int k_dims,
                              int r,
                              QueryMethodCache& cache) {
    switch (kind) {
        case Nn1LbMethod::MethodMv:
        case Nn1LbMethod::MethodPMv:
            build_mv_envelopes(q_t, m, k_dims, r, cache.q_l, cache.q_u);
            break;
        case Nn1LbMethod::MethodWebMv:
            build_mv_envelopes(q_t, m, k_dims, r, cache.q_l, cache.q_u);
            build_web_envelopes(q_t, m, k_dims, r, cache.l_uq, cache.u_lq);
            break;
        case Nn1LbMethod::MethodMc1:
            cache.mu_q.assign(m, 0.0);
            cache.sigma_q.assign(m, 0.0);
            mu_sigma(const_cast<std::vector<std::vector<double>>&>(query), cache.mu_q, cache.sigma_q);
            cache.mu_q_l.assign(m, 0.0);
            cache.mu_q_u.assign(m, 0.0);
            lower_upper_lemire(cache.mu_q, m, r, cache.mu_q_l, cache.mu_q_u);
            break;
        case Nn1LbMethod::MethodMc2:
            cache.mu_q.assign(m, 0.0);
            cache.sigma_q.assign(m, 0.0);
            mu_sigma(const_cast<std::vector<std::vector<double>>&>(query), cache.mu_q, cache.sigma_q);
            cache.mu_q_l.assign(m, 0.0);
            cache.mu_q_u.assign(m, 0.0);
            cache.sigma_q_l.assign(m, 0.0);
            cache.sigma_q_u.assign(m, 0.0);
            lower_upper_lemire(cache.mu_q, m, r, cache.mu_q_l, cache.mu_q_u);
            lower_upper_lemire(cache.sigma_q, m, r, cache.sigma_q_l, cache.sigma_q_u);
            break;
        case Nn1LbMethod::MethodMc3:
        case Nn1LbMethod::MethodMc3e:
            prepare_mc3_query(query, m, k_dims, r, cache);
            break;
        case Nn1LbMethod::MethodPc:
            cache.pc_bboxes = LB_PC_Preprocess(query, r, k_dims, 2.0);
            break;
        case Nn1LbMethod::MethodTi:
            break;
    }
}

bool method_needs_train_mu_sigma(Nn1LbMethod kind) {
    switch (kind) {
        case Nn1LbMethod::MethodMc1:
        case Nn1LbMethod::MethodMc2:
        case Nn1LbMethod::MethodMc3:
        case Nn1LbMethod::MethodMc3e:
            return true;
        default:
            return false;
    }
}

void prepare_all_trains_for_method(Nn1LbMethod kind,
                                   const std::vector<std::vector<std::vector<double>>>& train_set,
                                   int m,
                                   int k_dims,
                                   TrainMethodCache& cache) {
    if (!method_needs_train_mu_sigma(kind)) {
        return;
    }

    const int n_train = static_cast<int>(train_set.size());
    cache.mu_a.assign(n_train, std::vector<double>(m));
    cache.sigma_a.assign(n_train, std::vector<double>(m, 0.0));

    for (int j = 0; j < n_train; ++j) {
        if (kind == Nn1LbMethod::MethodMc1) {
            only_mu(const_cast<std::vector<std::vector<double>>&>(train_set[j]), cache.mu_a[j]);
        } else {
            mu_sigma(const_cast<std::vector<std::vector<double>>&>(train_set[j]), cache.mu_a[j],
                     cache.sigma_a[j]);
        }
    }
}

double compute_lb(Nn1LbMethod kind,
                  const std::vector<std::vector<double>>& query,
                  const std::vector<std::vector<double>>& train,
                  int train_idx,
                  int m,
                  int k_dims,
                  int r,
                  int ti_p,
                  double best_dist,
                  bool use_early_stopping,
                  QueryMethodCache& qcache,
                  const TrainMethodCache& tcache,
                  std::vector<std::vector<double>>& q_t,
                  std::vector<std::vector<double>>& a_t) {
    auto& q_mut = const_cast<std::vector<std::vector<double>>&>(query);
    auto& t_mut = const_cast<std::vector<std::vector<double>>&>(train);

    if (use_early_stopping) {
        switch (kind) {
            case Nn1LbMethod::MethodMv:
                return LB_MV_early_env(a_t, k_dims, m, qcache.q_l, qcache.q_u, best_dist);
            case Nn1LbMethod::MethodPMv:
                return LB_P_MV_early(a_t, q_t, qcache.q_l, qcache.q_u, r, best_dist);
            case Nn1LbMethod::MethodWebMv:
                return LB_WEB_MV_early(a_t, q_t, qcache.q_l, qcache.q_u, qcache.l_uq, qcache.u_lq, r,
                                       best_dist);
            case Nn1LbMethod::MethodMc1:
                return LB_MC1_early_cached(tcache.mu_a[train_idx], qcache.mu_q,
                                           qcache.mu_q_l, qcache.mu_q_u, k_dims, m, best_dist);
            case Nn1LbMethod::MethodMc2:
                return LB_MC2_early_cached(tcache.mu_a[train_idx], tcache.sigma_a[train_idx],
                                         qcache.mu_q, qcache.sigma_q,
                                         qcache.mu_q_l, qcache.mu_q_u,
                                         qcache.sigma_q_l, qcache.sigma_q_u,
                                         k_dims, m, best_dist);
            case Nn1LbMethod::MethodTi:
                return LB_TI_early(query, train, r, ti_p, best_dist);
            case Nn1LbMethod::MethodPc:
                return LB_PC_early(query, train, qcache.pc_bboxes, best_dist);
            case Nn1LbMethod::MethodMc3:
                return LB_MC3_early(t_mut, qcache.u_mu_q, qcache.l_mu_q, qcache.u_sigma_q, qcache.l_sigma_q,
                                    const_cast<std::vector<double>&>(tcache.mu_a[train_idx]),
                                    const_cast<std::vector<double>&>(tcache.sigma_a[train_idx]),
                                    qcache.dbl, qcache.dbr, qcache.dtl, qcache.dtr,
                                    qcache.t_cells, qcache.interval_mu, qcache.interval_sigma, qcache.segment_num,
                                    best_dist);
            case Nn1LbMethod::MethodMc3e: {
                std::vector<double> lb_record(m, 0.0);
                return LB_MC3E_early(t_mut, q_mut, a_t, q_t,
                                     qcache.u_mu_q, qcache.l_mu_q, qcache.u_sigma_q, qcache.l_sigma_q,
                                     const_cast<std::vector<double>&>(tcache.mu_a[train_idx]),
                                     const_cast<std::vector<double>&>(tcache.sigma_a[train_idx]),
                                     qcache.dbl, qcache.dbr, qcache.dtl, qcache.dtr,
                                     qcache.t_cells, qcache.interval_mu, qcache.interval_sigma, qcache.segment_num,
                                     lb_record, r, best_dist);
            }
        }
        return 0.0;
    }

    switch (kind) {
        case Nn1LbMethod::MethodMv:
            return LB_MV(q_mut, t_mut, k_dims, m, r);
        case Nn1LbMethod::MethodPMv:
            return LB_P_MV(a_t, q_t, qcache.q_l, qcache.q_u, r);
        case Nn1LbMethod::MethodWebMv:
            return LB_WEB_MV(a_t, q_t, qcache.q_l, qcache.q_u, qcache.l_uq, qcache.u_lq, r);
        case Nn1LbMethod::MethodMc1:
            return LB_MC1(q_mut, t_mut, const_cast<std::vector<double>&>(qcache.mu_q), k_dims, m, r);
        case Nn1LbMethod::MethodMc2:
            return LB_MC2(q_mut, t_mut, const_cast<std::vector<double>&>(qcache.mu_q),
                          const_cast<std::vector<double>&>(qcache.sigma_q), k_dims, m, r);
        case Nn1LbMethod::MethodTi:
            return LB_TI(query, train, r, ti_p);
        case Nn1LbMethod::MethodPc:
            return LB_PC(query, train, qcache.pc_bboxes);
        case Nn1LbMethod::MethodMc3:
            return LB_MC3(t_mut, qcache.u_mu_q, qcache.l_mu_q, qcache.u_sigma_q, qcache.l_sigma_q,
                        const_cast<std::vector<double>&>(tcache.mu_a[train_idx]),
                        const_cast<std::vector<double>&>(tcache.sigma_a[train_idx]),
                        qcache.dbl, qcache.dbr, qcache.dtl, qcache.dtr,
                        qcache.t_cells, qcache.interval_mu, qcache.interval_sigma, qcache.segment_num);
        case Nn1LbMethod::MethodMc3e: {
            std::vector<double> lb_record(m, 0.0);
            return LB_MC3E(t_mut, q_mut, a_t, q_t,
                           qcache.u_mu_q, qcache.l_mu_q, qcache.u_sigma_q, qcache.l_sigma_q,
                           const_cast<std::vector<double>&>(tcache.mu_a[train_idx]),
                           const_cast<std::vector<double>&>(tcache.sigma_a[train_idx]),
                           qcache.dbl, qcache.dbr, qcache.dtl, qcache.dtr,
                           qcache.t_cells, qcache.interval_mu, qcache.interval_sigma, qcache.segment_num,
                           lb_record, r);
        }
    }
    return 0.0;
}

constexpr const char* kDtwBaselineCacheFile = "..//dtw_baseline_times.txt";

bool load_dtw_baseline_from_cache(const std::string& dataset_name, double r_ratio, double& out_time) {
    std::ifstream in(kDtwBaselineCacheFile);
    if (!in.is_open()) {
        return false;
    }
    std::string header;
    if (!std::getline(in, header)) {
        return false;
    }
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        std::istringstream iss(line);
        std::string name;
        double cached_ratio = 0.0;
        int cached_r = 0;
        double cached_time = 0.0;
        if (!(iss >> name >> cached_ratio >> cached_r >> cached_time)) {
            continue;
        }
        if (name == dataset_name && std::abs(cached_ratio - r_ratio) < 1e-9) {
            out_time = cached_time;
            return true;
        }
    }
    return false;
}

void append_dtw_baseline_to_cache(const std::string& dataset_name, double r_ratio, int r, double time_sec) {
    const bool exists = std::ifstream(kDtwBaselineCacheFile).good();
    std::ofstream out(kDtwBaselineCacheFile, std::ios::app);
    if (!out.is_open()) {
        return;
    }
    if (!exists) {
        out << "Name\tr_ratio\tr\tDtwNn1Time\n";
    }
    out << dataset_name << '\t' << r_ratio << '\t' << r << '\t' << time_sec << '\n';
}

double run_dtw_only_nn1(const std::vector<std::vector<std::vector<double>>>& test_set,
                        const std::vector<std::vector<std::vector<double>>>& train_set,
                        int m,
                        int k_dims,
                        int r) {
    const int n_test = static_cast<int>(test_set.size());
    const int n_train = static_cast<int>(train_set.size());
    const clock_t start = clock();
    for (int i = 0; i < n_test; ++i) {
        double best_dist = std::numeric_limits<double>::infinity();
        for (int j = 0; j < n_train; ++j) {
            const double dist = dtw(const_cast<std::vector<std::vector<double>>&>(test_set[i]),
                                    const_cast<std::vector<std::vector<double>>&>(train_set[j]),
                                    k_dims, m, r);
            if (dist < best_dist) {
                best_dist = dist;
            }
        }
    }
    return elapsed_sec(start, clock());
}

double get_or_run_dtw_baseline(const std::vector<std::vector<std::vector<double>>>& test_set,
                               const std::vector<std::vector<std::vector<double>>>& train_set,
                               const std::string& dataset_name,
                               double r_ratio,
                               int r,
                               int m,
                               int k_dims) {
    double cached_time = 0.0;
    if (load_dtw_baseline_from_cache(dataset_name, r_ratio, cached_time)) {
        std::cout << "  [DTW baseline] loaded from cache: " << cached_time << "s" << std::endl;
        return cached_time;
    }
    std::cout << "  [DTW baseline] running pure DTW 1NN ..." << std::flush;
    const double time_sec = run_dtw_only_nn1(test_set, train_set, m, k_dims, r);
    std::cout << " " << time_sec << "s (saved to cache)" << std::endl;
    append_dtw_baseline_to_cache(dataset_name, r_ratio, r, time_sec);
    return time_sec;
}

std::vector<Nn1LbMethod> all_nn1_methods() {
    return {Nn1LbMethod::MethodMv,     Nn1LbMethod::MethodPMv,   Nn1LbMethod::MethodWebMv,
            Nn1LbMethod::MethodMc1,    Nn1LbMethod::MethodMc2,   Nn1LbMethod::MethodTi,
            Nn1LbMethod::MethodPc,     Nn1LbMethod::MethodMc3};
}

std::optional<Nn1LbMethod> parse_nn1_method_token(const std::string& raw) {
    std::string t = raw;
    for (char& c : t) {
        if (c >= 'a' && c <= 'z') {
            c = static_cast<char>(c - 'a' + 'A');
        }
    }
    if (t == "MV") {
        return Nn1LbMethod::MethodMv;
    }
    if (t == "P_MV" || t == "PMV") {
        return Nn1LbMethod::MethodPMv;
    }
    if (t == "WEB_MV" || t == "WEBMV") {
        return Nn1LbMethod::MethodWebMv;
    }
    if (t == "MC1") {
        return Nn1LbMethod::MethodMc1;
    }
    if (t == "MC2") {
        return Nn1LbMethod::MethodMc2;
    }
    if (t == "TI") {
        return Nn1LbMethod::MethodTi;
    }
    if (t == "PC") {
        return Nn1LbMethod::MethodPc;
    }
    if (t == "MC3") {
        return Nn1LbMethod::MethodMc3;
    }
    if (t == "MC3E") {
        return Nn1LbMethod::MethodMc3e;
    }
    return std::nullopt;
}

std::vector<Nn1LbMethod> resolve_nn1_methods(const std::vector<std::string>& filter) {
    if (filter.empty()) {
        return all_nn1_methods();
    }
    std::vector<Nn1LbMethod> out;
    for (const std::string& raw : filter) {
        std::stringstream ss(raw);
        std::string token;
        while (std::getline(ss, token, ',')) {
            if (token.empty()) {
                continue;
            }
            auto kind = parse_nn1_method_token(token);
            if (!kind) {
                std::cerr << "Unknown method in filter: " << token << std::endl;
                continue;
            }
            out.push_back(*kind);
        }
    }
    return out;
}

const char* lb_name(Nn1LbMethod kind) {
    switch (kind) {
        case Nn1LbMethod::MethodMv: return "MV";
        case Nn1LbMethod::MethodPMv: return "P_MV";
        case Nn1LbMethod::MethodWebMv: return "WEB_MV";
        case Nn1LbMethod::MethodMc1: return "MC1";
        case Nn1LbMethod::MethodMc2: return "MC2";
        case Nn1LbMethod::MethodTi: return "TI";
        case Nn1LbMethod::MethodPc: return "PC";
        case Nn1LbMethod::MethodMc3: return "MC3";
        case Nn1LbMethod::MethodMc3e: return "MC3E";
    }
    return "UNKNOWN";
}

Nn1LbStats run_one_method(Nn1LbMethod kind,
                          const std::vector<std::vector<std::vector<double>>>& test_set,
                          const std::vector<std::vector<std::vector<double>>>& train_set,
                          const OfflineTranspose& offline,
                          int m,
                          int k_dims,
                          int r,
                          int ti_p,
                          bool use_early_stopping) {
    const int n_test = static_cast<int>(test_set.size());
    const int n_train = static_cast<int>(train_set.size());

    Nn1LbStats stats;
    stats.name = lb_name(kind);

    clock_t prep_clock = 0;
    clock_t lb_clock = 0;
    clock_t dtw_clock = 0;

    TrainMethodCache train_cache;
    clock_t t0 = clock();
    prepare_all_trains_for_method(kind, train_set, m, k_dims, train_cache);
    prep_clock += clock() - t0;

    for (int i = 0; i < n_test; ++i) {
        QueryMethodCache qcache;
        t0 = clock();
        prepare_query_for_method(kind, test_set[i], offline.test_t[i], m, k_dims, r, qcache);
        prep_clock += clock() - t0;

        double best_dist = std::numeric_limits<double>::infinity();

        for (int j = 0; j < n_train; ++j) {
            clock_t lb_start = clock();
            const double lb = compute_lb(kind, test_set[i], train_set[j], j, m, k_dims, r, ti_p, best_dist,
                                         use_early_stopping, qcache, train_cache,
                                         const_cast<std::vector<std::vector<double>>&>(offline.test_t[i]),
                                         const_cast<std::vector<std::vector<double>>&>(offline.train_t[j]));
            lb_clock += clock() - lb_start;

            if (lb >= best_dist) {
                ++stats.pruned;
                continue;
            }

            clock_t dtw_start = clock();
            double dist = 0.0;
            if (use_early_stopping) {
                dist = dtw(const_cast<std::vector<std::vector<double>>&>(test_set[i]),
                           const_cast<std::vector<std::vector<double>>&>(train_set[j]),
                           k_dims, m, r, best_dist);
            } else {
                dist = dtw(const_cast<std::vector<std::vector<double>>&>(test_set[i]),
                           const_cast<std::vector<std::vector<double>>&>(train_set[j]),
                           k_dims, m, r);
            }
            dtw_clock += clock() - dtw_start;
            ++stats.dtw_calls;

            if (dist < best_dist) {
                best_dist = dist;
            }
        }
    }

    stats.prep_time = elapsed_sec(0, prep_clock);
    stats.lb_time = elapsed_sec(0, lb_clock);
    stats.dtw_time = elapsed_sec(0, dtw_clock);
    stats.total_time = stats.prep_time + stats.lb_time + stats.dtw_time;
    return stats;
}

void write_header_if_needed(std::ofstream& out, const std::string& path) {
    static std::unordered_set<std::string> header_written;
    if (header_written.count(path)) {
        return;
    }
    if (std::filesystem::exists(path) && std::filesystem::file_size(path) > 0) {
        header_written.insert(path);
        return;
    }
    out << "Name\tr_ratio\tr\tEarlyStop\tMethod\tPruneRate\tPruned\tDtwCalls\t"
           "DtwBaselineTime\tPrepTime\tLbTime\tDtwTime\tTotalTime\tSpeedup\n";
    header_written.insert(path);
}

void write_stats(std::ofstream& out,
                 const std::string& dataset_name,
                 double r_ratio,
                 int r,
                 bool use_early_stopping,
                 const Nn1LbStats& stats) {
    out << dataset_name << '\t' << r_ratio << '\t' << r << '\t' << (use_early_stopping ? 1 : 0) << '\t'
        << stats.name << '\t'
        << stats.prune_rate() << '\t' << stats.pruned << '\t' << stats.dtw_calls << '\t'
        << stats.dtw_baseline_time << '\t'
        << stats.prep_time << '\t' << stats.lb_time << '\t' << stats.dtw_time << '\t'
        << stats.total_time << '\t' << stats.speedup << '\n';
}

}  // namespace

void run_nn1_lb_test(const std::vector<std::vector<std::vector<double>>>& test_set,
                     const std::vector<std::vector<std::vector<double>>>& train_set,
                     const std::string& dataset_name,
                     const Nn1Config& config) {
    if (test_set.empty() || train_set.empty()) {
        return;
    }

    const int m = static_cast<int>(test_set[0].size());
    const int k_dims = static_cast<int>(test_set[0][0].size());
    const double r_ratio = config.r_ratio;
    const int r = static_cast<int>(std::ceil(r_ratio * m));
    const int ti_p = 5;
    const bool use_early_stopping = config.use_early_stopping;

    const std::vector<Nn1LbMethod> methods = resolve_nn1_methods(config.methods_filter);
    if (methods.empty()) {
        std::cerr << "No valid methods to run (check methods_filter)." << std::endl;
        return;
    }

    OfflineTranspose offline;
    offline_transpose_all(test_set, train_set, offline);

    if (config.result_file.empty()) {
        std::cerr << "result_file is empty; set Nn1Config.result_file in main before testing."
                  << std::endl;
        return;
    }
    std::ofstream out(config.result_file, std::ios::app);
    if (!out.is_open()) {
        std::cerr << "Failed to open result file: " << config.result_file << std::endl;
        return;
    }
    write_header_if_needed(out, config.result_file);

    std::cout << dataset_name << " r_ratio=" << r_ratio << " r=" << r
              << " early_stop=" << (use_early_stopping ? "on" : "off");
    if (!config.methods_filter.empty()) {
        std::cout << " methods=" << config.methods_filter[0];
    }
    std::cout << std::endl;
    const double dtw_baseline_time =
        get_or_run_dtw_baseline(test_set, train_set, dataset_name, r_ratio, r, m, k_dims);

    for (Nn1LbMethod kind : methods) {
        std::cout << "  method=" << lb_name(kind) << std::flush;
        Nn1LbStats stats = run_one_method(kind, test_set, train_set, offline, m, k_dims, r, ti_p,
                                          use_early_stopping);
        stats.dtw_baseline_time = dtw_baseline_time;
        stats.speedup = (stats.total_time > 0.0) ? (dtw_baseline_time / stats.total_time) : 0.0;
        write_stats(out, dataset_name, r_ratio, r, use_early_stopping, stats);
        std::cout << " prune_rate=" << stats.prune_rate()
                  << " prep=" << stats.prep_time << "s"
                  << " lb=" << stats.lb_time << "s"
                  << " dtw=" << stats.dtw_time << "s"
                  << " total=" << stats.total_time << "s"
                  << " speedup=" << stats.speedup
                  << std::endl;
    }

    out << '\n';
    out.close();
}
