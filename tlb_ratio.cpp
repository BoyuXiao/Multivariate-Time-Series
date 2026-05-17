#include "tlb_ratio.h"
#include "funcs.h"

#include <cmath>
#include <cctype>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <optional>
#include <sstream>
#include <unordered_set>
#include <vector>

namespace {

enum class TlbMethod {
    Mv,
    PMv,
    WebMv,
    Mc1,
    Mc2,
    Ti,
    Pc,
    Mc3,
};

struct QueryCache {
    std::vector<std::vector<double>> q_l, q_u, l_uq, u_lq;
    std::vector<double> mu_q, sigma_q;
    std::vector<double> u_mu_q, l_mu_q, u_sigma_q, l_sigma_q;
    std::vector<std::vector<std::vector<Cell>>> t_cells;
    std::vector<double> interval_mu, interval_sigma, dbl, dbr, dtl, dtr;
    std::vector<int> segment_num;
    std::vector<std::vector<std::vector<std::vector<double>>>> pc_bboxes;
};

struct TrainCache {
    std::vector<std::vector<double>> mu_a, sigma_a;
};

double elapsed_sec(clock_t start, clock_t end) {
    return static_cast<double>(end - start) / CLOCKS_PER_SEC;
}

void build_mv_envelopes(const std::vector<std::vector<double>>& series_t, int m, int k_dims, int r,
                        std::vector<std::vector<double>>& lower,
                        std::vector<std::vector<double>>& upper) {
    lower.assign(k_dims, std::vector<double>(m));
    upper.assign(k_dims, std::vector<double>(m));
    for (int k = 0; k < k_dims; ++k) {
        lower_upper_lemire(const_cast<std::vector<double>&>(series_t[k]), m, r, lower[k], upper[k]);
    }
}

void build_web_envelopes(const std::vector<std::vector<double>>& series_t, int m, int k_dims, int r,
                         std::vector<std::vector<double>>& l_uq,
                         std::vector<std::vector<double>>& u_lq) {
    l_uq.assign(k_dims, std::vector<double>(m));
    u_lq.assign(k_dims, std::vector<double>(m));
    for (int k = 0; k < k_dims; ++k) {
        lower_lemire(const_cast<std::vector<double>&>(series_t[k]), m, r, l_uq[k]);
        upper_lemire(const_cast<std::vector<double>&>(series_t[k]), m, r, u_lq[k]);
    }
}

void prepare_mc3_query(const std::vector<std::vector<double>>& query, int m, int k_dims, int r,
                       QueryCache& cache) {
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
                                 cache.mu_q, cache.sigma_q, cache.u_mu_q, cache.l_mu_q,
                                 cache.u_sigma_q, cache.l_sigma_q, cache.t_cells, cache.interval_mu,
                                 cache.interval_sigma, cache.segment_num, r);
    cache.dbl = dnn[0];
    cache.dbr = dnn[1];
    cache.dtl = dnn[2];
    cache.dtr = dnn[3];
}

void prepare_query(TlbMethod kind, const std::vector<std::vector<double>>& query,
                   const std::vector<std::vector<double>>& q_t, int m, int k_dims, int r,
                   double pc_q, QueryCache& cache) {
    switch (kind) {
        case TlbMethod::Mv:
        case TlbMethod::PMv:
            build_mv_envelopes(q_t, m, k_dims, r, cache.q_l, cache.q_u);
            break;
        case TlbMethod::WebMv:
            build_mv_envelopes(q_t, m, k_dims, r, cache.q_l, cache.q_u);
            build_web_envelopes(q_t, m, k_dims, r, cache.l_uq, cache.u_lq);
            break;
        case TlbMethod::Mc1:
            cache.mu_q.assign(m, 0.0);
            only_mu(const_cast<std::vector<std::vector<double>>&>(query), cache.mu_q);
            break;
        case TlbMethod::Mc2:
            cache.mu_q.assign(m, 0.0);
            cache.sigma_q.assign(m, 0.0);
            mu_sigma(const_cast<std::vector<std::vector<double>>&>(query), cache.mu_q, cache.sigma_q);
            break;
        case TlbMethod::Mc3:
            prepare_mc3_query(query, m, k_dims, r, cache);
            break;
        case TlbMethod::Pc:
            cache.pc_bboxes = LB_PC_Preprocess(query, r, k_dims, pc_q);
            break;
        case TlbMethod::Ti:
            break;
    }
}

bool needs_train_cache(TlbMethod kind) {
    return kind == TlbMethod::Mc1 || kind == TlbMethod::Mc2 || kind == TlbMethod::Mc3;
}

void prepare_trains(TlbMethod kind,
                    const std::vector<std::vector<std::vector<double>>>& train_set, int m,
                    int k_dims, TrainCache& cache) {
    if (!needs_train_cache(kind)) {
        return;
    }
    const int n_train = static_cast<int>(train_set.size());
    cache.mu_a.assign(n_train, std::vector<double>(m));
    cache.sigma_a.assign(n_train, std::vector<double>(m, 0.0));
    for (int j = 0; j < n_train; ++j) {
        if (kind == TlbMethod::Mc1) {
            only_mu(const_cast<std::vector<std::vector<double>>&>(train_set[j]), cache.mu_a[j]);
        } else {
            mu_sigma(const_cast<std::vector<std::vector<double>>&>(train_set[j]), cache.mu_a[j],
                     cache.sigma_a[j]);
        }
    }
}

double compute_lb_full(TlbMethod kind,
                       const std::vector<std::vector<double>>& query,
                       const std::vector<std::vector<double>>& train,
                       int train_idx,
                       int m,
                       int k_dims,
                       int r,
                       int ti_p,
                       QueryCache& qcache,
                       const TrainCache& tcache,
                       std::vector<std::vector<double>>& q_t,
                       std::vector<std::vector<double>>& a_t) {
    auto& q_mut = const_cast<std::vector<std::vector<double>>&>(query);
    auto& t_mut = const_cast<std::vector<std::vector<double>>&>(train);

    switch (kind) {
        case TlbMethod::Mv:
            return LB_MV(q_mut, t_mut, k_dims, m, r);
        case TlbMethod::PMv:
            return LB_P_MV(a_t, q_t, qcache.q_l, qcache.q_u, r);
        case TlbMethod::WebMv:
            return LB_WEB_MV(a_t, q_t, qcache.q_l, qcache.q_u, qcache.l_uq, qcache.u_lq, r);
        case TlbMethod::Mc1:
            return LB_MC1(q_mut, t_mut, const_cast<std::vector<double>&>(qcache.mu_q), k_dims, m, r);
        case TlbMethod::Mc2:
            return LB_MC2(q_mut, t_mut, const_cast<std::vector<double>&>(qcache.mu_q),
                          const_cast<std::vector<double>&>(qcache.sigma_q), k_dims, m, r);
        case TlbMethod::Ti:
            return LB_TI(query, train, r, ti_p);
        case TlbMethod::Pc:
            return LB_PC(q_mut, t_mut, qcache.pc_bboxes);
        case TlbMethod::Mc3:
            return LB_MC3(t_mut, qcache.u_mu_q, qcache.l_mu_q, qcache.u_sigma_q, qcache.l_sigma_q,
                          const_cast<std::vector<double>&>(tcache.mu_a[train_idx]),
                          const_cast<std::vector<double>&>(tcache.sigma_a[train_idx]), qcache.dbl,
                          qcache.dbr, qcache.dtl, qcache.dtr, qcache.t_cells, qcache.interval_mu,
                          qcache.interval_sigma, qcache.segment_num);
    }
    return 0.0;
}

std::string normalize_method_token(std::string s) {
    for (char& c : s) {
        if (c == ' ') {
            c = '_';
        }
        c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    }
    return s;
}

std::optional<TlbMethod> parse_method_token(const std::string& token) {
    const std::string t = normalize_method_token(token);
    if (t == "MV") {
        return TlbMethod::Mv;
    }
    if (t == "P_MV" || t == "PMV") {
        return TlbMethod::PMv;
    }
    if (t == "WEB_MV" || t == "WEBMV") {
        return TlbMethod::WebMv;
    }
    if (t == "MC1") {
        return TlbMethod::Mc1;
    }
    if (t == "MC2") {
        return TlbMethod::Mc2;
    }
    if (t == "TI") {
        return TlbMethod::Ti;
    }
    if (t == "PC") {
        return TlbMethod::Pc;
    }
    if (t == "MC3") {
        return TlbMethod::Mc3;
    }
    return std::nullopt;
}

std::vector<TlbMethod> all_tlb_methods() {
    return {TlbMethod::Mv, TlbMethod::PMv, TlbMethod::WebMv, TlbMethod::Mc1,
            TlbMethod::Mc2, TlbMethod::Ti, TlbMethod::Pc, TlbMethod::Mc3};
}

std::vector<TlbMethod> resolve_methods(const std::vector<std::string>& filter) {
    if (filter.empty()) {
        return all_tlb_methods();
    }
    std::vector<TlbMethod> out;
    for (const std::string& raw : filter) {
        std::stringstream ss(raw);
        std::string token;
        while (std::getline(ss, token, ',')) {
            if (token.empty()) {
                continue;
            }
            auto kind = parse_method_token(token);
            if (!kind) {
                std::cerr << "Unknown method in filter: " << token << std::endl;
                continue;
            }
            out.push_back(*kind);
        }
    }
    return out;
}

const char* method_name(TlbMethod kind) {
    switch (kind) {
        case TlbMethod::Mv:
            return "MV";
        case TlbMethod::PMv:
            return "P_MV";
        case TlbMethod::WebMv:
            return "WEB_MV";
        case TlbMethod::Mc1:
            return "MC1";
        case TlbMethod::Mc2:
            return "MC2";
        case TlbMethod::Ti:
            return "TI";
        case TlbMethod::Pc:
            return "PC";
        case TlbMethod::Mc3:
            return "MC3";
    }
    return "UNKNOWN";
}

void write_header_if_needed(std::ofstream& out, const std::string& path) {
    static std::unordered_set<std::string> done;
    if (done.count(path)) {
        return;
    }
    if (std::filesystem::exists(path) && std::filesystem::file_size(path) > 0) {
        done.insert(path);
        return;
    }
    out << "Name\tr_ratio\tr\tti_p\tpc_q\tMethod\tTLB_ratio\tPairs\tSkippedZeroDtw\t"
           "PrepTime\tLbTime\tDtwTime\n";
    done.insert(path);
}

TlbRatioStats run_one_method(
    TlbMethod kind,
    const std::vector<std::vector<std::vector<double>>>& test_set,
    const std::vector<std::vector<std::vector<double>>>& train_set,
    const std::vector<std::vector<std::vector<double>>>& test_t,
    const std::vector<std::vector<std::vector<double>>>& train_t,
    const std::vector<std::vector<double>>& dtw_mat,
    int m,
    int k_dims,
    int r,
    int ti_p,
    double pc_q) {
    const int n_test = static_cast<int>(test_set.size());
    const int n_train = static_cast<int>(train_set.size());
    const int total_pairs = n_test * n_train;

    TlbRatioStats stats;
    stats.name = method_name(kind);
    stats.pairs = total_pairs;

    TrainCache train_cache;
    clock_t t0 = clock();
    prepare_trains(kind, train_set, m, k_dims, train_cache);
    stats.prep_time = elapsed_sec(0, clock() - t0);

    double sum_ratio = 0.0;
    clock_t lb_clock = 0;

    for (int i = 0; i < n_test; ++i) {
        QueryCache qcache;
        t0 = clock();
        prepare_query(kind, test_set[i], test_t[i], m, k_dims, r, pc_q, qcache);
        stats.prep_time += elapsed_sec(0, clock() - t0);

        auto q_t_mut = const_cast<std::vector<std::vector<double>>&>(test_t[i]);
        for (int j = 0; j < n_train; ++j) {
            const double d = dtw_mat[i][j];
            if (d <= 1e-12) {
                ++stats.skipped_zero_dtw;
                continue;
            }
            auto a_t_mut = const_cast<std::vector<std::vector<double>>&>(train_t[j]);
            clock_t lb_start = clock();
            const double lb = compute_lb_full(kind, test_set[i], train_set[j], j, m, k_dims, r, ti_p,
                                              qcache, train_cache, q_t_mut, a_t_mut);
            lb_clock += clock() - lb_start;
            sum_ratio += lb / d;
        }
    }

    stats.lb_time = elapsed_sec(0, lb_clock);
    const int valid = total_pairs - stats.skipped_zero_dtw;
    stats.tlb_ratio = (valid > 0) ? (sum_ratio / valid) : 0.0;
    return stats;
}

}  // namespace

void run_tlb_ratio_test(const std::vector<std::vector<std::vector<double>>>& test_set,
                        const std::vector<std::vector<std::vector<double>>>& train_set,
                        const std::string& dataset_name,
                        const TlbRatioConfig& config) {
    if (test_set.empty() || train_set.empty()) {
        return;
    }

    const int m = static_cast<int>(test_set[0].size());
    const int k_dims = static_cast<int>(test_set[0][0].size());
    const int r = static_cast<int>(std::ceil(config.r_ratio * m));
    const int ti_p = config.ti_p;
    const double pc_q = config.pc_q;

    const std::vector<TlbMethod> methods = resolve_methods(config.methods_filter);
    if (methods.empty()) {
        std::cerr << "No valid methods to run (check methods_filter)." << std::endl;
        return;
    }

    const int n_test = static_cast<int>(test_set.size());
    const int n_train = static_cast<int>(train_set.size());

    std::vector<std::vector<std::vector<double>>> test_t(n_test);
    std::vector<std::vector<std::vector<double>>> train_t(n_train);
    for (int i = 0; i < n_test; ++i) {
        test_t[i] = transpose(test_set[i]);
    }
    for (int j = 0; j < n_train; ++j) {
        train_t[j] = transpose(train_set[j]);
    }

    std::vector<std::vector<double>> dtw_mat(n_test, std::vector<double>(n_train, 0.0));
    clock_t dtw_clock = clock();
    for (int i = 0; i < n_test; ++i) {
        for (int j = 0; j < n_train; ++j) {
            dtw_mat[i][j] =
                dtw(const_cast<std::vector<std::vector<double>>&>(test_set[i]),
                    const_cast<std::vector<std::vector<double>>&>(train_set[j]), k_dims, m, r);
        }
    }
    const double dtw_time = elapsed_sec(dtw_clock, clock());

    if (config.result_file.empty()) {
        std::cerr << "TlbRatioConfig.result_file is empty." << std::endl;
        return;
    }
    std::ofstream out(config.result_file, std::ios::app);
    if (!out.is_open()) {
        std::cerr << "Failed to open: " << config.result_file << std::endl;
        return;
    }
    write_header_if_needed(out, config.result_file);

    std::cout << dataset_name << " TLB ratio  r_ratio=" << config.r_ratio << " r=" << r
              << " ti_p=" << ti_p << " pc_q=" << pc_q << " pairs=" << (n_test * n_train)
              << " methods=";
    for (size_t i = 0; i < methods.size(); ++i) {
        if (i) {
            std::cout << ',';
        }
        std::cout << method_name(methods[i]);
    }
    std::cout << std::endl;
    std::cout << "  DTW(all pairs): " << dtw_time << "s" << std::endl;

    for (TlbMethod kind : methods) {
        TlbRatioStats stats =
            run_one_method(kind, test_set, train_set, test_t, train_t, dtw_mat, m, k_dims, r, ti_p,
                           pc_q);
        stats.dtw_time = dtw_time;
        std::cout << "  " << stats.name << " TLB=" << stats.tlb_ratio << " prep=" << stats.prep_time
                  << "s lb=" << stats.lb_time << "s" << std::endl;

        out << dataset_name << '\t' << config.r_ratio << '\t' << r << '\t' << ti_p << '\t' << pc_q
            << '\t' << stats.name << '\t' << stats.tlb_ratio << '\t' << stats.pairs << '\t'
            << stats.skipped_zero_dtw << '\t' << stats.prep_time << '\t' << stats.lb_time << '\t'
            << stats.dtw_time << '\n';
    }
    out << '\n';
    out.close();
}
