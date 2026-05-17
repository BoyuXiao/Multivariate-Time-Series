#ifndef TLB_RATIO_H
#define TLB_RATIO_H

#include <string>
#include <vector>

/** 经典 TLB 指标：对所有 (query, train) 对计算 mean(LB/DTW)，使用完整下界（非 1NN early stopping）。 */
struct TlbRatioStats {
    std::string name;
    double tlb_ratio = 0.0;  // mean(LB/DTW)
    int pairs = 0;
    int skipped_zero_dtw = 0;
    double prep_time = 0.0;
    double lb_time = 0.0;
    double dtw_time = 0.0;  // 全对 DTW（各方法共享，仅首行有意义）
};

struct TlbRatioConfig {
    double r_ratio = 0.01;
    int ti_p = 5;       // LB_TI 的 P
    double pc_q = 2.0;  // LB_PC 的 Q
    std::string result_file;
    /** 为空则跑全部方法；否则只跑列出的方法名，如 MC1、MV */
    std::vector<std::string> methods_filter;
};

void run_tlb_ratio_test(const std::vector<std::vector<std::vector<double>>>& test_set,
                        const std::vector<std::vector<std::vector<double>>>& train_set,
                        const std::string& dataset_name,
                        const TlbRatioConfig& config = TlbRatioConfig{});

#endif
