#ifndef NN1_LB_H
#define NN1_LB_H

#include <string>
#include <vector>

struct Nn1LbStats
{
    std::string name;
    int pruned = 0;
    int dtw_calls = 0;
    double prep_time = 0.0;  // 该方法专用下界预处理（不含 transpose）
    double lb_time = 0.0;    // 下界计算
    double dtw_time = 0.0;   // DTW 计算
    double total_time = 0.0; // prep + lb + dtw
    double dtw_baseline_time = 0.0;
    double speedup = 0.0;

    double prune_rate() const
    {
        const int comparisons = pruned + dtw_calls;
        return comparisons == 0 ? 0.0 : static_cast<double>(pruned) / comparisons;
    }
};

struct Nn1Config
{
    double r_ratio = 0.01;
    bool use_early_stopping = true;
    // 一次测试运行的结果 TSV（由 main 写入 res/nn1_YYYYMMDD_HHMMSS.txt）
    std::string result_file;
    // 逗号分隔方法名，如 "MC1,MC2,MC3"；空则跑全部默认方法
    std::vector<std::string> methods_filter;
};

void run_nn1_lb_test(const std::vector<std::vector<std::vector<double>>> &test_set,
                     const std::vector<std::vector<std::vector<double>>> &train_set,
                     const std::string &dataset_name,
                     const Nn1Config &config = Nn1Config{});

#endif
