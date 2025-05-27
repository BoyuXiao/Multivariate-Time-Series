
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <random>

double calculate_std_dev(const std::vector<double>& data) {
    if (data.empty()) return 0.0;
    double sum = 0.0;
    for (double val : data) {
        sum += val;
    }
    double mean = sum / data.size();

    double variance = 0.0;
    for (double val : data) {
        variance += (val - mean)*(val - mean);
    }
    variance /= data.size();
    return std::sqrt(variance);
}

// 计算度量值
double calculate_metric(const std::vector<std::vector<double>>& Q, const std::vector<int>& split_indices) {
    int m = Q.size();
    int K = Q[0].size();
    double total_metric = 0.0;
    for (int i = 0; i < m; ++i) {
        std::vector<double> I1, I2;
        for (int k = 0; k < K; ++k) {
            if (std::find(split_indices.begin(), split_indices.end(), k) != split_indices.end()) {
                I1.push_back(Q[i][k]);
            } else {
                I2.push_back(Q[i][k]);
            }
        }
        double sigma1 = calculate_std_dev(I1);
        double sigma2 = calculate_std_dev(I2);
        int size1 = I1.size();
        int size2 = I2.size();
        total_metric += size1 * sigma1*sigma1 + size2 * sigma2*sigma2;
    }
    return total_metric;
}

// 生成一个邻域解
std::vector<int> generate_neighbor(const std::vector<int>& split_indices, int K) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, K - 1);

    std::vector<int> neighbor = split_indices;
    int random_index = dis(gen);
    auto it = std::find(neighbor.begin(), neighbor.end(), random_index);
    if (it != neighbor.end()) {
        neighbor.erase(it);
    } else {
        neighbor.push_back(random_index);
    }
    return neighbor;
}

// 模拟退火算法
std::vector<int> simulated_annealing(const std::vector<std::vector<double>>& Q, double initial_temp, double cooling_rate, int max_iterations) {
    int K = Q[0].size();
    std::mt19937 gen(42);
    std::uniform_int_distribution<> dis(0, K - 1);

    // 随机初始化一个分割方案
    std::vector<int> current_split;
    for (int i = 0; i < K / 2; ++i) {
        int index;
        do {
            index = dis(gen);
        } while (std::find(current_split.begin(), current_split.end(), index) != current_split.end());
        current_split.push_back(index);
    }

    double current_metric = calculate_metric(Q, current_split);
    std::vector<int> best_split = current_split;
    double best_metric = current_metric;

    double temp = initial_temp;
    for (int iter = 0; iter < max_iterations; ++iter) {
        std::vector<int> neighbor = generate_neighbor(current_split, K);
        double neighbor_metric = calculate_metric(Q, neighbor);

        double delta = neighbor_metric - current_metric;
        if (delta < 0 || std::exp(-delta / temp) > ((double) rand() / RAND_MAX)) {
            current_split = neighbor;
            current_metric = neighbor_metric;
            if (current_metric < best_metric) {
                best_split = current_split;
                best_metric = current_metric;
            }
        }

        temp *= cooling_rate;
    }

    return best_split;
}

// 分割Q矩阵的函数
void splitQ_min_metric_heuristic(const std::vector<std::vector<double>>& Q,
                                 std::vector<std::vector<double>>& Q1,
                                 std::vector<std::vector<double>>& Q2,
                                 std::vector<int>& split_indices) {
    double initial_temp = 100.0;
    double cooling_rate = 0.99;
    int max_iterations = 1000;

    split_indices = simulated_annealing(Q, initial_temp, cooling_rate, max_iterations);

    int m = Q.size();
    int K = Q[0].size();
    for (int i = 0; i < m; ++i) {
        std::vector<double> q1_dim, q2_dim;
        for (int k = 0; k < K; ++k) {
            if (std::find(split_indices.begin(), split_indices.end(), k) != split_indices.end()) {
                q1_dim.push_back(Q[i][k]);
            } else {
                q2_dim.push_back(Q[i][k]);
            }
        }
        Q1[i] = q1_dim;
        Q2[i] = q2_dim;
    }
}