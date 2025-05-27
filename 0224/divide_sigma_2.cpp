
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>

// 计算一组数据的标准差
double calculate_std_dev2(const std::vector<double>& data) {
    if (data.empty()) return 0.0;
    double sum = 0.0;
    for (double val : data) {
        sum += val;
    }
    double mean = sum / data.size();

    double variance = 0.0;
    for (double val : data) {
        variance += std::pow(val - mean, 2);
    }
    variance /= data.size();
    return std::sqrt(variance);
}

// 计算度量值
double calculate_metric2(const std::vector<double>& I1, const std::vector<double>& I2) {
    double sigma1 = calculate_std_dev2(I1);
    double sigma2 = calculate_std_dev2(I2);
    int size1 = I1.size();
    int size2 = I2.size();
    return size1 * std::pow(sigma1, 2) + size2 * std::pow(sigma2, 2);
}

// 分割Q矩阵的函数
void splitQ_min_metric(const std::vector<std::vector<double>>& Q,
                       std::vector<std::vector<double>>& Q1,
                       std::vector<std::vector<double>>& Q2,
                       std::vector<int>& split_indices) {
    int m = Q.size();
    int K = Q[0].size();

    double min_metric = std::numeric_limits<double>::max();
    std::vector<int> best_split_indices;

    // 生成所有可能的分割方式
    for (int mask = 0; mask < (1 << K); ++mask) {
        std::vector<int> current_split_indices;
        for (int k = 0; k < K; ++k) {
            if (mask & (1 << k)) {
                current_split_indices.push_back(k);
            }
        }

        double total_metric = 0.0;
        for (int i = 0; i < m; ++i) {
            std::vector<double> I1, I2;
            for (int k = 0; k < K; ++k) {
                if (std::find(current_split_indices.begin(), current_split_indices.end(), k) != current_split_indices.end()) {
                    I1.push_back(Q[i][k]);
                } else {
                    I2.push_back(Q[i][k]);
                }
            }
            total_metric += calculate_metric2(I1, I2);
        }

        if (total_metric < min_metric) {
            min_metric = total_metric;
            best_split_indices = current_split_indices;
        }
    }

    split_indices = best_split_indices;

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