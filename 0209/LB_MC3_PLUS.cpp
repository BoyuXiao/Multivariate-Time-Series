// #include <iostream>
// #include <vector>
// #include <cmath>
// #include <algorithm>
// #include "funcs.h"
// using namespace std;
//
//
// // 拆分测试集 Q 为 Q1 和 Q2
// void splitQ(const std::vector<std::vector<double>>& Q,
//             std::vector<std::vector<double>>& Q1,
//             std::vector<std::vector<double>>& Q2,
//             std::vector<std::vector<int>>& split_indices,int r) {
//     int m = Q.size();    // 维度数量
//     int K = Q[0].size(); // 每个维度的变量数量
//     int split_size = K / 2;
//     vector<int> start(m);
//     vector<int> end(m);
//     Q1.clear();
//     Q2.clear();
//     split_indices.clear();
//
//     for (int i = 0; i < m; ++i) {
//         start[i] = MAX(0,i-r);
//         end[i] = MIN(m-1,i+r);
//     }
//
//     vector<pair<double, int>> means(K, {0, 0});
//     for (int i = 0; i < K; i++) {
//         means[i].second = i;
//     }
//     vector<int> count(K, 0);  // 用于记录每个 k 参与累加的元素个数
//     for (int i = 0; i < m; ++i) {
//         for (int j = start[i]; j <= end[i]; ++j) {
//             for (int k = 0; k < K; ++k) {
//                 means[k].first += Q[j][k];
//                 count[k]++;  // 累加元素个数
//             }
//         }
//     }
//     // 计算均值
//     for (int k = 0; k < K; ++k) {
//         if (count[k] != 0) {
//             means[k].first /= count[k];
//         }
//     }
//
//         // 按均值降序排序
//         sort(means.rbegin(), means.rend());
//
//         // 记录拆分索引
//         vector<int> indices1, indices2;
//         for (int j = 0; j < split_size; ++j) {
//             indices1.push_back(means[j].second);
//         }
//         for (int j = split_size; j < K; ++j) {
//             indices2.push_back(means[j].second);
//         }
//         split_indices.push_back(indices1); // 记录当前维度的拆分规则
//
//         // 构建 Q1 和 Q2 的当前维度
//         vector<double> q1_dim, q2_dim;
//         for (int idx : indices1) {
//             q1_dim.push_back(Q[i][idx]);
//         }
//         for (int idx : indices2) {
//             q2_dim.push_back(Q[i][idx]);
//         }
//
//         Q1.push_back(q1_dim);
//         Q2.push_back(q2_dim);
//     }
// }


#include <vector>
#include <algorithm>
#include <utility>
#include <cmath> // 用于 MAX 和 MIN 函数

// 辅助函数：计算一个窗口内的均值
double calculateWindowMean(const std::vector<std::vector<double>>& Q, int i, int start, int end, int j) {
    double sum = 0.0;
    int count = 0;
    for (int k = start; k <= end; ++k) {
        sum += Q[k][j]; // Q[k][j] 是第 k 个时间点的第 j 个变量
        count++;
    }
    return sum / count;
}

// 拆分测试集 Q 为 Q1 和 Q2
void splitQ(const std::vector<std::vector<double>>& Q,
            std::vector<std::vector<double>>& Q1,
            std::vector<std::vector<double>>& Q2,
            std::vector<std::vector<int>>& split_indices,
            int r) {
    int m = Q.size();    // 时间点数量
    int K = Q[0].size(); // 每个时间点的变量数量
    int split_size = K / 2;

    // 清空输出
    Q1.clear();
    Q2.clear();
    split_indices.clear();

    // 计算每个时间点 i 的 start 和 end
    std::vector<int> start(m), end(m);
    for (int i = 0; i < m; ++i) {
        start[i] = std::max(0, i - r);
        end[i] = std::min(m - 1, i + r);
    }

    // 对每个时间点 i 进行拆分
    for (int i = 0; i < m; ++i) {
        // 提取当前时间点 i 的窗口内的变量并计算均值
        std::vector<std::pair<double, int>> means;
        for (int j = 0; j < K; ++j) {
            double mean = calculateWindowMean(Q, i, start[i], end[i], j); // 计算窗口均值
            means.emplace_back(mean, j);
        }

        // 按均值降序排序
        std::sort(means.rbegin(), means.rend());

        // 记录拆分索引
        std::vector<int> indices1, indices2;
        for (int j = 0; j < split_size; ++j) {
            indices1.push_back(means[j].second);
        }
        for (int j = split_size; j < K; ++j) {
            indices2.push_back(means[j].second);
        }
        split_indices.push_back(indices1); // 记录当前时间点的拆分规则

        // 构建 Q1 和 Q2 的当前时间点
        std::vector<double> q1_dim, q2_dim;
        for (int idx : indices1) {
            q1_dim.push_back(Q[i][idx]);
        }
        for (int idx : indices2) {
            q2_dim.push_back(Q[i][idx]);
        }

        Q1.push_back(q1_dim);
        Q2.push_back(q2_dim);
    }
}

// 根据 Q 的拆分规则拆分训练集 A
void splitA(const std::vector<std::vector<double>>& A,
            const std::vector<std::vector<int>>& split_indices,
            std::vector<std::vector<double>>& A1,
            std::vector<std::vector<double>>& A2) {
    int m = A.size();    // 维度数量
    int K = A[0].size(); // 每个维度的变量数量

    // 清空输出
    A1.clear();
    A2.clear();

    // 对每个维度 i 进行拆分
    for (int i = 0; i < m; ++i) {
        const std::vector<int>& indices1 = split_indices[i]; // 获取当前维度的拆分规则

        // 构建 A1 和 A2 的当前维度
        std::vector<double> a1_dim, a2_dim;
        for (int j = 0; j < K; ++j) {
            if (std::find(indices1.begin(), indices1.end(), j) != indices1.end()) {
                a1_dim.push_back(A[i][j]); // 属于 Q1 的变量
            } else {
                a2_dim.push_back(A[i][j]); // 属于 Q2 的变量
            }
        }

        A1.push_back(a1_dim);
        A2.push_back(a2_dim);
    }
}