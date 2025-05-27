#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>
#include "funcs.h"

// 计算两点之间的距离，假设两点是向量
double distance(const std::vector<double>& p1, const std::vector<double>& p2) {
    double sum = 0.0;
    for (size_t i = 0; i < p1.size(); ++i) {
        double diff = p1[i] - p2[i];
        sum += diff * diff;
    }
    // return std::sqrt(sum);
    return sum;
}
double calculateLB_TI(const std::vector<std::vector<double>>& query, const std::vector<std::vector<double>>& candidate, int windowSize, int period) {
    // int qLen = static_cast<int>(query.size());
    // int cLen = static_cast<int>(candidate.size());
    // int dim = static_cast<int>(query[0].size());
    // double lb_ti = 0.0;
    //
    // std::vector<double> dxx(qLen - 1);
    // for (int i = 0; i < qLen - 1; ++i) {
    //     dxx[i] = distance(query[i], query[i + 1]);
    // }
    //
    // std::vector<std::vector<double>> upperBounds(qLen, std::vector<double>(2 * windowSize + 1, 0.0));
    // std::vector<std::vector<double>> lowerBounds(qLen, std::vector<double>(2 * windowSize + 1, 0.0));
    //
    // for (int t = 0; t < qLen; ++t) {
    //     int startIdx = (t > windowSize) ? 0 : windowSize - t;
    //     int lw = std::max(0, t - windowSize);
    //     int tp = std::min(t + windowSize + 1, cLen);
    //
    //     if (t % period == 0) {
    //         std::vector<double> dxyInit;
    //         for (int i = lw; i < tp; ++i) {
    //             dxyInit.push_back(distance(query[t], candidate[i]));
    //         }
    //         for (int i = 0; i < tp - lw; ++i) {
    //             upperBounds[t][startIdx + i] = dxyInit[i];
    //             lowerBounds[t][startIdx + i] = dxyInit[i];
    //         }
    //         lb_ti += *std::min_element(dxyInit.begin(), dxyInit.end());
    //     } else {
    //         int lr = (t < windowSize) ? 0 : t - windowSize;
    //         int ur = (cLen - 1 < t + windowSize) ? cLen - 1 : t + windowSize;
    //         double thisdxx = dxx[t - 1];
    //         int startIdx_lr = startIdx - lr + 1;
    //         int t_1 = t - 1;
    //         int idx = ur - lr - 1;
    //
    //         if (t + windowSize <= cLen - 1) {
    //             for (int i = lr; i < ur; ++i) {
    //                 upperBounds[t][startIdx + i - lr] = upperBounds[t_1][startIdx_lr + i] + thisdxx;
    //                 if (lowerBounds[t_1][startIdx_lr + i] > thisdxx) {
    //                     lowerBounds[t][startIdx + i - lr] = lowerBounds[t_1][startIdx_lr + i] - thisdxx;
    //                 } else if (thisdxx < upperBounds[t_1][startIdx_lr + i]) {
    //                     lowerBounds[t][startIdx + i - lr] = 0;
    //                 } else {
    //                     lowerBounds[t][startIdx + i - lr] = thisdxx - upperBounds[t_1][startIdx_lr + i];
    //                 }
    //             }
    //             double temp = distance(query[t], candidate[ur]);
    //             upperBounds[t][startIdx + idx + 1] = temp;
    //             lowerBounds[t][startIdx + idx + 1] = temp;
    //             lb_ti += *std::min_element(lowerBounds[t].begin() + startIdx, lowerBounds[t].begin() + startIdx + idx + 2);
    //         } else {
    //             for (int i = lr; i <= ur; ++i) {
    //                 upperBounds[t][startIdx + i - lr] = upperBounds[t_1][startIdx_lr + i] + thisdxx;
    //                 if (lowerBounds[t_1][startIdx_lr + i] > thisdxx) {
    //                     lowerBounds[t][startIdx + i - lr] = lowerBounds[t_1][startIdx_lr + i] - thisdxx;
    //                 } else if (thisdxx < upperBounds[t_1][startIdx_lr + i]) {
    //                     lowerBounds[t][startIdx + i - lr] = 0;
    //                 } else {
    //                     lowerBounds[t][startIdx + i - lr] = thisdxx - upperBounds[t_1][startIdx_lr + i];
    //                 }
    //             }
    //             lb_ti += *std::min_element(lowerBounds[t].begin() + startIdx, lowerBounds[t].begin() + startIdx + idx + 2);
    //         }
    //     }
    // }

    int qLen = static_cast<int>(query.size());
    int cLen = static_cast<int>(candidate.size());
    double lb_ti = 0.0;

    std::vector<double> dxx(qLen - 1);
    for (int i = 0; i < qLen - 1; ++i) {
        dxx[i] = distance(query[i], query[i + 1]);
    }

    std::vector<std::vector<double>> upperBounds(qLen, std::vector<double>(2 * windowSize + 1, 0.0));
    std::vector<std::vector<double>> lowerBounds(qLen, std::vector<double>(2 * windowSize + 1, 0.0));

    // 初始化第一个点的上下界
    int startIdx = windowSize;
    int lw = 0;
    int tp = std::min(windowSize + 1, cLen);
    std::vector<double> dxyInit;
    for (int i = lw; i < tp; ++i) {
        dxyInit.push_back(distance(query[0], candidate[i]));
    }
    for (int i = 0; i < tp - lw; ++i) {
        upperBounds[0][startIdx + i] = dxyInit[i];
        lowerBounds[0][startIdx + i] = dxyInit[i];
    }
    lb_ti += *std::min_element(dxyInit.begin(), dxyInit.end());

    for (int t = 1; t < qLen; ++t) {
        startIdx = (t > windowSize) ? 0 : windowSize - t;
        int lr = (t < windowSize) ? 0 : t - windowSize;
        int ur = (cLen - 1 < t + windowSize) ? cLen - 1 : t + windowSize;
        double thisdxx = dxx[t - 1];
        int startIdx_lr = startIdx - lr + 1;
        int t_1 = t - 1;
        int idx = ur - lr - 1;

        if (t + windowSize <= cLen - 1) {
            for (int i = lr; i < ur; ++i) {
                upperBounds[t][startIdx + i - lr] = upperBounds[t_1][startIdx_lr + i] + thisdxx;
                if (lowerBounds[t_1][startIdx_lr + i] > thisdxx) {
                    lowerBounds[t][startIdx + i - lr] = lowerBounds[t_1][startIdx_lr + i] - thisdxx;
                } else if (thisdxx < upperBounds[t_1][startIdx_lr + i]) {
                    lowerBounds[t][startIdx + i - lr] = 0;
                } else {
                    lowerBounds[t][startIdx + i - lr] = thisdxx - upperBounds[t_1][startIdx_lr + i];
                }
            }
            double temp = distance(query[t], candidate[ur]);
            upperBounds[t][startIdx + idx + 1] = temp;
            lowerBounds[t][startIdx + idx + 1] = temp;
            lb_ti += *std::min_element(lowerBounds[t].begin() + startIdx, lowerBounds[t].begin() + startIdx + idx + 2);
        } else {
            for (int i = lr; i <= ur; ++i) {
                upperBounds[t][startIdx + i - lr] = upperBounds[t_1][startIdx_lr + i] + thisdxx;
                if (lowerBounds[t_1][startIdx_lr + i] > thisdxx) {
                    lowerBounds[t][startIdx + i - lr] = lowerBounds[t_1][startIdx_lr + i] - thisdxx;
                } else if (thisdxx < upperBounds[t_1][startIdx_lr + i]) {
                    lowerBounds[t][startIdx + i - lr] = 0;
                } else {
                    lowerBounds[t][startIdx + i - lr] = thisdxx - upperBounds[t_1][startIdx_lr + i];
                }
            }
            lb_ti += *std::min_element(lowerBounds[t].begin() + startIdx, lowerBounds[t].begin() + startIdx + idx + 2);
        }
    }

    return sqrt(lb_ti);
}

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>

// // 计算两点之间的欧几里得距离
// double distance(const std::vector<double>& p1, const std::vector<double>& p2) {
//     double sum = 0.0;
//     for (size_t i = 0; i < p1.size(); ++i) {
//         double diff = p1[i] - p2[i];
//         sum += diff * diff;
//     }
//     // return std::sqrt(sum);
//     return sum;
// }
//
// // 计算 LB_TI 下界
// double calculateLB_TI(const std::vector<std::vector<double>>& C, const std::vector<std::vector<double>>& Q, int r) {
//     int cLen = static_cast<int>(C.size());
//     int qLen = static_cast<int>(Q.size());
//     double lb_ti = 0.0;
//
//     // 计算 C 和 Q 的相邻点距离
//     std::vector<double> delta_C(cLen - 1);
//     std::vector<double> delta_Q(qLen - 1);
//     for (int i = 0; i < cLen - 1; ++i) {
//         delta_C[i] = distance(C[i], C[i + 1]);
//     }
//     for (int i = 0; i < qLen - 1; ++i) {
//         delta_Q[i] = distance(Q[i], Q[i + 1]);
//     }
//
//     // 计算从 C[0] 到 Q[0:r] 的距离，存储到 prevLB 和 prevUB
//     std::vector<double> prevLB;
//     std::vector<double> prevUB;
//     for (int j = 0; j < std::min(r, qLen); ++j) {
//         double dist = distance(C[0], Q[j]);
//         prevLB.push_back(dist);
//         prevUB.push_back(dist);
//     }
//     lb_ti = *std::min_element(prevLB.begin(), prevLB.end());
//
//     for (int i = 1; i < cLen; ++i) {
//         std::vector<double> L;
//         std::vector<double> U;
//         const auto& c = C[i];
//
//         // 遍历 Q 中位于 c 的窗口内（不包括窗口顶部）的点
//         for (int j = std::max(0, i - r); j < std::min(i + r - 1, qLen); ++j) {
//             int prevIndex = j - (i - r);
//             // 计算上界
//             U.push_back(prevUB[prevIndex] + delta_C[i - 1]);
//             // 计算下界
//             if (prevLB[prevIndex] > delta_C[i - 1]) {
//                 L.push_back(prevLB[prevIndex] - delta_C[i - 1]);
//             } else if (delta_C[i - 1] < prevUB[prevIndex]) {
//                 L.push_back(0);
//             } else {
//                 L.push_back(delta_C[i - 1] - prevUB[prevIndex]);
//             }
//         }
//
//         // 计算窗口顶部点的边界
//         if (i + r - 1 < qLen) {
//             int prevIndex = i + r - 1 - (i - r);
//             // 计算窗口顶部点的下界
//             if (prevLB[prevIndex] > delta_Q[i + r - 1]) {
//                 L.push_back(prevLB[prevIndex] - delta_Q[i + r - 1]);
//             } else if (delta_Q[i + r - 1] < prevUB[prevIndex]) {
//                 L.push_back(0);
//             } else {
//                 L.push_back(delta_Q[i + r - 1] - prevUB[prevIndex]);
//             }
//             // 计算窗口顶部点的上界
//             U.push_back(prevUB[prevIndex] + delta_Q[i + r - 1]);
//         }
//
//         lb_ti += *std::min_element(L.begin(), L.end());
//         prevLB = L;
//         prevUB = U;
//     }
//
//     return sqrt(lb_ti);
// }