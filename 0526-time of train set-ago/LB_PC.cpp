#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>
#include <unordered_map>
#include "funcs.h"
#include <map>
#include <numeric>

double getLB_oneQR(const vector<vector<double>>& Q,
                  const vector<vector<double>>& R,
                  const vector<vector<vector<vector<double>>>>& bboxes,
                  double dist) {

    double LB_sum = 0;
    const int dim = Q[0].size(); // 获取维度

    // 遍历参考序列的每个点
    for(int idr = 0; idr < R.size(); ++idr) {
        const auto& r = R[idr];

        // 获取当前点对应的所有边界框
        const auto& cur_boxes = bboxes[idr];
        vector<double> bounds;

        // 遍历每个边界框
        for(const auto& box : cur_boxes) {
            const auto& l = box[0]; // 下界
            const auto& u = box[1]; // 上界

            double temp = 0;

            // 计算每个维度的距离贡献
            for(int d = 0; d < dim; ++d) {
                if(r[d] > u[d]) {
                    temp += pow(r[d] - u[d], 2);
                } else if(r[d] < l[d]) {
                    temp += pow(l[d] - r[d], 2);
                }
                // 在边界框内时不贡献距离
            }

            bounds.push_back(temp);
        }

        // 取最小边界贡献
        LB_sum += *min_element(bounds.begin(), bounds.end());

        // 提前终止检查
        if(sqrt(LB_sum) >= dist) {
            return sqrt(LB_sum);
        }
    }

    return sqrt(LB_sum);
}


vector<vector<vector<double>>> findboxes(const vector<vector<double>>& awindow,
                                        int K,
                                        double Q) {

    const int dims = awindow[0].size();
    vector<double> overall_ls(dims), overall_us(dims);
    vector<int> cells(dims);
    vector<double> celllens(dims);

    // Step 1: 计算每个维度的全局最小最大值
    for (int d = 0; d < dims; ++d) {
        double min_val = numeric_limits<double>::max();
        double max_val = -numeric_limits<double>::max();
        for (const auto& point : awindow) {
            if (point[d] < min_val) min_val = point[d];
            if (point[d] > max_val) max_val = point[d];
        }
        overall_ls[d] = min_val;
        overall_us[d] = max_val;
    }

    // Step 2: 计算分箱参数
    for (int d = 0; d < dims; ++d) {
        double range = overall_us[d] - overall_ls[d];
        cells[d] = 1 + static_cast<int>(range * Q);
        celllens[d] = range / cells[d] + 1e-8; // 添加小量防止除以零
    }

    // Step 3: 分配点到网格单元
    map<string, vector<vector<double>>> cellMembers;
    for (const auto& e : awindow) {
        vector<int> coords(dims);
        for (int d = 0; d < dims; ++d) {
            double pos = (e[d] - overall_ls[d]) / celllens[d];
            coords[d] = static_cast<int>(pos);
        }

        // 生成唯一字符串键
        ostringstream key;
        for (int c : coords) {
            key << c << ",";
        }

        cellMembers[key.str()].push_back(e);
    }

    // Step 4: 生成初始边界框
    vector<vector<vector<double>>> bboxes;
    for (const auto& [key, points] : cellMembers) {
        vector<double> l(dims, numeric_limits<double>::max());
        vector<double> u(dims, -numeric_limits<double>::max());

        for (const auto& p : points) {
            for (int d = 0; d < dims; ++d) {
                if (p[d] < l[d]) l[d] = p[d];
                if (p[d] > u[d]) u[d] = p[d];
            }
        }
        bboxes.push_back({l, u});
    }

    // Step 5: 合并多余边界框
    if (bboxes.size() > K && K > 0) {
        vector<vector<vector<double>>> remaining(bboxes.begin() + K - 1, bboxes.end());
        vector<double> combinedL(dims, numeric_limits<double>::max());
        vector<double> combinedU(dims, -numeric_limits<double>::max());

        for (const auto& b : remaining) {
            for (int d = 0; d < dims; ++d) {
                if (b[0][d] < combinedL[d]) combinedL[d] = b[0][d];
                if (b[1][d] > combinedU[d]) combinedU[d] = b[1][d];
            }
        }

        bboxes.resize(K - 1);
        bboxes.push_back({combinedL, combinedU});
    }

    return bboxes;
}


tuple<double, vector<vector<vector<vector<double>>>>> DTWwbbox(
    const vector<vector<double>>& s1,
    const vector<vector<double>>& s2,
    int windowSize,
    int K,
    double Q) {
    const int len1 = s1.size();
    const int len2 = s2.size();
    const int w = max(windowSize, abs(len1 - len2));

    // 初始化DTW表 (len1+1)x(len2+1) 避免边界检查
    vector<vector<double>> DTW(len1, vector<double>(len2, INF));
    vector<vector<vector<vector<double>>>> bboxes;

    // 设置初始条件
    DTW[0][0] = distance(s1[0], s2[0]);

    for (int i = 0; i < len1; ++i) {
        int left = std::max(0, i - w);
        int right = std::min(len2-1, i + w);
        for (int j = left; j <= right; ++j) {
            if (i == 0 && j == 0) continue;
            double dist = distance(s1[i], s2[j]);
            double min_val = INT_MAX;
            if (i > 0) min_val = std::min(min_val, DTW[i - 1][j]);
            if (j > 0) min_val = std::min(min_val, DTW[i][j - 1]);
            if (i > 0 && j > 0) min_val = std::min(min_val, DTW[i - 1][j - 1]);
            DTW[i][j] = dist + min_val;
        }
        bboxes.push_back(findboxes(s1, K, Q));
    }


    return {sqrt(DTW[len1][len2]), bboxes};
}

// tuple<double, int, int, double, int> DTWDistanceWindowLB_Ordered_LBPC_(
//     int K,
//     double Q,
//     const vector<vector<double>>& query,
//     const vector<vector<vector<double>>>& references,
//     int W,
//     double TH = 1.0)
// {
//     // 初始化统计量
//     int skip = 0;
//     int cluster_cals = 0;
//     double coretime = 0.0;
//
//     const int ql = query.size();
//     const int dim = query[0].size();
//     vector<vector<vector<double>>> bounds(ql);
//
//     // Step 1: 计算查询序列的边界框
//     for (int idx = 0; idx < ql; ++idx) {
//         int start = max(0, idx - W);
//         int end = min(ql-1, idx + W) + 1; // 包含右端点
//
//         vector<vector<double>> segment;
//         for (int i = start; i < end; ++i) {
//             segment.push_back(query[i]);
//         }
//
//         vector<double> l(dim, numeric_limits<double>::max());
//         vector<double> u(dim, -numeric_limits<double>::max());
//
//         for (const auto& p : segment) {
//             for (int d = 0; d < dim; ++d) {
//                 if (p[d] < l[d]) l[d] = p[d];
//                 if (p[d] > u[d]) u[d] = p[d];
//             }
//         }
//         bounds[idx] = {l, u};
//     }
//
//     // Step 2: 计算初步下界并排序
//     auto M0LBs = getLB_oneQ_qbox(query, references, bounds);
//
//     // 生成排序索引
//     vector<size_t> LBSortedIndex(references.size());
//     iota(LBSortedIndex.begin(), LBSortedIndex.end(), 0);
//     sort(LBSortedIndex.begin(), LBSortedIndex.end(),
//         [&M0LBs](size_t i1, size_t i2) { return M0LBs[i1] < M0LBs[i2]; });
//
//     int predId = LBSortedIndex[0];
//
//     // Step 3: 计算初始DTW和边界框
//     auto [dist, bboxes] = DTWwbbox(query, references[predId], W, K, Q);
//
//     // Step 4: 二次筛选过程
//     for (size_t x = 1; x < LBSortedIndex.size(); ++x) {
//         int thisrefid = LBSortedIndex[x];
//         double currentLB = M0LBs[thisrefid];
//
//         if (currentLB >= dist) {
//             skip += (LBSortedIndex.size() - x);
//             break;
//         }
//
//         if (currentLB >= dist - TH * dist) {
//             // 使用聚类下界
//             double c_lb = getLB_oneQR(query, references[thisrefid], bboxes, dist);
//             cluster_cals++;
//
//             if (c_lb < dist) {
//                 double dist2 = DTW_a(query, references[thisrefid], W, dist);
//                 if (dist2 < dist) {
//                     dist = dist2;
//                     predId = thisrefid;
//                 }
//             } else {
//                 skip++;
//             }
//         } else {
//             // 直接计算DTW
//             double dist2 = DTW_a(query, references[thisrefid], W, dist);
//             if (dist2 < dist) {
//                 dist = dist2;
//                 predId = thisrefid;
//             }
//         }
//     }
//
//     return make_tuple(dist, predId, skip, coretime, cluster_cals);
// }



vector<vector<vector<vector<double>>>> LB_PC_Preprocess(
    const vector<vector<double>>& s1,
    int windowSize,
    int K,
    double Q) {
    const int len1 = s1.size();

    vector<vector<vector<vector<double>>>> bboxes;
    for (int i = 0; i < len1; ++i) {
        int left = std::max(0, i - windowSize);
        int right = std::min(len1-1, i + windowSize);
        std::vector<std::vector<double>> sub_s1(s1.begin() + left, s1.begin() + right + 1);
        bboxes.push_back(findboxes(sub_s1, K, Q));
    }

    return bboxes;
}

double LB_PC(const vector<vector<double>>& Q,
                  const vector<vector<double>>& R,vector<vector<vector<vector<double>>>>& bboxes) {
    double LB_sum = 0;
    const int dim = Q[0].size(); // 获取维度

    // 遍历参考序列的每个点
    for(int idr = 0; idr < R.size(); ++idr) {
        const auto& r = R[idr];

        // 获取当前点对应的所有边界框
        const auto& cur_boxes = bboxes[idr];
        vector<double> bounds;

        // 遍历每个边界框
        for(const auto& box : cur_boxes) {
            const auto& l = box[0]; // 下界
            const auto& u = box[1]; // 上界

            double temp = 0;

            // 计算每个维度的距离贡献
            for(int d = 0; d < dim; ++d) {
                if(r[d] > u[d]) {
                    temp += pow(r[d] - u[d], 2);
                } else if(r[d] < l[d]) {
                    temp += pow(l[d] - r[d], 2);
                }
                // 在边界框内时不贡献距离
            }

            bounds.push_back(temp);
        }

        // 取最小边界贡献
        LB_sum += *min_element(bounds.begin(), bounds.end());
        // if (LB_sum>INF) {
        //     cout<<1;
        // }
    }

    return sqrt(LB_sum);
}