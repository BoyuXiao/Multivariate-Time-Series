#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
#include "funcs.h"
#include "lb_early.h"
#include <tuple>
#include <string>
#include <sstream>
#include <numeric>

double distance(const std::vector<double>& a, const std::vector<double>& b) {
    double sum = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        double diff = a[i] - b[i];
        sum += diff * diff;
    }
    return sum;
}

// Get the lower bounds between one query series X and one candidate series
vector<double> getLB_oneQ_qbox(const vector<vector<double>>& X,
                             const vector<vector<vector<double>>>& others,
                             const vector<pair<vector<double>, vector<double>>>& qbounds) {
    vector<double> lbs;
    size_t dim = X[0].size();


    //原版
    for (const auto& s2 : others) {
        double LB_sum = 0.0;
        for (size_t idy = 0; idy < s2.size(); ++idy) {
            const auto& l = qbounds[idy].first;
            const auto& u = qbounds[idy].second;
            double temp = 0.0;

            for (size_t idd = 0; idd < dim; ++idd) {
                double val = s2[idy][idd];
                if (val > u[idd]) temp += pow(val - u[idd], 2);
                else if (val < l[idd]) temp += pow(l[idd] - val, 2);
            }
            LB_sum += temp;
        }
        LB_sum =sqrt(LB_sum);
        // cout<<LB_sum<<endl;
        lbs.push_back(LB_sum);
    }

    //LB_MV版本
    //     for (const auto& s2 : others) {
    //     double LB_sum = 0.0;
    //     for (size_t idd = 0; idd < dim; ++idd) {
    //         for (size_t idy = 0; idy < s2.size(); ++idy) {
    //             const auto& l = qbounds[idy].first;
    //             const auto& u = qbounds[idy].second;
    //             double temp = 0.0;
    //             double val = s2[idy][idd];
    //             if (val > u[idd]) temp += pow(val - u[idd], 2);
    //             else if (val < l[idd]) temp += pow(l[idd] - val, 2);
    //             LB_sum += temp;
    //         }
    //     }
    //     LB_sum = sqrt(LB_sum);
    //     // cout<<LB_sum<<endl;
    //     lbs.push_back(LB_sum);
    // }


    return lbs;
}

// Compute the DTW distance between s1 and s2, and also the neighbor distances of s1
pair<double, vector<double>> DTWwnd(const std::vector<std::vector<double>>& s1, const std::vector<std::vector<double>>& s2, int windowSize) {
    const int len1 = s1.size();
    const int len2 = s2.size();
    const int w = max(windowSize, abs(len1 - len2));

    vector<vector<double>> DTW(len1, vector<double>(len2, INF));

    DTW[0][0] = distance(s1[0], s2[0]);

    vector<double> dxx;

    for (int i = 0; i < len1-1; ++i) {
        dxx.push_back(distance(s1[i+1], s1[i]));
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
    }

    // 处理最后一个元素
        const int i = len1-1;
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

    return {sqrt(DTW[len1-1][len2-1]), dxx};
}

// Compute the DTW distance between s1 and s2, and also the neighbor distances of s1
vector<double> DTWwnd_Pre(const std::vector<std::vector<double>>& s1, const std::vector<std::vector<double>>& s2) {
    const int len1 = s1.size();
    const int len2 = s2.size();

    vector<vector<double>> DTW(len1, vector<double>(len2, INF));

    DTW[0][0] = distance(s1[0], s2[0]);
    vector<double> dxx;
    for (int i = 0; i < len1-1; ++i) {
        dxx.push_back(distance(s1[i+1], s1[i]));
    }

    return  dxx;
}


// Same as tiBounds except that the true distances are calculated in every P samples of X
// And early abondoning is used.
double tiBounds_top_calP_list_comp_eb(
    const vector<vector<double>>& X,
    const vector<vector<double>>& Y,
    int W,
    int P,
    const vector<double>& dxx,
    double dist_threshold)
{
    const int Xlen = X.size();
    const int Ylen = Y.size();
    const int window_size = 2 * W + 1;

    // 初始化上下界矩阵
    vector<vector<double>> upperBounds(Xlen, vector<double>(window_size, INF));
    vector<vector<double>> lowerBounds(Xlen, vector<double>(window_size, INF));

    double lbrst = 0.0;

        for (int t = 0; t < Xlen; ++t) {
        // if (t==141) {
        //     cout << "lbrst = " << lbrst << endl;
        // }

        const int startIdx = (t > W) ? 0 : W - t;

        if (t % P == 0) {
            // 计算真实距离的case
            const int lw = max(0, t - W);
            const int tp = min(t + W + 1, Ylen);
            vector<double> dxyInit;

            // 计算初始距离
            for (int i = lw; i < tp; ++i) {
                dxyInit.push_back(distance(X[t], Y[i]));
            }

            // 填充到矩阵
            const int copy_len = tp - lw;
            fill(upperBounds[t].begin() + startIdx, upperBounds[t].begin() + startIdx + copy_len, INF);
            copy(dxyInit.begin(), dxyInit.end(), upperBounds[t].begin() + startIdx);
            copy(dxyInit.begin(), dxyInit.end(), lowerBounds[t].begin() + startIdx);

            // 更新累计下界
            auto min_val = *min_element(dxyInit.begin(), dxyInit.end());
            lbrst += min_val;
        } else {
            // 使用递推公式计算
            const int t_1 = t - 1;
            const int lr = max(0, t - W);
            const int ur = min(Ylen - 1, t + W);
            const double thisdxx = dxx[t_1];
            const int startIdx_lr = startIdx - lr + 1;

            // 计算有效窗口范围
            const int valid_len = ur - lr + 1;
            const int idx = ur - lr - 1;

            // 处理正常情况
            if (t + W <= Ylen - 1) {
                // 更新上界
                for (int i = lr,row=0; i < ur; ++i,row++) {
                    const int prev_col = startIdx_lr + i;
                    if (prev_col >= 0 && prev_col < window_size) {
                        upperBounds[t][startIdx + row] = upperBounds[t_1][prev_col] + thisdxx;
                    }
                }

                // 更新下界
                for (int i = lr,row=0; i < ur; ++i,row++) {
                    const int prev_col = startIdx_lr + i;
                    if (prev_col >= 0 && prev_col < window_size) {
                        const double prev_lower = lowerBounds[t_1][prev_col];
                        if (prev_lower > thisdxx) {
                            lowerBounds[t][startIdx + row] = prev_lower - thisdxx;
                        } else {
                            const double prev_upper = upperBounds[t_1][prev_col];
                            lowerBounds[t][startIdx + row] = (thisdxx < prev_upper) ? 0 : thisdxx - prev_upper;
                        }
                    }
                }

                // 处理最后一个点
                const double temp = distance(X[t], Y[ur]);
                upperBounds[t][startIdx + idx + 1] = temp;
                lowerBounds[t][startIdx + idx + 1] = temp;

                // 计算当前行最小值
                auto begin = lowerBounds[t].begin() + startIdx;
                auto end = begin + idx + 2;
                lbrst += *min_element(begin, end);
            } else {
                // 处理边界case
                for (int i = lr,row=0; i < ur+1; ++i,row++) {
                    const int prev_col = startIdx_lr + i;
                    if (prev_col >= 0 && prev_col < window_size) {
                        upperBounds[t][startIdx + row] = upperBounds[t_1][prev_col] + thisdxx;

                        const double prev_lower = lowerBounds[t_1][prev_col];
                        lowerBounds[t][startIdx + row] = (prev_lower > thisdxx) ?
                            prev_lower - thisdxx :
                            (thisdxx < upperBounds[t_1][prev_col]) ? 0 : thisdxx - upperBounds[t_1][prev_col];
                    }
                }

                // 计算当前行最小值
                auto begin = lowerBounds[t].begin() + startIdx;
                auto end = begin + idx + 2;
                lbrst += *min_element(begin, end);
            }
        }
            // 早期终止检查
            if (sqrt(lbrst) >= dist_threshold) {
                return sqrt(lbrst);
            }

    }





    return sqrt(lbrst);
}

double DTW_a(const vector<vector<double>>& s1,
             const vector<vector<double>>& s2,
             int windowSize,
             double bestdist) {
    const size_t n = s1.size();
    const size_t m = s2.size();
    const int w = max(windowSize, abs(static_cast<int>(n) - static_cast<int>(m)));

    // 初始化DTW矩阵
    vector<vector<double>> DTW(n, vector<double>(m, INF));
    DTW[0][0] = 0.0;

    for (size_t i = 0; i < n; ++i) {
        double current_min = INF;
        const int left = max(0, static_cast<int>(i) - w);
        const int right = min(static_cast<int>(m), static_cast<int>(i) + w + 1);

        for (int j = left; j < right; ++j) {
            const double cost = distance(s1[i], s2[j]);
            double min_prev = INF;

            if (i > 0 && j > 0) min_prev = min(min_prev, DTW[i-1][j-1]);
            if (i > 0)           min_prev = min(min_prev, DTW[i-1][j]);
            if (j > 0)           min_prev = min(min_prev, DTW[i][j-1]);

            DTW[i][j] = cost + (min_prev != INF ? min_prev : 0);
            current_min = min(current_min, DTW[i][j]);
        }

        // 早期终止检查
        if (current_min >= bestdist) {
            return current_min;
        }
    }

    return sqrt(DTW[n-1][m-1]);
}

//
// Compute the DTW distance between a query series and a reference series.
tuple< int, int> DTWDistanceWindowLB_Ordered_LBTI(
    const vector<vector<double>>& query,
    const vector<vector<vector<double>>>& references,
    int W,double threshold, int P) {

    int skips = 0, p_cals = 0;
    size_t ql = query.size();
    size_t dim = query[0].size();

    // 计算查询边界
    vector<pair<vector<double>, vector<double>>> qbounds;
    for (size_t idx = 0; idx < ql; ++idx) {
        int start = max(0, static_cast<int>(idx) - W);
        int end = min(static_cast<int>(ql), static_cast<int>(idx) + W + 1);

        vector<double> l(dim, INF);
        vector<double> u(dim, -INF);
        for (int i = start; i < end; ++i) {
            for (size_t d = 0; d < dim; ++d) {
                l[d] = min(l[d], query[i][d]);
                u[d] = max(u[d], query[i][d]);
            }
        }
        qbounds.emplace_back(l, u);
    }

    // 计算初始下界
    auto LBs = getLB_oneQ_qbox(query, references, qbounds);
    vector<size_t> LBSortedIndex(LBs.size());
    iota(LBSortedIndex.begin(), LBSortedIndex.end(), 0);
    sort(LBSortedIndex.begin(), LBSortedIndex.end(),
        [&LBs](size_t i1, size_t i2) { return LBs[i1] < LBs[i2]; });

    size_t predId = LBSortedIndex[0];
    auto [dist, dxx] = DTWwnd(query, references[predId], W);

    // 过滤处理
    for (size_t x = 1; x < LBSortedIndex.size(); ++x) {
        size_t thisrefid = LBSortedIndex[x];
        // cout<<x<<" "<<LBs[thisrefid]<<endl;
        if (LBs[thisrefid] >= threshold) {
            skips = LBs.size() - x;
            break;
        }

        if (LBs[thisrefid] >= 0.9 * threshold) {
            double p_lb = tiBounds_top_calP_list_comp_eb(query, references[thisrefid], W, P, dxx, threshold);
            // cout<<p_lb<<endl;
            if (p_lb < threshold) {
                auto new_dist = DTW_a(query, references[thisrefid], W,threshold);
                if (new_dist < threshold) {
                    p_cals++;
                }
            } else {
                skips = LBs.size() - x;
                break;
            }
        } else {
            auto new_dist = DTW_a(query, references[thisrefid], W,threshold);
            if (new_dist < threshold) {
                p_cals++;
            }
        }
    }

    return {skips, p_cals};
}

double LB_TI_early(const vector<vector<double>>& X, const vector<vector<double>>& Y, int W, int P,
                   double best_dist) {
    const int Xlen = static_cast<int>(X.size());
    const int Ylen = static_cast<int>(Y.size());
    const int window_size = 2 * W + 1;

    auto res = DTWwnd(X, Y, W);
    vector<double> dxx = res.second;

    vector<vector<double>> upperBounds(Xlen, vector<double>(window_size, INF));
    vector<vector<double>> lowerBounds(Xlen, vector<double>(window_size, INF));

    double lbrst = 0.0;
    const double thresh_sq = lb_threshold_sq(best_dist);

    for (int t = 0; t < Xlen; ++t) {
        const int startIdx = (t > W) ? 0 : W - t;

        if (t % P == 0) {
            const int lw = max(0, t - W);
            const int tp = min(t + W + 1, Ylen);
            vector<double> dxyInit;
            for (int i = lw; i < tp; ++i) {
                dxyInit.push_back(distance(X[t], Y[i]));
            }
            const int copy_len = tp - lw;
            fill(upperBounds[t].begin() + startIdx, upperBounds[t].begin() + startIdx + copy_len, INF);
            copy(dxyInit.begin(), dxyInit.end(), upperBounds[t].begin() + startIdx);
            copy(dxyInit.begin(), dxyInit.end(), lowerBounds[t].begin() + startIdx);
            auto min_val = *min_element(dxyInit.begin(), dxyInit.end());
            lbrst += min_val;
        } else {
            const int lr = max(0, t - W);
            const int ur = min(Ylen - 1, t + W);
            const double thisdxx = dxx[t - 1];
            const int startIdx_lr = startIdx - lr + 1;
            const int t_1 = t - 1;
            const int idx = ur - lr - 1;

            if (t + W <= Ylen - 1) {
                for (int i = lr, row = 0; i < ur; ++i, row++) {
                    const int prev_col = startIdx_lr + i;
                    upperBounds[t][startIdx + row] = upperBounds[t_1][prev_col] + thisdxx;
                }
                for (int i = lr, row = 0; i < ur; ++i, row++) {
                    const int prev_col = startIdx_lr + i;
                    const double prev_lower = lowerBounds[t_1][prev_col];
                    const double prev_upper = upperBounds[t_1][prev_col];
                    if (prev_lower > thisdxx) {
                        lowerBounds[t][startIdx + row] = prev_lower - thisdxx;
                    } else {
                        lowerBounds[t][startIdx + row] = (thisdxx < prev_upper) ? 0 : thisdxx - prev_upper;
                    }
                }
                const double temp = distance(X[t], Y[ur]);
                upperBounds[t][startIdx + idx + 1] = temp;
                lowerBounds[t][startIdx + idx + 1] = temp;
                auto begin = lowerBounds[t].begin() + startIdx;
                auto end = begin + idx + 2;
                lbrst += *min_element(begin, end);
            } else {
                for (int i = lr, row = 0; i < ur + 1; ++i, row++) {
                    const int prev_col = startIdx_lr + i;
                    upperBounds[t][startIdx + row] = upperBounds[t_1][prev_col] + thisdxx;
                    const double prev_lower = lowerBounds[t_1][prev_col];
                    lowerBounds[t][startIdx + row] = (prev_lower > thisdxx) ?
                        prev_lower - thisdxx :
                        (thisdxx < upperBounds[t_1][prev_col]) ? 0 : thisdxx - upperBounds[t_1][prev_col];
                }
                auto begin = lowerBounds[t].begin() + startIdx;
                auto end = begin + idx + 2;
                lbrst += *min_element(begin, end);
            }
        }

        if (lb_can_early_stop(best_dist) && lbrst >= thresh_sq) {
            return sqrt(lbrst);
        }
    }
    return sqrt(lbrst);
}

double LB_TI(const vector<vector<double>>& X,const vector<vector<double>>& Y,int W,int P) {
    const int Xlen = X.size();
    const int Ylen = Y.size();
    const int window_size = 2 * W + 1;

    // auto dxx = DTWwnd_Pre(X, Y);
    auto res = DTWwnd(X, Y,W);
    vector<double>dxx = res.second;

    // 初始化上下界矩阵
    vector<vector<double>> upperBounds(Xlen, vector<double>(window_size, INF));
    vector<vector<double>> lowerBounds(Xlen, vector<double>(window_size, INF));

    double lbrst = 0.0;

    for (int t = 0; t < Xlen; ++t) {
        const int startIdx = (t > W) ? 0 : W - t;

        if (t % P == 0) {
            // 计算真实距离的case
            const int lw = max(0, t - W);
            const int tp = min(t + W + 1, Ylen);
            vector<double> dxyInit;

            // 计算初始距离
            for (int i = lw; i < tp; ++i) {
                dxyInit.push_back(distance(X[t], Y[i]));
            }

            // 填充到矩阵
            const int copy_len = tp - lw;
            fill(upperBounds[t].begin() + startIdx, upperBounds[t].begin() + startIdx + copy_len, INF);
            copy(dxyInit.begin(), dxyInit.end(), upperBounds[t].begin() + startIdx);
            copy(dxyInit.begin(), dxyInit.end(), lowerBounds[t].begin() + startIdx);

            // 更新累计下界
            auto min_val = *min_element(dxyInit.begin(), dxyInit.end());
            lbrst += min_val;
        }
        else
            {
            // 使用递推公式计算
           
            const int lr = max(0, t - W);
            const int ur = min(Ylen - 1, t + W);
            const double thisdxx = dxx[t-1];
            const int startIdx_lr = startIdx - lr + 1;
            const int t_1 = t-1;
            const int idx = ur - lr - 1;

            // 处理正常情况
            if (t + W <= Ylen - 1) {
                // 更新上界
                for (int i = lr,row=0; i < ur; ++i,row++) {
                    const int prev_col = startIdx_lr + i;
                    // if (prev_col >= 0 && prev_col < window_size) {
                        upperBounds[t][startIdx + row] = upperBounds[t_1][prev_col] + thisdxx;
                    // }
                }

                // 更新下界
                for (int i = lr,row=0; i < ur; ++i,row++) {
                    const int prev_col = startIdx_lr + i;
                    const double prev_lower = lowerBounds[t_1][prev_col];
                    const double prev_upper = upperBounds[t_1][prev_col];
                    if (prev_lower > thisdxx) {
                        lowerBounds[t][startIdx + row] = prev_lower - thisdxx;
                    } else {
                        lowerBounds[t][startIdx + row] = (thisdxx < prev_upper) ? 0 : thisdxx - prev_upper;
                    }
                }

                // 处理最后一个点
                const double temp = distance(X[t], Y[ur]);
                upperBounds[t][startIdx + idx + 1] = temp;
                lowerBounds[t][startIdx + idx + 1] = temp;

                // 计算当前行最小值
                auto begin = lowerBounds[t].begin() + startIdx;
                auto end = begin + idx + 2;
                lbrst += *min_element(begin, end);
            } else {
                // 处理边界case
                for (int i = lr,row=0; i < ur+1; ++i,row++) {
                    const int prev_col = startIdx_lr + i;
                    upperBounds[t][startIdx + row] = upperBounds[t_1][prev_col] + thisdxx;
                    const double prev_lower = lowerBounds[t_1][prev_col];
                    lowerBounds[t][startIdx + row] = (prev_lower > thisdxx) ?
                        prev_lower - thisdxx :
                        (thisdxx < upperBounds[t_1][prev_col]) ? 0 : thisdxx - upperBounds[t_1][prev_col];
                }

                // 计算当前行最小值
                auto begin = lowerBounds[t].begin() + startIdx;
                auto end = begin + idx + 2;
                lbrst += *min_element(begin, end);
            }
        }
    }

    return sqrt(lbrst);
}
