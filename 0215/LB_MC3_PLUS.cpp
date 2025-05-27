#include <iostream>
#include <vector>
#include <algorithm>
#include "funcs.h"
#include <utility>
#include <random>
// --------------------------------------------------------------------
void splitQ_sort(const std::vector<std::vector<double>>& Q,
            std::vector<std::vector<double>>& Q1,
            std::vector<std::vector<double>>& Q2,
            vector<int>& split_indices,
            int r) {
    split_indices.clear();
    int m = Q.size();
    int K = Q[0].size();
    int split_size = K / 2;
    vector<int>start(m);
    vector<int>end(m);
    for (int i = 0; i < m; ++i) {
        start[i] = MAX(0,i-r);
        end[i] = MIN(m-1,i+r);
    }

    vector<double> dimension_sums(K, 0.0);

    for (int i=0;i<m;++i) {
        for (int j=start[i];j<=end[i];++j) {
            for (int k=0;k<K;++k) {
                dimension_sums[k] += Q[j][k];
            }
        }
    }

    vector<pair<double, int>> sorted_dimensions;
    for (int k = 0; k < K; ++k) {
        sorted_dimensions.emplace_back(dimension_sums[k], k);
    }
    sort(sorted_dimensions.rbegin(), sorted_dimensions.rend());

    for (int i = 0; i < split_size; ++i) {
        split_indices.push_back(sorted_dimensions[i].second);
    }

    for (int i = 0; i < m; ++i) {
        vector<double> q1_dim, q2_dim;
        for (int k = 0; k < K; ++k) {
            if (find(split_indices.begin(), split_indices.end(), k) != split_indices.end()) {
                q1_dim.push_back(Q[i][k]);
            } else {
                q2_dim.push_back(Q[i][k]);
            }
        }
        Q1[i] = q1_dim;
        Q2[i] = q2_dim;
    }
}

void splitA_sort(const std::vector<std::vector<double>>& A,
                 const std::vector<int>& split_indices,
                 std::vector<std::vector<double>>& A1,
                 std::vector<std::vector<double>>& A2) {
    int m = A.size();
    int K = A[0].size();

    for (int i = 0; i < m; ++i) {
        vector<double> a1_dim, a2_dim;
        for (int k = 0; k < K; ++k) {
            if (find(split_indices.begin(), split_indices.end(), k) != split_indices.end()) {
                a1_dim.push_back(A[i][k]);
            } else {
                a2_dim.push_back(A[i][k]);
            }
        }
        A1[i] = a1_dim;
        A2[i] = a2_dim;
    }
}
// --------------------------------------------------------------------
void split_half(const std::vector<std::vector<double>>& Q,
            std::vector<std::vector<double>>& Q1,
            std::vector<std::vector<double>>& Q2) {
    int m = Q.size();
    int K = Q[0].size();
    int split_size = K / 2;
    for (int i = 0; i < m; ++i) {
        Q1[i].assign(Q[i].begin(), Q[i].begin() + split_size);
        Q2[i].assign(Q[i].begin() + split_size, Q[i].end());
    }
}

// --------------------------------------------------------------------
// vector<vector<double>> compute_covariance_matrix(const vector<vector<double>>& Q) {
//     int m = Q.size();
//     int K = Q[0].size();
//     vector<vector<double>> cov_matrix(K, vector<double>(K, 0.0));
//
//     vector<double> means(K, 0.0);
//     for (int k = 0; k < K; ++k) {
//         double sum = 0.0;
//         for (int i = 0; i < m; ++i) {
//             sum += Q[i][k];
//         }
//         means[k] = sum / m;
//     }
//
//     for (int k1 = 0; k1 < K; ++k1) {
//         for (int k2 = 0; k2 < K; ++k2) {
//             double cov = 0.0;
//             for (int i = 0; i < m; ++i) {
//                 cov += (Q[i][k1] - means[k1]) * (Q[i][k2] - means[k2]);
//             }
//             cov_matrix[k1][k2] = cov / m;
//         }
//     }
//
//     return cov_matrix;
// }
//
//
// std::vector<int> kmeans_clustering(const std::vector<std::vector<double>>& data, int num_clusters, int max_iters = 100) {
//     int n = data.size();
//     int dim = data[0].size();
//
//     std::mt19937 gen(42); // 使用固定种子 42 初始化随机数生成器
//     std::uniform_int_distribution<> dis(0, n - 1); // 均匀分布，范围 [0, n-1]
//     std::vector<std::vector<double>> centers(num_clusters);
//     for (int i = 0; i < num_clusters; ++i) {
//         centers[i] = data[dis(gen)];
//     }
//
//     // 初始化聚类标签
//     std::vector<int> labels(n, 0);
//
//     for (int iter = 0; iter < max_iters; ++iter) {
//         // 分配每个数据点到最近的聚类中心
//         for (int i = 0; i < n; ++i) {
//             double min_dist = MAX_VAL;
//             for (int j = 0; j < num_clusters; ++j) {
//                 double dist = 0.0;
//                 for (int d = 0; d < dim; ++d) {
//                     dist += (data[i][d] - centers[j][d])*(data[i][d] - centers[j][d]);
//                 }
//                 if (dist < min_dist) {
//                     min_dist = dist;
//                     labels[i] = j;
//                 }
//             }
//         }
//
//         // 更新聚类中心
//         std::vector<std::vector<double>> new_centers(num_clusters, std::vector<double>(dim, 0.0));
//         std::vector<int> counts(num_clusters, 0);
//         for (int i = 0; i < n; ++i) {
//             int cluster = labels[i];
//             for (int d = 0; d < dim; ++d) {
//                 new_centers[cluster][d] += data[i][d];
//             }
//             counts[cluster]++;
//         }
//         for (int j = 0; j < num_clusters; ++j) {
//             if (counts[j] > 0) {
//                 for (int d = 0; d < dim; ++d) {
//                     centers[j][d] = new_centers[j][d] / counts[j];
//                 }
//             }
//         }
//     }
//
//     return labels;
// }

//修正版本--------------------------------------------------------------------------------------------------------------------------------------
vector<vector<double>> compute_covariance_matrix(const vector<vector<double>>& Q) {
    int m = Q.size();
    int K = Q[0].size();
    vector<vector<double>> cov_matrix(K, vector<double>(K, 0.0));

    vector<double> means(K, 0.0);
    for (int k = 0; k < K; ++k) {
        double sum = 0.0;
        for (int i = 0; i < m; ++i) {
            sum += Q[i][k];
        }
        means[k] = sum / m;
    }

    for (int k1 = 0; k1 < K; ++k1) {
        for (int k2 = 0; k2 < K; ++k2) {
            double cov = 0.0;
            for (int i = 0; i < m; ++i) {
                cov += (Q[i][k1] - means[k1]) * (Q[i][k2] - means[k2]);
            }
            cov_matrix[k1][k2] = cov / m;
        }
    }

    // //打印cov_matrix
    // for (int k1 = 0; k1 < K; ++k1) {
    //     for (int k2 = 0; k2 < K; ++k2) {
    //         cout << cov_matrix[k1][k2] << ", ";
    //     }
    //     cout << endl;
    // }
    // cout<<endl<<endl<<endl;
    return cov_matrix;
}


std::vector<int> kmedoids_clustering(const std::vector<std::vector<double>>& data, int num_clusters, int max_iters = 100) {
    int n = data.size();
    std::mt19937 gen(42);
    std::uniform_int_distribution<> dis(0, n - 1);
    std::vector<int> medoids;
    for (int i = 0; i < num_clusters; ++i) {
        int index;
        do {
            index = dis(gen);
        } while (std::find(medoids.begin(), medoids.end(), index) != medoids.end());
        medoids.push_back(index);
    }

    std::vector<int> labels(n, 0);

    for (int iter = 0; iter < max_iters; ++iter) {
        // 分配每个数据点到最近的 medoid
        for (int i = 0; i < n; ++i) {
            double min_dist = MAX_VAL;
            for (int j = 0; j < num_clusters; ++j) {
                double dist = data[i][medoids[j]];
                if (dist < min_dist) {
                    min_dist = dist;
                    labels[i] = j;
                }
            }
        }

        // 更新 medoids
        bool medoids_changed = false;
        for (int j = 0; j < num_clusters; ++j) {
            std::vector<int> cluster_points;
            for (int i = 0; i < n; ++i) {
                if (labels[i] == j) {
                    cluster_points.push_back(i);
                }
            }

            double min_cost = std::numeric_limits<double>::max();
            int new_medoid = medoids[j];
            for (int k : cluster_points) {
                double cost = 0.0;
                for (int l : cluster_points) {
                    cost += data[k][l];
                }
                if (cost < min_cost) {
                    min_cost = cost;
                    new_medoid = k;
                }
            }

            if (new_medoid != medoids[j]) {
                medoids[j] = new_medoid;
                medoids_changed = true;
            }
        }

        if (!medoids_changed) {
            break;
        }
    }

    return labels;
}

void splitQ_correlation(const std::vector<std::vector<double>>& Q,
                        std::vector<std::vector<double>>& Q1,
                        std::vector<std::vector<double>>& Q2,
                        std::vector<int>& split_indices,
                        int r) {
    int m = Q.size();    // 行数
    int K = Q[0].size(); // 列数（变量数）

    // 计算协方差矩阵
    std::vector<std::vector<double>> cov_matrix = compute_covariance_matrix(Q);

    // 将协方差矩阵作为特征进行二聚类
    std::vector<int> labels = kmedoids_clustering(cov_matrix, 2);

    // 根据聚类结果划分变量
    split_indices.clear();
    for (int k = 0; k < K; ++k) {
        if (labels[k] == 0) {
            split_indices.push_back(k); // 属于第一个聚类
        }
    }

    // 根据 split_indices 拆分 Q
    Q1.resize(m);
    Q2.resize(m);
    for (int i = 0; i < m; ++i) {
        std::vector<double> q1_dim, q2_dim;
        for (int k = 0; k < K; ++k) {
            if (std::find(split_indices.begin(), split_indices.end(), k) != split_indices.end()) {
                q1_dim.push_back(Q[i][k]); // 属于第一个聚类
            } else {
                q2_dim.push_back(Q[i][k]); // 属于第二个聚类
            }
        }
        Q1[i] = q1_dim;
        Q2[i] = q2_dim;
    }
}


// ---------------------------------------------------------

