#include <iostream>
#include <vector>
#include <algorithm>
#include "funcs.h"
#include <utility>
#include <random>
#include <unordered_set>
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


void splitQ_sort_strength(const std::vector<std::vector<double>>& Q,
            std::vector<std::vector<double>>& Q1,
            std::vector<std::vector<double>>& Q2,
            vector<int>& split_indices) {
    split_indices.clear();
    int m = Q.size();
    int K = Q[0].size();
    int split_size = K / 2;
    vector<double> dimension_sums(K, 0.0);

    for (int i=0;i<m;++i) {
        for (int k=0;k<K;++k) {
            dimension_sums[k] += Q[i][k];
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


void splitA_sort_single(const vector<double>& A,
                 const std::vector<int>& split_indices,
                 vector<double>& A1,
                 vector<double>& A2) {
    A1.clear();
    A2.clear();
    int K = A.size();
    for (int k = 0; k < K; ++k) {
        if (std::find(split_indices.begin(), split_indices.end(), k) != split_indices.end()) {
            A1.push_back(A[k]);
        } else {
            A2.push_back(A[k]);
        }
    }
}
// --------------------------------------------------------------------
void split_half(const std::vector<std::vector<double>>& Q,
            std::vector<std::vector<double>>& Q1,
            std::vector<std::vector<double>>& Q2) {
    int m = Q.size();
    int K = Q[0].size();
    int split_size = K / 2;
    Q1.resize(m);
    Q2.resize(m);
    for (int i = 0; i < m; ++i) {
        Q1[i].assign(Q[i].begin(), Q[i].begin() + split_size);
        Q2[i].assign(Q[i].begin() + split_size, Q[i].end());
    }
}

void split_first_and_second(const std::vector<std::vector<double>>& Q,
            std::vector<std::vector<double>>& Q1) {
    int m = Q.size();
    int split_size = 2;
    Q1.resize(m);
    for (int i = 0; i < m; ++i) {
        Q1[i].assign(Q[i].begin(), Q[i].begin() + split_size);
    }
}

void split_third(const std::vector<std::vector<double>>& Q,
                 std::vector<std::vector<double>>& Q1,
                 std::vector<std::vector<double>>& Q2,
                 std::vector<std::vector<double>>& Q3) {
    int m = Q.size();
    int K = Q[0].size();
    int split_size = K / 3;
    for (int i = 0; i < m; ++i) {
        Q1[i].assign(Q[i].begin(), Q[i].begin() + split_size);
        Q2[i].assign(Q[i].begin() + split_size, Q[i].begin() + 2 * split_size);
        Q3[i].assign(Q[i].begin() + 2 * split_size, Q[i].end());
    }
}
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
        // cout<<means[k]<<endl;
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

    //打印cov_matrix
    // for (int k1 = 0; k1 < K; ++k1) {
    //     for (int k2 = 0; k2 < K; ++k2) {
    //         cout << cov_matrix[k1][k2] << ", ";
    //     }
    //     cout << endl;
    // }
    // cout<<endl<<endl<<endl;
    return cov_matrix;
}


vector<int> kmedoids_clustering(const vector<vector<double>>& distance_matrix, int num_clusters, int max_iters = 100) {
    int K = distance_matrix.size();
    vector<int> medoids(num_clusters);
    vector<int> labels(K, 0);
    mt19937 gen(42);
    uniform_int_distribution<> dis(0, K-1);

    // 初始化medoids，确保不重复
    unordered_set<int> selected_medoids;
    medoids[0] = dis(gen);
    selected_medoids.insert(medoids[0]);
    for (int i = 1; i < num_clusters; ++i) {
        int candidate;
        do {
            candidate = dis(gen);
        } while (selected_medoids.count(candidate));
        medoids[i] = candidate;
        selected_medoids.insert(candidate);
    }

    bool changed;
    for (int iter = 0; iter < max_iters; ++iter) {
        changed = false;

        // 分配每个点到最近的medoid
        for (int i = 0; i < K; ++i) {
            double min_dist = MAX_VAL;
            int best_cluster = -1;
            for (int m = 0; m < num_clusters; ++m) {
                double dist = distance_matrix[i][medoids[m]];
                if (dist < min_dist) {
                    min_dist = dist;
                    best_cluster = m;
                }
            }
            if (labels[i] != best_cluster) {
                labels[i] = best_cluster;
                changed = true;
            }
        }

        if (!changed) break;

        // 更新medoids
        for (int m = 0; m < num_clusters; ++m) {
            vector<int> cluster_points;
            for (int i = 0; i < K; ++i) {
                if (labels[i] == m) {
                    cluster_points.push_back(i);
                }
            }

            double min_total_dist = MAX_VAL;
            int new_medoid = medoids[m];
            for (int candidate : cluster_points) {
                double total_dist = 0.0;
                for (int point : cluster_points) {
                    total_dist += distance_matrix[point][candidate];
                }
                if (total_dist < min_total_dist) {
                    min_total_dist = total_dist;
                    new_medoid = candidate;
                }
            }

            if (new_medoid != medoids[m]) {
                medoids[m] = new_medoid;
                changed = true;
            }
        }

        if (!changed) break;
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

    // 将协方差矩阵取负数
    std::vector<std::vector<double>> negative_cov_matrix(K, std::vector<double>(K));
    for (int i = 0; i < K; ++i) {
        for (int j = 0; j < K; ++j) {
            negative_cov_matrix[i][j] = -cov_matrix[i][j];
        }
    }

    // std::vector<int> labels = kmedoids_clustering(cov_matrix, 2);
    std::vector<int> labels = kmedoids_clustering(negative_cov_matrix, 2);

    split_indices.clear();
    for (int k = 0; k < K; ++k) {
        if (labels[k] == 0) {
            split_indices.push_back(k);
        }
    }

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

