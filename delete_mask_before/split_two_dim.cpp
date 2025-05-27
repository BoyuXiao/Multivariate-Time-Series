#include <vector>
#include "funcs.h"
void split_two_dim(const std::vector<std::vector<double>>& Q,
            vector<double>& Q1,
            vector<double>& Q2) {
    int m = Q.size();  // 行数
    Q1.resize(m);
    Q2.resize(m);
    for (int i = 0; i < m; ++i) {
        Q1[i] = Q[i][0];
        Q2[i] = Q[i][1];
    }
}