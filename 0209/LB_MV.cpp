#include <iostream>
#include <ostream>

#include "funcs.h"
using namespace std;
double LB_MV(vector<vector<double>>& Q, vector<vector<double>>& A, int K, int m, int r) {
    double lb = 0.0;
    for (int k = 0; k < K; ++k) {
        vector<double> Q_col;
        vector<double> A_col;
        for (int i = 0; i < Q.size(); ++i) {
            Q_col.push_back(Q[i][k]);
        }
        for (int i = 0; i < A.size(); ++i) {
            A_col.push_back(A[i][k]);
        }
        lb += LB_Keogh_2(Q_col, A_col, m, r);
    }
    return sqrt(lb);
}

double LB_MV_T(vector<vector<double>>& Q, vector<vector<double>>& A, int K, int m, int r) {
    double lb = 0.0;
    for (int k = 0; k < m; ++k) {
        lb += LB_Keogh_2(Q[k], A[k], K, r);
    }
    return sqrt(lb);
}
