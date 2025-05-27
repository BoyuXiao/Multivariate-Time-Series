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
        lb += LB_Keogh_2(A_col, Q_col, m, r);
    }
    cout<<sqrt(lb)<<endl;
    return sqrt(lb);
}

double LB_MV_test(vector<vector<double>>& Q, vector<vector<double>>& A, int K, int m, int r) {
    double lb = 0.0;
    for (int i=0;i<m;i++) {
        lb+=LB_Keogh_2(A[i], Q[i], K, r);
    }
    return sqrt(lb);
}

double LB_MV_T(vector<vector<double>>& Q, vector<vector<double>>& A, int K, int m, int r) {
    double lb = 0.0;
    for (int k = 0; k < m; ++k) {
        lb += LB_Keogh_2(A[k], Q[k], K, r);
    }
    // cout << lb << endl;
    return sqrt(lb);
}
