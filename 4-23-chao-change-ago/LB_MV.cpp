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
    // cout<<sqrt(lb)<<endl;
    return sqrt(lb);
}

double LB_MV_graph(vector<double>& Q, vector<vector<double>>& A, int K) {
    double lb = 0.0;
    int m = A.size();
    for (int k = 0; k < K; ++k) {
        vector<double> A_col;
        for (int i = 0; i < m; ++i) {
            A_col.push_back(A[i][k]);
        }
        auto max_item = max_element(A_col.begin(), A_col.end());
        auto min_item = min_element(A_col.begin(), A_col.end());
        double upperBound = *max_item;
        double lowerBound = *min_item;

        if (Q[k] < lowerBound) {
            lb+= dist(Q[k], lowerBound);
        }else if (Q[k] > upperBound) {
            lb+= dist(Q[k], upperBound);
        }
    }
    return lb;
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
    vector<double> lowerBound(K);
    vector<double> upperBound(K);
    for (int k = 0; k < m; ++k) {
        lower_upper_lemire(Q[k],K,r,lowerBound,upperBound);
        for (int i=0;i<K;++i) {
            if (A[k][i]<lowerBound[i]) {
                lb+=dist(A[k][i],lowerBound[i]);
            }else if (A[k][i]>upperBound[i]) {
                lb+=dist(A[k][i],upperBound[i]);
            }
        }
    }
    return sqrt(lb);
}


double LB_MV_graph_mask(vector<double>& Q, vector<vector<double>>& A, int K,int start,int end,vector<char>&mask) {
    double lb = 0.0;
    for (int k = 0; k < K; ++k) {
        vector<double> A_col;
        for (int i = start; i <= end; ++i) {
            A_col.push_back(A[i][k]);
        }
        auto max_item = max_element(A_col.begin(), A_col.end());
        auto min_item = min_element(A_col.begin(), A_col.end());
        double upperBound = *max_item;
        double lowerBound = *min_item;
        int max_index = start + distance(A_col.begin(), max_item);
        int min_index = start + distance(A_col.begin(), min_item);

        if (Q[k] < lowerBound) {
            lb+= dist(Q[k], lowerBound);
            mask[min_index]=1;
        }else if (Q[k] > upperBound) {
            lb+= dist(Q[k], upperBound);
            mask[max_index]=1;
        }
    }
    return lb;
}

//不是K次累加
// double LB_MV_double(vector<vector<double>>& A,vector<vector<double>>& Q,vector<double>&lb_record,vector<char>&mask,int r) {
//
//     int m = A.size();
//     int K = A[0].size();
//     double lb = 0.0;
//
//     for (int i=0;i<m;i++) {
//         double delta = 0;
//         int start = MAX(0,i-r);
//         int end = MIN(m-1,i+r);
//         delta = LB_MV_graph_mask(A[i],Q,K,start,end,mask);
//         lb_record[i] = delta;
//         lb += lb_record[i];
//     }
//
//       // // //计算右侧空缺的点
//       vector<double>max_val;
//       upper_lemire(lb_record,m,r,max_val);
//       for (int i=0;i<m;i++) {
//           if (mask[i])continue;
//           double delta = 0;
//           int start = MAX(0,i-r);
//           int end = MIN(m-1,i+r);
//           vector<vector<double>> A_part(A.begin() + start, A.begin() + end + 1);
//           delta = LB_MV_graph(Q[i],A_part,K);
//           lb += MAX(0.0,delta-max_val[i]);
//       }
//
//       return sqrt(lb);
// }

//K次累加 转置版本
double LB_MV_double_T(vector<vector<double>>& A,vector<vector<double>>& Q, int r) {
    int m = A[0].size();
    int K = A.size();
    double lb = 0.0;
    vector<pair<int,double>> Q_L(m),Q_U(m),A_L(m),A_U(m);
    vector<double>max_val(m);
    vector<double> lb_record(m,0);
    vector<char> mask(m,0);
    for (int k=0;k<K;k++) {
        lb_record.assign(m, 0);
        mask.assign(m, 0);
        lower_upper_lemire_index(Q[k],m,r,Q_L,Q_U);
        lower_upper_lemire_index(A[k],m,r,A_L,A_U);
        for (int i=0;i<m;i++) {
            if (A[k][i]<Q_L[i].second) {
                lb_record[i] = dist(A[k][i],Q_L[i].second);
                mask[Q_L[i].first] = 1;
            }else if (A[k][i]>Q_U[i].second) {
                lb_record[i] = dist(A[k][i],Q_U[i].second);
                mask[Q_U[i].first] = 1;
            }
            lb+=lb_record[i];
        }
        //二次计算
        upper_lemire(lb_record,m,r,max_val);
        for (int i=0;i<m;i++) {
            if (mask[i]) continue;
            double delta = 0;
            if (Q[k][i]<A_L[i].second) {
                delta = dist(Q[k][i],A_L[i].second);
            }else if (Q[k][i]>A_U[i].second) {
                delta = dist(Q[k][i],A_U[i].second);
            }
            lb+=MAX(0.0,delta-max_val[i]);
        }
    }

    return sqrt(lb);
}