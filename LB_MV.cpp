#include <iostream>
#include <ostream>
#include <fstream>
#include "funcs.h"
using namespace std;


double LB_KIM_NEW(const std::vector<double>& t, const std::vector<double>& q, int m) {
    double d = 0;
    d += dist(q[0], t[0]);
    d += dist(q[m-1], t[m-1]);
    double left_min = MAX_VAL;
    double right_min = MAX_VAL;
    double d10 = dist(t[1],q[0]);
    double d20 = dist(t[2],q[0]);
    double d01 = dist(t[0],q[1]);
    double d11 = dist(t[1],q[1]);
    double d21 = dist(t[2],q[1]);
    double d02 = dist(t[0],q[2]);
    double d12 = dist(t[1],q[2]);
    double d22 = dist(t[2],q[2]);
    left_min = min({d10+d20,d10+d21,d11+d21,d11+d22,d11+d12,d01+d02,d01+d12});
    d+=left_min;
    d10 = dist(t[m-2],q[m-1]);
    d20 = dist(t[m-3],q[m-1]);
    d01 = dist(t[m-1],q[m-2]);
    d11 = dist(t[m-2],q[m-2]);
    d21 = dist(t[m-3],q[m-2]);
    d02 = dist(t[m-1],q[m-3]);
    d12 = dist(t[m-2],q[m-3]);
    d22 = dist(t[m-3],q[m-3]);
    right_min = min({d10+d20,d10+d21,d11+d21,d11+d22,d11+d12,d01+d02,d01+d12});
    d+=right_min;
    return d;
}




double lbPetitjean_new(std::vector<double>& q,const std::vector<double>& uq, const std::vector<double>& lq, std::vector<double>& x,
                     const std::vector<double>& ux, const std::vector<double>& lx, int w) {
    double lbk2 = 0;
    int m = q.size();
    vector<double> p(m);

    lbk2 += LB_KIM_NEW(q, x,  m);
    for (int i = 3; i < m-3; i++) {
        double qi = q[i];
        if (qi > ux[i]) {
            lbk2 += dist(qi, ux[i]);
        } else if (qi < lx[i]) {
            lbk2 += dist(qi, lx[i]);
        }
    }

    for (int i=0;i<m;i++) {
        if (i<3||i>=m-3) {
            p[i]=q[i];
        }
        else if (q[i]>ux[i]) {
            p[i]=ux[i];
        }else if (q[i]<lx[i]) {
            p[i]=lx[i];
        }else p[i]=q[i];
    }
    vector<double>lp(m),up(m);
    lower_upper_lemire(p,m,w,lp,up);

    for (int ii=3;ii<m-3;ii++) {
        if(x[ii] > up[ii] && up[ii] > uq[ii]) {
            lbk2 += (dist(x[ii], uq[ii]) - dist(up[ii], uq[ii]));
        } else if(x[ii] < lp[ii] && lp[ii] < lq[ii]) {
            lbk2 += (dist(x[ii], lq[ii]) - dist(lp[ii], lq[ii]));
        } else if(x[ii] > up[ii]) {
            lbk2 += dist(x[ii], up[ii]);
        } else if(x[ii] < lp[ii]) {
            lbk2 += dist(x[ii], lp[ii]);
        }
    }


    return lbk2;
}

double lbPetitjean_new_early(std::vector<double>& q,const std::vector<double>& uq, const std::vector<double>& lq, std::vector<double>& x,
                     const std::vector<double>& ux, const std::vector<double>& lx, int w,double threshold,int &flag,double &lbk2) {
    int m = q.size();
    vector<double> p(m);

    lbk2 += LB_KIM_NEW(q, x,  m);
    for (int i = 3; i < m-3; i++) {
        double qi = q[i];
        if (qi > ux[i]) {
            lbk2 += dist(qi, ux[i]);
        } else if (qi < lx[i]) {
            lbk2 += dist(qi, lx[i]);
        }
        if (lbk2 > threshold) {flag = 1;return lbk2;}
    }
    if (lbk2 > threshold) {flag = 1;return lbk2;}
    for (int i=0;i<m;i++) {
        if (i<3||i>=m-3) {
            p[i]=q[i];
        }
        else if (q[i]>ux[i]) {
            p[i]=ux[i];
        }else if (q[i]<lx[i]) {
            p[i]=lx[i];
        }else p[i]=q[i];
    }
    vector<double>lp(m),up(m);
    lower_upper_lemire(p,m,w,lp,up);

    for (int ii=3;ii<m-3;ii++) {
        if(x[ii] > up[ii] && up[ii] > uq[ii]) {
            lbk2 += (dist(x[ii], uq[ii]) - dist(up[ii], uq[ii]));
        } else if(x[ii] < lp[ii] && lp[ii] < lq[ii]) {
            lbk2 += (dist(x[ii], lq[ii]) - dist(lp[ii], lq[ii]));
        } else if(x[ii] > up[ii]) {
            lbk2 += dist(x[ii], up[ii]);
        } else if(x[ii] < lp[ii]) {
            lbk2 += dist(x[ii], lp[ii]);
        }
        if (lbk2 > threshold) {flag = 1;return lbk2;}
    }

    return lbk2;
}

double lbPetitjean_new_early(std::vector<double>& q,const std::vector<double>& uq, const std::vector<double>& lq, std::vector<double>& x,
                     const std::vector<double>& ux, const std::vector<double>& lx, int w,double threshold,int &flag,double &lbk2,vector<double>&cb) {
    int m = q.size();
    vector<double> p(m);

    for (int i = 0; i < m; i++) {
        double qi = q[i];
        double temp=0;
        if (qi > ux[i]) {
            temp = dist(qi, ux[i]);
        } else if (qi < lx[i]) {
            temp = dist(qi, lx[i]);
        }
        lbk2+=temp;
        cb[i]+=temp;
        if (lbk2 > threshold) {flag = 1;return lbk2;}
    }
    for (int i=0;i<m;i++) {
        if (i<3||i>=m-3) {
            p[i]=q[i];
        }
        else if (q[i]>ux[i]) {
            p[i]=ux[i];
        }else if (q[i]<lx[i]) {
            p[i]=lx[i];
        }else p[i]=q[i];
    }
    vector<double>lp(m),up(m);
    lower_upper_lemire(p,m,w,lp,up);

    for (int ii=0;ii<m;ii++) {
        double temp = 0;
        if(x[ii] > up[ii] && up[ii] > uq[ii]) {
            temp = (dist(x[ii], uq[ii]) - dist(up[ii], uq[ii]));
        } else if(x[ii] < lp[ii] && lp[ii] < lq[ii]) {
            temp = (dist(x[ii], lq[ii]) - dist(lp[ii], lq[ii]));
        } else if(x[ii] > up[ii]) {
            temp = dist(x[ii], up[ii]);
        } else if(x[ii] < lp[ii]) {
            temp = dist(x[ii], lp[ii]);
        }
        lbk2+=temp;
        cb[ii]+=temp;
        if (lbk2 > threshold) {flag = 1;return lbk2;}
    }

    return lbk2;
}

double lbWebb(const std::vector<double>& q, const std::vector<double>& uq, const std::vector<double>& lq,
              const std::vector<double>& luq, const std::vector<double>& ulq, const std::vector<double>& t,
              const std::vector<double>& ut, const std::vector<double>& lt, const std::vector<double>& lut,
              const std::vector<double>& ult, int window) {
    double lb = 0;
    int m = q.size();
    lb += LB_KIM_NEW(t,q,m);
    int freeCountAbove = window;
    int freeCountBelow = window;

    // lbKeogh
    for (int i = 3; i < m-3; ++i) {
        if (q[i] > ut[i]) {
            lb += (q[i]-ut[i])*(q[i]-ut[i]);
            if (ut[i] >= ulq[i]) {
                freeCountBelow++;
            } else {
                freeCountBelow = 0;
            }
            freeCountAbove = 0;
        } else if (q[i] < lt[i]) {
            lb += (q[i]-lt[i])*(q[i]-lt[i]);
            if (lt[i] <= luq[i]) {
                freeCountAbove++;
            } else {
                freeCountAbove = 0;
            }
            freeCountBelow = 0;
        } else {
            freeCountAbove++;
            freeCountBelow++;
        }

        if (i >= window + 3) {
            int j = i - window;
            if (t[j] > uq[j]) {
                if (freeCountAbove > 2 * window) {
                    lb += (t[j]-uq[j])*(t[j]-uq[j]);
                } else {
                    if (t[j] > ult[j] && ult[j] >= uq[j]) {
                        lb += (t[j]-uq[j])*(t[j]-uq[j]) - (ult[j]-uq[j])*(ult[j]-uq[j]);
                    }
                }
            } else {
                if (t[j] < lq[j]) {
                    if (freeCountBelow > 2 * window) {
                        lb += (t[j]-lq[j])*(t[j]-lq[j]);
                    } else {
                        if (t[j] < lut[j] && lut[j] <= lq[j]) {
                            lb += (t[j]-lq[j])*(t[j]-lq[j]) - (lut[j]-lq[j])*(lut[j]-lq[j]);
                        }
                    }
                }
            }
        }
    }

    for (int j = m-3 - window; j < m-3; ++j) {
        if (t[j] > uq[j]) {
            if (j >= m-3 - freeCountAbove + window) {
                lb += (t[j]-uq[j])*(t[j]-uq[j]);
            } else {
                if (t[j] > ult[j] &&  ult[j] >= uq[j]) {
                    lb += (t[j]-uq[j])*(t[j]-uq[j]) - ( ult[j]-uq[j])*( ult[j]-uq[j]);
                }
            }
        } else {
            if (t[j] < lq[j]) {
                if (j >= m-3 - freeCountBelow + window) {
                    lb += (t[j]-lq[j])*(t[j]-lq[j]);
                } else {
                    if (t[j] < lut[j] && lut[j] <= lq[j]) {
                        lb += (t[j]-lq[j])*(t[j]-lq[j]) - (lut[j]-lq[j])*(lut[j]-lq[j]);
                    }
                }
            }
        }
    }

    return lb;
}



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

double LB_MV_T(vector<vector<double>>& lowerBound,vector<vector<double>>& upperBound, vector<vector<double>>& A, int K, int m) {
    double lb = 0.0;
    for (int k = 0; k < m; ++k) {
        for (int i=0;i<K;++i) {
            if (A[k][i]<lowerBound[k][i]) {
                lb+=dist(A[k][i],lowerBound[k][i]);
            }else if (A[k][i]>upperBound[k][i]) {
                lb+=dist(A[k][i],upperBound[k][i]);
            }
        }
    }
    return sqrt(lb);
}

void LB_MV_T_hot_map(vector<vector<double>>& lowerBound,vector<vector<double>>& upperBound, vector<vector<double>>& A, int m, int K) {
    vector<double> lb(m,0);
    for (int k = 0; k < K; ++k) {
        for (int i=0;i<m;++i) {
            if (A[k][i]<lowerBound[k][i]) {
                lb[i]+=dist(A[k][i],lowerBound[k][i]);
            }else if (A[k][i]>upperBound[k][i]) {
                lb[i]+=dist(A[k][i],upperBound[k][i]);
            }
        }
    }
    std::ofstream outputFile;
    outputFile.open("..//hot_map.txt",std::ios::app);
    for (int i=0;i<m;++i) {
        outputFile<<lb[i]<<" ";
    }
    outputFile << "\n";
    outputFile.close();

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
double LB_MV_double_T(vector<vector<double>>& A,vector<vector<double>>& Q, vector<vector<double>>& Q_L,vector<vector<double>>& Q_U,int r) {
    int m = A[0].size();
    int K = A.size();
    double lb = 0.0;
    vector<double>A_L(m),A_U(m);
    // vector<pair<int,double>>A_L(m),A_U(m);
    vector<double>max_val(m);
    vector<double> lb_record(m,0);
    for (int k=0;k<K;k++) {
        lb_record.assign(m, 0);
        for (int i=0;i<m;i++) {
            if (A[k][i]<Q_L[k][i]) {
                lb_record[i] = dist(A[k][i],Q_L[k][i]);
            }else if (A[k][i]>Q_U[k][i]) {
                lb_record[i] = dist(A[k][i],Q_U[k][i]);
            }
            lb+=lb_record[i];
        }

        //二次计算
        lower_upper_lemire(A[k],m,r,A_L,A_U);
        upper_lemire(lb_record,m,r,max_val);
        for (int i=0;i<m;i++) {
            double delta = 0;
            if (Q[k][i]<A_L[i]) {
                delta = dist(Q[k][i],A_L[i]);
            }else if (Q[k][i]>A_U[i]) {
                delta = dist(Q[k][i],A_U[i]);
            }
            lb+=MAX(0.0,delta-max_val[i]);
        }
    }

    return sqrt(lb);
}


double LB_MV_double_T_V2(vector<vector<double>>& A,vector<vector<double>>& Q,vector<vector<double>>& Q_L,vector<vector<double>>& Q_U,int r) {
    int m = A[0].size();
    int K = A.size();
    double lb = 0.0;
    double lb_temp = 0.0;
    vector<double>max_val(m);
    vector<double>A_L(m),A_U(m);
    vector<double> lb_record(m,0);
    vector<double> delta(m,0);
    for (int k=0;k<K;k++) {
        for (int i=0;i<m;i++) {
            lb_temp = 0;
            if (A[k][i]<Q_L[k][i]) {
                lb_temp = dist(A[k][i],Q_L[k][i]);
                lb_record[i] += lb_temp;
            }else if (A[k][i]>Q_U[k][i]) {
                lb_temp = dist(A[k][i],Q_U[k][i]);
                lb_record[i] += lb_temp;
            }
            lb+=lb_temp;
        }
        //准备二次计算的资源
        lower_upper_lemire(A[k],m,r,A_L,A_U);
        for (int i=0;i<m;i++) {
            if (Q[k][i]<A_L[i]) {
                delta[i] += dist(Q[k][i],A_L[i]);
            }else if (Q[k][i]>A_U[i]) {
                delta[i] += dist(Q[k][i],A_U[i]);
            }
        }
    }
    upper_lemire(lb_record,m,r,max_val);
    for (int i=0;i<m;i++) {
        lb+=MAX(0.0,delta[i]-max_val[i]);
    }
    return sqrt(lb);
}

double LB_MV_double_KIM_T(vector<vector<double>>& A,vector<vector<double>>& Q, vector<vector<double>>& Q_L,vector<vector<double>>& Q_U,int r) {
    int m = A[0].size();
    int K = A.size();
    double lb = 0.0;
    vector<double>A_L(m),A_U(m);
    // vector<pair<int,double>>A_L(m),A_U(m);
    vector<double>max_val(m);
    vector<double> lb_record(m,0);
    for (int k=0;k<K;k++) {
        lb_record.assign(m, 0);
        for (int i=3;i<m-3;i++) {
            if (A[k][i]<Q_L[k][i]) {
                lb_record[i] = dist(A[k][i],Q_L[k][i]);
            }else if (A[k][i]>Q_U[k][i]) {
                lb_record[i] = dist(A[k][i],Q_U[k][i]);
            }
            lb+=lb_record[i];
        }

        lb+=LB_KIM_NEW(Q[k], A[k], m);

        //二次计算
        lower_upper_lemire(A[k],m,r,A_L,A_U);
        upper_lemire(lb_record,m,r,max_val);
        for (int i=3;i<m-3;i++) {
            double delta = 0;
            if (Q[k][i]<A_L[i]) {
                delta = dist(Q[k][i],A_L[i]);
            }else if (Q[k][i]>A_U[i]) {
                delta = dist(Q[k][i],A_U[i]);
            }
            lb+=MAX(0.0,delta-max_val[i]);
        }
    }

    return sqrt(lb);
}

// double LB_MV_double_IDEA_T(vector<vector<double>>& A,vector<vector<double>>& Q, vector<vector<double>>& Q_L,vector<vector<double>>& Q_U,int r) {
//     int m = A[0].size();
//     int K = A.size();
//     double lb = 0.0;
//     vector<double> lb_record(m,0);
//     for (int k=0;k<K;k++) {
//         lb_record.assign(m, 0);
//         for (int i=3;i<m-3;i++) {
//             if (A[k][i]<Q_L[k][i]) {
//                 lb_record[i] = dist(A[k][i],Q_L[k][i]);
//             }else if (A[k][i]>Q_U[k][i]) {
//                 lb_record[i] = dist(A[k][i],Q_U[k][i]);
//             }
//             lb+=lb_record[i];
//         }
//         lb+=LB_KIM_NEW(Q[k], A[k], m);
//         for (int i=3;i<m-3;i++) {
//             const int start = MAX(0,i-r);
//             const int end = MIN(m-1,i+r);
//             double min_val = MAX_VAL;
//             for (int j=start;j<=end;j++) {
//                 double distance = dist(Q[k][i],A[k][j]) - lb_record[j];
//                 min_val = MIN(min_val,distance);
//             }
//             lb+=MAX(0.0,min_val);
//         }
//     }
//
//     return sqrt(lb);
// }

double LB_MV_double_IDEA_T(vector<vector<double>>& A,vector<vector<double>>& Q, vector<vector<double>>& Q_L,vector<vector<double>>& Q_U,int r) {
    int m = A[0].size();
    int K = A.size();
    double lb = 0.0;
    vector<double> lb_record(m,0);
    for (int k=0;k<K;k++) {
        lb_record.assign(m, 0);
        for (int i=0;i<m;i++) {
            if (A[k][i]<Q_L[k][i]) {
                lb_record[i] = dist(A[k][i],Q_L[k][i]);
            }else if (A[k][i]>Q_U[k][i]) {
                lb_record[i] = dist(A[k][i],Q_U[k][i]);
            }
            lb+=lb_record[i];
        }
        // lb+=LB_KIM_NEW(Q[k], A[k], m);
        for (int i=0;i<m;i++) {
            const int start = MAX(0,i-r);
            const int end = MIN(m-1,i+r);
            double min_val = MAX_VAL;
            for (int j=start;j<=end;j++) {
                double distance = dist(Q[k][i],A[k][j]) - lb_record[j];
                // min_val = MIN(min_val,distance);
                min_val = MAX(min_val,distance);
            }
            lb+=MAX(0.0,min_val);
        }
    }

    return sqrt(lb);
}


double LB_MV_double_IDEA_T_V2(vector<vector<double>>& A,vector<vector<double>>& Q, vector<vector<double>>& A_t,vector<vector<double>>& Q_t,vector<vector<double>>& Q_L,vector<vector<double>>& Q_U,int r) {
    int m = A[0].size();
    int K = A.size();
    double lb = 0.0;
    double lb_temp = 0.0;
    vector<double> lb_record(m,0);
    for (int k=0;k<K;k++) {
        for (int i=0;i<m;i++) {
            lb_temp = 0;
            if (A[k][i]<Q_L[k][i]) {
                lb_temp = dist(A[k][i],Q_L[k][i]);
                lb_record[i] += lb_temp;
            }else if (A[k][i]>Q_U[k][i]) {
                lb_temp = dist(A[k][i],Q_U[k][i]);
                lb_record[i] += lb_temp;
            }
            lb+=lb_temp;
        }
        // lb+=LB_KIM_NEW(Q[k], A[k], m);
    }
        for (int i=0;i<m;i++) {
            const int start = MAX(0,i-r);
            const int end = MIN(m-1,i+r);
            double min_val = -MAX_VAL;
            for (int j=start;j<=end;j++) {
                double distance = euclideanDistance_2(Q_t[i],A_t[j]) - lb_record[j];
                // min_val = MIN(min_val,distance);
                min_val = MAX(min_val,distance);
            }
            lb+=MAX(0.0,min_val);
        }
    return sqrt(lb);
}

// double LB_MV_double_IDEA_T_V2(vector<vector<double>>& A,vector<vector<double>>& Q, vector<vector<double>>& A_t,vector<vector<double>>& Q_t,vector<vector<double>>& Q_L,vector<vector<double>>& Q_U,int r) {
//     int m = A[0].size();
//     int K = A.size();
//     double lb = 0.0;
//     double lb_temp = 0.0;
//     vector<double> lb_record(m,0);
//     for (int k=0;k<K;k++) {
//         for (int i=3;i<m-3;i++) {
//             lb_temp = 0;
//             if (A[k][i]<Q_L[k][i]) {
//                 lb_temp = dist(A[k][i],Q_L[k][i]);
//                 lb_record[i] += lb_temp;
//             }else if (A[k][i]>Q_U[k][i]) {
//                 lb_temp = dist(A[k][i],Q_U[k][i]);
//                 lb_record[i] += lb_temp;
//             }
//             lb+=lb_temp;
//         }
//         lb+=LB_KIM_NEW(Q[k], A[k], m);
//     }
//     for (int i=3;i<m-3;i++) {
//         const int start = MAX(0,i-r);
//         const int end = MIN(m-1,i+r);
//         double min_val = MAX_VAL;
//         for (int j=start;j<=end;j++) {
//             double distance = euclideanDistance_2(Q_t[i],A_t[j]) - lb_record[j];
//             min_val = MIN(min_val,distance);
//         }
//         lb+=MAX(0.0,min_val);
//     }
//     return sqrt(lb);
// }

double LB_P_MV(vector<vector<double>>& A,vector<vector<double>>& Q, vector<vector<double>>& Q_L,vector<vector<double>>& Q_U,int r) {
    int m = A[0].size();
    int K = A.size();
    double lb = 0.0;
    vector<double>A_L(m),A_U(m);
    for (int k=0;k<K;k++) {
        lower_upper_lemire(A[k],m,r,A_L,A_U);
        lb+=lbPetitjean_new(A[k],A_U,A_L,Q[k],Q_U[k],Q_L[k],r);
    }
    return sqrt(lb);
}

void LB_P_MV_early(vector<vector<double>>& A,vector<vector<double>>& Q, vector<vector<double>>& Q_L,vector<vector<double>>& Q_U,int r,double threshold,int &flag,double &lb) {
    int m = A[0].size();
    int K = A.size();
    vector<double>A_L(m),A_U(m);
    for (int k=0;k<K;k++) {
        lower_upper_lemire(A[k],m,r,A_L,A_U);
        flag = 0 ;
        lbPetitjean_new_early(A[k],A_U,A_L,Q[k],Q_U[k],Q_L[k],r,threshold*threshold,flag,lb);
        if (flag)break;
    }
    // return sqrt(lb);
}

void LB_P_MV_early(vector<vector<double>>& A,vector<vector<double>>& Q, vector<vector<double>>& Q_L,vector<vector<double>>& Q_U,int r,double threshold,int &flag,double &lb,vector<double>&cb) {
    int m = A[0].size();
    int K = A.size();
    vector<double>A_L(m),A_U(m);
    for (int k=0;k<K;k++) {
        lower_upper_lemire(A[k],m,r,A_L,A_U);
        lbPetitjean_new_early(A[k],A_U,A_L,Q[k],Q_U[k],Q_L[k],r,threshold,flag,lb,cb);
        if (flag)break;
    }
    // return sqrt(lb);
}


double LB_WEB(const std::vector<double>& q, vector<double>& uq, vector<double>& lq,const std::vector<double>& t,vector<double>& ut, vector<double>& lt,const std::vector<double>& lut,
              const std::vector<double>& ult, int window) {
    int m = q.size();
    vector<double> luq(m),ulq(m);
    lower_lemire(uq,m,window,luq);
    upper_lemire(lq,m,window,ulq);
    return lbWebb(q,uq,lq,luq,ulq,t,ut,lt,lut,ult,window);

}



double LB_WEB_MV(vector<vector<double>>& A,vector<vector<double>>& Q, vector<vector<double>>& Q_L,vector<vector<double>>& Q_U,vector<vector<double>>& L_UQ,vector<vector<double>>& U_LQ,int r) {
    int m = A[0].size();
    int K = A.size();
    double lb = 0.0;
    vector<double>A_L(m),A_U(m);
    for (int k=0;k<K;k++) {
        lower_upper_lemire(A[k],m,r,A_L,A_U);
        lb+=LB_WEB(A[k],A_U,A_L,Q[k],Q_U[k],Q_L[k],L_UQ[k],U_LQ[k],r);
    }
    return sqrt(lb);
}