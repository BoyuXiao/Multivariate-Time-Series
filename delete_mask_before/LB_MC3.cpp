#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>

#include "funcs.h"
using namespace std;

int H(double X, double LQ,double UQ, double intervalLength,int lambda,int n) {
    if (X<=LQ)return 0;
    else if (X>=UQ)return lambda*n-1;
    else return static_cast<int>((X - LQ) / intervalLength);
}

double initializeT(std::vector<Bucket>& T, std::vector<double>& Q, double LQ, double UQ, int lambda, int n) {
    double intervalLength = (UQ - LQ) / (lambda * n);
    // 初始化所有的 Bucket
    for (int i = 0; i < lambda * n; ++i) {
        double left = LQ + i * intervalLength;
        double right = left + intervalLength;
        T.emplace_back(left, right);
    }
    int T_size = T.size();
    // 有序地将所有点插入到对应的 Bucket 中
    for (double val : Q) {
        int index = H(val, LQ, UQ, intervalLength, lambda, n);
        if (index >= 0 && index < T_size) {
            // 调用 insertData 函数将数据有序插入到对应的 Bucket 中
            T[index].dataList.push_back(val);
        }
    }

    //left
    for (int i = 0; i < T_size; ++i) {
        int j = i-1;
        if (j<0)continue;
        if (T[j].leftNearestDistance!=MAX_VAL) {
            T[i].leftNearestDistance = T[j].leftNearestDistance + intervalLength;
        }
        if (!T[j].dataList.empty()) {
            T[i].leftNearestDistance = T[i].left-T[j].max_val();
        }
    }

    //right
    for (int i = T_size-1; i >= 0; --i) {
        int j = i+1;
        if (j>=T_size)continue;
        if (T[j].rightNearestDistance!=MAX_VAL) {
            T[i].rightNearestDistance = T[j].rightNearestDistance + intervalLength;
        }
        if (!T[j].dataList.empty()) {
            T[i].rightNearestDistance=T[j].min_val()-T[i].right;
        }
    }

    return intervalLength;
}

double CornerDist(const vector<double>& X, const vector<double>& V, double dNN) {
    double dMu = (X[0] - V[0])*(X[0] - V[0]);
    double dSigma = (X[1] - V[1])*(X[1] - V[1]);
    return dMu + dSigma + dNN*dNN + 2 * dNN * sqrt(MIN(dMu, dSigma));
}


double CenterDist(double X, double LYi, double UYi, vector<Bucket>& T,double intervalLength,int lambda,int num) {
    double dist = 0.0;
    int t = H(X,LYi,UYi,intervalLength,lambda,num);
    if (T[t].dataList.size() >= 2) {
        return dist;
    } else {
        double leftDist=MAX_VAL;
        double rightDist=MAX_VAL;
        leftDist = MIN((X-T[t].left+T[t].leftNearestDistance)*(X-T[t].left+T[t].leftNearestDistance),leftDist);
        rightDist = MIN((T[t].right-X +T[t].rightNearestDistance)*(T[t].right-X +T[t].rightNearestDistance),rightDist);
        double epsilon = MIN(leftDist, rightDist);
        if (T[t].dataList.size()== 1)
            epsilon = MIN(epsilon, (X-T[t].dataList[0])*(X-T[t].dataList[0]));
        dist += epsilon;
    }
    // cout<<dist<<endl;
    return dist;
}


double EdgeDist(const vector<double>& X, const vector<double>& V, double dNN,bool MiddleLeftOrRight,bool TopOrBottomMiddle,double L_mu_Q, double U_mu_Q,
    double L_sigma_Q, double U_sigma_Q,vector<Bucket>& T_mu, vector<Bucket>& T_sigma,double intervalLength_mu,double intervalLength_sigma,int lambda,int num) {
    double dMu = (X[0] - V[0])*(X[0] - V[0]);
    double dSigma = (X[1] - V[1])*(X[1] - V[1]);
    if (MiddleLeftOrRight) {
        if (dNN*dNN <= dSigma) {
            return dMu+CenterDist(X[1],L_sigma_Q,U_sigma_Q,T_sigma,intervalLength_sigma,lambda,num);
        } else {
            return MAX(dMu+CenterDist(X[1],L_sigma_Q,U_sigma_Q,T_sigma,intervalLength_sigma,lambda,num),dNN*dNN + dMu + dSigma - 2 * dNN * sqrt(dSigma));
        }
    } else if (TopOrBottomMiddle) {
        if (dNN*dNN <= dMu) {
            return dSigma+CenterDist(X[0],L_mu_Q,U_mu_Q,T_mu,intervalLength_mu,lambda,num);
        } else {
            return MAX(dSigma+CenterDist(X[0],L_mu_Q,U_mu_Q,T_mu,intervalLength_mu,lambda,num),dNN*dNN + dMu + dSigma - 2 * dNN * sqrt(dMu));
        }
    }

    return -1;
}



vector<vector<double>> Preprocessing_MC3(vector<vector<double>>& Q,
    vector<double>& mu_Q,vector<double>& sigma_Q,
    vector<double>& U_mu_Q,vector<double>& L_mu_Q,
    vector<double>& U_sigma_Q,vector<double>& L_sigma_Q,
    vector<vector<Bucket>>& T_mu,vector<vector<Bucket>>& T_sigma,
    vector <double>& intervalLength_mu,vector <double>& intervalLength_sigma ,vector <int>& num,int r) {
    int lambda =3;
    int m = Q.size();
    lower_upper_lemire(mu_Q, m, r, L_mu_Q, U_mu_Q);
    lower_upper_lemire(sigma_Q, m, r, L_sigma_Q, U_sigma_Q);

    vector<vector<double>> result(4, vector<double>(m));


    double minDbl,minDbr,minDtl,minDtr;
    vector<int> start(m);
    vector<int> end(m);
    vector<double> Qt(2);
    vector<double> corner(2);
    double dist;

    for (int i = 0; i < m; ++i) {
        start[i] = MAX(0,i-r);
        end[i] = MIN(m-1,i+r);
    }
    //计算四个dNN预选值
    for (int i = 0; i < m; ++i) {
        minDbl = MAX_VAL;
        for (int t = start[i]; t <= end[i]; ++t) {
            Qt = {mu_Q[t], sigma_Q[t]};
            corner = {L_mu_Q[i], L_sigma_Q[i]};
            dist = euclideanDistance_2(Qt, corner,Qt.size());
            minDbl = MIN(minDbl, dist);
        }
        result[0][i] = minDbl;
    }
    for (int i = 0; i < m; ++i) {
        minDbr = MAX_VAL;
        for (int t = start[i]; t <= end[i]; ++t) {
            Qt = {mu_Q[t], sigma_Q[t]};
            corner = {U_mu_Q[i], L_sigma_Q[i]};
            dist = euclideanDistance_2(Qt, corner,Qt.size());
            minDbr = MIN(minDbr, dist);
        }
        result[1][i] = minDbr;
    }
    for (int i = 0; i < m; ++i) {
        minDtl = MAX_VAL;
        for (int t = start[i]; t <= end[i]; ++t) {
            Qt = {mu_Q[t], sigma_Q[t]};
            corner = {L_mu_Q[i], U_sigma_Q[i]};
            dist = euclideanDistance_2(Qt, corner,Qt.size());
            minDtl = MIN(minDtl, dist);
        }
        result[2][i] = minDtl;
    }
    for (int i = 0; i < m; ++i) {
        minDtr = MAX_VAL;
        for (int t = start[i]; t <= end[i]; ++t) {
            Qt = {mu_Q[t], sigma_Q[t]};
            corner = {U_mu_Q[i], U_sigma_Q[i]};
            dist = euclideanDistance_2(Qt, corner,Qt.size());
            minDtr = MIN(minDtr, dist);
        }
        result[3][i] = minDtr;
    }

    //计算Q的m个T表
    for (int i = 0; i < m; ++i) {
        num[i] = end[i]-start[i]+1;
        vector<double> mu_Q_part(mu_Q.begin() + start[i], mu_Q.begin() + end[i] + 1);
        vector<double> sigma_Q_part(sigma_Q.begin() + start[i], sigma_Q.begin() + end[i] + 1);
        intervalLength_mu[i] = initializeT(T_mu[i], mu_Q_part, L_mu_Q[i], U_mu_Q[i], lambda, num[i]);
        intervalLength_sigma[i] = initializeT(T_sigma[i], sigma_Q_part, L_sigma_Q[i], U_sigma_Q[i], lambda, num[i]);
    }

    for (int i = 0; i < result.size(); ++i) {
        for (int j = 0; j < result[i].size(); ++j) {
            result[i][j] = sqrt(result[i][j]);
        }
    }

    return result;
}



bool isCornerCase(vector<double>& Ai, double& L_mu_Q, double& U_mu_Q, double& L_sigma_Q, double& U_sigma_Q) {
    return (Ai[0] < L_mu_Q && Ai[1] < L_sigma_Q) ||
           (Ai[0] > U_mu_Q && Ai[1] < L_sigma_Q) ||
           (Ai[0] < L_mu_Q && Ai[1] > U_sigma_Q) ||
           (Ai[0] > U_mu_Q && Ai[1] > U_sigma_Q);
}

bool isEdgeCase(vector<double>& Ai, double& L_mu_Q, double& U_mu_Q, double& L_sigma_Q, double& U_sigma_Q) {
    return ((Ai[0] < L_mu_Q || Ai[0] > U_mu_Q) && (L_sigma_Q <= Ai[1] && Ai[1] <= U_sigma_Q)) ||
           ((L_mu_Q <= Ai[0] && Ai[0] <= U_mu_Q) && (Ai[1] < L_sigma_Q || Ai[1] > U_sigma_Q));
}


double LB_MC3(vector<vector<double>>& Q, vector<vector<double>>& A,
    vector<double>& mu_Q,vector<double>& sigma_Q,
    vector<double>& U_mu_Q,vector<double>& L_mu_Q,
    vector<double>& U_sigma_Q,vector<double>& L_sigma_Q,
    vector<double>& Dbl, vector<double>& Dbr,
    vector<double>& Dtl, vector<double>& Dtr,
    vector<vector<Bucket>>& T_mu,vector<vector<Bucket>>& T_sigma,
    vector <double>& intervalLength_mu,vector <double>& intervalLength_sigma ,vector <int>& num,
    int K)
{

    std::ofstream outputFile_data;
    outputFile_data.open("..//data.txt",std::ios::app);

    int m = A.size();
    int lambda  = 3;
    vector<double> mu_A(m);
    vector<double> sigma_A(m);
    // mu_sigma(A,K,m,mu_A,sigma_A);
    if (K==2) {
        split_two_dim(A,mu_A,sigma_A);
    }else {
        mu_sigma(A,mu_A,sigma_A);
    }
    double lb = 0.0;

    vector<double> Ai(2);
    vector<double> V(2);
    double dNN;
    vector<double> V_plus(2);
    vector<double> V_minus(2);
    double d_plus;
    double d_minus;
    bool MiddleLeftOrRight=false;
    bool TopOrBottomMiddle=false;


    for (int i = 0; i < m; ++i) {
        Ai = {mu_A[i],sigma_A[i]};
        if (isCornerCase(Ai, L_mu_Q[i], U_mu_Q[i], L_sigma_Q[i], U_sigma_Q[i])) {
            if (Ai[0] < L_mu_Q[i] && Ai[1] < L_sigma_Q[i]) {
                V = {L_mu_Q[i], L_sigma_Q[i]};
                dNN = Dbl[i];
            } else if (Ai[0] > U_mu_Q[i] && Ai[1] < L_sigma_Q[i]) {
                V = {U_mu_Q[i], L_sigma_Q[i]};
                dNN = Dbr[i];
            } else if (Ai[0] < L_mu_Q[i] && Ai[1] > U_sigma_Q[i]) {
                V = {L_mu_Q[i], U_sigma_Q[i]};
                dNN = Dtl[i];
            } else {
                V = {U_mu_Q[i], U_sigma_Q[i]};
                dNN = Dtr[i];
            }
            lb += CornerDist(Ai, V, dNN);
            // outputFile_data <<i<<" "<< CornerDist(Ai, V, dNN) << " " << "corner" << endl;
        }
        else if (isEdgeCase(Ai, L_mu_Q[i], U_mu_Q[i], L_sigma_Q[i], U_sigma_Q[i])) {
            MiddleLeftOrRight=false;
            TopOrBottomMiddle=false;
            if (Ai[0] < L_mu_Q[i] && L_sigma_Q[i] <= Ai[1] && Ai[1] <= U_sigma_Q[i]) {
                V_plus = {L_mu_Q[i], L_sigma_Q[i]};
                V_minus = {L_mu_Q[i], U_sigma_Q[i]};
                d_plus = Dbl[i];
                d_minus = Dtl[i];
                MiddleLeftOrRight = true;
            } else if (Ai[0] > U_mu_Q[i] && L_sigma_Q[i] <= Ai[1] && Ai[1] <= U_sigma_Q[i]) {
                V_plus = {U_mu_Q[i], L_sigma_Q[i]};
                V_minus = {U_mu_Q[i], U_sigma_Q[i]};
                d_plus = Dbr[i];
                d_minus = Dtr[i];
                MiddleLeftOrRight = true;
            } else if (L_mu_Q[i] <= Ai[0] && Ai[0] <= U_mu_Q[i] && Ai[1] < L_sigma_Q[i]) {
                V_plus = {L_mu_Q[i], L_sigma_Q[i]};
                V_minus = {U_mu_Q[i], L_sigma_Q[i]};
                d_plus = Dbl[i];
                d_minus = Dbr[i];
                TopOrBottomMiddle = true;
            } else {
                V_plus = {L_mu_Q[i], U_sigma_Q[i]};
                V_minus = {U_mu_Q[i], U_sigma_Q[i]};
                d_plus = Dtl[i];
                d_minus = Dtr[i];
                TopOrBottomMiddle = true;
            }
            double delta_plus = EdgeDist(Ai, V_plus, d_plus,MiddleLeftOrRight,TopOrBottomMiddle,L_mu_Q[i],U_mu_Q[i],L_sigma_Q[i],U_sigma_Q[i],T_mu[i],T_sigma[i],intervalLength_mu[i],intervalLength_sigma[i],lambda,num[i]);
            double delta_minus = EdgeDist(Ai, V_minus, d_minus,MiddleLeftOrRight,TopOrBottomMiddle,L_mu_Q[i],U_mu_Q[i],L_sigma_Q[i],U_sigma_Q[i],T_mu[i],T_sigma[i],intervalLength_mu[i],intervalLength_sigma[i],lambda,num[i]);
            lb += MIN(delta_plus, delta_minus);
            // outputFile_data <<i<<" "<< min(delta_plus, delta_minus) << " " << "edge" << endl;
            // edge++;
        } else {
            // 中心情况
            lb += (CenterDist(mu_A[i],L_mu_Q[i],U_mu_Q[i],T_mu[i],intervalLength_mu[i],lambda,num[i])
                +CenterDist(sigma_A[i],L_sigma_Q[i],U_sigma_Q[i],T_sigma[i],intervalLength_sigma[i],lambda,num[i]));
            //outputFile_data <<i<<" "<< (CenterDist(mu_A[i],L_mu_Q[i],U_mu_Q[i],T_mu[i],intervalLength_mu[i],lambda,num[i])
            //+CenterDist(sigma_A[i],L_sigma_Q[i],U_sigma_Q[i],T_sigma[i],intervalLength_sigma[i],lambda,num[i])) << " " << "center" << endl;
        }


    }
    if (K==2) {
        return sqrt(lb);
    }
    return sqrt(K*lb);
}
