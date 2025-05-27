
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <utility>

#include "funcs.h"
using namespace std;

//TODO min max 的宏定义 所有的max用10086

//新版哈希
pair<int,int> H(double mu_X,double sigma_X, double L_mu_Q,double U_mu_Q,double L_sigma_Q,double U_sigma_Q, double intervalLength_mu,double intervalLength_sigma ,int lambda,int n) {
    int numIntervals = lambda * static_cast<int>(std::sqrt(n));
    int muIndex;
    if (mu_X <= L_mu_Q) {
        muIndex = 0;
    } else if (mu_X >= U_mu_Q) {
        muIndex = numIntervals - 1;
    } else {
        muIndex = static_cast<int>((mu_X - L_mu_Q) / intervalLength_mu);
    }

    int sigmaIndex;
    if (sigma_X <= L_sigma_Q) {
        sigmaIndex = 0;
    } else if (sigma_X >= U_sigma_Q) {
        sigmaIndex = numIntervals - 1;
    } else {
        sigmaIndex = static_cast<int>((sigma_X - L_sigma_Q) / intervalLength_sigma);
    }

    return {sigmaIndex, muIndex};
}


//新版初始化
// pair<double,double> initializeT(std::vector<vector<Cell>>& T, std::pair<vector<double>,vector<double>>& Q, double L_mu_Q, double U_mu_Q, double L_sigma_Q, double U_sigma_Q, int lambda, int n) {
//     int numIntervals = lambda * static_cast<int>(std::sqrt(n));
//     double intervalLength_mu = (U_mu_Q - L_mu_Q) / numIntervals;
//     double intervalLength_sigma = (U_sigma_Q - L_sigma_Q) / numIntervals;
//
//     T.resize(numIntervals);
//     for (int i = 0; i < numIntervals; ++i) {
//         double sigmaBottom = L_sigma_Q + i * intervalLength_sigma;
//         double sigmaTop = sigmaBottom + intervalLength_sigma;
//         for (int j = 0; j < numIntervals; ++j) {
//             double muLeft = L_mu_Q + j * intervalLength_mu;
//             double muRight = muLeft + intervalLength_mu;
//             // 初始化线段顶点坐标
//             T[i].emplace_back(muLeft, muRight, sigmaBottom, sigmaTop);
//         }
//     }
//
//     for (int i = 0; i < n; ++i) {
//             double mu = Q.first[i];
//             double sigma = Q.second[i];
//             auto [muIndex, sigmaIndex] = H(mu, sigma, L_mu_Q, U_mu_Q, L_sigma_Q, U_sigma_Q,intervalLength_mu,intervalLength_sigma,lambda,n);
//             if (muIndex >= 0 && muIndex < numIntervals && sigmaIndex >= 0 && sigmaIndex < numIntervals) {
//                 T[muIndex][sigmaIndex].dataList.emplace_back(mu,sigma);
//             }
//         }
//
//     for (int i = 0; i < numIntervals; ++i) {
//         for (int j = 0; j < numIntervals; ++j) {
//             //Lleft:
//             for (int row=0;row< numIntervals;row++) {
//                 for (int col=0;col<j;col++) {
//                     if (row==i) {
//                         for (int k=0;k<T[row][col].dataList.size();k++) {
//                             T[i][j].haveleft=true;
//                             double dist = T[i][j].Lleft_bottom.first-T[row][col].dataList[k].first;
//                             T[i][j].minDistToLleft=MIN(T[i][j].minDistToLleft,dist);
//                         }
//                     }else if (row<i) {
//                         for (int k=0;k<T[row][col].dataList.size();k++) {
//                             T[i][j].haveleft=true;
//                             double dist = sqrt(dist(T[i][j].Lleft_bottom.first,T[row][col].dataList[k].first)+dist(T[i][j].Lleft_bottom.second,T[row][col].dataList[k].second));
//                             T[i][j].minDistToLleft=MIN(T[i][j].minDistToLleft,dist);
//                         }
//                     }else {
//                         for (int k=0;k<T[row][col].dataList.size();k++) {
//                             T[i][j].haveleft=true;
//                             double dist = sqrt(dist(T[i][j].Lleft_top.first,T[row][col].dataList[k].first)+dist(T[i][j].Lleft_top.second,T[row][col].dataList[k].second));
//                             T[i][j].minDistToLleft=MIN(T[i][j].minDistToLleft,dist);
//                         }
//                     }
//                     }
//                 }
//
//             //LRight:
//             for (int row=0;row<numIntervals;row++) {
//                 for (int col=j+1;col<numIntervals;col++) {
//                     if (row==i) {
//                         for (int k=0;k<T[row][col].dataList.size();k++) {
//                             T[i][j].haveright=true;
//                             double dist = T[row][col].dataList[k].first-T[i][j].Lright_bottom.first;
//                             T[i][j].minDistToLright=MIN(T[i][j].minDistToLright,dist);
//                         }
//                     }else if (row<i) {
//                         for (int k=0;k<T[row][col].dataList.size();k++) {
//                             T[i][j].haveright=true;
//                             double dist = sqrt(dist(T[i][j].Lright_bottom.first,T[row][col].dataList[k].first)+dist(T[i][j].Lright_bottom.second,T[row][col].dataList[k].second));
//                             T[i][j].minDistToLright=MIN(T[i][j].minDistToLright,dist);
//                         }
//                     }else {
//                         for (int k=0;k<T[row][col].dataList.size();k++) {
//                             T[i][j].haveright=true;
//                             double dist = sqrt(dist(T[i][j].Lright_top.first,T[row][col].dataList[k].first)+dist(T[i][j].Lright_top.second,T[row][col].dataList[k].second));
//                             T[i][j].minDistToLright=MIN(T[i][j].minDistToLright,dist);
//                         }
//                     }
//                 }
//             }
//
//             //LBottom:
//             for (int row=0;row<i;row++) {
//                 for (int col=0;col<numIntervals;col++) {
//                     if (col==j) {
//                         for (int k=0;k<T[row][col].dataList.size();k++) {
//                             T[i][j].havebottom=true;
//                             double dist = T[i][j].Lleft_bottom.second-T[row][col].dataList[k].second;
//                             T[i][j].minDistToLbottom=MIN(T[i][j].minDistToLbottom,dist);
//                         }
//                     }else if (col<j) {
//                         for (int k=0;k<T[row][col].dataList.size();k++) {
//                             T[i][j].havebottom=true;
//                             double dist = sqrt(dist(T[i][j].Lleft_bottom.first,T[row][col].dataList[k].first)+dist(T[i][j].Lleft_bottom.second,T[row][col].dataList[k].second));
//                             T[i][j].minDistToLbottom=MIN(T[i][j].minDistToLbottom,dist);
//                         }
//                     }else {
//                         for (int k=0;k<T[row][col].dataList.size();k++) {
//                             T[i][j].havebottom=true;
//                             double dist = sqrt(dist(T[i][j].Lright_bottom.first,T[row][col].dataList[k].first)+dist(T[i][j].Lright_bottom.second,T[row][col].dataList[k].second));
//                             T[i][j].minDistToLbottom=MIN(T[i][j].minDistToLbottom,dist);
//                         }
//                     }
//                 }
//             }
//
//             //LTop:
//             for (int row=i+1;row<numIntervals;row++) {
//                 for (int col=0;col<numIntervals;col++) {
//                     if (col==j) {
//                         for (int k=0;k<T[row][col].dataList.size();k++) {
//                             T[i][j].havetop=true;
//                             double dist = T[row][col].dataList[k].second-T[i][j].Lleft_top.second;
//                             T[i][j].minDistToLtop=MIN(T[i][j].minDistToLtop,dist);
//                         }
//                     }else if (col<j) {
//                         for (int k=0;k<T[row][col].dataList.size();k++) {
//                             T[i][j].havetop=true;
//                             double dist = sqrt(dist(T[i][j].Lleft_top.first,T[row][col].dataList[k].first)+dist(T[i][j].Lleft_top.second,T[row][col].dataList[k].second));
//                             T[i][j].minDistToLtop=MIN(T[i][j].minDistToLtop,dist);
//                         }
//                     }else {
//                         for (int k=0;k<T[row][col].dataList.size();k++) {
//                             T[i][j].havetop=true;
//                             double dist = sqrt(dist(T[i][j].Lright_top.first,T[row][col].dataList[k].first)+dist(T[i][j].Lright_top.second,T[row][col].dataList[k].second));
//                             T[i][j].minDistToLtop=MIN(T[i][j].minDistToLtop,dist);
//                         }
//                     }
//                 }
//             }
//             }
//         }
//
//     return {intervalLength_mu, intervalLength_sigma};
// }


//优化版本
pair<double,double> initializeT(std::vector<vector<Cell>>& T, std::pair<vector<double>,vector<double>>& Q, double L_mu_Q, double U_mu_Q, double L_sigma_Q, double U_sigma_Q, int lambda, int n) {
    int numIntervals = lambda * static_cast<int>(std::sqrt(n));
    double intervalLength_mu = (U_mu_Q - L_mu_Q) / numIntervals;
    double intervalLength_sigma = (U_sigma_Q - L_sigma_Q) / numIntervals;

    T.resize(numIntervals);
    for (int i = 0; i < numIntervals; ++i) {
        double sigmaBottom = L_sigma_Q + i * intervalLength_sigma;
        double sigmaTop = sigmaBottom + intervalLength_sigma;
        for (int j = 0; j < numIntervals; ++j) {
            double muLeft = L_mu_Q + j * intervalLength_mu;
            double muRight = muLeft + intervalLength_mu;
            T[i].emplace_back(muLeft, muRight, sigmaBottom, sigmaTop);
        }
    }

    for (int i = 0; i < n; ++i) {
            double mu = Q.first[i];
            double sigma = Q.second[i];
            auto [sigmaIndex, muIndex] = H(mu, sigma, L_mu_Q, U_mu_Q, L_sigma_Q, U_sigma_Q,intervalLength_mu,intervalLength_sigma,lambda,n);
            if (muIndex >= 0 && muIndex < numIntervals && sigmaIndex >= 0 && sigmaIndex < numIntervals) {
                T[sigmaIndex][muIndex].dataList.emplace_back(mu,sigma);
            }
        }

    //Lleft:
    for (int i = 0; i < numIntervals; ++i) {
        for (int j = 0; j < numIntervals; ++j) {
            //继承上一个
            if (T[i][MAX(j-1,0)].minDistToLleft!=MAX_VAL) {
                double mu = T[i][MAX(j-1,0)].left_nearest_mu;
                double sigma = T[i][MAX(j-1,0)].left_nearest_sigma;
                auto [sigmaIndex, muIndex] = H(mu,sigma, L_mu_Q, U_mu_Q, L_sigma_Q, U_sigma_Q,intervalLength_mu,intervalLength_sigma,lambda,n);
                if (sigmaIndex==i) {
                    T[i][j].minDistToLleft = T[i][j].Lleft_bottom.first-mu;
                }else if (sigmaIndex<i) {
                    T[i][j].minDistToLleft = sqrt(dist(T[i][j].Lleft_bottom.first,mu)+dist(T[i][j].Lleft_bottom.second,sigma));
                }else {
                    T[i][j].minDistToLleft = sqrt(dist(T[i][j].Lleft_top.first,mu)+dist(T[i][j].Lleft_top.second,sigma));
                }
                T[i][j].left_nearest_mu = mu;
                T[i][j].left_nearest_sigma = sigma;
            }
            //left
            for (int row=0;row< numIntervals;row++) {
                int col = j-1;
                if (col<0)continue;
                    for (int k=0;k<T[row][col].dataList.size();k++) {
                        double dist;
                        double mu = T[row][col].dataList[k].first;
                        double sigma = T[row][col].dataList[k].second;
                        if (row==i) {
                            dist = T[i][j].Lleft_bottom.first-mu;
                        }else if (row<i) {
                            dist = sqrt(dist(T[i][j].Lleft_bottom.first,mu)+dist(T[i][j].Lleft_bottom.second,sigma));
                        }else {
                            dist = sqrt(dist(T[i][j].Lleft_top.first,mu)+dist(T[i][j].Lleft_top.second,sigma));
                        }
                        if (dist < T[i][j].minDistToLleft) {
                            T[i][j].minDistToLleft = dist;
                            T[i][j].left_nearest_mu = mu;
                            T[i][j].left_nearest_sigma = sigma;
                        }
                    }
                }
            }
        }
    //LRight:
    for (int i = 0; i < numIntervals; ++i) {
        for (int j = numIntervals-1; j >= 0; --j) {
            // 继承上一个
            if (T[i][MIN(j + 1, numIntervals - 1)].minDistToLright != MAX_VAL) {
                double mu = T[i][MIN(j + 1, numIntervals - 1)].right_nearest_mu;
                double sigma = T[i][MIN(j + 1, numIntervals - 1)].right_nearest_sigma;
                auto [sigmaIndex, muIndex] = H(mu, sigma, L_mu_Q, U_mu_Q, L_sigma_Q, U_sigma_Q, intervalLength_mu, intervalLength_sigma, lambda, n);
                if (sigmaIndex == i) {
                    T[i][j].minDistToLright = mu - T[i][j].Lright_bottom.first;
                } else if (sigmaIndex < i) {
                    T[i][j].minDistToLright = sqrt(dist(T[i][j].Lright_bottom.first, mu) + dist(T[i][j].Lright_bottom.second, sigma));
                } else {
                    T[i][j].minDistToLright = sqrt(dist(T[i][j].Lright_top.first, mu) + dist(T[i][j].Lright_top.second, sigma));
                }
                T[i][j].right_nearest_mu = mu;
                T[i][j].right_nearest_sigma = sigma;
            }
            //right
            for (int row=0;row<numIntervals;row++) {
                int col = j+1;
                if (col>=numIntervals)continue;
                    for (int k = 0; k < T[row][col].dataList.size(); k++) {
                        double dist;
                        double mu = T[row][col].dataList[k].first;
                        double sigma = T[row][col].dataList[k].second;
                        if (row == i) {
                            dist = mu - T[i][j].Lright_bottom.first;
                        } else if (row < i) {
                            dist = sqrt(dist(T[i][j].Lright_bottom.first, mu) + dist(T[i][j].Lright_bottom.second, sigma));
                        } else {
                            dist = sqrt(dist(T[i][j].Lright_top.first, mu) + dist(T[i][j].Lright_top.second, sigma));
                        }
                        if (dist < T[i][j].minDistToLright) {
                            T[i][j].minDistToLright = dist;
                            T[i][j].right_nearest_mu = mu;
                            T[i][j].right_nearest_sigma = sigma;
                        }
                    }
            }
        }
    }
    //LBottom:
    for (int i = 0; i < numIntervals; ++i) {
        for (int j = 0; j < numIntervals; ++j) {
            // 继承上一个
            if (T[MAX(i - 1, 0)][j].minDistToLbottom != MAX_VAL) {
                double mu = T[MAX(i - 1, 0)][j].bottom_nearest_mu;
                double sigma = T[MAX(i - 1, 0)][j].bottom_nearest_sigma;
                auto [sigmaIndex, muIndex] = H(mu, sigma, L_mu_Q, U_mu_Q, L_sigma_Q, U_sigma_Q, intervalLength_mu, intervalLength_sigma, lambda, n);
                if (sigmaIndex == i) {
                    T[i][j].minDistToLbottom = T[i][j].Lleft_bottom.second - sigma;
                } else if (sigmaIndex < i) {
                    T[i][j].minDistToLbottom = sqrt(dist(T[i][j].Lleft_bottom.first, mu) + dist(T[i][j].Lleft_bottom.second, sigma));
                } else {
                    T[i][j].minDistToLbottom = sqrt(dist(T[i][j].Lright_bottom.first, mu) + dist(T[i][j].Lright_bottom.second, sigma));
                }
                T[i][j].bottom_nearest_mu = mu;
                T[i][j].bottom_nearest_sigma = sigma;
            }
            //bottom
            int row = i-1;
            if (row<0)continue;
            for (int col=0;col<numIntervals;col++) {
                for (int k = 0; k < T[row][col].dataList.size(); k++) {
                    double dist;
                    double mu = T[row][col].dataList[k].first;
                    double sigma = T[row][col].dataList[k].second;
                    if (col == j) {
                        dist = T[i][j].Lleft_bottom.second - sigma;
                    } else if (col < j) {
                        dist = sqrt(dist(T[i][j].Lleft_bottom.first, mu) + dist(T[i][j].Lleft_bottom.second, sigma));
                    } else {
                        dist = sqrt(dist(T[i][j].Lright_bottom.first, mu) + dist(T[i][j].Lright_bottom.second, sigma));
                    }
                    if (dist < T[i][j].minDistToLbottom) {
                        T[i][j].minDistToLbottom = dist;
                        T[i][j].bottom_nearest_mu = mu;
                        T[i][j].bottom_nearest_sigma = sigma;
                    }
                }
            }
        }
    }
    //LTop:
    for (int i = numIntervals-1; i >= 0; --i) {
        for (int j = 0; j < numIntervals; ++j) {
            // 继承上一个
            if (T[MIN(i + 1, numIntervals - 1)][j].minDistToLtop != MAX_VAL) {
                double mu = T[MIN(i + 1, numIntervals - 1)][j].top_nearest_mu;
                double sigma = T[MIN(i + 1, numIntervals - 1)][j].top_nearest_sigma;
                auto [sigmaIndex, muIndex] = H(mu, sigma, L_mu_Q, U_mu_Q, L_sigma_Q, U_sigma_Q, intervalLength_mu, intervalLength_sigma, lambda, n);
                if (sigmaIndex == i) {
                    T[i][j].minDistToLtop = sigma - T[i][j].Lleft_top.second;
                } else if (sigmaIndex < i) {
                    T[i][j].minDistToLtop = sqrt(dist(T[i][j].Lleft_top.first, mu) + dist(T[i][j].Lleft_top.second, sigma));
                } else {
                    T[i][j].minDistToLtop = sqrt(dist(T[i][j].Lright_top.first, mu) + dist(T[i][j].Lright_top.second, sigma));
                }
                T[i][j].top_nearest_mu = mu;
                T[i][j].top_nearest_sigma = sigma;
            }
            //top
            int row = i+1;
            if (row>=numIntervals)continue;
            for (int col=0;col<numIntervals;col++) {
                for (int k = 0; k < T[row][col].dataList.size(); k++) {
                    double dist;
                    double mu = T[row][col].dataList[k].first;
                    double sigma = T[row][col].dataList[k].second;
                    if (col == j) {
                        dist = sigma - T[i][j].Lleft_top.second;
                    } else if (col < j) {
                        dist = sqrt(dist(T[i][j].Lleft_top.first, mu) + dist(T[i][j].Lleft_top.second, sigma));
                    } else {
                        dist = sqrt(dist(T[i][j].Lright_top.first, mu) + dist(T[i][j].Lright_top.second, sigma));
                    }
                    if (dist < T[i][j].minDistToLtop) {
                        T[i][j].minDistToLtop = dist;
                        T[i][j].top_nearest_mu = mu;
                        T[i][j].top_nearest_sigma = sigma;
                    }
                }
            }
        }
    }
    return {intervalLength_mu, intervalLength_sigma};
}


double CornerDistE(const vector<double>& X, const vector<double>& V, double dNN) {
    double dMu = (X[0] - V[0])*(X[0] - V[0]);
    double dSigma = (X[1] - V[1])*(X[1] - V[1]);
    return dMu + dSigma + dNN*dNN + 2 * dNN * sqrt(MIN(dMu, dSigma));
}


double CenterDist2D(double mu_X, double sigma_X, double L_mu_Q, double U_mu_Q, double L_sigma_Q, double U_sigma_Q, std::vector<vector<Cell>>& T,double intervalLength_mu,double intervalLength_sigma,int lambda,int num) {
    double dist = 0.0;
    auto [row,col] = H(mu_X,sigma_X,L_mu_Q,U_mu_Q,L_sigma_Q,U_sigma_Q,intervalLength_mu,intervalLength_sigma,lambda,num);
    if (T[row][col].dataList.size() >= 2) {
        return dist;
    } else {
        double leftDist=MAX_VAL;
        double rightDist=MAX_VAL;
        double bottomDist=MAX_VAL;
        double topDist=MAX_VAL;
        leftDist = MIN(leftDist,T[row][col].minDistToLleft+mu_X-T[row][col].Lleft_bottom.first)*(T[row][col].minDistToLleft+mu_X-T[row][col].Lleft_bottom.first);
        rightDist = MIN(rightDist,T[row][col].minDistToLright+T[row][col].Lright_bottom.first-mu_X)*(T[row][col].minDistToLright+T[row][col].Lright_bottom.first-mu_X);
        bottomDist = MIN(bottomDist,T[row][col].minDistToLbottom+sigma_X-T[row][col].Lleft_bottom.second)*(T[row][col].minDistToLbottom+sigma_X-T[row][col].Lleft_bottom.second);
        topDist = MIN(topDist,T[row][col].minDistToLtop+T[row][col].Lleft_top.second-sigma_X)*(T[row][col].minDistToLtop+T[row][col].Lleft_top.second-sigma_X);
        double epsilon = MIN(MIN(leftDist, rightDist),MIN(bottomDist,topDist));
        if (T[row][col].dataList.size()== 1)
            epsilon = MIN(epsilon, dist(mu_X,T[row][col].dataList[0].first)+dist(sigma_X,T[row][col].dataList[0].second));
        dist += epsilon;
    }
    return dist;
}


double EdgeDistE(const vector<double>& X, const vector<double>& V, double dNN,bool MiddleLeftOrRight,bool TopOrBottomMiddle,double L_mu_Q, double U_mu_Q,
    double L_sigma_Q, double U_sigma_Q,vector<vector<Cell>>& T,double intervalLength_mu,double intervalLength_sigma,int lambda,int num) {
    double dMu = (X[0] - V[0])*(X[0] - V[0]);
    double dSigma = (X[1] - V[1])*(X[1] - V[1]);
    if (MiddleLeftOrRight) {
        if (dNN*dNN <= dSigma) {
            return dMu + CenterDist2D(V[0],X[1],L_mu_Q,U_mu_Q,L_sigma_Q,U_sigma_Q,T,intervalLength_mu,intervalLength_sigma,lambda,num);
        } else {
            // return dNN*dNN + dMu + dSigma - 2 * dNN * sqrt(dSigma);
            return MAX(dMu + CenterDist2D(V[0],X[1],L_mu_Q,U_mu_Q,L_sigma_Q,U_sigma_Q,T,intervalLength_mu,intervalLength_sigma,lambda,num),dNN*dNN + dMu + dSigma - 2 * dNN * sqrt(dSigma));
        }
    } else if (TopOrBottomMiddle) {
        if (dNN*dNN <= dMu) {
            return dSigma + CenterDist2D(X[0],V[1],L_mu_Q,U_mu_Q,L_sigma_Q,U_sigma_Q,T,intervalLength_mu,intervalLength_sigma,lambda,num);
            // return dSigma+CenterDist(X[0],L_mu_Q,U_mu_Q,T_mu,intervalLength_mu,lambda,num);
        } else {
            // return dNN*dNN + dMu + dSigma - 2 * dNN * sqrt(dMu);
            return MAX(dSigma + CenterDist2D(X[0],V[1],L_mu_Q,U_mu_Q,L_sigma_Q,U_sigma_Q,T,intervalLength_mu,intervalLength_sigma,lambda,num),dNN*dNN + dMu + dSigma - 2 * dNN * sqrt(dMu));
        }
    }

    return -1;
}



vector<vector<double>> Preprocessing_MC3E(vector<vector<double>>& Q,
    vector<double>& mu_Q,vector<double>& sigma_Q,
    vector<double>& U_mu_Q,vector<double>& L_mu_Q,
    vector<double>& U_sigma_Q,vector<double>& L_sigma_Q,
    vector<vector<vector<Cell>>>& T,
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



    std::pair<double, double> intervalLength;
    //计算Q的m个T表
    for (int i = 0; i < m; ++i) {
        num[i] = end[i]-start[i]+1;
        vector<double> mu_Q_part(mu_Q.begin() + start[i], mu_Q.begin() + end[i] + 1);
        vector<double> sigma_Q_part(sigma_Q.begin() + start[i], sigma_Q.begin() + end[i] + 1);
        pair<vector<double>,vector<double>> Q_part = {mu_Q_part,sigma_Q_part};
        intervalLength = initializeT(T[i], Q_part, L_mu_Q[i], U_mu_Q[i], L_sigma_Q[i],U_sigma_Q[i],lambda, num[i]);
        intervalLength_mu[i]=intervalLength.first;
        intervalLength_sigma[i]=intervalLength.second;
    }

    for (int i = 0; i < result.size(); ++i) {
        for (int j = 0; j < result[i].size(); ++j) {
            result[i][j] = sqrt(result[i][j]);
        }
    }

    return result;
}



bool isCornerCaseE(vector<double>& Ai, double& L_mu_Q, double& U_mu_Q, double& L_sigma_Q, double& U_sigma_Q) {
    return (Ai[0] < L_mu_Q && Ai[1] < L_sigma_Q) ||
           (Ai[0] > U_mu_Q && Ai[1] < L_sigma_Q) ||
           (Ai[0] < L_mu_Q && Ai[1] > U_sigma_Q) ||
           (Ai[0] > U_mu_Q && Ai[1] > U_sigma_Q);
}

bool isEdgeCaseE(vector<double>& Ai, double& L_mu_Q, double& U_mu_Q, double& L_sigma_Q, double& U_sigma_Q) {
    return ((Ai[0] < L_mu_Q || Ai[0] > U_mu_Q) && (L_sigma_Q <= Ai[1] && Ai[1] <= U_sigma_Q)) ||
           ((L_mu_Q <= Ai[0] && Ai[0] <= U_mu_Q) && (Ai[1] < L_sigma_Q || Ai[1] > U_sigma_Q));
}

double LB_MC3E(vector<vector<double>>& Q, vector<vector<double>>& A,
    vector<double>& mu_Q,vector<double>& sigma_Q,
    vector<double>& U_mu_Q,vector<double>& L_mu_Q,
    vector<double>& U_sigma_Q,vector<double>& L_sigma_Q,
    vector<double>& Dbl, vector<double>& Dbr,
    vector<double>& Dtl, vector<double>& Dtr,
    vector<vector<vector<Cell>>>& T,
    vector <double>& intervalLength_mu,vector <double>& intervalLength_sigma ,vector <int>& num,
    int K)
{

    // std::ofstream outputFile_data;
    // outputFile_data.open("..//data.txt",std::ios::app);

    int m = A.size();
    int lambda  = 3;
    vector<double> mu_A(m);
    vector<double> sigma_A(m);
    mu_sigma(A,K,m,mu_A,sigma_A);
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
        if (isCornerCaseE(Ai, L_mu_Q[i], U_mu_Q[i], L_sigma_Q[i], U_sigma_Q[i])) {
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
            lb += CornerDistE(Ai, V, dNN);
            // outputFile_data <<i<<" "<< CornerDistE(Ai, V, dNN) << " " << "corner" << endl;
        }
        else if (isEdgeCaseE(Ai, L_mu_Q[i], U_mu_Q[i], L_sigma_Q[i], U_sigma_Q[i])) {
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
            double delta_plus = EdgeDistE(Ai, V_plus, d_plus,MiddleLeftOrRight,TopOrBottomMiddle,L_mu_Q[i],U_mu_Q[i],L_sigma_Q[i],U_sigma_Q[i],T[i],intervalLength_mu[i],intervalLength_sigma[i],lambda,num[i]);
            double delta_minus = EdgeDistE(Ai, V_minus, d_minus,MiddleLeftOrRight,TopOrBottomMiddle,L_mu_Q[i],U_mu_Q[i],L_sigma_Q[i],U_sigma_Q[i],T[i],intervalLength_mu[i],intervalLength_sigma[i],lambda,num[i]);
            lb += MIN(delta_plus, delta_minus);
            // outputFile_data <<i<<" "<< MIN(delta_plus, delta_minus) << " " << "edge" << endl;
        } else {
            // 中心情况
            lb += CenterDist2D(mu_A[i],sigma_A[i],L_mu_Q[i],U_mu_Q[i],L_sigma_Q[i],U_sigma_Q[i],T[i],intervalLength_mu[i],intervalLength_sigma[i],lambda,num[i]);
            // outputFile_data <<i<<" "<< CenterDist2D(mu_A[i],sigma_A[i],L_mu_Q[i],U_mu_Q[i],L_sigma_Q[i],U_sigma_Q[i],T[i],intervalLength_mu[i],intervalLength_sigma[i],lambda,num[i]) << " " << "center" << endl;
        }

    }

    return sqrt(K*lb);
}
