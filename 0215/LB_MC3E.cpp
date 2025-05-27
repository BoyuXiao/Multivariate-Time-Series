
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <utility>
#include "funcs.h"
using namespace std;

// 优化版本
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
      T[0][0].numIntervals_sigma = numIntervals;
      T[0][0].numIntervals_mu = numIntervals;
      intervalLength_mu = 1 / intervalLength_mu;
      intervalLength_sigma = 1 / intervalLength_sigma;
      for (int i = 0; i < n; ++i) {
              double mu = Q.first[i];
              double sigma = Q.second[i];
          // auto [sigmaIndex, muIndex] = H(mu, sigma, L_mu_Q, U_mu_Q, L_sigma_Q, U_sigma_Q,intervalLength_mu,intervalLength_sigma,lambda,n);
          int muIndex;
          int sigmaIndex;

          if (mu!=U_mu_Q) muIndex = static_cast<int>((mu - L_mu_Q) * intervalLength_mu);
          else muIndex = numIntervals -1;
          if (sigma!=U_sigma_Q) sigmaIndex = static_cast<int>((sigma - L_sigma_Q) * intervalLength_sigma);
          else sigmaIndex = numIntervals -1;
          // if (muIndex >= 0 && muIndex < numIntervals && sigmaIndex >= 0 && sigmaIndex < numIntervals) {
              T[sigmaIndex][muIndex].dataList.emplace_back(mu,sigma);
          // }
     }

      //leftbottom point and left,bottom edge
      for (int i = 0; i < numIntervals; ++i) {
          for (int j = 0; j < numIntervals; ++j) {
              //继承两个点
              int row = MAX(i-1,0);
              int col = MAX(j-1,0);
              if (T[row][col].left_top_nearest_dist1!=MAX_VAL) {
                  double mu = T[row][col].left_top_nearest_mu1;
                  double sigma = T[row][col].left_top_nearest_sigma1;
                  double dist = sqrt(dist(T[i][j].Lleft_bottom.first,mu)+dist(T[i][j].Lleft_bottom.second,sigma));
                  T[i][j].left_bottom_nearest_dist1 = dist;
                  T[i][j].left_bottom_nearest_mu1 = mu;
                  T[i][j].left_bottom_nearest_sigma1 = sigma;
                  if (i>=1&&j>=1) {
                      T[i-1][j-1].right_top_nearest_dist1 = dist;
                      T[i-1][j-1].right_top_nearest_mu1 = mu;
                      T[i-1][j-1].right_top_nearest_sigma1 = sigma;
                      T[i-1][j].left_top_nearest_dist1 = dist;
                      T[i-1][j].left_top_nearest_mu1 = mu;
                      T[i-1][j].left_top_nearest_sigma1 = sigma;
                      T[i][j-1].right_bottom_nearest_dist1 = dist;
                      T[i][j-1].right_bottom_nearest_mu1 = mu;
                      T[i][j-1].right_bottom_nearest_sigma1 = sigma;
                  }else if (i>=1) {
                      T[i-1][j].left_top_nearest_dist1 = dist;
                      T[i-1][j].left_top_nearest_mu1 = mu;
                      T[i-1][j].left_top_nearest_sigma1 = sigma;
                  }else if (j>=1) {
                      T[i][j-1].right_bottom_nearest_dist1 = dist;
                      T[i][j-1].right_bottom_nearest_mu1 = mu;
                      T[i][j-1].right_bottom_nearest_sigma1 = sigma;
                  }
              }
              if (T[row][col].right_bottom_nearest_dist1!=MAX_VAL) {
                  double mu = T[row][col].right_bottom_nearest_mu1;
                  double sigma = T[row][col].right_bottom_nearest_sigma1;
                  double dist = sqrt(dist(T[i][j].Lleft_bottom.first,mu)+dist(T[i][j].Lleft_bottom.second,sigma));
                  if (dist < T[i][j].left_bottom_nearest_dist1) {
                      T[i][j].left_bottom_nearest_dist1 = dist;
                      T[i][j].left_bottom_nearest_mu1 = mu;
                      T[i][j].left_bottom_nearest_sigma1 = sigma;
                      if (i>=1&&j>=1) {
                          T[i-1][j-1].right_top_nearest_dist1 = dist;
                          T[i-1][j-1].right_top_nearest_mu1 = mu;
                          T[i-1][j-1].right_top_nearest_sigma1 = sigma;
                          T[i-1][j].left_top_nearest_dist1 = dist;
                          T[i-1][j].left_top_nearest_mu1 = mu;
                          T[i-1][j].left_top_nearest_sigma1 = sigma;
                          T[i][j-1].right_bottom_nearest_dist1 = dist;
                          T[i][j-1].right_bottom_nearest_mu1 = mu;
                          T[i][j-1].right_bottom_nearest_sigma1 = sigma;
                      }else if (i>=1) {
                          T[i-1][j].left_top_nearest_dist1 = dist;
                          T[i-1][j].left_top_nearest_mu1 = mu;
                          T[i-1][j].left_top_nearest_sigma1 = sigma;
                      }else if (j>=1) {
                          T[i][j-1].right_bottom_nearest_dist1 = dist;
                          T[i][j-1].right_bottom_nearest_mu1 = mu;
                          T[i][j-1].right_bottom_nearest_sigma1 = sigma;
                      }
                  }
              }
              //遍历一个格子来给点赋值
              for (int k=0;k<T[row][col].dataList.size();k++) {
                  if (i<1||j<1) break;
                  double mu = T[row][col].dataList[k].first;
                  double sigma = T[row][col].dataList[k].second;
                  double dist = sqrt(dist(T[i][j].Lleft_bottom.first,mu)+dist(T[i][j].Lleft_bottom.second,sigma));
                  if (dist < T[i][j].left_bottom_nearest_dist1) {
                      T[i][j].left_bottom_nearest_dist1 = dist;
                      T[i][j].left_bottom_nearest_mu1 = mu;
                      T[i][j].left_bottom_nearest_sigma1 = sigma;
                      if (i>=1&&j>=1) {
                          T[i-1][j-1].right_top_nearest_dist1 = dist;
                          T[i-1][j-1].right_top_nearest_mu1 = mu;
                          T[i-1][j-1].right_top_nearest_sigma1 = sigma;
                          T[i-1][j].left_top_nearest_dist1 = dist;
                          T[i-1][j].left_top_nearest_mu1 = mu;
                          T[i-1][j].left_top_nearest_sigma1 = sigma;
                          T[i][j-1].right_bottom_nearest_dist1 = dist;
                          T[i][j-1].right_bottom_nearest_mu1 = mu;
                          T[i][j-1].right_bottom_nearest_sigma1 = sigma;
                      }else if (i>=1) {
                          T[i-1][j].left_top_nearest_dist1 = dist;
                          T[i-1][j].left_top_nearest_mu1 = mu;
                          T[i-1][j].left_top_nearest_sigma1 = sigma;
                      }else if (j>=1) {
                          T[i][j-1].right_bottom_nearest_dist1 = dist;
                          T[i][j-1].right_bottom_nearest_mu1 = mu;
                          T[i][j-1].right_bottom_nearest_sigma1 = sigma;
                      }
                  }
              }
              //遍历一个格子来给left edge赋值
              row = i;
              col = MAX(0,j-1);
              if (T[row][col].dataList.size()!=0) {
                  for (int k=0;k<T[row][col].dataList.size();k++) {
                      if (j<1)break;
                      double mu = T[row][col].dataList[k].first;
                      double sigma = T[row][col].dataList[k].second;
                      double dist = T[i][j].Lleft_bottom.first - mu;
                      if (dist < T[i][j].minDistToLleft) {
                          T[i][j].minDistToLleft = dist;
                          T[i][j].left_nearest_mu = mu;
                          T[i][j].left_nearest_sigma = sigma;
                      }
                  }
              }else {
                  if (T[row][col].minDistToLleft!=MAX_VAL) {
                      double mu = T[row][col].left_nearest_mu;
                      double sigma = T[row][col].left_nearest_sigma;
                      double dist = T[i][j].Lleft_bottom.first - mu;
                      T[i][j].minDistToLleft = dist;
                      T[i][j].left_nearest_mu = mu;
                      T[i][j].left_nearest_sigma = sigma;
                  }
              }
              row = MAX(i-1,0);
              col = j;
              //遍历一个格子来给bottom edge赋值
              if (T[row][col].dataList.size()!=0) {
                  for (int k=0;k<T[row][col].dataList.size();k++) {
                      if (i<1)break;
                      double mu = T[row][col].dataList[k].first;
                      double sigma = T[row][col].dataList[k].second;
                      double dist = T[i][j].Lright_bottom.second - sigma;
                      if (dist < T[i][j].minDistToLbottom) {
                          T[i][j].minDistToLbottom = dist;
                          T[i][j].bottom_nearest_mu = mu;
                          T[i][j].bottom_nearest_sigma = sigma;
                      }
                  }
              }else {
                  if (T[row][col].minDistToLbottom!=MAX_VAL) {
                      double mu = T[row][col].bottom_nearest_mu;
                      double sigma = T[row][col].bottom_nearest_sigma;
                      double dist = T[i][j].Lright_bottom.second - sigma;
                      T[i][j].minDistToLbottom = dist;
                      T[i][j].bottom_nearest_mu = mu;
                      T[i][j].bottom_nearest_sigma = sigma;
                  }
              }
          }
      }

      //righttop point and right,top edge
      for (int i = numIntervals - 1; i >= 0; --i) {
          for (int j = numIntervals - 1; j >= 0; --j) {
              //继承两个点
              int row = MIN(i+1,numIntervals-1);
              int col  = MIN(j+1,numIntervals - 1);
              if (T[row][col].left_top_nearest_dist2!=MAX_VAL) {
                  double mu = T[row][col].left_top_nearest_mu2;
                  double sigma = T[row][col].left_top_nearest_sigma2;
                  double dist = sqrt(dist(T[i][j].Lright_top.first,mu)+dist(T[i][j].Lright_top.second,sigma));
                  T[i][j].right_top_nearest_dist2 = dist;
                  T[i][j].right_top_nearest_mu2 = mu;
                  T[i][j].right_top_nearest_sigma2 = sigma;
                  if (i+2<=numIntervals&&j+2<=numIntervals) {
                      T[i+1][j+1].left_bottom_nearest_dist2 = dist;
                      T[i+1][j+1].left_bottom_nearest_mu2 = mu;
                      T[i+1][j+1].left_bottom_nearest_sigma2 = sigma;

                      T[i][j+1].left_top_nearest_dist2 = dist;
                      T[i][j+1].left_top_nearest_mu2 = mu;
                      T[i][j+1].left_top_nearest_sigma2 = sigma;

                      T[i+1][j].right_bottom_nearest_dist2 = dist;
                      T[i+1][j].right_bottom_nearest_mu2 = mu;
                      T[i+1][j].right_bottom_nearest_sigma2 = sigma;
                  }else if (j+2<=numIntervals) {
                      T[i][j+1].left_top_nearest_dist2 = dist;
                      T[i][j+1].left_top_nearest_mu2 = mu;
                      T[i][j+1].left_top_nearest_sigma2 = sigma;
                  }else if (i+2<=numIntervals) {
                      T[i+1][j].right_bottom_nearest_dist2 = dist;
                      T[i+1][j].right_bottom_nearest_mu2 = mu;
                      T[i+1][j].right_bottom_nearest_sigma2 = sigma;
                  }
              }
              if (T[row][col].right_bottom_nearest_dist2!=MAX_VAL) {
                  double mu = T[row][col].right_bottom_nearest_mu2;
                  double sigma = T[row][col].right_bottom_nearest_sigma2;
                  double dist = sqrt(dist(T[i][j].Lright_top.first,mu)+dist(T[i][j].Lright_top.second,sigma));
                  if (dist < T[i][j].right_top_nearest_dist2) {
                      T[i][j].right_top_nearest_dist2 = dist;
                      T[i][j].right_top_nearest_mu2 = mu;
                      T[i][j].right_top_nearest_sigma2 = sigma;
                      if (i+2<=numIntervals&&j+2<=numIntervals) {
                          T[i+1][j+1].left_bottom_nearest_dist2 = dist;
                          T[i+1][j+1].left_bottom_nearest_mu2 = mu;
                          T[i+1][j+1].left_bottom_nearest_sigma2 = sigma;

                          T[i][j+1].left_top_nearest_dist2 = dist;
                          T[i][j+1].left_top_nearest_mu2 = mu;
                          T[i][j+1].left_top_nearest_sigma2 = sigma;

                          T[i+1][j].right_bottom_nearest_dist2 = dist;
                          T[i+1][j].right_bottom_nearest_mu2 = mu;
                          T[i+1][j].right_bottom_nearest_sigma2 = sigma;
                      }else if (j+2<=numIntervals) {
                          T[i][j+1].left_top_nearest_dist2 = dist;
                          T[i][j+1].left_top_nearest_mu2 = mu;
                          T[i][j+1].left_top_nearest_sigma2 = sigma;
                      }else if (i+2<=numIntervals) {
                          T[i+1][j].right_bottom_nearest_dist2 = dist;
                          T[i+1][j].right_bottom_nearest_mu2 = mu;
                          T[i+1][j].right_bottom_nearest_sigma2 = sigma;
                      }
                  }
              }
              //遍历一个格子来给点赋值
              for (int k=0;k<T[row][col].dataList.size();k++) {
                  if (i+2>numIntervals||j+2>numIntervals)break;
                  double mu = T[row][col].dataList[k].first;
                  double sigma = T[row][col].dataList[k].second;
                  double dist = sqrt(dist(T[i][j].Lright_top.first,mu)+dist(T[i][j].Lright_top.second,sigma));
                  if (dist < T[i][j].right_top_nearest_dist2) {
                      T[i][j].right_top_nearest_dist2 = dist;
                      T[i][j].right_top_nearest_mu2 = mu;
                      T[i][j].right_top_nearest_sigma2 = sigma;
                      if (i+2<=numIntervals&&j+2<=numIntervals) {
                          T[i+1][j+1].left_bottom_nearest_dist2 = dist;
                          T[i+1][j+1].left_bottom_nearest_mu2 = mu;
                          T[i+1][j+1].left_bottom_nearest_sigma2 = sigma;

                          T[i][j+1].left_top_nearest_dist2 = dist;
                          T[i][j+1].left_top_nearest_mu2 = mu;
                          T[i][j+1].left_top_nearest_sigma2 = sigma;

                          T[i+1][j].right_bottom_nearest_dist2 = dist;
                          T[i+1][j].right_bottom_nearest_mu2 = mu;
                          T[i+1][j].right_bottom_nearest_sigma2 = sigma;
                      }else if (j+2<=numIntervals) {
                          T[i][j+1].left_top_nearest_dist2 = dist;
                          T[i][j+1].left_top_nearest_mu2 = mu;
                          T[i][j+1].left_top_nearest_sigma2 = sigma;
                      }else if (i+2<=numIntervals) {
                          T[i+1][j].right_bottom_nearest_dist2 = dist;
                          T[i+1][j].right_bottom_nearest_mu2 = mu;
                          T[i+1][j].right_bottom_nearest_sigma2 = sigma;
                      }
                  }
              }
              //遍历一个格子来给right edge赋值
              row = i;
              col = MIN(j+1,numIntervals - 1);
              if (T[row][col].dataList.size()!=0) {
                  for (int k=0;k<T[row][col].dataList.size();k++) {
                      if (j+2>numIntervals)break;
                      double mu = T[row][col].dataList[k].first;
                      double sigma = T[row][col].dataList[k].second;
                      double dist = mu - T[i][j].Lright_bottom.first;
                      if (dist < T[i][j].minDistToLright) {
                          T[i][j].minDistToLright = dist;
                          T[i][j].right_nearest_mu = mu;
                          T[i][j].right_nearest_sigma = sigma;
                      }
                  }
              }else {
                  if (T[row][col].minDistToLright!=MAX_VAL) {
                      double mu = T[row][col].right_nearest_mu;
                      double sigma = T[row][col].right_nearest_sigma;
                      double dist = mu - T[i][j].Lright_bottom.first;
                      T[i][j].minDistToLright = dist;
                      T[i][j].right_nearest_mu = mu;
                      T[i][j].right_nearest_sigma = sigma;
                  }
              }
              row = MIN(i+1,numIntervals-1);
              col = j;
              //遍历一个格子来给top edge赋值
              if (T[row][col].dataList.size()!=0) {
                  for (int k=0;k<T[row][col].dataList.size();k++) {
                      if (i+2>numIntervals)break;
                      double mu = T[row][col].dataList[k].first;
                      double sigma = T[row][col].dataList[k].second;
                      double dist = sigma - T[i][j].Lright_top.second;
                      if (dist < T[i][j].minDistToLtop) {
                          T[i][j].minDistToLtop = dist;
                          T[i][j].top_nearest_mu = mu;
                          T[i][j].top_nearest_sigma = sigma;
                      }
                  }
              }else {
                  if (T[row][col].minDistToLtop!=MAX_VAL) {
                      double mu = T[row][col].top_nearest_mu;
                      double sigma = T[row][col].top_nearest_sigma;
                      double dist = sigma - T[i][j].Lright_top.second;
                      T[i][j].minDistToLtop = dist;
                      T[i][j].top_nearest_mu = mu;
                      T[i][j].top_nearest_sigma = sigma;
                  }
              }
          }
      }


       //rightbottom point
       for (int i = 0; i < numIntervals; ++i) {
           for (int j = numIntervals - 1; j >= 0; --j) {
               //继承两个点
               int row = MAX(i-1,0);
               int col = MIN(j+1,numIntervals - 1);
               if (T[row][col].left_bottom_nearest_dist3!=MAX_VAL) {
                   double mu = T[row][col].left_bottom_nearest_mu3;
                   double sigma = T[row][col].left_bottom_nearest_sigma3;
                   double dist = sqrt(dist(T[i][j].Lright_bottom.first,mu)+dist(T[i][j].Lright_bottom.second,sigma));
                   T[i][j].right_bottom_nearest_dist3 = dist;
                   T[i][j].right_bottom_nearest_mu3 = mu;
                   T[i][j].right_bottom_nearest_sigma3 = sigma;
                   if (i>=1&&j+2<=numIntervals) {
                       T[i-1][j+1].left_top_nearest_dist3 = dist;
                       T[i-1][j+1].left_top_nearest_mu3 = mu;
                       T[i-1][j+1].left_top_nearest_sigma3 = sigma;
                       T[i-1][j].right_top_nearest_dist3 = dist;
                       T[i-1][j].right_top_nearest_mu3 = mu;
                       T[i-1][j].right_top_nearest_sigma3 = sigma;
                       T[i][j+1].left_bottom_nearest_dist3 = dist;
                       T[i][j+1].left_bottom_nearest_mu3 = mu;
                       T[i][j+1].left_bottom_nearest_sigma3 = sigma;
                   }else if (i>=1) {
                       T[i-1][j].right_top_nearest_dist3 = dist;
                       T[i-1][j].right_top_nearest_mu3 = mu;
                       T[i-1][j].right_top_nearest_sigma3 = sigma;
                   }else if (j+2<=numIntervals) {
                       T[i][j+1].left_bottom_nearest_dist3 = dist;
                       T[i][j+1].left_bottom_nearest_mu3 = mu;
                       T[i][j+1].left_bottom_nearest_sigma3 = sigma;
                   }
               }
               if (T[row][col].right_top_nearest_dist3!=MAX_VAL) {
                   double mu = T[row][col].right_top_nearest_mu3;
                   double sigma = T[row][col].right_top_nearest_sigma3;
                   double dist = sqrt(dist(T[i][j].Lright_bottom.first,mu)+dist(T[i][j].Lright_bottom.second,sigma));
                   if (dist < T[i][j].right_bottom_nearest_dist3) {
                       T[i][j].right_bottom_nearest_dist3 = dist;
                       T[i][j].right_bottom_nearest_mu3 = mu;
                       T[i][j].right_bottom_nearest_sigma3 = sigma;
                       if (i>=1&&j+2<=numIntervals) {
                           T[i-1][j+1].left_top_nearest_dist3 = dist;
                           T[i-1][j+1].left_top_nearest_mu3 = mu;
                           T[i-1][j+1].left_top_nearest_sigma3 = sigma;
                           T[i-1][j].right_top_nearest_dist3 = dist;
                           T[i-1][j].right_top_nearest_mu3 = mu;
                           T[i-1][j].right_top_nearest_sigma3 = sigma;
                           T[i][j+1].left_bottom_nearest_dist3 = dist;
                           T[i][j+1].left_bottom_nearest_mu3 = mu;
                           T[i][j+1].left_bottom_nearest_sigma3 = sigma;
                       }else if (i>=1) {
                           T[i-1][j].right_top_nearest_dist3 = dist;
                           T[i-1][j].right_top_nearest_mu3 = mu;
                           T[i-1][j].right_top_nearest_sigma3 = sigma;
                       }else if (j+2<=numIntervals) {
                           T[i][j+1].left_bottom_nearest_dist3 = dist;
                           T[i][j+1].left_bottom_nearest_mu3 = mu;
                           T[i][j+1].left_bottom_nearest_sigma3 = sigma;
                       }
                   }
               }
               //遍历一个格子来给点赋值
               for (int k=0;k<T[row][col].dataList.size();k++) {
                   if (i<1||j+2>numIntervals)break;
                   double mu = T[row][col].dataList[k].first;
                   double sigma = T[row][col].dataList[k].second;
                   double dist = sqrt(dist(T[i][j].Lright_bottom.first,mu)+dist(T[i][j].Lright_bottom.second,sigma));
                   if (dist < T[i][j].right_bottom_nearest_dist3) {
                       T[i][j].right_bottom_nearest_dist3 = dist;
                       T[i][j].right_bottom_nearest_mu3 = mu;
                       T[i][j].right_bottom_nearest_sigma3 = sigma;
                       if (i>=1&&j+2<=numIntervals) {
                           T[i-1][j+1].left_top_nearest_dist3 = dist;
                           T[i-1][j+1].left_top_nearest_mu3 = mu;
                           T[i-1][j+1].left_top_nearest_sigma3 = sigma;
                           T[i-1][j].right_top_nearest_dist3 = dist;
                           T[i-1][j].right_top_nearest_mu3 = mu;
                           T[i-1][j].right_top_nearest_sigma3 = sigma;
                           T[i][j+1].left_bottom_nearest_dist3 = dist;
                           T[i][j+1].left_bottom_nearest_mu3 = mu;
                           T[i][j+1].left_bottom_nearest_sigma3 = sigma;
                       }else if (i>=1) {
                           T[i-1][j].right_top_nearest_dist3 = dist;
                           T[i-1][j].right_top_nearest_mu3 = mu;
                           T[i-1][j].right_top_nearest_sigma3 = sigma;
                       }else if (j+2<=numIntervals) {
                           T[i][j+1].left_bottom_nearest_dist3 = dist;
                           T[i][j+1].left_bottom_nearest_mu3 = mu;
                           T[i][j+1].left_bottom_nearest_sigma3 = sigma;
                       }
                   }
               }
           }
       }

       //lefttop point
       for (int i = numIntervals - 1; i >= 0; --i) {
           for (int j = 0; j < numIntervals; ++j) {
               //继承两个点
               int row = MIN(i+1,numIntervals-1);
               int col = MAX(j-1,0);
               if (T[row][col].left_bottom_nearest_dist4!=MAX_VAL) {
                   double mu = T[row][col].left_bottom_nearest_mu4;
                   double sigma = T[row][col].left_bottom_nearest_sigma4;
                   double dist = sqrt(dist(T[i][j].Lleft_top.first,mu)+dist(T[i][j].Lleft_top.second,sigma));
                   T[i][j].left_top_nearest_dist4 = dist;
                   T[i][j].left_top_nearest_mu4 = mu;
                   T[i][j].left_top_nearest_sigma4 = sigma;
                   if (j>=1&&i+2<=numIntervals) {
                       T[i+1][j-1].right_bottom_nearest_dist4 = dist;
                       T[i+1][j-1].right_bottom_nearest_mu4 = mu;
                       T[i+1][j-1].right_bottom_nearest_sigma4 = sigma;
                   }else if (j>=1) {
                       T[i][j-1].right_top_nearest_dist4 = dist;
                       T[i][j-1].right_top_nearest_mu4 = mu;
                       T[i][j-1].right_top_nearest_sigma4 = sigma;
                   }else if (i+2<=numIntervals) {
                       T[i+1][j].left_bottom_nearest_dist4 = dist;
                       T[i+1][j].left_bottom_nearest_mu4 = mu;
                       T[i+1][j].left_bottom_nearest_sigma4 = sigma;
                   }
               }
               if (T[row][col].right_top_nearest_dist4!=MAX_VAL) {
                   double mu = T[row][col].right_top_nearest_mu4;
                   double sigma = T[row][col].right_top_nearest_sigma4;
                   double dist = sqrt(dist(T[i][j].Lleft_top.first,mu)+dist(T[i][j].Lleft_top.second,sigma));
                   if (dist < T[i][j].left_top_nearest_dist4) {
                       T[i][j].left_top_nearest_dist4 = dist;
                       T[i][j].left_top_nearest_mu4 = mu;
                       T[i][j].left_top_nearest_sigma4 = sigma;
                       if (j>=1&&i+2<=numIntervals) {
                           T[i+1][j-1].right_bottom_nearest_dist4 = dist;
                           T[i+1][j-1].right_bottom_nearest_mu4 = mu;
                           T[i+1][j-1].right_bottom_nearest_sigma4 = sigma;
                       }else if (j>=1) {
                           T[i][j-1].right_top_nearest_dist4 = dist;
                           T[i][j-1].right_top_nearest_mu4 = mu;
                           T[i][j-1].right_top_nearest_sigma4 = sigma;
                       }else if (i+2<=numIntervals) {
                           T[i+1][j].left_bottom_nearest_dist4 = dist;
                           T[i+1][j].left_bottom_nearest_mu4 = mu;
                           T[i+1][j].left_bottom_nearest_sigma4 = sigma;
                       }
                   }
               }
               //遍历一个格子来给点赋值
               for (int k=0;k<T[row][col].dataList.size();k++) {
                   if (j<1||i+2>numIntervals)break;
                   double mu = T[row][col].dataList[k].first;
                   double sigma = T[row][col].dataList[k].second;
                   double dist = sqrt(dist(T[i][j].Lleft_top.first,mu)+dist(T[i][j].Lleft_top.second,sigma));
                   if (dist < T[i][j].left_top_nearest_dist4) {
                       T[i][j].left_top_nearest_dist4 = dist;
                       T[i][j].left_top_nearest_mu4 = mu;
                       T[i][j].left_top_nearest_sigma4 = sigma;
                       if (j>=1&&i+2<=numIntervals) {
                           T[i+1][j-1].right_bottom_nearest_dist4 = dist;
                           T[i+1][j-1].right_bottom_nearest_mu4 = mu;
                           T[i+1][j-1].right_bottom_nearest_sigma4 = sigma;
                       }else if (j>=1) {
                           T[i][j-1].right_top_nearest_dist4 = dist;
                           T[i][j-1].right_top_nearest_mu4 = mu;
                           T[i][j-1].right_top_nearest_sigma4 = sigma;
                       }else if (i+2<=numIntervals) {
                           T[i+1][j].left_bottom_nearest_dist4 = dist;
                           T[i+1][j].left_bottom_nearest_mu4 = mu;
                           T[i+1][j].left_bottom_nearest_sigma4 = sigma;
                       }
                   }
               }
           }
       }

       //update
       for (int i = 0; i < numIntervals; ++i) {
           for (int j = 0; j < numIntervals; ++j) {
               T[i][j].minDistToLleft = MIN(MIN(T[i][j].minDistToLleft,T[i][j].left_bottom_nearest_dist1),T[i][j].left_top_nearest_dist4);
               T[i][j].minDistToLright = MIN(MIN(T[i][j].minDistToLright,T[i][j].right_bottom_nearest_dist3),T[i][j].right_top_nearest_dist2);
               T[i][j].minDistToLbottom = MIN(MIN(T[i][j].minDistToLbottom,T[i][j].left_bottom_nearest_dist1),T[i][j].right_bottom_nearest_dist3);
               T[i][j].minDistToLtop = MIN(MIN(T[i][j].minDistToLtop,T[i][j].left_top_nearest_dist4),T[i][j].right_top_nearest_dist2);
           }
       }
       return {intervalLength_mu, intervalLength_sigma};
   }
//
// //

pair<double,double> initializeT_divide(std::vector<vector<Cell>>& T, std::pair<vector<double>,vector<double>>& Q, double L_mu_Q, double U_mu_Q, double L_sigma_Q, double U_sigma_Q, int lambda, int n) {
    int numIntervals = 2 * lambda * static_cast<int>(std::sqrt(n));
    int numIntervals_mu ;
    int numIntervals_sigma ;
    double range_mu = U_mu_Q - L_mu_Q;
    double range_sigma = U_sigma_Q - L_sigma_Q;
    if (range_mu == 0 && range_sigma == 0) {
        numIntervals_mu = 1;
        numIntervals_sigma = 1;
    }
    else if (range_mu > range_sigma) {
        numIntervals_mu = static_cast<int>(numIntervals * (range_mu / (range_mu + range_sigma)));
        numIntervals_sigma = numIntervals - numIntervals_mu;
    } else {
        numIntervals_sigma = static_cast<int>(numIntervals * (range_sigma / (range_mu + range_sigma)));
        numIntervals_mu = numIntervals - numIntervals_sigma;
    }
    double intervalLength_mu = (U_mu_Q - L_mu_Q) / numIntervals_mu;
    double intervalLength_sigma = (U_sigma_Q - L_sigma_Q) / numIntervals_sigma;
    // cout<<range_mu<<" "<<range_sigma<<""<<numIntervals_sigma<<" "<<numIntervals_mu<<" "<<intervalLength_sigma<<" "<<intervalLength_mu<<endl;
    T.resize(numIntervals_sigma);
    for (int i = 0; i < numIntervals_sigma; ++i) {
        double sigmaBottom = L_sigma_Q + i * intervalLength_sigma;
        double sigmaTop = sigmaBottom + intervalLength_sigma;
        for (int j = 0; j < numIntervals_mu; ++j) {
            double muLeft = L_mu_Q + j * intervalLength_mu;
            double muRight = muLeft + intervalLength_mu;
            T[i].emplace_back(muLeft, muRight, sigmaBottom, sigmaTop);
        }
    }
      T[0][0].numIntervals_sigma = numIntervals_sigma;
      T[0][0].numIntervals_mu = numIntervals_mu;
      intervalLength_mu = 1 / intervalLength_mu;
      intervalLength_sigma = 1 / intervalLength_sigma;

    for (int i = 0; i < n; ++i) {
        double mu = Q.first[i];
        double sigma = Q.second[i];
        int muIndex;
        int sigmaIndex;
        if (mu!=U_mu_Q) muIndex = static_cast<int>((mu - L_mu_Q) * intervalLength_mu);
        else muIndex = numIntervals_mu -1;
        if (sigma!=U_sigma_Q) sigmaIndex = static_cast<int>((sigma - L_sigma_Q) * intervalLength_sigma);
        else sigmaIndex = numIntervals_sigma -1;
        // if (muIndex >= 0 && muIndex < numIntervals_mu && sigmaIndex >= 0 && sigmaIndex < numIntervals_sigma) {
            T[sigmaIndex][muIndex].dataList.emplace_back(mu,sigma);
        // }
    }
    //leftbottom point and left,bottom edge
    for (int i = 0; i < numIntervals_sigma; ++i) {
        for (int j = 0; j < numIntervals_mu; ++j) {
            //继承两个点
        int row = MAX(i-1,0);
        int col = MAX(j-1,0);
        if (T[row][col].left_top_nearest_dist1!=MAX_VAL) {
            double mu = T[row][col].left_top_nearest_mu1;
            double sigma = T[row][col].left_top_nearest_sigma1;
            double dist = sqrt(dist(T[i][j].Lleft_bottom.first,mu)+dist(T[i][j].Lleft_bottom.second,sigma));
            T[i][j].left_bottom_nearest_dist1 = dist;
            T[i][j].left_bottom_nearest_mu1 = mu;
            T[i][j].left_bottom_nearest_sigma1 = sigma;
            if (i>=1&&j>=1) {
                T[i-1][j-1].right_top_nearest_dist1 = dist;
                T[i-1][j-1].right_top_nearest_mu1 = mu;
                T[i-1][j-1].right_top_nearest_sigma1 = sigma;
                T[i-1][j].left_top_nearest_dist1 = dist;
                T[i-1][j].left_top_nearest_mu1 = mu;
                T[i-1][j].left_top_nearest_sigma1 = sigma;
                T[i][j-1].right_bottom_nearest_dist1 = dist;
                T[i][j-1].right_bottom_nearest_mu1 = mu;
                T[i][j-1].right_bottom_nearest_sigma1 = sigma;
            }else if (i>=1) {
                T[i-1][j].left_top_nearest_dist1 = dist;
                T[i-1][j].left_top_nearest_mu1 = mu;
                T[i-1][j].left_top_nearest_sigma1 = sigma;
            }else if (j>=1) {
                T[i][j-1].right_bottom_nearest_dist1 = dist;
                T[i][j-1].right_bottom_nearest_mu1 = mu;
                T[i][j-1].right_bottom_nearest_sigma1 = sigma;
            }
        }
        if (T[row][col].right_bottom_nearest_dist1!=MAX_VAL) {
            double mu = T[row][col].right_bottom_nearest_mu1;
            double sigma = T[row][col].right_bottom_nearest_sigma1;
            double dist = sqrt(dist(T[i][j].Lleft_bottom.first,mu)+dist(T[i][j].Lleft_bottom.second,sigma));
            if (dist < T[i][j].left_bottom_nearest_dist1) {
                T[i][j].left_bottom_nearest_dist1 = dist;
                T[i][j].left_bottom_nearest_mu1 = mu;
                T[i][j].left_bottom_nearest_sigma1 = sigma;
                if (i>=1&&j>=1) {
                    T[i-1][j-1].right_top_nearest_dist1 = dist;
                    T[i-1][j-1].right_top_nearest_mu1 = mu;
                    T[i-1][j-1].right_top_nearest_sigma1 = sigma;
                    T[i-1][j].left_top_nearest_dist1 = dist;
                    T[i-1][j].left_top_nearest_mu1 = mu;
                    T[i-1][j].left_top_nearest_sigma1 = sigma;
                    T[i][j-1].right_bottom_nearest_dist1 = dist;
                    T[i][j-1].right_bottom_nearest_mu1 = mu;
                    T[i][j-1].right_bottom_nearest_sigma1 = sigma;
                }else if (i>=1) {
                    T[i-1][j].left_top_nearest_dist1 = dist;
                    T[i-1][j].left_top_nearest_mu1 = mu;
                    T[i-1][j].left_top_nearest_sigma1 = sigma;
                }else if (j>=1) {
                    T[i][j-1].right_bottom_nearest_dist1 = dist;
                    T[i][j-1].right_bottom_nearest_mu1 = mu;
                    T[i][j-1].right_bottom_nearest_sigma1 = sigma;
                }
            }
        }
        //遍历一个格子来给点赋值
        for (int k=0;k<T[row][col].dataList.size();k++) {
            if (i<1||j<1) break;
            double mu = T[row][col].dataList[k].first;
            double sigma = T[row][col].dataList[k].second;
            double dist = sqrt(dist(T[i][j].Lleft_bottom.first,mu)+dist(T[i][j].Lleft_bottom.second,sigma));
            if (dist < T[i][j].left_bottom_nearest_dist1) {
                T[i][j].left_bottom_nearest_dist1 = dist;
                T[i][j].left_bottom_nearest_mu1 = mu;
                T[i][j].left_bottom_nearest_sigma1 = sigma;
                if (i>=1&&j>=1) {
                    T[i-1][j-1].right_top_nearest_dist1 = dist;
                    T[i-1][j-1].right_top_nearest_mu1 = mu;
                    T[i-1][j-1].right_top_nearest_sigma1 = sigma;
                    T[i-1][j].left_top_nearest_dist1 = dist;
                    T[i-1][j].left_top_nearest_mu1 = mu;
                    T[i-1][j].left_top_nearest_sigma1 = sigma;
                    T[i][j-1].right_bottom_nearest_dist1 = dist;
                    T[i][j-1].right_bottom_nearest_mu1 = mu;
                    T[i][j-1].right_bottom_nearest_sigma1 = sigma;
                }else if (i>=1) {
                    T[i-1][j].left_top_nearest_dist1 = dist;
                    T[i-1][j].left_top_nearest_mu1 = mu;
                    T[i-1][j].left_top_nearest_sigma1 = sigma;
                }else if (j>=1) {
                    T[i][j-1].right_bottom_nearest_dist1 = dist;
                    T[i][j-1].right_bottom_nearest_mu1 = mu;
                    T[i][j-1].right_bottom_nearest_sigma1 = sigma;
                }
            }
        }
        //遍历一个格子来给left edge赋值
        row = i;
        col = MAX(0,j-1);
        if (T[row][col].dataList.size()!=0) {
            for (int k=0;k<T[row][col].dataList.size();k++) {
                if (j<1)break;
                double mu = T[row][col].dataList[k].first;
                double sigma = T[row][col].dataList[k].second;
                double dist = T[i][j].Lleft_bottom.first - mu;
                if (dist < T[i][j].minDistToLleft) {
                    T[i][j].minDistToLleft = dist;
                    T[i][j].left_nearest_mu = mu;
                    T[i][j].left_nearest_sigma = sigma;
                }
            }
        }else {
            if (T[row][col].minDistToLleft!=MAX_VAL) {
                double mu = T[row][col].left_nearest_mu;
                double sigma = T[row][col].left_nearest_sigma;
                double dist = T[i][j].Lleft_bottom.first - mu;
                T[i][j].minDistToLleft = dist;
                T[i][j].left_nearest_mu = mu;
                T[i][j].left_nearest_sigma = sigma;
            }
        }
        row = MAX(i-1,0);
        col = j;
        //遍历一个格子来给bottom edge赋值
        if (T[row][col].dataList.size()!=0) {
            for (int k=0;k<T[row][col].dataList.size();k++) {
                if (i<1)break;
                double mu = T[row][col].dataList[k].first;
                double sigma = T[row][col].dataList[k].second;
                double dist = T[i][j].Lright_bottom.second - sigma;
                if (dist < T[i][j].minDistToLbottom) {
                    T[i][j].minDistToLbottom = dist;
                    T[i][j].bottom_nearest_mu = mu;
                    T[i][j].bottom_nearest_sigma = sigma;
                }
            }
        }else {
            if (T[row][col].minDistToLbottom!=MAX_VAL) {
                double mu = T[row][col].bottom_nearest_mu;
                double sigma = T[row][col].bottom_nearest_sigma;
                double dist = T[i][j].Lright_bottom.second - sigma;
                T[i][j].minDistToLbottom = dist;
                T[i][j].bottom_nearest_mu = mu;
                T[i][j].bottom_nearest_sigma = sigma;
            }
        }
        }
    }

    //righttop point and right,top edge
    for (int i = numIntervals_sigma - 1; i >= 0; --i) {
        for (int j = numIntervals_mu - 1; j >= 0; --j) {
            //继承两个点
            int row = MIN(i+1,numIntervals_sigma-1);
            int col = MIN(j+1,numIntervals_mu - 1);
            if (T[row][col].left_top_nearest_dist2!=MAX_VAL) {
                double mu = T[row][col].left_top_nearest_mu2;
                double sigma = T[row][col].left_top_nearest_sigma2;
                double dist = sqrt(dist(T[i][j].Lright_top.first,mu)+dist(T[i][j].Lright_top.second,sigma));
                T[i][j].right_top_nearest_dist2 = dist;
                T[i][j].right_top_nearest_mu2 = mu;
                T[i][j].right_top_nearest_sigma2 = sigma;
                if (i+2<=numIntervals_sigma&&j+2<=numIntervals_mu) {
                    T[i+1][j+1].left_bottom_nearest_dist2 = dist;
                    T[i+1][j+1].left_bottom_nearest_mu2 = mu;
                    T[i+1][j+1].left_bottom_nearest_sigma2 = sigma;

                    T[i][j+1].left_top_nearest_dist2 = dist;
                    T[i][j+1].left_top_nearest_mu2 = mu;
                    T[i][j+1].left_top_nearest_sigma2 = sigma;

                    T[i+1][j].right_bottom_nearest_dist2 = dist;
                    T[i+1][j].right_bottom_nearest_mu2 = mu;
                    T[i+1][j].right_bottom_nearest_sigma2 = sigma;
                }else if (j+2<=numIntervals_mu) {
                    T[i][j+1].left_top_nearest_dist2 = dist;
                    T[i][j+1].left_top_nearest_mu2 = mu;
                    T[i][j+1].left_top_nearest_sigma2 = sigma;
                }else if (i+2<=numIntervals_sigma) {
                    T[i+1][j].right_bottom_nearest_dist2 = dist;
                    T[i+1][j].right_bottom_nearest_mu2 = mu;
                    T[i+1][j].right_bottom_nearest_sigma2 = sigma;
                }
            }
            if (T[row][col].right_bottom_nearest_dist2!=MAX_VAL) {
                double mu = T[row][col].right_bottom_nearest_mu2;
                double sigma = T[row][col].right_bottom_nearest_sigma2;
                double dist = sqrt(dist(T[i][j].Lright_top.first,mu)+dist(T[i][j].Lright_top.second,sigma));
                if (dist < T[i][j].right_top_nearest_dist2) {
                    T[i][j].right_top_nearest_dist2 = dist;
                    T[i][j].right_top_nearest_mu2 = mu;
                    T[i][j].right_top_nearest_sigma2 = sigma;
                    if (i+2<=numIntervals_sigma&&j+2<=numIntervals_mu) {
                        T[i+1][j+1].left_bottom_nearest_dist2 = dist;
                        T[i+1][j+1].left_bottom_nearest_mu2 = mu;
                        T[i+1][j+1].left_bottom_nearest_sigma2 = sigma;

                        T[i][j+1].left_top_nearest_dist2 = dist;
                        T[i][j+1].left_top_nearest_mu2 = mu;
                        T[i][j+1].left_top_nearest_sigma2 = sigma;

                        T[i+1][j].right_bottom_nearest_dist2 = dist;
                        T[i+1][j].right_bottom_nearest_mu2 = mu;
                        T[i+1][j].right_bottom_nearest_sigma2 = sigma;
                    }else if (j+2<=numIntervals_mu) {
                        T[i][j+1].left_top_nearest_dist2 = dist;
                        T[i][j+1].left_top_nearest_mu2 = mu;
                        T[i][j+1].left_top_nearest_sigma2 = sigma;
                    }else if (i+2<=numIntervals_sigma) {
                        T[i+1][j].right_bottom_nearest_dist2 = dist;
                        T[i+1][j].right_bottom_nearest_mu2 = mu;
                        T[i+1][j].right_bottom_nearest_sigma2 = sigma;
                    }
                }
            }
            //遍历一个格子来给点赋值
            for (int k=0;k<T[row][col].dataList.size();k++) {
                if (i+2>numIntervals_sigma||j+2>numIntervals_mu)break;
                double mu = T[row][col].dataList[k].first;
                double sigma = T[row][col].dataList[k].second;
                double dist = sqrt(dist(T[i][j].Lright_top.first,mu)+dist(T[i][j].Lright_top.second,sigma));
                if (dist < T[i][j].right_top_nearest_dist2) {
                    T[i][j].right_top_nearest_dist2 = dist;
                    T[i][j].right_top_nearest_mu2 = mu;
                    T[i][j].right_top_nearest_sigma2 = sigma;
                    if (i+2<=numIntervals_sigma&&j+2<=numIntervals_mu) {
                        T[i+1][j+1].left_bottom_nearest_dist2 = dist;
                        T[i+1][j+1].left_bottom_nearest_mu2 = mu;
                        T[i+1][j+1].left_bottom_nearest_sigma2 = sigma;

                        T[i][j+1].left_top_nearest_dist2 = dist;
                        T[i][j+1].left_top_nearest_mu2 = mu;
                        T[i][j+1].left_top_nearest_sigma2 = sigma;

                        T[i+1][j].right_bottom_nearest_dist2 = dist;
                        T[i+1][j].right_bottom_nearest_mu2 = mu;
                        T[i+1][j].right_bottom_nearest_sigma2 = sigma;
                    }else if (j+2<=numIntervals_mu) {
                        T[i][j+1].left_top_nearest_dist2 = dist;
                        T[i][j+1].left_top_nearest_mu2 = mu;
                        T[i][j+1].left_top_nearest_sigma2 = sigma;
                    }else if (i+2<=numIntervals_sigma) {
                        T[i+1][j].right_bottom_nearest_dist2 = dist;
                        T[i+1][j].right_bottom_nearest_mu2 = mu;
                        T[i+1][j].right_bottom_nearest_sigma2 = sigma;
                    }
                }
            }
            row = i;
            col = MIN(j+1,numIntervals_mu - 1);
            //遍历一个格子来给right edge赋值
            if (T[row][col].dataList.size()!=0) {
                for (int k=0;k<T[row][col].dataList.size();k++) {
                    if (j+2>numIntervals_mu)break;
                    double mu = T[row][col].dataList[k].first;
                    double sigma = T[row][col].dataList[k].second;
                    double dist = mu - T[i][j].Lright_bottom.first;
                    if (dist < T[i][j].minDistToLright) {
                        T[i][j].minDistToLright = dist;
                        T[i][j].right_nearest_mu = mu;
                        T[i][j].right_nearest_sigma = sigma;
                    }
                }
            }else {
                if (T[row][col].minDistToLright!=MAX_VAL) {
                    double mu = T[row][col].right_nearest_mu;
                    double sigma = T[row][col].right_nearest_sigma;
                    double dist = mu - T[i][j].Lright_bottom.first;
                    T[i][j].minDistToLright = dist;
                    T[i][j].right_nearest_mu = mu;
                    T[i][j].right_nearest_sigma = sigma;
                }
            }
            row = MIN(i+1,numIntervals_sigma-1);
            col = j;
            //遍历一个格子来给top edge赋值
            if (T[row][col].dataList.size()!=0) {
                for (int k=0;k<T[row][col].dataList.size();k++) {
                    if (i+2>numIntervals_sigma)break;
                    double mu = T[row][col].dataList[k].first;
                    double sigma = T[row][col].dataList[k].second;
                    double dist = sigma - T[i][j].Lright_top.second;
                    if (dist < T[i][j].minDistToLtop) {
                        T[i][j].minDistToLtop = dist;
                        T[i][j].top_nearest_mu = mu;
                        T[i][j].top_nearest_sigma = sigma;
                    }
                }
            }else {
                if (T[row][col].minDistToLtop!=MAX_VAL) {
                    double mu = T[row][col].top_nearest_mu;
                    double sigma = T[row][col].top_nearest_sigma;
                    double dist = sigma - T[i][j].Lright_top.second;
                    T[i][j].minDistToLtop = dist;
                    T[i][j].top_nearest_mu = mu;
                    T[i][j].top_nearest_sigma = sigma;
                }
            }
        }
    }


    //rightbottom point
    for (int i = 0; i < numIntervals_sigma; ++i) {
        for (int j = numIntervals_mu - 1; j >= 0; --j) {
            //继承两个点
            int row = MAX(i-1,0);
            int col = MIN(j+1,numIntervals_mu - 1);
            if (T[row][col].left_bottom_nearest_dist3!=MAX_VAL) {
                double mu = T[row][col].left_bottom_nearest_mu3;
                double sigma = T[row][col].left_bottom_nearest_sigma3;
                double dist = sqrt(dist(T[i][j].Lright_bottom.first,mu)+dist(T[i][j].Lright_bottom.second,sigma));
                T[i][j].right_bottom_nearest_dist3 = dist;
                T[i][j].right_bottom_nearest_mu3 = mu;
                T[i][j].right_bottom_nearest_sigma3 = sigma;
                if (i>=1&&j+2<=numIntervals_mu) {
                    T[i-1][j+1].left_top_nearest_dist3 = dist;
                    T[i-1][j+1].left_top_nearest_mu3 = mu;
                    T[i-1][j+1].left_top_nearest_sigma3 = sigma;
                    T[i-1][j].right_top_nearest_dist3 = dist;
                    T[i-1][j].right_top_nearest_mu3 = mu;
                    T[i-1][j].right_top_nearest_sigma3 = sigma;
                    T[i][j+1].left_bottom_nearest_dist3 = dist;
                    T[i][j+1].left_bottom_nearest_mu3 = mu;
                    T[i][j+1].left_bottom_nearest_sigma3 = sigma;
                }else if (i>=1) {
                    T[i-1][j].right_top_nearest_dist3 = dist;
                    T[i-1][j].right_top_nearest_mu3 = mu;
                    T[i-1][j].right_top_nearest_sigma3 = sigma;
                }else if (j+2<=numIntervals_mu) {
                    T[i][j+1].left_bottom_nearest_dist3 = dist;
                    T[i][j+1].left_bottom_nearest_mu3 = mu;
                    T[i][j+1].left_bottom_nearest_sigma3 = sigma;
                }
            }
            if (T[row][col].right_top_nearest_dist3!=MAX_VAL) {
                double mu = T[row][col].right_top_nearest_mu3;
                double sigma = T[row][col].right_top_nearest_sigma3;
                double dist = sqrt(dist(T[i][j].Lright_bottom.first,mu)+dist(T[i][j].Lright_bottom.second,sigma));
                if (dist < T[i][j].right_bottom_nearest_dist3) {
                    T[i][j].right_bottom_nearest_dist3 = dist;
                    T[i][j].right_bottom_nearest_mu3 = mu;
                    T[i][j].right_bottom_nearest_sigma3 = sigma;
                    if (i>=1&&j+2<=numIntervals_mu) {
                        T[i-1][j+1].left_top_nearest_dist3 = dist;
                        T[i-1][j+1].left_top_nearest_mu3 = mu;
                        T[i-1][j+1].left_top_nearest_sigma3 = sigma;
                        T[i-1][j].right_top_nearest_dist3 = dist;
                        T[i-1][j].right_top_nearest_mu3 = mu;
                        T[i-1][j].right_top_nearest_sigma3 = sigma;
                        T[i][j+1].left_bottom_nearest_dist3 = dist;
                        T[i][j+1].left_bottom_nearest_mu3 = mu;
                        T[i][j+1].left_bottom_nearest_sigma3 = sigma;
                    }else if (i>=1) {
                        T[i-1][j].right_top_nearest_dist3 = dist;
                        T[i-1][j].right_top_nearest_mu3 = mu;
                        T[i-1][j].right_top_nearest_sigma3 = sigma;
                    }else if (j+2<=numIntervals_mu) {
                        T[i][j+1].left_bottom_nearest_dist3 = dist;
                        T[i][j+1].left_bottom_nearest_mu3 = mu;
                        T[i][j+1].left_bottom_nearest_sigma3 = sigma;
                    }
                }
            }
            //遍历一个格子来给点赋值
            for (int k=0;k<T[row][col].dataList.size();k++) {
                if (i<1||j+2>numIntervals_mu)break;
                double mu = T[row][col].dataList[k].first;
                double sigma = T[row][col].dataList[k].second;
                double dist = sqrt(dist(T[i][j].Lright_bottom.first,mu)+dist(T[i][j].Lright_bottom.second,sigma));
                if (dist < T[i][j].right_bottom_nearest_dist3) {
                    T[i][j].right_bottom_nearest_dist3 = dist;
                    T[i][j].right_bottom_nearest_mu3 = mu;
                    T[i][j].right_bottom_nearest_sigma3 = sigma;
                    if (i>=1&&j+2<=numIntervals_mu) {
                        T[i-1][j+1].left_top_nearest_dist3 = dist;
                        T[i-1][j+1].left_top_nearest_mu3 = mu;
                        T[i-1][j+1].left_top_nearest_sigma3 = sigma;
                        T[i-1][j].right_top_nearest_dist3 = dist;
                        T[i-1][j].right_top_nearest_mu3 = mu;
                        T[i-1][j].right_top_nearest_sigma3 = sigma;
                        T[i][j+1].left_bottom_nearest_dist3 = dist;
                        T[i][j+1].left_bottom_nearest_mu3 = mu;
                        T[i][j+1].left_bottom_nearest_sigma3 = sigma;
                    }else if (i>=1) {
                        T[i-1][j].right_top_nearest_dist3 = dist;
                        T[i-1][j].right_top_nearest_mu3 = mu;
                        T[i-1][j].right_top_nearest_sigma3 = sigma;
                    }else if (j+2<=numIntervals_mu) {
                        T[i][j+1].left_bottom_nearest_dist3 = dist;
                        T[i][j+1].left_bottom_nearest_mu3 = mu;
                        T[i][j+1].left_bottom_nearest_sigma3 = sigma;
                    }
                }
            }
        }
    }

    //lefttop point
    for (int i = numIntervals_sigma - 1; i >= 0; --i) {
        for (int j = 0; j < numIntervals_mu; ++j) {
            int row = MIN(i+1,numIntervals_sigma-1);
            int col = MAX(j-1,0);
            //继承两个点
            if (T[row][col].left_bottom_nearest_dist4!=MAX_VAL) {
                double mu = T[row][col].left_bottom_nearest_mu4;
                double sigma = T[row][col].left_bottom_nearest_sigma4;
                double dist = sqrt(dist(T[i][j].Lleft_top.first,mu)+dist(T[i][j].Lleft_top.second,sigma));
                T[i][j].left_top_nearest_dist4 = dist;
                T[i][j].left_top_nearest_mu4 = mu;
                T[i][j].left_top_nearest_sigma4 = sigma;
                if (j>=1&&i+2<=numIntervals_sigma) {
                    T[i+1][j-1].right_bottom_nearest_dist4 = dist;
                    T[i+1][j-1].right_bottom_nearest_mu4 = mu;
                    T[i+1][j-1].right_bottom_nearest_sigma4 = sigma;
                }else if (j>=1) {
                    T[i][j-1].right_top_nearest_dist4 = dist;
                    T[i][j-1].right_top_nearest_mu4 = mu;
                    T[i][j-1].right_top_nearest_sigma4 = sigma;
                }else if (i+2<=numIntervals_sigma) {
                    T[i+1][j].left_bottom_nearest_dist4 = dist;
                    T[i+1][j].left_bottom_nearest_mu4 = mu;
                    T[i+1][j].left_bottom_nearest_sigma4 = sigma;
                }
            }
            if (T[row][col].right_top_nearest_dist4!=MAX_VAL) {
                double mu = T[row][col].right_top_nearest_mu4;
                double sigma = T[row][col].right_top_nearest_sigma4;
                double dist = sqrt(dist(T[i][j].Lleft_top.first,mu)+dist(T[i][j].Lleft_top.second,sigma));
                if (dist < T[i][j].left_top_nearest_dist4) {
                    T[i][j].left_top_nearest_dist4 = dist;
                    T[i][j].left_top_nearest_mu4 = mu;
                    T[i][j].left_top_nearest_sigma4 = sigma;
                    if (j>=1&&i+2<=numIntervals_sigma) {
                        T[i+1][j-1].right_bottom_nearest_dist4 = dist;
                        T[i+1][j-1].right_bottom_nearest_mu4 = mu;
                        T[i+1][j-1].right_bottom_nearest_sigma4 = sigma;
                    }else if (j>=1) {
                        T[i][j-1].right_top_nearest_dist4 = dist;
                        T[i][j-1].right_top_nearest_mu4 = mu;
                        T[i][j-1].right_top_nearest_sigma4 = sigma;
                    }else if (i+2<=numIntervals_sigma) {
                        T[i+1][j].left_bottom_nearest_dist4 = dist;
                        T[i+1][j].left_bottom_nearest_mu4 = mu;
                        T[i+1][j].left_bottom_nearest_sigma4 = sigma;
                    }
                }
            }
            //遍历一个格子来给点赋值
            for (int k=0;k<T[row][col].dataList.size();k++) {
                if (j<1||i+2>numIntervals_sigma)break;
                double mu = T[row][col].dataList[k].first;
                double sigma = T[row][col].dataList[k].second;
                double dist = sqrt(dist(T[i][j].Lleft_top.first,mu)+dist(T[i][j].Lleft_top.second,sigma));
                if (dist < T[i][j].left_top_nearest_dist4) {
                    T[i][j].left_top_nearest_dist4 = dist;
                    T[i][j].left_top_nearest_mu4 = mu;
                    T[i][j].left_top_nearest_sigma4 = sigma;
                    if (j>=1&&i+2<=numIntervals_sigma) {
                        T[i+1][j-1].right_bottom_nearest_dist4 = dist;
                        T[i+1][j-1].right_bottom_nearest_mu4 = mu;
                        T[i+1][j-1].right_bottom_nearest_sigma4 = sigma;
                    }else if (j>=1) {
                        T[i][j-1].right_top_nearest_dist4 = dist;
                        T[i][j-1].right_top_nearest_mu4 = mu;
                        T[i][j-1].right_top_nearest_sigma4 = sigma;
                    }else if (i+2<=numIntervals_sigma) {
                        T[i+1][j].left_bottom_nearest_dist4 = dist;
                        T[i+1][j].left_bottom_nearest_mu4 = mu;
                        T[i+1][j].left_bottom_nearest_sigma4 = sigma;
                    }
                }
            }
        }
    }

    //update
    for (int i = 0; i < numIntervals_sigma; ++i) {
        for (int j = 0; j < numIntervals_mu; ++j) {
            T[i][j].minDistToLleft = MIN(MIN(T[i][j].minDistToLleft,T[i][j].left_bottom_nearest_dist1),T[i][j].left_top_nearest_dist4);
            T[i][j].minDistToLright = MIN(MIN(T[i][j].minDistToLright,T[i][j].right_bottom_nearest_dist3),T[i][j].right_top_nearest_dist2);
            T[i][j].minDistToLbottom = MIN(MIN(T[i][j].minDistToLbottom,T[i][j].left_bottom_nearest_dist1),T[i][j].right_bottom_nearest_dist3);
            T[i][j].minDistToLtop = MIN(MIN(T[i][j].minDistToLtop,T[i][j].left_top_nearest_dist4),T[i][j].right_top_nearest_dist2);
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
     int col;
     int row;
     if (mu_X!=U_mu_Q) col = static_cast<int>((mu_X - L_mu_Q) * intervalLength_mu);
     else col = T[0][0].numIntervals_mu -1;
     if (sigma_X!=U_sigma_Q) row = static_cast<int>((sigma_X - L_sigma_Q) * intervalLength_sigma);
     else row = T[0][0].numIntervals_sigma -1;
    if (T[row][col].dataList.size() >= 2) {
        return dist;
    } else {
        double leftDist = (T[row][col].minDistToLleft+mu_X-T[row][col].Lleft_bottom.first)*(T[row][col].minDistToLleft+mu_X-T[row][col].Lleft_bottom.first);
        double rightDist = (T[row][col].minDistToLright+T[row][col].Lright_bottom.first-mu_X)*(T[row][col].minDistToLright+T[row][col].Lright_bottom.first-mu_X);
        double bottomDist = (T[row][col].minDistToLbottom+sigma_X-T[row][col].Lleft_bottom.second)*(T[row][col].minDistToLbottom+sigma_X-T[row][col].Lleft_bottom.second);
        double topDist = (T[row][col].minDistToLtop+T[row][col].Lleft_top.second-sigma_X)*(T[row][col].minDistToLtop+T[row][col].Lleft_top.second-sigma_X);
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


vector<vector<double>> Preprocessing_MC3E_DIVIDE(vector<vector<double>>& Q,
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
        intervalLength = initializeT_divide(T[i], Q_part, L_mu_Q[i], U_mu_Q[i], L_sigma_Q[i],U_sigma_Q[i],lambda, num[i]);
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
    return ((Ai[0] <= L_mu_Q || Ai[0] >= U_mu_Q) && (L_sigma_Q <= Ai[1] && Ai[1] <= U_sigma_Q)) ||
           ((L_mu_Q <= Ai[0] && Ai[0] <= U_mu_Q) && (Ai[1] <= L_sigma_Q || Ai[1] >= U_sigma_Q));
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
            if (Ai[0] <= L_mu_Q[i] && L_sigma_Q[i] <= Ai[1] && Ai[1] <= U_sigma_Q[i]) {
                V_plus = {L_mu_Q[i], L_sigma_Q[i]};
                V_minus = {L_mu_Q[i], U_sigma_Q[i]};
                d_plus = Dbl[i];
                d_minus = Dtl[i];
                MiddleLeftOrRight = true;
            } else if (Ai[0] >= U_mu_Q[i] && L_sigma_Q[i] <= Ai[1] && Ai[1] <= U_sigma_Q[i]) {
                V_plus = {U_mu_Q[i], L_sigma_Q[i]};
                V_minus = {U_mu_Q[i], U_sigma_Q[i]};
                d_plus = Dbr[i];
                d_minus = Dtr[i];
                MiddleLeftOrRight = true;
            } else if (L_mu_Q[i] <= Ai[0] && Ai[0] <= U_mu_Q[i] && Ai[1] <= L_sigma_Q[i]) {
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
