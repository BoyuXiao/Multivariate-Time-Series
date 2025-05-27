#include <vector>
#include "funcs.h"
using namespace std;

double LB_MC2(vector<vector<double>>& Q, vector<vector<double>>& A,vector<double>& mu_Q, vector<double>& sigma_Q, int K, int m, int r) {
    vector<double> mu_A(m);
    vector<double> sigma_A(m);
    // mu_sigma(A, K, m, mu_A, sigma_A);
    if (K==2) {
        split_two_dim(A,mu_A,sigma_A);
    }else {
        mu_sigma(A,mu_A,sigma_A);
    }
    if (K==2) {
        double result = LB_Keogh_2(mu_A, mu_Q, m, r) + K * LB_Keogh_2(sigma_A, sigma_Q, m, r);
        return sqrt(result);
    }
    double result = K * LB_Keogh_2(mu_A, mu_Q, m, r) + K * LB_Keogh_2(sigma_A, sigma_Q, m, r);
    return sqrt(result);
}







