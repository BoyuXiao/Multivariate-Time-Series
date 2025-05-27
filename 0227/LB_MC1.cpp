#include <vector>
#include "funcs.h"
using namespace std;
//TODO time
double LB_MC1(vector<vector<double>>& Q, vector<vector<double>>& A,vector<double>& mu_Q, int K, int m, int r) {
    vector<double> mu_A(m);

    // only_mu(A, K, m, mu_A);
    if (K==2) {
        split_two_dim(A,mu_A,mu_A);
    }else {
        mu_sigma(A,K,m,mu_A,mu_A);
    }
    if (K==2) {
        double result = LB_Keogh_2(mu_A, mu_Q, m, r);
        return sqrt(result);
    }
    double result = K * LB_Keogh_2(mu_A, mu_Q, m, r);
    return sqrt(result);
}







