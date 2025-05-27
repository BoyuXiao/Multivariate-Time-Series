#include <vector>
#include "funcs.h"
using namespace std;

double LB_MC1(vector<vector<double>>& Q, vector<vector<double>>& A,vector<double>& mu_Q, int K, int m, int r) {
    vector<double> mu_A(m);

    only_mu(A, K, m, mu_A);

    double result = K * LB_Keogh_2(mu_A, mu_Q, m, r);
    return sqrt(result);
}







