#include <vector>
#include "funcs.h"
using namespace std;

double LB_MC_NAIVE1(vector<vector<double>> &Q, vector<vector<double>> &A, vector<double> &mu_Q,
                    int K, int m, int r)
{
    vector<double> mu_A(m);
    only_mu(A, mu_A);

    double result = 0.0;
    for (int i = 0; i < m; ++i)
    {
        double min_val = MAX_VAL;
        for (int j = MAX(0, i - r); j <= MIN(m - 1, i + r); ++j)
        {
            min_val = MIN(min_val, (mu_A[i] - mu_Q[j]) * (mu_A[i] - mu_Q[j]));
        }
        result += min_val;
    }
    return sqrt(static_cast<double>(K) * result);
}
