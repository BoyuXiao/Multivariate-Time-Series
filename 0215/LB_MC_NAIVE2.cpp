#include <vector>
#include "funcs.h"
#include<fstream>
using namespace std;

double LB_MC_NAIVE2(vector<vector<double>>& Q, vector<vector<double>>& A,vector<double>& mu_Q, vector<double>& sigma_Q, int K, int m, int r) {
    vector<double> mu_A(m);
    vector<double> sigma_A(m);
    mu_sigma(A,K,m,mu_A,sigma_A);
    // std::ofstream outputFile_data;
    // outputFile_data.open("..//data.txt",std::ios::app);
    double result = 0.0;
    for(int i=0;i<m;++i) {
        double min_val = MAX_VAL;
        for (int j = MAX(0, i - r); j <= MIN(m - 1, i + r); ++j) {
            min_val = MIN(min_val, K*((mu_A[i]-mu_Q[j])*(mu_A[i]-mu_Q[j])+(sigma_A[i]-sigma_Q[j])*(sigma_A[i]-sigma_Q[j])));
        }
        result+=min_val;
        // outputFile_data <<i<<" "<< min_val/K << " " << "naive" << endl;
    }
    // outputFile_data.close();
    return sqrt(result);
}







