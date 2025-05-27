#include "funcs.h"
#include <iostream>
#include <limits>
using namespace std;
double LB_Keogh_2(vector<double>& A, vector<double>& Q, int m, int r) {
    double lb = 0.0;
    vector<double> lowerBound(m);
    vector<double> upperBound(m);
    lower_upper_lemire(Q,m,r,lowerBound,upperBound);
    for (int i = 0; i < m; ++i) {
        if (A[i] < lowerBound[i]) {
            lb += (A[i] - lowerBound[i]) * (A[i] - lowerBound[i]);
        } else if (A[i] > upperBound[i]) {
            lb += (A[i] - upperBound[i]) * (A[i] - upperBound[i]);
        }
    }
    return lb;
}


double LB_Keogh_old(vector<double>& Q, vector<double>& A, int m, int r) {
    double lb = 0.0;
    for (int i = 0; i < m; ++i) {
        double max_val = -numeric_limits<double>::infinity();
        double min_val = numeric_limits<double>::infinity();
        for (int j = MAX(0, i - r); j <= MIN(m - 1, i + r); ++j) {
            max_val = MAX(max_val, Q[j]);
            min_val = MIN(min_val, Q[j]);
        }
        double lowerBound = min_val;
        double upperBound = max_val;
        if (A[i] < lowerBound) {
            lb += (A[i] - lowerBound) * (A[i] - lowerBound);
        } else if (A[i] > upperBound) {
            lb += (A[i] - upperBound) * (A[i] - upperBound);
        }
    }
    return lb;
}









