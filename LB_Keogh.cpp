#include "funcs.h"
#include "lb_early.h"
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

double LB_Keogh_2_early(vector<double>& A, vector<double>& Q, int m, int r, double best_dist) {
    const double thresh_sq = lb_threshold_sq(best_dist);
    double lb = 0.0;
    vector<double> lowerBound(m);
    vector<double> upperBound(m);
    lower_upper_lemire(Q, m, r, lowerBound, upperBound);
    for (int i = 0; i < m; ++i) {
        if (A[i] < lowerBound[i]) {
            lb += (A[i] - lowerBound[i]) * (A[i] - lowerBound[i]);
        } else if (A[i] > upperBound[i]) {
            lb += (A[i] - upperBound[i]) * (A[i] - upperBound[i]);
        }
        if (lb_can_early_stop(best_dist) && lb >= thresh_sq) {
            return lb;
        }
    }
    return lb;
}

double LB_Keogh_2_env_early(vector<double>& A,
                          const vector<double>& q_lower,
                          const vector<double>& q_upper,
                          int m,
                          double best_dist) {
    const double thresh_sq = lb_threshold_sq(best_dist);
    double lb = 0.0;
    for (int i = 0; i < m; ++i) {
        if (A[i] < q_lower[i]) {
            lb += (A[i] - q_lower[i]) * (A[i] - q_lower[i]);
        } else if (A[i] > q_upper[i]) {
            lb += (A[i] - q_upper[i]) * (A[i] - q_upper[i]);
        }
        if (lb_can_early_stop(best_dist) && lb >= thresh_sq) {
            return lb;
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









