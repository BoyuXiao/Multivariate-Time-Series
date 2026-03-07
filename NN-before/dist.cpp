#include "funcs.h"

double euclideanDistance_2(vector<double>& point1, vector<double>& point2) {
    double distance = 0.0;
    int K = point1.size();
    for (int i = 0; i < K; ++i) {
        double diff = point1[i] - point2[i];
        distance += diff * diff;
    }
    return distance;
}


