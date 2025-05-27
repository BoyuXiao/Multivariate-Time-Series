#include <iostream>
#include <vector>
#include "funcs.h"

vector<vector<double>> transpose(const std::vector<std::vector<double>>& matrix) {
    if (matrix.empty()) {
        return {};
    }
    vector<vector<double>> result(matrix[0].size(), vector<double>(matrix.size()));
    for (size_t i = 0; i < matrix.size(); ++i) {
        for (size_t j = 0; j < matrix[0].size(); ++j) {
            result[j][i] = matrix[i][j];
        }
    }
    return result;
}

// int main() {}

