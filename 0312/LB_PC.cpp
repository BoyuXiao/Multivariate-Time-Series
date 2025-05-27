#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>

// 计算两点之间的欧几里得距离
double euclideanDistance(const std::vector<double>& p1, const std::vector<double>& p2) {
    double sum = 0.0;
    for (size_t i = 0; i < p1.size(); ++i) {
        double diff = p1[i] - p2[i];
        sum += diff * diff;
    }
    // return std::sqrt(sum);
    return sum;
}

// 查找一个窗口内的边界框
std::vector<std::pair<std::vector<double>, std::vector<double>>> findBoundingBoxes(
    const std::vector<std::vector<double>>& window, int K, int Q) {
    std::vector<std::vector<std::vector<double>>> cellMembers;
    std::vector<std::pair<std::vector<double>, std::vector<double>>> boundingBoxes;
    size_t dim = window[0].size();

    // 找到每个维度的最小值和最大值
    std::vector<double> minValues(dim, std::numeric_limits<double>::max());
    std::vector<double> maxValues(dim, std::numeric_limits<double>::lowest());
    for (const auto& point : window) {
        for (size_t i = 0; i < dim; ++i) {
            minValues[i] = std::min(minValues[i], point[i]);
            maxValues[i] = std::max(maxValues[i], point[i]);
        }
    }

    // 计算每个维度的单元格数量和长度
    std::vector<int> numCells(dim);
    std::vector<double> cellLengths(dim);
    for (size_t i = 0; i < dim; ++i) {
        numCells[i] = 1 + static_cast<int>((maxValues[i] - minValues[i]) * Q);
        cellLengths[i] = (maxValues[i] - minValues[i]) / numCells[i] + 1e-8;
    }

    // 将点分配到单元格中
    for (const auto& point : window) {
        std::vector<int> cellIndices(dim);
        for (size_t i = 0; i < dim; ++i) {
            cellIndices[i] = static_cast<int>((point[i] - minValues[i]) / cellLengths[i]);
        }
        size_t cellId = 0;
        for (size_t i = 0; i < dim; ++i) {
            cellId = cellId * numCells[i] + cellIndices[i];
        }
        if (cellMembers.size() <= cellId) {
            cellMembers.resize(cellId + 1);
        }
        cellMembers[cellId].push_back(point);
    }

    // 为每个非空单元格创建边界框
    for (const auto& group : cellMembers) {
        if (!group.empty()) {
            std::vector<double> minBounds(dim, std::numeric_limits<double>::max());
            std::vector<double> maxBounds(dim, std::numeric_limits<double>::lowest());
            for (const auto& point : group) {
                for (size_t i = 0; i < dim; ++i) {
                    minBounds[i] = std::min(minBounds[i], point[i]);
                    maxBounds[i] = std::max(maxBounds[i], point[i]);
                }
            }
            boundingBoxes.emplace_back(minBounds, maxBounds);
        }
    }

    // 如果边界框数量超过K，合并多余的边界框
    if (static_cast<int>(boundingBoxes.size()) > K) {
        std::vector<double> combinedMin(dim, std::numeric_limits<double>::max());
        std::vector<double> combinedMax(dim, std::numeric_limits<double>::lowest());
        for (int i = K - 1; i < static_cast<int>(boundingBoxes.size()); ++i) {
            for (size_t j = 0; j < dim; ++j) {
                combinedMin[j] = std::min(combinedMin[j], boundingBoxes[i].first[j]);
                combinedMax[j] = std::max(combinedMax[j], boundingBoxes[i].second[j]);
            }
        }
        boundingBoxes.resize(K - 1);
        boundingBoxes.emplace_back(combinedMin, combinedMax);
    }

    return boundingBoxes;
}

// 计算单个查询序列与参考序列之间的下界
double calculateLBPCLowerBound(
    const std::vector<std::vector<double>>& query,
    const std::vector<std::vector<double>>& reference,
    int windowSize, int K, int Q) {
    size_t queryLength = query.size();
    size_t referenceLength = reference.size();
    int w = std::max(windowSize, static_cast<int>(std::abs(static_cast<int>(queryLength) - static_cast<int>(referenceLength))));

    double lowerBound = 0.0;
    for (size_t i = 0; i < queryLength; ++i) {
        int left = std::max(0, static_cast<int>(i) - w);
        int right = std::min(static_cast<int>(referenceLength), static_cast<int>(i) + w);
        std::vector<std::vector<double>> window;
        for (int j = left; j < right; ++j) {
            window.push_back(reference[j]);
        }
        auto boundingBoxes = findBoundingBoxes(window, K, Q);

        double minDist = std::numeric_limits<double>::max();
        for (const auto& box : boundingBoxes) {
            double dist = 0.0;
            for (size_t d = 0; d < query[i].size(); ++d) {
                if (query[i][d] < box.first[d]) {
                    dist += std::pow(box.first[d] - query[i][d], 2);
                } else if (query[i][d] > box.second[d]) {
                    dist += std::pow(query[i][d] - box.second[d], 2);
                }
            }
            // dist = std::sqrt(dist);
            minDist = std::min(minDist, dist);
        }
        lowerBound += minDist;
    }
    return sqrt(lowerBound);
}