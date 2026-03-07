#ifndef BUCKET_H
#define BUCKET_H
#include <algorithm>
// 定义 Bucket 结构体
struct Bucket {
    double left;  // Bucket 的左边界
    double right; // Bucket 的右边界
    std::vector<double> dataList; // 用于有序存储插入到该 Bucket 中的数据
    double leftNearestDistance; // 左边界左侧最近点到左边界的距离
    double rightNearestDistance; // 右边界右侧最近点到右边界的距离

    // 构造函数，初始化成员变量
    Bucket(double l = 0, double r = 0) : left(l), right(r),
                                         leftNearestDistance(MAX_VAL),
                                         rightNearestDistance(MAX_VAL){}

    double min_val() {
        auto minIt = std::min_element(dataList.begin(), dataList.end());
        return *minIt;
    }
    double max_val() {
        auto max_it = std::max_element(dataList.begin(), dataList.end());
        return *max_it;
    }
};
#endif