#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

#include "funcs.h"

using namespace std;

vector<vector<vector<double>>> readCSV(const string& filePath) {
    vector<vector<vector<double>>> data;
    ifstream file(filePath);
    cout << "Reading CSV file:" << filePath << endl;
    if (!file.is_open()) {
        cerr << "open failed" << endl;
        return data;
    }

    string line;
    vector<vector<double>> currentSequence;
    bool emptyLineFound = false;  // 标记是否找到空行
    while (getline(file, line)) {
        if (line.empty()) {
            if (!currentSequence.empty()) {
                data.push_back(currentSequence);
                currentSequence.clear();
            }
            emptyLineFound = true;
            continue;
        }
        emptyLineFound = false;
        stringstream ss(line);
        string value;
        vector<double> row;
        while (getline(ss, value, ',')) {
            try {
                row.push_back(stod(value));
            }
            catch (const invalid_argument& e) {
                cerr << "failed to convert " << value << endl;
            }
        }
        if (!emptyLineFound) {
            currentSequence.push_back(row);
        }
    }
    // 处理文件结束时，若还有未添加的序列
    if (!currentSequence.empty()) {
        data.push_back(currentSequence);
    }
    file.close();

    return data;
}


//
// #include <iostream>
// #include <fstream>
// #include <sstream>
// #include <vector>
// #include <string>
//
// #include "funcs.h"
//
// using namespace std;
//
// vector<vector<vector<double>>> readCSV(const string& filePath) {
//     vector<vector<vector<double>>> data;
//     ifstream file(filePath);
//     cout << "Reading CSV file:" << filePath << endl;
//     if (!file.is_open()) {
//         cerr << "open failed" << endl;
//         return data;
//     }
//
//     string line;
//     vector<vector<double>> currentSequence;
//     bool emptyLineFound = false;  // 标记是否找到空行
//     int lineCount = 0;  // 记录当前组的行数
//
//     while (getline(file, line)) {
//         if (line.empty()) {
//             if (!currentSequence.empty()) {
//                 data.push_back(currentSequence);
//                 currentSequence.clear();
//             }
//             emptyLineFound = true;
//             lineCount = 0;  // 遇到空行，重置行数计数器
//             continue;
//         }
//         emptyLineFound = false;
//
//         if (lineCount < 60) {
//             stringstream ss(line);
//             string value;
//             vector<double> row;
//             while (getline(ss, value, ',')) {
//                 try {
//                     row.push_back(stod(value));
//                 }
//                 catch (const invalid_argument& e) {
//                     cerr << "failed to convert " << value << endl;
//                 }
//             }
//             currentSequence.push_back(row);
//             lineCount++;
//         }
//     }
//
//     // 处理文件结束时，若还有未添加的序列
//     if (!currentSequence.empty()) {
//         data.push_back(currentSequence);
//     }
//     file.close();
//
//     return data;
// }