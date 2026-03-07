#include <iostream>
#include <fstream>
#include <vector>
#include "funcs.h"
#include "def.h"
using namespace std;


int main(int argc, char* argv[]) {
    vector<string> testFilePaths = {data1_test, data2_test, data3_test, data4_test, data5_test, data6_test, data7_test, data8_test,
        data9_test, data10_test, data11_test, data12_test, data13_test, data14_test, data15_test,
        data16_test, data17_test, data18_test, data19_test, data20_test, data21_test, data22_test,
        data23_test, data24_test, data25_test, data26_test
    };

    vector<string> trainFilePaths = {data1_train, data2_train, data3_train, data4_train, data5_train, data6_train, data7_train, data8_train,
        data9_train, data10_train, data11_train, data12_train, data13_train, data14_train, data15_train,
        data16_train, data17_train, data18_train, data19_train, data20_train, data21_train, data22_train,
        data23_train, data24_train, data25_train, data26_train
    };

    vector<std::string> dataNames = {
        data1_name, data2_name, data3_name, data4_name, data5_name, data6_name, data7_name, data8_name,
        data9_name, data10_name, data11_name, data12_name, data13_name, data14_name, data15_name,
        data16_name, data17_name, data18_name, data19_name, data20_name, data21_name, data22_name,
        data23_name, data24_name, data25_name, data26_name
    };

    if (argc != 3) {
        cerr << "Usage: " << argv[0] << " start_index end_index" << endl;
        return 1;
    }

    int start_index = atoi(argv[1]);
    int end_index = atoi(argv[2]);

    if (start_index < 0 || start_index >= testFilePaths.size() || end_index < 0 || end_index >= testFilePaths.size() || start_index > end_index) {
        cerr << "Invalid start or end index." << endl;
        return 1;
    }

    for (int i = start_index; i <= end_index; i++) {
        if (i==5||i==-4||i==16||i==9||i==18||i==20)continue;
        string filePath_test = testFilePaths[i];
        string filePath_train = trainFilePaths[i];
        string file_name = dataNames[i];
        vector<vector<vector<double>>> data_test = readCSV(filePath_test);
        vector<vector<vector<double>>> data_train = readCSV(filePath_train);
        std::ofstream outputFile;
        outputFile.open("..//output.txt",std::ios::app);
        if (outputFile.is_open()) {
            outputFile << "Name\tTime_DTW\tTLB_MV\tTime_MV\tpruning rate_MV\tTLB_MC1\tTime_MC1\tpruning rate_MC1\tTLB_MCNaive1\tTime_MCNaive1\tpruning rate_MCNaive1\tTLB_MC2\tTime_MC2\tpruning rate_MC2\tTLB_MCNaive2\tTime_MCNaive2\tpruning rate_MCNaive2\tTLB_MC3\tTime_MC3\tpruning rate_MC3\tTLB_MC3E\tTime_MC3E\tpruning rate_MC3E\tTLB_MC3_PLUS\tTime_MC3_PLUS\tpruning rate_MC3_PLUS\tTLB_MC3_PLUS_SORT\tTime_MC3_PLUS_SORT\tpruning rate_MC3_PLUS_SORT\tTLB_MC3_PLUS_COV\tTime_MC3_PLUS_COV\tpruning rate_MC3_PLUS_COV\tTLB_MC3_PLUS_SIGMA\tTime_MC3_PLUS_SIGMA\tpruning rate_MC3_PLUS_SIGMA\tTLB_PCA\tTime_PCA\tpruning rate_PCA\tTLB_PCA4\tTime_PCA4\tpruning rate_PCA4\n";
        } else {
            std::cerr << "Failed to open the output file." << std::endl;
        }
        outputFile.close();

        TLB(data_test,data_train,file_name);
    }
    return 0;
}