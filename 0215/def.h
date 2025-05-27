
#ifndef SCALE_DTW_DEF_H
#define SCALE_DTW_DEF_H
#include <stdlib.h>


#define data1_train "..//data//UEA//train//ArticularyWordRecognition_TRAIN_out.csv"
#define data1_test "..//data//UEA//test//ArticularyWordRecognition_TEST_out.csv"
#define data1_name "ArticularyWordRecognition"


#define data2_train "..//data//UEA//train//AtrialFibrillation_TRAIN_out.csv"
#define data2_test "..//data//UEA//test//AtrialFibrillation_TEST_out.csv"
#define data2_name "AtrialFibrillation"


#define data3_train "..//data//UEA//train//BasicMotions_TRAIN_out.csv"
#define data3_test "..//data//UEA//test//BasicMotions_TEST_out.csv"
#define data3_name "BasicMotions"


#define data4_train "..//data//UEA//train//Cricket_TRAIN_out.csv"
#define data4_test "..//data//UEA//test//Cricket_TEST_out.csv"
#define data4_name "Cricket"


#define data5_train "..//data//UEA//train//DuckDuckGeese_TRAIN_out.csv"
#define data5_test "..//data//UEA//test//DuckDuckGeese_TEST_out.csv"
#define data5_name "DuckDuckGeese"


#define data6_train "..//data//UEA//train//EigenWorms_TRAIN_out.csv"
#define data6_test "..//data//UEA//test//EigenWorms_TEST_out.csv"
#define data6_name "EigenWorms"


#define data7_train "..//data//UEA//train//Epilepsy_TRAIN_out.csv"
#define data7_test "..//data//UEA//test//Epilepsy_TEST_out.csv"
#define data7_name "Epilepsy"


#define data8_train "..//data//UEA//train//ERing_TRAIN_out.csv"
#define data8_test "..//data//UEA//test//ERing_TEST_out.csv"
#define data8_name "ERing"


#define data9_train "..//data//UEA//train//EthanolConcentration_TRAIN_out.csv"
#define data9_test "..//data//UEA//test//EthanolConcentration_TEST_out.csv"
#define data9_name "EthanolConcentration"


#define data10_train "..//data//UEA//train//FaceDetection_TRAIN_out.csv"
#define data10_test "..//data//UEA//test//FaceDetection_TEST_out.csv"
#define data10_name "FaceDetection"


#define data11_train "..//data//UEA//train//FingerMovements_TRAIN_out.csv"
#define data11_test "..//data//UEA//test//FingerMovements_TEST_out.csv"
#define data11_name "FingerMovements"


#define data12_train "..//data//UEA//train//HandMovementDirection_TRAIN_out.csv"
#define data12_test "..//data//UEA//test//HandMovementDirection_TEST_out.csv"
#define data12_name "HandMovementDirection"


#define data13_train "..//data//UEA//train//Handwriting_TRAIN_out.csv"
#define data13_test "..//data//UEA//test//Handwriting_TEST_out.csv"
#define data13_name "Handwriting"


#define data14_train "..//data//UEA//train//Heartbeat_TRAIN_out.csv"
#define data14_test "..//data//UEA//test//Heartbeat_TEST_out.csv"
#define data14_name "Heartbeat"


#define data15_train "..//data//UEA//train//Libras_TRAIN_out.csv"
#define data15_test "..//data//UEA//test//Libras_TEST_out.csv"
#define data15_name "Libras"


#define data16_train "..//data//UEA//train//LSST_TRAIN_out.csv"
#define data16_test "..//data//UEA//test//LSST_TEST_out.csv"
#define data16_name "LSST"


#define data17_train "..//data//UEA//train//MotorImagery_TRAIN_out.csv"
#define data17_test "..//data//UEA//test//MotorImagery_TEST_out.csv"
#define data17_name "MotorImagery"


#define data18_train "..//data//UEA//train//NATOPS_TRAIN_out.csv"
#define data18_test "..//data//UEA//test//NATOPS_TEST_out.csv"
#define data18_name "NATOPS"


#define data19_train "..//data//UEA//train//PEMS-SF_TRAIN_out.csv"
#define data19_test "..//data//UEA//test//PEMS-SF_TEST_out.csv"
#define data19_name "PEMS-SF"


#define data20_train "..//data//UEA//train//PenDigits_TRAIN_out.csv"
#define data20_test "..//data//UEA//test//PenDigits_TEST_out.csv"
#define data20_name "PenDigits"


#define data21_train "..//data//UEA//train//PhonemeSpectra_TRAIN_out.csv"
#define data21_test "..//data//UEA//test//PhonemeSpectra_TEST_out.csv"
#define data21_name "PhonemeSpectra"


#define data22_train "..//data//UEA//train//RacketSports_TRAIN_out.csv"
#define data22_test "..//data//UEA//test//RacketSports_TEST_out.csv"
#define data22_name "RacketSports"


#define data23_train "..//data//UEA//train//SelfRegulationSCP1_TRAIN_out.csv"
#define data23_test "..//data//UEA//test//SelfRegulationSCP1_TEST_out.csv"
#define data23_name "SelfRegulationSCP1"


#define data24_train "..//data//UEA//train//SelfRegulationSCP2_TRAIN_out.csv"
#define data24_test "..//data//UEA//test//SelfRegulationSCP2_TEST_out.csv"
#define data24_name "SelfRegulationSCP2"


#define data25_train "..//data//UEA//train//StandWalkJump_TRAIN_out.csv"
#define data25_test "..//data//UEA//test//StandWalkJump_TEST_out.csv"
#define data25_name "StandWalkJump"


#define data26_train "..//data//UEA//train//UWaveGestureLibrary_TRAIN_out.csv"
#define data26_test "..//data//UEA//test//UWaveGestureLibrary_TEST_out.csv"
#define data26_name "UWaveGestureLibrary"


using namespace std;

#define MV 0
#define MC1 0
#define MC2 0
#define MC3 0
#define MC3E 0
#define DIVIDE 1
#define PLUS_HALF 0
#define PLUS_SORT 0
#define PLUS_COV 1
#define PCA 0
#define NAIVE1 0
#define NAIVE2 0


#define dist(x,y) ((x-y)*(x-y))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))

#define INF 1e9
#define MAX_VAL 1008610086


#endif 
