cmake_minimum_required(VERSION 3.30)
project(TLB)

set(CMAKE_CXX_STANDARD 20)

add_executable(TLB main.cpp
        dtw.cpp
        dist.cpp
        LB_MV.cpp
        LB_MC1.cpp
        LB_Keogh.cpp
        LB_MC2.cpp
        LB_MC_NAIVE1.cpp
        LB_MC_NAIVE2.cpp
        mu_sigma.cpp
        TLB.cpp
        read_data.cpp
        lower_upper.cpp
        only_mu.cpp
        LB_MC3.cpp
        transpose.cpp
        bucket.h
        LB_MC3E.cpp
        cell.h
        LB_MC3_PLUS.cpp
)
