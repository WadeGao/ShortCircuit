/*
 * @Author: your name
 * @Date: 2021-04-08 19:18:41
 * @LastEditTime: 2021-04-09 10:29:20
 * @LastEditors: Please set LastEditors
 * @Description: In User Settings Edit
 * @FilePath: /VSCode/Common.h
 */
#pragma once
#include "Equipment.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cmath>
#include <list>
#include <map>
#include <omp.h>
#include <set>
#include <thread>
#include <tuple>
#include <utility>
#include <vector>

constexpr auto PI = 3.14159265358979323846;
constexpr auto epsilon = 1e-6;
constexpr auto Imaginer = cf{0, 1};
const cf alpha = {cosf(PI * 2 / 3), sinf(PI * 2 / 3)};
const auto maxThreadsNum = std::thread::hardware_concurrency();
const Eigen::MatrixXcf St = (Eigen::MatrixXcf(3, 3) << 1, 1, 1, alpha * alpha, alpha, 1, alpha, alpha * alpha, 1).finished(); //�������(1 / 3 * S)����

using ThreeSequenceData = struct ThreeSequenceData
{
    Eigen::MatrixXf lineData1;
    Eigen::MatrixXf lineData2;
    Eigen::MatrixXf lineData0;
};

using TripVecType = std::vector<Eigen::Triplet<cf>>;
