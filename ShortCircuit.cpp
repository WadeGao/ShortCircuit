/*
 * @Author: your name
 * @Date: 2021-04-08 19:54:54
 * @LastEditTime: 2021-04-09 13:06:19
 * @LastEditors: Please set LastEditors
 * @Description: In User Settings Edit
 * @FilePath: /VSCode/ShortCircuit.cpp
 */
#define EIGEN_USE_MKL_ALL

#include "Common.h"
#include "Grid.h"
#include <iostream>

int main()
{
    Eigen::MatrixXf mat(3, 5);

    mat << 2, 3, 0.024, 0.065, 0.032,
        2, 4, 0.030, 0.080, 0.040,
        4, 3, 0.018, 0.050, 0.026;

    Grid my_grid(5, {mat, mat, mat}, {{{1, 2, 1.00}, cf{0, 0.105}}, {{5, 4, 1.00}, cf{0, 0.184}}}, {}, {{1, 120, cf{0, 0.15}, 120}, {5, 120, cf{0, 0.22}, 120}});

    const auto ret = my_grid.SymmetricShortCircuit(3, cf(0.0f, 0.0f), 1);
    std::cout << "各点电压分布: " << std::endl << std::get<0>(ret) << std::endl << std::endl;

    std::cout << "各点电流分布: " << std::endl;
    const auto &I_list = std::get<1>(ret);
    for(const auto &iter : I_list)
    {
        std::cout << "I" << iter.first.first << iter.first.second << ": " << std::abs(iter.second) << " A" << std::endl;
    }

    return 0;
}

/*
P121:6-1

Eigen::MatrixXf mat(3, 5);

    mat << 2, 3, 0.024, 0.065, 0.032,
        2, 4, 0.030, 0.080, 0.040,
        4, 3, 0.018, 0.050, 0.026;

    Grid my_grid(5, {mat, mat, mat}, {{{1, 2, 1.05}, cf{0, 0.105}}, {{5, 4, 0.96}, cf{0, 0.184}}}, {}, {{1, 120, cf{0, 0.15}, 120}, {5, 120, cf{0, 0.22}, 120}});
*/

/*
P84:4-1

mat <<
    3, 4, 0.000, 0.43554, 0.018515,
    3, 5, 0.000, 0.29036, 0.012343,
    4, 5, 0.000, 0.25407, 0.010800;

    Transformer2 tran1(1, 3, 1, 120, 10.5, 120), tran2(2, 4, 1, 60, 10.5, 120);
    Generator gene1(1, 120, 0.23, 120), gene2(2, 60, 0.14, 120);

    Grid my_grid(5, {mat, mat, mat}, {}, {tran1, tran2}, {gene1, gene2});
*/

/*
P125:6-3
Eigen::MatrixXf mat(3, 5);

    mat << 2, 3, 0.024, 0.065, 0.032,
        2, 4, 0.030, 0.080, 0.040,
        4, 3, 0.018, 0.050, 0.026;

    Grid my_grid(5, {mat, mat, mat}, {{{1, 2, 1.00}, cf{0, 0.105}}, {{5, 4, 1.00}, cf{0, 0.184}}}, {}, {{1, 120, cf{0, 0.15}, 120}, {5, 120, cf{0, 0.22}, 120}});

    const auto ret = my_grid.SymmetricShortCircuit(3, cf(0.0f, 0.0f), 1);
    std::cout << "各点电压分布: " << std::endl << std::get<0>(ret) << std::endl << std::endl;

    std::cout << "各点电流分布: " << std::endl;
    const auto &I_list = std::get<1>(ret);
    for(const auto &iter : I_list)
    {
        std::cout << "I" << iter.first.first << iter.first.second << ": " << std::abs(iter.second) << " A" << std::endl;
    }
*/
