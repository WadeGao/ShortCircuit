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

    /*mat << 1, 2, 0.000, 0.105, 0.000, 1.05,
        2, 3, 0.024, 0.065, 0.032, 1.00,
        2, 4, 0.030, 0.080, 0.040, 1.00,
        4, 3, 0.018, 0.050, 0.026, 1.00,
        5, 4, 0.000, 0.184, 0.000, 0.96;*/

    /*mat <<
        2, 3, 0.024, 0.065, 0.032,
        2, 4, 0.030, 0.080, 0.040,
        4, 3, 0.018, 0.050, 0.026;*/
    //IdealTransformer2 tran1(5, 4, 0.96), tran2(1, 2, 1.05);
    //Grid my_grid(5, {mat, mat, mat}, {{tran1, cf{0, 0.184}}, {tran2, cf{0, 0.105}}}, {}, {});

    mat <<
    3, 4, 0.000, 0.43554, 0.018515,
    3, 5, 0.000, 0.29036, 0.012343,
    4, 5, 0.000, 0.25407, 0.010800;

    Transformer2 tran1(1, 3, 1, 120, 10.5, 120), tran2(2, 4, 1, 60, 10.5, 120);
    Generator gene1(1, 120, 0.23, 120), gene2(2, 60, 0.14, 120);

    Grid my_grid(5, {mat, mat, mat}, {}, {tran1, tran2}, {gene1, gene2});
    std::cout << my_grid.getYx(SEQUENCE::POSITIVE) << std::endl;
    return 0;
}
