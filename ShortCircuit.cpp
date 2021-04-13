/*
 * @Author: your name
 * @Date: 2021-04-08 19:54:54
 * @LastEditTime: 2021-04-11 21:56:56
 * @LastEditors: Please set LastEditors
 * @Description: In User Settings Edit
 * @FilePath: /VSCode/ShortCircuit.cpp
 */

//-lgomp -lpthread
#define EIGEN_USE_MKL_ALL

#include "Common.h"
#include "Database.h"
#include "Grid.h"
#include "readMySQL.h"
#include <cstdlib>
#include <iostream>

int main()
{
    DataFetcher util("wadegao.tpddns.net", "root", "Grid");
    const auto &LineData = util.getLineData();
    const auto &IdealTransformer2Data = util.getIdealTransWithReactanceList();
    const auto &Transformer2Data = util.getTransformer2List(120);
    const auto &GeneratorData = util.getGeneratorList(120);

    Grid my_grid(5, {LineData, LineData, LineData}, IdealTransformer2Data, Transformer2Data, GeneratorData);

    const auto ret = my_grid.SymmetricShortCircuit(3, cf(0.0f, 0.0f), 1);

    std::cout << "各点电压分布: " << std::endl << std::get<0>(ret) << std::endl << std::endl;

    std::cout << "各点电流分布: " << std::endl;
    const auto &I_list = std::get<1>(ret);
    for(const auto &iter : I_list)
    {
        std::cout << "I" << iter.first.first << iter.first.second << ": " << std::abs(iter.second)  << std::endl;
    }

    return 0;
}
/*
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
