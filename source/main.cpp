/*
 * @Author: your name
 * @Date: 2021-04-08 19:54:54
 * @LastEditTime: 2021-05-07 10:37:44
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
#include <ctime>
#include <iostream>

int main()
{
    DataFetcher util("wadegao.tpddns.cn", "root", "IEEE300");
    const auto &LineData = util.getLineTripletList();
    const auto &IdealTransformer2Data = util.getIdealTransWithTripletReactanceList();
    //const auto &Transformer2Data = util.getTransformer2List(120);
    const auto &GeneratorData = util.getGeneratorTripletList(120);
    const auto &NodeData = util.getNodeTripletList();

    Grid my_grid(300, LineData, NodeData, IdealTransformer2Data, {}, GeneratorData);

    const auto ret = my_grid.SymmetricShortCircuit(13, cf(0.0f, 0.0f), 1);
    auto task = [&my_grid](NodeType node, const cf &z, const DeviceArgType u) -> lllReturnType { return my_grid.SymmetricShortCircuit(node, z, u); };

    auto startTime = clock();
    auto results = my_grid.lllWholeGridScan();

    size_t len = 0;
    for (auto &&res : results)
        len += std::get<1>(res.get()).size();

    std::cout << "len = " << len << std::endl;
    printf("time: %lfs\n", (double)(clock() - startTime) / CLOCKS_PER_SEC);

    startTime = clock();
    size_t ans{0};
    for (int i = 1; i < 300; i++)
    {
        const auto &ret = my_grid.SymmetricShortCircuit(i, cf(0.0f, 0.0f), 1);
        const auto &I_list = std::get<1>(ret);
        ans += I_list.size();
    }
    std::cout << "ans = " << ans << std::endl;
    printf("time: %lfs\n", (double)(clock() - startTime) / CLOCKS_PER_SEC);

    /*for (auto &&result : results)
    {
        auto ret = result.get();
        std::cout << "各点电压分布: " << std::endl
                  << std::get<0>(ret) << std::endl
                  << std::endl;

        std::cout << "各点电流分布: " << std::endl;
        const auto &I_list = std::get<1>(ret);
        for (const auto &iter : I_list)
            std::cout << "I[" << iter.first.first + 1 << "][" << iter.first.second + 1 << "]: (" << iter.second.real() << ", " << iter.second.imag() << ")" << std::endl;
    }*/

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
