/*
 * @Author: your name
 * @Date: 2021-04-11 21:50:11
 * @LastEditTime: 2021-04-11 21:55:51
 * @LastEditors: your name
 * @Description: In User Settings Edit
 * @FilePath: /VSCode/readMySQL.cpp
 */
#include "readMySQL.h"

DataFetcher::DataFetcher(const char *host, const char *user, const char *db, unsigned int port)
{
    if (!this->db.ConnectMySQL(host, user, db, port))
    {
        fprintf(stdout, "Failed to connect to database\n");
        exit(-1);
    }
}

std::vector<std::vector<float>> DataFetcher::getData(const std::string &query)
{
    const auto &TableSize = this->db.getQuerySize(query);
    //TODO: 这里，只要矩阵维度参数有一个是0，就不行，我已经定位问题了
    if(!std::get<0>(TableSize))
        return {};

    std::vector<std::vector<float>> ret(std::get<0>(TableSize), std::vector(std::get<1>(TableSize), 0.00f));
    const DatabaseTableType &Data = this->db.ReadMySQL(query);

#pragma omp parallel for
    for (decltype(Data.size()) row = 0; row < Data.size(); row++)
    {
        const auto &curRowData = Data.at(row);
        for (decltype(curRowData.size()) col = 0; col < curRowData.size(); col++)
            ret[row][col] = atof(curRowData.at(col).c_str());
    }

    return ret;
}

Eigen::MatrixXf DataFetcher::getLineData()
{
    const auto &ret = this->getData(queryLine);
    auto Rows = ret.size(), Cols = ret.at(0).size();
    Eigen::MatrixXf LineData = Eigen::MatrixXf::Zero(Rows, Cols);

#pragma omp parallel for
    for(decltype(Rows) i = 0; i < Rows; i++){
        const auto &curRow = ret.at(i);
        for(decltype(Cols) j = 0; j < Cols; j++)
            LineData(i, j) = curRow.at(j);
    }
    return LineData;
}

std::vector<std::pair<IdealTransformer2, cf>> DataFetcher::getIdealTransWithReactanceList()
{
    const auto &rawData = this->getData(queryTransformer);
    if(!rawData.size()) return {};
    std::vector<std::pair<IdealTransformer2, cf>> ret(rawData.size(), {{1, 2, 1.00}, {0, 0}});

#pragma omp parallel for
    for (decltype(ret.size()) i = 0; i < ret.size(); i++)
    {
        IdealTransformer2 thisTrans(rawData.at(i).at(0), rawData.at(i).at(1), rawData.at(i).at(5));
        cf priZ{rawData.at(i).at(2), rawData.at(i).at(3)};
        ret[i] = std::make_pair(thisTrans, priZ);
    }
    return ret;
}

/*
std::vector<Transformer2> DataFetcher::getTransformer2List(const DeviceArgType SB)
{
    const auto &rawData = this->getData("SELECT * FROM Transformer2");
    if(!rawData.size()) return {};
    std::vector<Transformer2> ret(rawData.size(), {1, 2, 1, 1, 1, 1});

#pragma omp parallel for
    for (decltype(ret.size()) i = 0; i < ret.size(); i++)
        ret[i] = Transformer2(rawData.at(i).at(0), rawData.at(i).at(1), rawData.at(i).at(2), rawData.at(i).at(3), rawData.at(i).at(4), SB);

    return ret;
}
*/

std::vector<Generator> DataFetcher::getGeneratorList(const DeviceArgType SB)
{
    const auto &rawData = this->getData(queryGenerator);
    if(!rawData.size()) return {};
    std::vector<Generator> ret(rawData.size(), {1, 1, {0, 0}, 1});

#pragma omp parallel for
    for (decltype(ret.size()) i = 0; i < ret.size(); i++){
        const auto &curRowData = rawData.at(i);
        ret[i] = Generator(curRowData.at(0), std::abs(cf(curRowData.at(1), curRowData.at(2))), cf(0.0f, 0.0f), SB);
    }
    return ret;
}

std::vector<std::tuple<NodeType, float, float>> DataFetcher::getNodeArgList()
{
    const auto rawData = this->getData(queryNode);
    if(!rawData.size())    return {};
    std::vector<std::tuple<NodeType, float, float>> ret(rawData.size(), {0, 0.0f, 0.0f});
#pragma omp parallel for
    for (decltype(ret.size()) i = 0; i < ret.size(); i++){
        const auto &curRowData = rawData.at(i);
        ret[i] = {curRowData.at(0), curRowData.at(1), curRowData.at(2)};
    }
    return ret;
}
