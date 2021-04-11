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

Eigen::MatrixXf DataFetcher::getData(const char *TableName)
{
    const auto &TableSize = this->db.getTableSize(TableName);
    //TODO: 这里，只要矩阵维度参数有一个是0，就不行，我已经定位问题了
    Eigen::MatrixXf ret(std::get<0>(TableSize), std::get<1>(TableSize));
    const auto &Data = this->db.ReadMySQL("SELECT * FROM Line;");

#pragma omp parallel for
    for (decltype(Data.size()) row = 0; row < Data.size(); row++)
    {
        const auto &curRowData = Data.at(row);
        for (decltype(curRowData.size()) col = 0; col < curRowData.size(); col++)
            ret(row, col) = atof(curRowData.at(col).c_str());
    }
    std::cout << __LINE__ << std::endl;
    return ret;
}

Eigen::MatrixXf DataFetcher::getLineData()
{
    const auto &ret = this->getData("Line");
    return ret;
}


std::vector<std::pair<IdealTransformer2, cf>> DataFetcher::getIdealTransWithReactanceList()
{

    const auto &rawData = this->getData("IdealTransformer2");
    if(!rawData.rows()) return {};
    std::vector<std::pair<IdealTransformer2, cf>> ret(rawData.rows(), {{1, 2, 1.00}, {0, 0}});
    //ret.resize(rawData.rows());

#pragma omp parallel for
    for (decltype(ret.size()) i = 0; i < ret.size(); i++)
    {
        IdealTransformer2 thisTrans(rawData(i, 0), rawData(i, 1), rawData(i, 2));
        cf priZ{rawData(i, 3), rawData(i, 4)};
        ret[i] = std::make_pair(thisTrans, priZ);
    }
    return ret;
}


std::vector<Transformer2> DataFetcher::getTransformer2List(const DeviceArgType SB)
{
    const auto &rawData = this->getData("Transformer2");
    if(!rawData.rows()) return {};
    std::vector<Transformer2> ret(rawData.rows(), {1, 2, 1, 1, 1, 1});
    //ret.resize(rawData.rows());

#pragma omp parallel for
    for (decltype(ret.size()) i = 0; i < ret.size(); i++)
        ret[i] = Transformer2(rawData(i, 0), rawData(i, 1), rawData(i, 2), rawData(i, 3), rawData(i, 4), SB);

    return ret;
}

std::vector<Generator> DataFetcher::getGeneratorList(const DeviceArgType SB)
{
    const auto &rawData = this->getData("Generator");
    if(!rawData.rows()) return {};
    std::vector<Generator> ret(rawData.rows(), {1, 1, {0, 0}, 1});

#pragma omp parallel for
    for (decltype(ret.size()) i = 0; i < ret.size(); i++)
        ret[i] = Generator(rawData(i, 0), rawData(i, 1), cf{rawData(i, 2), rawData(i, 3)}, SB);

    return ret;
}
