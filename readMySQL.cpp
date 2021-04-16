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

std::vector<Eigen::Triplet<cf>> DataFetcher::getLineTripletList()
{
    const auto &rawData = this->getData(queryLine);
    auto Rows = rawData.size();
    if(!Rows)    return {};

    std::vector<Eigen::Triplet<cf>> ret(3 * Rows, {0, 0, cf(0, 0)});
#pragma omp parallel for
    for(decltype(Rows) i = 0; i < Rows; i++)
    {
        const auto &curRow = rawData.at(i);

        const int from = curRow.at(0) - 1, to = curRow.at(1) - 1;
        const auto &y = (cf(1, 0) / cf(curRow.at(2), curRow.at(3)));
        const auto &B = cf(0, curRow.at(4));
        const auto Y_self = y + B / cf(2, 0);

        ret[3 * i] = Eigen::Triplet<cf>{std::max(from, to), std::min(from, to), cf(0, 0) - y};
        ret[3 * i + 1] = Eigen::Triplet<cf>{from, from, Y_self};
        ret[3 * i + 2] = Eigen::Triplet<cf>{to, to, Y_self};
    }
    return ret;
}

std::vector<Eigen::Triplet<cf>> DataFetcher::getIdealTransWithTripletReactanceList()
{
    const auto &rawData = this->getData(queryTransformer);
    if(!rawData.size())    return {};
    std::vector<Eigen::Triplet<cf>> ret(3 * rawData.size(), {0, 0, cf(0, 0)});

#pragma omp parallel for
    for (decltype(rawData.size()) i = 0; i < rawData.size(); i++)
    {
        const auto &curRow = rawData.at(i);

        const int pNode = curRow.at(0) - 1, qNode = curRow.at(1) - 1;
        const auto &k = curRow.at(4);
        const auto &y = cf{1, 0} / cf{curRow.at(2), curRow.at(3)};

        ret[3 * i] = Eigen::Triplet<cf>{std::max(pNode, qNode), std::min(pNode, qNode), cf(0, 0) - y / k};
        ret[3 * i + 1] = Eigen::Triplet<cf>{pNode, pNode, y};
        ret[3 * i + 2] = Eigen::Triplet<cf>{qNode, qNode, y / (k * k)};

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

std::vector<Eigen::Triplet<cf>> DataFetcher::getGeneratorTripletList(const DeviceArgType SB)
{
    const auto &rawData = this->getData(queryGenerator);
    if(!rawData.size())    return {};
    std::vector<Eigen::Triplet<cf>> ret(rawData.size(), {0, 0, cf(0, 0)});

#pragma omp parallel for
    for (decltype(ret.size()) i = 0; i < ret.size(); i++){
        const auto &curRowData = rawData.at(i);
        // Node(Node_), Sn(Sn_), xd_(__xd), Xd(__xd * SB_ / Sn_)
        const int node = curRowData.at(0) - 1;
        const auto Sn = std::abs(cf(curRowData.at(1), curRowData.at(2)));
        const auto __xd = cf(0.0f, 0.0f);
        const auto Xd = __xd * SB / Sn;
        const auto y = (std::abs(Xd) > epsilon) ? (cf{1, 0} / Xd) : cf(0, 0);
        ret[i] = Eigen::Triplet<cf>{node, node, y};
    }

    return ret;
}

std::vector<Eigen::Triplet<cf>> DataFetcher::getNodeTripletList()
{
    const auto rawData = this->getData(queryNode);
    if(!rawData.size())    return {};
    std::vector<Eigen::Triplet<cf>> ret(rawData.size(), {0, 0, cf(0, 0)});
#pragma omp parallel for
    for (decltype(ret.size()) i = 0; i < ret.size(); i++){
        const auto &curRowData = rawData.at(i);
        const int &node = curRowData.at(0);
        const auto Gs_add_Bs = cf(0.0f, curRowData.at(1) + curRowData.at(2));
        ret[i] = Eigen::Triplet<cf>{node, node, Gs_add_Bs};
    }
    return ret;
}
