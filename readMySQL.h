/*
 * @Author: your name
 * @Date: 2021-04-11 21:43:56
 * @LastEditTime: 2021-04-11 21:56:12
 * @LastEditors: your name
 * @Description: In User Settings Edit
 * @FilePath: /VSCode/readMySQL.h
 */
#pragma once
#include "Common.h"
#include "Database.h"
#include "Equipment.h"
#include "Generator.h"
#include "Transformer.h"

class DataFetcher
{
private:
    Database db;
    std::vector<std::vector<float>> getData(const std::string &query);

public:
    DataFetcher(const char *host, const char *user, const char *db, unsigned int port = 3306);
    Eigen::MatrixXf getLineData();
    std::vector<std::pair<IdealTransformer2, cf>> getIdealTransWithReactanceList();
    //std::vector<Transformer2> getTransformer2List(const DeviceArgType SB);
    std::vector<Generator> getGeneratorList(const DeviceArgType SB);
    std::vector<std::tuple<NodeType, float, float>> getNodeArgList();
};

const std::string queryLine = "select fbus, tbus, r, x, b from branch where status = 1 and ratio = 0;";
const std::string queryTransformer = "select tbus, fbus, ratio, r, x from branch where status = 1 and ratio <> 0;";
const std::string queryGenerator = "SELECT bus, mBase as Sn from generator where status = 1;";
const std::string queryNode = "SELECT bus_i, Gs, Bs FROM bus WHERE Gs <> 0 OR Bs <> 0;";
