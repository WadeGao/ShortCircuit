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
    std::vector<std::vector<float>> getData(const char *TableName);

public:
    DataFetcher(const char *host, const char *user, const char *db, unsigned int port = 3306);
    Eigen::MatrixXf getLineData();
    std::vector<std::pair<IdealTransformer2, cf>> getIdealTransWithReactanceList();
    std::vector<Transformer2> getTransformer2List(const DeviceArgType SB);
    std::vector<Generator> getGeneratorList(const DeviceArgType SB);
};
