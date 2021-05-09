/*
 * @Author: your name
 * @Date: 2021-04-11 21:43:56
 * @LastEditTime: 2021-05-09 10:34:41
 * @LastEditors: Please set LastEditors
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
    std::vector<Eigen::Triplet<cf>> getLineTripletList();
    std::vector<Eigen::Triplet<cf>> getIdealTransWithTripletReactanceList();
    //std::vector<Transformer2> getTransformer2List(const DeviceArgType SB);
    std::vector<Eigen::Triplet<cf>> getGeneratorTripletList(const DeviceArgType SB);
    std::vector<Eigen::Triplet<cf>> getNodeTripletList();
};
/*
const std::string queryLine = "select tapBusNo, ZbusNo, R, X, B from branch where transformerFinalTurnsRatio = 0;";
const std::string queryTransformer = "select tapBusNo, ZbusNo, R, X, transformerFinalTurnsRatio from branch where transformerFinalTurnsRatio <> 0;";
const std::string queryGenerator = "select bus, generationMW, generationMVAR from bus where (generationMW <> 0 or generationMVAR <> 0);";
const std::string queryNode = "select bus, G, B FROM bus WHERE (G <> 0 or B <> 0) and (generationMW = 0 and generationMVAR = 0);";
*/
const std::string queryLine = "select fbus, tbus, r, x, b from branch where ratio = 0;";
const std::string queryTransformer = "select fbus, tbus, r, x, ratio from branch where ratio <> 0 and status = 1;";
const std::string queryGenerator = "select bus, Pg, Qg from bus where (Pg <> 0 or Qg <> 0) and status = 1;";
const std::string queryNode = "select bus_i, Gs, Bs FROM bus WHERE (Gs <> 0 or Bs <> 0);";
