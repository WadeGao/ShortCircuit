/*
 * @Author: your name
 * @Date: 2021-04-08 19:08:26
 * @LastEditTime: 2021-04-26 09:48:05
 * @LastEditors: Please set LastEditors
 * @Description: In User Settings Edit
 * @FilePath: /VSCode/Grid.h
 */
#pragma once
#include "Common.h"
#include "Generator.h"
#include "ThreadPool.hpp"
#include "Transformer.h"

enum class SEQUENCE
{
    POSITIVE = 1,
    NEGATIVE,
    ZERO
};

using lllReturnType = std::tuple<Eigen::VectorXcf, std::list<std::pair<socketType, cf>>>;

class Grid
{
private:
    NodeType NodeNum{0};
    //ThreadPool myPool;
    Eigen::MatrixXcf Z1, Z2, Z0;
    Eigen::SparseMatrix<cf> Y1, Y2, Y0;
    //socket->Y
    std::map<socketType, cf> SocketData1{}, SocketData2{}, SocketData0{};

    static void wrapper(const TripVecType &triList, std::map<socketType, cf> &sockMap);
    void getBusVoltageAndCurrent(const Eigen::VectorXcf &If, const NodeType faultNode);

public:
    Grid(const NodeType node_, const TripVecType &lineData, const TripVecType &nodeList, const TripVecType &idealTransList, const TripVecType &transList, const TripVecType &geneList);
    ~Grid();

    Eigen::MatrixXcf getYx(const SEQUENCE whichSeq) const;
    Eigen::MatrixXcf getZx(const SEQUENCE whichSeq) const;

    lllReturnType SymmetricShortCircuit(const NodeType shortPoint, const cf &Zf, const DeviceArgType Uav);

    void lgShortCircuit(const NodeType faultNode, const cf &Zf);

    void llShortCircuit(const NodeType faultNode, const cf &Zf);

    void llgShortCircuit(const NodeType faultNode, const cf &Zf);
};
