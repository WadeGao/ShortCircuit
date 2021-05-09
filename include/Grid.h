/*
 * @Author: your name
 * @Date: 2021-04-08 19:08:26
 * @LastEditTime: 2021-05-09 20:47:01
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
using SequenceVoltType = std::tuple<Eigen::MatrixXcf, Eigen::MatrixXcf, Eigen::MatrixXcf>;
using SequenceCurrentType = std::map<std::pair<int, int>, std::vector<cf>>;

class Grid
{
private:
    ThreadPool myPool;

    NodeType NodeNum{0};
    //ThreadPool myPool;c
    Eigen::MatrixXcf Z1, Z2, Z0;
    Eigen::SparseMatrix<cf> Y1, Y2, Y0;
    //socket->Y
    std::map<socketType, cf> SocketData1{}; //, SocketData2{}, SocketData0{};

    static void wrapper(const TripVecType &triList, std::map<socketType, cf> &sockMap);
    //void getBusVoltageAndCurrent(const Eigen::VectorXcf &If, const NodeType faultNode);

public:
    Grid(const NodeType node_, const TripVecType &lineData, const TripVecType &nodeList, const TripVecType &idealTransList, const TripVecType &transList, const TripVecType &geneList);
    ~Grid();

    Eigen::MatrixXcf getYx(const SEQUENCE whichSeq) const;
    Eigen::MatrixXcf getZx(const SEQUENCE whichSeq) const;

    lllReturnType SymmetricShortCircuit(const NodeType shortPoint, const cf &Zf, const DeviceArgType Uav);

    cf lgShortCircuit(const NodeType faultNode, const cf &Zf);

    cf llShortCircuit(const NodeType faultNode, const cf &Zf);

    cf llgShortCircuit(const NodeType faultNode, const cf &Zf);

    std::list<std::future<lllReturnType>> lllWholeGridScan();
    std::list<std::future<cf>> llgWholeGridScan();
    std::list<std::future<cf>> lgWholeGridScan();
    std::list<std::future<cf>> llWholeGridScan();
    SequenceVoltType getSequenceVolt(const Eigen::VectorXcf &If, const NodeType faultNode);
    Eigen::MatrixXf getBusVoltage(const Eigen::VectorXcf &If, const NodeType faultNode);
    SequenceCurrentType getBusCurrent(const Eigen::VectorXcf &If, const NodeType faultNode);
    NodeType getNodeNum() const;
};
