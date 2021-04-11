/*
 * @Author: your name
 * @Date: 2021-04-08 19:08:26
 * @LastEditTime: 2021-04-11 21:56:03
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

class Grid
{
private:
    NodeType NodeNum{0};
    ThreadPool myPool;
    Eigen::MatrixXcf Z1, Z2, Z0;
    Eigen::MatrixXcf Y1, Y2, Y0;
    //socket->{Y, B}
    std::map<socketType, std::pair<cf, cf>> SocketData1{}, SocketData2{}, SocketData0{};

    void setYxFromSheet(const Eigen::MatrixXf &line_data_sheet, Eigen::MatrixXcf &Y, std::map<socketType, std::pair<cf, cf>> &socketData);

    void adjustIdealTransformer2_PrimarySideReactanceRatio(const std::vector<Transformer2> &transList);

    void mountTransformer2(const std::vector<Transformer2> &transList);

    void mountIdealTransformer2_PrimarySideReactance(const std::vector<std::pair<IdealTransformer2, cf>> &idealTransList);

    void mountGenerator(const std::vector<Generator> &geneList);

    void getBusVoltageAndCurrent(const Eigen::VectorXcf &If, const NodeType faultNode);

public:
    Grid(const NodeType node_, const ThreeSequenceData &data, const std::vector<std::pair<IdealTransformer2, cf>> &idealTransformList, const std::vector<Transformer2> &transList, const std::vector<Generator> &geneList);
    ~Grid();

    Eigen::MatrixXcf getYx(const SEQUENCE whichSeq) const;
    Eigen::MatrixXcf getZx(const SEQUENCE whichSeq) const;

    std::tuple<Eigen::VectorXcf, std::list<std::pair<socketType, cf>>> SymmetricShortCircuit(const NodeType shortPoint, const cf &Zf, const DeviceArgType Uav);

    void lgShortCircuit(const NodeType faultNode, const cf &Zf);

    void llShortCircuit(const NodeType faultNode, const cf &Zf);

    void llgShortCircuit(const NodeType faultNode, const cf &Zf);
};
