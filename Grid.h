#pragma once
#include "Common.h"
#include "Transformer.h"
#include "ThreadPool.hpp"

enum class SEQUENCE
{
    POSITIVE = 1,
    NEGATIVE,
    ZERO
};

class Grid
{
private:
    ThreadPool myPool;
    Eigen::MatrixXcf Z1, Z2, Z0;
    Eigen::MatrixXcf Y1, Y2, Y0;
    //socket->{Y, B}
    std::map<std::pair<NodeType, NodeType>, std::pair<cf, cf>> SocketData1{}, SocketData2{}, SocketData0{};
    NodeType NodesNum{ 0 };

    //通过数据集获得线路原始导纳矩阵
    Eigen::MatrixXcf setYxFromSheet(const Eigen::MatrixXf& line_data_sheet, std::map<std::pair<NodeType, NodeType>, std::pair<cf, cf>>& socketData);

    //调整变压器变比，同时对节点导纳矩阵作出修改
    void adjustTransformerRatio(const std::vector<Transformer2> &transList);

    //对称短路计算
    void SymmetricShortCircuit(const Eigen::MatrixXcf& AdmittanceMatrix, const NodeType shortPoint, const DeviceArgType Sb, const DeviceArgType Uav);

    //单相短路计算
    void lgShortCircuit(const NodeType faultNode, const cf& Zf);

    //两相短路计算
    void llShortCircuit(const NodeType faultNode, const cf& Zf);

    //两相对地短路计算
    void llgShortCircuit(const NodeType faultNode, const cf& Zf);

    //非对称短路时计算线路电压电流
    void getBusVoltageAndCurrent(const Eigen::VectorXcf& If, const NodeType faultNode);

public:
    Grid(const ThreeSequenceData& data, const std::vector<Transformer2> &transformList);
    ~Grid();

    //getter方法
    Eigen::MatrixXcf getYx(const SEQUENCE whichSeq) const;
    Eigen::MatrixXcf getZx(const SEQUENCE whichSeq) const;
};

