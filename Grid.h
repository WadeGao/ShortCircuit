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

    //ͨ�����ݼ������·ԭʼ���ɾ���
    Eigen::MatrixXcf setYxFromSheet(const Eigen::MatrixXf& line_data_sheet, std::map<std::pair<NodeType, NodeType>, std::pair<cf, cf>>& socketData);

    //������ѹ����ȣ�ͬʱ�Խڵ㵼�ɾ��������޸�
    void adjustTransformerRatio(const std::vector<Transformer2> &transList);

    //�Գƶ�·����
    void SymmetricShortCircuit(const Eigen::MatrixXcf& AdmittanceMatrix, const NodeType shortPoint, const DeviceArgType Sb, const DeviceArgType Uav);

    //�����·����
    void lgShortCircuit(const NodeType faultNode, const cf& Zf);

    //�����·����
    void llShortCircuit(const NodeType faultNode, const cf& Zf);

    //����Եض�·����
    void llgShortCircuit(const NodeType faultNode, const cf& Zf);

    //�ǶԳƶ�·ʱ������·��ѹ����
    void getBusVoltageAndCurrent(const Eigen::VectorXcf& If, const NodeType faultNode);

public:
    Grid(const ThreeSequenceData& data, const std::vector<Transformer2> &transformList);
    ~Grid();

    //getter����
    Eigen::MatrixXcf getYx(const SEQUENCE whichSeq) const;
    Eigen::MatrixXcf getZx(const SEQUENCE whichSeq) const;
};

