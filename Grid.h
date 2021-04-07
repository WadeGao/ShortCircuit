#pragma once
#include "Common.h"
#include "Transformer.h"
#include "ThreadPool.hpp"

enum class SEQUENCE {
	POSITIVE = 1,
	NEGATIVE,
	ZERO
};

class Grid {
private:
	ThreadPool myPool;
	Eigen::MatrixXcf Z1, Z2, Z0;
	Eigen::MatrixXcf Y1, Y2, Y0;
	//socket->{Y, B}
	std::map<std::pair<size_t, size_t>, std::pair<cf, cf>> SocketData1{}, SocketData2{}, SocketData0{};
	size_t NodesNum{ 0 };

	//通过阻抗矩阵获得导纳矩阵
	//static Eigen::MatrixXcf getAdmittanceMatrixFromReactanceMatrix(const Eigen::MatrixXcf& ReactanceMatrix);

	//通过数据集获得导纳矩阵
	Eigen::MatrixXcf getAdmittanceMatrixBySheet(const Eigen::MatrixXf& line_data_sheet, std::map<std::pair<size_t, size_t>, std::pair<cf, cf>>& socketData);

	//对称短路计算
	void SymmetricShortCircuit(const Eigen::MatrixXcf& AdmittanceMatrix, const size_t shortPoint, const float Sb, const float Uav);

	//单相短路计算
	void lgShortCircuit(const size_t faultNode, const cf& Zf);

	//两相短路计算
	void llShortCircuit(const size_t faultNode, const cf& Zf);

	//两相对地短路计算
	void llgShortCircuit(const size_t faultNode, const cf& Zf);

    //非对称短路时计算线路电压电流
	void getBusVoltageAndCurrent(const Eigen::VectorXcf& If, const size_t faultNode);

	//调整变压器变比，同时对节点导纳矩阵作出修改
	void adjustTransformerRatio(const std::list<Transformer2> &transList);

public:
	Grid(const ThreeSequenceData& data, const std::list<Transformer2>& transformList);
	~Grid();
	
	Eigen::MatrixXcf getAdmittanceMatrix(const SEQUENCE whichSeq);
};

