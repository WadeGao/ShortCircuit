#pragma once
#include "Common.h"
#include "ThreadPool.hpp"

class Grid {
private:
	ThreadPool myPool;
	Eigen::MatrixXcf Z1, Z2, Z0;
	Eigen::MatrixXcf Y1, Y2, Y0;
	size_t NodesNum{ 0 };

	//通过阻抗矩阵获得导纳矩阵
	static Eigen::MatrixXcf getAdmittanceMatrixFromReactanceMatrix(const Eigen::MatrixXcf& ReactanceMatrix);

	//对称短路计算
	void SymmetricShortCircuit(const Eigen::MatrixXcf& AdmittanceMatrix, const size_t shortPoint, const float Sb, const float Uav);

	//单相短路计算
	void lgShortCircuit(const size_t faultNode, const std::complex<float>& Zf);

	//两相短路计算
	void llShortCircuit(const size_t faultNode, const std::complex<float>& Zf);

	//两相对地短路计算
	void llgShortCircuit(const size_t faultNode, const std::complex<float>& Zf);

    //非对称短路时计算线路电压电流
	void getBusVoltageAndCurrent(const Eigen::VectorXcf& If, const size_t faultNode);

	//调整变压器变比，同时对节点导纳矩阵作出修改
	void adjustTransformerRatio(const size_t PrimaryWindingNode, const size_t SecondaryWindingNode, const float newRatio, const std::complex<float> &z);

public:
	Grid(const Eigen::MatrixXcf& Z1_, const Eigen::MatrixXcf& Z2_, const Eigen::MatrixXcf& Z0_);
	~Grid();

	//通过数据集获得导纳矩阵
	static Eigen::MatrixXcf getAdmittanceMatrixBySheet(const Eigen::MatrixXf& line_data_sheet);
};
