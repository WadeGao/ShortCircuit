#pragma once
#include "Common.h"
#include "ThreadPool.hpp"

class Grid {
private:
	ThreadPool myPool;
	Eigen::MatrixXcf Z1, Z2, Z0;
	Eigen::MatrixXcf Y1, Y2, Y0;
	size_t NodesNum{ 0 };

	//ͨ�����ݼ���õ��ɾ���
	static Eigen::MatrixXcf getAdmittanceMatrixBySheet(const Eigen::MatrixXf& line_data_sheet);

	//ͨ���迹�����õ��ɾ���
	static Eigen::MatrixXcf getAdmittanceMatrixFromReactanceMatrix(const Eigen::MatrixXcf& ReactanceMatrix);

	//�Գƶ�·����
	void SymmetricShortCircuit(const Eigen::MatrixXcf& AdmittanceMatrix, const size_t shortPoint, const float Sb, const float Uav);

	//�����·����
	void lgShortCircuit(const size_t faultNode, const std::complex<float>& Zf);

	//�����·����
	void llShortCircuit(const size_t faultNode, const std::complex<float>& Zf);

	//����Եض�·����
	void llgShortCircuit(const size_t faultNode, const std::complex<float>& Zf);

	void getBusVoltageAndCurrent(const Eigen::VectorXcf& If, const size_t faultNode);

public:
	Grid(const Eigen::MatrixXcf& Z1_, const Eigen::MatrixXcf& Z2_, const Eigen::MatrixXcf& Z0_);
	~Grid();
};