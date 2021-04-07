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

	//ͨ���迹�����õ��ɾ���
	//static Eigen::MatrixXcf getAdmittanceMatrixFromReactanceMatrix(const Eigen::MatrixXcf& ReactanceMatrix);

	//ͨ�����ݼ���õ��ɾ���
	Eigen::MatrixXcf getAdmittanceMatrixBySheet(const Eigen::MatrixXf& line_data_sheet, std::map<std::pair<size_t, size_t>, std::pair<cf, cf>>& socketData);

	//�Գƶ�·����
	void SymmetricShortCircuit(const Eigen::MatrixXcf& AdmittanceMatrix, const size_t shortPoint, const float Sb, const float Uav);

	//�����·����
	void lgShortCircuit(const size_t faultNode, const cf& Zf);

	//�����·����
	void llShortCircuit(const size_t faultNode, const cf& Zf);

	//����Եض�·����
	void llgShortCircuit(const size_t faultNode, const cf& Zf);

    //�ǶԳƶ�·ʱ������·��ѹ����
	void getBusVoltageAndCurrent(const Eigen::VectorXcf& If, const size_t faultNode);

	//������ѹ����ȣ�ͬʱ�Խڵ㵼�ɾ��������޸�
	void adjustTransformerRatio(const std::list<Transformer2> &transList);

public:
	Grid(const ThreeSequenceData& data, const std::list<Transformer2>& transformList);
	~Grid();
	
	Eigen::MatrixXcf getAdmittanceMatrix(const SEQUENCE whichSeq);
};

