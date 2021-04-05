#include "Grid.h"
#include <iostream>

//ͨ�����ݼ���õ��ɾ���
Eigen::MatrixXcf Grid::getAdmittanceMatrixBySheet(const Eigen::MatrixXf& line_data_sheet)
{
	//dataSheet: https://img-blog.csdn.net/20180616204918263?watermark/2/text/aHR0cHM6Ly9ibG9nLmNzZG4ubmV0L3FxXzMyNDEyNzU5/font/5a6L5L2T/fontsize/400/fill/I0JBQkFCMA==/dissolve/70
	const Eigen::VectorXi inNode = line_data_sheet.col(0).cast<int>().array();
	const Eigen::VectorXi outNode = line_data_sheet.col(1).cast<int>().array();
	const Eigen::VectorXcf Z = line_data_sheet.col(2) + Imaginer * line_data_sheet.col(3);
	const Eigen::VectorXcf B = line_data_sheet.col(4) * Imaginer;
	const Eigen::VectorXf K = line_data_sheet.col(5);

	auto branchNum = inNode.size();
	auto busNum = std::max(inNode.maxCoeff(), outNode.maxCoeff());

	Eigen::VectorXcf y = Eigen::VectorXcf(branchNum).setOnes().array() / Z.array();
	Eigen::MatrixXcf Y = Eigen::MatrixXcf::Zero(busNum, busNum);


#pragma omp parallel for
	for (decltype(branchNum) i = 0; i < branchNum; i++)
	{
		auto from = inNode(i) - 1, to = outNode(i) - 1;
		//Y(from, to) -= K(i) * y(i);
		//Y(to, from) = Y(from, to);
		//Y(from, from) += K(i) * K(i) * y(i) + Imaginer * B(i);
		//Y(to, to) += y(i) + Imaginer * B(i);

		Y(from, to) -= y(i) / K(i);
		Y(to, from) = Y(from, to);
		Y(to, to) += y(i) / (K(i) * K(i)) + (B(i) / std::complex<float>{2, 0});
		Y(from, from) += y(i) + (B(i) / std::complex<float>{2, 0});
	}
	return Y;
}

//ͨ���迹����Zdx��õ��ɾ���
Eigen::MatrixXcf Grid::getAdmittanceMatrixFromReactanceMatrix(const Eigen::MatrixXcf& Zdx)
{
	//dataSheet: https://blog.csdn.net/qq_32412759/article/details/80715617
	Eigen::VectorXi inNode = Zdx.col(0).real().cast<int>().array();
	Eigen::VectorXi outNode = Zdx.col(1).real().cast<int>().array();
	Eigen::MatrixXcf yd = Eigen::VectorXcf(inNode.size()).setOnes().array() / Zdx.col(2).array();
	auto busNum = std::max(inNode.maxCoeff(), outNode.maxCoeff());
	Eigen::MatrixXcf Y = Eigen::MatrixXcf::Zero(busNum, busNum);

	//�������map����Ϊ��pair��Ϊunordered_map��keyʱ����Ҫhash��������map����Ҫ
	std::map<std::pair<int, int>, std::complex<float>> socket_yd{};
	std::unordered_map<int, std::unordered_set<int>> DirectConnectNodePair{};

#pragma omp parallel for
	for (size_t i = 0; i < inNode.size(); i++)
	{
		//��ʵ����ʡһ��ռ�ģ��ÿռ任ʱ��ɣ�����Ҳ�ǶԳƾ�����
		socket_yd.insert({ {inNode(i), outNode(i)}, yd(i) });
		socket_yd.insert({ {outNode(i), inNode(i)}, yd(i) });
		DirectConnectNodePair[inNode(i)].insert(outNode(i));
		DirectConnectNodePair[outNode(i)].insert(inNode(i));
	}

#pragma omp parallel for
	for (decltype(busNum) i = 0; i < busNum; i++)
	{
		for (decltype(i) j = 0; j < i; j++)
		{
			auto iter = socket_yd.find({ i + 1, j + 1 });
			if (iter != socket_yd.end())
			{
				//TODO:���ﵽ���Ǽӻ��Ǽ�
				//Y(i, j) += iter->second;
				Y(i, j) -= iter->second;
				Y(j, i) = Y(i, j);
			}
		}
	}

#pragma omp parallel for
	for (decltype(busNum) i = 0; i < busNum; i++)
	{
		auto iter = DirectConnectNodePair.find(i + 1);
		for (const auto& hashSetIter : iter->second)
		{
			auto MapIter = socket_yd.find({ i + 1, hashSetIter });
			Y(i, i) += MapIter->second;
		}
	}
	return Y;
}

//TODO:�Գƶ�·����
void Grid::SymmetricShortCircuit(const Eigen::MatrixXcf& AdmittanceMatrix, const size_t shortPoint, const float Sb, const float Uav)
{
	const auto busNum = AdmittanceMatrix.rows();
	Eigen::VectorXf Is = Eigen::VectorXf::Zero(busNum, 1);
}

//�����·����
void Grid::lgShortCircuit(const size_t faultNode, const std::complex<float>& Zf)
{
	auto faultNode2 = faultNode, faultNode0 = faultNode;

	std::complex<float> Ifa1 = std::complex<float>{ 1, 0 } / (Z1(faultNode - 1, faultNode - 1) + Z2(faultNode2 - 1, faultNode2 - 1) + Z0(faultNode0 - 1, faultNode0 - 1) + Zf * std::complex<float>{0, 3});
	Eigen::VectorXcf If(3, 1);
	//auto Ifa2{ Ifa1 }, Ifa0{ Ifa1 };
	//If << Ifa1, Ifa2, Ifa0;
	If << Ifa1, Ifa1, Ifa1;
	Eigen::VectorXcf If_abc = St * If;

	fprintf(stdout, "���ϵ�abc���·����: \n");
	for (int i = 0; i < If_abc.size(); i++)
		fprintf(stdout, "%f\n", abs(If_abc(i)));

	auto Is = If_abc(0);
	std::cout << "���ϵ���: " << Is << std::endl;

	this->getBusVoltageAndCurrent(If, faultNode);
}

//�����·����
void Grid::llShortCircuit(const size_t faultNode, const std::complex<float>& Zf)
{
	auto faultNode2 = faultNode;

	Eigen::VectorXcf If(3, 1);
	std::complex<float> Ifa1 = std::complex<float>{ 1, 0 } / (Z1(faultNode - 1, faultNode - 1) + Z2(faultNode2 - 1, faultNode2 - 1) + Zf * std::complex<float>{0, 1});
	std::complex<float> Ifa2{ -Ifa1.real(), -Ifa1.imag() }, Ifa0{ 0, 0 };
	If << Ifa1, Ifa2, Ifa0;
	Eigen::VectorXcf If_abc = St * If;

	fprintf(stdout, "���ϵ�abc���·����: \n");
	for (int i = 0; i < If_abc.size(); i++)
		fprintf(stdout, "%f\n", abs(If_abc(i)));

	auto Is = If_abc(2);
	std::cout << "���ϵ���: " << Is << std::endl;

	this->getBusVoltageAndCurrent(If, faultNode);
}

//����Եض�·����
void Grid::llgShortCircuit(const size_t faultNode, const std::complex<float>& Zf)
{
	auto faultNode2 = faultNode, faultNode0 = faultNode;

	Eigen::VectorXcf If(3, 1);
	std::complex<float> Ifa1 = std::complex<float>{ 1, 0 } / (Z1(faultNode - 1, faultNode - 1) + Z2(faultNode2 - 1, faultNode2 - 1) * (Z0(faultNode0 - 1, faultNode0 - 1) + Zf * std::complex<float>{0, 3}) / (Z2(faultNode2 - 1, faultNode2 - 1) + Z0(faultNode0 - 1, faultNode0 - 1) + Zf * std::complex<float>{0, 3}));
	std::complex<float> Ifa2 = -(std::complex<float>{1, 0} - Z1(faultNode - 1, faultNode - 1) * Ifa1) / Z2(faultNode2 - 1, faultNode2 - 1);
	std::complex<float> Ifa0 = -(std::complex<float>{1, 0} - Z1(faultNode - 1, faultNode - 1) * Ifa1) / (Z0(faultNode2 - 1, faultNode2 - 1) + Zf * std::complex<float>{0, 3});

	If << Ifa1, Ifa2, Ifa0;
	Eigen::VectorXcf If_abc = St * If;

	fprintf(stdout, "���ϵ�abc���·����: \n");
	for (int i = 0; i < If_abc.size(); i++)
		fprintf(stdout, "%f\n", abs(If_abc(i)));

	auto Is = If_abc(1) + If_abc(2);
	std::cout << "���ϵ���: " << Is << std::endl;

	this->getBusVoltageAndCurrent(If, faultNode);
}

void Grid::getBusVoltageAndCurrent(const Eigen::VectorXcf& If, const size_t faultNode)
{
	const auto branchNum = this->Y1.rows();
	Eigen::MatrixXcf U1 = Eigen::MatrixXcf(branchNum, 1);
	auto U2{ U1 }, U0{ U1 };
	Eigen::MatrixXf Uabc = Eigen::MatrixXf(branchNum, 3);

#pragma omp parallel for
	for (size_t i = 0; i < branchNum; i++) {
		U1(i) = std::complex<float>{ 1,0 } - this->Z1(faultNode - 1, i) * If(0);
		U2(i) = -(this->Z2(faultNode - 1, i) * If(1));
		U0(i) = -(this->Z0(faultNode - 1, i) * If(2));
		Uabc.row(i) = (St * ((Eigen::MatrixXcf(3, 1) << U1(i), U2(i), U0(i)).finished())).cwiseAbs();
	}

	for (size_t i = 0; i < branchNum; i++)
		std::cout << "�ڵ�" << i + 1 << "��abc���·��ѹ: " << Uabc.row(i) << std::endl;

	std::map < std::pair<int, int>, std::tuple<std::complex<float>, std::complex<float>, std::complex<float>>> Isx;

#pragma omp parallel for
	for (size_t i = 0; i < branchNum; i++) {
		for (decltype(i) j = 1; j < i; j++) {
			std::pair<int, int> curSocket = std::make_pair(i, j);
			std::complex<float> Is1{ 0,0 }, Is2{ 0,0 }, Is0{ 0,0 };
			if (abs(this->Y1(i, j)) > epsilon)	Is1 = -(U1(i) - U1(j)) * this->Y1(i, j);
			if (abs(this->Y2(i, j)) > epsilon)	Is2 = -(U2(i) - U2(j)) * this->Y2(i, j);
			if (abs(this->Y0(i, j)) > epsilon)	Is0 = -(U0(i) - U0(j)) * this->Y0(i, j);
			Isx.insert({ curSocket,{Is1,Is2,Is0} });
		}
	}

	for (const auto& iter : Isx) {
		auto& thisTuple = iter.second;
		Eigen::MatrixXcf I = St * ((Eigen::MatrixXcf(3, 1) << std::get<0>(thisTuple), std::get<1>(thisTuple), std::get<2>(thisTuple)).finished());
		std::cout << "�ڵ�" << iter.first.first << "��ڵ�" << iter.first.second << "���abc���·����: " << std::endl;
		std::cout << I << std::endl;
	}

}

Grid::Grid(const Eigen::MatrixXcf& Y1_, const Eigen::MatrixXcf& Y2_, const Eigen::MatrixXcf& Y0_) : myPool(maxThreadsNum), Y1(Y1_), Y2(Y2_), Y0(Y0_)
{
	omp_set_num_threads(std::thread::hardware_concurrency());
	std::vector<std::function<void()>> Y_Z_interConvert;
	Y_Z_interConvert.reserve(3);

	Y_Z_interConvert.emplace_back([this]()->void { this->Z1 = this->Y1.inverse(); });
	Y_Z_interConvert.emplace_back([this]()->void { this->Z2 = this->Y2.inverse(); });
	Y_Z_interConvert.emplace_back([this]()->void { this->Z0 = this->Y0.inverse(); });

	for (const auto& iter : Y_Z_interConvert)
		this->myPool.enqueue(iter);

}

Grid::~Grid()
{
}


