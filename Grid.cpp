#include "Grid.h"
#include <iostream>

Grid::Grid(const ThreeSequenceData& data, const std::list<Transformer2>& transformList) : myPool(maxThreadsNum)
{
	omp_set_num_threads(std::thread::hardware_concurrency());

	/*
	std::vector<std::function<void()>> RawLineWithoutGeneratorAndTransformer;
	RawLineWithoutGeneratorAndTransformer.reserve(3);
	RawLineWithoutGeneratorAndTransformer.emplace_back([this, &Z1_data_sheet]()->void { this->Y1 = Grid::getAdmittanceMatrixBySheet(Z1_data_sheet); });
	RawLineWithoutGeneratorAndTransformer.emplace_back([this, &Z2_data_sheet]()->void { this->Y2 = Grid::getAdmittanceMatrixBySheet(Z2_data_sheet); });
	RawLineWithoutGeneratorAndTransformer.emplace_back([this, &Z0_data_sheet]()->void { this->Y0 = Grid::getAdmittanceMatrixBySheet(Z0_data_sheet); });

	for (const auto& iter : RawLineWithoutGeneratorAndTransformer)
		this->myPool.enqueue(iter);
	*/

	this->NodesNum = std::max(data.lineData1.col(0).maxCoeff(), data.lineData1.col(1).maxCoeff());

	this->Y1 = this->getAdmittanceMatrixBySheet(data.lineData1, this->SocketData1);
	this->Y2 = this->getAdmittanceMatrixBySheet(data.lineData2, this->SocketData2);
	this->Y0 = this->getAdmittanceMatrixBySheet(data.lineData0, this->SocketData0);

	this->adjustTransformerRatio(transformList);

	this->Z1 = this->Y1.inverse();
	this->Z2 = this->Y2.inverse();
	this->Z0 = this->Y0.inverse();
}

Grid::~Grid()
{
}

//通过数据集获得导纳矩阵
Eigen::MatrixXcf Grid::getAdmittanceMatrixBySheet(const Eigen::MatrixXf& line_data_sheet, std::map<std::pair<size_t, size_t>, std::pair<cf, cf>>& socketData)
{
	//dataSheet: https://img-blog.csdn.net/20180616204918263?watermark/2/text/aHR0cHM6Ly9ibG9nLmNzZG4ubmV0L3FxXzMyNDEyNzU5/font/5a6L5L2T/fontsize/400/fill/I0JBQkFCMA==/dissolve/70
	//但是没有上面变压器变比那一列，因为大部分线路的变比都是1

	const Eigen::VectorXi inNode = line_data_sheet.col(0).cast<int>().array();
	const Eigen::VectorXi outNode = line_data_sheet.col(1).cast<int>().array();
	const Eigen::VectorXcf Z = line_data_sheet.col(2) + Imaginer * line_data_sheet.col(3);
	const Eigen::VectorXcf B = line_data_sheet.col(4) * Imaginer;

	auto branchNum = inNode.size();
	auto busNum = std::max(inNode.maxCoeff(), outNode.maxCoeff());

	Eigen::VectorXcf y = Eigen::VectorXcf(branchNum).setOnes().array() / Z.array();
	Eigen::MatrixXcf Y = Eigen::MatrixXcf::Zero(busNum, busNum);


#pragma omp parallel for
	for (decltype(branchNum) i = 0; i < branchNum; i++)
	{
		const size_t from = inNode(i) - 1, to = outNode(i) - 1;
		/*
		Y(from, to) -= y(i) / K(i);
		Y(to, from) = Y(from, to);
		Y(from, from) += y(i) + (B(i) / cf{2, 0});
		Y(to, to) += y(i) / (K(i) * K(i)) + (B(i) / cf{2, 0});
		*/

		Y(from, to) -= y(i);
		Y(to, from) = Y(from, to);
		Y(from, from) += y(i) + (B(i) / cf{ 2, 0 });
		Y(to, to) += y(i) + (B(i) / cf{ 2, 0 });

		const auto line_data = std::make_pair(y(i), B(i));
		socketData.insert({ {inNode(i), outNode(i)}, line_data });
		socketData.insert({ {outNode(i), inNode(i)}, line_data });
	}

	return Y;
}
/*
//通过阻抗矩阵Zdx获得导纳矩阵
Eigen::MatrixXcf Grid::getAdmittanceMatrixFromReactanceMatrix(const Eigen::MatrixXcf& Zdx)
{
	//dataSheet: https://blog.csdn.net/qq_32412759/article/details/80715617
	Eigen::VectorXi inNode = Zdx.col(0).real().cast<int>().array();
	Eigen::VectorXi outNode = Zdx.col(1).real().cast<int>().array();
	Eigen::MatrixXcf yd = Eigen::VectorXcf(inNode.size()).setOnes().array() / Zdx.col(2).array();
	auto busNum = std::max(inNode.maxCoeff(), outNode.maxCoeff());
	Eigen::MatrixXcf Y = Eigen::MatrixXcf::Zero(busNum, busNum);

	//这里采用map，因为当pair作为unordered_map的key时，需要hash函数，而map不需要
	std::map<std::pair<int, int>, cf> socket_yd{};
	std::unordered_map<int, std::unordered_set<int>> DirectConnectNodePair{};

#pragma omp parallel for
	for (decltype(inNode.size()) i = 0; i < inNode.size(); i++)
	{
		//其实可以省一半空间的，用空间换时间吧，反正也是对称矩阵嘛
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
				//TODO:这里到底是加还是减
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
*/

//TODO:对称短路计算
void Grid::SymmetricShortCircuit(const Eigen::MatrixXcf& AdmittanceMatrix, const size_t shortPoint, const float Sb, const float Uav)
{
	const auto busNum = AdmittanceMatrix.rows();
	Eigen::VectorXf Is = Eigen::VectorXf::Zero(busNum, 1);
}

//单相短路计算
void Grid::lgShortCircuit(const size_t faultNode, const cf& Zf)
{
	auto faultNode2 = faultNode, faultNode0 = faultNode;

	cf Ifa1 = cf{ 1, 0 } / (Z1(faultNode - 1, faultNode - 1) + Z2(faultNode2 - 1, faultNode2 - 1) + Z0(faultNode0 - 1, faultNode0 - 1) + Zf * cf{ 0, 3 });
	Eigen::VectorXcf If(3, 1);
	//auto Ifa2{ Ifa1 }, Ifa0{ Ifa1 };
	//If << Ifa1, Ifa2, Ifa0;
	If << Ifa1, Ifa1, Ifa1;
	Eigen::VectorXcf If_abc = St * If;

	fprintf(stdout, "故障点abc相短路电流: \n");
	for (int i = 0; i < If_abc.size(); i++)
		fprintf(stdout, "%f\n", abs(If_abc(i)));

	auto Is = If_abc(0);
	std::cout << "故障电流: " << Is << std::endl;

	this->getBusVoltageAndCurrent(If, faultNode);
}

//两相短路计算
void Grid::llShortCircuit(const size_t faultNode, const cf& Zf)
{
	auto faultNode2 = faultNode;

	Eigen::VectorXcf If(3, 1);
	cf Ifa1 = cf{ 1, 0 } / (Z1(faultNode - 1, faultNode - 1) + Z2(faultNode2 - 1, faultNode2 - 1) + Zf * cf{ 0, 1 });
	cf Ifa2{ -Ifa1.real(), -Ifa1.imag() }, Ifa0{ 0, 0 };
	If << Ifa1, Ifa2, Ifa0;
	Eigen::VectorXcf If_abc = St * If;

	fprintf(stdout, "故障点abc相短路电流: \n");
	for (int i = 0; i < If_abc.size(); i++)
		fprintf(stdout, "%f\n", abs(If_abc(i)));

	auto Is = If_abc(2);
	std::cout << "故障电流: " << Is << std::endl;

	this->getBusVoltageAndCurrent(If, faultNode);
}

//两相对地短路计算
void Grid::llgShortCircuit(const size_t faultNode, const cf& Zf)
{
	auto faultNode2 = faultNode, faultNode0 = faultNode;

	Eigen::VectorXcf If(3, 1);
	cf Ifa1 = cf{ 1, 0 } / (Z1(faultNode - 1, faultNode - 1) + Z2(faultNode2 - 1, faultNode2 - 1) * (Z0(faultNode0 - 1, faultNode0 - 1) + Zf * cf{ 0, 3 }) / (Z2(faultNode2 - 1, faultNode2 - 1) + Z0(faultNode0 - 1, faultNode0 - 1) + Zf * cf{ 0, 3 }));
	cf Ifa2 = -(cf{ 1, 0 } - Z1(faultNode - 1, faultNode - 1) * Ifa1) / Z2(faultNode2 - 1, faultNode2 - 1);
	cf Ifa0 = -(cf{ 1, 0 } - Z1(faultNode - 1, faultNode - 1) * Ifa1) / (Z0(faultNode2 - 1, faultNode2 - 1) + Zf * cf{ 0, 3 });

	If << Ifa1, Ifa2, Ifa0;
	Eigen::VectorXcf If_abc = St * If;

	fprintf(stdout, "故障点abc相短路电流: \n");
	for (int i = 0; i < If_abc.size(); i++)
		fprintf(stdout, "%f\n", abs(If_abc(i)));

	auto Is = If_abc(1) + If_abc(2);
	std::cout << "故障电流: " << Is << std::endl;

	this->getBusVoltageAndCurrent(If, faultNode);
}

void Grid::getBusVoltageAndCurrent(const Eigen::VectorXcf& If, const size_t faultNode)
{
	auto branchNum = this->Y1.rows();
	Eigen::MatrixXcf U1 = Eigen::MatrixXcf(branchNum, 1);
	auto U2{ U1 }, U0{ U1 };
	Eigen::MatrixXf Uabc = Eigen::MatrixXf(branchNum, 3);

#pragma omp parallel for
	for (decltype(branchNum) i = 0; i < branchNum; i++) {
		U1(i) = cf{ 1,0 } - this->Z1(faultNode - 1, i) * If(0);
		U2(i) = -(this->Z2(faultNode - 1, i) * If(1));
		U0(i) = -(this->Z0(faultNode - 1, i) * If(2));
		Uabc.row(i) = (St * ((Eigen::MatrixXcf(3, 1) << U1(i), U2(i), U0(i)).finished())).cwiseAbs();
	}

	for (size_t i = 0; i < branchNum; i++)
		std::cout << "节点" << i + 1 << "的abc相短路电压: " << Uabc.row(i) << std::endl;

	std::map < std::pair<int, int>, std::tuple<cf, cf, cf>> Isx;

#pragma omp parallel for
	for (size_t i = 0; i < branchNum; i++) {
		for (decltype(i) j = 1; j < i; j++) {
			std::pair<int, int> curSocket = std::make_pair(i, j);
			cf Is1{ 0,0 }, Is2{ 0,0 }, Is0{ 0,0 };
			if (abs(this->Y1(i, j)) > epsilon)  Is1 = -(U1(i) - U1(j)) * this->Y1(i, j);
			if (abs(this->Y2(i, j)) > epsilon)	Is2 = -(U2(i) - U2(j)) * this->Y2(i, j);
			if (abs(this->Y0(i, j)) > epsilon)	Is0 = -(U0(i) - U0(j)) * this->Y0(i, j);
			Isx.insert({ curSocket,{Is1,Is2,Is0} });
		}
	}

	for (const auto& iter : Isx) {
		auto& thisTuple = iter.second;
		Eigen::MatrixXcf I = St * ((Eigen::MatrixXcf(3, 1) << std::get<0>(thisTuple), std::get<1>(thisTuple), std::get<2>(thisTuple)).finished());
		std::cout << "节点" << iter.first.first << "与节点" << iter.first.second << "间的abc相短路电流: " << std::endl;
		std::cout << I << std::endl;
	}

}

void Grid::adjustTransformerRatio(const std::list<Transformer2>& transList)
{

#pragma omp parallel for
	for (const auto& transformer : transList) {
		const auto primaryNode = transformer.getPrimaryNode();
		const auto secondNode = transformer.getSecondaryNode();

		const auto y = this->SocketData1.find({ primaryNode ,secondNode })->second.first;

		const cf delta_Y_tt = (1 / (transformer.getRatio() * transformer.getRatio()) - 1) * y;
		const cf delta_Y_ft = y * (1 - 1 / transformer.getRatio());

		this->Y1(secondNode, secondNode) += delta_Y_tt;
		this->Y1(primaryNode, secondNode) += delta_Y_ft;
		this->Y1(secondNode, primaryNode) = this->Y1(primaryNode, secondNode);

		this->Y2(secondNode, secondNode) += delta_Y_tt;
		this->Y2(primaryNode, secondNode) += delta_Y_ft;
		this->Y2(secondNode, primaryNode) = this->Y2(primaryNode, secondNode);

		this->Y0(secondNode, secondNode) += delta_Y_tt;
		this->Y0(primaryNode, secondNode) += delta_Y_ft;
		this->Y0(secondNode, primaryNode) = this->Y0(primaryNode, secondNode);
	}
}

Eigen::MatrixXcf Grid::getAdmittanceMatrix(const SEQUENCE whichSeq){
	switch (whichSeq)
	{
	case SEQUENCE::POSITIVE:
		return this->Y1;
	case SEQUENCE::NEGATIVE:
		return this->Y2;
	case SEQUENCE::ZERO:
		return this->Y0;
	default:
		break;
	}
	return Eigen::MatrixXcf(0, 0);
}