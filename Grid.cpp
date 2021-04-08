#include "Grid.h"
#include <iostream>

Grid::Grid(const ThreeSequenceData& data, const std::vector<Transformer2>& transformList) : myPool(maxThreadsNum)
{
    omp_set_num_threads(std::thread::hardware_concurrency());

    this->NodesNum = std::max(data.lineData1.col(0).maxCoeff(), data.lineData1.col(1).maxCoeff());

    this->Y1 = this->setYxFromSheet(data.lineData1, this->SocketData1);
    this->Y2 = this->setYxFromSheet(data.lineData2, this->SocketData2);
    this->Y0 = this->setYxFromSheet(data.lineData0, this->SocketData0);

    this->adjustTransformerRatio(transformList);

    this->Z1 = this->Y1.inverse();
    this->Z2 = this->Y2.inverse();
    this->Z0 = this->Y0.inverse();
}

Grid::~Grid()
{
}

//ͨ�����ݼ���õ��ɾ���
Eigen::MatrixXcf Grid::setYxFromSheet(const Eigen::MatrixXf &line_data_sheet, std::map<std::pair<size_t, size_t>, std::pair<cf, cf>> &socketData)
{
    //dataSheet: https://img-blog.csdn.net/20180616204918263?watermark/2/text/aHR0cHM6Ly9ibG9nLmNzZG4ubmV0L3FxXzMyNDEyNzU5/font/5a6L5L2T/fontsize/400/fill/I0JBQkFCMA==/dissolve/70
    //����û�������ѹ�������һ�У���Ϊ�󲿷���·�ı�ȶ���1

    const Eigen::VectorXi inNode = line_data_sheet.col(0).cast<int>().array();
    const Eigen::VectorXi outNode = line_data_sheet.col(1).cast<int>().array();
    const Eigen::VectorXcf Z = line_data_sheet.col(2) + Imaginer * line_data_sheet.col(3);
    const Eigen::VectorXcf B = line_data_sheet.col(4) * Imaginer;

    auto branchNum = inNode.size();
    const auto busNum = std::max(inNode.maxCoeff(), outNode.maxCoeff());

    Eigen::VectorXcf y = Eigen::VectorXcf(branchNum).setOnes().array() / Z.array();
    Eigen::MatrixXcf Y = Eigen::MatrixXcf::Zero(busNum, busNum);


//#pragma omp parallel for
//����̶߳�ͬһ������ͬʱд���̲߳���ȫ
    for (decltype(branchNum) i = 0; i < branchNum; i++)
    {
        const size_t from = inNode(i) - 1, to = outNode(i) - 1;

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

//TODO:�Գƶ�·����
void Grid::SymmetricShortCircuit(const Eigen::MatrixXcf& AdmittanceMatrix, const size_t shortPoint, const float Sb, const float Uav)
{
    const auto busNum = AdmittanceMatrix.rows();
    Eigen::VectorXf Is = Eigen::VectorXf::Zero(busNum, 1);
}

//�����·����
void Grid::lgShortCircuit(const size_t faultNode, const cf& Zf)
{
    const auto faultNode2 = faultNode, faultNode0 = faultNode;

    const cf Ifa1 = cf{ 1, 0 } / (Z1(faultNode - 1, faultNode - 1) + Z2(faultNode2 - 1, faultNode2 - 1) + Z0(faultNode0 - 1, faultNode0 - 1) + Zf * cf{ 0, 3 });
    Eigen::VectorXcf If(3, 1);
    //auto Ifa2{ Ifa1 }, Ifa0{ Ifa1 };
    //If << Ifa1, Ifa2, Ifa0;
    If << Ifa1, Ifa1, Ifa1;
    const Eigen::VectorXcf If_abc = St * If;

    fprintf(stdout, "���ϵ�abc���·����: \n");
    for (decltype(If_abc.size()) i = 0; i < If_abc.size(); i++)
        fprintf(stdout, "%f\n", abs(If_abc(i)));

    const auto Is = If_abc(0);
    std::cout << "���ϵ���: " << Is << std::endl;

    this->getBusVoltageAndCurrent(If, faultNode);
}

//�����·����
void Grid::llShortCircuit(const size_t faultNode, const cf& Zf)
{
    const auto faultNode2 = faultNode;

    Eigen::VectorXcf If(3, 1);
    const cf Ifa1 = cf{ 1, 0 } / (Z1(faultNode - 1, faultNode - 1) + Z2(faultNode2 - 1, faultNode2 - 1) + Zf * cf{ 0, 1 });
    const cf Ifa2{ -Ifa1.real(), -Ifa1.imag() }, Ifa0{ 0, 0 };
    If << Ifa1, Ifa2, Ifa0;
    const Eigen::VectorXcf If_abc = St * If;

    fprintf(stdout, "���ϵ�abc���·����: \n");
    for (decltype(If_abc.size()) i = 0; i < If_abc.size(); i++)
        fprintf(stdout, "%f\n", abs(If_abc(i)));

    const auto Is = If_abc(2);
    std::cout << "���ϵ���: " << Is << std::endl;

    this->getBusVoltageAndCurrent(If, faultNode);
}

//����Եض�·����
void Grid::llgShortCircuit(const size_t faultNode, const cf& Zf)
{
    const auto faultNode2 = faultNode, faultNode0 = faultNode;

    Eigen::VectorXcf If(3, 1);
    const cf Ifa1 = cf{ 1, 0 } / (Z1(faultNode - 1, faultNode - 1) + Z2(faultNode2 - 1, faultNode2 - 1) * (Z0(faultNode0 - 1, faultNode0 - 1) + Zf * cf{ 0, 3 }) / (Z2(faultNode2 - 1, faultNode2 - 1) + Z0(faultNode0 - 1, faultNode0 - 1) + Zf * cf{ 0, 3 }));
    const cf Ifa2 = -(cf{ 1, 0 } - Z1(faultNode - 1, faultNode - 1) * Ifa1) / Z2(faultNode2 - 1, faultNode2 - 1);
    const cf Ifa0 = -(cf{ 1, 0 } - Z1(faultNode - 1, faultNode - 1) * Ifa1) / (Z0(faultNode2 - 1, faultNode2 - 1) + Zf * cf{ 0, 3 });

    If << Ifa1, Ifa2, Ifa0;
    const Eigen::VectorXcf If_abc = St * If;

    fprintf(stdout, "���ϵ�abc���·����: \n");
    for (int i = 0; i < If_abc.size(); i++)
        fprintf(stdout, "%f\n", abs(If_abc(i)));

    const auto Is = If_abc(1) + If_abc(2);
    std::cout << "���ϵ���: " << Is << std::endl;

    this->getBusVoltageAndCurrent(If, faultNode);
}

void Grid::getBusVoltageAndCurrent(const Eigen::VectorXcf& If, const size_t faultNode)
{
    auto branchNum = this->Y1.rows();
    Eigen::MatrixXcf U1 = Eigen::MatrixXcf(branchNum, 1);
    auto U2{ U1 }, U0{ U1 };
    Eigen::MatrixXf Uabc = Eigen::MatrixXf(branchNum, 3);

    #pragma omp parallel for
    for (decltype(branchNum) i = 0; i < branchNum; i++)
    {
        U1(i) = cf{ 1,0 } - this->Z1(faultNode - 1, i) * If(0);
        U2(i) = -(this->Z2(faultNode - 1, i) * If(1));
        U0(i) = -(this->Z0(faultNode - 1, i) * If(2));
        Uabc.row(i) = (St * ((Eigen::MatrixXcf(3, 1) << U1(i), U2(i), U0(i)).finished())).cwiseAbs();
    }

    for (decltype(branchNum) i = 0; i < branchNum; i++)
        std::cout << "�ڵ�" << i + 1 << "��abc���·��ѹ: " << Uabc.row(i) << std::endl;

    std::map < std::pair<int, int>, std::tuple<cf, cf, cf>> Isx;

//#pragma omp parallel for
    for (decltype(branchNum) i = 0; i < branchNum; i++)
    {
        for (decltype(i) j = 1; j < i; j++)
        {
            std::pair<int, int> curSocket = std::make_pair(i, j);
            cf Is1{ 0,0 }, Is2{ 0,0 }, Is0{ 0,0 };
            if (abs(this->Y1(i, j)) > epsilon)  Is1 = -(U1(i) - U1(j)) * this->Y1(i, j);
            if (abs(this->Y2(i, j)) > epsilon)	Is2 = -(U2(i) - U2(j)) * this->Y2(i, j);
            if (abs(this->Y0(i, j)) > epsilon)	Is0 = -(U0(i) - U0(j)) * this->Y0(i, j);
            Isx.insert({ curSocket,{Is1,Is2,Is0} });
        }
    }

    for (const auto &iter : Isx)
    {
        const auto& thisTuple = iter.second;
        Eigen::MatrixXcf I = St * ((Eigen::MatrixXcf(3, 1) << std::get<0>(thisTuple), std::get<1>(thisTuple), std::get<2>(thisTuple)).finished());
        std::cout << "�ڵ�" << iter.first.first << "��ڵ�" << iter.first.second << "���abc���·����: " << std::endl;
        std::cout << I << std::endl;
    }
}

void Grid::adjustTransformerRatio(const std::vector<Transformer2>& transList)
{
    #pragma omp parallel for
    for (decltype(transList.size()) i = 0; i < transList.size(); i++)
    {
        const auto &transformer = transList.at(i);
        const auto primaryNode = transformer.getPrimaryNode() - 1;
        const auto secondNode = transformer.getSecondaryNode() - 1;
        const auto k = transformer.getRatio();

        const auto &y1 = this->SocketData1.find({ transformer.getPrimaryNode(), transformer.getSecondaryNode() })->second.first;
        const auto &y2 = this->SocketData2.find({ transformer.getPrimaryNode(), transformer.getSecondaryNode() })->second.first;
        const auto &y0 = this->SocketData0.find({ transformer.getPrimaryNode(), transformer.getSecondaryNode() })->second.first;

        const cf delta_Y1_tt = y1 / (k * k) - y1;
        const cf delta_Y1_ft = y1 - y1 / k;

        const cf delta_Y2_tt = y2 / (k * k) - y2;
        const cf delta_Y2_ft = y2 - y2 / k;

        const cf delta_Y0_tt = y0 / (k * k) - y0;
        const cf delta_Y0_ft = y0 - y0 / k;

        this->Y1(secondNode, secondNode) += delta_Y1_tt;
        this->Y1(primaryNode, secondNode) += delta_Y1_ft;
        this->Y1(secondNode, primaryNode) = this->Y1(primaryNode, secondNode);

        this->Y2(secondNode, secondNode) += delta_Y2_tt;
        this->Y2(primaryNode, secondNode) += delta_Y2_ft;
        this->Y2(secondNode, primaryNode) = this->Y2(primaryNode, secondNode);

        this->Y0(secondNode, secondNode) += delta_Y0_tt;
        this->Y0(primaryNode, secondNode) += delta_Y0_ft;
        this->Y0(secondNode, primaryNode) = this->Y0(primaryNode, secondNode);
    }
}

Eigen::MatrixXcf Grid::getYx(const SEQUENCE whichSeq) const
{
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

Eigen::MatrixXcf Grid::getZx(const SEQUENCE whichSeq) const
{
    switch (whichSeq)
    {
    case SEQUENCE::POSITIVE:
        return this->Z1;
    case SEQUENCE::NEGATIVE:
        return this->Z2;
    case SEQUENCE::ZERO:
        return this->Z0;
    default:
        break;
    }
    return Eigen::MatrixXcf(0, 0);
}
