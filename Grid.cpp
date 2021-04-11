#include "Grid.h"
#include <iostream>

Grid::Grid(const NodeType node_, const ThreeSequenceData &data, const std::vector<std::pair<IdealTransformer2, cf>> &idealTransformList, const std::vector<Transformer2> &transList, const std::vector<Generator> &geneList) : NodeNum(node_), myPool(maxThreadsNum)
{
    omp_set_num_threads(std::thread::hardware_concurrency());

    this->Y2 = this->Y0 = this->Y1 = Eigen::MatrixXcf::Zero(this->NodeNum, this->NodeNum);

    this->setYxFromSheet(data.lineData1, this->Y1, this->SocketData1);
    this->setYxFromSheet(data.lineData2, this->Y2, this->SocketData2);
    this->setYxFromSheet(data.lineData0, this->Y0, this->SocketData0);

    this->mountIdealTransformer2_PrimarySideReactance(idealTransformList);
    this->mountTransformer2(transList);
    this->mountGenerator(geneList);

    this->Z1 = this->Y1.inverse();
    /*this->Z2 = this->Y2.inverse();
    this->Z0 = this->Y0.inverse();*/
}

Grid::~Grid()
{
}

//通过数据集获得导纳矩阵
void Grid::setYxFromSheet(const Eigen::MatrixXf &line_data_sheet, Eigen::MatrixXcf &Y, std::map<socketType, std::pair<cf, cf>> &socketData)
{
    //dataSheet: https://img-blog.csdn.net/20180616204918263?watermark/2/text/aHR0cHM6Ly9ibG9nLmNzZG4ubmV0L3FxXzMyNDEyNzU5/font/5a6L5L2T/fontsize/400/fill/I0JBQkFCMA==/dissolve/70
    //但是没有上面变压器变比那一列，因为大部分线路的变比都是1

    const Eigen::VectorXi inNode = line_data_sheet.col(0).cast<int>().array();
    const Eigen::VectorXi outNode = line_data_sheet.col(1).cast<int>().array();
    const Eigen::VectorXcf Z = line_data_sheet.col(2) + Imaginer * line_data_sheet.col(3);
    const Eigen::VectorXcf B = line_data_sheet.col(4) * Imaginer;

    auto branchNum = inNode.size();
    Eigen::VectorXcf y = Eigen::VectorXcf(branchNum).setOnes().array() / Z.array();

    //#pragma omp parallel for
    //多个线程对同一个容器同时写，线程不安全
    for (decltype(branchNum) i = 0; i < branchNum; i++)
    {
        //TODO: 应该完善直接接地的电抗
        //如果两点最小的那个为零，说明这是一个节点接地电抗
        if (!std::min(inNode(i), outNode(i)))
        {
            const auto insertedNode = std::max(inNode(i), outNode(i)) - 1;
            Y(insertedNode, insertedNode) += y(i);
            continue;
        }
        const NodeType from = inNode(i) - 1, to = outNode(i) - 1;

        Y(from, to) -= y(i);
        Y(to, from) = Y(from, to);
        Y(from, from) += y(i) + (B(i) / cf{2, 0});
        Y(to, to) += y(i) + (B(i) / cf{2, 0});

        const auto line_data_y_B = std::make_pair(y(i), B(i));
        socketData.insert(
            {{{inNode(i), outNode(i)}, line_data_y_B},
             {{outNode(i), inNode(i)}, line_data_y_B}});
    }
}

//TODO:对称短路计算
std::tuple<Eigen::VectorXcf, std::list<std::pair<socketType, cf>>> Grid::SymmetricShortCircuit(const NodeType shortPoint, const cf &Zf, const DeviceArgType Uav)
{
    const Eigen::VectorXcf &Z = this->Z1.col(shortPoint - 1);
    Eigen::VectorXcf Ui = Eigen::VectorXcf(this->NodeNum).setOnes();
    const auto If = cf(Uav, 0) / (Z(shortPoint - 1) + Zf);
    Ui = Ui.array() - Z.array() * If;
    Ui(shortPoint - 1) = cf{0.0f, 0.0f};

    std::set<socketType> pqNodeSum{};
    std::list<std::pair<socketType, cf>> shortCurrent{};

    //TODO:这里把变压器当成了k=1的变压器，或者说根本没有变压器
    const auto &sockData = this->SocketData1;
    for (const auto &iter : sockData)
    {
        const auto &socket = iter.first;
        //const auto sum = socket.first + socket.second;
        if (pqNodeSum.find(socket) != pqNodeSum.end())
            continue;

        const auto &y = iter.second.first;

        pqNodeSum.insert(
            {{socket.second, socket.first},
             socket});

        const auto Ipq = (Ui(socket.first - 1) - Ui(socket.second - 1)) * y;
        shortCurrent.emplace_back(std::make_pair(socket, Ipq));
    }
    return std::make_tuple(Ui, shortCurrent);
}

//单相短路计算
void Grid::lgShortCircuit(const NodeType faultNode, const cf &Zf)
{
    const auto faultNode2 = faultNode, faultNode0 = faultNode;

    const cf Ifa1 = cf{1, 0} / (Z1(faultNode - 1, faultNode - 1) + Z2(faultNode2 - 1, faultNode2 - 1) + Z0(faultNode0 - 1, faultNode0 - 1) + Zf * cf{3, 0});
    Eigen::VectorXcf If(3, 1);
    //auto Ifa2{ Ifa1 }, Ifa0{ Ifa1 };
    //If << Ifa1, Ifa2, Ifa0;
    If << Ifa1, Ifa1, Ifa1;
    const Eigen::VectorXcf If_abc = St * If;

    fprintf(stdout, "故障点abc相短路电流: \n");
    for (decltype(If_abc.size()) i = 0; i < If_abc.size(); i++)
        fprintf(stdout, "%f\n", abs(If_abc(i)));

    const auto Is = If_abc(0);
    std::cout << "故障电流: " << Is << std::endl;

    this->getBusVoltageAndCurrent(If, faultNode);
}

//两相短路计算
void Grid::llShortCircuit(const NodeType faultNode, const cf &Zf)
{
    const auto faultNode2 = faultNode;

    Eigen::VectorXcf If(3, 1);
    const cf Ifa1 = cf{1, 0} / (Z1(faultNode - 1, faultNode - 1) + Z2(faultNode2 - 1, faultNode2 - 1) + Zf);
    const cf Ifa2{-Ifa1.real(), -Ifa1.imag()}, Ifa0{0, 0};
    If << Ifa1, Ifa2, Ifa0;
    const Eigen::VectorXcf If_abc = St * If;

    fprintf(stdout, "故障点abc相短路电流: \n");
    for (decltype(If_abc.size()) i = 0; i < If_abc.size(); i++)
        fprintf(stdout, "%f\n", abs(If_abc(i)));

    const auto Is = If_abc(2);
    std::cout << "故障电流: " << Is << std::endl;

    this->getBusVoltageAndCurrent(If, faultNode);
}

//两相对地短路计算
void Grid::llgShortCircuit(const NodeType faultNode, const cf &Zf)
{
    const auto faultNode2 = faultNode, faultNode0 = faultNode;

    Eigen::VectorXcf If(3, 1);
    const cf Ifa1 = cf{1, 0} / (Z1(faultNode - 1, faultNode - 1) + Z2(faultNode2 - 1, faultNode2 - 1) * (Z0(faultNode0 - 1, faultNode0 - 1) + Zf * cf{3, 0}) / (Z2(faultNode2 - 1, faultNode2 - 1) + Z0(faultNode0 - 1, faultNode0 - 1) + Zf * cf{3, 0}));
    const cf Ifa2 = -(cf{1, 0} - Z1(faultNode - 1, faultNode - 1) * Ifa1) / Z2(faultNode2 - 1, faultNode2 - 1);
    const cf Ifa0 = -(cf{1, 0} - Z1(faultNode - 1, faultNode - 1) * Ifa1) / (Z0(faultNode2 - 1, faultNode2 - 1) + Zf * cf{0, 3});

    If << Ifa1, Ifa2, Ifa0;
    const Eigen::VectorXcf If_abc = St * If;

    fprintf(stdout, "故障点abc相短路电流: \n");
    for (int i = 0; i < If_abc.size(); i++)
        fprintf(stdout, "%f\n", abs(If_abc(i)));

    const auto Is = If_abc(1) + If_abc(2);
    std::cout << "故障电流: " << Is << std::endl;

    this->getBusVoltageAndCurrent(If, faultNode);
}

void Grid::getBusVoltageAndCurrent(const Eigen::VectorXcf &If, const NodeType faultNode)
{
    auto branchNum = this->Y1.rows();
    Eigen::MatrixXcf U1 = Eigen::MatrixXcf(branchNum, 1);
    auto U2{U1}, U0{U1};
    Eigen::MatrixXf Uabc = Eigen::MatrixXf(branchNum, 3);

#pragma omp parallel for
    for (decltype(branchNum) i = 0; i < branchNum; i++)
    {
        U1(i) = cf{1, 0} - this->Z1(faultNode - 1, i) * If(0);
        U2(i) = -(this->Z2(faultNode - 1, i) * If(1));
        U0(i) = -(this->Z0(faultNode - 1, i) * If(2));
        Uabc.row(i) = (St * ((Eigen::MatrixXcf(3, 1) << U1(i), U2(i), U0(i)).finished())).cwiseAbs();
    }

    for (decltype(branchNum) i = 0; i < branchNum; i++)
        std::cout << "节点" << i + 1 << "的abc相短路电压: " << Uabc.row(i) << std::endl;

    std::map<std::pair<int, int>, std::tuple<cf, cf, cf>> Isx{};

    //#pragma omp parallel for
    for (decltype(branchNum) i = 0; i < branchNum; i++)
    {
        for (decltype(i) j = 1; j < i; j++)
        {
            std::pair<int, int> curSocket = std::make_pair(i, j);
            cf Is1{0, 0}, Is2{0, 0}, Is0{0, 0};
            if (std::abs(this->Y1(i, j)) > epsilon)
                Is1 = -(U1(i) - U1(j)) * this->Y1(i, j);
            if (std::abs(this->Y2(i, j)) > epsilon)
                Is2 = -(U2(i) - U2(j)) * this->Y2(i, j);
            if (std::abs(this->Y0(i, j)) > epsilon)
                Is0 = -(U0(i) - U0(j)) * this->Y0(i, j);
            Isx.insert({curSocket, {Is1, Is2, Is0}});
        }
    }

    for (const auto &iter : Isx)
    {
        const auto &thisTuple = iter.second;
        Eigen::MatrixXcf I = St * ((Eigen::MatrixXcf(3, 1) << std::get<0>(thisTuple), std::get<1>(thisTuple), std::get<2>(thisTuple)).finished());
        std::cout << "节点" << iter.first.first << "与节点" << iter.first.second << "间的abc相短路电流: " << std::endl;
        std::cout << I << std::endl;
    }
}

void Grid::adjustIdealTransformer2_PrimarySideReactanceRatio(const std::vector<Transformer2> &transList)
{
#pragma omp parallel for
    for (decltype(transList.size()) i = 0; i < transList.size(); i++)
    {
        const auto &transformer = transList.at(i);
        const auto primaryNode = transformer.getPrimaryNode() - 1;
        const auto secondNode = transformer.getSecondaryNode() - 1;
        const auto k = transformer.getRatio();

        const auto &y1 = this->SocketData1.find({transformer.getPrimaryNode(), transformer.getSecondaryNode()})->second.first;
        const auto &y2 = this->SocketData2.find({transformer.getPrimaryNode(), transformer.getSecondaryNode()})->second.first;
        const auto &y0 = this->SocketData0.find({transformer.getPrimaryNode(), transformer.getSecondaryNode()})->second.first;

        const auto delta_Y1_tt = y1 / (k * k) - y1;
        const auto delta_Y1_ft = y1 - y1 / k;

        const auto delta_Y2_tt = y2 / (k * k) - y2;
        const auto delta_Y2_ft = y2 - y2 / k;

        const auto delta_Y0_tt = y0 / (k * k) - y0;
        const auto delta_Y0_ft = y0 - y0 / k;

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

void Grid::mountIdealTransformer2_PrimarySideReactance(const std::vector<std::pair<IdealTransformer2, cf>> &idealTransList)
{
    for (decltype(idealTransList.size()) i = 0; i < idealTransList.size(); i++)
    {
        const auto &thisTrans = idealTransList.at(i).first;
        const auto y = cf{1, 0} / idealTransList.at(i).second;
        const auto pNode = thisTrans.getPrimaryNode() - 1;
        const auto qNode = thisTrans.getSecondaryNode() - 1;
        const auto k = thisTrans.getRatio();

        this->Y1(qNode, pNode) -= y / k;
        this->Y1(pNode, qNode) = this->Y1(qNode, pNode);
        this->Y1(pNode, pNode) += y;
        this->Y1(qNode, qNode) += y / (k * k);

        //TODO: 要不要+1
        this->SocketData1.insert(
            {{{pNode + 1, qNode + 1}, {y, cf(0.0f, 0.0f)}},
             {{qNode + 1, pNode + 1}, {y, cf(0.0f, 0.0f)}}});
    }
}

void Grid::mountTransformer2(const std::vector<Transformer2> &transList)
{
#pragma omp parallel for
    for (decltype(transList.size()) i = 0; i < transList.size(); i++)
    {
        const auto &thisTrans = transList.at(i);
        const auto y = cf{0, -1 / thisTrans.getXd()};
        const auto pNode = thisTrans.getPrimaryNode() - 1;
        const auto qNode = thisTrans.getSecondaryNode() - 1;

        //TODO: 这里到底是+还是-
        //this->Y1(qNode, pNode) += y;
        this->Y1(qNode, pNode) -= y;
        this->Y1(pNode, qNode) = this->Y1(qNode, pNode);
        this->Y1(pNode, pNode) += y;
        this->Y1(qNode, qNode) += y;
    }
}

void Grid::mountGenerator(const std::vector<Generator> &geneList)
{
#pragma omp parallel for
    for (decltype(geneList.size()) i = 0; i < geneList.size(); i++)
    {
        const auto &gene = geneList.at(i);
        const auto node = gene.getNode() - 1;
        const auto y = cf{1, 0} / gene.getXd();
        this->Y1(node, node) += y;
    }
}
