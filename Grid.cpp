#include "Grid.h"
#include <iostream>

Grid::Grid(const NodeType node_, const TripVecType &lineData, const TripVecType &nodeList, const TripVecType &idealTransList, const TripVecType &transList, const TripVecType &geneList) : NodeNum(node_), myPool(maxThreadsNum)
{
    //omp_set_num_threads(std::thread::hardware_concurrency());

    this->Y2 = this->Y0 = this->Y1 = Eigen::SparseMatrix<cf>(this->NodeNum, this->NodeNum);

    Grid::wrapper(lineData, this->SocketData1);
    Grid::wrapper(nodeList, this->SocketData1);
    Grid::wrapper(idealTransList, this->SocketData1);
    Grid::wrapper(transList, this->SocketData1);
    Grid::wrapper(geneList, this->SocketData1);

    TripVecType initList;
    //initList.insert
    for(const auto &iter: this->SocketData1)
        initList.emplace_back(Eigen::Triplet<cf>(iter.first.first, iter.first.second, iter.second));

    this->Y1.setFromTriplets(initList.begin(), initList.end());

    //在这里做了对比试验，同样的数据集，采用稀疏矩阵也可以保证形成相同的Y矩阵
    //std::cout << this->Y1.selfadjointView<Eigen::Lower>() << std::endl;

    Eigen::SparseMatrix<cf> E(this->Y1.rows(), this->Y1.cols());
    E.setIdentity();
    std::cout << E << std::endl;
    Eigen::SparseLU<Eigen::SparseMatrix<cf>> solver;
    //这里有bug: 比如(8, -3)翻到上面去就变成了(8, 3)导致Z矩阵计算错误
	solver.compute(this->Y1.selfadjointView<Eigen::Lower>());
    this->Z1 = solver.solve(E);

    //std::cout << this->Z1 << std::endl;
    /*this->Z2 = this->Y2.inverse();
    this->Z0 = this->Y0.inverse();*/

}

Grid::~Grid()
{
}

void Grid::wrapper(const TripVecType &triList, std::map<socketType, cf> &sockMap)
{
//#pragma omp parallel for
    for(decltype(triList.size()) i = 0; i < triList.size(); i++)
    {
        const auto &thisTriplet = triList.at(i);
        const socketType &sock = {thisTriplet.row(), thisTriplet.col()};
        if(sockMap.find(sock) == sockMap.end())
            sockMap.insert({sock, thisTriplet.value()});
        else
            sockMap[sock] += thisTriplet.value();
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

    std::list<std::pair<socketType, cf>> shortCurrent{};

    //TODO:这里把变压器当成了k=1的变压器，或者说根本没有变压器
    for (const auto &iter : this->SocketData1)
    {
        const auto &socket = iter.first;
        const auto &y = iter.second;

        if(socket.first != socket.second)
        {
            const auto Ipq = (Ui(socket.first) - Ui(socket.second)) * y;
            shortCurrent.emplace_back(std::make_pair(socket, Ipq));
        }
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

    using SpIterType = Eigen::SparseMatrix<cf>::InnerIterator;
    for (decltype(this->Y1.outerSize()) k = 0; k < this->Y1.outerSize(); ++k)
    {
        cf Is1{0, 0}, Is2{0, 0}, Is0{0, 0};

        for (SpIterType it1(this->Y1, k), it2(this->Y2, k), it0(this->Y0, k); it1 && it2 && it0; ++it1, ++it2, ++it0)
		{
            Is1 = -(U1(it1.row()) - U1(it1.col())) * it1.value();
            Is2 = -(U2(it2.row()) - U2(it2.col())) * it2.value();
            Is0 = -(U0(it0.row()) - U0(it0.col())) * it0.value();
            const std::pair<int, int> &curSocket = std::make_pair(it0.row(), it0.col());
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

