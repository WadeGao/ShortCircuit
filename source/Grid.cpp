#define EIGEN_USE_MKL_ALL
#include "Grid.h"
#include <iostream>

Grid::Grid(const NodeType node_, const TripVecType &lineData, const TripVecType &nodeList, const TripVecType &idealTransList, const TripVecType &transList, const TripVecType &geneList) : NodeNum(node_), myPool(maxThreadsNum)
{
    //omp_set_num_threads(std::thread::hardware_concurrency());
    this->Y2 = this->Y0 = this->Y1 = Eigen::SparseMatrix<cf>(node_, node_);

    Grid::wrapper(lineData, this->SocketData1);
    Grid::wrapper(nodeList, this->SocketData1);
    Grid::wrapper(idealTransList, this->SocketData1);
    Grid::wrapper(transList, this->SocketData1);
    Grid::wrapper(geneList, this->SocketData1);

    {
        std::list<Eigen::Triplet<cf>> initList{};
        for (const auto &iter : this->SocketData1)
        {
            const auto &sock = iter.first;
            initList.emplace_back(Eigen::Triplet<cf>(sock.first, sock.second, iter.second));
            //不这么搞，就(8, 3) -> (8, -3)
            if (sock.first != sock.second)
                initList.emplace_back(Eigen::Triplet<cf>(sock.second, sock.first, iter.second));
        }
        //在这里做了对比试验，同样的数据集，采用稀疏矩阵也可以保证形成与密集情况下相同的Y矩阵
        this->Y1.setFromTriplets(initList.begin(), initList.end());
        this->Y2 = this->Y1;
    }

    Eigen::SparseMatrix<cf> E(node_, node_);
    E.setIdentity();
    Eigen::SparseLU<Eigen::SparseMatrix<cf>> solver;
    solver.compute(this->Y1);
    //这里，通过大量的对比试验。证明密集矩阵的inverse()方法不精确
    //下面的语句是精确的求逆方法。
    this->Z0 = this->Z2 = this->Z1 = solver.solve(E);
}

Grid::~Grid()
{
}

void Grid::wrapper(const TripVecType &triList, std::map<socketType, cf> &sockMap)
{
    for (decltype(triList.size()) i = 0; i < triList.size(); i++)
    {
        const auto &thisTriplet = triList.at(i);
        const socketType &sock = {thisTriplet.row(), thisTriplet.col()};
        if (sockMap.find(sock) == sockMap.end())
            sockMap.insert({sock, thisTriplet.value()});
        else
            sockMap[sock] += thisTriplet.value();
    }
}

//TODO:对称短路计算
lllReturnType Grid::SymmetricShortCircuit(const NodeType shortPoint, const cf &Zf, const DeviceArgType Uav)
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

        if (socket.first != socket.second)
        {
            const auto Ipq = (Ui(socket.first) - Ui(socket.second)) * y;
            shortCurrent.emplace_back(std::make_pair(socket, Ipq));
        }
    }
    return lllReturnType(Ui, shortCurrent); //std::make_tuple(Ui, shortCurrent);
}

//单相短路计算
cf Grid::lgShortCircuit(const NodeType faultNode, const cf &Zf)
{
    const auto faultNode2 = faultNode, faultNode0 = faultNode;

    const cf Ifa1 = cf{1, 0} / (Z1(faultNode - 1, faultNode - 1) + Z2(faultNode2 - 1, faultNode2 - 1) + Z0(faultNode0 - 1, faultNode0 - 1) + Zf * cf{3, 0});
    Eigen::VectorXcf If(3, 1);
    //auto Ifa2{ Ifa1 }, Ifa0{ Ifa1 };
    //If << Ifa1, Ifa2, Ifa0;
    If << Ifa1, Ifa1, Ifa1;
    const Eigen::VectorXcf If_abc = St * If;

    /*fprintf(stdout, "故障点abc相短路电流: \n");
    for (decltype(If_abc.size()) i = 0; i < If_abc.size(); i++)
        fprintf(stdout, "%f\n", abs(If_abc(i)));*/

    const auto Is = If_abc(0);
    //std::cout << "故障电流: " << Is << std::endl;

    return Is;
    //this->getBusVoltageAndCurrent(If, faultNode);
}

//两相短路计算
cf Grid::llShortCircuit(const NodeType faultNode, const cf &Zf)
{
    const auto faultNode2 = faultNode;

    Eigen::VectorXcf If(3, 1);
    const cf Ifa1 = cf{1, 0} / (Z1(faultNode - 1, faultNode - 1) + Z2(faultNode2 - 1, faultNode2 - 1) + Zf);
    const cf Ifa2{-Ifa1.real(), -Ifa1.imag()}, Ifa0{0, 0};
    If << Ifa1, Ifa2, Ifa0;
    const Eigen::VectorXcf If_abc = St * If;

    /*fprintf(stdout, "故障点abc相短路电流: \n");
    for (decltype(If_abc.size()) i = 0; i < If_abc.size(); i++)
        fprintf(stdout, "%f\n", abs(If_abc(i)));*/

    const auto Is = If_abc(2);
    //std::cout << "故障电流: " << Is << std::endl;
    return Is;
    //this->getBusVoltageAndCurrent(If, faultNode);
}

//两相对地短路计算
cf Grid::llgShortCircuit(const NodeType faultNode, const cf &Zf)
{
    const auto faultNode2 = faultNode, faultNode0 = faultNode;

    Eigen::VectorXcf If(3, 1);
    const cf Ifa1 = cf{1, 0} / (Z1(faultNode - 1, faultNode - 1) + Z2(faultNode2 - 1, faultNode2 - 1) * (Z0(faultNode0 - 1, faultNode0 - 1) + Zf * cf{3, 0}) / (Z2(faultNode2 - 1, faultNode2 - 1) + Z0(faultNode0 - 1, faultNode0 - 1) + Zf * cf{3, 0}));
    const cf Ifa2 = -(cf{1, 0} - Z1(faultNode - 1, faultNode - 1) * Ifa1) / Z2(faultNode2 - 1, faultNode2 - 1);
    const cf Ifa0 = -(cf{1, 0} - Z1(faultNode - 1, faultNode - 1) * Ifa1) / (Z0(faultNode2 - 1, faultNode2 - 1) + Zf * cf{0, 3});

    If << Ifa1, Ifa2, Ifa0;
    const Eigen::VectorXcf If_abc = St * If;

    /*fprintf(stdout, "故障点abc相短路电流: \n");
    for (int i = 0; i < If_abc.size(); i++)
        fprintf(stdout, "%f\n", abs(If_abc(i)));*/

    const auto Is = If_abc(1) + If_abc(2);
    //std::cout << "故障电流: " << Is << std::endl;
    return Is;
    //this->getBusVoltageAndCurrent(If, faultNode);
}

SequenceVoltType Grid::getSequenceVolt(const Eigen::VectorXcf &If, const NodeType faultNode)
{
    auto branchNum = this->Y1.rows();
    Eigen::MatrixXcf U1 = Eigen::MatrixXcf(branchNum, 1);
    auto U2{U1}, U0{U1};
#pragma omp parallel for
    for (decltype(branchNum) i = 0; i < branchNum; i++)
    {
        U1(i) = cf{1, 0} - this->Z1(faultNode - 1, i) * If(0);
        U2(i) = -(this->Z2(faultNode - 1, i) * If(1));
        U0(i) = -(this->Z0(faultNode - 1, i) * If(2));
    }

    return {U1, U2, U0};
}

Eigen::MatrixXf Grid::getBusVoltage(const Eigen::VectorXcf &If, const NodeType faultNode)
{
    const auto &seqVolt = this->getSequenceVolt(If, faultNode);
    auto branchNum = this->Y1.rows();
    Eigen::MatrixXf Uabc = Eigen::MatrixXf(branchNum, 3);
    auto &U1 = std::get<0>(seqVolt), &U2 = std::get<1>(seqVolt), &U0 = std::get<2>(seqVolt);

#pragma omp parallel for
    for (decltype(branchNum) i = 0; i < branchNum; i++)
        Uabc.row(i) = (St * ((Eigen::MatrixXcf(3, 1) << U1(i), U2(i), U0(i)).finished())).cwiseAbs();

    /*for (decltype(branchNum) i = 0; i < branchNum; i++)
        std::cout << "节点" << i + 1 << "的abc相短路电压: " << Uabc.row(i) << std::endl;*/
    return Uabc;
}

SequenceCurrentType Grid::getBusCurrent(const Eigen::VectorXcf &If, const NodeType faultNode)
{
    auto branchNum = this->Y1.rows();
    const auto &seqVolt = this->getSequenceVolt(If, faultNode);
    auto &U1 = std::get<0>(seqVolt), &U2 = std::get<1>(seqVolt), &U0 = std::get<2>(seqVolt);

    std::map<std::pair<int, int>, cf> Ib1{}, Ib2{}, Ib0{};

    using SpIterType = Eigen::SparseMatrix<cf>::InnerIterator;

    auto task1 = [&Ib1, &U1, this]() -> void {
        for (decltype(this->Y1.outerSize()) k = 0; k < this->Y1.outerSize(); ++k)
            for (SpIterType it1(this->Y1, k); it1; ++it1)
                Ib1.insert({{it1.row(), it1.col()}, -(U1(it1.row()) - U1(it1.col())) * it1.value()});
    };
    auto task2 = [&Ib2, &U2, this]() -> void {
        for (decltype(this->Y2.outerSize()) k = 0; k < this->Y2.outerSize(); ++k)
            for (SpIterType it2(this->Y1, k); it2; ++it2)
                Ib2.insert({{it2.row(), it2.col()}, -(U2(it2.row()) - U2(it2.col())) * it2.value()});
    };
    auto task0 = [&Ib0, &U0, this]() -> void {
        for (decltype(this->Y0.outerSize()) k = 0; k < this->Y0.outerSize(); ++k)
            for (SpIterType it0(this->Y0, k); it0; ++it0)
                Ib0.insert({{it0.row(), it0.col()}, -(U0(it0.row()) - U0(it0.col())) * it0.value()});
    };

    std::thread calcIbx[3]{std::thread(task0), std::thread(task1), std::thread(task2)};
    for (auto &th : calcIbx)
        th.join();

    SequenceCurrentType Isx{};

    for (const auto &iter1 : Ib1)
    {
        std::vector<cf> ins(3, cf{0, 0});
        ins.at(0) = iter1.second;
        Isx.insert({iter1.first, ins});
    }

    for (const auto &iter2 : Ib2)
    {
        if (Isx.find(iter2.first) == Isx.end())
        {
            std::vector<cf> ins(3, cf{0, 0});
            ins.at(1) = iter2.second;
            Isx.insert({iter2.first, ins});
        }
        else
            Isx[iter2.first].at(1) = iter2.second;
    }

    for (const auto &iter0 : Ib0)
    {
        if (Isx.find(iter0.first) == Isx.end())
        {
            std::vector<cf> ins(3, cf{0, 0});
            ins.at(2) = iter0.second;
            Isx.insert({iter0.first, ins});
        }
        else
            Isx[iter0.first].at(2) = iter0.second;
    }

    /*
    for (const auto &iter : Isx)
    {
        const auto &thisTuple = iter.second;
        Eigen::MatrixXcf I = St * ((Eigen::MatrixXcf(3, 1) << thisTuple.at(0), thisTuple.at(1), thisTuple.at(2)).finished());
        std::cout << "节点" << iter.first.first << "与节点" << iter.first.second << "间的abc相短路电流: " << std::endl;
        std::cout << I << std::endl;
    }
    */

    return Isx;
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

std::list<std::future<lllReturnType>> Grid::lllWholeGridScan()
{
    std::list<std::future<lllReturnType>> results;

    auto task = [this](const NodeType node, const cf &z, const DeviceArgType u) -> lllReturnType { return this->SymmetricShortCircuit(node, z, u); };
    for (int a = 1; a <= this->getNodeNum(); a++)
        for (decltype(this->NodeNum) i = 1; i <= this->NodeNum; i++)
            results.emplace_back(this->myPool.enqueue(task, i, cf(0.0f, 0.0f), 1));

    return results;
}

std::list<std::future<cf>> Grid::llgWholeGridScan()
{
    std::list<std::future<cf>> results;

    auto task = [this](const NodeType node, const cf &z) -> cf { return this->llgShortCircuit(node, z); };
    for (int a = 1; a <= this->getNodeNum(); a++)
        for (decltype(this->NodeNum) i = 1; i <= this->NodeNum; i++)
            results.emplace_back(this->myPool.enqueue(task, i, cf(0.0f, 0.0f)));

    return results;
}

std::list<std::future<cf>> Grid::llWholeGridScan()
{
    std::list<std::future<cf>> results;

    auto task = [this](const NodeType node, const cf &z) -> cf { return this->llShortCircuit(node, z); };
    for (int a = 1; a <= this->getNodeNum(); a++)
        for (decltype(this->NodeNum) i = 1; i <= this->NodeNum; i++)
            results.emplace_back(this->myPool.enqueue(task, i, cf(0.0f, 0.0f)));

    return results;
}

std::list<std::future<cf>> Grid::lgWholeGridScan()
{
    std::list<std::future<cf>> results;

    auto task = [this](const NodeType node, const cf &z) -> cf { return this->lgShortCircuit(node, z); };
    for (int a = 1; a <= this->getNodeNum(); a++)
        for (decltype(this->NodeNum) i = 1; i <= this->NodeNum; i++)
            results.emplace_back(this->myPool.enqueue(task, i, cf(0.0f, 0.0f)));

    return results;
}
NodeType Grid::getNodeNum() const { return this->NodeNum; }