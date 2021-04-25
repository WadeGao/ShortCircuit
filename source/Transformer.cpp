/*
 * @Author: your name
 * @Date: 2021-04-08 19:52:30
 * @LastEditTime: 2021-04-11 21:56:21
 * @LastEditors: your name
 * @Description: In User Settings Edit
 * @FilePath: /VSCode/Transformer.cpp
 */
#include "Transformer.h"

IdealTransformer2::IdealTransformer2(const NodeType primary_node_, const NodeType secondary_node_, const DeviceArgType ratio_) : PrimaryNode(primary_node_), SecondaryNode(secondary_node_), Ratio(ratio_) {}

void IdealTransformer2::setRatio(const DeviceArgType newRatio)
{
    this->Ratio = newRatio;
}

NodeType IdealTransformer2::getPrimaryNode() const
{
    return this->PrimaryNode;
}

NodeType IdealTransformer2::getSecondaryNode() const
{
    return this->SecondaryNode;
}

DeviceArgType IdealTransformer2::getRatio() const
{
    return this->Ratio;
}

Transformer2::Transformer2(const NodeType primary_node_, const NodeType secondary_node_, const DeviceArgType ratio_, const DeviceArgType Sn_, const DeviceArgType Vs_, const DeviceArgType SB_)
    : IdealTransformer2(primary_node_, secondary_node_, ratio_), Sn(Sn_), Vs_percent(Vs_), Xd(Vs_ * SB_ / 100 / Sn_) {}

DeviceArgType Transformer2::getXd() const
{
    return this->Xd;
}

DeviceArgType Transformer2::getSn() const
{
    return this->Sn;
}

DeviceArgType Transformer2::getVs() const
{
    return this->Vs_percent;
}
