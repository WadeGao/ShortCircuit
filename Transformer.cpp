#include "Transformer.h"

Transformer2::Transformer2(const NodeType primary_node_, const NodeType secondary_node_, const DeviceArgType ratio_) : PrimaryNode(primary_node_), SecondaryNode(secondary_node_), Ratio(ratio_) {}

Transformer2::Transformer2(const NodeType primary_node_, const NodeType secondary_node_, const DeviceArgType ratio_, const DeviceArgType Sn_, const DeviceArgType Vs_, const DeviceArgType SB_)
    : PrimaryNode(primary_node_), SecondaryNode(secondary_node_), Ratio(ratio_), Sn(Sn_), Vs_percent(Vs_), Xd(Vs_ * SB_ / 100 / Sn_) {}

void Transformer2::setRatio(const DeviceArgType newRatio)
{
    this->Ratio = newRatio;
}

NodeType Transformer2::getPrimaryNode() const
{
    return this->PrimaryNode;
}

NodeType Transformer2::getSecondaryNode() const
{
    return this->SecondaryNode;
}

DeviceArgType Transformer2::getRatio() const
{
    return this->Ratio;
}

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
