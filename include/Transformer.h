#pragma once
#include "Equipment.h"

class IdealTransformer2
{
private:
    NodeType PrimaryNode{0}, SecondaryNode{0};
    DeviceArgType Ratio{-1.0f};

public:
    IdealTransformer2(const NodeType primary_node_, const NodeType secondary_node_, const DeviceArgType ratio_);

    void setRatio(const DeviceArgType newRatio);

    NodeType getPrimaryNode() const;
    NodeType getSecondaryNode() const;
    DeviceArgType getRatio() const;
};

class Transformer2 : public IdealTransformer2
{
private:
    DeviceArgType Sn{-1.0f};
    DeviceArgType Vs_percent{-1.0f};
    DeviceArgType Xd{-1.0f};

public:
    Transformer2(const NodeType primary_node_, const NodeType secondary_node_, const DeviceArgType ratio_, const DeviceArgType Sn_, const DeviceArgType Vs_, const DeviceArgType SB_);

    DeviceArgType getSn() const;
    DeviceArgType getXd() const;
    DeviceArgType getVs() const;
};
