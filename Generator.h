#pragma once
#include "Equipment.h"

class Generator{
private:
    NodeType Node{0};
    DeviceArgType Sn{-1.0f};
    cf xd_;
    cf Xd;

public:
    Generator(const NodeType Node_, const DeviceArgType Sn_, const cf __xd, const DeviceArgType SB_);

    NodeType getNode() const;
    DeviceArgType getSn() const;
    cf getXd() const;
};
