#pragma once
#include "Equipment.h"

class Generator{
private:
    NodeType Node{0};
    DeviceArgType Sn{-1.0f};
    DeviceArgType xd_{-1.0f};
    DeviceArgType Xd{-1.0f};

public:
    Generator(const NodeType Node_, const DeviceArgType Sn_, const DeviceArgType __xd, const DeviceArgType SB_);

    NodeType getNode() const;
    DeviceArgType getSn() const;
    DeviceArgType getXd() const;
};
