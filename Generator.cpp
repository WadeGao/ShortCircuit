#include "Generator.h"

Generator::Generator(const NodeType Node_, const DeviceArgType Sn_, const DeviceArgType __xd, const DeviceArgType SB_) : Node(Node_), Sn(Sn_), xd_(__xd), Xd(__xd * SB_ / Sn_) {}

NodeType Generator::getNode() const
{
    return this->Node;
}

DeviceArgType Generator::getSn() const
{
    return this->Sn;
}

DeviceArgType Generator::getXd() const
{
    return this->Xd;
}
