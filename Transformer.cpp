#include "Transformer.h"

Transformer2::Transformer2(const size_t primary_node_, const size_t secondary_node_, const float ratio_) : PrimaryNode(primary_node_), SecondaryNode(secondary_node_), Ratio(ratio_) {}

Transformer2::Transformer2(const Transformer2& trans) : PrimaryNode(trans.PrimaryNode), SecondaryNode(trans.SecondaryNode), Ratio(trans.Ratio) {}

void Transformer2::setRatio(const float newRatio)
{
    this->Ratio = newRatio;
}

float Transformer2::getRatio() const
{
    return this->Ratio;
}

size_t  Transformer2::getPrimaryNode() const
{
    return this->PrimaryNode;
}
size_t Transformer2::getSecondaryNode() const
{
    return this->SecondaryNode;
}
