#include "Transformer.h"

Transformer2::Transformer2(const size_t primary_node_, const size_t secondary_node_, const float ratio_):PrimaryNode(primary_node_), SecondaryNode(secondary_node_), Ratio(ratio_)
{
}

void Transformer2::adjustRatio(const float newRatio)
{
    this->Ratio = newRatio;
}

float Transformer2::getRatio() const
{
    return this->Ratio;
}
