#pragma once
#include <cstddef>
class Transformer2
{
private:
    size_t PrimaryNode{ 0 }, SecondaryNode{ 0 };
    float Ratio{ 1 };

public:
    Transformer2(const size_t primary_node_, const size_t secondary_node_, const float ratio_);
    Transformer2(const Transformer2 &trans);
    size_t getPrimaryNode() const;
    size_t getSecondaryNode()const;
    void setRatio(const float newRatio);
    float getRatio() const;
};
