/*
 * @Author: your name
 * @Date: 2021-04-09 20:35:04
 * @LastEditTime: 2021-04-11 21:56:41
 * @LastEditors: your name
 * @Description: In User Settings Edit
 * @FilePath: /VSCode/Generator.h
 */
#pragma once
#include "Equipment.h"

class Generator
{
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
