#pragma once

#include "Common.h"
#include <iostream>
#include <mysql/mysql.h>
#include <string>

using DatabaseTableType = std::vector<std::vector<std::string>>;

class Database
{
private:
    MYSQL mysql{};

public:
    Database();
    //断开数据库连接
    ~Database();
    //连接到数据库
    bool ConnectMySQL(const char *host, const char *user, const char *db, unsigned int port = 3306);
    //从数据库中创建表、更新表、 删除表
    bool CUD_MySQL(const std::string &query);
    //返回数据库中表的集合
    std::vector<std::string> GetTables();
    //读取数据库内容
    DatabaseTableType ReadMySQL(const std::string &query);

    std::tuple<size_t, size_t> getTableSize(const std::string &tableName);
};

