#pragma once

#include <iostream>
#include <string>
#include <list>
#include <mysql/mysql.h>

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
    std::list<std::string> GetTables();
    //返回数据库中表的集合到屏幕
    void GetTables2Screen();
    //读取数据库内容
    std::list<std::list<std::string>> ReadMySQL(const std::string &query);
    //读取数据库内容到屏幕
    void ReadMySQL2Screen(const std::string &query);
};
