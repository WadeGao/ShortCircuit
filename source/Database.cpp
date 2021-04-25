/*
 * @Author: your name
 * @Date: 2021-04-25 21:59:06
 * @LastEditTime: 2021-04-25 22:16:27
 * @LastEditors: your name
 * @Description: In User Settings Edit
 * @FilePath: /ShortCircuit/source/Database.cpp
 */
#include "Database.h"

Database::Database()
{
    mysql_init(&mysql);
}

Database::~Database()
{
    mysql_close(&mysql);
    //fprintf(stdout, "Database Closed\n");
}

bool Database::ConnectMySQL(const char *host, const char *user, const char *db, unsigned int port)
{
    char pwd[32]{0};
    uint8_t failedTimes = 0;
    while (true)
    {
        fprintf(stdout, "Please input Database's password: ");

        std::cin.getline(pwd, sizeof(pwd));

        if (!mysql_real_connect(&this->mysql, host, user, pwd, db, port, nullptr, CLIENT_FOUND_ROWS))
        {
            fprintf(stderr, "Connection Failed: %s\n\n", mysql_error(&this->mysql));
            if (++failedTimes >= 3)
                return false;
        }
        else
            break;
    }
    return true;
}

bool Database::CUD_MySQL(const std::string &query)
{
    const auto res = mysql_query(&this->mysql, query.c_str());
    const bool isSuccessful = (res == 0);
    isSuccessful ? fprintf(stdout, "%s\tQuery Successfully!\n", query.c_str()) : fprintf(stderr, "%s\tSyntax Error!\n", query.c_str());
    return isSuccessful;
}

std::vector<std::string> Database::GetTables()
{
    std::vector<std::string> allTables;
    CUD_MySQL("SHOW TABLES;");

    MYSQL_RES *result = mysql_store_result(&this->mysql);
    MYSQL_ROW row;
    while ((row = mysql_fetch_row(result)))
        allTables.emplace_back(row[0]);

    return allTables;
}

DatabaseTableType Database::ReadMySQL(const std::string &query)
{
    if (mysql_query(&this->mysql, query.c_str()))
    {
        std::vector<std::vector<std::string>> ret;
        return ret;
    }
    MYSQL_RES *result = mysql_store_result(&this->mysql);

    auto row_count = mysql_num_rows(result);
    auto field_count = mysql_num_fields(result);

    MYSQL_ROW row;
    decltype(row_count) cur_line = 0;

    //这里，初始化指定大小，避免push_back导致内存重新分配引起的迭代器失效，进而避免多线程读写的问题
    //下面代码使用下标访问，赋值操作不是插入新元素，而是修改旧元素
    std::vector<std::vector<std::string>> res(row_count, std::vector<std::string>(field_count, ""));
    while ((row = mysql_fetch_row(result)))
    {
#pragma omp parallel for
        for (decltype(field_count) i = 0; i < field_count; i++)
            res[cur_line][i] = row[i];
        cur_line++;
    }
    mysql_free_result(result);
    return res;
}

std::tuple<size_t, size_t> Database::getQuerySize(const std::string &query)
{
    //const std::string query = "SELECT * FROM " + tableName;
    if (mysql_query(&this->mysql, query.c_str()))
        return {0, 0};

    MYSQL_RES *result = mysql_store_result(&this->mysql);
    return {mysql_num_rows(result), mysql_num_fields(result)};
}
