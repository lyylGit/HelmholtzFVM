#pragma once

#include <string>
#include <functional>
#include <math.h>

using namespace std;

class staticFunctions
{

private:
    /* data */
public:
    staticFunctions() {};
    ~staticFunctions() {};

    std::function<void(string)> fnShowPrompt;
    string m_strErr;

    void ShowPrompt(string s)
    {
        if (0 != fnShowPrompt)
        fnShowPrompt(s);
    }

    // 清理内存
    template <class T2>
    void _deleteArray(T2 *&a)
    {
        if (0 != a)
            delete[] a;
        a = 0;
    };

    template <class T2>
    void _deleteArray(T2 **&a, int n)
    {
        if (0 != a)
        {
            if (n > 0)
                for (int i = 0; i < n; i++)
                    delete[] (a[i]);
            delete[] a;
        }
        a = 0;
    };
    // 字符串处理
    string LTrim(const string &str)
    {
        if (str == "")
            return "";
        return str.substr(str.find_first_not_of(" \n\r\t"));
    }
    string RTrim(const string &str)
    {
        if (str == "")
            return "";
        return str.substr(0, str.find_last_not_of(" \n\r\t") + 1);
    }
    string Trim(const string &str)
    {
        if (str == "")
            return "";
        return LTrim(RTrim(str));
    }
    string ToLower(string str)
    {
        transform(str.begin(), str.end(), str.begin(), ::tolower);
        return str;
    }
    vector<string> split(const string &s, const string &seperator)
    {
        vector<string> result;
        typedef string::size_type string_size;
        string_size i = 0;

        while (i != s.size())
        {
            // 找到字符串中首个不等于分隔符的字母；
            int flag = 0;
            while (i != s.size() && flag == 0)
            {
                flag = 1;
                for (string_size x = 0; x < seperator.size(); ++x)
                    if (s[i] == seperator[x])
                    {
                        ++i;
                        flag = 0;
                        break;
                    }
            }

            // 找到又一个分隔符，将两个分隔符之间的字符串取出；
            flag = 0;
            string_size j = i;
            while (j != s.size() && flag == 0)
            {
                for (string_size x = 0; x < seperator.size(); ++x)
                    if (s[j] == seperator[x])
                    {
                        flag = 1;
                        break;
                    }
                if (flag == 0)
                    ++j;
            }
            if (i != j)
            {
                result.push_back(s.substr(i, j - i));
                i = j;
            }
        }
        return result;
    }

    // 截取浮点数
    vector<double> getMathValueDouble(string s)
    {
        vector<double> ret;
        auto sa = split(s, " ");
        try
        {
            for (unsigned int i = 0; i < sa.size(); i++)
                ret.push_back(stod(sa[i]));
        }
        catch (...)
        {
        }
        return ret;
    }
    // 截取16进制整型
    vector<long> getMathValueHex(string s)
    {
        vector<long> ret;
        auto sa = split(s, " ");
        try
        {
            for (unsigned int i = 0; i < sa.size(); i++)
                ret.push_back(stol(sa[i], 0, 16));
        }
        catch (...)
        {
        }
        return ret;
    }
    int o2i(string s)
    {
        return stoi(s);
    }

    // 判断字符串s中是否含有taget标志
    bool isTarget(string s, string taget)
    {
        if (s.length() < taget.length())
            return false;
        int n = taget.length();
        for (int i = 0; i < n; i++)
        {
            if (s[i] != taget[i])
                return false;
        }
        return true;
    }
    // 取字符串s中tag的内容
    string getContentString(string s, string tag)
    {
        s = s.substr(tag.length(), s.length() - tag.length());
        // 找到第一个)位置
        int i = s.find_first_of(')');
        if (i > 0)
            s = s.substr(0, i + 1);
        s = Trim(s);
        return s;
    }

    // 计算两点间(p,p0)的距离点p长度为h的中间点
    template <class T2>
    T2 *_calcuPoint(T2 *nodes, long p, long p0, T2 h)
    {
        T2 *ret = new T2[3];
        // p,p0的长度
        T2 len = 0;
        for (long i = 0; i < 3; i++)
        {
            T2 a = (nodes[p * 3 + i] - nodes[p0 * 3 + i]);
            len += a * a;
        }
        len = sqrt(len);

        T2 r = h / len;
        // T2 A[3];
        for (long i = 0; i < 3; i++)
        {
            ret[i] = nodes[p * 3 + i] + r * (nodes[p0 * 3 + i] - nodes[p * 3 + i]);
        }
        return ret;
    };
};
