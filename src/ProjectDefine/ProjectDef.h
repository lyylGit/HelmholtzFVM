#pragma once
#include <iostream>

#include <cstring>
#include <functional>

#include "ProjectDefType.h"
#include "staticFunctions.h"

#include "json.hpp"

#include "mpParser.h"
#include "mpDefines.h"

using namespace std;

using json = nlohmann::json;
using namespace mup;

// 使用json定义工程
class ProjectConfig
{

    map<string, defCaseType> DefCaseType = {
        {"HelmFreq", defCaseType::HelmFreq},
        {"HelmTime", defCaseType::HelmTime},
        {"HelmModal", defCaseType::HelmModal},
        {"LEE", defCaseType::LEE},
    };
    map<string, defMeshValueType> DefMeshValueType = {
        {"float", defMeshValueType::Float},
        {"double", defMeshValueType::Double},
    };
    map<string, defMeshFileType> DefMeshFileType = {
        {"fluent_txt", defMeshFileType::FluentTxt},
        {"fluent_bin", defMeshFileType::FluentBin},
    };
    map<string, defFieldType> DefFieldType = {
        {"uniform", defFieldType::Uniform},
        {"field", defFieldType::Field},
    };

    // 运算符分割
    vector<char> split = {'+', '-', '/', '*', '^', '(', ')', ' '};
    // 保留变量
    vector<string> keyParam = {"Hz"};

    string confFileName;

public:
    // 结构
    string caseName, meshFile, outPath;
    // bool outpathclean=false;
    double mm = 1.0;
    defCaseType caseType;
    defMeshValueType tCase;
    defMeshFileType meshFileType;

    double freqFrom, freqTo, freqDelta;

    defFieldValueDefine rho, p0, gamma;
    defHelmModal helmModal;
    defHelmFreq helmFreq;

    vector<shared_ptr<defFTF>> ftfs;
    vector<shared_ptr<defBC>> bcs;

    vector<shared_ptr<EqDef>> eqDefs;

    //
    staticFunctions SF;
    map<string, double> param;

public:
    ProjectConfig(string fn) : confFileName(fn)
    {
    }

    bool parseProj()
    {
        ifstream fi(confFileName);
        if (!fi.good())
            return false;
        json data = json::parse(fi, nullptr, true, true);
        // 项目名称
        caseName = data["casename"].get<string>();
        // 项目类型
        string s = data["casetype"].get<string>();
        if (DefCaseType.count(s) <= 0)
        {
            cout << "casetype: " << s << " 无定义\n";
            return false;
        }
        caseType = DefCaseType[s];
        // 网格数据类型
        s = data["tcase"].get<string>();
        if (DefMeshValueType.count(s) <= 0)
        {
            cout << "casetype: " << s << " 无定义\n";
            return false;
        }
        tCase = DefMeshValueType[s];
        // 网格文件类型
        s = data["meshfiletype"].get<string>();
        if (DefMeshFileType.count(s) <= 0)
        {
            cout << "casetype: " << s << " 无定义\n";
            return false;
        }
        meshFileType = DefMeshFileType[s];

        if (data.contains("mm"))
            mm = data["mm"].get<double>();

        meshFile = data["meshfile"].get<string>();
        outPath = data["outpath"].get<string>();

        if(!data["outpathclean"].empty())
        {
            if(data["outpathclean"].get<bool>())
            {
                /* 该代码不安全，取消使用
                std::string command = "rm -rf " + std::string(outPath) + "/*";  // 构造删除命令
                int result = std::system(command.c_str());  // 调用系统命令
                if (result == 0) 
                {
                    std::cout << "成功删除文件夹：" << outPath << std::endl;
                } 
                else 
                {
                    std::cout << "删除文件夹失败：" << outPath << std::endl;
                }
                */
                std::cout << "删除文件夹失败：" << outPath << std::endl;
            }
        }

        // 频率；
        auto freq = data["freqrange"];
        if (!freq.empty() && freq.size() >= 3)
        {
            freqFrom = freq.at(0).get<double>();
            freqTo = freq.at(1).get<double>();
            freqDelta = freq.at(2).get<double>();
        }
        else
            cout << "freqrange: " << s << " 无定义\n";

        // feast;
        auto fst = data["HelmModal"];
        if (!fst.empty())
        {
            helmModal.output = fst["output"].get<int>();
            helmModal.err = fst["err"].get<int>();
            helmModal.criteria = fst["criteria"].get<int>();
            helmModal.center =(freqFrom+freqTo)/2.0;
            helmModal.r = (freqTo-freqFrom)/2.0;
            helmModal.M0 = fst["M0"].get<int>();
            helmModal.it = fst["it"].get<int>();
        }
        // HelmFreq
        auto hfreq = data["HelmFreq"];
        if (!hfreq.empty())
        {
            auto J = hfreq["source"];
            if (!J.empty() && !parseSource(helmFreq.sources, J))
            {
                cout << "解析源项出错\n";
                return false;
            }
            J = hfreq["output"];
            if (!J.empty() && !parseOutput(helmFreq.outputs, J))
            {
                cout << "解析源项出错\n";
                return false;
            }
        }

        // 变量
        param.clear();
        auto _para = data["param"];
        for (int i = 0; i < _para.size(); i++)
        {
            auto item = _para.at(i);
            if (item.contains("name") && item.contains("value"))
                param["$" + item["name"].get<string>()] = item["value"].get<double>();
        }
        // 场
        auto fields = data["fields"];
        for (int i = 0; i < fields.size(); i++)
        {
            auto item = fields.at(i);
            string name = item["name"].get<string>();
            if (name == "rho" && !parseField(rho, item))
            {
                cout << "解析 rho 出错\n";
                return false;
            }
            if (name == "p0" && !parseField(p0, item))
            {
                cout << "解析 p0 出错\n";
                return false;
            }
            if (name == "gamma" && !parseField(gamma, item))
            {
                cout << "解析 gamma 出错\n";
                return false;
            }
        }
        // ftfs
        auto f = data["ftf"];
        if (!f.empty() && !parseFTFs(ftfs, f))
        {
            cout << "解析 ftf 出错\n";
            return false;
        }
        //bcs
        f = data["bcs"];
        if (!f.empty() && !parseBCs(bcs, f))
        {
            cout << "解析 bcs 出错\n";
            return false;
        }
        // 方程
        auto eq = data["equations"];
        for (int i = 0; i < eq.size(); i++)
        {
            auto item = eq.at(i);

            auto eqdef = make_shared<EqDef>();
            eqDefs.push_back(eqdef);

            eqdef->name = item["name"].get<string>();
        }

        return true;
    }
    
    bool parseSource(vector<shared_ptr<defSource>> &F, json &J)
    {
        for (int i = 0; i < J.size(); i++)
        {
            auto K = J.at(i);
            auto source = make_shared<defSource>();
            F.push_back(source);
            source->type = K["type"].get<string>();
            auto pos = K["pos"];
            if (pos.is_array() && pos.size() >= 3)
            {
                source->pos.x = pos.at(0).get<double>();
                source->pos.y = pos.at(1).get<double>();
                source->pos.z = pos.at(2).get<double>();
            }
            if (!K["q"].empty())
                source->q = K["q"].get<double>();
            else
            {
                cout << "source 没有定义 q\n";
                return false;
            }
        }

        return true;
    }
    bool parseOutput(vector<shared_ptr<defOutput>> &F, json &J)
    {
        for (int i = 0; i < J.size(); i++)
        {
            auto K = J.at(i);
            auto source = make_shared<defOutput>();
            F.push_back(source);
            source->filename = K["filename"].get<string>();
            auto pos = K["mic"];
            if (pos.is_array() && pos.size() >= 3)
            {
                source->mic.x = pos.at(0).get<double>();
                source->mic.y = pos.at(1).get<double>();
                source->mic.z = pos.at(2).get<double>();
            }
            if (!K["value"].empty())
                source->value = K["value"].get<string>();
            else
            {
                cout << "output 没有定义 value\n";
                return false;
            }
        }
        return true;
    }
    bool parseField(defFieldValueDefine &F, json &J)
    {
        string s = J["type"];
        if (DefFieldType.count(s) <= 0)
        {
            cout << "解析 field type: " << s << " 出错\n";
            return false;
        }
        F.type = DefFieldType[s];
        switch (F.type)
        {
        case defFieldType::Uniform:
        {
            if (J["value"].is_number())
                F.value = J["value"].get<double>();
            else
                F.value = calcu(J["value"].get<string>());
            break;
        }
        case defFieldType::Field:
        {
            auto vv = J["value"];
            F.filename = vv["filename"].get<string>();
            auto secs = vv["sections"];
            if (!secs.empty())
            {
                for (int i = 0; i < secs.size(); i++)
                {
                    auto sec = secs.at(i);
                    auto def = make_shared<defSection>();
                    def->dimension = sec["dimension"].get<string>();
                    F.sections.push_back(def);
                    auto ss = sec["segments"];
                    if (ss.empty())
                    {
                        cout << "定义segments 出错\n";
                        return false;
                    }
                    for (int j = 0; j < ss.size(); j++)
                    {
                        s = ss.at(j).get<string>();
                        def->segments.push_back(getSegment(s));
                    }
                }
            }
            break;
        }
        default:
        {
            cout << "定义segments 出错\n";
            return false;
        }
        }
        return true;
    }
    tuple<double, double, int, int, double> getSegment(string s)
    {
        auto sa = SF.split(s, ":");
        if (sa.size() != 2)
        {
            cout << "定义segments 出错\n";
            return make_tuple(-1, -1, -1, -1, -1);
        }
        string s0 = SF.Trim(sa[0]), s1 = SF.Trim(sa[1]);
        double a, b, c;
        int d = -1, e = -1;

        if (!getSegmentRange(s0, a, b, d, e))
            return make_tuple(-1, -1, -1, -1, -1);
        c = calcu(s1);

        return make_tuple(a, b, d, e, c);
    }
    bool getSegmentRange(string s0, double &from, double &to, int &d, int &e)
    {
        char cc = s0[0];
        if (cc == '(')
            d = 0;
        else if (cc == '[')
            d = 1;
        s0 = s0.substr(1);
        cc = s0[s0.length() - 1];
        if (cc == ')')
            e = 0;
        else if (cc == ']')
            e = 1;
        if (d == -1 || e == -1)
        {
            cout << "定义segments 范围 " << s0 << " 出错\n";
            return false;
        }
        s0 = s0.substr(0, s0.length() - 1);
        auto sa = SF.split(s0, "|");
        if (sa.size() != 2)
        {
            cout << "定义segments 范围 " << s0 << " 出错\n";
            return false;
        }
        from = calcu(sa[0]);
        to = calcu(sa[1]);
        return true;
    }
    double calcu(string s)
    {
        if (s.find_first_of("$") >= 0)
        {
            vector<string> left;
            s = replacePara(s, param, left);
            if (left.size() > 0)
            {
                cout << "定义segments:" + s + " 出错\n";
                return 0;
            }
        }

        ParserX parser;
        parser.SetExpr(s.c_str());
        return parser.Eval().GetFloat();
    }
    complex<double> calcuZ(string s)
    {
        if (s.find_first_of("$") >= 0)
        {
            vector<string> left;
            s = replacePara(s, param, left);
            if (left.size() > 0)
            {
                cout << "定义segments:" + s + " 出错\n";
                return 0;
            }
        }

        ParserX parser;
        parser.SetExpr(s.c_str());
        return parser.Eval().GetComplex();
    }
    //  complex<double> calcuZ(string s,vector<tuple<string,double>>paraList)
    // {

    //     ParserX parser;
    //     parser.SetExpr(s.c_str());
    //     parser.DefineVar(
    //     return parser.Eval().GetComplex();
    // }
    bool parseFTFs(vector<shared_ptr<defFTF>> &F, json &J)
    {
        for (int i = 0; i < J.size(); i++)
        {
            auto ftf = make_shared<defFTF>();
            F.push_back(ftf);
            auto f = J.at(i);
            ftf->order = f["order"].get<int>();
            auto g = f["posref"];
            if (g.empty())
            {
                cout << "解析 FTF 出错,没有" << "posref\n";
                return false;
            }
            ftf->posref.x = g.at(0).get<double>();
            ftf->posref.y = g.at(1).get<double>();
            ftf->posref.z = g.at(2).get<double>();

            g = f["nref"];
            if (g.empty())
            {
                cout << "解析 FTF 出错,没有" << "nref\n";
                return false;
            }
            ftf->nref.x = g.at(0).get<double>();
            ftf->nref.y = g.at(1).get<double>();
            ftf->nref.z = g.at(2).get<double>();

            g = f["heatsection"];
            if (g.empty())
            {
                cout << "解析 FTF 出错,没有" << "heatsection \n";
                return false;
            }
            ftf->heatsection.dimension = g["dimension"].get<string>();
            auto section = g["segments"];
            if (section.empty())
            {
                cout << "解析 FTF 出错,没有" << "segments \n";
                return false;
            }
            for (int k = 0; k < section.size(); k++)
            {
                string s0 = section.at(i).get<string>();
                ftf->heatsection.segments.push_back(getSegment(s0));
            }
            auto model = f["model"];
            if (model.empty())
            {
                cout << "解析 FTF 出错,没有" << "model \n";
                return false;
            }
            ftf->ftfmodel.type = model["type"].get<string>();
            if (ftf->ftfmodel.type == "nt")
            {
                if (model["tao"].is_number())
                    ftf->ftfmodel.tao = model["tao"].get<double>();
                else
                    ftf->ftfmodel.tao = calcu(model["tao"].get<string>());

                if (model["n"].is_number())
                    ftf->ftfmodel.n = model["n"].get<double>();
                else
                    ftf->ftfmodel.n = calcu(model["n"].get<string>());
            }
            else
            {
                cout << "未支持的ftf model\n";
                return false;
            }
        }
        return true;
    }
    bool parseBCs(vector<shared_ptr<defBC>> &F, json &J)
    {
        if (!J.empty())
        {
            for (int j = 0; j < J.size(); j++)
            {
                auto bc = J.at(j);
                auto defbc = make_shared<defBC>();
                F.push_back(defbc);
                if (bc.contains("id"))
                {
                    auto ids=bc["id"];
                    if (ids.is_array())
                    {
                        for (int k = 0; k < ids.size(); k++)
                            defbc->ids.push_back(ids.at(k).get<int>());
                    }
                }
                if (!bc.contains("type"))
                {
                    cout << "bcs 定义缺失 type\n";
                    return false;
                }
                defbc->type = bc["type"].get<string>();

                if (bc.contains("order"))
                    defbc->order = bc["order"].get<int>();

                if (bc.contains("grad"))
                {
                    auto g = bc["grad"];
                    if (g.size() >= 2)
                        defbc->v = g.at(0).get<double>() + 1i * g.at(1).get<double>();
                    else if (g.size() >= 1)
                        defbc->v = g.at(0).get<double>();
                }

                if (bc.contains("Z"))
                {
                    vector<string> left;
                    string formula = bc["Z"].get<string>();
                    formula = replacePara(formula, param, defbc->param);
                    if (defbc->param.size() <= 0)
                        defbc->v = calcuZ(formula);
                    else
                        defbc->formula = formula;
                }
            }
        }
        return true;
    }
    string replacePara(string formula, map<string, double> &param, vector<string> &left)
    {
        string ret = "";
        bool b = false;
        string p = "";
        for (int i = 0; i < formula.length(); i++)
        {
            auto c = formula[i];
            if (b)
            {
                if (find(split.begin(), split.end(), c) != split.end())
                {
                    if (param.count(p) > 0)
                    {
                        ret += to_string(param[p]) + c;
                    }
                    else
                    {
                        // 去掉了$符号
                        ret += p.substr(1) + c;
                        left.push_back(p.substr(1));
                    }

                    p = "";
                    b = false;
                }
                else
                    p += c;
            }
            else if (c == '$')
            {
                p += c;
                b = true;
            }
            else
                ret += c;
        }
        if (b)
        {
            if (param.count(p) > 0)
            {
                ret += to_string(param[p]);
            }
            else
            {
                // 去掉了$符号
                ret += p.substr(1);
                left.push_back(p.substr(1));
            }
        }
        return ret;
    }
};