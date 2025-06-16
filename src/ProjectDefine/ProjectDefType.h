#pragma once

#include <string>
#include <vector>
#include <tuple>
#include "vector3d.h"
#include <complex>
#include <memory>

using namespace std;
// 定义计算问题的一些常数

enum class defCaseType
{
    HelmFreq,
    HelmTime,
    HelmModal,
    LEE
};

enum class defMeshValueType
{
    Float,
    Double
};
enum class defMeshFileType
{
    FluentTxt,
    FluentBin,
};

enum class defFieldType
{
    Uniform,
    Field
};
class defHelmModal
{
public:
    int output = 1, err = 4, criteria = 0, M0 = 180, it = 0;
    double center, r;
};
class defSource
{
public:
    string type;
    vector3D<double> pos;
    double q;
};
class defOutput
{
public:
    string filename;
    vector3D<double> mic;
    string value;
};
class defHelmFreq
{
public:
    vector<shared_ptr<defSource>> sources;
    vector<shared_ptr<defOutput>> outputs;
};

class defSection
{
public:
    string dimension;
    // from-to-tag-tag-val
    // tag=0,表示开，tag=1 表示闭
    vector<tuple<double, double, int, int, double>> segments;
};
// 变量定义
class defFieldValueDefine
{
public:
    defFieldType type;
    double value;
    vector<shared_ptr<defSection>> sections;
    string filename;
};
class defFTFModel
{
public:
    string type;
    double n;
    double tao;
};
class defFTF
{
public:
    int order;
    vector3D<double> posref, nref;
    defSection heatsection;
    defFTFModel ftfmodel;
};

class defBC
{
public:
    vector<int> ids;
    string type = "";
    int order = -1;

    string formula = "";
    vector<string> param;

    complex<double> v;
};
// 方程定义
class EqDef
{
public:
    string name;
};
