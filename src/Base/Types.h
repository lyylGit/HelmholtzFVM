#pragma once

#define PI 3.141592653589793238462643383279502884197169399375105820974944592307816406L

enum class VariableType
{
    PhysicalProperty, //物性参数
    Variable,         //待求解变量
    CAL,              //需要根据公式进行计算
    NUM,              //常数
};

//数据类型,标量,矢量或张量
enum  class FieldType
{
    _Scalar, // 1
    _Vector, // 1*3
    _Tensor, // 3*3
};

enum  class DiscreteMethodSchemeType
{
    //时间方向离散
    Explicit,       //显式
    Implicit,       //全隐式
    CrankNicholson, // C-N混合格式,

    //面值离散
    Central, //中心差分格式
    Upwind,  //上风格式

};
enum  class FaceValueSchemeType
{
    //面值离散
    Central,   //中心差分格式 
    Upwind,    //上风格式
    Corrected, // Rhie-Chow 修正
};
enum  class DiscreteSchemeTransient
{
    BackwardEuler, //一阶
    AdamsMoulton,  //二阶
};
//对网格控制体梯度离散方法
enum  class DiscreteVolumeGradSchemeType
{
    //高斯定理
    GreenGauss,
    //最小二乘
    LeastSquare,
};
//对面法向梯度的离散方法
enum  class DiscreteFaceGradSchemeType
{
    //无修正
    NoCorrection,
    //非正交修正
    NonOrthogonalityCorrection,
    //倾斜修正
    SkewnessCorrection
};
enum class DataStoreOn
{
    //存储于点
    onVertex,
    //存储于面上
    onFace,
    //存储于体上
    onCell,
};

// 方程离散方法定义
class EquationConfigue
{
public:
    // 时间步长
    double _DeltaT;
    // 时间步长
    double _DeltaTao;

    // 松弛因子
    double _Sor = 1e-2;

    // 离散选项
    // 瞬态项
    DiscreteSchemeTransient dmsTransient;
    // 时间方向离散格式
    DiscreteMethodSchemeType dmsTemporal;
    // 交界面计算
    DiscreteMethodSchemeType dmsFaceValue;
    // 对流项交界面值
    DiscreteMethodSchemeType dmsFaceValueConvection;
    // Grad计算
    DiscreteVolumeGradSchemeType dmsGrad;
    // 面法向梯度计算
    DiscreteFaceGradSchemeType dmsFaceGrad;
    // 源项计算
    DiscreteMethodSchemeType dmsSource;
};


enum class BoundaryConditionType
{
    // 第一类边界条件,给定值, Dirichlet
    Value,
    // 第二类边界条件,给定梯度, Neumann
    Gradient,
    //第三类边界条件:边界上 A+Bϕ+Γ∇ϕ=0
    Mixed,
    //边界是阻抗，且是特征值问题
    MixedModal,

    // 单元格的边界信息仅为处理方便
    // 单元格为内部,
    // Inner,
    // 单元格有多个边界
    // Multi,

    // 使用center,有些中间变量在进行通量计算时,没有边界条件,直接采用单元中心值
    UseCenter,
    //
    // 未定义
    Null,
};

enum class PhysicsBCType
{
    Wall,     // 固壁
    Inlet,    // 进口
    Outlet,   // 出口
    Symmetry, // 对称

    // 单元格为内部,
    Inner,
    // 单元格有多个边界
    Multi,
};
// 面属性,与Fluent一致
enum class FaceType
{
    Interior = 2,         // 内部面
    Wall = 3,             // 固壁
    PressurInlet = 4,     // pressure-inlet, inlet-vent, intake-fan
    PressureOutlet = 5,   // pressure-outlet, exhaust-fan, outlet-vent
    Symmetry = 7,         // symmetry
    PeriodicShadow = 8,   // periodic-shadow
    PressureFarField = 9, // pressure-far-field
    VelocityInlet = 10,   // velocity-inlet
    Periodic = 12,        // periodic
    Fan = 14,             // fan, porous-jump, radiator
    MassFlowInlet = 20,   // mass-flow-inlet
    Interface = 24,       // interface
    Parent = 31,          // parent(hanging node)
    Outflow = 36,         // outflow
    Axis = 37,            // axis
};