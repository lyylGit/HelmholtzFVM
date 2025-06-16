#pragma once
#include "CaseBase.h"
#include "complex.h"
#include "MeshOctree.h"
#include <fstream>

using namespace std;

template <class TCase>
class CaseHelmholtz : public CaseBase<TCase>
{
    // 求解Helmholtz方程

private:
    /* data */

    long nCells; // 网格数
    long nFaces; // 面元数

public:
    TCase c0 = 340;      // 声速
    TCase rho_0 = 1.225; // 密度
    // TCase omega;         // 当前计算频率
    TCase freq = 600; // Hz
    // TCase Sx = 1.1, Sy = 0.2, Sz = 0.35; // 声源位置，单位m
    TCase Sx = 0.10, Sy = 0.2, Sz = 0.8;              // 声源位置，单位m
    double q = 0.0001;                                //* 1e3;             // 声源强度
    TCase Mx = 1.2, My = 0.5, Mz = 0.4;               // 检测位置
    TCase freqFrom = 50, freqTo = 400, freqDelta = 5; // 计算模态时使用
public:
    // 场
    string tag_F_P = "Field_P",
           tag_F_Rho = "Field_Density",
           tag_F_Q = "Field_Q",
           tag_E_H = "EQ_Helmholtz";

public:
    CaseHelmholtz(Mesh3<TCase> *_pMesh)
        : CaseBase<TCase>(_pMesh)
    {
        // 网格信息
        nCells = this->_Mesh->nCells;
        nFaces = this->_Mesh->nFaces; // 面元数
    }
    ~CaseHelmholtz() {}

    // 方便使用
    FieldSimple<TCase, complex<TCase>> *getF_P()
    {
        return (FieldSimple<TCase, complex<TCase>> *)(this->GetFieldByName(tag_F_P));
    }
    FieldSimple<TCase, double> *getF_Rho()
    {
        return (FieldSimple<TCase, double> *)(this->GetFieldByName(tag_F_Rho));
    }
    FieldSimple<TCase, double> *getF_Q()
    {
        return (FieldSimple<TCase, double> *)(this->GetFieldByName(tag_F_Q));
    }
    EquationHelmholtz<TCase, complex<TCase>> *getE_H()
    {
        return (EquationHelmholtz<TCase, complex<TCase>> *)(this->GetEqByName(tag_E_H));
    }

public:
    // 重载 定义Case
    void Define()
    {

        this->Name = "Helmholtz方程 计算Case";
        this->Stability = true;  // 计算稳态方程
        this->_bUseChem = false; // 无反应
        this->_DeltaT = 1;       // 1e-3;

        // 复数压力场
        function<void()> _p_updata = bind(&CaseHelmholtz::fn_p_updata, this);

        this->fields[tag_F_P] = make_shared<FieldSimple<TCase, complex<TCase>>>(this, tag_F_P, _p_updata);
        // 密度场，可以假定为常数，此处设置为场
        // 类型为double，不需要更新
        auto Rho = make_shared<FieldSimple<TCase, double>>(this, tag_F_Rho, nullptr);
        this->fields[tag_F_Rho] = Rho;

        // 声源项
        function<void()> _q_updata = bind(&CaseHelmholtz::fn_q_updata, this);
        auto Q = make_shared<FieldSimple<TCase, double>>(this, tag_F_Q, _q_updata);
        this->fields[tag_F_Q] = Q;

        InitFields();

        // 定义方程
        this->equations[tag_E_H] = make_shared<EquationHelmholtz<TCase, complex<TCase>>>(this);
    }

    // 初始化场
    void InitFields()
    {

        // 密度场,设置为恒定值
        auto F = getF_Rho();
        F->SetCellValue(rho_0, true);

        // 压力场
        auto P = getF_P();
        complex<TCase> p0(0, 0);
        P->SetCellValue(p0, true);

        // 声源强度项
        auto Q = getF_Q();

        // 设置初始声压场
        // 声源位置
        vector3D<TCase> QP(Sx, Sy, Sz);
        TCase delta = 0.015;
        auto d1 = 1.0 / (sqrt(2.0 * 3.14159) * delta);
        auto d2 = -1.0 / (2 * delta * delta);
        // lamda
        auto calcuQ = [delta, d1, d2, &QP](vector3D<TCase> &pos) -> TCase
        {
            TCase f = d1 * exp(d2 * (pos - QP).norm2());
            return f;
        };

        // 按面布局，然后按体
        for (long i = 0; i < nFaces; i++)
        {
            auto fcenter = this->_Mesh->faceCenter[i];
            (*Q->fFaceValues)[i] = calcuQ(fcenter);
        }
        Q->fFaceValues->Unit();
        Q->fFaceValues->MultiInplace(q);

        for (long i = 0; i < nCells; i++)
        {
            auto ccenter = this->_Mesh->cellCenter[i];
            (*Q->fValues)[i] = calcuQ(ccenter);
        }
        Q->fValues->Unit();
        Q->fValues->MultiInplace(q);

        //
        // this->_Mesh->template outputTecPlotFEM<double>("d:/qq.plt", "qq", Q->fValues->data);

        // 压力场边界条件，四周均为固壁
        for (long bcID : this->_Mesh->FaceBD)
        {
            shared_ptr<BC<TCase, complex<TCase>>> bc = make_shared<BC<TCase, complex<TCase>>>(P->Name, bcID);
            bc->InitInfo(this->_Mesh->faceSections[bcID]);

            if (bcID == 11)
            {
                // 声阻抗，2阶模型
                shared_ptr<Impedance<TCase, complex<TCase>>> imp = make_shared<Impedance<TCase, complex<TCase>>>(3);

                // 计算系数.直接计算模态时使用
                imp->CalcuImpCoeffs = [](Impedance<TCase, complex<TCase>> *_imp) -> void {
                };
                // 计算具体值
                imp->CalcuImpValues = [this](double fHz, Impedance<TCase, complex<TCase>> *_imp, void *_mem) -> void
                {
                    // 均匀阻抗分布
                    /*  TCase k0 = 2 * PI * fHz / c0;
                      TCase omega = fHz * 2.0 * PI;
                      complex<TCase> Z = rho_0 * c0 / (1.0 + 1.0 / (complex<TCase>(0, 1) * k0)) * 1e4;
                      cout<<Z<<endl;
                      auto gamma = Z / (complex<TCase>(0, 1) * omega * rho_0);

                      auto mem = (MemSegmentBC<TCase, complex<TCase>> *)_mem;
                      mem->SetValues(gamma);
                     */
                    // 穿孔板
                    // Z=ρ0c0[7.337e-3*(1+72.23Ma)+i*2.2245e-5(1+51t)(1+204d)f]/σ
                    TCase Ma = 0;
                    TCase t = 1e-3;       // 穿孔板厚度
                    TCase d = 2e-3;       // 穿孔孔径
                    TCase sigma = 0.0314; // 穿孔率
                    // fHz = 757.0;
                    complex<TCase> Z = rho_0 * c0 / sigma * (7.337e-3 * (1 + 72.23 * Ma) + 1i * 2.2245e-5 * (1 + 51 * t) * (1 + 204 * d) * fHz);
                    // cout << Z << " \t" << abs(Z) << endl;
                    // auto gm = (Z - rho_0 * c0) / (Z + rho_0 * c0);
                    // auto al = 1 - abs(gm*gm);
                    // cout << "alpha=" << al << endl;

                    TCase omega = fHz * 2.0 * PI;
                    auto gamma = Z / (1i * omega * rho_0);
                    // cout << gamma << endl;
                    auto mem = (MemSegmentBC<TCase, complex<TCase>> *)_mem;
                    mem->SetValues(gamma);
                };

                bc->_BCType = BoundaryConditionType::Mixed;
                bc->_Field = P;
                // 边界上：p+Z/(iωρ0)∇(p)=0
                // A=0;
                // B=1;
                // Γ=Z/(iωρ0)
                // TCase k0 = 2 * PI * freq / c0;
                // TCase omega = freq * 2.0 * PI;
                // complex<TCase> Z = rho_0 * c0 / (1.0 + 1.0 / (complex<TCase>(0, 1) * k0)) * 1e4;

                // auto gamma = Z / (complex<TCase>(0, 1) * omega * rho_0);
                bc->SetValueA3(0);
                bc->SetValueB3(1);
                bc->SetValueGamma3(imp);
                this->_BCsNeedCalcu.push_back(bc);
            }
            else
            {
                // bc->_BCType = BoundaryConditionType::Value;
                bc->_BCType = BoundaryConditionType::Gradient;
                bc->SetValue12(complex<TCase>(0, 0));

                bc->_Field = P;
            }
            bc->Init();

            this->ShowPrompt(bc->Name);

            this->_BCs.push_back(bc);
        }
    }
    // 迭代
    bool Iterate()
    {
        // 输出监测点
        // MeshOctree MO(this->_Mesh);
        // auto ip = MO.GetClosetCell(Mx, My, Mz);
        // auto sp = MO.GetClosetCell(Sx, Sy, Sz);
        // 输出结果
        // ofstream sw("d:/freq.txt", ios::out);
        long sp = 11058;
        // long sp = 66632;

        auto eqH = getE_H();
        double k = 162.5;
        // for (int k = freqFrom; k < freqTo; k += freqDelta)
        //  for (int i =0; i < 100; i ++)
        {
            // cout << "iterations=" << k << endl;

            this->freq = k;
            // 更新需要计算的边界阻抗
            for (auto _bcc : this->_BCsNeedCalcu)
            {
                BC<TCase, complex<TCase>> *_bc = (BC<TCase, complex<TCase>> *)(_bcc.get());
                _bc->CalcuImpValues(k);
            }

            eqH->Discrete(0);

            if (!eqH->IterateEquation())
                return false;

            auto P = getF_P();
            // 输出 tecplot

            double *pp = new double[this->_Mesh->nCells];
            for (long i = 0; i < this->_Mesh->nCells; i++)
            {
                pp[i] = abs(P->fValues->data[i]); //.imag();//.real();
            }
            this->_Mesh->template outputTecPlotFEM<double>("d:/pxx_" + to_string(k) + "Hz.plt", "pa", pp);
            delete[] pp;

            auto db = 20.0 * log10(abs(P->fValues->data[sp]) / (2.0e-5));
            cout << sp << "\t freq: " << k << "\tdb=" << db << "\tp=" << P->fValues->data[sp] << endl;

            P->SetCellValue(complex<TCase>(0, 0), true);
        }

        // sw.close();
        return true;
    }

    // 更新函数，压力场
    void fn_p_updata()
    {
    }
    // 更新声源项
    void fn_q_updata()
    {
    }
};