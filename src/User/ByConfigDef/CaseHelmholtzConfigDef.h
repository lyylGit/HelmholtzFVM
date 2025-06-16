#pragma once
#include "CaseBase.h"
#include "complex.h"
#include "ProjectDef.h"

using namespace std;

template <class TCase>
class CaseHelmholtzConfigDef : public CaseBase<TCase>
{
    // 求解Helmholtz方程

public:
    vector3D<TCase> nref = vector3D<TCase>(1.0, 0, 0);       // 参考点流动法向
    long refId = -1;                                         // 参考点
    vector<long> spIds;                                      // 检测点
    vector<string> spOutputfile;                             // 检测点数据
    shared_ptr<FTFSource<TCase, complex<TCase>>> heatSource; // 释热源项部分

public:
    // 场
    string tag_F_P = "Field_P",
           tag_F_Rho = "Field_Density",
           tag_F_Q = "Field_Q",
           tag_E_H = "EQ_Helmholtz";

public:
    CaseHelmholtzConfigDef(Mesh3<TCase> *_pMesh)
        : CaseBase<TCase>(_pMesh)
    {
    }
    ~CaseHelmholtzConfigDef() {}

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
    EquationHelmholtzConfigDef<TCase, complex<TCase>> *getE_H()
    {
        return (EquationHelmholtzConfigDef<TCase, complex<TCase>> *)(this->GetEqByName(tag_E_H));
    }

public:
    // 重载 定义Case
    void Define()
    {
        if (this->_ProjectConfig == NULL)
        {
            this->_Mesh->ShowPrompt("没有定义配置");
            return;
        }

        this->Name = this->_ProjectConfig->caseName;
        this->Stability = true;  // 计算稳态方程
        this->_bUseChem = false; // 无反应
        this->_DeltaT = 1;       // 1e-3;

        // 频率范围:
        this->freqFrom = this->_ProjectConfig->freqFrom;
        this->freqTo = this->_ProjectConfig->freqTo;
        this->freqDelta = this->_ProjectConfig->freqDelta;

        // 复数压力场
        function<void()> _p_updata = bind(&CaseHelmholtzConfigDef::fn_p_updata, this);

        this->fields[tag_F_P] = make_shared<FieldSimple<TCase, complex<TCase>>>(this, tag_F_P, _p_updata);
        // 密度的倒数场，可以假定为常数，此处设置为场
        // 类型为double，不需要更新
        auto Rho = make_shared<FieldSimple<TCase, double>>(this, tag_F_Rho, nullptr);
        this->fields[tag_F_Rho] = Rho;

        // 声源项
        function<void()> _q_updata = bind(&CaseHelmholtzConfigDef::fn_q_updata, this);
        auto Q = make_shared<FieldSimple<TCase, double>>(this, tag_F_Q, _q_updata);
        this->fields[tag_F_Q] = Q;

        InitFields();

        // 定义方程
        this->equations[tag_E_H] = make_shared<EquationHelmholtzConfigDef<TCase, complex<TCase>>>(this);
    }

    // 初始化场
    void InitFields()
    {
        MemSegment<TCase, TCase> rho_f(this, "", DataStoreOn::onCell);
        this->InitField(rho_f, this->_ProjectConfig->rho);
        rho_f.DividedInplace(1.0); // 1.0/rho_f

        auto F = getF_Rho();
        F->SetCellValue(rho_f, true);

        // this->_Mesh->template outputTecPlotFEM<double>(this->_ProjectConfig->outPath + "rr.plt", "cccc", F->fValues->data);

        // 压力场
        // MemSegment<TCase, complex<TCase>> p0(this, "", DataStoreOn::onCell);
        // this->InitField(p0, this->_ProjectConfig->p0);
        auto P = getF_P();
        P->SetCellValue(complex<TCase>(0, 0), true);

        MeshOctree MO(this->_Mesh);

        if (this->_ProjectConfig->helmFreq.outputs.size() > 0)
        {
            for (auto op : this->_ProjectConfig->helmFreq.outputs)
            {
                spIds.push_back(MO.GetClosetCell(op->mic));
                spOutputfile.push_back(op->filename);
            }
        }

        // 声源项
        auto Q = getF_Q();
        this->InitSources(Q);
        // 释热源项：ftf
        this->InitHeatSources(MO, heatSource, nref, refId);

        //
        // 边界条件
        for (long bcID : this->_Mesh->FaceBD)
        {
            shared_ptr<BC<TCase, complex<TCase>>> bc = make_shared<BC<TCase, complex<TCase>>>(P->Name, bcID);
            bc->InitInfo(this->_Mesh->faceSections[bcID]);

            auto def = this->getBCDef(bcID);
            if (def != NULL)
            {
                this->BuildBCbyDef(bc, def, F);
            }
            bc->_Field = P;

            bc->Init();

            this->ShowPrompt(bc->Name);

            this->_BCs.push_back(bc);
        }
    }
    // 迭代
    bool Iterate()
    {
        auto eqH = getE_H();

        if (this->_ProjectConfig->caseType == defCaseType::HelmFreq)
        {

            for (int k = this->freqFrom; k < this->freqTo; k += this->freqDelta)
            {
                cout << "iterations freq=" << k << endl;

                this->freqCurrent = k;
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
                    // auto real= P->fValues->data[i].real();
                    int m = P->fValues->data[i].real() >= 0 ? 1 : -1;
                    pp[i] = abs(P->fValues->data[i]) * m; //*real/abs(real);
                }
                this->_Mesh->template outputTecPlotFEM<double>(this->_ProjectConfig->outPath + "pxx_" + to_string(k) + "Hz.plt", "pa", pp);
                delete[] pp;
                for (int kk = 0; kk < spIds.size(); kk++)
                {
                    auto sp = spIds[kk];
                    auto db = 20.0 * log10(abs(P->fValues->data[sp]) / (2.0e-5));
                    cout << sp << "\t freq : " << k << "\t db = " << db << "\t p = " << P->fValues->data[sp] << endl;
                    ofstream sw(this->_ProjectConfig->outPath + spOutputfile[kk], ios::out | ios::app);
                    sw << k << "\t" << db << endl;
                    sw.close();
                }

                P->SetCellValue(complex<TCase>(0, 0), true);
            }
        }
        else
            eqH->Discrete(0);

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