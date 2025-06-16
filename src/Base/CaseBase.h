#pragma once

#include <memory>
#include <functional>
#include "staticFunctions.h"
#include "Mesh3.h"
#include "PhysicsBC.h"
#include "IName.h"
#include "Types.h"
#include "ProjectDef.h"
#include "vector3d.h"
#include <functional>

using namespace std;

template <class TMesh, class TField>
class FieldBase;
template <class TMesh, class TField>
class BC;
template <class TMesh, class TField>
class EquationBase;

template <class TMesh, class TField>
class FTFSource;

template <class TMesh>
class MeshOctree;

template <class TMesh, class TField>
class MemSegment;

template <class T>
class PolyFit;

template <class TMesh>
class CaseBase : public staticFunctions
{
private:
public:
    // 网格
    Mesh3<TMesh> *_Mesh;
    // 名称
    string Name;
    // 是否稳态
    bool Stability = true;

    // 多组分的内容
    // 是否使用Chemicallib
    bool _bUseChem = false;
    // 如果使用chemicalib
    // 热力学数据文件，输运数据，反应机理文件
    string _fnThermoData = "";
    string _fnTransData = "";
    string _fnReactMech = "";
    // 混合物构成,最后一个组分不计算，因此最后一个组分通常为大组分
    vector<string> _mixtureNames;
    // ClassMixtures _mixtures = null;
    // //组分
    // map<string, ClassSpecies> _htSpecies = null;
    // //反应
    // ClassReactions _reactionPath = null;
    // 时间控制
    double _DeltaT = 1e-5;
    double _DeltaTao = 1;

    // 相
    vector<string> FluidPhase,
        SolidPhase;
    // 场
    // map<string, FieldBase<TMesh,TField> *> fields;
    map<string, shared_ptr<void>> fields;
    // 组分场索引
    //  map<string, FLOAT.Concentration> htYFields;
    //  map<string, FLOAT.Wdot> htWdot ;

    // 方程
    //  map<string, EquationBase<TMesh,TField>*> equations ;
    map<string, shared_ptr<void>> equations;
    // //组分方程
    // vector<EquationBase> YEquations ;
    // 边界条件
    vector<shared_ptr<BCBase>> _BCs;
    // 需要计算的边界条件
    vector<shared_ptr<BCBase>> _BCsNeedCalcu;
    // 物理边界条件
    vector<shared_ptr<PhysicsBCBase>> PhysicsBCs;

    // 事件
    //  vector<Event> Events;

    // 配置
    ProjectConfig *_ProjectConfig = NULL;

    // TMesh freqFrom = 50, freqTo = 400, freqDelta = 5; // 计算模态时使用，拟合时使用
    TMesh freqFrom = 50, freqTo = 400, freqDelta = 5; // 计算模态时使用，拟合时使用
    TMesh freqCurrent = 0;                            // 在模式为HelmFreq时，当前计算频率

    long nCells; // 网格数
    long nFaces; // 面元数

public:
    CaseBase(Mesh3<TMesh> *_pMesh)
    {
        _Mesh = _pMesh;
        nCells = _Mesh->nCells;
        nFaces = _Mesh->nFaces; // 面元数
    }
    ~CaseBase() {};

public:
    // 根据名称获取场Field： void *
    void *GetFieldByName(string tag)
    {
        if (fields.count(tag) > 0)
        {
            shared_ptr<void> pp = this->fields[tag];
            return pp.get();
        }
        else
            return nullptr;
    }

    // 根据名称获取方程Eq： void *
    void *GetEqByName(string tag)
    {
        shared_ptr<void> pp = this->equations[tag];
        return pp.get();
    }
    // 获取该Case的f场的边界条件
    template <class TField>
    vector<shared_ptr<void>> GetBoundaryConditions(FieldBase<TMesh, TField> *f)
    {
        vector<shared_ptr<void>> ret;

        size_t nbcs = _BCs.size();
        for (size_t j = 0; j < nbcs; j++)
        {
            auto bc = _BCs[j];
            if (bc->Name == f->Name)
                ret.push_back(bc);
        }
        return ret;
    }

    // 获取该Case的f场的bcID的边界条件
    template <class TField>
    void *GetBoundaryConditions(const FieldBase<TMesh, TField> *f, long bcID)
    {
        size_t nbcs = _BCs.size();
        for (size_t j = 0; j < nbcs; j++)
        {
            auto bc = _BCs[j];
            if (bc->Name == f->Name && bc->ID == bcID)
            {
                return bc.get();
            }
        }
        return nullptr;
    }
    template <class TField>
    void *GetBoundaryConditionOfFace(FieldBase<TMesh, TField> *f, long faceID)
    {
        // 首先确定faceID所在的bcID
        long bcID = -1;
        map<int, tuple<long, long, int, int>>::iterator it = _Mesh->faceSections.begin();
        map<int, tuple<long, long, int, int>>::iterator itEnd = _Mesh->faceSections.end();
        while (it != itEnd)
        {
            int id = it->first;
            auto fs = it->second;
            int istart = get<0>(fs), iend = get<1>(fs);
            if (faceID >= istart && faceID <= iend)
            {
                bcID = id;
                break;
            }
        }

        return GetBoundaryConditions(f, bcID);
    }
    // 物理边界
    void *GetBC(long bcID)
    {
        size_t nbcs = PhysicsBCs.size();
        for (size_t j = 0; j < nbcs; j++)
        {
            auto bc = dynamic_pointer_cast<PhysicsBCBase>(PhysicsBCs[j]);

            if (bc->ID == bcID)
                return bc.get();
        }

        return NULL;
    }

    defBC *getBCDef(long bcID)
    {
        defBC *def = NULL;
        defBC *normal = NULL;
        for (int k = 0; k < _ProjectConfig->bcs.size(); k++)
        {
            auto bc = _ProjectConfig->bcs[k];
            if (find(bc->ids.begin(), bc->ids.end(), bcID) != bc->ids.end())
            {
                def = bc.get();
                break;
            }
            else if (bc->ids.size() <= 0)
                normal = bc.get();
        }

        if (def != NULL)
            return def;
        else
            return normal;
    }
    template <class TField>
    void BuildBCbyDef(shared_ptr<BC<TMesh, TField>> bc, defBC *def, FieldBase<TMesh, TMesh> *rho1)
    {
        if (def->type == "MixedModal")
        {
            long bfidFrom = bc->start;
            long bfidTo = bc->end;
            shared_ptr<Impedance<TMesh, TField>> imp = make_shared<Impedance<TMesh, TField>>(def->order + 1, bc->length, bc->start);

            // 计算系数
            imp->CalcuImpCoeffs = [this, rho1, def, bfidFrom, bfidTo](Impedance<TMesh, TField> *_imp) -> void
            {
                stopWatch();
                double minRSqure = 1e30;
                // 计算拟合用数据对
                int npts = (int)((freqTo - freqFrom) / freqDelta + 1);
                TField *x = new TField[npts];
                TField *zz = new TField[npts];
                {
                    int i = 0;
                    // 定值
                    if (def->param.size() <= 0)
                        for (TMesh fHz = freqFrom; fHz <= freqTo; fHz += freqDelta)
                        {
                            TMesh omg = 2 * PI * fHz;
                            zz[i] = def->v; // 采用Z
                            x[i] = omg;
                            i++;
                        }
                    else
                    {
                        ParserX parser;
                        parser.SetExpr(def->formula.c_str());
                        Value aHz(freqFrom);
                        parser.DefineVar("Hz", Variable(&aHz));
                        
                        for (TMesh fHz = freqFrom; fHz <= freqTo; fHz += freqDelta)
                        {
                            aHz = fHz;
                            TMesh omg = 2 * PI * fHz;
                            zz[i] = parser.Eval().GetComplex(); // 采用Z
                            x[i] = omg;
                            i++;
                        }
                    }
                }
                // #pragma omp parallel for
                for (long faceId = bfidFrom; faceId <= bfidTo; faceId++)
                {
                    auto rho_1 = 1.0 / (*rho1->fFaceValues)[faceId];
                    TField *y = new TField[npts];
                    int i = 0;
                    for (TMesh fHz = freqFrom; fHz <= freqTo; fHz += freqDelta)
                    {
                        auto omg = x[i];
                        TField Z = zz[i];
                        TMesh delta = this->_Mesh->dcf[faceId];
                        y[i] = omg * rho_1 / (1i * Z - omg * rho_1 * delta);
                        i++;
                    }
                    // 拟合
                    PolyFit fit(_imp->N - 1, npts, x, y);
                    fit.Fit();

                    // 系数送入imp中
                    for (int j = 0; j < _imp->N; j++)
                    {
                        (_imp->A[j])[faceId] = fit.coefbeta[j];
                    }

                    minRSqure = min(minRSqure, fit.R2);
                    delete[] y;
                }
                // sw.close();
                delete[] x;
                delete[] zz;
                cout << "imp:" << minRSqure << endl;
                cout << stopWatch() << endl;
            };
            // 计算具体值,计算模态时不需要实现
            imp->CalcuImpValues = [this, rho1, def, bfidFrom, bfidTo](double fHz, Impedance<TMesh, TField> *_imp, void *_mem) -> void
            {
                TField Z = 0;
                // 定值
                if (def->param.size() <= 0)
                    Z = def->v; // 采用Z
                else
                {
                    ParserX parser;
                    parser.SetExpr(def->formula.c_str());
                    Value aHz(fHz);
                    parser.DefineVar("Hz", Variable(&aHz));

                    Z = parser.Eval().GetComplex(); // 采用Z
                }

                auto mem = (MemSegmentBC<TMesh, TField> *)_mem;
                for (long faceId = bfidFrom; faceId <= bfidTo; faceId++)
                {
                    auto rho_1 = 1.0 / (*rho1->fFaceValues)[faceId];
                    TMesh omega = fHz * 2.0 * PI;
                    (*mem)[faceId] = Z / (1i * omega * rho_1);
                }
            };
            if (_ProjectConfig->caseType == defCaseType::HelmModal)
            {
                // 现在就可以计算系数
                imp->CalcuImpCoeffs(imp.get());
                // 计算模态的问题
                bc->_BCType = BoundaryConditionType::MixedModal;
            }
            else // 扫频计算时，此为定值，因此为mixed
            {
                bc->_BCType = BoundaryConditionType::Mixed;
                bc->SetValueA3(0);
                bc->SetValueB3(1);
            }
            bc->SetValueGamma3(imp);
            this->_BCsNeedCalcu.push_back(bc);
        }
        else if (def->type == "Mixed")
        {
            bc->_BCType = BoundaryConditionType::Mixed;
            throw "not finished";

            /* TMesh k0 = 2 * PI * freq / c0;
             TMesh omega = freq * 2.0 * PI;
             complex<TMesh> Z = rho_0 * c0 / (1.0 + 1.0 / (complex<TMesh>(0, 1) * k0)); //* 1e4;

             auto gamma = Z / (complex<TMesh>(0, 1) * omega * rho_0);
             bc->SetValueA3(0);
             bc->SetValueB3(1);
             bc->SetValueGamma3(gamma);*/
        }
        else if (def->type == "Gradient")
        {
            bc->_BCType = BoundaryConditionType::Gradient;
            bc->SetValue12(def->v);
        }
        else if (def->type == "Value")
        {
            bc->_BCType = BoundaryConditionType::Value;
            bc->SetValue12(def->v);
        }
        else
            throw "边界类型 " + def->type + " 未定义";
    }

    // 通过config定义场
    template <class TField>
    void InitField(MemSegment<TMesh, TField> &f, defFieldValueDefine &def)
    {

        if (def.type == defFieldType::Uniform)
        {
            f.SetValues(def.value);
        }
        else if (def.type == defFieldType::Field)
        {
            // 从文件读取
            if (def.filename != "")
            {
            }
            else if (def.sections.size() > 0) // 分段设置
            {
                for (int k = 0; k < def.sections.size(); k++)
                {
                    auto sec = def.sections[k];
                    InitFieldSegments(f.data.get(), *(sec.get()));
                }
            }
        }
    }
    template <class TField>
    void InitFieldSegments(TField *f, defSection &sec)
    {
        long nCells = _Mesh->nCells;

        function<void(TField &, vector3D<TMesh> &)> v;
        string dim = sec.dimension;

        for (int m = 0; m < sec.segments.size(); m++)
        {
            auto seg = sec.segments[m];
            double from = get<0>(seg), to = get<1>(seg), val = get<4>(seg);
            int a = get<2>(seg), b = get<3>(seg);

            if (dim == "x")
            {
                if (a == 0)
                {
                    if (b == 0)
                        v = [this, from, to, val](TField &fv, vector3D<TMesh> &cell) -> void
                        {
                            if (cell.x > from && cell.x < to)
                                fv = val;
                        };
                    else if (b == 1)
                        v = [this, from, to, val](TField &fv, vector3D<TMesh> &cell) -> void
                        {
                            if (cell.x > from && cell.x <= to)
                                fv = val;
                        };
                }
                else if (a == 1)
                {
                    if (b == 0)
                        v = [this, from, to, val](TField &fv, vector3D<TMesh> &cell) -> void
                        {
                            if (cell.x > from && cell.x <= to)
                                fv = val;
                        };
                    else if (b == 1)
                        v = [this, from, to, val](TField &fv, vector3D<TMesh> &cell) -> void
                        {
                            if (cell.x >= from && cell.x <= to)
                                fv = val;
                        };
                }
            }
            else if (dim == "y")
            {
                if (a == 0)
                {
                    if (b == 0)
                        v = [this, from, to, val](TField &fv, vector3D<TMesh> &cell) -> void
                        {
                            if (cell.z > from && cell.z < to)
                                fv = val;
                        };
                    else if (b == 1)
                        v = [this, from, to, val](TField &fv, vector3D<TMesh> &cell) -> void
                        {
                            if (cell.z > from && cell.z <= to)
                                fv = val;
                        };
                }
                else if (a == 1)
                {
                    if (b == 0)
                        v = [this, from, to, val](TField &fv, vector3D<TMesh> &cell) -> void
                        {
                            if (cell.z > from && cell.z <= to)
                                fv = val;
                        };
                    else if (b == 1)
                        v = [this, from, to, val](TField &fv, vector3D<TMesh> &cell) -> void
                        {
                            if (cell.z >= from && cell.z <= to)
                                fv = val;
                        };
                }
            }
            else if (dim == "z")
            {
                if (a == 0)
                {
                    if (b == 0)
                        v = [this, from, to, val](TField &fv, vector3D<TMesh> &cell) -> void
                        {
                            if (cell.y > from && cell.y < to)
                                fv = val;
                        };
                    else if (b == 1)
                        v = [this, from, to, val](TField &fv, vector3D<TMesh> &cell) -> void
                        {
                            if (cell.y > from && cell.y <= to)
                                fv = val;
                        };
                }
                else if (a == 1)
                {
                    if (b == 0)
                        v = [this, from, to, val](TField &fv, vector3D<TMesh> &cell) -> void
                        {
                            if (cell.y > from && cell.y <= to)
                                fv = val;
                        };
                    else if (b == 1)
                        v = [this, from, to, val](TField &fv, vector3D<TMesh> &cell) -> void
                        {
                            if (cell.y >= from && cell.y <= to)
                                fv = val;
                        };
                }
            }
            for (long i = 0; i < nCells; i++)
            {
                v(f[i], _Mesh->cellCenter[i]);
            }
        }
    }
    template <class TField>
    void InitSources(FieldBase<TMesh, TField> *Q)
    {
        if (this->_ProjectConfig->caseType == defCaseType::HelmFreq && this->_ProjectConfig->helmFreq.sources.size() > 0)
        {
            for (auto source : this->_ProjectConfig->helmFreq.sources)
            {
                // 用高斯分布，暂时只支持一个。
                if (source->type == "point")
                {
                    TMesh delta = 0.015;
                    auto d1 = 1.0 / (sqrt(2.0 * PI) * delta);
                    auto d2 = -1.0 / (2 * delta * delta);
                    // lamda
                    auto calcuQ = [delta, d1, d2, &source](vector3D<TMesh> &pos) -> TMesh
                    {
                        TMesh f = d1 * exp(d2 * (pos - source->pos).norm2());
                        return f;
                    };

                    // 按面布局，然后按体
                    for (long i = 0; i < nFaces; i++)
                    {
                        auto fcenter = this->_Mesh->faceCenter[i];
                        (*Q->fFaceValues)[i] = calcuQ(fcenter);
                    }
                    Q->fFaceValues->Unit();
                    Q->fFaceValues->MultiInplace(source->q);

                    for (long i = 0; i < nCells; i++)
                    {
                        auto ccenter = this->_Mesh->cellCenter[i];
                        (*Q->fValues)[i] = calcuQ(ccenter);
                    }
                    Q->fValues->Unit();
                    Q->fValues->MultiInplace(source->q);
                }
            }
        }
    }

    void InitHeatSources(MeshOctree<TMesh> &MO, shared_ptr<FTFSource<TMesh, complex<TMesh>>> &heatSource, vector3D<TMesh> &nref, long &refId)
    {
        // 释热,暂支持一个响应
        if (this->_ProjectConfig->ftfs.size() > 0)
        {
            for (int i = 0; i < this->_ProjectConfig->ftfs.size(); i++)
            {
                auto def = this->_ProjectConfig->ftfs[i];
                auto posref = def->posref;
                nref = def->nref;

                // 参考点的cell id
                refId = MO.GetClosetCell(posref);

                // 释热源项
                heatSource = make_shared<FTFSource<TMesh, complex<TMesh>>>(def->order + 1, this); // 2阶拟合
                heatSource->CalcuFTFSourceCoeffs = [this, def](FTFSource<TMesh, complex<TMesh>> *ftf) -> void
                {
                    stopWatch();
                    double minRSqure = 1e30;

                    // n-τ模型
                    TMesh Nux = def->ftfmodel.n;   // 增益
                    TMesh tao = def->ftfmodel.tao; // 1e-4;  // 时滞 s

                    // 计算拟合用数据对,假设与空间无关
                    int npts = (int)((this->freqTo - this->freqFrom) / this->freqDelta + 1);
                    complex<TMesh> *x = new complex<TMesh>[npts];
                    complex<TMesh> *y = new complex<TMesh>[npts];
                    {
                        int i = 0;
                        for (TMesh fHz = this->freqFrom; fHz <= this->freqTo; fHz += this->freqDelta)
                        {
                            TMesh omg = 2 * PI * fHz;
                            complex<TMesh> q = Nux * exp(1i * omg * tao);
                            y[i] = q;
                            x[i] = omg;
                            i++;
                        }
                    }

                    // 拟合
                    PolyFit fit(ftf->N - 1, npts, x, y);
                    fit.Fit();

                    minRSqure = min(minRSqure, fit.R2);

                    delete[] y;
                    delete[] x;

                    TMesh *qq = new TMesh[nCells];
                    this->InitFieldSegments(qq, def->heatsection);

                    for (long i = 0; i < nCells; i++)
                    {
                        // 在火焰区赋值
                        // 系数送入FTF中
                        for (int j = 0; j < ftf->N; j++)
                        {
                            (ftf->A[j])[i] = fit.coefbeta[j] * qq[i];
                        }
                    }
                    // this->_Mesh->template outputTecPlotFEM<double>(this->_ProjectConfig->outPath + "pp.plt", "cccc", qq);
                    delete[] qq;

                    cout << "q:" << minRSqure << endl;
                    cout << stopWatch() << endl;
                };
                heatSource->CalcuFTFSourceValues = [this, def](double fHz, FTFSource<TMesh, complex<TMesh>> *ftf, void *_mem) -> void
                {
                    // n-τ模型
                    TMesh Nux = def->ftfmodel.n;   // 增益
                    TMesh tao = def->ftfmodel.tao; // 1e-4;  // 时滞 s

                    TMesh omg = 2 * PI * fHz;
                    complex<TMesh> q = Nux * exp(1i * omg * tao);

                    MemSegment<TMesh, complex<TMesh>> *mem = (MemSegment<TMesh, complex<TMesh>> *)_mem;
                    complex<TMesh> *qq = new complex<TMesh>[nCells];
                    this->InitFieldSegments(qq, def->heatsection);

                    for (long i = 0; i < nCells; i++)
                    {
                        qq[i] *= q;
                    }
                    mem->SetValues(qq);
                    // this->_Mesh->template outputTecPlotFEM<double>(this->_ProjectConfig->outPath + "pp.plt", "cccc", qq);
                    delete[] qq;
                };
                // 计算释热源项的拟合系数
                if (this->_ProjectConfig->caseType == defCaseType::HelmModal)
                    heatSource->CalcuFTFSourceCoeffs(heatSource.get());
            }
        }
    }

public:
    // 定义Case,由子类实现
    virtual void Define() {};
    // 迭代求解
    virtual bool Iterate() { return true; };
};
