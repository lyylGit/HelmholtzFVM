#pragma once

#include "staticFunctions.h"
#include <map>
#include <set>

#ifdef MPI
#include <metis.h>
#endif

// 为了方便，需要交换的cellid进行偏移
#define CELLID_SHIFT(X) (-X - 111)
// #include "Cell.h"
using namespace std;

// 网格单元类型
enum elementType
{
    mixed = 0,
    triangular = 1,
    tetrahedral = 2,
    quadrilateral = 3,
    hexahedral = 4,
    pyramid = 5,
    wedge = 6,
    polyhedral
};
// 节点信息
template <class T>
struct VertexInfo
{
    // 围绕此节点的cell编号,
    // 该面元分摊的距离到此节点
    map<long, T> distPart;

    // 总长度
    T totalDist = 0;
    T D()
    {
        if (totalDist == 0)
        {
            for (auto part : distPart)
                totalDist += part.second;
        }
        return totalDist;
    }
    // v是按cell编号构成的内存
    T calcuValue(T *v)
    {
        T ret = 0;
        D();
        for (auto part : distPart)
            ret += part.second / totalDist * v[part.first];
        return ret;
    }
};

// extern double m_NearZero ;

template <class T>
class Mesh3 : public staticFunctions
{
public:
    // 网格数据:以下内容为读取
    // 维数：3
    int nDimensions;
    // 节点数,面元数,网格单元数
    long nNodes, nFaces, nCells;
    // 节点，三维坐标
    shared_ptr<vector3D<T>[]> _Nodes;
    // 面元 [nFaces][pt]
    shared_ptr<vector<long>[]> _Faces;
    // 网格单元的面组成
    shared_ptr<vector<long>[]> _Cells;
    // 网格单元的点组成，输出tecplot用；
    // shared_ptr<set<long>[]> _CellPts;
    // 该面元的前向和后向单元序号
    shared_ptr<long[]> neighborCell; // F-->prevCell
    shared_ptr<long[]> ownerCell;    // C  -->nextCell

    // 以下内容为计算

    // 面元面积矢量
    shared_ptr<vector3D<T>[]> surfaceVector;
    // 面元面积大小
    shared_ptr<T[]> faceArea;
    // 面元中心 nFaces*3
    shared_ptr<vector3D<T>[]> faceCenter;
    // Cell中心 nCells*3
    shared_ptr<vector3D<T>[]> cellCenter;
    shared_ptr<T[]> cellVolume; // 网格体积
    // cell中每个面的法向系数；
    shared_ptr<map<long, short>[]> cellFaceSign;
    // 每个面的ownerCell的法向符号
    //  shared_ptr<short[]> ownerFaceSign;
    //  边界面到其网格中心的距离
    shared_ptr<T[]> gc;            // dFf/dFC
    shared_ptr<vector3D<T>[]> ECF; // 单位矢量C->F,
    shared_ptr<T[]> dcf;           // 距离

    // 记录那些是内部面、边界面。
    vector<long> inFaces, bdFaces;
    // 分区时交界面,两侧单元需要交换数据
    vector<long> interFaces;
    // 需要交换的cell:partionID, local face id, global cell id
    map<int, map<long, long>> RecivCells;
    map<int, vector<long>> SendCells;

    // 面集合的id,from,to,bcType,faceType;
    map<int, tuple<long, long, int, int>> faceSections;
    vector<int> faceIceSections;
    // cell的id,from,to,bctype,facetype;
    map<int, tuple<long, long, int, int>> cellSections;
    map<int, string> zoneSections;

    // 边界面: from--to--id: 从---到--groupid
    // id 对应到faceSection 中的id
    vector<int> FaceBD, FaceInner;

    // 网格倍数,程序计算采用m，
    // 如果网格尺寸是mm，则需要乘以倍数: mm=1e-3, mm2=mm*mm;
    T mm = 1;

    const long NReadBytes = 1024 * 5; // 一次读取的数据长度
                                      // 用户需要读取的数据
    map<string, T *> cellData;        // 网格数据
    map<string, T *> faceData;        // 面元数据
    // 输出按点的tecplot文件使用。
    shared_ptr<VertexInfo<T>[]> vertex2Cell;

private:
    vector<string> vertexStructuresFEM, cellStructuresFEM, vertexStructuresCellCenter; // tecplot 输出时的网格结构

private:
    // T mm2 = 1;
    void Init()
    {
        Clear();

        nDimensions = -1;
        nNodes = -1;
        nFaces = -1;
        nCells = -1;

        // mm2 = mm * mm; // 面积
    }
    void Clear()
    {
        faceSections.clear();
        zoneSections.clear();
    }

public:
    Mesh3(/* args */)
    {
    }
    // 按面切分时产生的子网格
    Mesh3(long nfaces)
    {
        nDimensions = 3;
        nFaces = nfaces;
        _Faces = make_shared<vector<long>[]>(nFaces);
        neighborCell = make_shared<long[]>(nFaces);
        ownerCell = make_shared<long[]>(nFaces);
    }
    ~Mesh3() {}

    void setFacesNum(long nfaces)
    {
        nDimensions = 3;
        nFaces = nfaces;
        _Faces = make_shared<vector<long>[]>(nFaces);
        neighborCell = make_shared<long[]>(nFaces);
        ownerCell = make_shared<long[]>(nFaces);
    }

public:
    // 读取msh文本格式的网格信息
    bool ReadMesh_Fluent_Mesh_Text(string mshFileName)
    {
        ifstream sr(mshFileName, ios::in);
        if (!sr.good())
        {
            m_strErr = mshFileName + " 文件不存在";
            return false;
        }

        // 信息复位
        Init();

        // 格式的标志
        string tagDim = "(2 ",
               tagNode = "(10 ",
               tagFace = "(13 ",
               tagCell = "(12 ",
               tagZoneTye = "(45 ",
               tagZoneSection = "(39 ";
        bool bDim = false,
             bNode = false,
             bFace = false,
             bCell = false;
        bool bSuccess = true;

        string line, s;
        while (getline(sr, line))
        {
            if (line == "")
                continue;
            line = Trim(line);
            if (line == "")
                continue;
            if (!bDim && isTarget(line, tagDim))
            {
                s = getContentString(line, tagDim);
                nDimensions = o2i(s);
                if (nDimensions != 3)
                {
                    m_strErr = "不是三维问题,暂不能处理";
                    bSuccess = false;
                    break;
                }
                bDim = true;
                ShowPrompt("维数:" + to_string(nDimensions));
                continue;
            }
            if (!bNode && isTarget(line, tagNode))
            {
                s = getContentString(line, tagNode);
                s = s.substr(1);
                s = s.substr(0, s.length() - 1);
                auto sa = split(s, " ");
                if (sa.size() >= 4)
                {
                    nNodes = stoi(sa[2], 0, 16);
                    _Nodes = make_shared<vector3D<T>[]>(nNodes);
                }
                bNode = true;

                ShowPrompt("节点数:" + to_string(nNodes));
                continue;
            }
            if (!bFace && isTarget(line, tagFace))
            {
                s = getContentString(line, tagFace);
                s = s.substr(1);
                s = s.substr(0, s.length() - 1);
                auto sa = split(s, " ");
                if (sa.size() >= 4)
                {
                    nFaces = stoi(sa[2], 0, 16);
                    _Faces = make_shared<vector<long>[]>(nFaces);
                    neighborCell = make_shared<long[]>(nFaces);
                    ownerCell = make_shared<long[]>(nFaces);
                }
                bFace = true;

                ShowPrompt("面元数:" + to_string(nFaces));
                continue;
            }
            if (!bCell && isTarget(line, tagCell))
            {
                s = getContentString(line, tagCell);
                s = s.substr(1);
                s = s.substr(0, s.length() - 1);
                auto sa = split(s, " ");
                if (sa.size() >= 4)
                {
                    nCells = stoi(sa[2], 0, 16);
                }
                bCell = true;

                ShowPrompt("网格数:" + to_string(nCells));
                continue;
            }
            // 节点,可能会有多段
            if (isTarget(line, tagNode))
            {
                s = getContentString(line, tagNode);
                s = s.substr(1);
                s = s.substr(0, s.length() - 1);
                auto sa = split(s, " ");
                if (sa.size() >= 4)
                {
                    int istart = stoi(sa[1], 0, 16) - 1;
                    int iend = stoi(sa[2], 0, 16) - 1;

                    ShowPrompt("正在读取 节点:" + s);
                    // 读取
                    for (int j = istart; j <= iend; j++)
                    {
                        // 每行是3个数

                        getline(sr, line);
                        line = Trim(line);
                        if (line == "")
                        {
                            j--;
                            continue;
                        }

                        auto f = getMathValueDouble(line);
                        if (f.size() < 3)
                        {
                            j--;
                            continue;
                        }
                        _Nodes[j] = vector3D<T>(f[0] * mm, f[1] * mm, f[2] * mm);
                    }
                }
                continue;
            }
            // 面
            if (isTarget(line, tagFace))
            {
                s = getContentString(line, tagFace);
                s = s.substr(1);
                s = s.substr(0, s.length() - 1);
                auto sa = split(s, " ");
                if (sa.size() >= 5)
                {
                    ShowPrompt("正在读取 面元:" + s);
                    int id = stoi(sa[0], 0, 16);
                    long istart = stol(sa[1], 0, 16) - 1;
                    long iend = stol(sa[2], 0, 16) - 1;
                    int bcType = stoi(sa[3], 0, 16);
                    int facetype = stoi(sa[4]);
                    // 读取
                    if (facetype == 0)
                    {
                        ShowPrompt("混合网格id:" + to_string(id));
                    }
                    if (facetype == 5)
                    {
                        ShowPrompt("多边形网格id:" + to_string(id));
                    }
                    if (bcType != 2)
                    {
                        ShowPrompt("边界id:" + to_string(id));
                        FaceBD.push_back(id);
                    }
                    else
                        FaceInner.push_back(id);
                    // 存储面集合信息
                    faceSections[id] = make_tuple(istart, iend, bcType, facetype);

                    for (long j = istart; j <= iend; j++)
                    {
                        getline(sr, line);
                        line = Trim(line);
                        if (line == "")
                        {
                            j--;
                            continue;
                        }

                        auto f = getMathValueHex(line);
                        // 一行内容为: pt1 pt2 pt3 .... prev next
                        // 如果是混合或多边形,则前面多一个节点数:nd pt1 pt2 pt3 .... prev next
                        //  I为起点
                        int I = 0;
                        if (facetype == 0 || facetype == 5)
                            I = 1;
                        unsigned int nodecount = f.size() - 2 - I;
                        // 面至少由3点组成
                        if (nodecount < 3)
                        {
                            j--;
                            continue;
                        }

                        vector<long> pts;
                        for (unsigned int k = I; k < nodecount + I; k++)
                        {
                            // 原始序号从1开始，改为从0开始
                            pts.push_back(f[k] - 1);
                        }
                        _Faces[j] = pts;
                        // 该面的前指和后指网格单元序号
                        neighborCell[j] = f[nodecount + I] - 1;
                        ownerCell[j] = f[nodecount + I + 1] - 1;
                    }
                }
                continue;
            }
            if (isTarget(line, tagCell))
            {
                s = getContentString(line, tagCell);
                s = s.substr(1);
                s = s.substr(0, s.length() - 1);
                auto sa = split(s, " ");
                if (sa.size() >= 5)
                {
                    int id = stoi(sa[0], 0, 16); // - 1;
                    long istart = stol(sa[1], 0, 16) - 1;
                    long iend = stol(sa[2], 0, 16) - 1;
                    int eletype = stoi(sa[4]);
                    int bcType = stoi(sa[3], 0, 16);

                    // 存储cell信息
                    cellSections[id] = make_tuple(istart, iend, bcType, eletype);
                }
                continue;
            }
        }

        sr.close();
        return bSuccess;
    }
    bool ReadMesh_Fluent_Mesh_Binary(string mshFileName)
    {
        ifstream sr(mshFileName, ios::in | ios::binary);
        if (!sr.good())
        {
            m_strErr = mshFileName + " 文件不存在";
            return false;
        }

        // 信息复位
        Init();

        // 格式的标志
        int tagDim = 2,
            tagNode = 10, tagNodeDouble = 3010, tagNodeSingle = 2010,
            tagFace = 13, tagFaceInt = 2013,
            tagCell = 12, tagCellInt = 2012,
            tagZoneTye = 45,
            tagZoneSection = 39;
        bool bDim = false,
             bNode = false,
             bFace = false,
             bCell = false;
        bool bSuccess = true;
        string endSection = "End of Binary Section   ";
        auto seg = ReadLine(sr); // ReadSegment(sr, '(', ')');
        while (seg.size() > 0)
        {
            // ProcessSegment(seg);
            int splitPos = -1, tag = -1;
            if (getTag(seg, ' ', splitPos, tag))
            {
                if (!bDim && tag == tagDim)
                {
                    nDimensions = o2i(string(seg.begin() + splitPos + 1, seg.end()));
                    if (nDimensions != 3)
                    {
                        m_strErr = "不是三维问题,暂不能处理";
                        bSuccess = false;
                        break;
                    }
                    bDim = true;
                    ShowPrompt("维数:" + to_string(nDimensions));
                }
                else if (!bNode && tag == tagNode)
                {
                    auto s = getStringContent(seg, splitPos);
                    auto sa = split(s, " ");
                    if (sa.size() >= 4)
                    {
                        nNodes = stoi(sa[2], 0, 16);
                        _Nodes = make_shared<vector3D<T>[]>(nNodes);
                    }
                    bNode = true;

                    ShowPrompt("节点数:" + to_string(nNodes));
                }
                else if (!bFace && tag == tagFace)
                {
                    auto s = getStringContent(seg, splitPos);
                    auto sa = split(s, " ");
                    if (sa.size() >= 4)
                    {
                        nFaces = stoi(sa[2], 0, 16);

                        _Faces = make_shared<vector<long>[]>(nFaces);
                        neighborCell = make_shared<long[]>(nFaces);
                        ownerCell = make_shared<long[]>(nFaces);
                    }
                    bFace = true;

                    ShowPrompt("面元数:" + to_string(nFaces));
                }
                else if (!bCell && tag == tagCell)
                {
                    auto s = getStringContent(seg, splitPos);
                    auto sa = split(s, " ");
                    if (sa.size() >= 4)
                    {
                        nCells = stoi(sa[2], 0, 16);
                    }
                    bCell = true;

                    ShowPrompt("网格数:" + to_string(nCells));
                }
                // 节点,可能会有多段
                else if (tag == tagNodeDouble) // double 型数据
                {
                    auto s = getStringContent(seg, splitPos);
                    auto sa = split(s, " ");
                    if (sa.size() >= 4)
                    {
                        int istart = stoi(sa[1], 0, 16) - 1;
                        int iend = stoi(sa[2], 0, 16) - 1;

                        ShowPrompt("正在读取 节点:" + s);
                        char read[1];
                        // 应是'('
                        sr.read(read, 1);
                        long npts = (iend - istart + 1);
                        double *temp = new double[npts * 3];
                        sr.read((char *)temp, npts * 3 * sizeof(double));
                        double *p = temp;
                        for (long k = istart; k <= iend; k++)
                        {
                            // long L=k-istart;
                            // auto p = temp+(3 * L);
                            _Nodes[k] = vector3D<T>(p[0] * mm, p[1] * mm, p[2] * mm);
                            p += 3;
                        }

                        // 应是')'
                        seg = ReadLine(sr);
                        // 应是 End of Binary....
                        seg = ReadLine(sr);
                    }
                }
                // 面
                else if (tag == tagFaceInt)
                {
                    auto s = getStringContent(seg, splitPos);
                    auto sa = split(s, " ");
                    if (sa.size() >= 5)
                    {
                        ShowPrompt("正在读取 面元:" + s);
                        int id = stoi(sa[0], 0, 16);
                        long istart = stol(sa[1], 0, 16) - 1;
                        long iend = stol(sa[2], 0, 16) - 1;
                        int bcType = stoi(sa[3], 0, 16);
                        int facetype = stoi(sa[4]);
                        // 读取
                        if (facetype == 0)
                        {
                            ShowPrompt("混合网格id:" + to_string(id));
                        }
                        if (facetype == 5)
                        {
                            ShowPrompt("多边形网格id:" + to_string(id));
                        }
                        if (bcType != 2)
                        {
                            ShowPrompt("边界id:" + to_string(id));
                            FaceBD.push_back(id);
                        }
                        else
                            FaceInner.push_back(id);
                        // 存储面集合信息
                        faceSections[id] = make_tuple(istart, iend, bcType, facetype);

                        char read[1];
                        // 应是'('
                        sr.read(read, 1);
                        seg = ReadToEnd(sr, endSection);
                        int nodecount = facetype;
                        //  I为起点
                        int I = 0;
                        // 数据内存
                        int *L = (int *)seg.data();
                        if (facetype == 0 || facetype == 5)
                        {
                            // 如果是混合或多边形,则前面多一个节点数:nd pt1 pt2 pt3 .... prev next
                            // 尚未测试
                            // throw("not test");
                            I = 1;
                            int J = 0;
                            for (long j = istart; j <= iend; j++)
                            {
                                nodecount = L[J];
                                vector<long> pts;
                                for (int k = I; k < nodecount + I; k++)
                                {
                                    // 原始序号从1开始，改为从0开始
                                    J++;
                                    pts.push_back(L[J] - 1);
                                }
                                _Faces[j] = pts;
                                // 该面的前指和后指网格单元序号
                                J++;
                                neighborCell[j] = L[J] - 1;
                                J++;
                                ownerCell[j] = L[J] - 1;
                                J++;
                            }
                        }
                        else // 一行内容为: pt1 pt2 pt3 .... prev next
                        {
                            for (long j = istart; j <= iend; j++)
                            {
                                long _start = (j - istart) * (nodecount + I + 2);
                                vector<long> pts;
                                for (int k = I; k < nodecount + I; k++)
                                {
                                    // 原始序号从1开始，改为从0开始
                                    pts.push_back(L[_start + k] - 1);
                                }
                                _Faces[j] = pts;
                                // 该面的前指和后指网格单元序号
                                neighborCell[j] = L[_start + nodecount] - 1;
                                ownerCell[j] = L[_start + nodecount + 1] - 1;
                            }
                        }
                    }
                }
                else if (tag == tagCellInt)
                {
                    auto s = getStringContent(seg, splitPos);
                    auto sa = split(s, " ");
                    if (sa.size() >= 5)
                    {
                        ShowPrompt("正在读取 网格单元:" + s);
                        int id = stoi(sa[0], 0, 16); // - 1;
                        long istart = stol(sa[1], 0, 16) - 1;
                        long iend = stol(sa[2], 0, 16) - 1;
                        int eletype = stoi(sa[4]);
                        int bcType = stoi(sa[3], 0, 16);

                        // 存储cell信息
                        cellSections[id] = make_tuple(istart, iend, bcType, eletype);
                        char read[1];
                        // 应是'('
                        sr.read(read, 1);
                        seg = ReadToEnd(sr, endSection);
                    }
                }
                else if (tag == tagZoneSection)
                {
                    auto s = getStringContent(seg, splitPos);
                    auto sa = split(s, " ");
                    if (sa.size() >= 3)
                    {
                        string s2 = ToLower(sa[2]);
                        ShowPrompt(s2);
                    }

                    // 保存信息
                    zoneSections[stoi(sa[0])] = sa[1] + " " + sa[2];
                }
                else if (tag > 2000)
                {
                    GoToEnd(sr, endSection);
                }
            }

            seg.clear();

            seg = ReadLine(sr);
        }
        seg.clear();
        sr.close();

        return true;
    }
    // 保存为msh文本格式文件
    void SaveMesh_Fluent_Mesh_Text(string mshFileName)
    {
        ofstream sw(mshFileName, ios::out);
        // id
        int ID = 0;
        // 三维问题
        sw << "(2 3)" << endl;
        // 总网格数
        sw << "(10 (0 1 " << hex << nNodes << " 0 3))" << endl;
        // 开始网格内容
        ID++;
        // 节点坐标
        sw << "(10 (" << hex << ID << " 1 " << hex << nNodes << " 1 3)" << endl;
        sw << "(" << endl;
        for (long i = 0; i < nNodes; i++)
            sw << _Nodes[i].x / mm << " " << _Nodes[i].y / mm << " " << _Nodes[i].z / mm << endl;
        sw << "))" << endl;
        // cell
        ID++;
        sw << "(12 (0 1 " << hex << nCells << " 0 0 ))" << endl;
        /* sw << "(12 (" << hex << ID << " 1 " << hex << nCells << " 1 0)(" << endl;
         int n = nCells / 10;
         for (long i = 0; i < n; i++)
         {
             for (int j = 0; j < 10; j++)
                 sw << cellElementType[i * 10 + j] << "  ";
             sw << endl;
         }
         for (long i = n * 10; i < nCells; i++)
             sw << cellElementType[i] << "  ";
         sw << endl;
         sw << "))" << endl;
 */
        // faces
        vector<int> usedID;
        sw << "(13 (0 1 " << hex << nFaces << " 0 0))" << endl;
        map<int, tuple<long, long, int, int>>::iterator it = faceSections.begin();
        map<int, tuple<long, long, int, int>>::iterator itEnd = faceSections.end();
        while (it != itEnd)
        {
            int id = it->first;
            usedID.push_back(id);

            auto fs = it->second;
            int istart = get<0>(fs), iend = get<1>(fs), bcType = get<2>(fs), faceType = get<3>(fs);
            sw << "(13 (" << hex << id << "  " << hex << istart + 1 << " " << hex << iend + 1 << " " << bcType << " " << faceType << ")(" << endl;
            if (faceType != 0)
                for (long i = istart; i <= iend; i++)
                {
                    long n = _Faces[i].size();
                    for (long j = 0; j < n; j++)
                        sw << hex << _Faces[i][j] + 1 << " ";
                    sw << hex << neighborCell[i] + 1 << " " << hex << ownerCell[i] + 1 << endl;
                }
            else
                for (long i = istart; i <= iend; i++)
                {
                    long n = _Faces[i].size();
                    sw << n << "  ";
                    for (long j = 0; j < n; j++)
                        sw << hex << _Faces[i][j] + 1 << " ";
                    sw << hex << neighborCell[i] + 1 << " " << hex << ownerCell[i] + 1 << endl;
                }
            sw << "))" << endl;
            it++;
        }

        for (size_t i = 0; i < usedID.size(); i++)
        {
            auto id = usedID[i];
            auto content = zoneSections[id];
            sw << "(39 (" << dec << id << " " << content << ")())" << endl;
        }

        sw.close();
    }

private:
    void ProcessSegment(vector<char> &seg)
    {
        int splitPos = -1, tag = -1;
        if (getTag(seg, ' ', splitPos, tag))
        {
        }

        seg.clear();
    }
    // 对形如（### ******）的数据块取###,认为是整数
    bool getTag(vector<char> &seg, char split, int &splitPos, int &tag, int _start = 1)
    {
        splitPos = -1;
        tag = -1;
        for (size_t j = _start; j < seg.size(); j++)
            if (seg[j] == split && j > 0)
            {
                splitPos = j;
                try
                {
                    tag = stoi(string(seg.begin() + _start, seg.begin() + splitPos));
                }
                catch (...)
                {
                    return false;
                }
                return true;
            }
        return false;
    }
    string getStringContent(vector<char> &seg, int &start)
    {
        int _s = -1, _e = -1;
        for (size_t i = start; i < seg.size(); i++)
        {
            if (_s < 0 && seg[i] == '(')
                _s = i;
            else if (_e < 0 && seg[i] == ')')
            {
                _e = i;
                break;
            }
        }
        if (_s < 0)
            _s = 0;
        if (_e < 0)
            _e = seg.size() - 1;
        start = _e;
        auto s = string(seg.begin() + _s + 1, seg.begin() + _e);
        return s;
    }
    bool getSegmentData(vector<char> &seg, int &start, int &end)
    {
        end = -1;
        for (int i = seg.size() - 1; i > 0; i--)
            if (seg[i] == ')')
            {
                end = i;
                break;
            }
        if (end < 0)
            return false;

        int _start = -1;
        for (int i = start; i < seg.size(); i++)
            if (seg[i] == '(')
            {
                _start = i;
                break;
            }
        if (_start < 0)
            return false;

        start = _start;
        return true;
    }

    vector<char> strToByte(string s)
    {
        vector<char> ret;
        char *bytes = new char[s.length()];
        memcpy(bytes, s.data(), s.length());
        ret.insert(ret.end(), bytes, bytes + s.length());
        delete[] bytes;
        return ret;
    }
    long FindPos(char *data, char find, size_t iStart = 0)
    {
        long j = 0;
        int len = strlen(data) + iStart;
        for (j = iStart; j < len; j++)
            if (data[j] == find)
                return j;
        return -1;
    }
    long FindPos(char *data, string find, size_t iStart, size_t len)
    {
        size_t j = 0;
        size_t lenFind = find.length();
        for (j = iStart; j < len; j++)
        {
            // 找到起点
            if (data[j] == find[0])
            {
                bool b = true;
                if (j + lenFind < len)
                {
                    for (size_t i = 1; i < lenFind; i++)
                    {
                        if (find[i] != data[j + i])
                        {
                            j = j + i;
                            b = false;
                            break;
                        }
                    }
                    if (b)
                        return j;
                }
                else
                {
                    for (size_t i = 1; i < len - j; i++)
                    {
                        if (find[i] != data[j + i])
                        {
                            j = j + i;
                            b = false;
                            break;
                        }
                    }
                    if (b)
                        return j;
                }
            }
        }
        return -1;
    }
    long FindPos(vector<char> data, string find, size_t iStart, size_t len)
    {
        size_t j = 0;
        size_t lenFind = find.length();
        for (j = iStart; j < len; j++)
        {
            // 找到起点
            if (data[j] == find[0])
            {
                bool b = true;
                if (j + lenFind < len)
                {
                    for (size_t i = 1; i < lenFind; i++)
                    {
                        if (find[i] != data[j + i])
                        {
                            j = j + i;
                            b = false;
                            break;
                        }
                    }
                    if (b)
                        return j;
                }
                else
                {
                    for (size_t i = 1; i < len - j; i++)
                    {
                        if (find[i] != data[j + i])
                        {
                            j = j + i;
                            b = false;
                            break;
                        }
                    }
                    if (b)
                        return j;
                }
            }
        }
        return -1;
    }
    string FindTarget(vector<char> data, map<string, vector<char>> tagDataTitle) //, vector<char> isDataTitle)
    {
        bool isData = true;
        for (auto it = tagDataTitle.begin(); it != tagDataTitle.end(); it++)
        {
            string key = it->first;
            vector<char> title = it->second;
            if (data.size() < title.size())
                continue;

            int n = title.size();
            isData = true;
            for (int i = 0; i < n; i++)
            {
                if (title[i] != data[i])
                {
                    isData = false;
                    break;
                }
            }
            if (isData)
                return key;
        }
        return "";
    }
    string trim(string str)
    {
        const char *typeOfWhitespaces = " tnrfv";
        str.erase(str.find_last_not_of(typeOfWhitespaces) + 1);
        str.erase(0, str.find_first_not_of(typeOfWhitespaces));
        return str;
    }
    vector<string> stringSplit(const string &str, char delim)
    {
        stringstream ss(str);
        string item;
        vector<string> elems;
        while (getline(ss, item, delim))
        {
            if (!item.empty())
            {
                elems.push_back(trim(item));
            }
        }
        return elems;
    }
    // 读取到标志,内容不保存
    void GoToEnd(ifstream &br, string end)
    {

        while (br)
        {
            char read[NReadBytes];
            br.read(read, NReadBytes);
            size_t n = br.gcount();
            long j = FindPos(read, end, 0, n);
            if (j >= 0 && j + end.length() <= n)
            {
                br.seekg(j + end.length() - n - 1, ios_base::cur);
                return;
            }
            if (j > 0)
                br.seekg(j - n - 1, ios_base::cur);
            string tag0 = "(0 \"";
            j = FindPos(read, end, 0, n);
            if (j >= 0 && j + end.length() <= n)
            {
                br.seekg(j + end.length() - n - 1, ios_base::cur);
                return;
            }
            if (j > 0)
                br.seekg(j - n - 1, ios_base::cur);
        }
    }
    void GoToEnd(ifstream &br, string end, vector<char> read)
    {
        size_t n = read.size();
        long j = FindPos(read, end, 0, n);
        if (j >= 0 && j + end.length() <= n)
        {
            br.seekg(j + end.length() - n - 1, ios_base::cur);
            return;
        }
        if (j > 0)
            br.seekg(j - n - 1, ios_base::cur);
        string tag0 = "(0 \"";
        j = FindPos(read, end, 0, n);
        if (j >= 0 && j + end.length() <= n)
        {
            br.seekg(j + end.length() - n - 1, ios_base::cur);
            return;
        }
        if (j > 0)
            br.seekg(j - n - 1, ios_base::cur);
    }
    // 一次读取从()内的内容
    vector<char> ReadSegment(ifstream &br, char start, char end, bool bnextchar = true)
    {
        vector<char> bits;
        bits.reserve(NReadBytes * 3);
        bool good = false; // 是否读到结尾了？
        int b = 0;
        while (br)
        {
            char read[NReadBytes];
            br.read(read, NReadBytes);
            size_t n = br.gcount();
            for (size_t j = 0; j < n; j++)
            {
                if (read[j] == end && b > 0)
                {
                    char next;
                    if (bnextchar)
                    {
                        if (j == n - 1)
                        {
                            br.read(&next, 1);
                            if (br.gcount() <= 0)
                                break;
                        }
                        else
                            next = read[j + 1];
                    }
                    else
                        next = end;

                    if (next == ')' || next == 0x20 || next == 0xa || next == '\t')
                    {
                        b--;
                        if (b == 0)
                        {
                            if (bnextchar)
                                br.seekg(j - n + 2, ios_base::cur);
                            else
                                br.seekg(j - n + 1, ios_base::cur);

                            good = true;
                            break;
                        }
                    }
                }
                if (b > 0)
                {
                    bits.push_back(read[j]);
                }
                if (read[j] == start && b < 2)
                {
                    char prev = 0x20;
                    if (j > 0)
                        prev = read[j - 1];
                    if (prev == 0x20 || prev == 0xa || prev == '\t')
                        b++;
                }
            }
            if (good)
                break;
        }

        return bits;
    }
    vector<char> ReadLine(ifstream &br)
    {
        vector<char> bits;
        while (br)
        {
            char read[NReadBytes];
            br.read(read, NReadBytes);
            size_t n = br.gcount();
            long j = FindPos(read, 0x0a);
            if (j < 0)
                j = FindPos(read, ')');
            if (bits.size() > 0 && j >= 0)
            {
                br.seekg(j - n + 1, ios_base::cur);
                bits.insert(bits.end(), read, read + j);
                break;
            }
            else if (bits.size() == 0 && j == 0)
            {
                j = FindPos(read, 0x0a, 1);
                if (j < 0)
                    j = FindPos(read, ')', 1);
                if (j > 0)
                {
                    br.seekg(j - n + 1, ios_base::cur);
                    bits.insert(bits.end(), read + 1, read + j);
                    break;
                }
            }
            else if (j > 0)
            {
                br.seekg(j - n + 1, ios_base::cur);
                bits.insert(bits.end(), read, read + j);
                break;
            }

            bits.insert(bits.end(), read, read + n);
        }

        return bits;
    }
    vector<char> ReadToEnd(ifstream &br, string end)
    {
        vector<char> bits;
        bits.reserve(3 * NReadBytes);

        long endLen = end.length();
        while (br)
        {
            char read[NReadBytes];
            br.read(read, NReadBytes);
            long n = (long)br.gcount();
            long j = FindPos(read, end, 0, n);
            if (j >= 0 && j + endLen <= n)
            {
                br.seekg(j + endLen - n - 1, ios_base::cur);
                bits.insert(bits.end(), read, read + j);
                return bits;
            }
            if (j > 0)
            {
                br.seekg(j - n - 1, ios_base::cur);
                bits.insert(bits.end(), read, read + j);
            }
            else
                bits.insert(bits.end(), read, read + n);
        }
        return bits;
    }
    // 清除数据
    void ClearData()
    {
        size_t n = cellData.size();
        if (n > 0)
        {
            for (auto it = cellData.begin(); it != cellData.end(); it++)
            {
                delete[] (it->second);
            }
        }
        cellData.clear();

        n = faceData.size();
        if (n > 0)
        {
            for (auto it = faceData.begin(); it != faceData.end(); it++)
            {
                delete[] (it->second);
            }
        }
        faceData.clear();
    }

    // 读取数据
    bool ReadData_Fluent_Data_Binary(string datFileName, map<string, string> datMap, bool bClearData = true)
    {
        if (bClearData)
            ClearData();

        size_t n = datMap.size();
        if (n <= 0)
            return false;
        // 分配内存
        // auto isDataTitle = strToByte("(0 \"SV_");
        string tagDataStartDouble = "(3300 (",
               tagDataStart3314 = "(3314 (",
               tagSV = "(0 \"SV_",
               tag0 = "(0 \"";
        // string end = "End of Binary Section   3300))";
        string end = "End of Binary Section   ";
        map<string, vector<char>> tagDataTitle; // 二进制标记
        for (auto it = datMap.begin(); it != datMap.end(); it++)
        {
            // ShowPrompt(it->first);
            cellData[it->second] = new T[nCells];
            faceData[it->second] = new T[nFaces];
            tagDataTitle[it->first] = strToByte("(0 \"" + it->first + ",");
        }
        // 开始读取
        bool bSuccess = true;
        T *curData = NULL;
        ifstream sr(datFileName, ios::in | ios::binary);
        if (!sr.good())
        {
            m_strErr = datFileName + " 文件不存在";
            return false;
        }

        GoToEnd(sr, "Data:\")");

        while (sr)
        {
            vector<char> bytes = ReadLine(sr);
            if (bytes.size() == 0)
                continue;
            string key = FindTarget(bytes, tagDataTitle);
            if (key != "")
            {
                key = datMap[key];
                string s(bytes.begin(), bytes.end());
                // ShowPrompt(" readed " + s);
                vector<string> sa = stringSplit(s, ',');
                s = sa.back();
                if (s.find("cells:") != s.npos)
                {
                    curData = cellData[key];
                    // 是cell的值,后面应跟着 (3300
                    bytes.clear();
                    bytes = ReadLine(sr);
                    if (bytes.size() == 0)
                        continue;
                    //
                    string s1(bytes.begin(), bytes.end());
                    if (isTarget(s1, tagDataStartDouble))
                    {
                        s = getContentString(s1, tagDataStartDouble);
                        s = s.substr(0, s.length() - 1);
                        sa = stringSplit(s, ' ');
                        if (sa.size() >= 6)
                        {
                            long istart = stol(sa[5]) - 1,
                                 iend = stol(sa[6]) - 1;
                            // 读取到curData
                            n = (iend - istart + 1) * sizeof(double);
                            char read[1];
                            // 应是'('
                            sr.read(read, 1);
                            sr.read((char *)curData, n);

                            ShowPrompt(key + " readed " + s);
                        }
                    }
                    else if (isTarget(s1, tagDataStart3314))
                    {
                        s = getContentString(s1, tagDataStart3314);
                        s = s.substr(0, s.length() - 1);
                        sa = stringSplit(s, ' ');
                        if (sa.size() >= 6)
                        {
                            long istart = stol(sa[5]) - 1,
                                 iend = stol(sa[6]) - 1;
                            // 读取到curData
                            n = (iend - istart + 1) * sizeof(double);
                            char read[1];
                            // 应是'('
                            sr.read(read, 1);
                            sr.read((char *)curData, n);

                            ShowPrompt(key + " readed " + s);
                        }
                    }
                }
                else if (s.find("faces:") != s.npos) // 面上的值
                {
                    curData = faceData[key];
                    bytes.clear();
                    bytes = ReadLine(sr);
                    if (bytes.size() == 0)
                        continue;
                    //
                    string s1(bytes.begin(), bytes.end());
                    if (isTarget(s1, tagDataStartDouble))
                    {
                        s = getContentString(s1, tagDataStartDouble);
                        s = s.substr(0, s.length() - 1);
                        sa = stringSplit(s, ' ');
                        if (sa.size() >= 6)
                        {
                            long istart = stol(sa[5]) - 1,
                                 iend = stol(sa[6]) - 1;
                            // 应是'('+double
                            n = (iend - istart + 1) * sizeof(double);
                            size_t m = istart * sizeof(double);
                            char read[1];
                            // 应是'('
                            sr.read(read, 1);
                            sr.read((char *)curData + m, n);

                            ShowPrompt(key + " readed " + s);
                            // sr.seekg(n, ios_base::cur);
                        }
                    }
                    else if (isTarget(s1, tagDataStart3314))
                    {
                        s = getContentString(s1, tagDataStart3314);
                        s = s.substr(0, s.length() - 1);
                        sa = stringSplit(s, ' ');
                        if (sa.size() >= 6)
                        {
                            long istart = stol(sa[5]) - 1,
                                 iend = stol(sa[6]) - 1;
                            // 应是'('+double
                            n = (iend - istart + 1) * sizeof(double);
                            size_t m = istart * sizeof(double);
                            char read[1];
                            // 应是'('
                            sr.read(read, 1);
                            sr.read((char *)curData + m, n);

                            ShowPrompt(key + " readed " + s);
                            // sr.seekg(n, ios_base::cur);
                        }
                    }
                }
                GoToEnd(sr, end);
            }
            else // 不是要的数据
            {
                string s(bytes.begin(), bytes.end());
                if (isTarget(s, tagSV) || isTarget(s, tag0))
                {
                    int bug = 0;
                    if (s == "(0 \"Residuals:\")")
                    {
                        bug = 1;
                    }
                    ShowPrompt(" readed " + s);
                    bytes.clear();
                    bytes = ReadLine(sr);
                    if (bytes.size() == 0)
                        continue;
                    //
                    string s1(bytes.begin(), bytes.end());
                    if (isTarget(s1, tagDataStartDouble))
                    {
                        s = getContentString(s1, tagDataStartDouble);
                        s = s.substr(0, s.length() - 1);
                        vector<string> sa = stringSplit(s, ' ');
                        if (sa.size() >= 6)
                        {
                            long istart = stol(sa[5]) - 1,
                                 iend = stol(sa[6]) - 1;
                            // 矢量个数?
                            long vn = stol(sa[2]);
                            // 应是'('+double
                            n = (iend - istart + 1) * sizeof(double) * vn + 1;
                            sr.seekg(n, ios_base::cur);
                        }
                        GoToEnd(sr, end);
                    }
                    else if (isTarget(s1, tagDataStart3314))
                    {
                        s = getContentString(s1, tagDataStart3314);
                        s = s.substr(0, s.length() - 1);
                        vector<string> sa = stringSplit(s, ' ');
                        if (sa.size() >= 6)
                        {
                            long istart = stol(sa[5]) - 1,
                                 iend = stol(sa[6]) - 1;
                            // 矢量个数?
                            long vn = stol(sa[2]);
                            // 应是'('+double
                            n = (iend - istart + 1) * sizeof(double) * vn + 1;
                            sr.seekg(n, ios_base::cur);
                        }
                        GoToEnd(sr, end);
                    }
                    else
                    {
                        key = FindTarget(bytes, tagDataTitle);
                        if (key != "")
                        {
                            // 回退
                            sr.seekg(-bytes.size() - 1, ios_base::cur);
                            // continue;
                        }
                        else
                            GoToEnd(sr, end);
                    }
                }
                else
                    GoToEnd(sr, end, bytes);
            }

            bytes.clear();
        }

        sr.close();
        return true;
    }
    // 保存为data二进制格式文件
    void SaveData_Fluent_Data_Binary(string datFileName)
    {
        string tagDataStartDouble = "(3300 (",
               tagSV = "(0 \"SV_";
        string end = "End of Binary Section   3300))\n\n";

        ofstream sw(datFileName, ios::out | ios::binary);
        string sdata = "(0 \"Data:\")\n";
        sw.write(sdata.data(), sdata.length());

        size_t n = cellData.size();
        if (n > 0)
        {
            for (auto it = cellData.begin(); it != cellData.end(); it++)
            {
                // 写title
                string key(it->first);
                string s = tagSV + key + ",domain X, zoneX, X cells:\")\n" + tagDataStartDouble + "15 15 1 0 1 1 " + to_string(nCells) + ")\n(";
                sw.write(s.data(), s.length());
                // 写数据段
                sw.write(reinterpret_cast<const char *>(it->second), sizeof(double) * nCells);
                sw.write(end.data(), end.length());
            }
        }
        n = faceData.size();
        if (n > 0)
        {
            for (auto it = faceData.begin(); it != faceData.end(); it++)
            {
                // 写title
                string key(it->first);
                string s = tagSV + key + ",domain X, zoneX, X faces:\")\n" + tagDataStartDouble + "15 15 1 0 1 1 " + to_string(nFaces) + ")\n(";
                sw.write(s.data(), s.length());
                // 写数据段
                sw.write(reinterpret_cast<const char *>(it->second), sizeof(double) * nFaces);
                sw.write(end.data(), end.length());
            }
        }

        sw.close();
    }

public:
    void ClearUnusedPoints()
    {
        auto all = new VertexInfo<T>[nNodes];
        for (long i = 0; i < nFaces; i++)
        {
            auto f = _Faces[i];
            // 该面点数
            auto npts = f.size();
            if (ownerCell[i] >= 0)
            {
                for (long k = 0; k < npts; k++)
                {
                    all[f[k]].distPart[ownerCell[i]] = 0.0;
                }
            }
            if (neighborCell[i] >= 0)
            {
                for (long k = 0; k < npts; k++)
                {

                    all[f[k]].distPart[neighborCell[i]] = 0.0;
                }
            }
        }
        // 搜索没有用到的点
        vector<long> unusedPts;

        for (long i = 0; i < nNodes; i++)
        {
            if (all[i].distPart.size() <= 0)
            {
                unusedPts.push_back(i);
            }
        }
        delete[] all;

        int n = unusedPts.size();
        if (n <= 0)
            return;
        // 开始清理
        long start = -1, end = -1;
        // start -- --- subtraction
        vector<tuple<long, long>> minus;
        for (int i = n - 1; i >= 0; i--)
        {
            // 清理
            if (start == -1)
            {
                start = unusedPts[i];
                end = start;
                continue;
            }
            if (end - unusedPts[i] == 1)
            {
                end -= 1;
            }
            else
            {
                // start---end 为连续段
                minus.push_back(make_tuple(start + 1, start - end + 1));
                start = -1;
            }
        }
        // 剩余
        minus.push_back(make_tuple(start + 1, start - end + 1));

        // 调整面中的序号
        n = minus.size();
        for (int i = n - 2; i >= 0; i--)
        {
            auto item = minus[i];
            get<1>(item) = get<1>(item) + get<1>(minus[i + 1]);
        }
        // 生成map
        //  old-->new
        map<long, long> change, change2;

        start = nNodes - 1;
        for (int i = n - 1; i >= 0; i--)
        {
            end = get<0>(minus[i]);
            auto m = get<1>(minus[i]);
            for (long j = end; j <= start; j++)
            {
                change[j] = j - m;
                change2[j - m] = j;
            }
        }
        for (long i = 0; i < nFaces; i++)
        {
            // 该面点数
            auto npts = _Faces[i].size();
            for (long k = 0; k < npts; k++)
            {
                _Faces[i][k] = change[_Faces[i][k]];
            }
        }

        nNodes = nNodes - unusedPts.size();
        // 前移
        for (long i = 0; i < nNodes; i++)
        {
            _Nodes[i] = _Nodes[change2[i]];
        }
    }
    // 消除vect中重复单元并排序
    void removeDups(vector<long> &vect)
    {
        sort(vect.begin(), vect.end());
        vect.erase(unique(vect.begin(), vect.end()), vect.end());
        // return;
        // int u = 0;

        // for (int i = 1; i < vect.size();)
        // {
        //     if (vect[u] == vect[i])
        //     {
        //         ++i;
        //     }
        //     else
        //     {
        //         ++u;
        //         swap(vect[u], vect[i]);
        //         ++i;
        //     }
        // }
        // vect.erase(vect.begin() + u + 1, vect.end());
    }
    void MeshPartion(int _nParts)
    {
#ifdef MPI

        // 对网格进行分区,以cell为基础
        // 先拼接cell
        MergeCellsForPartion();

        // eptr: 面指针位置；nFaces+1
        // eind: 面元索引；
        vector<idx_t> eptr, eind;
        eptr.reserve(nCells + 1);
        eind.reserve(nCells * 4); // 假定为4面体

        idx_t ptr = 0;
        eptr.push_back(ptr);
        for (long i = 0; i < nCells; i++)
        {
            auto npts = _Cells[i].size();
            ptr += npts;
            eptr.push_back(ptr);
            eind.insert(end(eind), begin(_Cells[i]), end(_Cells[i]));
        }

        idx_t objval, ne = nCells, nn = nFaces, nParts = _nParts;
        vector<idx_t> epart(nCells, 0), npart(nFaces, 0);
        int ret = METIS_PartMeshNodal(&ne, &nn, eptr.data(), eind.data(), NULL, NULL, &nParts, NULL, NULL, &objval, epart.data(), npart.data());

        map<int, vector<long>> cellMap;
        for (int i = 0; i < _nParts; i++)
        {
            cellMap[i].reserve(nCells / _nParts + 1);
        }
        for (long i = 0; i < nCells; i++)
        {
            cellMap[epart[i]].push_back(i);
        }

        // 形成mesh
        shared_ptr<Mesh3<T>[]> meshs = make_shared<Mesh3<T>[]>(_nParts);
        // 分区时，交界面出的单元序号映射:global id--(partid-local id)
        map<long, tuple<int, long>> mapInterCells;
        // 分区时，每个分区中需要接收的cell id ： local face id -- global cell id
        shared_ptr<map<long, long>[]> RecivCells = make_shared<map<long, long>[]>(_nParts);
        for (int i = 0; i < _nParts; i++)
        {
            // 统计用到的node和face id
            vector<long> nodeSet, faceSet;
            map<long, long> nodeOld2New, faceOld2New, cellOld2New;
            long ncells = cellMap[i].size();
            long j = 0;
            for (j = 0; j < ncells; j++)
            {
                long oldCellID = cellMap[i][j];
                cellOld2New[oldCellID] = j;

                faceSet.insert(end(faceSet), begin(_Cells[oldCellID]), end(_Cells[oldCellID]));
            }
            removeDups(faceSet);

            long nfaces = faceSet.size();
            Mesh3<T> *subMesh = &(meshs[i]);
            subMesh->setFacesNum(nfaces);
            subMesh->nCells = ncells;
            // 面序号排序
            for (j = 0; j < nfaces; j++)
            {
                long oldFaceID = faceSet[j];
                faceOld2New[oldFaceID] = j;

                nodeSet.insert(end(nodeSet), begin(_Faces[oldFaceID]), end(_Faces[oldFaceID]));
            }
            removeDups(nodeSet);

            // 点序号排序
            long nnodes = nodeSet.size();
            subMesh->nNodes = nnodes;
            subMesh->_Nodes = make_shared<vector3D<T>[]>(nnodes);
            for (j = 0; j < nnodes; j++)
            {
                long oldNodeID = nodeSet[j];
                nodeOld2New[oldNodeID] = j;
                subMesh->_Nodes[j].x = _Nodes[oldNodeID].x;
                subMesh->_Nodes[j].y = _Nodes[oldNodeID].y;
                subMesh->_Nodes[j].z = _Nodes[oldNodeID].z;
            }

            // 更新面序号
            for (j = 0; j < nfaces; j++)
            {
                long oldFaceID = faceSet[j];
                int np = _Faces[oldFaceID].size();
                for (long k = 0; k < np; k++)
                    subMesh->_Faces[j].push_back(nodeOld2New[_Faces[oldFaceID][k]]);
                // 面的cell
                int interType = 0;
                if (ownerCell[oldFaceID] >= 0)
                {
                    if (cellOld2New.count(ownerCell[oldFaceID]) > 0)
                        subMesh->ownerCell[j] = cellOld2New[ownerCell[oldFaceID]];
                    else
                    {
                        subMesh->ownerCell[j] = CELLID_SHIFT(ownerCell[oldFaceID]);
                        // 此面为交界面
                        subMesh->interFaces.push_back(j);
                        RecivCells[i][j] = ownerCell[oldFaceID];
                        interType = 1;
                    }
                }
                else
                    subMesh->ownerCell[j] = -1;

                if (neighborCell[oldFaceID] >= 0)
                {
                    if (cellOld2New.count(neighborCell[oldFaceID]) > 0)
                        subMesh->neighborCell[j] = cellOld2New[neighborCell[oldFaceID]];
                    else
                    {
                        subMesh->neighborCell[j] = CELLID_SHIFT(neighborCell[oldFaceID]);
                        // 此面为交界面
                        subMesh->interFaces.push_back(j);
                        RecivCells[i][j] = neighborCell[oldFaceID];
                        interType = 2;
                    }
                }
                else
                    subMesh->neighborCell[j] = -1;

                // 是否为交界面
                switch (interType)
                {
                case 1:
                    mapInterCells[neighborCell[oldFaceID]] = make_tuple(i, cellOld2New[neighborCell[oldFaceID]]);
                    break;
                case 2:
                    mapInterCells[ownerCell[oldFaceID]] = make_tuple(i, cellOld2New[ownerCell[oldFaceID]]);
                    break;
                default:
                    break;
                }
            }
        }

        for (int i = 0; i < _nParts; i++)
        {
            Mesh3<T> *subMesh = &(meshs[i]);
            // 更新其数据交换
            for (auto it = RecivCells[i].begin(); it != RecivCells[i].end(); it++)
            {
                long faceid = it->first, cellid = it->second;
                // 找到cellid所在的单元
                auto partionid = get<0>(mapInterCells[cellid]);
                auto localcellid = get<1>(mapInterCells[cellid]);

                subMesh->RecivCells[partionid][faceid] = localcellid;

                // 对方要发送
                meshs[partionid].SendCells[i].push_back(localcellid);
            }

            subMesh->ConstructCells();

            // 输出测试
            // string fn = "D:\\" + to_string(i) + ".plt";
            // T *data = new T[subMesh->nCells];
            // for (long j = 0; j < subMesh->nCells; j++)
            //     data[j] = (T)i;
            // subMesh->template outputTecPlotFEM<T>(fn, "pp", data);
            // delete data;
            //
        }

#endif
    }

    // 从面拼接cells
    void MergeCellsForPartion()
    {
        _Cells = make_shared<vector<long>[]>(nCells);
        for (long i = 0; i < nFaces; i++)
        {
            if (ownerCell[i] >= 0)
                _Cells[ownerCell[i]].push_back(i);

            if (neighborCell[i] >= 0)
                _Cells[neighborCell[i]].push_back(i);
        }
    }
    // 计算网格参数
    bool ConstructCells()
    {
        ClearUnusedPoints();

        gc = make_shared<T[]>(nFaces);
        dcf = make_shared<T[]>(nFaces);

        faceArea = make_shared<T[]>(nFaces);
        surfaceVector = make_shared<vector3D<T>[]>(nFaces);
        faceCenter = make_shared<vector3D<T>[]>(nFaces);
        ECF = make_shared<vector3D<T>[]>(nFaces);

        _Cells = make_shared<vector<long>[]>(nCells);
        cellCenter = make_shared<vector3D<T>[]>(nCells);
        cellVolume = make_shared<T[]>(nCells);

        cellFaceSign = make_shared<map<long, short>[]>(nCells);

        vertex2Cell = make_shared<VertexInfo<T>[]>(nNodes);

        for (long i = 0; i < nFaces; i++)
        {
            auto f = _Faces[i];
            // 该面点数
            auto npts = f.size();
            // 计算几何中心
            vector3D<T> centre = _Nodes[f[0]];
            // 该面的结构
            for (long k = 1; k < npts; k++)
                centre += _Nodes[f[k]];
            // 几何中心
            centre /= (T)npts;

            // 计算面积向量和面积中心
            vector3D<T> centroid;
            vector3D<T> SF;
            T area = 0.0;

            for (long k = 0; k < npts; k++)
            {
                if (f[k] == 0)
                {
                    int bug = 1;
                }
                auto p2 = _Nodes[f[k]];
                auto p3 = _Nodes[f[(k + 1) % npts]];

                auto lcenter = (centre + p2 + p3) / 3;
                auto lSf = 0.5 * cross(p2 - centre, p3 - centre);
                auto larea = lSf.norm();

                centroid += larea * lcenter;
                SF += lSf;
                area += larea;

                // 统计点集
                // if (ownerCell[i] >= 0)
                //     _CellPts[ownerCell[i]].insert(f[k]);
                // if (neighborCell[i] >= 0)
                //     _CellPts[neighborCell[i]].insert(f[k]);
            }
            T ownerFaceSign = 1;
            if (ownerCell[i] >= 0)
            {
                _Cells[ownerCell[i]].push_back(i);
                cellFaceSign[ownerCell[i]][i] = 1.0;
                ownerFaceSign = 1;
            }
            else
            {
                // 交换
                ownerCell[i] = neighborCell[i];
                neighborCell[i] = -1;
                ownerFaceSign = -1;

                _Cells[ownerCell[i]].push_back(i);
                // 因为 ownerFacesign修正了，导致surfaceVector反映了正确的法向，
                // 所以在该cell中，该面的符号已经正确
                cellFaceSign[ownerCell[i]][i] = 1.0;
            }

            if (neighborCell[i] >= 0)
            {
                _Cells[neighborCell[i]].push_back(i);
                cellFaceSign[neighborCell[i]][i] = -1.0;
            }

            faceCenter[i] = centroid / area;       //* mm;
            surfaceVector[i] = SF * ownerFaceSign; // * mm2;
            faceArea[i] = area;                    //* mm2;
        }

        // 计算单元体积
        for (long i = 0; i < nCells; i++)
        {
            auto c = _Cells[i];
            // 该cell的面元数
            auto nfcs = c.size();
            // 计算几何中心
            auto centre = faceCenter[c[0]];
            for (long k = 1; k < nfcs; k++)
                centre += faceCenter[c[k]];
            centre /= (T)nfcs;

            // 计算体积和体积中心
            vector3D<T> centroid;
            T vol = 0.0;
            for (long k = 0; k < nfcs; k++)
            {
                // auto f2 = _Faces[c[k]];
                auto sign = cellFaceSign[i][c[k]];
                auto SF = surfaceVector[c[k]] * sign;
                auto Cf = faceCenter[c[k]] - centre;
                // 分体积
                auto lvol = SF * Cf / 3.0;
                auto lcenter = 0.75 * faceCenter[c[k]] + 0.25 * centre;
                centroid += lvol * lcenter;
                vol += lvol;
            }
            // 由于几何前面已经转换，此处不需要转换
            cellVolume[i] = vol;
            cellCenter[i] = centroid / vol;
        }
        // // 计算内部面的gc
        for (long i = 0; i < nFaces; i++)
        {
            if (ownerCell[i] >= 0 && neighborCell[i] >= 0)
            {
                auto dFf = cellCenter[neighborCell[i]] - faceCenter[i];
                ECF[i] = cellCenter[neighborCell[i]] - cellCenter[ownerCell[i]];

                dcf[i] = ECF[i].norm();
                gc[i] = dFf.norm() / dcf[i];

                ECF[i] /= dcf[i];
                // 内部面
                inFaces.push_back(i);

                // 记录点与cell间的关系
                auto f = _Faces[i];
                // 该面点数
                auto npts = f.size();
                for (long k = 0; k < npts; k++)
                {
                    vertex2Cell[f[k]].distPart[ownerCell[i]] = 0.0;
                    vertex2Cell[f[k]].distPart[neighborCell[i]] = 0.0;
                }
            }
            else
            {
                ECF[i] = faceCenter[i] - cellCenter[ownerCell[i]];
                dcf[i] = ECF[i].norm();
                gc[i] = 1;

                ECF[i] /= dcf[i];

                bdFaces.push_back(i);

                // 记录点与cell间的关系
                auto f = _Faces[i];
                // 该面点数
                auto npts = f.size();
                for (long k = 0; k < npts; k++)
                {
                    vertex2Cell[f[k]].distPart[ownerCell[i]] = 0.0;
                }
            }
        }
        //从bdFace中清除交界面的信息，只保留真实边界信息。
        if (interFaces.size() > 0)
            bdFaces.erase(remove_if(bdFaces.begin(), bdFaces.end(),
                                    [this](const long &id)
                                    {
                                        auto it = find(interFaces.begin(), interFaces.end(), id);
                                        return (it != interFaces.end());
                                    }),
                          bdFaces.end());
        // 计算点与cellCenter间的距离
        for (long i = 0; i < nNodes; i++)
        {
            for (auto part : vertex2Cell[i].distPart)
                vertex2Cell[i].distPart[part.first] = (_Nodes[i] - cellCenter[part.first]).norm();
        }

        return true;
    }

    // 输出tecplot data 数据；isOnVertex 该数据是定义在点上吗？
    template <class TOther>
    void outputTecPlotFEM(string fileName, string title, vector<shared_ptr<TOther[]>> &data)
    {
        ofstream sw(fileName, ios::out);
        sw << "VARIABLES = \"X\", \"Y\", \"Z\"," << title << endl;
        string ss("4");
        for (size_t i = 1; i < data.size(); i++)
            ss += "," + to_string(i + 4);
        sw << "ZONE  NODES=" << nNodes << ", ELEMENTS=" << nCells << ", DATAPACKING=BLOCK,ZONETYPE= FETETRAHEDRON , VARLOCATION=([" << ss << "]=CELLCENTERED)" << endl;
        // 输出点

        // x
        for (long i = 0; i < nNodes; i++)
            sw << _Nodes[i].x / mm << endl;
        // y
        for (long i = 0; i < nNodes; i++)
            sw << _Nodes[i].y / mm << endl;
        // z
        for (long i = 0; i < nNodes; i++)
            sw << _Nodes[i].z / mm << endl;

        // 输出结果
        for (size_t i = 0; i < data.size(); i++)
        {
            TOther *f = data[i]->fValues.data.get();

            for (long j = 0; j < nCells; j++)
                sw << f[j] << " " << endl;
        }

        for (long k = 0; k < nCells; k++)
        {
            // cell结构,4个点
            map<long, int> pts;
            // 面元
            auto fs = _Cells[k];
            for (int j = 0; j < fs.size(); j++)
            {
                // 点
                auto f = _Faces[fs[j]];
                for (long p : f)
                    pts[p] = 1;
            }
            auto ib = pts.begin();
            auto ie = pts.end();
            while (ib != ie)
            {
                sw << ib->first + 1 << " ";
                ib++;
            }
            sw << endl;
        }
        sw.close();
    }

    template <class TOther>
    void outputTecPlotFEM(string fileName, string title, std::shared_ptr<TOther[]> &data)
    {
        outputTecPlotFEM(fileName, title, data.get());
    }

    template <class TOther>
    void outputTecPlotFEM(string fileName, string title, TOther *data)
    {
        // cout << "save data " << fileName << endl;
        // stopWatch();

        ofstream sw(fileName, ios::out);

        sw << "VARIABLES = \"X\", \"Y\", \"Z\"," << title << "\n";
        string ss("4");
        sw << "ZONE  NODES=" << nNodes << ", ELEMENTS=" << nCells << ", DATAPACKING=BLOCK,ZONETYPE= FETETRAHEDRON , VARLOCATION=([" << ss << "]=CELLCENTERED)" << "\n";

        // 输出点
        if (vertexStructuresFEM.size() <= 0)
        {
            vertexStructuresFEM.reserve(nNodes * 3);
            // x
            for (long i = 0; i < nNodes; i++)
                vertexStructuresFEM.push_back(to_string(_Nodes[i].x / mm));
            // y
            for (long i = 0; i < nNodes; i++)
                vertexStructuresFEM.push_back(to_string(_Nodes[i].y / mm));
            // z
            for (long i = 0; i < nNodes; i++)
                vertexStructuresFEM.push_back(to_string(_Nodes[i].z / mm));
            // cout << stopWatch() << "end\n";
        }
        // 输出结果
        vector<string> dataLines;
        dataLines.reserve(nCells);
        {
            for (long j = 0; j < nCells; j++)
                dataLines.push_back(to_string(data[j]));
        }
        // cout << stopWatch() << "end\n";
        if (cellStructuresFEM.size() <= 0)
        {
            cellStructuresFEM.reserve(nCells);
            set<long> pts;
            for (long k = 0; k < nCells; k++)
            {
                // cell结构,4个点
                pts.clear();
                auto fs = _Cells[k];
                for (int j = 0; j < fs.size(); j++)
                {
                    // 点
                    auto f = _Faces[fs[j]];
                    for (long p : f)
                        pts.insert(p);
                }

                string s;
                for (auto it = pts.begin(); it != pts.end(); it++)
                    s += to_string(*it + 1) + " ";
                cellStructuresFEM.push_back(s);
            }
        }

        // cout << stopWatch() << "end\n";
        for (auto s : vertexStructuresFEM)
            sw << s << "\n";
        for (auto s : dataLines)
            sw << s << "\n";
        for (auto s : cellStructuresFEM)
            sw << s << "\n";

        sw.close();

        // cout << stopWatch() << "end\n";
    }
    template <class TOther>
    void outputTecPlotCellCenter(string fileName, string title, TOther *data)
    {
        // cout << "save data " << fileName << endl;
        // stopWatch();

        ofstream sw(fileName, ios::out);

        sw << "VARIABLES = \"X\", \"Y\", \"Z\",\"" << title << "\"\n";
        string ss("4");
        sw << "ZONE  NODES=" << nNodes << ", ELEMENTS=" << nCells << ", F=FEPOINT,ET= TETRAHEDRON  " << "\n";

        // 输出点
        if (vertexStructuresCellCenter.size() <= 0)
        {
            vertexStructuresCellCenter.reserve(nNodes);
            for (long i = 0; i < nNodes; i++)
            {
                vertexStructuresCellCenter.push_back(to_string(_Nodes[i].x / mm) + " " + to_string(_Nodes[i].y / mm) + " " + to_string(_Nodes[i].z / mm) + " ");
            }
        }
        vector<string> dataLines;
        dataLines.reserve(nNodes);
        for (long i = 0; i < nNodes; i++)
        {
            ostringstream stream;stream << vertex2Cell[i].calcuValue(data);
            dataLines.push_back(stream.str());
        }
        // cout << stopWatch() << "end\n";
        if (cellStructuresFEM.size() <= 0)
        {
            cellStructuresFEM.reserve(nCells);
            set<long> pts;
            for (long k = 0; k < nCells; k++)
            {
                // cell结构,4个点
                pts.clear();
                auto fs = _Cells[k];
                for (int j = 0; j < fs.size(); j++)
                {
                    // 点
                    auto f = _Faces[fs[j]];
                    for (long p : f)
                        pts.insert(p);
                }

                string s;
                for (auto it = pts.begin(); it != pts.end(); it++)
                    s += to_string(*it + 1) + " ";
                cellStructuresFEM.push_back(s);
            }
        }

        // cout << stopWatch() << "end\n";
        for (long i = 0; i < nNodes; i++)
            sw << vertexStructuresCellCenter[i] << dataLines[i] << "\n";
        for (auto s : cellStructuresFEM)
            sw << s << "\n";

        sw.close();

        // cout << stopWatch() << "end\n";
    }
};
