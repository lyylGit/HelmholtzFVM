#include <iostream>
#include <vector>
#include <string>
#include <cstring>
#include <fstream>
#include <sstream>

#include "ProjectDef.h"

#include "Mesh"
#include "Case"

void ShowPrompt1(string s)
{
    cout << s << endl; 
}


bool calcuProjDouble(ProjectConfig &conf)
{
    Mesh3<double> mesh; 
    // 网格中单位为mm，需要转换为m 
    mesh.mm = conf.mm;
    mesh.fnShowPrompt = ShowPrompt1;
    if (conf.meshFileType == defMeshFileType::FluentBin)
        mesh.ReadMesh_Fluent_Mesh_Binary(conf.meshFile);
    else if (conf.meshFileType == defMeshFileType::FluentTxt)
        mesh.ReadMesh_Fluent_Mesh_Text(conf.meshFile);

    // // 处理网格信息  
    mesh.ConstructCells();

    CaseHelmholtzConfigDef<double> helmholtzCase(&mesh);
    helmholtzCase._ProjectConfig=&conf;
    helmholtzCase.fnShowPrompt = ShowPrompt1;
    helmholtzCase.Define(); 

    if (helmholtzCase.Iterate())
        ShowPrompt1("ok");  
    else 
        ShowPrompt1("err");   

    return true;
}
bool calcuProjFloat(ProjectConfig &conf)
{ 
 
    return true;
}

bool SetupProject(ProjectConfig &conf) 
{
    // 构造计算方案
    if (conf.tCase == defMeshValueType::Double)
        return calcuProjDouble(conf);
    else if (conf.tCase == defMeshValueType::Float)
        return calcuProjFloat(conf);

    return false; 
}
