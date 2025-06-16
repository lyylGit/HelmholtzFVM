#include <iostream>
#include <vector>
#include <string>
#include <cstring>
#include <fstream>
#include <sstream>

#include "ProjectDef.h"
// #define CFDSolver

#include "Base/Mesh"
#include "Base/Case"

#include <complex>
// #include "EqMKL.h"
//  #include "EigenFEAST.h"

using namespace std;

void ShowPrompt(string s)
{
   cout << s << endl;
}
void createProjdef(string config_path);

int main(int argc, char* argv[])
{
   if (argc < 2) {
        std::cerr << "请提供配置文件路径作为参数" << std::endl;
        return -1;
    }
   
   string config_path = argv[1];
   createProjdef(config_path);
   return 1;
}

bool SetupProject(ProjectConfig &conf);

void createProjdef(string config_path)
{
   ProjectConfig conf(config_path);  
   if (!conf.parseProj())
      return;

   SetupProject(conf);

}
