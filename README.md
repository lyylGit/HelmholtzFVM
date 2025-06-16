# HelmholtzFVM 求解器

[![License](https://img.shields.io/badge/License-HIT-orange.svg)](LICENSE)  

> 基于有限体积法（FVM）的频域亥姆霍兹方程求解器，支持频率相关阻抗边界下的声学模态分析/ A frequency-dependent impedance boundary acoustic modal solver based on the finite volume method

---

## 目录
1. [核心功能](#核心功能)  
2. [安装指南](#安装指南)  
3. [测试算例](#测试算例)  
4. [程序I/O规范](#程序io规范)  
5. [贡献指南](#贡献指南)  
6. [开源许可证](#开源许可证)  

---

## 1. 核心功能
求解频域亥姆霍兹方程：

	$$ \nabla\left ( \frac{1}{\rho_{0}} \nabla \hat{p}(x,\omega)\right ) +\frac{\omega ^2}{\gamma p_{0}}\hat{p}(x,\omega)=0 $$

- **频率相关阻抗边界条件**（声衬/共振腔等）声学模态求解  
- 输出声压场、模态分布
- 网格类型支持 fluent mesh 文本格式与二进制格式; 支持Fluent Meshing生成的三角形、多面体、以及六面体-多面体混合网格。

---

## 2. 安装指南

### 依赖项
开发平台： Ubuntu 18.04

| 工具库   | 版本         | 功能说明            |
| -------- | ------------ | ------------------- |
| IntelMKL | 2024.2.1.100 | 高度优化的数学函数库 |
| FEAST | 4.0 | FEAST 算法是一个通用的特征值求解器 |
| Eigen | 3.4.0        | 线性代数求解器接口  |
| muParserX |  | 公式解析  |
| nlohmann/json |  | 配置文件解析 |

### 编译步骤

下载
```bash
# 克隆仓库并进入目录
git clone https://github.com/lyylGit/HelmholtzFVM.git
cd HelmholtzFVM
```

编译，Makefile 文件如下：

```makefile
CXX =g++ # g++-13

# define any compile-time flags 
CXXFLAGS	:= -std=c++20   -fopenmp -w -g

# define library paths in addition to /usr/lib
LFLAGS = -lfeast -lgfortran -lgomp -lpthread -lm -lmkl_rt -lmetis -lGKlib -lmuparserx

# define output directory
OUTPUT	:= /YOUR/PATH/output

# define source directory
SRC		:= application 
# application/Base 

ifeq ($(MKLROOT),)
MKLROOT := /opt/intel/oneapi/mkl/2024.2
endif

ifeq ($(FEASTROOT),)
FEASTROOT := /YOUR/PATH/FEAST/4.0
endif

# define include directory
INCLUDE	:= src/Base
INCLUDE	+= src/Octree
INCLUDE	+= src/ProjectDefine
INCLUDE	+= src/User/CFD
INCLUDE	+= src/User/Modal
INCLUDE	+= src/User/ByConfigDef 
INCLUDE	+= $(MKLROOT)/include
INCLUDE	+= $(FEASTROOT)/include

INCLUDE	+= /usr/local/include/nlohmann
INCLUDE	+= /usr/local/include/muparserx
INCLUDE	+= /usr/include/eigen3/


# define lib directory
LIB := lib \
		$(MKLROOT)/lib/intel64 \
		$(FEASTROOT)/lib/x64 \
		/YOUR/PATH/local/lib \
		/usr/local/lib \
#LIB += ../muparserx/build

#define blaze
ifeq ($(OS),Windows_NT)
MAIN	:= HelmholtzFVM.exe
SOURCEDIRS	:= $(SRC) 
INCLUDEDIRS	:= $(INCLUDE) $(EIGEN)
LIBDIRS		:= $(LIB)
FIXPATH = $(subst /,\,$1)
RM			:= del /q /f
MD	:= mkdir
else
MAIN	:= HelmholtzFVM
SOURCEDIRS	:= $(shell find $(SRC) -type d)
INCLUDEDIRS	:= $(shell find $(INCLUDE) -type d)
LIBDIRS		:= $(LIB)
FIXPATH = $1
RM = rm -f
MD	:= mkdir -p
endif

# define any directories containing header files other than /usr/include
INCLUDES	:= $(patsubst %,-I%, $(INCLUDEDIRS:%/=%))

# define the C libs
LIBS		:= $(patsubst %,-L%, $(LIBDIRS:%/=%))

# define the C source files
SOURCES		:= $(wildcard $(patsubst %,%/*.cpp, $(SOURCEDIRS)))

# define the C object files 
OBJECTS		:= $(SOURCES:.cpp=.o)

#
# The following part of the makefile is generic; it can be used to 
# build any executable just by changing the definitions above and by
# deleting dependencies appended to the file from 'make depend'
#

OUTPUTMAIN	:= $(call FIXPATH,$(OUTPUT)/$(MAIN))

all: $(OUTPUT) $(MAIN)
	@echo Executing 'all' complete!

$(OUTPUT):
	$(MD) $(OUTPUT)

$(MAIN): $(OBJECTS) 
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(OUTPUTMAIN) $(OBJECTS) $(LFLAGS) $(LIBS)

# this is a suffix replacement rule for building .o's from .c's
# it uses automatic variables $<: the name of the prerequisite of
# the rule(a .c file) and $@: the name of the target of the rule (a .o file) 
# (see the gnu make manual section about automatic variables)
.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $<  -o $@

.PHONY: clean
clean:
	$(RM) $(OUTPUTMAIN)
	$(RM) $(call FIXPATH,$(OBJECTS))
	@echo Cleanup complete!

run: all
	./$(OUTPUTMAIN)
	@echo Executing 'run: all' complete!
```

---

## 3 测试算例

矩形管道固壁边界条件：

```bash
./YOUR/PATH/HelmholtzFVM /YOUR/PATH/Case/test_all_rigid_bc.json
```

网格类型兼容性测试：

```bash
./YOUR/PATH/HelmholtzFVM /YOUR/PATH/Case/test_all_rigid_meshtype.json 
```

阻抗边界测试：

```bash
./YOUR/PATH/HelmholtzFVM /YOUR/PATH/Case/test_fdi_bc.json
```

---

## 4 程序I/O规范

### 输入文件

josn文件，参考测试算例。josn文件除定义了casename、网格文件路径、求解器类型等基本信息外，也包括边界条件定义。

```json
"//": "================脚本定义中用的变量,$** 来引用",
	"param": [
		{
			"name": "rho_0",
			"value": 1.225,
			"memo": "密度1"
		},
		{
			"name": "rho_1",
			"value": 0.3,
			"memo": "密度2"
		},
		{
			"name": "c0",
			"value": 347,
			"memo": "声速1"
		},
		{
			"name": "c1",
			"value": 694.36,
			"memo": "声速2"
		},
		{
			"//": "================计算声阻抗用"
		},
		{
			"name": "sigma",
			"value": -1000,
			"memo": "阻抗边界计算参数"
		},
		{
			"$Hz": "================保留，表示计算中引用系统中的频率"
		}
	],
	"//": "================定义初始场",
	"fields": [
		{
			"name": "rho",
			"type": "uniform",
			"value": "$rho_0"
		},
		{
			"name": "p0",
			"type": "uniform",
			"value": 1e5
		},
		{
			"name": "gamma",
			"type": "uniform",
			"value": 1.44
		}
	],
	"bcs": [
		{
			"//": "================网格文件的边界ID",
			"id": [
				38
			],
			"type": "MixedModal",
			"//": "================计算模态时，拟合阶数",
			"order": 2,
			"//": "================阻抗计算",
			"Z": "$rho_1*$c1 / (-1 + 1i * $sigma /( 2 *  3.14159 * $Hz))"
		},
		{
			"type": "Gradient",
			"//": "================梯度为0+0i，固壁",
			"grad": [
				0,
				0
			]
		}
	],
	"//": "================方程变量,uniform：均匀介质，field：输入密度场",
	"equations": [
		{
			"name": "HelmholtzEq"
		}
	],
```



### 输入文件

| 求解器类型 | 输出文件                               |
| ---------- | -------------------------------------- |
| HelmFreq   | 监测点位置声压（文本文件）以及不同频率下的声压分布Tecplot格式） |
| HelmModal   | 模态主频与声压分布（Tecplot格式） |

