
CXX =g++ # g++-13

# define any compile-time flags 
CXXFLAGS	:= -std=c++20   -fopenmp -w -g

# define library paths in addition to /usr/lib
LFLAGS = -lfeast -lgfortran -lgomp -lpthread -lm -lmkl_rt -lmetis -lGKlib -lmuparserx

# define output directory
OUTPUT	:= output

# define source directory
SRC		:= application 
# application/Base 

ifeq ($(MKLROOT),)
MKLROOT := /opt/intel/oneapi/mkl/2024.2
endif

ifeq ($(FEASTROOT),)
FEASTROOT := your/path/FEAST/4.0
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

INCLUDE	+= your/path/local/include 
INCLUDE	+= your/path/METIS/build/xinclude  
INCLUDE	+= /usr/local/include/nlohmann
INCLUDE	+= /usr/local/include/muparserx

INCLUDE	+= /usr/include/eigen3/


# define lib directory
LIB := lib \
		$(MKLROOT)/lib/intel64 \
		$(FEASTROOT)/lib/x64 \
		your/path/local/lib \
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