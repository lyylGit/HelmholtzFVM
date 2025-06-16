#include <mkl.h>
#include "stdio.h"
using namespace std;

void testMKL()
{
    MKLVersion mkl_version;
    mkl_get_version(&mkl_version);
    printf("You are using oneMKL %d.%d\n", mkl_version.MajorVersion, mkl_version.UpdateVersion);
}