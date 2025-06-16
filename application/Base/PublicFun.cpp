#include "publicFun.h"


#include <vector>
#include <algorithm>
#include "math.h"
#include <chrono>

double m_NearZero = 1e-20;

string stopWatch()
{
    static auto start = std::chrono::high_resolution_clock::now();
    auto finish = std::chrono::high_resolution_clock::now(); // 结束计时
    auto delta= std::chrono::duration_cast<std::chrono::milliseconds>(finish - start).count();
    start=finish;
    return "耗时为:"+to_string(delta)+"ms";
}
