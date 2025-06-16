#pragma once
#include <stdint.h>
#include <vector>
#include <unordered_map>
#include <numeric>

// 稀疏矩阵求解器
//https://stackoverflow.com/questions/20511347/a-good-hash-function-for-a-vector/72073933#72073933
// 定义结构
using index_type_int32 = std::vector<uint32_t>;

struct index_hash
{
    // std::size_t operator()(index_type const &i) const noexcept
    // {
    //     // 确定是2个坐标
    //     std::size_t r = i[0];
    //     int x = i[1];
    //     r += std::hash<int>()(x) + 0x9e3779b9 + (r << 6) + (r >> 2);
    //     return r;
    // }
    std::size_t operator()(index_type_int32 const &vec) const
    {
        std::size_t seed = vec.size();
        for (auto x : vec)
        {
            x = ((x >> 16) ^ x) * 0x45d9f3b;
            x = ((x >> 16) ^ x) * 0x45d9f3b;
            x = (x >> 16) ^ x;
            seed ^= x + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

template <typename T>
using sparse_array = std::unordered_map<index_type_int32, T, index_hash>;

//64bit
using index_type_int64 = std::vector<uint64_t>;

struct index_hash_64
{
    std::size_t operator()(index_type_int64 const &vec) const
    {
        std::size_t seed = vec.size();
        for (auto x : vec)
        {
            x = (x ^ (x >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
            x = (x ^ (x >> 27)) * UINT64_C(0x94d049bb133111eb);
            x = x ^ (x >> 31);
            seed ^= x + 0x9e3779b97f4a7c55 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

template <typename T>
using sparse_array_64 = std::unordered_map<index_type_int64, T, index_hash_64>;
