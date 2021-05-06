#pragma once

#if defined(__CUDACC__) || defined(__CUDABE__)
#else
#include <vector>
#endif

struct CSGPrimSpec
{
    const float*    aabb ; 
    const unsigned* sbtIndexOffset ;  
    unsigned        num_prim ; 
    unsigned        stride_in_bytes ; 
    bool            device ; 

#if defined(__CUDACC__) || defined(__CUDABE__)
#else
    void downloadDump(const char* msg="CSGPrimSpec::downloadDump") const ; 
    void gather(std::vector<float>& out) const ;
    static void Dump(std::vector<float>& out);
    void dump(const char* msg="CSGPrimSpec::Dump") const ; 
#endif
};

