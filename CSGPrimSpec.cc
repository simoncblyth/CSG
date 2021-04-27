
#include <iostream>
#include <iomanip>
#include <cassert>
#include <cstring>

#include "CSGPrimSpec.h"


void CSGPrimSpec::gather(std::vector<float>& out) const 
{
    assert( device == false ); 
    unsigned size_in_floats = 6 ; 
    out.resize( num_prim*size_in_floats ); 

    unsigned stride_in_floats = stride_in_bytes/sizeof(float) ; 
    for(unsigned i=0 ; i < num_prim ; i++) 
    {   
        float* dst = out.data() + size_in_floats*i ;   
        const float* src = aabb + stride_in_floats*i ;   
        memcpy(dst, src,  sizeof(float)*size_in_floats );  
    }   
}

void CSGPrimSpec::Dump(std::vector<float>& out)  // static 
{
     std::cout << " gather " << out.size() << std::endl ; 
     for(unsigned i=0 ; i < out.size() ; i++) 
     {    
         if(i % 6 == 0) std::cout << std::endl ; 
         std::cout << std::setw(10) << out[i] << " " ; 
     } 
     std::cout << std::endl ; 
}


void CSGPrimSpec::dump(const char* msg) const 
{
    assert( stride_in_bytes % sizeof(float) == 0 ); 
    unsigned stride_in_floats = stride_in_bytes/sizeof(float) ; 
    std::cout 
        << msg 
        << " num_prim " << num_prim 
        << " stride_in_bytes " << stride_in_bytes 
        << " stride_in_floats " << stride_in_floats 
        << std::endl 
        ; 

    for(unsigned i=0 ; i < num_prim ; i++)
    {   
        std::cout 
            << " i " << std::setw(4) << i 
            << " sbtIndexOffset " << std::setw(4) << *(sbtIndexOffset + i*stride_in_floats)   
            ; 
        for(unsigned j=0 ; j < 6 ; j++)  
            std::cout << std::setw(10) << std::fixed << std::setprecision(3) << *(aabb + i*stride_in_floats + j ) << " "  ;   
        std::cout << std::endl ; 
    }   
}

