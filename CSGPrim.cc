
#if defined(__CUDACC__) || defined(__CUDABE__)
#else

#include "sutil_vec_math.h"
#include "qat4.h"
#include "CSGPrim.h"

#include <vector>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cassert>
#include <cstring>


std::string CSGPrim::desc() const 
{  
    std::stringstream ss ; 
    ss 
      << "CSGPrim"
      << " mn " << mn() 
      << " mx " << mx() 
      << " numNode/sbt/node/tran/planOffset " 
      << std::setw(3) << numNode() 
      << std::setw(3) << sbtIndexOffset() 
      << std::setw(3) << nodeOffset()
      << std::setw(3) << tranOffset()
      << std::setw(3) << planOffset()
      ;
    
    std::string s = ss.str(); 
    return s ; 
}


CSGPrimSpec CSGPrim::MakeSpec( const CSGPrim* prim0,  unsigned primIdx, unsigned numPrim ) // static 
{
    const CSGPrim* prim = prim0 + primIdx ; 

    CSGPrimSpec ps ; 
    ps.aabb = prim->AABB() ; 
    ps.sbtIndexOffset = prim->sbtIndexOffsetPtr() ;  
    ps.num_prim = numPrim ; 
    ps.stride_in_bytes = sizeof(CSGPrim); 
    return ps ; 
}

#endif


