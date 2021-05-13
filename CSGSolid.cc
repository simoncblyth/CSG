
#include "sutil_vec_math.h"
#include "CSGSolid.h"

#if defined(__CUDACC__) || defined(__CUDABE__)
#else

#include <sstream>
#include <iostream>
#include <iomanip>
#include <cstring>



CSGSolid CSGSolid::Make( const char* label_, int numPrim_, int primOffset_ )
{
    CSGSolid so = {} ; 

    strncpy( so.label, label_, sizeof(label) );
    so.numPrim = numPrim_ ; 
    so.primOffset = primOffset_ ; 
    so.type = STANDARD_SOLID ;  
    so.center_extent = make_float4(0.f, 0.f, 0.f, 0.f) ;  // changed later 

    return so ; 
}

std::string CSGSolid::desc() const 
{
    std::string label8(label, 8); 
    std::stringstream ss ; 
    ss << "CSGSolid " 
       << std::setw(9) << label8 
       << " primNum/Offset " 
       << std::setw(5) << numPrim 
       << std::setw(5) << primOffset
       << " ce " << center_extent
       ; 

    if( type == ONE_PRIM_SOLID )
    {
       ss << " ONE_PRIM_SOLID "
          << " origin_solidIdx " << origin_solidIdx 
          << " origin_primIdxRel " << origin_primIdxRel 
          << " origin_nodeIdxRel " << origin_nodeIdxRel 
          ;
    }
    std::string s = ss.str(); 
    return s ; 
}

std::string CSGSolid::MakeLabel(char type, unsigned idx )
{
    std::stringstream ss ; 
    //ss << type << std::setfill('0') << std::setw(3) << idx ; 
    ss << type << idx ; 
    std::string s = ss.str();  
    return s ; 
}



#endif

