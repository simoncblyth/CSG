
#include "sutil_vec_math.h"
#include "CSGSolid.h"

#if defined(__CUDACC__) || defined(__CUDABE__)
#else

#include <sstream>
#include <iostream>
#include <iomanip>
#include <cstring>


std::string CSGSolid::desc() const 
{
    std::stringstream ss ; 
    ss << "CSGSolid " 
       << std::setw(30) << label 
       << " numPrim " << std::setw(3) << numPrim 
       << " primOffset " << std::setw(3) << primOffset
       << " center_extent " << center_extent
       ; 
    std::string s = ss.str(); 
    return s ; 
}

#endif

