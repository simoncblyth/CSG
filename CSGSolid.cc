
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
    std::string label4(label, 4); 
    std::stringstream ss ; 
    ss << "CSGSolid " 
       << std::setw(5) << label4 
       << " primNum/Offset " 
       << std::setw(5) << numPrim 
       << std::setw(5) << primOffset
       << " ce " << center_extent
       ; 
    std::string s = ss.str(); 
    return s ; 
}

#endif

