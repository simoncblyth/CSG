#pragma once

#include "CSGEnum.h"

#if defined(__CUDACC__) || defined(__CUDABE__)
#else
#include <string>
#endif


/**
CSGSolid
==========

Currently this is not uploaded to GPU, but still coding like 
it might be in future, ie keep it simple, no refs, 
128 bit (16 byte) alignment

**/

struct CSGSolid   // Composite shape 
{
    char        label[16] ; 

    int         numPrim ; 
    int         primOffset ;
    int         type = STANDARD_SOLID ;       
    int         padding ;  

    float4      center_extent ; 


#if defined(__CUDACC__) || defined(__CUDABE__)
#else
    static CSGSolid Make( const char* label_, int numPrim_, int primOffset_=-1 ); 
    static std::string MakeLabel(const char* typ0, unsigned idx0, char delim='_' );  
    static std::string MakeLabel(char typ0, unsigned idx0 );  
    static std::string MakeLabel(char typ0, unsigned idx0, char typ1, unsigned idx1 );  
    static std::string MakeLabel(char typ0, unsigned idx0, char typ1, unsigned idx1, char typ2, unsigned idx2 );  
    std::string desc() const ; 
#endif

};


