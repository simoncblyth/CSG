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
    char        label[8] ; 
    int         numPrim ; 
    int         primOffset ;

    float4      center_extent ; 

    int         type = STANDARD_SOLID ;       
    int         origin_solidIdx  ;   // these are used to identify extra debugging solids 
    int         origin_primIdxRel ; 
    int         origin_nodeIdxRel ; 


#if defined(__CUDACC__) || defined(__CUDABE__)
#else
    static CSGSolid Make( const char* label_, int numPrim_, int primOffset_=-1 ); 
    static std::string MakeLabel(char typ0, unsigned idx0 );  
    static std::string MakeLabel(char typ0, unsigned idx0, char typ1, unsigned idx1 );  
    static std::string MakeLabel(char typ0, unsigned idx0, char typ1, unsigned idx1, char typ2, unsigned idx2 );  
    std::string desc() const ; 
#endif

};


