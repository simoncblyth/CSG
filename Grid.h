#pragma once

#include <array>
#include <vector>
#include <glm/glm.hpp>

/**
Grid
======

Input high level specification example. 
Should not be part of the CSG model.

Currently using glm::mat4 for the transforms but as 
there are no rotations there is not much reason for this. 
Using qat4.h would avoid the Sys.h, due to inherent multi-typing. 

Actually is makes more sense for the foundry to hold
the transforms and for Grid+Geo to just act as an example interface 
to populating the foundry.

The foundry is central to the CSG model.

**/

struct Grid
{ 
    unsigned               ias_idx ; 
    unsigned               num_solid ; 
    float                  gridscale ; 

    std::array<int,9>      grid ; 
    std::vector<unsigned>  solid_modulo ;  
    std::vector<unsigned>  solid_single ;  
    std::vector<glm::mat4> trs ;  

    Grid(unsigned ias_idx, unsigned num_solid);

    const float4 center_extent() const ;
    std::string desc() const ;

    void init(); 
    void write(const char* base, const char* rel, unsigned idx ) const ;
};


