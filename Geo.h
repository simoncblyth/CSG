#pragma once
#include <array>
#include <vector>
#include <string>


#include "CSGPrimSpec.h"

struct CSGSolid ; 
struct CSGPrim ; 
struct CSGFoundry ; 

struct Grid ; 


/**
Can Foundry replace Geo ?

* no foundry is for now just for making local solids 
  and holding the constituent Node, Prim etc..

Geo is for high level definition of specific examples of geometry, 
it is fine as a test of the CSG model but it definitely does not belong 
in the model package.

**/

struct Geo
{
    Geo(CSGFoundry* foundry_);

    void init();
    void init_sphere_containing_grid_of_spheres(float& tminf, float& tmaxf, unsigned layers);
    void init_parade(float& tminf, float& tmaxf );
    void init_layered(const char* name, float& tminf, float& tmaxf, unsigned layers);
    void init_clustered(const char* name, float& tminf, float& tmaxf );

    void init(const char* name, float& tminf, float& tmaxf);
    std::string desc() const ;

    unsigned getNumSolid() const ; 
    unsigned getNumPrim() const ; 

    const CSGSolid* getSolid(unsigned solidIdx) const ; 
    CSGPrimSpec getPrimSpec(unsigned solidIdx) const ;
    const CSGPrim* getPrim(unsigned primIdx) const ; 

    unsigned getNumGrid() const ; 
    const Grid*  getGrid(unsigned gridIdx) const ; 
    const Grid*  getGrid_(int gridIdx) const ; 

    void addGrid(const Grid* grid) ;

    void write(const char* prefix) const ; 
    void setCenterExtent(const float4& center_extent); 
    float4 getCenterExtent() const ; 

    float tmin = 0.f ; 
    float tmax = 1e16f ; 
    float4 center_extent = {0.f, 0.f, 0.f, 100.f} ; 

    CSGFoundry*               foundry ; 
    std::vector<const Grid*>  grids ; 
    const char*               top ;  

};


