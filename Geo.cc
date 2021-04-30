#include <iostream>
#include <iomanip>
#include <cstring>
#include <vector_types.h>

#include "NP.hh"
#include "sutil_vec_math.h"

#include "CSGPrim.h"
#include "CSGFoundry.h"

#include "Util.h"
#include "Geo.h"
#include "Grid.h"

Geo::Geo(CSGFoundry* foundry_)
    :
    foundry(foundry_),
    top("i0")
{
    init();
}

void Geo::init()
{
    float tminf(0.1) ; 
    float tmaxf(10000.f) ; 

    std::string geometry = Util::GetEValue<std::string>("GEOMETRY", "sphere_containing_grid_of_spheres" ); 
    unsigned layers = Util::GetEValue<unsigned>("LAYERS", 1) ; 

    std::cout
        << "Geo::init"
        << " geometry " << geometry
        << " layers " << layers 
        << std::endl 
        ;    

    if(strcmp(geometry.c_str(), "sphere_containing_grid_of_spheres") == 0)
    {
        init_sphere_containing_grid_of_spheres(tminf, tmaxf, layers );
    }
    else if(strcmp(geometry.c_str(), "parade") == 0)
    {
        init_parade(tminf, tmaxf);
    }
    else if(Util::StartsWith(geometry.c_str(), "clustered_"))
    {
        init_clustered( geometry.c_str() + strlen("clustered_"), tminf, tmaxf ); 
    }
    else if(strcmp(geometry.c_str(), "layered_sphere") == 0 )
    {
        init_layered("sphere", tminf, tmaxf, layers);
    }
    else if(strcmp(geometry.c_str(), "layered_zsphere") == 0 )
    {
        init_layered("zsphere", tminf, tmaxf, layers);
    }
    else
    {
        init(geometry.c_str(), tminf, tmaxf); 
    }

    const float4 ce  = getCenterExtent(); 
    tmin = ce.w*tminf ; 
    tmax = ce.w*tmaxf ; 
    std::cout 
        << "Geo::init" 
        << " top_extent " << ce.w  
        << " tminf " << tminf 
        << " tmin " << tmin 
        << " tmaxf " << tmaxf 
        << " tmax " << tmax 
        << std::endl 
        ; 

    float e_tminf = Util::GetEValue<float>("TMIN", -1.0) ; 
    if(e_tminf > 0.f )
    {
        tmin = ce.w*e_tminf ; 
        std::cout << "Geo::init e_tminf TMIN " << e_tminf << " override tmin " << tmin << std::endl ; 
    }
    
    float e_tmaxf = Util::GetEValue<float>("TMAX", -1.0) ; 
    if(e_tmaxf > 0.f )
    {
        tmax = ce.w*e_tmaxf ; 
        std::cout << "Geo::init e_tmaxf TMAX " << e_tmaxf << " override tmax " << tmax << std::endl ; 
    }
}

/**
Geo::init_sphere_containing_grid_of_spheres
---------------------------------------------

A cube of side 1 (halfside 0.5) has diagonal sqrt(3):1.7320508075688772 
that will fit inside a sphere of diameter sqrt(3) (radius sqrt(3)/2 : 0.86602540378443)
Container sphere "extent" needs to be sqrt(3) larger than the grid extent.

**/

void Geo::init_sphere_containing_grid_of_spheres(float& tminf, float& tmaxf, unsigned layers )
{
    std::cout << "Geo::init_sphere_containing_grid_of_spheres : layers " << layers << std::endl ; 

    unsigned ias_idx = 0 ; 
    unsigned num_solid = 3 ; 
    const float4 ce = Grid::AddInstances(foundry, ias_idx, num_solid) ; 

    float big_radius = float(ce.w)*sqrtf(3.f) ;
    std::cout << " big_radius " << big_radius << std::endl ; 

    foundry->makeLayered("sphere", 0.7f, layers ); 
    foundry->makeLayered("sphere", 1.0f, layers ); 
    CSGSolid* so = foundry->makeLayered("sphere", big_radius, 1 ); 

    top = strdup("i0") ; 

    setCenterExtent(so->center_extent); 

    tminf = 0.75f ; 
    tmaxf = 10000.f ; 
}

void Geo::init_parade(float& tminf, float& tmaxf )
{
    foundry->makeDemoSolids();  
    unsigned num_solid = foundry->getNumSolid() ; 

    unsigned ias_idx = 0 ; 
    const float4 ce = Grid::AddInstances( foundry, ias_idx, num_solid ); 

    top = strdup("i0") ; 
    setCenterExtent(ce); 

    tminf = 0.75f ;   
    tmaxf = 10000.f ; 
}

/**
Geo::init_clustered
--------------------

Aiming to test a GAS containing multiple spread (non-concentric) 
placements of the same type of single node Prim.  
Will need to assign appropriate node transforms and get those applied 
to the bbox at Prim+Node(?) level.

**/
void Geo::init_clustered(const char* name, float& tminf, float& tmaxf )
{
    float unit = Util::GetEValue<float>("CLUSTERUNIT", 200.f ); 
    std::string clusterspec = Util::GetEValue<std::string>("CLUSTERSPEC","-1:2:1,-1:2:1,-1:2:1") ; 

    std::cout 
        << "Geo::init_clustered " << name 
        << " clusterspec " << clusterspec 
        << " unit " << unit 
        << std::endl
        ; 

    bool inbox = false ; 
    std::array<int,9> cl ; 
    Util::ParseGridSpec(cl, clusterspec.c_str()); // string parsed into array of 9 ints 
    CSGSolid* so = foundry->makeClustered(name, cl[0],cl[1],cl[2],cl[3],cl[4],cl[5],cl[6],cl[7],cl[8], unit, inbox ); 
    std::cout << "Geo::init_layered" << name << " so.center_extent " << so->center_extent << std::endl ; 

    setCenterExtent(so->center_extent); 
    top = strdup("g0") ; 

    tminf = 1.60f ;  
    tmaxf = 10000.f ; 
}

void Geo::init_layered(const char* name, float& tminf, float& tmaxf, unsigned layers)
{
    CSGSolid* so = foundry->makeLayered(name, 100.f, layers ); 
    std::cout << "Geo::init_layered" << name << " so.center_extent " << so->center_extent << std::endl ; 
    setCenterExtent(so->center_extent); 
    top = strdup("g0") ; 

    tminf = 1.60f ; 
    tmaxf = 10000.f ; 
}

void Geo::init(const char* name, float& tminf, float& tmaxf )
{
    CSGSolid* so = foundry->make(name) ; 
    std::cout << "Geo::init" << name << " so.center_extent " << so->center_extent << std::endl ; 
    setCenterExtent(so->center_extent); 
    top = strdup("g0") ; 

    tminf = 1.60f ;   //  hmm depends on viewpoint, aiming to cut into the sphere with the tmin
    tmaxf = 10000.f ; 
}


std::string Geo::desc() const
{
    std::stringstream ss ; 
    ss << "Geo " ;
    std::string s = ss.str(); 
    return s ; 
}

void Geo::setCenterExtent(const float4&  center_extent_){ center_extent = center_extent_ ;  }
float4 Geo::getCenterExtent() const  { return center_extent ;  }


unsigned        Geo::getNumSolid() const {                           return foundry->getNumSolid() ;      }
unsigned        Geo::getNumPrim() const {                            return foundry->getNumPrim() ;      }
const CSGSolid*    Geo::getSolid(         unsigned solidIdx) const { return foundry->getSolid(solidIdx);  }
CSGPrimSpec        Geo::getPrimSpec(      unsigned solidIdx) const { return foundry->getPrimSpec(solidIdx);  }
const CSGPrim*     Geo::getPrim(          unsigned primIdx) const  { return foundry->getPrim(primIdx);  }



void Geo::write(const char* dir) const 
{
    std::cout << "Geo::write " << dir << std::endl ;  

    NP uspec("<u4", 5); 
    unsigned* u = uspec.values<unsigned>() ; 
    *(u+0) = getNumSolid(); 
    *(u+1) = 0 ; //getNumGrid(); 
    *(u+2) = 0 ; //InstanceId::ins_bits ; 
    *(u+3) = 0 ; //InstanceId::gas_bits ; 
    *(u+4) = Util::Encode4(top) ; 
    uspec.save(dir, "uspec.npy"); 

    NP fspec("<f4", 4); 
    float* f = fspec.values<float>() ; 
    *(f+0) = center_extent.x ; 
    *(f+1) = center_extent.y ; 
    *(f+2) = center_extent.z ; 
    *(f+3) = center_extent.w ; 
    fspec.save(dir, "fspec.npy"); 

    // with foundry write everything at once, not individual solids/grids
    foundry->write(dir, "foundry"); 

}
