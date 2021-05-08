#include "SSys.hh"
#include "sutil_vec_math.h"
#include "CSGFoundry.h"

#include "OPTICKS_LOG.hh"

int main(int argc, char** argv)
{
    unsigned mesh_idx = argc > 1 ? std::atoi(argv[1]) : 130 ; 

    OPTICKS_LOG(argc, argv); 

    const char* cfbase = SSys::getenvvar("CFBASE", "$TMP/CSG_GGeo" );  
    CSGFoundry* fd = CSGFoundry::Load(cfbase, "CSGFoundry"); 

    LOG(info) << "foundry " << fd->desc() ; 
    fd->summary(); 


    std::vector<CSGPrim> prim ; 
    fd->getMeshPrim(prim, mesh_idx ); 

    LOG(info) 
        << " mesh_idx " << mesh_idx 
        << " prim.size " << prim.size()
        ; 
    
    for(unsigned i=0 ; i < prim.size() ; i++)
    {
        const CSGPrim& pr = prim[i]; 
        std::cout << pr.desc() << std::endl ; 
    }


    return 0 ; 
}



