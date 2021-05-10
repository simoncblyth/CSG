#include "SSys.hh"
#include "OPTICKS_LOG.hh"

#include "sutil_vec_math.h"
#include "CSGFoundry.h"
#include "DemoGeo.h"


int main(int argc, char** argv)
{
    OPTICKS_LOG(argc, argv); 

    CSGFoundry foundry ;
    DemoGeo dg(&foundry) ;  

    const char* cfbase = SSys::getenvvar("CFBASE", "$TMP/CSGDemoTest" );
    const char* rel = "CSGFoundry" ; 

    foundry.write(cfbase, rel );    // expects existing directory $CFBASE/CSGFoundry 

    CSGFoundry* fd = CSGFoundry::Load(cfbase, rel);  // load foundary and check identical bytes
    assert( 0 == CSGFoundry::Compare(&foundry, fd ) );  


    unsigned ias_idx = 0u ; 
    unsigned long long emm = 0ull ;
    LOG(info) << "descInst" << std::endl << fd->descInst(ias_idx, emm);  

    AABB bb = fd->iasBB(ias_idx, emm ); 
    float4 ce = bb.center_extent() ;  

    LOG(info) << "bb:" << bb.desc() ; 
    LOG(info) << "ce:" << ce ; 

    std::vector<float3> corners ; 
    AABB::cube_corners(corners, ce ); 
    for(int i=0 ; i < int(corners.size()) ; i++) LOG(info) << corners[i] ;  


    return 0 ; 
}
