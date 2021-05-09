#include "SSys.hh"
#include "sutil_vec_math.h"
#include "CSGFoundry.h"

#include "OPTICKS_LOG.hh"


int test_getGlobalCenterExtent(const CSGFoundry& fd,  float4& gce, int midx, int mord, int iidx )
{
    std::vector<CSGPrim> prim ; 
    fd.getMeshPrim(prim, midx ); // collect prim matching the MIDX 

    LOG(info) 
        << " MIDX " << midx 
        << " prim.size " << prim.size()
        ; 
    
    for(unsigned i=0 ; i < prim.size() ; i++)
    {
        const CSGPrim& pr = prim[i]; 
        std::cout << pr.desc() << std::endl ; 
    }

    bool mord_in_range = mord < prim.size() ; 
    LOG(info)  
        << " midx " << midx
        << " mord " << mord 
        << " prim.size " << prim.size()
        << " mord_in_range " << mord_in_range 
        ;   

    if(!mord_in_range) return 1 ; 


    // first find the MORD-inal prim which has MIDX for its lvIdx
    const CSGPrim& lpr = prim[mord] ; 
    const float4 local_ce = lpr.ce() ; 

    // use the prim to lookup indices for  the solid and prim 
    unsigned repeatIdx = lpr.repeatIdx(); 
    unsigned primIdx = lpr.primIdx(); 
    unsigned gas_idx = repeatIdx ; 

    // collect the instances 
    std::vector<qat4> inst ; 
    fd.getInstanceTransformsGAS(inst, gas_idx ); 

    bool iidx_in_range = iidx < inst.size(); 

    LOG(info) 
        << " repeatIdx " << repeatIdx
        << " primIdx " << primIdx
        << " inst.size " << inst.size()
        << " iidx " << iidx
        << " iidx_in_range " << iidx_in_range 
        << " local_ce " << local_ce 
        ; 

    if(!iidx_in_range) return 2 ; 

    qat4 q(inst[iidx].cdata());   // copy the instance
    unsigned ins_idx, gas_idx2, ias_idx ; 
    q.getIdentity(ins_idx, gas_idx2, ias_idx )  ;
    q.clearIdentity();           // clear before doing any transforming 
    assert( gas_idx == gas_idx2 ); 

    CSGPrim gpr = {} ; 
    CSGPrim::copy(gpr, lpr); 
    q.transform_aabb_inplace( gpr.AABB_() ); 

    LOG(info) 
        << " q " << q 
        << " ins_idx " << ins_idx
        << " ias_idx " << ias_idx
        ; 

    LOG(info) 
        << " gpr " << gpr.desc()
        ; 


    return 0 ; 
}



void test_CE(const CSGFoundry* fd, int midx, int mord, int iidx)
{
    float4 gce0 = make_float4( 0.f, 0.f, 0.f, 0.f ); 
    int rc0 = test_getGlobalCenterExtent(*fd, gce0, midx, mord, iidx ) ; 
    assert( rc0 == 0 ); 
    LOG(info) << " gce0 " << gce0 ; 

    float4 gce1 = make_float4( 0.f, 0.f, 0.f, 0.f ); 
    int rc1 = fd->getGlobalCenterExtent(gce1, midx, mord, iidx) ;
    assert( rc1 == 0 ); 
    LOG(info) << " gce1 " << gce1 ; 


    float4 gce2 = make_float4( 0.f, 0.f, 0.f, 0.f ); 
    int rc2 = fd->getCenterExtent(gce2, midx, mord, iidx) ;
    assert( rc2 == 0 ); 
    LOG(info) << " gce2 " << gce2 ; 
}



int main(int argc, char** argv)
{
    OPTICKS_LOG(argc, argv); 

    CSGFoundry* fd = CSGFoundry::Load(SSys::getenvvar("CFBASE", "$TMP/CSG_GGeo" ), "CSGFoundry"); 
    LOG(info) << "foundry " << fd->desc() ; 
    //fd->summary(); 

    const char* moi = SSys::getenvvar("MOI", "sWorld:0:0"); 
    int midx, mord, iidx ; 
    fd->parseMOI(midx, mord, iidx,  moi ); 

    LOG(info) 
        << " MOI " << moi 
        << " midx " << midx 
        << " mord " << mord 
        << " iidx " << iidx
        ; 

    test_CE(fd, midx, mord, iidx); 



    return 0 ; 
}




