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

    const char* cfbase = SSys::getenvvar("CFBASE", "/tmp" );  
    const char* rel = "CSGFoundry" ; 

    foundry.write(cfbase, rel );    // expects existing directory $CFBASE/CSGFoundry 

    CSGFoundry* fd = CSGFoundry::Load(cfbase, rel);  // load foundary and check identical bytes
    assert( 0 == CSGFoundry::Compare(&foundry, fd ) );  


    return 0 ; 
}
