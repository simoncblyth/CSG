// ./ScanTest.sh

#include <vector>
#include <cassert>
#include <iostream>

#include "sutil_vec_math.h"
#include "CSGFoundry.h"
#include "CSGSolid.h"


#include "Scan.h"
#include "Geo.h"
#include "Util.h"
#include "View.h"


void test_Foundry_Scan()
{
    const char* dir = "/tmp/ScanTest_scans" ; 
    CSGFoundry fd ;  

    //fd.makeDemoSolids(); 
    fd.makeEllipsoid(); 
    //const char* name = "sphe" ; 
    //fd.makeClustered(name, 0,1,1, 0,1,1,  0,2,1, 100. );  

    unsigned numSolid = fd.getNumSolid() ; 
    std::cout << "numSolid " << numSolid << std::endl ; 

    for(unsigned i=0 ; i < numSolid ; i++)
    {
        const CSGSolid* solid = fd.getSolid(i); 

        Scan sc(dir, &fd, solid); 
        sc.axis_scan(); 
        sc.rectangle_scan(); 
        sc.circle_scan(); 
    }
}


int main(int argc, char** argv)
{
    CSGFoundry foundry ;  
    Geo geo(&foundry) ; 

    unsigned width = 1280u ; 
    unsigned height = 720u ; 

    glm::vec4 eye_model ; 
    Util::GetEVec(eye_model, "EYE", "-1.0,-1.0,1.0,1.0"); 

    const float4 gce = geo.getCenterExtent() ;   
    glm::vec4 ce(gce.x,gce.y,gce.z, gce.w*1.4f );   // defines the center-extent of the region to view

    View view = {} ; 
    view.update(eye_model, ce, width, height) ; 
    view.dump("View::dump"); 
    view.save("/tmp");  

    const char* dir = "/tmp/ScanTest_scans" ; 

    const CSGSolid* solid0 = foundry.getSolid(0); 

    Scan scan(dir, &foundry, solid0 ); 
    scan.axis_scan() ; 
    scan.rectangle_scan() ; 
    scan.circle_scan() ; 

    return 0 ;  
}
