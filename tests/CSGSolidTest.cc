// ./CSGSolidTest.sh

#include "sutil_vec_math.h"

#include "CSGSolid.h"
#include "NP.hh"
#include <iostream>

int main(int argc, char** argv)
{
    const char* path = argc > 1 ? argv[1] : "/tmp/CSGSolidTest.npy" ; 

    CSGSolid r =  CSGSolid::Make("red"  ,    1, 0 ) ; r.center_extent = {1.f, 1.f, 1.f, 1.f } ;
    CSGSolid g =  CSGSolid::Make("green",    1, 1 ) ; r.center_extent = {2.f, 2.f, 2.f, 2.f } ;
    CSGSolid b =  CSGSolid::Make("blue",     1, 2 ) ; r.center_extent = {2.f, 2.f, 2.f, 2.f } ;
    CSGSolid c =  CSGSolid::Make("cyan",     1, 3 ) ; r.center_extent = {3.f, 3.f, 3.f, 3.f } ;
    CSGSolid m =  CSGSolid::Make("magenta",  1, 4 ) ; r.center_extent = {4.f, 4.f, 4.f, 4.f } ;
    CSGSolid y =  CSGSolid::Make("yellow",   1, 5 ) ; r.center_extent = {5.f, 5.f, 5.f, 5.f } ;


    std::vector<CSGSolid> so ; 

    so.push_back(r); 
    so.push_back(g); 
    so.push_back(b); 
    so.push_back(c); 
    so.push_back(m); 
    so.push_back(y); 

    std::cout << "sizeof(CSGSolid)" << sizeof(CSGSolid) << std::endl ; 
    assert( sizeof(float) == sizeof(int));
    assert( sizeof(CSGSolid) == 8*sizeof(float) ); 

    unsigned num_items = so.size(); 
    assert( num_items == 6 );
    unsigned num_values = sizeof(CSGSolid)/sizeof(int) ; 
    assert( num_values == 8 ); 

    NP::Write( path, (int*)so.data(), num_items, num_values ) ;
    return 0 ; 
}
