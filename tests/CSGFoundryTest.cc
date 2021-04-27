// ./CSGFoundryTest.sh

#include <iostream>
#include <cassert>

#include "sutil_vec_math.h"
#include "CSGFoundry.h"


void test_layered()
{
    CSGFoundry fd ;  

    CSGSolid* s0 = fd.makeLayered("sphere", 100.f, 10 ); 
    CSGSolid* s1 = fd.makeLayered("sphere", 1000.f, 10 ); 
    CSGSolid* s2 = fd.makeLayered("sphere", 50.f, 5 ); 
    CSGSolid* s3 = fd.makeSphere() ; 

    fd.dump(); 

    assert( fd.getSolidIdx(s0) == 0 ); 
    assert( fd.getSolidIdx(s1) == 1 ); 
    assert( fd.getSolidIdx(s2) == 2 ); 
    assert( fd.getSolidIdx(s3) == 3 ); 

    fd.write("/tmp", "FoundryTest_" ); 
}

void test_PrimSpec()
{
    CSGFoundry fd ; 
    fd.makeDemoSolids(); 
    for(unsigned i = 0 ; i < fd.solid.size() ; i++ )
    {
        unsigned solidIdx = i ; 
        std::cout << "solidIdx " << solidIdx << std::endl ; 
        CSGPrimSpec ps = fd.getPrimSpec(solidIdx);
        ps.dump(""); 
    }

    std::string bmap = fd.getBashMap(); 
    std::cout  << bmap << std::endl ; 
}

void test_addTran()
{
    CSGFoundry fd ; 
    const Tran<double>* tr = Tran<double>::make_translate( 100., 200., 300. ) ; 
    unsigned idx = fd.addTran( *tr );   // this idx is 0-based 
    std::cout << "test_addTran idx " << idx << std::endl ; 
    assert( idx == 0u );   
    const qat4* t = fd.getTran(idx) ; 
    const qat4* v = fd.getItra(idx) ; 
 
    std::cout << "idx " << idx << std::endl ; 
    std::cout << "t" << *t << std::endl ; 
    std::cout << "v" << *v << std::endl ; 
}

void test_makeClustered()
{
    std::cout << "[test_makeClustered" << std::endl ; 
    CSGFoundry fd ; 
    bool inbox = false ; 
    fd.makeClustered("sphe", -1,2,1, -1,2,1, -1,2,1, 1000., inbox ); 
    fd.dumpPrim(0); 
    std::cout << "]test_makeClustered" << std::endl ; 
}





int main(int argc, char** argv)
{
    test_layered(); 
    test_PrimSpec(); 
    test_addTran(); 
    test_makeClustered(); 

    return 0 ; 
}
