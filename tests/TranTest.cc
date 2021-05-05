#include <iostream>
#include "Tran.h"

int main()
{
    const Tran<float>* i = Tran<float>::make_translate( 0.f, 0.f, 0.f ) ; 
    const Tran<float>* t = Tran<float>::make_translate( 1.f, 2.f, 3.f ) ; 
    const Tran<float>* s = Tran<float>::make_scale(     1.f, 2.f, 3.f ) ; 
    const Tran<float>* r = Tran<float>::make_rotate(    0.f, 0.f, 1.f, 45.f ) ; 

    const Tran<float>* ts0 = Tran<float>::product( t, s, false ); 
    const Tran<float>* ts1 = Tran<float>::product( t, s, true ); 

    std::cout << "i   " << i->brief() << std::endl ; 
    std::cout << "t   " << t->brief() << std::endl ; 
    std::cout << "s   " << s->brief() << std::endl ; 
    std::cout << "r   " << r->brief() << std::endl ; 
    std::cout << "ts0 " << ts0->brief() << std::endl ; 
    std::cout << "ts1 " << ts1->brief() << std::endl ; 

    return 0 ; 
}
