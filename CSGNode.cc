#if defined(__CUDACC__) || defined(__CUDABE__)
#else

#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector_types.h>

#include "sutil_vec_math.h"
#include "OpticksCSG.h"
#include "CSGNode.h"


const float CSGNode::UNBOUNDED_DEFAULT_EXTENT = 100.f ; 


std::string CSGNode::desc() const 
{
    const float* aabb = AABB(); 
    std::stringstream ss ; 
    ss
       << "CSGNode "
       << CSG::Tag((OpticksCSG_t)typecode())
       << " aabb: "
       ;    

    for(int i=0 ; i < 6 ; i++ ) ss << std::setw(7) << std::fixed << std::setprecision(1) << *(aabb+i) << " " ; 
    std::string s = ss.str();
    return s ; 
}


void CSGNode::Dump(const CSGNode* n_, unsigned ni, const char* label)
{
    std::cout << "CSGNode::Dump ni " << ni << " " ; 
    if(label) std::cout << label ;  
    std::cout << std::endl ; 

    for(unsigned i=0 ; i < ni ; i++)
    {
        const CSGNode* n = n_ + i ;        
        std::cout << "(" << i << ")" << std::endl ; 
        std::cout 
            << " node.q0.f.xyzw ( " 
            << std::setw(10) << std::fixed << std::setprecision(3) << n->q0.f.x  
            << std::setw(10) << std::fixed << std::setprecision(3) << n->q0.f.y  
            << std::setw(10) << std::fixed << std::setprecision(3) << n->q0.f.z  
            << std::setw(10) << std::fixed << std::setprecision(3) << n->q0.f.w
            << " ) " 
            << std::endl 
            << " node.q1.f.xyzw ( " 
            << std::setw(10) << std::fixed << std::setprecision(3) << n->q1.f.x  
            << std::setw(10) << std::fixed << std::setprecision(3) << n->q1.f.y  
            << std::setw(10) << std::fixed << std::setprecision(3) << n->q1.f.z  
            << std::setw(10) << std::fixed << std::setprecision(3) << n->q1.f.w
            << " ) " 
            << std::endl 
            << " node.q2.f.xyzw ( " 
            << std::setw(10) << std::fixed << std::setprecision(3) << n->q2.f.x  
            << std::setw(10) << std::fixed << std::setprecision(3) << n->q2.f.y  
            << std::setw(10) << std::fixed << std::setprecision(3) << n->q2.f.z  
            << std::setw(10) << std::fixed << std::setprecision(3) << n->q2.f.w
            << " ) " 
            << std::endl 
            << " node.q3.f.xyzw ( " 
            << std::setw(10) << std::fixed << std::setprecision(3) << n->q3.f.x  
            << std::setw(10) << std::fixed << std::setprecision(3) << n->q3.f.y  
            << std::setw(10) << std::fixed << std::setprecision(3) << n->q3.f.z  
            << std::setw(10) << std::fixed << std::setprecision(3) << n->q3.f.w
            << " ) " 
            << std::endl 
            ;

        std::cout 
            << " node.q0.i.xyzw ( " 
            << std::setw(10) << n->q0.i.x  
            << std::setw(10) << n->q0.i.y  
            << std::setw(10) << n->q0.i.z  
            << std::setw(10) << n->q0.i.w
            << " ) " 
            << std::endl 
            << " node.q1.i.xyzw ( " 
            << std::setw(10) << n->q1.i.x  
            << std::setw(10) << n->q1.i.y  
            << std::setw(10) << n->q1.i.z  
            << std::setw(10) << n->q1.i.w
            << " ) " 
            << std::endl 
            << " node.q2.i.xyzw ( " 
            << std::setw(10) << n->q2.i.x  
            << std::setw(10) << n->q2.i.y  
            << std::setw(10) << n->q2.i.z  
            << std::setw(10) << n->q2.i.w
            << " ) " 
            << std::endl 
            << " node.q3.i.xyzw ( " 
            << std::setw(10) << n->q3.i.x  
            << std::setw(10) << n->q3.i.y  
            << std::setw(10) << n->q3.i.z  
            << std::setw(10) << n->q3.i.w
            << " ) " 
            << std::endl 
            ;
    }
}



CSGNode CSGNode::Union(){         return CSGNode::BooleanOperator('U') ; }  // static
CSGNode CSGNode::Intersection(){  return CSGNode::BooleanOperator('I') ; }
CSGNode CSGNode::Difference(){    return CSGNode::BooleanOperator('D') ; }
CSGNode CSGNode::BooleanOperator(char op)   // static 
{
    CSGNode nd = {} ;
    nd.setTypecode(CSG::BooleanOperator(op)) ; 
    return nd ; 
}

void CSGNode::setAABBLocal()
{
    unsigned tc = typecode(); 
    if(tc == CSG_SPHERE)
    {
        float r = radius();  
        setAABB(  -r, -r, -r,  r, r, r  ); 
    }
    else if(tc == CSG_ZSPHERE)
    {
        float r = radius();  
        float z1_ = z1();  
        float z2_ = z2();  
        assert( z2_ > z1_ ); 
        setAABB(  -r, -r, z1_,  r, r, z2_  ); 
    }
    else if( tc == CSG_CONE )
    {
        float r1, z1, r2, z2, a, b ; 
        getParam( r1, z1, r2, z2, a, b ); 
        float rmax = fmaxf(r1, r2) ;
        setAABB( -rmax, -rmax, z1, rmax, rmax, z2 ); 
    }
    else if( tc == CSG_BOX3 )
    {
        float fx, fy, fz, a, b, c ; 
        getParam( fx, fy, fz, a, b, c ); 
        setAABB( -fx*0.5f , -fy*0.5f, -fz*0.5f, fx*0.5f , fy*0.5f, fz*0.5f );   
    }
    else if( tc == CSG_CYLINDER )
    {
        float px, py, a, radius, z1, z2 ; 
        getParam( px, py, a, radius, z1, z2) ; 
        setAABB( px-radius, py-radius, z1, px+radius, py+radius, z2 );   
    }
    else if( tc == CSG_DISC )
    {
        float px, py, ir, r, z1, z2 ; 
        getParam( px, py, ir, r, z1, z2 ); 
        setAABB( px - r , py - r , z1, px + r, py + r, z2 ); 
    }
    else if( tc == CSG_HYPERBOLOID )
    {
        float r0, zf, z1, z2, a, b ; 
        getParam(r0, zf, z1, z2, a, b ) ; 

        assert( z1 < z2 ); 
        const float rr0 = r0*r0 ; 
        const float z1s = z1/zf ; 
        const float z2s = z2/zf ; 

        const float rr1 = rr0 * ( z1s*z1s + 1.f ) ;
        const float rr2 = rr0 * ( z2s*z2s + 1.f ) ;
        const float rmx = sqrtf(fmaxf( rr1, rr2 )) ; 

        setAABB(  -rmx,  -rmx,  z1,  rmx, rmx, z2 ); 
    }
    else if( tc == CSG_UNION || tc == CSG_INTERSECTION || tc == CSG_DIFFERENCE )
    {
        setAABB( 0.f );  
    }
    else
    {
        setAABB( UNBOUNDED_DEFAULT_EXTENT ); 
    }
}

CSGNode CSGNode::Sphere(float radius)  // static
{
    assert( radius > 0.f); 
    CSGNode nd = {} ;
    nd.setParam( 0.f, 0.f, 0.f, radius,  0.f,  0.f ); 
    nd.setAABB(  -radius, -radius, -radius,  radius, radius, radius  ); 
    nd.setTypecode(CSG_SPHERE) ; 
    return nd ;
}
CSGNode CSGNode::ZSphere(float radius, float z1, float z2)  // static
{
    assert( radius > 0.f); 
    assert( z2 > z1 ); 
    CSGNode nd = {} ;
    nd.setParam( 0.f, 0.f, 0.f, radius, z1, z2 ); 
    nd.setAABB(  -radius, -radius, z1,  radius, radius, z2  ); 
    nd.setTypecode(CSG_ZSPHERE) ; 
    return nd ;
}
CSGNode CSGNode::Cone(float r1, float z1, float r2, float z2)  // static
{
    assert( z2 > z1 ); 
    float rmax = fmaxf(r1, r2) ;
    CSGNode nd = {} ;
    nd.setParam( r1, z1, r2, z2, 0.f, 0.f ) ;
    nd.setAABB( -rmax, -rmax, z1, rmax, rmax, z2 ); 
    nd.setTypecode(CSG_CONE) ; 
    return nd ; 
}
CSGNode CSGNode::Box3(float fullside)  // static 
{
    return Box3(fullside, fullside, fullside); 
}
CSGNode CSGNode::Box3(float fx, float fy, float fz )  // static 
{
    assert( fx > 0.f ); 
    assert( fy > 0.f ); 
    assert( fz > 0.f ); 

    CSGNode nd = {} ;
    nd.setParam( fx, fy, fz, 0.f, 0.f, 0.f ); 
    nd.setAABB( -fx*0.5f , -fy*0.5f, -fz*0.5f, fx*0.5f , fy*0.5f, fz*0.5f );   
    nd.setTypecode(CSG_BOX3) ; 
    return nd ; 
}
CSGNode CSGNode::Cylinder(float px, float py, float radius, float z1, float z2)
{
    CSGNode nd = {} ; 
    nd.setParam( px, py, 0.f, radius, z1, z2)  ; 
    nd.setAABB( px-radius, py-radius, z1, px+radius, py+radius, z2 );   
    nd.setTypecode(CSG_CYLINDER); 
    return nd ; 
} 
CSGNode CSGNode::Disc(float px, float py, float ir, float r, float z1, float z2)
{
    CSGNode nd = {} ;
    nd.setParam( px, py, ir, r, z1, z2 ); 
    nd.setAABB( px - r , py - r , z1, px + r, py + r, z2 ); 
    nd.setTypecode(CSG_DISC); 
    return nd ; 
} 

CSGNode CSGNode::Hyperboloid(float r0, float zf, float z1, float z2) // static
{
    assert( z1 < z2 ); 
    const float rr0 = r0*r0 ; 
    const float z1s = z1/zf ; 
    const float z2s = z2/zf ; 

    const float rr1 = rr0 * ( z1s*z1s + 1.f ) ;
    const float rr2 = rr0 * ( z2s*z2s + 1.f ) ;
    const float rmx = sqrtf(fmaxf( rr1, rr2 )) ; 

    CSGNode nd = {} ;
    nd.setParam(r0, zf, z1, z2, 0.f, 0.f ) ; 
    nd.setAABB(  -rmx,  -rmx,  z1,  rmx, rmx, z2 ); 
    nd.setTypecode(CSG_HYPERBOLOID) ; 
    return nd ; 
}


CSGNode CSGNode::Plane(float nx, float ny, float nz, float d)
{
    CSGNode nd = {} ;
    nd.setParam(nx, ny, nz, d, 0.f, 0.f ) ;
    nd.setTypecode(CSG_PLANE) ; 
    nd.setAABB( UNBOUNDED_DEFAULT_EXTENT ); 
    return nd ; 
}

CSGNode CSGNode::Slab(float nx, float ny, float nz, float d1, float d2 )
{
    CSGNode nd = {} ;
    nd.setParam( nx, ny, nz, 0.f, d1, d2 ); 
    nd.setTypecode(CSG_SLAB) ; 
    nd.setAABB( UNBOUNDED_DEFAULT_EXTENT ); 
    return nd ; 
}



/**
CSGNode::Make
------------

Only the first four chars of the name are used to select the type of node.

**/

CSGNode CSGNode::Make(const char* name) // static
{
    if(strncmp(name, "sphe", 4) == 0) return CSGNode::Sphere(100.f) ; 
    if(strncmp(name, "zsph", 4) == 0) return CSGNode::ZSphere(100.f, -50.f, 50.f) ; 
    if(strncmp(name, "cone", 4) == 0) return CSGNode::Cone(150.f, -150.f, 50.f, -50.f) ; 
    if(strncmp(name, "hype", 4) == 0) return CSGNode::Hyperboloid(100.f, 50.f, -50.f, 50.f) ; 
    if(strncmp(name, "box3", 4) == 0) return CSGNode::Box3(50.f, 100.f, 150.f) ; 
    if(strncmp(name, "plan", 4) == 0) return CSGNode::Plane(1.f, 0.f, 0.f, 0.f) ; 
    if(strncmp(name, "slab", 4) == 0) return CSGNode::Slab(1.f, 0.f, 0.f, -10.f, 10.f ) ; 
    if(strncmp(name, "cyli", 4) == 0) return CSGNode::Cylinder(0.f, 0.f, 100.f, -50.f, 50.f ) ; 
    if(strncmp(name, "disc", 4) == 0) return CSGNode::Disc(    0.f, 0.f, 50.f, 100.f, -2.f, 2.f ) ; 
    assert(0); 
    return CSGNode::Sphere(1.0); 
}

CSGNode CSGNode::Make(unsigned typecode, const float* param6, const float* aabb ) // static
{
    CSGNode nd = {} ;
    nd.setTypecode(typecode) ; 
    if(param6) nd.setParam( param6 );  
    if(aabb)    nd.setAABB( aabb );  
    return nd ; 
}



#endif

