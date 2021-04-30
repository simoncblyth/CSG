#include <iostream>
#include <iomanip>
#include <array>

#include <glm/glm.hpp>
#include <glm/gtx/string_cast.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "sutil_vec_math.h"
#include "OpticksCSG.h"
#include "CSGSolid.h"
#include "CU.h"
#include "NP.hh"
#include "CSGFoundry.h"
#include "AABB.h"


CSGFoundry::CSGFoundry()
    :
    imax(1000),
    d_solid(nullptr),
    d_prim(nullptr),
    d_node(nullptr),
    d_plan(nullptr),
    d_tran(nullptr),
    d_itra(nullptr)
{
    init(); 
}

void CSGFoundry::init()
{
    // without sufficient reserved the vectors may reallocate on any push_back invalidating prior pointers 
    solid.reserve(imax); 
    prim.reserve(imax); 
    node.reserve(imax); 
    plan.reserve(imax); 
    tran.reserve(imax); 
    itra.reserve(imax); 
}

void CSGFoundry::makeDemoSolids()
{
    makeSphere(); 
    makeZSphere(); 
    makeCone(); 
    makeHyperboloid(); 
    makeBox3(); 
    makePlane(); 
    makeSlab(); 
    makeCylinder() ; 
    makeDisc(); 
    makeConvexPolyhedronCube(); 
    makeConvexPolyhedronTetrahedron(); 
    makeEllipsoid(); 
    makeUnionBoxSphere();
    makeIntersectionBoxSphere();
    makeDifferenceBoxSphere();
    makeRotatedCylinder();
    dump(); 
}

std::string CSGFoundry::getBashMap()
{
    unsigned num_solids = getNumSolid(); 
    std::stringstream ss ; 
    ss << "#CSGFoundry::getBashMap " << std::endl ; 
    ss << "declare -A map " << std::endl ;  
    for(unsigned i=0 ; i < num_solids ; i++)
    {
        const CSGSolid* so = getSolid(i); 
        ss << "map[" << i << "]=" << so->label << std::endl ; 
    }
    return ss.str(); 
}


void CSGFoundry::dump() const 
{
    for(unsigned idx=0 ; idx < solid.size() ; idx++) dumpSolid(idx); 
}

void CSGFoundry::dumpSolid(unsigned solidIdx) const 
{
    std::cout << "CSGFoundry::dumpSolid " << solidIdx << std::endl ; 

    const CSGSolid* so = solid.data() + solidIdx ; 
    std::cout << so->desc() << std::endl ; 

    for(unsigned primIdx=so->primOffset ; primIdx < so->primOffset+so->numPrim ; primIdx++)
    {
        const CSGPrim* pr = prim.data() + primIdx ; 
        std::cout 
            << " primIdx " << std::setw(3) << primIdx << " "
            << pr->desc() 
            << std::endl 
            ; 

        for(unsigned nodeIdx=pr->nodeOffset() ; nodeIdx < pr->nodeOffset()+pr->numNode() ; nodeIdx++)
        {
            const CSGNode* nd = node.data() + nodeIdx ; 
            std::cout << nd->desc() << std::endl ; 
        }
    } 
}

void CSGFoundry::dumpPrim(unsigned solidIdx) const 
{
    std::cout << "CSGFoundry::dumpPrim" << solidIdx << std::endl ; 
    const CSGSolid* so = solid.data() + solidIdx ; 
    std::cout << so->desc() << std::endl ; 
    for(unsigned primIdx=so->primOffset ; primIdx < so->primOffset+so->numPrim ; primIdx++)
    {
        const CSGPrim* pr = prim.data() + primIdx ; 
        std::cout 
            << " primIdx " << std::setw(3) << primIdx << " "
            << pr->desc() 
            << std::endl 
            ; 
    } 
}

void CSGFoundry::dumpNode(unsigned solidIdx) const 
{
    std::cout << "CSGFoundry::dumpNode" << solidIdx << std::endl ; 
    const CSGSolid* so = solid.data() + solidIdx ; 
    std::cout << so->desc() << std::endl ; 
    for(unsigned primIdx=so->primOffset ; primIdx < so->primOffset+so->numPrim ; primIdx++)
    {
        const CSGPrim* pr = prim.data() + primIdx ; 
        for(unsigned nodeIdx=pr->nodeOffset() ; nodeIdx < pr->nodeOffset()+pr->numNode() ; nodeIdx++)
        {
            const CSGNode* nd = node.data() + nodeIdx ; 
            std::cout << nd->desc() << std::endl ; 
        }
    } 
}




/**
CSGFoundry::getPrimSpec
----------------------

Provides the specification to access the AABB and sbtIndexOffset of all CSGPrim 
of a CSGSolid.  The specification includes pointers, counts and stride.

NB PrimAABB is distinct from NodeAABB. Cannot directly use NodeAABB 
because the number of nodes for each prim (node tree) varies meaning 
that the strides are irregular. 

**/

CSGPrimSpec CSGFoundry::getPrimSpec(unsigned solidIdx) const 
{
    CSGPrimSpec ps = d_prim ? getPrimSpecDevice(solidIdx) : getPrimSpecHost(solidIdx) ; 
    if(ps.device == false) std::cout << "CSGFoundry::getPrimSpec WARNING using host PrimSpec " << std::endl ; 
    return ps ; 
}
CSGPrimSpec CSGFoundry::getPrimSpecHost(unsigned solidIdx) const 
{
    const CSGSolid* so = solid.data() + solidIdx ; 
    CSGPrimSpec ps = CSGPrim::MakeSpec( prim.data(),  so->primOffset, so->numPrim ); ; 
    ps.device = false ; 
    return ps ; 
}
CSGPrimSpec CSGFoundry::getPrimSpecDevice(unsigned solidIdx) const 
{
    assert( d_prim ); 
    const CSGSolid* so = solid.data() + solidIdx ; 
    CSGPrimSpec ps = CSGPrim::MakeSpec( d_prim,  so->primOffset, so->numPrim ); ; 
    ps.device = true ; 
    return ps ; 
}


unsigned CSGFoundry::getNumSolid() const { return solid.size(); } 
unsigned CSGFoundry::getNumPrim() const  { return prim.size();  } 
unsigned CSGFoundry::getNumNode() const  { return node.size(); }
unsigned CSGFoundry::getNumPlan() const  { return plan.size(); }
unsigned CSGFoundry::getNumTran() const  { return tran.size(); }
unsigned CSGFoundry::getNumItra() const  { return itra.size(); }
unsigned CSGFoundry::getNumInst() const  { return inst.size(); }

const CSGSolid*  CSGFoundry::getSolid(unsigned solidIdx) const { return solidIdx < solid.size() ? solid.data() + solidIdx  : nullptr ; } 
const CSGPrim*   CSGFoundry::getPrim(unsigned primIdx)   const { return primIdx  < prim.size()  ? prim.data()  + primIdx  : nullptr ; } 
const CSGNode*   CSGFoundry::getNode(unsigned nodeIdx)   const { return nodeIdx  < node.size()  ? node.data()  + nodeIdx  : nullptr ; }  
const float4*    CSGFoundry::getPlan(unsigned planIdx)   const { return planIdx  < plan.size()  ? plan.data()  + planIdx  : nullptr ; }
const qat4*      CSGFoundry::getTran(unsigned tranIdx)   const { return tranIdx  < tran.size()  ? tran.data()  + tranIdx  : nullptr ; }
const qat4*      CSGFoundry::getItra(unsigned itraIdx)   const { return itraIdx  < itra.size()  ? itra.data()  + itraIdx  : nullptr ; }
const qat4*      CSGFoundry::getInst(unsigned instIdx)   const { return instIdx  < inst.size()  ? inst.data()  + instIdx  : nullptr ; }


const CSGSolid*  CSGFoundry::getSolid_(int solidIdx_) const { 
    unsigned solidIdx = solidIdx_ < 0 ? unsigned(solid.size() + solidIdx_) : unsigned(solidIdx_)  ;   // -ve counts from end
    return getSolid(solidIdx); 
}   

const CSGSolid* CSGFoundry::getSolidByName(const char* name) const  // caution stored labels truncated to 4 char 
{
    unsigned missing = ~0u ; 
    unsigned idx = missing ; 
    for(unsigned i=0 ; i < solid.size() ; i++) if(strcmp(solid[i].label, name) == 0) idx = i ;  
    assert( idx != missing ); 
    return getSolid(idx) ; 
}

/**
CSGFoundry::getSolidIdx
----------------------

Without sufficient reserve allocation this is unreliable as pointers go stale on reallocations.

**/

unsigned CSGFoundry::getSolidIdx(const CSGSolid* so) const 
{
    unsigned idx = ~0u ; 
    for(unsigned i=0 ; i < solid.size() ; i++) 
    {
       const CSGSolid* s = solid.data() + i ; 
       std::cout << " i " << i << " s " << s << " so " << so << std::endl ; 
       if(s == so) idx = i ;  
    } 
    assert( idx != ~0u ); 
    return idx ; 
}




CSGSolid* CSGFoundry::make(const char* name)
{
    CSGSolid* so = nullptr ; 
    if(     strcmp(name, "sphe") == 0) so = makeSphere(name) ;
    else if(strcmp(name, "zsph") == 0) so = makeZSphere(name) ;
    else if(strcmp(name, "cone") == 0) so = makeCone(name) ;
    else if(strcmp(name, "hype") == 0) so = makeHyperboloid(name) ;
    else if(strcmp(name, "box3") == 0) so = makeBox3(name) ;
    else if(strcmp(name, "plan") == 0) so = makePlane(name) ;
    else if(strcmp(name, "slab") == 0) so = makeSlab(name) ;
    else if(strcmp(name, "cyli") == 0) so = makeCylinder(name) ;
    else if(strcmp(name, "disc") == 0) so = makeDisc(name) ;
    else if(strcmp(name, "vcub") == 0) so = makeConvexPolyhedronCube(name) ;
    else if(strcmp(name, "vtet") == 0) so = makeConvexPolyhedronTetrahedron(name) ;
    else if(strcmp(name, "elli") == 0) so = makeEllipsoid(name) ;
    else if(strcmp(name, "ubsp") == 0) so = makeUnionBoxSphere(name) ;
    else if(strcmp(name, "ibsp") == 0) so = makeIntersectionBoxSphere(name) ;
    else if(strcmp(name, "dbsp") == 0) so = makeDifferenceBoxSphere(name) ;
    else if(strcmp(name, "rcyl") == 0) so = makeRotatedCylinder(name) ;
    assert( so ); 
    return so ;  
}

/**
CSGFoundry::addNode
--------------------

**/

CSGNode* CSGFoundry::addNode(CSGNode nd, const std::vector<float4>* pl )
{
    unsigned num_planes = pl ? pl->size() : 0 ; 
    if(num_planes > 0)
    {
        nd.setTypecode(CSG_CONVEXPOLYHEDRON) ; 
        nd.setPlaneIdx(plan.size());    
        nd.setPlaneNum(num_planes);    
        for(unsigned i=0 ; i < num_planes ; i++) addPlan((*pl)[i]);  
    }

    unsigned idx = node.size() ;  
    assert( idx < imax ); 
    node.push_back(nd); 
    return node.data() + idx ; 
}

float4* CSGFoundry::addPlan(const float4& pl )
{
    unsigned idx = plan.size(); 
    assert( idx < imax ); 
    plan.push_back(pl); 
    return plan.data() + idx ; 
}

template<typename T>
unsigned CSGFoundry::addTran( const Tran<T>& tr  )
{
    qat4 t(glm::value_ptr(tr.t));  // narrowing when T=double
    qat4 v(glm::value_ptr(tr.v)); 

    unsigned idx = tran.size(); 

    /*
    std::cout << "CSGFoundry::addTran idx " << idx << std::endl ; 
    std::cout << " t " << t << std::endl ; 
    std::cout << " v " << v << std::endl ; 
    */

    assert( tran.size() == itra.size()) ; 
    tran.push_back(t); 
    itra.push_back(v); 

    return idx ;  
}


template unsigned CSGFoundry::addTran<float>(const Tran<float>& ) ;
template unsigned CSGFoundry::addTran<double>(const Tran<double>& ) ;


CSGNode* CSGFoundry::addNodes(const std::vector<CSGNode>& nds )
{
    unsigned idx = node.size() ; 
    for(unsigned i=0 ; i < nds.size() ; i++) 
    {
        const CSGNode& nd = nds[i]; 
        idx = node.size() ;  
        assert( idx < imax ); 
        node.push_back(nd); 
    }
    return node.data() + idx ; 
}

/**
CSGFoundry::addPrim
------------------

Offsets counts for  node, tran and plan are 
persisted into the CSGPrim. 
Thus must addPrim prior to adding any node, 
tran or plan needed for a prim.

**/

CSGPrim* CSGFoundry::addPrim(int num_node)  
{
    CSGPrim pr = {} ;
    pr.setNumNode(num_node) ; 
    pr.setNodeOffset(node.size()); 
    pr.setTranOffset(tran.size()); 
    pr.setPlanOffset(plan.size()); 
    pr.setSbtIndexOffset(0) ; 

    unsigned primIdx = prim.size(); 
    assert( primIdx < imax ); 
    prim.push_back(pr); 
    return prim.data() + primIdx ; 
}


CSGSolid* CSGFoundry::addSolid(unsigned num_prim, const char* label )
{
    unsigned idx = solid.size(); 
    assert( idx < imax ); 

    unsigned primOffset = prim.size(); 
    CSGSolid so = {} ; 
    memcpy( so.label, label, 4 ); 
    so.numPrim = num_prim ; 
    so.primOffset = primOffset ; 
    so.center_extent = make_float4(0.f, 0.f, 0.f, 100.f) ;  // should be overwritten 

    solid.push_back(so); 
    return solid.data() + idx  ; 
}



/**
Foundary::makeLayered
----------------------------

Once have transforms working can generalize to any shape. 
But prior to that just do layering for sphere for adiabatic transition
from Shape to CSGFoundry/CSGSolid.

NB Each layer is a separate CSGPrim with a single CSGNode 

**/

CSGSolid* CSGFoundry::makeLayered(const char* label, float outer_radius, unsigned layers )
{
    std::vector<float> radii ;
    for(unsigned i=0 ; i < layers ; i++) radii.push_back(outer_radius*float(layers-i)/float(layers)) ; 

    unsigned numPrim = layers ; 
    CSGSolid* so = addSolid(numPrim, label); 
    so->center_extent = make_float4( 0.f, 0.f, 0.f, outer_radius ) ; 

    for(unsigned i=0 ; i < numPrim ; i++)
    {
        unsigned numNode = 1 ; 
        CSGPrim* p = addPrim(numNode); 
        float radius = radii[i]; 

        CSGNode* n = nullptr ; 

        if(strcmp(label, "sphere") == 0)
        {
            n = addNode(CSGNode::Sphere(radius)); 
        }
        else if(strcmp(label, "zsphere") == 0)
        {
            n = addNode(CSGNode::ZSphere(radius, -radius/2.f , radius/2.f )); 
        }
        else
        {
            assert( 0 && "layered only implemented for sphere and zsphere currently" ); 
        } 

        p->setSbtIndexOffset(i) ; 
        p->setAABB( n->AABB() ); 
    }
    return so ; 
}

CSGSolid* CSGFoundry::makeClustered(const char* name,  int i0, int i1, int is, int j0, int j1, int js, int k0, int k1, int ks, double unit, bool inbox ) 
{

    /*
    std::array<int, 9> cluster = {i0,i1,is,j0,j1,js,k0,k1,ks} ; 
    int3 imn = make_int3(i0,j0,k0); 
    int3 imx = make_int3(i1,j1,k1);
    float3 mn = make_float3( float(imn.x), float(imn.y), float(imn.z))*float(unit) ;  
    float3 mx = make_float3( float(imx.x), float(imx.y), float(imx.z))*float(unit) ;  
    AABB cbb = {mn, mx } ; 
    float4 cluster_ce = cbb.center_extent(); 
    */
    
    unsigned numPrim = inbox ? 1 : 0 ; 
    for(int i=i0 ; i < i1 ; i+=is ) 
    for(int j=j0 ; j < j1 ; j+=js ) 
    for(int k=k0 ; k < k1 ; k+=ks ) 
    {
        //std::cout << std::setw(2) << numPrim << " (i,j,k) " << "(" << i << "," << j << "," << k << ") " << std::endl ; 
        numPrim += 1 ; 
    }


       
    std::cout 
        << "CSGFoundry::makeClustered " 
        << " name " << name  
        << " numPrim " << numPrim 
        << " inbox " << inbox
        << std::endl
        ;  

    CSGSolid* so = addSolid(numPrim, name);
    unsigned idx = 0 ; 

    AABB bb = {} ; 
 
    for(int i=i0 ; i < i1 ; i+=is ) 
    for(int j=j0 ; j < j1 ; j+=js ) 
    for(int k=k0 ; k < k1 ; k+=ks ) 
    {
        unsigned numNode = 1 ; 
        CSGPrim* p = addPrim(numNode); 
        CSGNode* n = addNode(CSGNode::Make(name)) ;
    
        const Tran<double>* translate = Tran<double>::make_translate( double(i)*unit, double(j)*unit, double(k)*unit ); 
        unsigned transform_idx = 1 + addTran(*translate);      // 1-based idx, 0 meaning None
        n->setTransform(transform_idx); 
        const qat4* t = getTran(transform_idx-1u) ; 

        t->transform_aabb_inplace( n->AABB() ); 

        bb.include_aabb( n->AABB() ); 

        p->setSbtIndexOffset(idx) ; 
        p->setAABB( n->AABB() );

        DumpAABB("p->AABB() aft setup", p->AABB() ); 
        
        std::cout << " idx " << idx << " transform_idx " << transform_idx << std::endl ; 
 
        idx += 1 ; 
    }


    if(inbox)
    {
        float4 ce = bb.center_extent(); 
        float fullside = ce.w*2.f ; 

        const Tran<float>* to_center = Tran<float>::make_translate( float(ce.x), float(ce.y), float(ce.z) ); 
        unsigned transform_idx = 1 + addTran(*to_center);  // 1-based idx, 0 meaning None
        const qat4* t = getTran(transform_idx-1u) ; 

        CSGPrim* p = addPrim(1); 
        CSGNode bx = CSGNode::Box3(fullside) ;
        CSGNode* n = addNode(CSGNode::Box3(fullside)); 
        n->setTransform(transform_idx); 
        t->transform_aabb_inplace( n->AABB() ); 

        p->setSbtIndexOffset(idx); 
        p->setAABB( n->AABB() );
        
        idx += 1 ; 
    }


    so->center_extent = bb.center_extent() ;   // contains AABB of all CSGPrim 
    std::cout << " so->center_extent " << so->center_extent << std::endl ; 
    return so ; 
}

void CSGFoundry::DumpAABB(const char* msg, const float* aabb) // static 
{
    int w = 4 ; 
    std::cout << msg << " " ; 
    std::cout << " | " ; 
    for(int l=0 ; l < 3 ; l++) std::cout << std::setw(w) << *(aabb+l) << " " ; 
    std::cout << " | " ; 
    for(int l=0 ; l < 3 ; l++) std::cout << std::setw(w) << *(aabb+l+3) << " " ; 
    std::cout << " | " ; 
    for(int l=0 ; l < 3 ; l++) std::cout << std::setw(w) << *(aabb+l+3) - *(aabb+l)  << " " ; 
    std::cout << std::endl ; 
}



/**
CSGFoundry::makeSolid11 makes 1-CSGPrim with 1-Node
------------------------------------------------
**/

CSGSolid* CSGFoundry::makeSolid11(const char* label, CSGNode nd, const std::vector<float4>* pl  ) 
{
    unsigned numPrim = 1 ; 
    CSGSolid* so = addSolid(numPrim, label);

    unsigned numNode = 1 ; 
    CSGPrim* p = addPrim(numNode); 
    CSGNode* n = addNode(nd, pl ); 
    p->setAABB( n->AABB() ); 

    float extent = p->extent(); 
    if(extent == 0.f )
        std::cout << "CSGFoundry::makeSolid11 : FATAL : " << label << " : got zero extent " << std::endl ; 
    assert( extent > 0.f ); 

    AABB bb = AABB::Make( p->AABB() ); 
    so->center_extent = bb.center_extent()  ; 
    std::cout << "CSGFoundry::makeSolid11 so.label " << so->label << " so.center_extent " << so->center_extent << std::endl ; 
    return so ; 
}

CSGSolid* CSGFoundry::makeBooleanBoxSphere( const char* label, char op_, float radius, float fullside )
{
    CSGNode op = CSGNode::BooleanOperator(op_); 
    CSGNode bx = CSGNode::Box3(fullside) ; 
    CSGNode sp = CSGNode::Sphere(radius); 

    unsigned numPrim = 1 ; 
    CSGSolid* so = addSolid(numPrim, label);

    unsigned numNode = 3 ; 
    CSGPrim* p = addPrim(numNode); 

    addNode(op); 
    addNode(bx); 
    addNode(sp); 
     
    AABB bb = {} ;
    bb.include_aabb( bx.AABB() );     // assumes any transforms have been applied to the Node AABB
    bb.include_aabb( sp.AABB() ); 
    p->setAABB( bb.data() );  

    so->center_extent = bb.center_extent()  ; 
    std::cout << "CSGFoundry::makeUnionBoxSphere so.label " << so->label << " so.center_extent " << so->center_extent << std::endl ; 
    return so ; 
}

CSGSolid* CSGFoundry::makeUnionBoxSphere( const char* label, float radius, float fullside ){
    return makeBooleanBoxSphere(label, 'U', radius, fullside ); 
}
CSGSolid* CSGFoundry::makeIntersectionBoxSphere( const char* label, float radius, float fullside ){
    return makeBooleanBoxSphere(label, 'I', radius, fullside ); 
}
CSGSolid* CSGFoundry::makeDifferenceBoxSphere( const char* label, float radius, float fullside ){
    return makeBooleanBoxSphere(label, 'D', radius, fullside ); 
}


CSGSolid* CSGFoundry::makeSphere(const char* label, float radius)
{
    CSGNode nd = CSGNode::Sphere(radius); 
    return makeSolid11(label, nd ); 
}
CSGSolid* CSGFoundry::makeEllipsoid(  const char* label, float rx, float ry, float rz )
{
    CSGNode nd = CSGNode::Sphere(rx);

    double dx = double(rx) ; 
    double dy = double(ry) ; 
    double dz = double(rz) ; 
 
    double sx = double(1.) ; 
    double sy = dy/dx ; 
    double sz = dz/dx ; 

    const Tran<double>* tr = Tran<double>::make_scale(sx, sy, sz ); 

    unsigned idx = 1 + addTran(*tr);      // 1-based idx, 0 meaning None
    std::cout << "CSGFoundry::makeEllipsoid " << *tr << std::endl ;

    nd.setTransform(idx); 
    return makeSolid11(label, nd ); 
}


CSGSolid* CSGFoundry::makeRotatedCylinder(const char* label, float px, float py, float radius, float z1, float z2, float ax, float ay, float az, float angle_deg )
{
    CSGNode nd = CSGNode::Cylinder( px, py, radius, z1, z2 ); 
    const Tran<float>* tr = Tran<float>::make_rotate(ax, ay, az, angle_deg ); 
    unsigned idx = 1 + addTran(*tr);      // 1-based idx, 0 meaning None
    std::cout << "CSGFoundry::makeRotatedCylinder " << *tr << std::endl ;
    nd.setTransform(idx); 
    return makeSolid11(label, nd ); 
}







CSGSolid* CSGFoundry::makeZSphere(const char* label, float radius, float z1, float z2)
{
    CSGNode nd = CSGNode::ZSphere(radius, z1, z2); 
    return makeSolid11(label, nd ); 
}

CSGSolid* CSGFoundry::makeCone(const char* label, float r1, float z1, float r2, float z2)
{
    CSGNode nd = CSGNode::Cone(r1, z1, r2, z2 ); 
    return makeSolid11(label, nd ); 
}

CSGSolid* CSGFoundry::makeHyperboloid(const char* label, float r0, float zf, float z1, float z2)
{
    CSGNode nd = CSGNode::Hyperboloid( r0, zf, z1, z2 ); 
    return makeSolid11(label, nd ); 
}

CSGSolid* CSGFoundry::makeBox3(const char* label, float fx, float fy, float fz )
{
    CSGNode nd = CSGNode::Box3(fx, fy, fz); 
    return makeSolid11(label, nd ); 
}

CSGSolid* CSGFoundry::makePlane(const char* label, float nx, float ny, float nz, float d)
{
    CSGNode nd = CSGNode::Plane(nx, ny, nz, d ); 
    return makeSolid11(label, nd ); 
}

CSGSolid* CSGFoundry::makeSlab(const char* label, float nx, float ny, float nz, float d1, float d2 )
{
    CSGNode nd = CSGNode::Slab( nx, ny, nz, d1, d1 ); 
    return makeSolid11(label, nd ); 
}

CSGSolid* CSGFoundry::makeCylinder(const char* label, float px, float py, float radius, float z1, float z2)
{
    CSGNode nd = CSGNode::Cylinder( px, py, radius, z1, z2 ); 
    return makeSolid11(label, nd ); 
}


CSGSolid* CSGFoundry::makeDisc(const char* label, float px, float py, float ir, float r, float z1, float z2)
{
    CSGNode nd = CSGNode::Disc(px, py, ir, r, z1, z2 ); 
    return makeSolid11(label, nd ); 
}


float4 CSGFoundry::TriPlane( const std::vector<float3>& v, unsigned i, unsigned j, unsigned k )  // static 
{
    // normal for plane through v[i] v[j] v[k]
    float3 ij = v[j] - v[i] ; 
    float3 ik = v[k] - v[i] ; 
    float3 n = normalize(cross(ij, ik )) ;
    float di = dot( n, v[i] ) ;
    float dj = dot( n, v[j] ) ;
    float dk = dot( n, v[k] ) ;
    //std::cout << " di " << di << " dj " << dj << " dk " << dk << " n (" << n.x << "," << n.y << "," << n.z << ")" << std::endl ; 
    float4 plane = make_float4( n, di ) ; 
    return plane ;  
}

CSGSolid* CSGFoundry::makeConvexPolyhedronCube(const char* label, float extent)
{
    float hx = extent ; 
    float hy = extent/2.f ; 
    float hz = extent/3.f ; 

    std::vector<float4> pl ; 
    pl.push_back( make_float4(  1.f,  0.f,  0.f, hx ) ); 
    pl.push_back( make_float4( -1.f,  0.f,  0.f, hx ) ); 
    pl.push_back( make_float4(  0.f,  1.f,  0.f, hy ) ); 
    pl.push_back( make_float4(  0.f, -1.f,  0.f, hy ) ); 
    pl.push_back( make_float4(  0.f,  0.f,  1.f, hz ) ); 
    pl.push_back( make_float4(  0.f,  0.f, -1.f, hz ) );

    CSGNode nd = {} ;
    nd.setAABB(-hx, -hy, -hz, hx, hy, hz); 
    return makeSolid11(label, nd, &pl ); 
}


/*  
     https://en.wikipedia.org/wiki/Tetrahedron

       0:(1,1,1)
       1:(1,−1,−1)
       2:(−1,1,−1) 
       3:(−1,−1,1)

                              (1,1,1)
                 +-----------0
                /|          /| 
     (-1,-1,1) / |         / |
              3-----------+  |
              |  |        |  |
              |  |        |  |
   (-1,1,-1)..|..2--------|--+
              | /         | /
              |/          |/
              +-----------1
                          (1,-1,-1)      

          Faces (right-hand-rule oriented outwards normals)
                0-1-2
                1-3-2
                3-0-2
                0-3-1

         z  y
         | /
         |/
         +---> x
*/

CSGSolid* CSGFoundry::makeConvexPolyhedronTetrahedron(const char* label, float extent)
{
    //extent = 100.f*sqrt(3); 
    float s = extent ; 

    std::vector<float3> vtx ; 
    vtx.push_back(make_float3( s, s, s));  
    vtx.push_back(make_float3( s,-s,-s)); 
    vtx.push_back(make_float3(-s, s,-s)); 
    vtx.push_back(make_float3(-s,-s, s)); 

    std::vector<float4> pl ; 
    pl.push_back(TriPlane(vtx, 0, 1, 2)) ;  
    pl.push_back(TriPlane(vtx, 1, 3, 2)) ;  
    pl.push_back(TriPlane(vtx, 3, 0, 2)) ;  
    pl.push_back(TriPlane(vtx, 0, 3, 1)) ;  

    //for(unsigned i=0 ; i < pl.size() ; i++) std::cout << " pl (" << pl[i].x << "," << pl[i].y << "," << pl[i].z << "," << pl[i].w << ") " << std::endl ;
 
    CSGNode nd = {} ;
    nd.setAABB(extent); 
    return makeSolid11(label, nd, &pl ); 
}








void CSGFoundry::write(const char* base, const char* rel) const 
{
    std::stringstream ss ;   
    ss << base << "/" << rel ; 
    std::string dir = ss.str();   

    std::cout << "CSGFoundry::write " << dir << std::endl ; 

    NP::Write(dir.c_str(), "solid.npy",  (int*)solid.data(),  solid.size(),  4, 2 ); 
    NP::Write(dir.c_str(), "prim.npy",   (float*)prim.data(), prim.size(),   4, 4 ); 
    NP::Write(dir.c_str(), "node.npy",   (float*)node.data(), node.size(),   4, 4 ); 
    NP::Write(dir.c_str(), "plan.npy",   (float*)plan.data(), plan.size(),   4    ); 
    NP::Write(dir.c_str(), "tran.npy",   (float*)tran.data(), tran.size(),   4, 4 ); 
    NP::Write(dir.c_str(), "itra.npy",   (float*)itra.data(), itra.size(),   4, 4 ); 

    unsigned num_inst = inst.size(); 
    if( num_inst > 0 )
    NP::Write(dir.c_str(), "inst.npy",   (float*)inst.data(), inst.size(),   4, 4 ); 
}

void CSGFoundry::upload()
{
    std::cout 
        << "CSGFoundry::upload"
        << " num_solid " << solid.size()
        << " num_prim " << prim.size()
        << " num_node " << node.size()
        << " num_plan " << plan.size()
        << " num_tran " << tran.size()
        << " num_itra " << itra.size()
        << " num_inst " << inst.size()
        << std::endl
        ;

    assert( tran.size() == itra.size() ); 

    d_solid = solid.size() > 0 ? CU::UploadArray<CSGSolid>(solid.data(), solid.size() ) : nullptr ; 
    d_prim = prim.size() > 0 ? CU::UploadArray<CSGPrim>(prim.data(), prim.size() ) : nullptr ; 
    d_node = node.size() > 0 ? CU::UploadArray<CSGNode>(node.data(), node.size() ) : nullptr ; 
    d_plan = plan.size() > 0 ? CU::UploadArray<float4>(plan.data(), plan.size() ) : nullptr ; 
    d_tran = tran.size() > 0 ? CU::UploadArray<qat4>(tran.data(), tran.size() ) : nullptr ; 
    d_itra = itra.size() > 0 ? CU::UploadArray<qat4>(itra.data(), itra.size() ) : nullptr ; 

    // note that d_solid and d_tran are not actually used on GPU currently 

    inst_find_unique(); 
}

void CSGFoundry::inst_find_unique()
{
    qat4::find_unique( inst, ins, gas, ias ); 
    std::cout 
        << "CSGFoundry::inst_find_unique"
        << " num_inst " << inst.size()
        << " ins " << ins.size()
        << " gas " << gas.size()
        << " ias " << ias.size()
        << std::endl
        ;
}

unsigned CSGFoundry::getNumUniqueIAS() const
{
    return ias.size(); 
}
unsigned CSGFoundry::getNumUniqueGAS() const
{
    return gas.size(); 
}
unsigned CSGFoundry::getNumUniqueINS() const
{
    return ins.size(); 
}

unsigned CSGFoundry::getNumInstancesIAS(unsigned ias_idx) const
{
    return qat4::count_ias(inst, ias_idx );  
}

void CSGFoundry::getInstanceTransformsIAS(std::vector<qat4>& ias_inst, unsigned ias_idx )
{
    qat4::select_instances_ias(inst, ias_inst, ias_idx ) ;
}




