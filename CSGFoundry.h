#pragma once

#include <vector>
#include <glm/glm.hpp>

#include "CSGSolid.h"
#include "CSGPrim.h"
#include "CSGNode.h"
#include "Quad.h"
#include "qat4.h"
#include "Tran.h"

/**
CSGFoundry
============

* CSGSolids contain one or more CSGPrim  (CSGPrim would correspond to Geant4 G4VSolid)
* CSGPrim contain one or more CSGNode    (CSGNode are CSG constituent nodes) 

**/


struct CSGFoundry
{
    static const unsigned IMAX ; 
    static CSGFoundry* Load(const char* base, const char* rel);
    static CSGFoundry* Load(const char* dir );
    static int Compare(const CSGFoundry* a , const CSGFoundry* b ); 

    template<typename T>
    static int CompareVec( const char* name, const std::vector<T>& a, const std::vector<T>& b );

    static int CompareBytes(const void* a, const void* b, unsigned num_bytes);


    CSGFoundry();
    void init(); 

    void makeDemoSolids() ;
    std::string desc() const ;
    void summary(const char* msg="CSGFoundry::summary") const ;


    void dump() const ;
    void dumpSolid(unsigned solidIdx ) const ;
    void dumpPrim(unsigned solidIdx ) const ;
    void dumpNode(unsigned solidIdx ) const ;

    CSGPrimSpec getPrimSpec(       unsigned solidIdx) const ;
    CSGPrimSpec getPrimSpecHost(   unsigned solidIdx) const ;
    CSGPrimSpec getPrimSpecDevice( unsigned solidIdx) const ;
    void        checkPrimSpec(     unsigned solidIdx) const ;
    void        checkPrimSpec() const ;


    const CSGSolid*   getSolidByName(const char* name) const ;
    const CSGSolid*   getSolid_(int solidIdx) const ;   // -ve counts from back 
    unsigned          getSolidIdx(const CSGSolid* so) const ; 

    unsigned getNumSolid() const; 
    unsigned getNumPrim() const;   
    unsigned getNumNode() const;
    unsigned getNumPlan() const;
    unsigned getNumTran() const;
    unsigned getNumItra() const;
    unsigned getNumInst() const;

    const CSGSolid*   getSolid(unsigned solidIdx) const ;  
    const CSGPrim*    getPrim(unsigned primIdx) const ;    
    const CSGNode*    getNode(unsigned nodeIdx) const ;
    const float4*     getPlan(unsigned planIdx) const ;
    const qat4*       getTran(unsigned tranIdx) const ;
    const qat4*       getItra(unsigned itraIdx) const ;
    const qat4*       getInst(unsigned instIdx) const ;

    CSGSolid* addSolid(unsigned num_prim, const char* label );
    CSGPrim*  addPrim(int num_node) ;
    CSGNode*  addNode(CSGNode nd, const std::vector<float4>* pl=nullptr );
    CSGNode*  addNodes(const std::vector<CSGNode>& nds );
    float4*   addPlan(const float4& pl );

    template<typename T> unsigned addTran( const Tran<T>& tr  );

    CSGSolid* make(const char* name); 
    CSGSolid* makeLayered( const char* label, float outer_radius, unsigned layers ) ;
    CSGSolid* makeClustered(const char* name,  int i0, int i1, int is, int j0, int j1, int js, int k0, int k1, int ks, double unit, bool inbox ) ;

    CSGSolid* makeSolid11(const char* label, CSGNode nd, const std::vector<float4>* pl=nullptr  );

    CSGSolid* makeBooleanBoxSphere( const char* label, char op, float radius, float fullside ) ;

    CSGSolid* makeUnionBoxSphere(        const char* label="ubsp", float radius=100.f, float fullside=150.f );
    CSGSolid* makeIntersectionBoxSphere( const char* label="ibsp", float radius=100.f, float fullside=150.f );
    CSGSolid* makeDifferenceBoxSphere(   const char* label="dbsp", float radius=100.f, float fullside=150.f );

    CSGSolid* makeSphere(     const char* label="sphe", float r=100.f ); 
    CSGSolid* makeEllipsoid(  const char* label="elli", float rx=100.f, float ry=100.f, float rz=50.f ); 

    CSGSolid* makeRotatedCylinder(const char* label="rcyl", float px=0.f, float py=0.f, float radius=100.f, float z1=-50.f, float z2=50.f, float ax=1.f, float ay=0.f, float az=0.f, float angle_deg=45.f  );

    CSGSolid* makeZSphere(    const char* label="zsph", float r=100.f,  float z1=-50.f , float z2=50.f ); 
    CSGSolid* makeCone(       const char* label="cone", float r1=300.f, float z1=-300.f, float r2=100.f,   float z2=-100.f ); 
    CSGSolid* makeHyperboloid(const char* label="hype", float r0=100.f, float zf=50.f,   float z1=-50.f,   float z2=50.f );
    CSGSolid* makeBox3(       const char* label="box3", float fx=100.f, float fy=200.f,  float fz=300.f );
    CSGSolid* makePlane(      const char* label="plan", float nx=1.0f,  float ny=0.f,    float nz=0.f,     float d=0.f );
    CSGSolid* makeSlab(       const char* label="slab", float nx=1.0f,  float ny=0.f,    float nz=0.f,     float d1=-10.f, float d2=10.f );
    CSGSolid* makeCylinder(   const char* label="cyli", float px=0.f,   float py=0.f,    float r=100.f,    float z1=-50.f, float z2=50.f );
    CSGSolid* makeDisc(       const char* label="disc", float px=0.f,   float py=0.f,    float ir=50.f,    float r=100.f,  float z1=-2.f, float z2=2.f);

    CSGSolid* makeConvexPolyhedronCube(       const char* label="vcub", float extent=100.f );
    CSGSolid* makeConvexPolyhedronTetrahedron(const char* label="vtet", float extent=100.f);

    static void DumpAABB(const char* msg, const float* aabb); 

    static float4 TriPlane( const std::vector<float3>& v, unsigned i, unsigned j, unsigned k );

    void write(const char* dir) const ;
    void write(const char* base, const char* rel) const ;
    void load( const char* base, const char* rel ) ; 
    void load( const char* dir ) ; 

    template<typename T> void loadArray( std::vector<T>& vec, const char* dir, const char* name ); 

    void upload();
    void inst_find_unique(); 

    unsigned getNumUniqueIAS() const ;
    unsigned getNumUniqueGAS() const ;
    unsigned getNumUniqueINS() const ;
    unsigned getNumInstancesIAS(unsigned ias_idx) const ;
    void     getInstanceTransformsIAS(std::vector<qat4>& ias_inst, unsigned ias_idx ) const ;

    std::vector<CSGSolid>  solid ;   
    std::vector<CSGPrim>   prim ; 
    std::vector<CSGNode>   node ; 
    std::vector<float4>    plan ; 
    std::vector<qat4>      tran ;  
    std::vector<qat4>      itra ;  
    std::vector<qat4>      inst ;  

    CSGSolid*   d_solid ; 
    CSGPrim*    d_prim ; 
    CSGNode*    d_node ; 
    float4*  d_plan ; 
    qat4*    d_tran ; 
    qat4*    d_itra ; 


    std::vector<unsigned>  ins ; 
    std::vector<unsigned>  gas ; 
    std::vector<unsigned>  ias ; 



};


