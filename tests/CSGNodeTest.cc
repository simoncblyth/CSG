// ./CSGNodeTest.sh 

#include <vector>
#include <iomanip>
#include <iostream>
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include "sutil_vec_math.h"
#include "CSGNode.h"
#include "Sys.h"

void test_copy()
{
   glm::mat4 m0(1.f); 
   glm::mat4 m1(2.f); 
   glm::mat4 m2(3.f); 

   m0[0][3] = Sys::int_as_float(42); 
   m1[0][3] = Sys::int_as_float(52); 
   m2[0][3] = Sys::int_as_float(62); 

   std::vector<glm::mat4> node ; 
   node.push_back(m0); 
   node.push_back(m1); 
   node.push_back(m2); 

   std::vector<CSGNode> node_(3) ; 

   memcpy( node_.data(), node.data(), sizeof(CSGNode)*node_.size() );  

   CSGNode* n_ = node_.data(); 
   CSGNode::Dump( n_, node_.size(), "CSGNodeTest" );  
}


void test_zero()
{

   CSGNode nd = {} ; 
   assert( nd.gtransformIdx() == 0u );  
   assert( nd.complement() == false );  

   unsigned tr = 42u ; 
   nd.setTransform( tr ); 

   assert( nd.gtransformIdx() == tr ); 
   assert( nd.complement() == false ); 

   nd.setComplement(true); 
   assert( nd.gtransformIdx() == tr ); 
   assert( nd.complement() == true ); 

   std::cout << nd.desc() << std::endl ; 
}


void test_sphere()
{
   CSGNode nd = CSGNode::Sphere(100.f); 
   std::cout << nd.desc() << std::endl ; 
}


int main(int argc, char** argv)
{
   std::cout << argv[0] << std::endl ; 

   test_zero(); 
   test_sphere(); 
   test_copy();  

   return 0 ; 
}
