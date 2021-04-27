#include <iostream>
#include "sutil_vec_math.h"    

#include "cuda_runtime.h"
#include "CUDA_CHECK.h"

#include "CSGSolid.h"
#include "CSGPrim.h"
#include "CSGNode.h"
#include "Quad.h"
#include "qat4.h"

#include "CU.h"

/**
CU::UploadArray
----------------

Allocate on device and copy from host to device

**/
template <typename T>
T* CU::UploadArray(const T* array, unsigned num_items ) // static
{
    std::cout << "CU::UploadArray num_items " << num_items << std::endl ; 
    T* d_array = nullptr ; 
    CUDA_CHECK( cudaMalloc(reinterpret_cast<void**>( &d_array ), num_items*sizeof(T) ));
    CUDA_CHECK( cudaMemcpy(reinterpret_cast<void*>( d_array ), array, sizeof(T)*num_items, cudaMemcpyHostToDevice ));
    return d_array ; 
}

/**
CU::UploadArray  
----------------

Allocate on host and copy from device to host 

**/

template <typename T>
T* CU::DownloadArray(const T* d_array, unsigned num_items ) // static
{
    std::cout << "CU::DownloadArray num_items " << num_items << std::endl ; 
    T* array = new T[num_items] ;  
    CUDA_CHECK( cudaMemcpy( array, d_array, sizeof(T)*num_items, cudaMemcpyDeviceToHost ));
    return array ; 
}



template float* CU::UploadArray<float>(const float* array, unsigned num_items) ;
template float* CU::DownloadArray<float>(const float* d_array, unsigned num_items) ;

template unsigned* CU::UploadArray<unsigned>(const unsigned* array, unsigned num_items) ;
template unsigned* CU::DownloadArray<unsigned>(const unsigned* d_array, unsigned num_items) ;

template float4* CU::UploadArray<float4>(const float4* array, unsigned num_items) ;
template float4* CU::DownloadArray<float4>(const float4* d_array, unsigned num_items) ;

template CSGNode* CU::UploadArray<CSGNode>(const CSGNode* d_array, unsigned num_items) ;
template CSGNode* CU::DownloadArray<CSGNode>(const CSGNode* d_array, unsigned num_items) ;

template quad4* CU::UploadArray<quad4>(const quad4* d_array, unsigned num_items) ;
template quad4* CU::DownloadArray<quad4>(const quad4* d_array, unsigned num_items) ;

template qat4* CU::UploadArray<qat4>(const qat4* d_array, unsigned num_items) ;
template qat4* CU::DownloadArray<qat4>(const qat4* d_array, unsigned num_items) ;

template CSGPrim* CU::UploadArray<CSGPrim>(const CSGPrim* d_array, unsigned num_items) ;
template CSGPrim* CU::DownloadArray<CSGPrim>(const CSGPrim* d_array, unsigned num_items) ;

template CSGSolid* CU::UploadArray<CSGSolid>(const CSGSolid* d_array, unsigned num_items) ;
template CSGSolid* CU::DownloadArray<CSGSolid>(const CSGSolid* d_array, unsigned num_items) ;



