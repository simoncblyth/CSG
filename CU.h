#pragma once
#include <vector>

struct CU
{
    template <typename T>
    static T* UploadArray(const T* array, unsigned num_items ) ; 

    template <typename T>
    static T* DownloadArray(const T* array, unsigned num_items ) ; 



    template <typename T>
    static T* UploadVec(const std::vector<T>& vec);

    template <typename T>
    static void DownloadVec(std::vector<T>& vec, const T* d_array, unsigned num_items);

};
