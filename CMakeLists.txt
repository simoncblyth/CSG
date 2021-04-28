cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
set(name CSG)
project(${name})

set(CMAKE_CXX_STANDARD 11) 
set(CMAKE_CXX_STANDARD_REQUIRED YES)


include(GNUInstallDirs) 

find_package(CUDA REQUIRED) 



include_directories(
    ${CUDA_INCLUDE_DIRS}
    ${CMAKE_INSTALL_PREFIX}/externals/glm/glm
    ${CMAKE_CURRENT_SOURCE_DIR} 
)


set(SOURCES)
set(HEADERS)


list(APPEND SOURCES
    CSGNode.cc
    CU.cc
    CSGPrim.cc
    CSGPrimSpec.cc
    CSGSolid.cc

    Tran.cc
    CSGFoundry.cc
)

list(APPEND HEADERS
    sutil_vec_math.h
    Quad.h 
    OpticksCSG.h
    CSGNode.h 

    qat4.h
    CSGPrimSpec.h
    CSGPrim.h
    CSGSolid.h

    Tran.h
    Util.h
    AABB.h
    CSGFoundry.h

)

message(STATUS "CMAKE_INSTALL_PREFIX     :${CMAKE_INSTALL_PREFIX}")
message(STATUS "CMAKE_INSTALL_INCLUDEDIR :${CMAKE_INSTALL_INCLUDEDIR}")
message(STATUS "CMAKE_INSTALL_LIBDIR     :${CMAKE_INSTALL_LIBDIR}")


add_library( ${name}  SHARED ${SOURCES} ${HEADERS} )

target_link_libraries(${name} 
     PUBLIC 
         ${CUDA_LIBRARIES}
)



install(FILES ${HEADERS}  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(TARGETS ${name}  DESTINATION ${CMAKE_INSTALL_LIBDIR})

add_subdirectory(tests)


#[=[

Sys.cc
Geo.cc
Grid.cc
Scan.cc
Tran.cc
Util.cc
View.cc



#]=]





