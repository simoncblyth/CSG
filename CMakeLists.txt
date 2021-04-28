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

    Scan.cc
    View.cc
    Util.cc
    Grid.cc
    Geo.cc
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
    AABB.h
    CSGFoundry.h

    Scan.h
    View.h

    Util.h
    Grid.h
    Geo.h
)


#[=[

#]=]

message(STATUS "CMAKE_INSTALL_PREFIX     : ${CMAKE_INSTALL_PREFIX}")
message(STATUS "CMAKE_INSTALL_INCLUDEDIR : ${CMAKE_INSTALL_INCLUDEDIR}")
message(STATUS "CMAKE_INSTALL_LIBDIR     : ${CMAKE_INSTALL_LIBDIR}")
message(STATUS "CMAKE_BINARY_DIR         : ${CMAKE_BINARY_DIR}")


add_library( ${name}  SHARED ${SOURCES} ${HEADERS} )
target_link_libraries(${name} 
     PUBLIC 
         ${CUDA_LIBRARIES}
)

target_include_directories(${name}
     PUBLIC
           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)



install(FILES ${HEADERS}  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
add_subdirectory(tests)


set(CSG_INCLUDE_DIRS "include")



#[=[
https://cmake.org/cmake/help/git-master/guide/tutorial/index.html#adding-export-configuration-step-11

https://cmake.org/cmake/help/git-stage/guide/importing-exporting/index.html
#]=]


set(CONFIG_DEST ${CMAKE_INSTALL_LIBDIR}/cmake/CSG) 


install(TARGETS ${name}  DESTINATION ${CMAKE_INSTALL_LIBDIR} EXPORT CSGTargets )

install(EXPORT CSGTargets
        FILE CSGTargets.cmake
        DESTINATION ${CONFIG_DEST}
)

include(CMakePackageConfigHelpers)

configure_package_config_file ( 
       "CSGConfig.cmake.in" "${CMAKE_BINARY_DIR}/CSGConfig.cmake"
       INSTALL_DESTINATION "${CONFIG_DEST}"
       PATH_VARS CSG_INCLUDE_DIRS 
)

install(FILES ${CMAKE_BINARY_DIR}/CSGConfig.cmake DESTINATION ${CONFIG_DEST})








