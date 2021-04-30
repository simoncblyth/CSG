#!/bin/bash 

sdir=$(pwd)
name=CSG

export CSG_PREFIX=/usr/local/csg
export PATH=$CSG_PREFIX/lib:$PATH

case $(uname) in
  Darwin) var=DYLD_LIBRARY_PATH ; lib="lib" ;;
  Linux)  var=LD_LIBRARY_PATH ; lib="lib64"  ;;
esac
export $var=$CSG_PREFIX/$lib  


msg="=== $0 :"
echo $msg sdir $sdir 
echo $msg name $name 
echo $msg CSG_PREFIX $CSG_PREFIX

echo $msg PATH
echo $PATH | tr ":" "\n"






#geometry=parade
#geometry=sphere_containing_grid_of_spheres
#geometry=layered_sphere
#geometry=layered_zsphere
#geometry=clustered_sphere
#geometry=sphe # 0
geometry=zsph # 1 
#geometry=cone # 2
#geometry=hype # 3
#geometry=box3 # 4 
#geometry=plan # 5 
#geometry=slab # 6  
#geometry=cyli # 7
#geometry=disc # 8 
#geometry=vcub # 9
#geometry=vtet # 10
#geometry=elli # 11
#geometry=ubsp # 12 
#geometry=ibsp # 13 
#geometry=dbsp # 14
#geometry=rcyl  # 15



export GEOMETRY=${GEOMETRY:-$geometry}

