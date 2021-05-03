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



