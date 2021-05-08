#!/bin/bash 

msg="=== $BASH_SOURCE :"
sdir=$(pwd)
name=$(basename $sdir)

bdir=/tmp/$USER/opticks/$name/build 
echo $msg bdir $bdir 

rm -rf $bdir && mkdir -p $bdir 
[ ! -d $bdir ] && exit 1

cd $bdir && pwd 

echo $msg CMAKE_PREFIX_PATH 
echo $CMAKE_PREFIX_PATH | tr ":" "\n"


[ -z "$OPTICKS_PREFIX" ] && echo $msg MISSING OPTICKS_PREFIX && exit 1

cmake $sdir \
     -DCMAKE_BUILD_TYPE=Debug \
     -DOPTICKS_PREFIX=$OPTICKS_PREFIX \
     -DCMAKE_MODULE_PATH=${OPTICKS_HOME}/cmake/Modules \
     -DCMAKE_INSTALL_PREFIX=${OPTICKS_PREFIX}


make
[ $? -ne 0 ] && echo $0 : make FAIL && exit 1
make install   
[ $? -ne 0 ] && echo $0 : install FAIL && exit 2

exit 0

