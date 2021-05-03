#!/bin/bash 

source ./env.sh 

bdir=$CSG_PREFIX/build 
echo $msg bdir $bdir 

rm -rf $bdir && mkdir -p $bdir 
[ ! -d $bdir ] && exit 1

cd $bdir && pwd 

echo $BASH_SOURCE CMAKE_PREFIX_PATH 
echo $CMAKE_PREFIX_PATH | tr ":" "\n"

cmake $sdir \
     -DCMAKE_BUILD_TYPE=Debug \
     -DOPTICKS_PREFIX=$OPTICKS_PREFIX \
     -DCMAKE_INSTALL_PREFIX=$CSG_PREFIX


rm -rf   $CSG_PREFIX/lib
mkdir -p $CSG_PREFIX/lib 

make
[ $? -ne 0 ] && echo $0 : make FAIL && exit 1
make install   
[ $? -ne 0 ] && echo $0 : install FAIL && exit 2

exit 0

