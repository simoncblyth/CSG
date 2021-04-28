#!/bin/bash -l 

source ./env.sh 

bins="CSGNodeTest CSGPrimTest CSGSolidTest CSGFoundryTest ScanTest"
for bin in $bins ; do
   echo $msg $(which $bin)
   $bin
   [ $? -ne 0 ] && echo runtime error from $bin && exit 1
done

exit 0

