#!/bin/bash -l 

source ./env.sh 

bins="CSGNodeTest CSGPrimTest CSGSolidTest CSGFoundryTest"
for bin in $bins ; do
   echo $msg $(which $bin)
   $bin
done



