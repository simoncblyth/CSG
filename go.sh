#!/bin/bash -l


source ./build.sh 
[ $? -ne 0 ] && exit 2

source ./run.sh 
[ $? -ne 0 ] && exit 3

exit 0 

