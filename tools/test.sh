#!/bin/bash

#rmove the all project files
rm *.o

#parameter
CPP_NAME="$1"
SRC_FILE="$CPP_NAME.cpp"
OUT_FILE="$CPP_NAME.o"

#command
g++ -g -std=c++11 -Wall -Wextra -pedantic $SRC_FILE -o $OUT_FILE

if [ -n "$2" ]
then
    if [ -n "$3" ] 
    then
        ./$OUT_FILE $2 $3
    else
        ./$OUT_FILE $2 
    fi
else
    ./$OUT_FILE
fi
