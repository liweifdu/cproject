#!/bin/bash

#rmove the all project files
rm *.o

#parameter
CPP_NAME="$1"
SRC_FILE="$CPP_NAME.cpp"
OUT_FILE="$CPP_NAME.o"

#command
g++ -g -std=c++11 -Wall -Wextra -pedantic $SRC_FILE -o $OUT_FILE
./$OUT_FILE
if [ -n "$2" ]
then
    echo Hello $1, glad to meet you
else
    ./$OUT_FILE
fi
