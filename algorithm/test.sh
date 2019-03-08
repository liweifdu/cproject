#!/bin/bash

rm *.o

g++ -g -std=c++11 $1 -o "$1.o"
./$1.o
