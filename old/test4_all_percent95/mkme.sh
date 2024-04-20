#!/bin/bash
rm boot.exe
#clang++ --std=c++17 boot.cxx -o boot.exe -g -O0 -Wall -Wextra -Wshadow -I/Users/37348458/bootstrap/c++ -I/Users/37348458/randutils
clang++ --std=c++17 boot.cxx -o boot.exe -O3 -I/Users/37348458/bootstrap/c++ -I/Users/37348458/randutils
