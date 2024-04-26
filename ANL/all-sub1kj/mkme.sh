#!/bin/bash
rm boot.exe
#clang++ --std=c++17 boot.cxx -o boot.exe -g -O0 -Wall -Wextra -Wshadow -I/Users/37348458/bootstrap/c++ -I/Users/37348458/randutils
clang++ --std=c++17 boot.cxx -o boot.exe -O3 -fopenmp -I/Users/37348458/bootstrap -I/Users/37348458/randutils -lomp
#clang++ --std=c++17 boot.cxx -o boot.exe -O0 -g -D_GLIBCXX_DEBUG -fopenmp -fsanitize-address-use-after-return=always -fsanitize-address-use-after-scope -I/Users/37348458/bootstrap -I/Users/37348458/randutils -lomp

