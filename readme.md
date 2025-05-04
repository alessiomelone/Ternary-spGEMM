# ASL_SparseGEMM

A project template for ASL Spring 2025.

## For C++ Implementation, only edit comp.cpp:

g++ -O3 cpp_impl/main.cpp cpp_impl/comp.cpp cpp_impl/perf.cpp -o cpp_impl/SparseGEMM.out -DPMU

sudo ./cpp_impl/SparseGEMM.out -M 32 -K 1024 -N 4096 -s 4

## For C Implementation:

#### For x86 systems

gcc -O3 ./c_impl/sparse_format.c ./c_impl/TestImpl.c -o ./c_impl/TestImpl.out

#### For ARM systems (e.g., Apple M1/M2)

gcc -DPMU -arch arm64 -O3 ./c_impl/sparse_format.c ./c_impl/TestImpl.c -o ./c_impl/TestImpl.out

sudo ./c_impl/TestImpl.out -M 32 -K 1024 -N 4096 -s 4
