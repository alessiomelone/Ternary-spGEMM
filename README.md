# ASL_SparseGEMM

A project template for ASL Spring 2025.

## For C Implementation:

#### For x86 systems

gcc -O3 ./c_impl/sparse_format.c ./c_impl/TestImpl.c -o ./c_impl/TestImpl.out

#### For ARM systems (e.g., Apple M1/M2)

gcc -DPMU -arch arm64 -O3 ./c_impl/sparse_format.c ./c_impl/TestImpl.c -o ./c_impl/TestImpl.out

sudo ./c_impl/TestImpl.out -M 32 -K 1024 -N 4096 -s 4

## For C++ Implementation:

#### For x86 systems

g++ -O3 SparseGEMM.cpp -o SparseGEMM.out

#### For ARM systems (e.g., Apple M1/M2)

g++ -arch arm64 -O3 SparseGEMM.cpp -o SparseGEMM.out -DPMU

sudo ./SparseGEMM.out -M 32 -K 1024 -N 4096 -s 4