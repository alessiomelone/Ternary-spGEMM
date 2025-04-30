# ASL_SparseGEMM

A project template for ASL Spring 2025.

## For C Implementation:

#### For x86 systems

gcc -O3 c_impl/ref.c -o c_impl/ref.out

#### For ARM systems (e.g., Apple M1/M2)

gcc -arch arm64 -O3 c_impl/ref.c -o c_impl/ref.out

## For C++ Implementation:

#### For x86 systems

g++ -O3 SparseGEMM.cpp -o SparseGEMM.out

#### For ARM systems (e.g., Apple M1/M2)

g++ -arch arm64 -O3 SparseGEMM.cpp -o SparseGEMM.out -DPMU

./SparseGEMM.out -M 16 -K 1024 -N 4096 -s 16