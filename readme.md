# ASL_SparseGEMM

A project template for ASL Spring 2025.

## C++ Implementation

Edit only `comp.cpp`:

```bash
g++ -O2 -march=native -mtune=native -fstrict-aliasing -DNDEBUG cpp_impl/main.cpp cpp_impl/comp.cpp cpp_impl/perf.cpp -o cpp_impl/SparseGEMM.out -DPMU
sudo ./cpp_impl/SparseGEMM.out -M 32 -K 1024 -N 4096 -s 4
```

## C Implementation

### For x86 systems

```bash
gcc -O3 ./c_impl/sparse_format.c ./c_impl/TestImpl.c -o ./c_impl/TestImpl.out
```

### For ARM systems (e.g., Apple M1/M2)

```bash
gcc -DPMU -arch arm64 -O3 ./c_impl/sparse_format.c ./c_impl/TestImpl.c -o ./c_impl/TestImpl.out
sudo ./c_impl/TestImpl.out -M 32 -K 1024 -N 4096 -s 4
```

---

## Considerations

### Flops of base implementation:

$MN\left(1 + \frac{K}{\text{nonZero}}\right)$

### Data size:

- **X**: $MK \cdot \text{sizeof}(T)$ — One row at a time is accessed then never used again ($K \cdot \text{sizeof}(T)$)
- **B**: $MN \cdot \text{sizeof}(T)$ — In the inner loop, one row is accessed then never used again ($N \cdot \text{sizeof}(T)$)
- **Y**: $MN \cdot \text{sizeof}(T)$ — In the inner loop, one row is accessed then never used again ($N \cdot \text{sizeof}(T)$)

**Total**: $M(K + 2N) \cdot \text{sizeof}(T)$


### Format:

Considering $W$ is $K \times N$

#### TCSC:

- **CSP and CSN**, total: $2(N + 1) \cdot \text{sizeof(int)} = 2N \cdot \text{sizeof(int)}$

TODO: Can it be reduced?  
If we think we just need the information about how many values are being read in that iteration (max of $n_{\text{rows}}$ values, so $K$), we could use $\log_2(K)$ bits for each value in this vector.

- **RIP and RIN**, total: $\frac{KN}{\text{nonZero}} \cdot \text{sizeof(int)}$

TODO: Can it be reduced?  
We only need to know in which row the values are located. We could also use $\log_2(K)$ bits.

**Total TCSC**: $N\left(\frac{K}{\text{nonZero}} + 2\right) \cdot \text{sizeof(int)}$

---

### Total, uncompressed:

Assuming $T = \text{float}$, $\text{sizeof(float)} = 2$ bytes, $\text{sizeof(int)} = 2$ bytes:

$2(MK + 2MN + 2N + \frac{KN}{\text{nonZero}})$

---

## Operational Intensity

1. First approach: compress TCSC, see how the code looks like — this would allow keeping the structure of the code  
2. Second approach: normal CSC, compressed values vector (only 1s and -1s, 8 bits for 5 values)  
   We need to try to code and think about how the code would look like
