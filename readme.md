# ASL_SparseGEMM

A project template for ASL Spring 2025.

## C++ Implementation

To register a new computational function:
1.  Define your templated function in `cpp_impl/comp.h`. This header is where all templated computational kernels should reside.
2.  In `cpp_impl/main.cpp`, register your function with the driver. This involves:
    *   Creating a lambda function that calls your templated function with the specific template arguments (e.g., `float`, specific unroll factors).
    *   Passing this lambda and a descriptive name string to the `add_function` utility.

To add a new data structure:
1. Create a `.h` file for your data structure (e.g., `my_data_structure.h`).
2. Include this new header in `cpp_impl/common.h` to make it accessible throughout the project.
3. You can then use this data structure within your algorithm definitions in `cpp_impl/comp.h`.

To build the C++ implementation, navigate to the `cpp_impl` directory and run `make`:
```bash
cd cpp_impl
make
```

To run the compiled program (example):
```bash
sudo ./SparseGEMM.out -M 32 -K 1024 -N 4096 -s 4
```

To clean the build (remove the executable):
```bash
make clean
```

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
