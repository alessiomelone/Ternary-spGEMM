# ASL_SparseGEMM

A project template for ASL Spring 2025.

## C++ Implementation

### Building the Project

To build the C++ implementation, run:
```bash
make
```

This will create an executable named `sparseGEMM` in the root directory.

### Running the Program

To run the compiled program (example):
```bash
sudo ./sparseGEMM.out -M 32 -K 1024 -N 4096 -s 4
```

### Registering New Functions

To register a new computational function:
1. Define your templated function in `cpp_impl/comp.h`. This header is where all templated computational kernels should reside.
2. In `cpp_impl/main.cpp`, register your function with the driver:
   - Create a lambda function that calls your templated function with specific template arguments (e.g., `float`, specific unroll factors)
   - Pass this lambda and a descriptive name string to the `add_function` utility

Example:
```cpp
// In main.cpp
add_function("MyFunction", [](const Matrix& A, const Matrix& B, Matrix& C) {
    my_templated_function<float, 4>(A, B, C);  // Example with float type and unroll factor 4
});
```

### Adding New Data Structures

To add a new data structure:
1. Create a new header file in `cpp_impl/data_structures/` (e.g., `MyDataStructure.h`)
2. Implement the `DataStructureInterface` interface defined in `cpp_impl/data_structures/DataStructureInterface.hpp`
3. Include your new header in `cpp_impl/common.h` to make it accessible throughout the project
4. You can then use this data structure within your algorithm definitions in `cpp_impl/comp.h`

Example data structure implementation:
```cpp
// MyDataStructure.h
class MyDataStructure : public DataStructureInterface {
public:
    void init(const int* matrix, int rows, int cols) override;
    std::vector<int> getVectorRepresentation(size_t expected_rows, size_t expected_cols) override;
    int getNumRows() const override;
    int getNumCols() const override;
    // ... other methods
};
```

### Testing Data Structures

To test a data structure:
```bash
g++ -Icpp_impl/data_structures -o test_data_structure.out cpp_impl/test_data_structure.cpp && ./test_data_structure.out
```


## Performance Considerations

### Flops of Base Implementation
$MN\left(1 + \frac{K}{\text{nonZero}}\right)$

### Memory Usage

- **X**: $MK \cdot \text{sizeof}(T)$ — One row at a time is accessed then never used again ($K \cdot \text{sizeof}(T)$)
- **B**: $MN \cdot \text{sizeof}(T)$ — In the inner loop, one row is accessed then never used again ($N \cdot \text{sizeof}(T)$)
- **Y**: $MN \cdot \text{sizeof}(T)$ — In the inner loop, one row is accessed then never used again ($N \cdot \text{sizeof}(T)$)

**Total**: $M(K + 2N) \cdot \text{sizeof}(T)$

### Data Structure Formats

#### TCSC Format
- **CSP and CSN**: $2(N + 1) \cdot \text{sizeof(int)} = 2N \cdot \text{sizeof(int)}$
- **RIP and RIN**: $\frac{KN}{\text{nonZero}} \cdot \text{sizeof(int)}$

**Total TCSC**: $N\left(\frac{K}{\text{nonZero}} + 2\right) \cdot \text{sizeof(int)}$

### Total Memory Usage (Uncompressed)
Assuming $T = \text{float}$, $\text{sizeof(float)} = 4$ bytes, $\text{sizeof(int)} = 4$ bytes:

$4(MK + 2MN + 2N + \frac{KN}{\text{nonZero}})$ bytes

## Optimization Approaches

1. Compress TCSC format while maintaining code structure
2. Use normal CSC with compressed values vector (1s and -1s, 8 bits for 5 values)
