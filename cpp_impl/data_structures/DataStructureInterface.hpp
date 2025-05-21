#pragma once
#include <vector>
#include <cstddef>

class DataStructureInterface {
public:
    // Initialize data structure with a raw pointer to an int vector plus dimensions
    virtual void init(const int* matrix, int rows, int cols) = 0;

    // Return a vector<int> representation given requested dimensions
    virtual std::vector<int> getVectorRepresentation(std::size_t rows, std::size_t cols) = 0;
};