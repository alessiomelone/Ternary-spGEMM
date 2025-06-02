#pragma once
#include <vector>
// If we want to use test_data_structure.cpp for debugging, we need to implement this interface.
class DataStructureInterface
{
public:
    virtual ~DataStructureInterface() = default;

    // Initialize the data structure from a raw matrix
    virtual void init(const int *matrix, int rows, int cols) = 0;

    // Convert back to a dense matrix representation
    virtual std::vector<int> getVectorRepresentation(size_t rows, size_t cols) = 0;
};