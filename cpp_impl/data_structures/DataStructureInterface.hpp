#pragma once
#include <vector>

class DataStructureInterface
{
public:
    virtual ~DataStructureInterface() = default;

    // Initialize the data structure from a raw matrix
    virtual void init(const int *matrix, int rows, int cols) = 0;

    // Convert back to a dense matrix representation
    virtual std::vector<int> getVectorRepresentation(size_t expected_rows, size_t expected_cols) = 0;

    // Get dimensions
    virtual int getNumRows() const = 0;
    virtual int getNumCols() const = 0;

    // Debug/print functionality
    virtual void printVars() = 0;
};