#include "BlockedTCSC.h"
#include <iostream>

void printVector(const std::vector<int> &vec, const std::string &name)
{
    std::cout << name << ": ";
    for (int val : vec)
    {
        std::cout << val << " ";
    }
    std::cout << std::endl;
}

int main()
{
    // Example matrix (4x4) with B=2
    // 1  0 -1  0
    // 0  1  0 -1
    // 1 -1  0  0
    // 0  0  1 -1
    int matrix[] = {
        1, 1, -1, 1,
        1, 1, 1, -1,
        1, -1, 0, 0,
        0, 0, 1, -1};

    BlockedTCSC<2> tcsc(matrix, 4, 4);

    std::cout << "Testing BlockedTCSC with B=2, K=4, N=4\n\n";

    printVector(tcsc.col_start_pos, "col_start_pos");
    printVector(tcsc.col_start_neg, "col_start_neg");
    printVector(tcsc.row_index_pos, "row_index_pos");
    printVector(tcsc.row_index_neg, "row_index_neg");

    return 0;
}