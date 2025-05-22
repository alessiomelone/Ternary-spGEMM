#pragma once

#include <string>
#include <functional>
#include "data_structures/CompressedCSC.h"
#include "data_structures/TCSRMatrix.h"
#include "data_structures/TCSCMatrix.h"
#include "data_structures/SparseFormatCSC.h"
#include "data_structures/SparseFormatCSR.h"

using comp_func = std::function<void(float *X, float *B, float *Y, int M, int N, int K)>;

void add_function(comp_func f, std::string name);