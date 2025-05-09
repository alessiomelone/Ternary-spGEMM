#pragma once

#include <string>
#include <functional>
#include "SparseGEMM.h"


using comp_func = std::function<void(float *X, float *B, float *Y, int M, int N, int K)>;

void add_function(comp_func f, std::string name);