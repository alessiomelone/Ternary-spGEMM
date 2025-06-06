#pragma once

#include <string>
#include <functional>
#include "data_structures/TCSC.h"
#include "data_structures/TCSR.h"
#include "data_structures/BlockedTCSC.h"
#include "data_structures/InterleavedBlockedTCSC.h"
#include "data_structures/InterleavedTCSC.h"
#include "data_structures/VectorTCSC.h"

using comp_func = std::function<void(float *X, float *B, float *Y, int M, int N, int K)>;
using comp_func_prelu = std::function<void(float *X, float *B, float *alpha, float *Y, int M, int N, int K)>;

void add_function(comp_func f, std::string name);
void add_prelu_function(comp_func_prelu f, std::string name);