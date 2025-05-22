#pragma once

#include <string>
#include <functional>
#include "data_structures/CompressedCSC.h"
#include "data_structures/TCSC.h"
#include "data_structures/TCSR.h"
#include "data_structures/BaseTCSC.h"
#include "data_structures/BaseTCSR.h"

using comp_func = std::function<void(float *X, float *B, float *Y, int M, int N, int K)>;

void add_function(comp_func f, std::string name);