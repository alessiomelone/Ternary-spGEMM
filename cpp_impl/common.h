#pragma once

#include <string>
#include <functional>
#include "data_structures/CompressedCSC.h"
#include "data_structures/ICSC.h"
#include "data_structures/ICSR.h"
#include "data_structures/BaseTCSC.h"
#include "data_structures/BaseTCSR.h"
#include "data_structures/BlockedTCSC.h"
#include "data_structures/BlockedTCSC_interleaved.h"
#include "data_structures/InterleavedTCSC.h"
#include "data_structures/InterleavedTCSCPadding.h"
#include "data_structures/InterleavedTCSC_baraq.h"

using comp_func = std::function<void(double *X, double *B, double *Y, int M, int N, int K)>;

void add_function(comp_func f, std::string name);