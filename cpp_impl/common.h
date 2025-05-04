#pragma once

#include <string>
#include "SparseGEMM.h"

typedef void (*comp_func)(float *X, int *col_start_pos, int *col_start_neg, int *row_index_pos, int *row_index_neg, float *b, float *Y, int M, int N, int K);
void add_function(comp_func f, std::string name);