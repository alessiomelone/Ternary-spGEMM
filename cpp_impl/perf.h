typedef void (*comp_func)(float *X, int *col_start_pos, int *col_start_neg, int *row_index_pos, int *row_index_neg, float *b, float *Y, int M, int N, int K);

double perf_test(comp_func f, int M, int K, int N, int nonZero);
