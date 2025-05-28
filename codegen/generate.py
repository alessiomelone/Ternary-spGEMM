# void BaseCSC(T *X, const BaseTCSC &W_csc, T *b, T *Y, int M, int N, int K)
# {
#     const int *col_start_pos = W_csc.col_start_pos.data();
#     const int *col_start_neg = W_csc.col_start_neg.data();
#     const int *row_index_pos = W_csc.row_index_pos.data();
#     const int *row_index_neg = W_csc.row_index_neg.data();

#     for (int m = 0; m < M; m++)
#     {
#         for (int n = 0; n < N; n++)
#         {
#             T y_val = 0;

#             // Process positive values
#             for (int k = col_start_pos[n]; k < col_start_pos[n + 1]; k++)
#             {
#                 T x_val = X[m * K + row_index_pos[k]];
#                 y_val += x_val;
#             }

#             // Process negative values
#             for (int k = col_start_neg[n]; k < col_start_neg[n + 1]; k++)
#             {
#                 T x_val = X[m * K + row_index_neg[k]];
#                 y_val -= x_val;
#             }

#             Y[m * N + n] = y_val + b[n];
#         }
#     }
# }

from M_16_K_1024_N_4096_s_16 import W, X, Y_expected, Y_actual, b, col_start_pos, col_start_neg, row_index_pos, row_index_neg, M, N, K, s
import numpy as np

def base_csc_flat(X_flat, col_start_pos, col_start_neg,
                  row_index_pos, row_index_neg, b,
                  M, N, K):
    """
    X_flat:         length M*K flat list or 1D ndarray
    col_start_pos:  length N+1 list of ints
    col_start_neg:  length N+1 list of ints
    row_index_pos:  list of ints (pos branch)
    row_index_neg:  list of ints (neg branch)
    b:              length N list or 1D ndarray
    M, N, K:        ints

    Returns Y_flat: length M*N 1D ndarray of float32
    """
    # ensure numpy arrays for fast indexing and correct dtype
    X = np.asarray(X_flat, dtype=np.float32)
    b = np.asarray(b, dtype=np.float32)

    Y = np.empty(M * N, dtype=np.float32)

    for m in range(M):
        base_x = m * K
        base_y = m * N
        for n in range(N):
            y_val = 0.0

            # positive contributions
            start_p = col_start_pos[n]
            end_p   = col_start_pos[n+1]
            for idx in range(start_p, end_p):
                row = row_index_pos[idx]
                y_val += X[base_x + row]

            # negative contributions
            start_n = col_start_neg[n]
            end_n   = col_start_neg[n+1]
            for idx in range(start_n, end_n):
                row = row_index_neg[idx]
                y_val -= X[base_x + row]

            Y[base_y + n] = y_val + b[n]

    return Y

def gen_code(X, col_start_pos, col_start_neg,
             row_index_pos, row_index_neg, b,
             M, N, K, log_filename):
    """
    Writes the exact sequence of ops (+X, -X, bias add) to a file,
    and computes Y using only Python lists, no numpy.
    """
    Y = [0.0] * (M * N)

    with open(log_filename, 'w') as f:
        f.write("inline void GenCSC(float *X, float *b, float *Y) {")
        f.write("float y_val;")
        for m in range(M):
            base_x = m * K
            base_y = m * N
            for n in range(N):
                # y_val = 0.0
                f.write("y_val = 0;")

                for idx in range(col_start_pos[n], col_start_pos[n+1]):
                    row = row_index_pos[idx]
                    f.write(f"y_val += X[{base_x + row}];\n")
                    # y_val += X[base_x + row]

                for idx in range(col_start_neg[n], col_start_neg[n+1]):
                    row = row_index_neg[idx]
                    f.write(f"y_val -= X[{base_x + row}];\n")
                    # y_val -= X[base_x + row]

                Y_idx = base_y + n
                f.write(f"Y[{Y_idx}] = y_val + b[{n}];\n")
                # Y[Y_idx] = y_val + b[n]
        f.write("}")
    # return Y

log_filename = "codegen/generated.h"
gen_code(X, col_start_pos, col_start_neg, row_index_pos, row_index_neg, b, M, N, K,
         log_filename)

# Y_py = base_csc_flat(X, col_start_pos, col_start_neg,
#                      row_index_pos, row_index_neg, b,
#                      M, N, K)

# print("max abs error:", np.max(np.abs(Y_py - Y_actual)))
# print("allclose?   ", np.allclose(Y_py, Y_actual, atol=1e-6))
