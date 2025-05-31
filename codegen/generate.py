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

from M_1_K_512_N_2048_s_16 import W, X, Y_expected, Y_actual, b, col_start_pos, col_start_neg, row_index_pos, row_index_neg, M, N, K, s
import numpy as np

def gen_code(X, col_start_pos, col_start_neg, row_index_pos, row_index_neg, b, M, N, K, log_filename):
    """
    Writes optimized sequence of ops (+X, -X, bias add) using batching and tree reductions.
    """
    BATCH_SIZE = 8
    with open(log_filename, 'w') as f:
        f.write("inline void GenCSC(float *X, float *b, float *Y) {\n")
        for m in range(M):
            base_x = m * K
            base_y = m * N
            for n in range(N):
                # collect positive and negative indices
                start_p = col_start_pos[n]
                end_p = col_start_pos[n+1]
                pos_indices = row_index_pos[start_p:end_p]
                start_n = col_start_neg[n]
                end_n = col_start_neg[n+1]
                neg_indices = row_index_neg[start_n:end_n]
                # build list of entries: (sign, index)
                entries = [(+1, idx) for idx in pos_indices] + [(-1, idx) for idx in neg_indices]
                # split into batches
                num_entries = len(entries)
                num_batches = (num_entries + BATCH_SIZE - 1) // BATCH_SIZE
                partial_sums = []
                for batch_id in range(num_batches):
                    batch_entries = entries[batch_id * BATCH_SIZE:(batch_id + 1) * BATCH_SIZE]
                    regs = []
                    # load into registers r0..r{len(batch_entries)-1}
                    for j, (sign, idx) in enumerate(batch_entries):
                        reg = f"r_{m}_{n}_{batch_id}_{j}"
                        if sign > 0:
                            f.write(f"    float {reg} = X[{base_x + idx}];\n")
                        else:
                            f.write(f"    float {reg} = -X[{base_x + idx}];\n")
                        regs.append(reg)
                    # tree-reduce within batch
                    tmp_count = 0
                    cur_regs = regs.copy()
                    while len(cur_regs) > 1:
                        new_regs = []
                        for i in range(0, len(cur_regs), 2):
                            if i+1 < len(cur_regs):
                                tname = f"t_{m}_{n}_{batch_id}_{tmp_count}"
                                f.write(f"    float {tname} = {cur_regs[i]} + {cur_regs[i+1]};\n")
                                new_regs.append(tname)
                                tmp_count += 1
                            else:
                                new_regs.append(cur_regs[i])
                        cur_regs = new_regs
                    # cur_regs[0] is the partial sum for this batch
                    part = f"partial{m}_{n}_{batch_id}"
                    f.write(f"    float {part} = {cur_regs[0]};\n")
                    partial_sums.append(part)
                # final tree reduction of partial sums
                cur_regs = partial_sums.copy()
                tmp_count = 0
                while len(cur_regs) > 1:
                    new_regs = []
                    for i in range(0, len(cur_regs), 2):
                        if i+1 < len(cur_regs):
                            tname = f"tt_{m}_{n}_{tmp_count}"
                            f.write(f"    float {tname} = {cur_regs[i]} + {cur_regs[i+1]};\n")
                            new_regs.append(tname)
                            tmp_count += 1
                        else:
                            new_regs.append(cur_regs[i])
                    cur_regs = new_regs
                final_sum = cur_regs[0] if cur_regs else "0.0f"
                # add bias and write Y
                f.write(f"    Y[{base_y + n}] = {final_sum} + b[{n}];\n")
        f.write("}\n")

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

log_filename = "codegen/generated.h"
gen_code(X, col_start_pos, col_start_neg, row_index_pos, row_index_neg, b, M, N, K,
         log_filename)

# Y_py = base_csc_flat(X, col_start_pos, col_start_neg,
#                      row_index_pos, row_index_neg, b,
#                      M, N, K)

# print("max abs error:", np.max(np.abs(Y_py - Y_actual)))
# print("allclose?   ", np.allclose(Y_py, Y_actual, atol=1e-6))
