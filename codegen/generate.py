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
    # writes optimized sequence of ops using batching and tree reductions.
    BATCH_SIZE = 8
    with open(log_filename, 'w') as f:
        f.write("inline void GenCSC(float * __restrict__ X, float *__restrict__ b, float *__restrict__ Y) {\n")
        for m in range(M):
            base_x = m * K
            base_y = m * N
            for n in range(N):
                start_p = col_start_pos[n]
                end_p = col_start_pos[n+1]
                pos_indices = row_index_pos[start_p:end_p]
                start_n = col_start_neg[n]
                end_n = col_start_neg[n+1]
                neg_indices = row_index_neg[start_n:end_n]
                entries = [(+1, idx) for idx in pos_indices] + [(-1, idx) for idx in neg_indices]
                num_entries = len(entries)
                num_batches = (num_entries + BATCH_SIZE - 1) // BATCH_SIZE
                partial_sums = []
                batch_id = 0
                while batch_id < num_batches:
                    curr_entries = entries[batch_id * BATCH_SIZE : min((batch_id + 1) * BATCH_SIZE, num_entries)]
                    next_exists  = (batch_id + 1) < num_batches
                    next_entries = entries[(batch_id + 1) * BATCH_SIZE : (batch_id + 2) * BATCH_SIZE] if next_exists else []

                    half_curr = len(curr_entries) // 2
                    half_next = len(next_entries) // 2

                    for j, (sign, idx) in enumerate(curr_entries[:half_curr]):
                        reg = f"r_{m}_{n}_{batch_id}_{j}"
                        if sign > 0:
                            f.write(f"    float {reg} = X[{base_x + idx}];\n")
                        else:
                            f.write(f"    float {reg} = -X[{base_x + idx}];\n")

                    for j, (sign, idx) in enumerate(next_entries[:half_next]):
                        reg = f"r_{m}_{n}_{batch_id+1}_{j}"
                        if sign > 0:
                            f.write(f"    float {reg} = X[{base_x + idx}];\n")
                        else:
                            f.write(f"    float {reg} = -X[{base_x + idx}];\n")

                    for j, (sign, idx) in enumerate(curr_entries[half_curr:], start=half_curr):
                        reg = f"r_{m}_{n}_{batch_id}_{j}"
                        if sign > 0:
                            f.write(f"    float {reg} = X[{base_x + idx}];\n")
                        else:
                            f.write(f"    float {reg} = -X[{base_x + idx}];\n")

                    curr_regs  = [f"r_{m}_{n}_{batch_id}_{j}"    for j in range(len(curr_entries))]
                    tmp_count  = 0
                    while len(curr_regs) > 1:
                        new_regs = []
                        for i in range(0, len(curr_regs), 2):
                            if i + 1 < len(curr_regs):
                                tname = f"t_{m}_{n}_{batch_id}_{tmp_count}"
                                f.write(f"    float {tname} = {curr_regs[i]} + {curr_regs[i+1]};\n")
                                new_regs.append(tname)
                                tmp_count += 1
                            else:
                                new_regs.append(curr_regs[i])
                        curr_regs = new_regs
                    part_cur = f"partial{m}_{n}_{batch_id}"
                    f.write(f"    float {part_cur} = {curr_regs[0] if curr_regs else '0.0f'};\n")
                    partial_sums.append(part_cur)

                    for j, (sign, idx) in enumerate(next_entries[half_next:], start=half_next):
                        reg = f"r_{m}_{n}_{batch_id+1}_{j}"
                        if sign > 0:
                            f.write(f"    float {reg} = X[{base_x + idx}];\n")
                        else:
                            f.write(f"    float {reg} = -X[{base_x + idx}];\n")

                    if next_exists:
                        next_regs = [f"r_{m}_{n}_{batch_id+1}_{j}" for j in range(len(next_entries))]
                        tmp_count = 0
                        while len(next_regs) > 1:
                            new_regs = []
                            for i in range(0, len(next_regs), 2):
                                if i + 1 < len(next_regs):
                                    tname = f"t_{m}_{n}_{batch_id+1}_{tmp_count}"
                                    f.write(f"    float {tname} = {next_regs[i]} + {next_regs[i+1]};\n")
                                    new_regs.append(tname)
                                    tmp_count += 1
                                else:
                                    new_regs.append(next_regs[i])
                            next_regs = new_regs
                        part_next = f"partial{m}_{n}_{batch_id+1}"
                        f.write(f"    float {part_next} = {next_regs[0] if next_regs else '0.0f'};\n")
                        partial_sums.append(part_next)

                    batch_id += 2
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
                f.write(f"    Y[{base_y + n}] = {final_sum} + b[{n}];\n")
        f.write("}\n")

def base_csc_flat(X_flat, col_start_pos, col_start_neg,
                  row_index_pos, row_index_neg, b,
                  M, N, K):
    X = np.asarray(X_flat, dtype=np.float32)
    b = np.asarray(b, dtype=np.float32)

    Y = np.empty(M * N, dtype=np.float32)

    for m in range(M):
        base_x = m * K
        base_y = m * N
        for n in range(N):
            y_val = 0.0

            start_p = col_start_pos[n]
            end_p   = col_start_pos[n+1]
            for idx in range(start_p, end_p):
                row = row_index_pos[idx]
                y_val += X[base_x + row]

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
