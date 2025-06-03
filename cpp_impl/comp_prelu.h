#ifndef COMP_H_PRELU
#define COMP_H_PRELU

#include "common.h"
#include <iostream>
#include <arm_neon.h>

#ifdef INSTRUMENTATION_RUN
extern long long flops;
#endif

template <typename T>
void BaseTCSC_PreLU(T *X, const TCSC &W_csc, T *b, T *alpha, T *Y, int M, int N, int K)
{
#ifdef INSTRUMENTATION_RUN
    flops = 0;
    ds_size = W_csc.getDataStructureSize();
#endif
    const int *col_start_pos = W_csc.col_start_pos.data();
    const int *col_start_neg = W_csc.col_start_neg.data();
    const int *row_index_pos = W_csc.row_index_pos.data();
    const int *row_index_neg = W_csc.row_index_neg.data();

    for (int m = 0; m < M; m++)
    {
        for (int n = 0; n < N; n++)
        {
            T y_val = 0;

            // Process positive values
            for (int k = col_start_pos[n]; k < col_start_pos[n + 1]; k++)
            {
                T x_val = X[m * K + row_index_pos[k]];
                y_val += x_val;
#ifdef INSTRUMENTATION_RUN
                flops++;
#endif
            }

            // Process negative values
            for (int k = col_start_neg[n]; k < col_start_neg[n + 1]; k++)
            {
                T x_val = X[m * K + row_index_neg[k]];
                y_val -= x_val;
#ifdef INSTRUMENTATION_RUN
                flops++;
#endif
            }

            // Apply bias
            y_val = y_val + b[n];
#ifdef INSTRUMENTATION_RUN
            flops++;
#endif

            // Apply PReLU activation: f(x) = x if x > 0, f(x) = α * x if x ≤ 0
            if (y_val > 0)
            {
                Y[m * N + n] = y_val;
            }
            else
            {
                Y[m * N + n] = alpha[n] * y_val;
#ifdef INSTRUMENTATION_RUN
                flops++;
#endif
            }
        }
    }
}

template <typename T>
void NeonTCSCHorizontalAdvanced_PReLU(float *X, const VectorTCSC &W_csc, float *b, float *a, float *Y, int M, int N, int K)
{
#ifdef INSTRUMENTATION_RUN
    flops = 0;
    ds_size = W_csc.getDataStructureSize();
#endif
    const int *row_index_pos = W_csc.row_index_pos.data();
    const int *row_index_neg = W_csc.row_index_neg.data();
    const int *cap_every_four = W_csc.cap_every_four.data();

    for (int m = 0; m < M; m++)
    {
        float* X_row_m = X + m * K;
        X_row_m[-1] = 0;

        int cap_idx = 0;
        int k3_end = 0;
        for (int n = 0; n < N; n += 4)
        {
            float32x4_t y_vec0 = vdupq_n_f32(0.0f);
            float32x4_t y_vec1 = vdupq_n_f32(0.0f);
            float32x4_t y_vec2 = vdupq_n_f32(0.0f);
            float32x4_t y_vec3 = vdupq_n_f32(0.0f);

            int cap = cap_every_four[cap_idx++];
            int k0 = k3_end;
            int k0_end = k0 + cap;
            int k1 = k0 + cap;
            int k2 = k1 + cap;
            int k3 = k2 + cap;
            k3_end = k3 + cap;

            for (int i = k0; i < k0_end; i += 4) {
                float32x4_t x_vals000 = { X_row_m[row_index_pos[k0]], X_row_m[row_index_pos[k0+1]], X_row_m[row_index_pos[k0+2]], X_row_m[row_index_pos[k0+3]] };
                float32x4_t x_vals001 = { X_row_m[row_index_neg[k0]], X_row_m[row_index_neg[k0+1]], X_row_m[row_index_neg[k0+2]], X_row_m[row_index_neg[k0+3]] };
                float32x4_t x_vals010 = { X_row_m[row_index_pos[k1]], X_row_m[row_index_pos[k1+1]], X_row_m[row_index_pos[k1+2]], X_row_m[row_index_pos[k1+3]] };
                float32x4_t x_vals011 = { X_row_m[row_index_neg[k1]], X_row_m[row_index_neg[k1+1]], X_row_m[row_index_neg[k1+2]], X_row_m[row_index_neg[k1+3]] };
                float32x4_t x_vals100 = { X_row_m[row_index_pos[k2]], X_row_m[row_index_pos[k2+1]], X_row_m[row_index_pos[k2+2]], X_row_m[row_index_pos[k2+3]] };
                float32x4_t x_vals101 = { X_row_m[row_index_neg[k2]], X_row_m[row_index_neg[k2+1]], X_row_m[row_index_neg[k2+2]], X_row_m[row_index_neg[k2+3]] };
                float32x4_t x_vals110 = { X_row_m[row_index_pos[k3]], X_row_m[row_index_pos[k3+1]], X_row_m[row_index_pos[k3+2]], X_row_m[row_index_pos[k3+3]] };
                float32x4_t x_vals111 = { X_row_m[row_index_neg[k3]], X_row_m[row_index_neg[k3+1]], X_row_m[row_index_neg[k3+2]], X_row_m[row_index_neg[k3+3]] };

                float32x4_t vec_tmp0 = vsubq_f32(x_vals000, x_vals001);
                float32x4_t vec_tmp1 = vsubq_f32(x_vals010, x_vals011);
                float32x4_t vec_tmp2 = vsubq_f32(x_vals100, x_vals101);
                float32x4_t vec_tmp3 = vsubq_f32(x_vals110, x_vals111);

                k0 += 4;
                k1 += 4;
                k2 += 4;
                k3 += 4;
#ifdef INSTRUMENTATION_RUN
                flops+= 32;
#endif
                y_vec0 = vaddq_f32(y_vec0, vec_tmp0);
                y_vec1 = vaddq_f32(y_vec1, vec_tmp1);
                y_vec2 = vaddq_f32(y_vec2, vec_tmp2);
                y_vec3 = vaddq_f32(y_vec3, vec_tmp3);
            }

#ifdef INSTRUMENTATION_RUN
                flops+= 16;
#endif
            // Horizontal add
            float y0 = vaddvq_f32(y_vec0);
            float y1 = vaddvq_f32(y_vec1);
            float y2 = vaddvq_f32(y_vec2);
            float y3 = vaddvq_f32(y_vec3);

#ifdef INSTRUMENTATION_RUN
                flops+= 12;
#endif
            float32x4_t b_vec = vld1q_f32(b + n);
            float32x4_t a_vec = vld1q_f32(a + n);
            float32x4_t y_res = { y0, y1, y2, y3 };
            float32x4_t biased_vec = vaddq_f32(y_res, b_vec);
            float32x4_t neg_path = vmulq_f32(biased_vec, a_vec);
            uint32x4_t mask = vcgtq_f32(biased_vec, vdupq_n_f32(0.0f));
            float32x4_t result_vec = vbslq_f32(mask, biased_vec, neg_path);
            vst1q_f32(Y + m * N + n, result_vec);
            // Store final results
            // Y[m * N + n]     = y0 + b[n];
            // Y[m * N + n + 1] = y1 + b[n + 1];
            // Y[m * N + n + 2] = y2 + b[n + 2];
            // Y[m * N + n + 3] = y3 + b[n + 3];
            // // Code below doesn't add much in performance smh
            // float32x4_t y_res = { y0, y1, y2, y3 };
            // float32x4_t b_vals = vld1q_f32(b + n);
            // float32x4_t store_in_y = vaddq_f32(y_res, b_vals);
            // vst1q_f32(Y + m * N + n, store_in_y);

        }
    }
}

template <typename T>
void NeonTCSCVertical_PReLU(float *X, const VectorTCSC &W_csc, float *b, float *a, float *Y, int M, int N, int K)
{
#ifdef INSTRUMENTATION_RUN
    flops = 0;
    ds_size = W_csc.getDataStructureSize();
#endif
    const int *row_index_pos = W_csc.row_index_pos.data();
    const int *row_index_neg = W_csc.row_index_neg.data();
    const int *cap_every_four = W_csc.cap_every_four.data();

    for (int m = 0; m < M; m++)
    {
        float *X_row_m = X + m * K;
        X_row_m[-1] = 0;

        int cap_idx = 0;
        int k3_end = 0;
        for (int n = 0; n < N; n += 4)
        {
            float32x4_t y_vec0 = vdupq_n_f32(0.0f);
            float32x4_t y_vec1 = vdupq_n_f32(0.0f);

            int cap = cap_every_four[cap_idx++];
            int k0 = k3_end;
            int k0_end = k0 + cap;
            int k1 = k0 + cap;
            int k2 = k1 + cap;
            int k3 = k2 + cap;
            k3_end = k3 + cap;

            for (int i = 0; i < cap; i += 4)
            {
                float32x4_t x_vals000 = {X_row_m[row_index_pos[k0]],
                                         X_row_m[row_index_pos[k1]],
                                         X_row_m[row_index_pos[k2]],
                                         X_row_m[row_index_pos[k3]]};
                float32x4_t x_vals001 = {X_row_m[row_index_neg[k0]],
                                         X_row_m[row_index_neg[k1]],
                                         X_row_m[row_index_neg[k2]],
                                         X_row_m[row_index_neg[k3]]};
                float32x4_t x_vals010 = {X_row_m[row_index_pos[k0 + 1]],
                                         X_row_m[row_index_pos[k1 + 1]],
                                         X_row_m[row_index_pos[k2 + 1]],
                                         X_row_m[row_index_pos[k3 + 1]]};
                float32x4_t x_vals011 = {X_row_m[row_index_neg[k0 + 1]],
                                         X_row_m[row_index_neg[k1 + 1]],
                                         X_row_m[row_index_neg[k2 + 1]],
                                         X_row_m[row_index_neg[k3 + 1]]};
                float32x4_t x_vals100 = {X_row_m[row_index_pos[k0 + 2]],
                                         X_row_m[row_index_pos[k1 + 2]],
                                         X_row_m[row_index_pos[k2 + 2]],
                                         X_row_m[row_index_pos[k3 + 2]]};
                float32x4_t x_vals101 = {X_row_m[row_index_neg[k0 + 2]],
                                         X_row_m[row_index_neg[k1 + 2]],
                                         X_row_m[row_index_neg[k2 + 2]],
                                         X_row_m[row_index_neg[k3 + 2]]};
                float32x4_t x_vals110 = {X_row_m[row_index_pos[k0 + 3]],
                                         X_row_m[row_index_pos[k1 + 3]],
                                         X_row_m[row_index_pos[k2 + 3]],
                                         X_row_m[row_index_pos[k3 + 3]]};
                float32x4_t x_vals111 = {X_row_m[row_index_neg[k0 + 3]],
                                         X_row_m[row_index_neg[k1 + 3]],
                                         X_row_m[row_index_neg[k2 + 3]],
                                         X_row_m[row_index_neg[k3 + 3]]};
                k0 += 4;
                k1 += 4;
                k2 += 4;
                k3 += 4;
#ifdef INSTRUMENTATION_RUN
                flops+= 32;
#endif
                float32x4_t vec_tmp0 = vsubq_f32(x_vals000, x_vals001);
                float32x4_t vec_tmp1 = vsubq_f32(x_vals010, x_vals011);
                float32x4_t vec_tmp2 = vsubq_f32(x_vals100, x_vals101);
                float32x4_t vec_tmp3 = vsubq_f32(x_vals110, x_vals111);
                float32x4_t vec_tmp4 = vaddq_f32(vec_tmp0, vec_tmp1);
                float32x4_t vec_tmp5 = vaddq_f32(vec_tmp2, vec_tmp3);
                y_vec0 = vaddq_f32(y_vec0, vec_tmp4);
                y_vec1 = vaddq_f32(y_vec1, vec_tmp5);
            }

#ifdef INSTRUMENTATION_RUN
                flops+= 12;
#endif
            float32x4_t b_vec = vld1q_f32(b + n);
            float32x4_t a_vec = vld1q_f32(a + n);
            float32x4_t y_res = vaddq_f32(y_vec0, y_vec1);
            float32x4_t biased_vec = vaddq_f32(y_res, b_vec);
            float32x4_t neg_path = vmulq_f32(biased_vec, a_vec);
            uint32x4_t mask = vcgtq_f32(biased_vec, vdupq_n_f32(0.0f));
            float32x4_t result_vec = vbslq_f32(mask, biased_vec, neg_path);
            vst1q_f32(Y + m * N + n, result_vec);

            // float32x4_t b_vals = vld1q_f32(b + n);
            // float32x4_t y_vec_res = vaddq_f32(y_vec0, y_vec1);
            // float32x4_t store_in_y = vaddq_f32(y_vec_res, b_vals);
            // vst1q_f32(Y + m * N + n, store_in_y);
        }
    }
}


#endif // COMP_H_PRELU
