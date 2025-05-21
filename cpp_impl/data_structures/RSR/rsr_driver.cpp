
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include "naive.h"
#include "rsr.h"
#include "rsrpp.h"
#include "utils.h"
#include "rsr_driver.h"

using namespace std;

static pair<vector<vector<int>>,vector<vector<int>>> convertSparseToTwoBinaryMatrices(vector<int> W_raw, int K, int N);
static vector<vector<int>> generateBinaryMatrix2(int k);

static void print_vector_of_ints(vector<int> &a) {
    for (int i = 0; i < a.size(); i++)
        cout << a[i] << ",";
    cout << "EOV" << endl;
}

#if 1
void MMPlusB(float *X_arg, RSR& rsr, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg) {
    int n = rsr.n;
    int perm_size = rsr.permutations_size;
    auto& bin_k = rsr.bin_k;
    auto& seg_size = rsr.seg_size;
    int k = rsr.k;
    auto& indices_Bplus = rsr.indices_Bplus;
    auto& indices_Bminus = rsr.indices_Bminus;
    auto& sz_and_idx = rsr.size_and_index;
    for (int row = 0; row < M_arg; row++) {
        float *v = X_arg + row * K_arg;
        float *res = Y_arg + row * N_arg;
        int index = 0, index1 = 0;
        int sz_and_idx_index = 0;
        for (int i = 0; i < perm_size; i++) {
            for (int j = 0; j < seg_size; j++) {
                float acc = 0;
                int plus_size = sz_and_idx[sz_and_idx_index++];
                for (int p = 0; p < plus_size; p++) {
                    acc += v[sz_and_idx[sz_and_idx_index++]];
                }
                int minus_size = sz_and_idx[sz_and_idx_index++];
                for (int p = 0; p < minus_size; p++) {
                    acc -= v[sz_and_idx[sz_and_idx_index++]];
                }
                if (acc != 0) {
                    for (int l=0; l < k; l++) {
                        res[i*k + l] += acc * bin_k[j][l];
                    }
                }

                // auto& local_indices = indices_Bplus[index++];
                // auto& local_indices_min = indices_Bminus[index1++];

                // int lisz = local_indices.size();
                // int lisz_min = local_indices_min.size();
                // if (lisz + lisz_min != 0) {
                //     float acc = 0.0;
                //     for (int l = 0; l < lisz; l++) {
                //         acc += v[local_indices[l]];
                //     }
                //     for (int l = 0; l < lisz_min; l++) {
                //         acc -= v[local_indices_min[l]];
                //     }
                //     for (int l=0; l < k; l++) {
                //         res[i*k + l] += acc * bin_k[j][l];
                //     }
                // }
            }
        }

        for (int i = 0; i < n; i++) {
            res[i] += B_arg[i];
        }
    }
}
#else
void MMPlusB(float *X_arg, RSR& rsr, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg) {
    int n = rsr.n;
    int perm_size = rsr.permutations_size;
    auto& bin_k = rsr.bin_k;
    auto& seg_size = rsr.seg_size;
    int k = rsr.k;
    auto& indices_Bplus = rsr.indices_Bplus;
    auto& indices_Bminus = rsr.indices_Bminus;
    for (int row = 0; row < M_arg; row++) {
        float *v = X_arg + row * K_arg;
        float *res = Y_arg + row * N_arg;
        int index = 0, index1 = 0;
        for (int i = 0; i < perm_size; i++) {
            for (int j = 0; j < seg_size; j++) {
                auto& local_indices = indices_Bplus[index++];
                auto& local_indices_min = indices_Bminus[index1++];

                int lisz = local_indices.size();
                int lisz_min = local_indices_min.size();
                if (lisz + lisz_min != 0) {
                    float acc = 0.0;
                    for (int l = 0; l < lisz; l++) {
                        acc += v[local_indices[l]];
                    }
                    for (int l = 0; l < lisz_min; l++) {
                        acc -= v[local_indices_min[l]];
                    }
                    for (int l=0; l < k; l++) {
                        res[i*k + l] += acc * bin_k[j][l];
                    }
                }
            }
        }

        for (int i = 0; i < n; i++) {
            res[i] += B_arg[i];
        }
    }
}
#endif
RSR::RSR(vector<int> W_raw, int K, int N) {
    int n = K; // > N ? K : N; // TODO : Experiment with that
    int k = static_cast<int>(ceil(log2(n) - log2(log2(n))));
    auto pr = convertSparseToTwoBinaryMatrices(W_raw, K, N);
    auto pre_Bminus = preprocess(pr.first, k);
    auto pre_Bplus = preprocess(pr.second, k);
    this->bin_k = generateBinaryMatrix2(k);
    this->k = k;
    this->n = pre_Bminus.first[0].size();
    this->seg_size = pre_Bminus.second[0].size();
    this->permutations_size = pre_Bminus.first.size();
    
    auto& permutations = pre_Bminus.first;
    auto& segments = pre_Bminus.second;
    vector<int> sizes_Bminus(segments.size());
    int sum_minus = 0;
    // vector<vector<int>> start_end_idx(permutations.size(), vector<int>(2 * segments.size()));
    for (int i = 0; i < permutations.size(); i++) {
        auto& segment = segments[i];
        auto& permutation = permutations[i];
        int start, end;
        // Each block
        for (int j = 0; j < segment.size(); j++) {
            start = segment[j];
            if (j < segment.size() - 1) {
                end = segment[j + 1];
            } else {
                end = n;
            }
            // Segmented sum
            int id = end - start;
            vector<int> local_indices(id > 0 ? id : 0);
            int k = 0;
            for (int index = start; index < end; index++) {
                local_indices[k++] = permutation[index];
            }
            this->indices_Bminus.push_back(local_indices);
            sizes_Bminus.push_back(k);
            sum_minus += k;
        }
    }
    permutations = pre_Bplus.first;
    segments = pre_Bplus.second;
    int sum_plus = 0;
    vector<int> sizes_Bplus(segments.size());
    // vector<vector<int>> start_end_idx(permutations.size(), vector<int>(2 * segments.size()));
    for (int i = 0; i < permutations.size(); i++) {
        auto& segment = segments[i];
        auto& permutation = permutations[i];
        int start, end;
        // Each block
        for (int j = 0; j < segment.size(); j++) {
            start = segment[j];
            if (j < segment.size() - 1) {
                end = segment[j + 1];
            } else {
                end = n;
            }
            // Segmented sum
            int id = end - start;
            vector<int> local_indices(id > 0 ? id : 0);
            int k = 0;
            for (int index = start; index < end; index++) {
                local_indices[k++] = permutation[index];
            }
            this->indices_Bplus.push_back(local_indices);
            sizes_Bplus.push_back(k);
            sum_plus += k;
        }
    }

    int max_size = sizes_Bminus.size() + sizes_Bplus.size() + sum_plus + sum_minus;
    vector<int> sz_and_idx(max_size);
    int i = 0;
    int idx_minus_ctr = 0, idx_plus_ctr = 0;
    int siz_minus_ctr = 0, siz_plus_ctr = 0;
    while (i < max_size) {
        int sz_Bplus = sizes_Bplus[siz_plus_ctr++];
        auto& idx_Bplus = this->indices_Bplus[idx_plus_ctr++];
        sz_and_idx[i++] = sz_Bplus;
        for (int j = 0; j < sz_Bplus; j++) {
            sz_and_idx[i++] = idx_Bplus[j];
        }
        int sz_Bminus = sizes_Bminus[siz_minus_ctr++];
        auto& idx_Bminus = this->indices_Bminus[idx_minus_ctr++];
        sz_and_idx[i++] = sz_Bminus;
        for (int j = 0; j < sz_Bminus; j++) {
            sz_and_idx[i++] = idx_Bminus[j];
        }

        if (idx_minus_ctr > this->indices_Bminus.size() ||
            idx_plus_ctr > this->indices_Bplus.size()   ||
            siz_plus_ctr > sizes_Bplus.size()           ||
            siz_minus_ctr > sizes_Bminus.size())
        { break; }
    }

    cout << "\n\n [*] Printing sizes_Bplus vector:\n\n" << endl; 
    print_vector_of_ints(sizes_Bplus);
    cout << "\n\n [*] Printing indices_Bplus vectors:\n\n" << endl; 
    for (int i = 0; i < this->indices_Bplus.size(); i++)
    {
        cout << " -- Vector " << endl;
        print_vector_of_ints(this->indices_Bplus[i]);
    }
    
    cout << "\n\n [*] Printing sizes_Bminus vector:\n\n" << endl; 
    print_vector_of_ints(sizes_Bminus);
    cout << "\n\n [*] Printing indices_Bminus vectors:\n\n" << endl; 
    for (int i = 0; i < this->indices_Bminus.size(); i++)
    {
        cout << " -- Vector " << endl;
        print_vector_of_ints(this->indices_Bminus[i]);
    }


    this->size_and_index = sz_and_idx;
    cout << "\n\n [*] Printing size_and_index vector:\n\n" << endl; 
    print_vector_of_ints(sz_and_idx);
    this->powK = pow(2, k);
}

static pair<vector<vector<int>>,vector<vector<int>>> convertSparseToTwoBinaryMatrices(vector<int> W_raw, int K, int N) {
    vector<vector<int>> Bminus(K, vector<int>(N));
    vector<vector<int>> Bplus(K, vector<int>(N));
    for (size_t i = 0; i < K; i++) {
        for (size_t j = 0; j < N; j++) {
            int index = i * N + j;
            if (W_raw[index] == 1) {
                Bplus[i][j] = 1;
                Bminus[i][j] = 0;
            }
            else if (W_raw[index] == -1) {
                Bplus[i][j] = 0;
                Bminus[i][j] = 1;
            }
            else {
                Bplus[i][j] = 0;
                Bminus[i][j] = 0;
            }
        }
    }
    return make_pair(Bminus, Bplus);
}

static vector<vector<int>> generateBinaryMatrix2(int k) {
    int rows = pow(2, k);  // 2^k rows
    vector<vector<int>> matrix(rows, vector<int>(k, 0));  // Initialize matrix with 0s

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < k; ++j) {
            // Generate the binary value for each position
            matrix[i][k - j - 1] = (i >> j) & 1;  // Extract the j-th bit from i
        }
    }

    return matrix;
}