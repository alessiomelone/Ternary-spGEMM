#pragma once
#include <vector>
#include <iostream>

using namespace std;


class RSR
{
public:
<<<<<<< HEAD
<<<<<<< HEAD
    int k;
    int *Bmin_ind2;
    int *Bplu_ind2;
    int perm0_size;
    int ssize;
    int perm_size;
    int us_buf_size_f;
    float *us_buf;
    float *result;
\
    RSR(vector<int> W_raw, int K, int N);
    ~RSR();
};

void MMPlusB(float *X_arg, RSR& rsr, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg);
=======
    pair<vector<vector<int>>, vector<vector<int>>> pre_Bminus;
    pair<vector<vector<int>>, vector<vector<int>>> pre_Bplus;
    vector<vector<int>> bin_k;
    int k;    
    RSR(vector<int> W_raw, int K, int N);
};

void MMPlusB(float *X_arg, const RSR& rsr, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg);
>>>>>>> e9a5eca (Add RSR data structure)
=======
    // pair<vector<vector<int>>, vector<vector<int>>> pre_Bminus;
    // pair<vector<vector<int>>, vector<vector<int>>> pre_Bplus;
    vector<vector<int>> bin_k;
    int k;
    int n;
    int n_preBminus;
    int n_preBplus;
    vector<vector<int>> indices_Bminus;
    vector<vector<int>> indices_Bplus;
    vector<int> segment_sizes_Bminus;
    vector<int> segment_sizes_Bplus;
    int permutations_size;
    float powK;
    vector<vector<vector<float>>> us1;
    vector<vector<vector<float>>> us2;
    RSR(vector<int> W_raw, int K, int N);
};

void MMPlusB(float *X_arg, RSR& rsr, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg);
>>>>>>> b3ff305 (Remove const from RSR API)


