#pragma once
#include <vector>
#include <iostream>

using namespace std;


class RSR
{
public:
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


