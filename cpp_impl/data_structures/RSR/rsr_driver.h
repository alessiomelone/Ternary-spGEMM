#pragma once
#include <vector>
#include <iostream>

using namespace std;


class RSR
{
public:
    pair<vector<vector<int>>, vector<vector<int>>> pre_Bminus;
    pair<vector<vector<int>>, vector<vector<int>>> pre_Bplus;
    vector<vector<int>> bin_k;
    int k;    
    RSR(vector<int> W_raw, int K, int N);
};

void MMPlusB(float *X_arg, const RSR& rsr, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg);


