#pragma once
#include<Vector>;
#include<iostream>;
using namespace std;

class SparseFormat {
public: 
	vector<int> col_start_pos;
	vector<int> col_start_neg;
	vector<int> row_index_pos;
	vector<int> row_index_neg;

	SparseFormat(int* matrix, int K, int N) {
		int column_start_pos = 0;
		int column_start_neg = 0;
		for (int n = 0; n < N; n++) {
			this->col_start_pos.push_back(column_start_pos);
			this->col_start_neg.push_back(column_start_neg);
			for (int k = 0; k < K; k++) {
				if (matrix[k * N + n] >= 1) {
					column_start_pos++;
					this->row_index_pos.push_back(k);
				}
				else if (matrix[k * N + n] <= -1) {
					column_start_neg++;
					this->row_index_neg.push_back(k);
				}
			}
		}
		this->col_start_pos.push_back(column_start_pos);
		this->col_start_neg.push_back(column_start_neg);
	}
};

template <typename T>
void sparseGEMM(T* X, SparseFormat W, T* b, T* Y, int M, int N, int K) {
	for (int m = 0; m < M; m++) {
		for (int n = 0; n < N; n++) {
			T y = 0;
			for (int k = W.col_start_pos[n]; k < W.col_start_pos[n + 1]; k++) {
				y += X[m * K + W.row_index_pos[k]];
			}
			for (int k = W.col_start_neg[n]; k < W.col_start_neg[n + 1]; k++) {
				y -= X[m * K + W.row_index_neg[k]];
			}
			Y[m * N + n] = y + b[n];
		}
	}
}

template <typename T>
void GEMM(T* X, T* W, T* b, T* Y, int M, int N, int K) {
	for (int m = 0; m < M; m++) {
		for (int n = 0; n < N; n++) {
			T y = 0;
			for (int k = 0; k < K; k++) {
				y += X[m * K + k] * W[k * N + n];
			}
			Y[m * N + n] = y + b[n];
		}
	}
}

template <typename T> 
bool compare_results(T* result, T* groundTruth, int H, int W) {
	for (int h = 0; h < H; h++) {
		for (int w = 0; w < W; w++) {
			int i = h * W + w;
			if (abs(result[i] - groundTruth[i]) > 10e-6) {
				cout << "Error at: H=" << h << ", W=" << w << ", result=" << result[i] << ", groundTruth=" << groundTruth[i] << endl;
				return false;
			}
		}
	}

	return true;
}