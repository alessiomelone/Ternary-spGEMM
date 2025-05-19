#pragma once
#include <vector>
#include <iostream> // Used by some template functions

#include <cstdlib>
#include <chrono>
#include <random>
#include <cstdint> // For int8_t if you define new structs using it

using namespace std; // Avoid in headers generally, but keeping for consistency with original

class SparseFormat
{
public:
	vector<int> col_start_pos;
	vector<int> col_start_neg;
	vector<int> row_index_pos;
	vector<int> row_index_neg;

	SparseFormat(int *matrix, int K, int N)
	{
		// ... (existing constructor logic is fine) ...
		int column_start_pos_val = 0; // Renamed for clarity
		int column_start_neg_val = 0; // Renamed for clarity
		for (int n_idx = 0; n_idx < N; n_idx++) // Renamed n
		{
			this->col_start_pos.push_back(column_start_pos_val);
			this->col_start_neg.push_back(column_start_neg_val);
			for (int k_idx = 0; k_idx < K; k_idx++) // Renamed k
			{
				if (matrix[k_idx * N + n_idx] >= 1)
				{
					column_start_pos_val++;
					this->row_index_pos.push_back(k_idx);
				}
				else if (matrix[k_idx * N + n_idx] <= -1)
				{
					column_start_neg_val++;
					this->row_index_neg.push_back(k_idx);
				}
			}
		}
		this->col_start_pos.push_back(column_start_pos_val);
		this->col_start_neg.push_back(column_start_neg_val);
	}
};

// Example of where you would define a new, "radically different" data structure
/*
struct CustomMixedTypeFormat {
    std::vector<int> main_indices;
    std::vector<int8_t> value_types;
    int K_dim; // Store necessary dimensions if needed by the algorithms
    int N_dim;

    // Constructor would take appropriate data to populate these vectors
    CustomMixedTypeFormat(int K, int N, const std::vector<int>& raw_indices, const std::vector<int8_t>& raw_types)
        : main_indices(raw_indices), value_types(raw_types), K_dim(K), N_dim(N) {}

    // Or a constructor that processes raw data like SparseFormat does
    CustomMixedTypeFormat(int* some_custom_matrix_format, int K, int N) : K_dim(K), N_dim(N) {
        // Logic to parse some_custom_matrix_format and populate main_indices and value_types
    }
};
*/

template <typename T>
vector<T> initX(int LEN, int Range)
{
	vector<T> X(LEN, 0);
	mt19937 generator(static_cast<unsigned int>(time(0)));
	uniform_int_distribution<int> range(-Range, Range);
	for (int i = 0; i < LEN; i++)
	{
		X[i] = range(generator);
	}
	return X;
};

template <typename T>
vector<T> generateSparseMatrix(int H, int W, int nonZero, bool uniformDistribution)
{
	vector<T> y = vector<T>(H * W, 0);
	if (uniformDistribution)
	{
		for (int h = 0; h < H; h++)
		{
			for (int w = 0; w < W; w += nonZero * 2)
			{
				// Assign +1, -1 to each 2 x nonZero slots
				int randomA = rand() % nonZero * 2;
				int randomB = rand() % nonZero * 2;
				y[h * W + w + randomA] = 1;
				while (randomA == randomB)
				{
					randomB = rand() % nonZero * 2;
				}
				y[h * W + w + randomB] = -1;
			}
		}
	}
	else
	{
		mt19937 generator(static_cast<unsigned int>(time(0)));
		uniform_int_distribution<int> range(0, W - 1);
		uniform_int_distribution<int> variRange(0, int(W / nonZero / 20 + 1)); // The variation among different columns
		for (int h = 0; h < H; h++)
		{
			int posVari = variRange(generator);
			int limitPos = (W / nonZero) / 2 + posVari;
			int limitNeg = (W / nonZero) / 2 - posVari;

			// Assign +1 to W / nonZero / 2 places
			int count = 0;
			while (count < limitPos)
			{
				int randomA = range(generator);
				if (y[h * W + randomA] == 0)
				{
					y[h * W + randomA] = 1;
					count++;
				}
			}

			// Assign -1 to W / nonZero / 2 places
			count = 0;
			while (count < limitNeg)
			{
				int randomA = range(generator);
				if (y[h * W + randomA] == 0)
				{
					y[h * W + randomA] = -1;
					count++;
				}
			}
		}
	}

	return y;
}

template <typename T>
void GEMM(T *X, T *W, T *b, T *Y, int M, int N, int K)
{
#pragma omp parallel for
	for (int m = 0; m < M; m++)
	{
		for (int n = 0; n < N; n++)
		{
			T y = 0;
			for (int k = 0; k < K; k++)
			{
				y += X[m * K + k] * W[k * N + n];
			}
			Y[m * N + n] = y + b[n];
		}
	}
}

template <typename T>
bool compare_results(T *result, T *groundTruth, int H, int W)
{
	for (int h = 0; h < H; h++)
	{
		for (int w = 0; w < W; w++)
		{
			int i = h * W + w;
			if (abs(result[i] - groundTruth[i]) > 10e-6)
			{
				cout << "Error at: H=" << h << ", W=" << w << ", result=" << result[i] << ", groundTruth=" << groundTruth[i] << endl;
				return false;
			}
		}
	}

	return true;
}