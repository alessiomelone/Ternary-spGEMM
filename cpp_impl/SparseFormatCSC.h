#pragma once
#include <vector>

class SparseFormatCSC
{
public:
	std::vector<int> col_start_pos;
	std::vector<int> col_start_neg;
	std::vector<int> row_index_pos;
	std::vector<int> row_index_neg;

	SparseFormatCSC(int *matrix, int K, int N)
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