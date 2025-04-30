#include"SparseGEMM.h"

int main() {
    // 10 Basic Test Cases: Sparsity = 50% (nonZeroFactor = 2), general sparsity
    vector<tuple<int, int, int>> basicTestCases = {
    {   1,  512,  2048},
    {   1, 1024,  4096},
    {   1, 2048,  8192},
    {   1, 4096, 16384},
    { 256,  512,  2048},
    { 256, 1024,  4096},
    { 256, 2048,  8192},
    { 256, 4096, 16384},
    {1000, 4096, 16384},
    {4000, 4096, 16384}
    };

    for (const auto& [M_ROW, K_LEN, N_COL] : basicTestCases) {
        cout << "M=" << M_ROW << ", K=" << K_LEN << ", N=" << N_COL << ", Sparsity = 50%, General Sparsity." << endl;
        vector<float> X = initX<float>(M_ROW * K_LEN, 512);
        vector<int> W = generateSparseMatrix<int>(K_LEN, N_COL, 2, false);
        vector<float> W_FP32(W.begin(), W.end());
        vector<float> Y(M_ROW * N_COL, 0);
        vector<float> B(N_COL, 2);
        vector<float> refY(M_ROW * N_COL, 0);
        SparseFormat sf = SparseFormat(W.data(), K_LEN, N_COL);
        sparseGEMM(X.data(), sf.col_start_pos.data(), sf.col_start_neg.data(), sf.row_index_pos.data(), sf.row_index_neg.data(), B.data(), Y.data(), M_ROW, N_COL, K_LEN);
        GEMM(X.data(), W_FP32.data(), B.data(), refY.data(), M_ROW, N_COL, K_LEN);
        if (compare_results(Y.data(), refY.data(), M_ROW, N_COL)) {
            cout << "Test case passed!" << endl;
        }
    }

    // 8x8x4=256 full test cases
    int TEST_CASES = 8;
    int M[] = { 1, 16, 64, 256, 1000, 4000, 16000, 64000 };
    int K[] = { 512, 1024, 2048,  4096, 2048, 4096, 8192, 16384 };
    int N[] = { 2048, 4096, 8192, 16384,  512, 1024, 2048,  4096 };
    int nonZero[] = { 2, 4, 8, 16 }; // 1/2, 1/4, 1/8, 1/16 non-zero values. For example, nonZero = 1/4 means 25% of total weights are +1/-1, equal to 75% sparsity

    // Below are full test cases, you can run those test cases in the final testing stage.
    // Remove those test case that has too long runtime (e.g., 1 function call > 1 minute)
    for (int m = 0; m < TEST_CASES; m++) {
        for (int n = 0; n < TEST_CASES; n++) {
            for (int s = 0; s < TEST_CASES / 2; s++) {
                cout << "M=" << M[m] << ", N=" << N[n] << ", K=" << K[n] << ", nonZeroFactor=" << nonZero[s] << endl;
                vector<float> X = initX<float>(M[m] * K[n], 512);
                vector<int> W = generateSparseMatrix<int>(K[n], N[n], nonZero[s], false);
                vector<float> W_FP32(W.begin(), W.end());
                vector<float> Y(M[m] * N[n], 0);
                vector<float> B(N[n], 2);
                vector<float> refY(M[m] * N[n], 0);
                SparseFormat sf = SparseFormat(W.data(), K[n], N[n]);
                sparseGEMM(X.data(), sf.col_start_pos.data(), sf.col_start_neg.data(), sf.row_index_pos.data(), sf.row_index_neg.data(), B.data(), Y.data(), M[m], N[n], K[n]);
                GEMM(X.data(), W_FP32.data(), B.data(), refY.data(), M[m], N[n], K[n]);
                if (compare_results(Y.data(), refY.data(), M[m], N[n])) {
                    cout << "Test case passed!" << endl;
                }
            }            
        }
    }
    return 0;
}