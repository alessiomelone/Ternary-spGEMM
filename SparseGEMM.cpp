#include "SparseGEMM.h"

vector<int> generateSparseMatrix(int H, int W, int nonZero, bool uniformDistribution) {
    vector<int> y = vector<int>(H * W, 0);
    if (uniformDistribution) {
        for (int h = 0; h < H; h++) {
            for (int w = 0; w < W; w += nonZero * 2) {
                // Assign +1, -1 to each 2 x nonZero slots
                int randomA = rand() % nonZero * 2;
                int randomB = rand() % nonZero * 2;
                y[w + randomA] = 1;
                while (randomA==randomB) {
                    randomB = rand() % nonZero * 2;
                }
                y[w + randomB] = 1;
            }
        }
    }
    else {
        for (int h = 0; h < H; h++) {
            // Assign +1 to W / nonZero / 2 places
            int count = 0;
            int limit = (W / nonZero) / 2 - 1;
            while (count < limit) {
                int randomA = rand() % W;
                if (y[randomA] == 0) {
                    y[randomA] = 1;
                    count++;
                }
            }

            // Assign -1 to W / nonZero / 2 places
            count = 0;
            while (count < (W / nonZero / 2 - 1)) {
                int randomA = rand() % W;
                if (y[randomA] == 0) {
                    y[randomA] = -1;
                    count++;
                }
            }
        }
    }

    return y;
}

int main() {
    int TEST_CASES = 8;
    int M[] = { 1, 16, 64, 256, 1000, 4000, 16000, 64000 };
    int K[] = {  512, 1024, 2048,  4096, 2048, 4096, 8192, 16384 };
    int N[] = { 2048, 4096, 8192, 16384,  512, 1024, 2048,  4096 };
    int nonZero[] = { 2, 4, 8, 16 }; // 1/2, 1/4, 1/8, 1/16 non-zero values

    vector<int> X = generateSparseMatrix(M[7], K[7], 4, true);
    cout << "X initialized."<< endl;
    vector<int> W = generateSparseMatrix(K[7], N[7], nonZero[3], true);
    cout << "W initialized." << endl;

    for (int m = 0; m < TEST_CASES; m++) {
        for (int n = 0; n < TEST_CASES; n++) {
            cout << "M=" << M[m] << ", N=" << N[n] << ", K=" << K[n] << endl;
            vector<int> Y = vector<int>(M[m] * N[n], 0);
            vector<int> B = vector<int>(N[n], 2);
            vector<int> refY = vector<int>(M[m] * N[n], 0);
            sparseGEMM(X.data(), SparseFormat(W.data(), K[n], N[n]), B.data(), Y.data(), M[m], N[n], K[n]);
            GEMM(X.data(), W.data(), B.data(), refY.data(), M[m], N[n], K[n]);
            if (compare_results(Y.data(), refY.data(), M[m], N[n])) {
                cout << "Test case passed!" << endl;
            }
        }
    }
}