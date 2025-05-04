#include <list>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

#include "perf.h"
#include "common.h"

using namespace std;

#define NR 32
#define CYCLES_REQUIRED 1e8
#define REP 10
#define EPS (1e-3)

void register_functions();
double perf_test(comp_func f, int M, int K, int N, int nonZero);
void add_function(comp_func f, string name);

vector<comp_func> userFuncs;
vector<string> funcNames;
int numFuncs = 0;

/*
 * Registers a user function to be tested by the driver program. Registers a
 * string description of the function as well
 */
void add_function(comp_func f, string name)
{
    userFuncs.push_back(f);
    funcNames.emplace_back(name);
    numFuncs++;
}

int main(int argc, char **argv)
{
    cout << "Starting program. ";
    double perf;
    int i;

    int M = 0, K = 0, N = 0, nonZero = 0;

    // Check all arguments passed
    if (argc < 9)
    {
        fprintf(stderr, "Usage: %s -M <int> -K <int> -N <int> -s <int>\n", argv[0]);
        return 1;
    }

    M = atoi(argv[2]);
    K = atoi(argv[4]);
    N = atoi(argv[6]);
    nonZero = atoi(argv[8]);

    // Check for valid parameters
    if (M <= 0 || K <= 0 || N <= 0 || nonZero <= 0)
    {
        fprintf(stderr, "ERROR: All dimensions must be positive integers.\n");
        fprintf(stderr, "Usage: %s -M <int> -K <int> -N <int> -s <int>\n", argv[0]);
        return 1;
    }

    register_functions();

    if (numFuncs == 0)
    {
        cout << endl;
        cout << "No functions registered - nothing for driver to do" << endl;
        cout << "Register functions by calling register_func(f, name)" << endl;
        cout << "in register_funcs()" << endl;

        return 0;
    }
    cout << numFuncs << " functions registered." << endl;

    vector<float> X = initX<float>(M * K, 512);
    vector<int> W = generateSparseMatrix<int>(K, N, nonZero, false);
    vector<float> W_FP32(W.begin(), W.end());
    vector<float> B(N, 2);
    vector<float> Y(M * N, 0);
    vector<float> refY(M * N, 0);
    SparseFormat sf = SparseFormat(W.data(), K, N);
    GEMM(X.data(), W_FP32.data(), B.data(), refY.data(), M, N, K);

    for (i = 0; i < numFuncs; i++)
    {
        fill(Y.begin(), Y.end(), 0);
        comp_func func = userFuncs[i];
        func(X.data(),
             sf.col_start_pos.data(), sf.col_start_neg.data(),
             sf.row_index_pos.data(), sf.row_index_neg.data(),
             B.data(), Y.data(), M, N, K);
        if (compare_results(Y.data(), refY.data(), M, N))
        {
            cout << "Test case passed!" << endl;
        }
        else
        {
            cout << "Test case failed!" << endl;
        }
    }

    // destroy all the matrices
    X.clear();
    W.clear();
    W_FP32.clear();
    B.clear();
    Y.clear();
    refY.clear();

    for (i = 0; i < numFuncs; i++)
    {
        perf = perf_test(userFuncs[i], M, K, N, nonZero);
        cout << endl
             << "Running: " << funcNames[i] << endl;
        cout << perf << " cycles" << endl;
    }

    return 0;
}