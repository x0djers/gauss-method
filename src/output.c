#include <stdio.h>

#include <gauss/output.h>
#include <matrix/output.h>
#include <gauss/algorithms.h>

void printGaussMatrix(const GaussMatrixResult gaussMatrix) {
    if (gaussMatrix.matrix.matrix) {
        printMatrix(gaussMatrix.matrix, outputToStd, stdout);
    }
}

void printGaussSolutions(const GaussSolutions gaussSolutions) {
    if (gaussSolutions.values) {
        for (size_t iter = 0; iter < gaussSolutions.count; iter++) {
            printf("x%lu = %.5f\n", iter, gaussSolutions.values[iter]);
        }
    }
}