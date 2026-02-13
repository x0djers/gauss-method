#include <math.h>

#include <gauss/solver.h>

GaussSolutions solveLinearSystemByGauss(const MatrixOutcome coeffsMatrix) {
    const GaussMatrixResult gaussMatrix = getGaussMatrix(coeffsMatrix);
    const GaussSolutions gaussSolutions = getGaussSolutions(gaussMatrix);

    return gaussSolutions;
}

MatrixDeterminant calcDeterminantByGauss(const MatrixOutcome coeffsMatrix) {
    MatrixDeterminant determinant = {0, NONE_ERROR};

    const GaussMatrixResult gaussMatrix = getGaussMatrix(coeffsMatrix);
    determinant.errorCode = gaussMatrix.matrix.errorCode;

    if (determinant.errorCode == NONE_ERROR) {
        determinant.determinant = pow(-1, gaussMatrix.stepsCount);
        for (size_t iter = 0; iter < gaussMatrix.matrix.matrix->cols; ++iter) {
            determinant.determinant *= gaussMatrix.matrix.matrix->data[iter][iter];
        }
    }

    return determinant;
}