#include <gauss/solver.h>

GaussSolutions solveLinearSystemByGauss(const MatrixOutcome coeffsMatrix) {
    const GaussMatrixResult gaussMatrix = getGaussMatrix(coeffsMatrix);
    const GaussSolutions gaussSolutions = getGaussSolutions(gaussMatrix);

    return gaussSolutions;
}
