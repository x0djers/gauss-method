#include <math.h>
#include <stdio.h>

#include <gauss/solver.h>
#include <gauss/output.h>

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

MatrixOutcome mergeMatrices(const MatrixOutcome A, const MatrixOutcome B) {
    MatrixOutcome result = {NULL, NONE_ERROR};
    if (A.matrix->rows == B.matrix->rows) {
        const size_t colsCount = A.matrix->cols + B.matrix->cols;
        result  = createMatrix(A.matrix->rows,colsCount);
        for (size_t rowIter = 0; rowIter < A.matrix->rows; rowIter++) {
            for (size_t colIter = 0; colIter < colsCount; colIter++) {
                if (colIter <  A.matrix->cols) {
                    result.matrix->data[rowIter][colIter] =
                        A.matrix->data[rowIter][colIter];
                } else {
                    result.matrix->data[rowIter][colIter] =
                        B.matrix->data[rowIter][colIter - A.matrix->cols];
                }
            }
        }
    } else {
        result.errorCode = SIZE_MISMATCH_ERROR;
    }

    return result;
}

MatrixOutcome getIdentityMatrix(const size_t rowsCount) {
    const MatrixOutcome identityMatrix = createMatrix(rowsCount, rowsCount);
    for (size_t rowIter = 0; rowIter < rowsCount; rowIter++) {
        for (size_t colIter = 0; colIter < rowsCount; colIter++) {
            identityMatrix.matrix->data[rowIter][colIter] =
                rowIter == colIter ? 1 : 0;
        }
    }

    return identityMatrix;
}

MatrixOutcome getColumnByIndex(const MatrixOutcome matrix, const size_t index) {
    const MatrixOutcome column = createMatrix(matrix.matrix->rows, 1);
    for (size_t rowIter = 0; rowIter < matrix.matrix->rows; rowIter++) {
        column.matrix->data[rowIter][0] = matrix.matrix->data[rowIter][index];
    }

    return column;
}

void addSolutionToInverseMatrix(const MatrixOutcome inverseMatrix,
                                const GaussSolutions solution,
                                const size_t index) {
    for (size_t rowIter = 0;
         rowIter < inverseMatrix.matrix->rows;
         rowIter++) {
        inverseMatrix.matrix->data[rowIter][index] = solution.values[rowIter];
    }
}

MatrixOutcome getInverseMatrixByGauss(const MatrixOutcome coeffsMatrix) {
    MatrixOutcome inverseMatrix = {NULL, NONE_ERROR};
    if (coeffsMatrix.matrix) {
        inverseMatrix = createMatrix(coeffsMatrix.matrix->rows,
                                     coeffsMatrix.matrix->cols);
        MatrixOutcome identityMatrix =
            getIdentityMatrix(coeffsMatrix.matrix->rows);
        MatrixOutcome supplementedMatrix =
            mergeMatrices(coeffsMatrix,identityMatrix);

        GaussMatrixResult coeffsGaussMatrix = getGaussMatrix(coeffsMatrix);
        GaussMatrixResult supplGaussMatrix =
            getGaussMatrix(supplementedMatrix);

        for (size_t colIter = coeffsMatrix.matrix->cols;
             colIter < supplementedMatrix.matrix->cols;
             colIter++) {
            MatrixOutcome solColumn = getColumnByIndex(supplGaussMatrix.matrix,
                                                       colIter);
            GaussMatrixResult mergedMatrix = {
                .matrix = mergeMatrices(coeffsGaussMatrix.matrix,
                                        solColumn),
                .stepsCount = 0
            };

            GaussSolutions solution = getGaussSolutions(mergedMatrix);

            addSolutionToInverseMatrix(inverseMatrix,
                                       solution,
                                       colIter - coeffsMatrix.matrix->cols);

            freeMatrixOutcome(&solColumn);
            freeMatrixOutcome(&mergedMatrix.matrix);

            free(solution.values);
            solution.values = NULL;
        }

        freeMatrixOutcome(&identityMatrix);
        freeMatrixOutcome(&supplementedMatrix);

        freeMatrixOutcome(&coeffsGaussMatrix.matrix);
        freeMatrixOutcome(&supplGaussMatrix.matrix);
    } else {
        inverseMatrix.errorCode = NULL_POINTER_ERROR;
    }

    return inverseMatrix;
}