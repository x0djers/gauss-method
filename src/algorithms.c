#include <stdio.h>
#include <string.h>
#include <tgmath.h>

#include <gauss/algorithms.h>

size_t findRowWithMaxAbsInColumn(const MatrixOutcome matrix,
                                 const size_t startRow,
                                 const size_t targetColumn) {
    size_t targetRow = startRow;
    MATRIX_TYPE currentMax = fabs(matrix.matrix->data[targetRow][targetColumn]);

    for (size_t rowIndex = startRow; rowIndex < matrix.matrix->rows; rowIndex++) {
        const MATRIX_TYPE currentAbsNumber = fabs(
            matrix.matrix->data[rowIndex][targetColumn]
            );
        if (currentAbsNumber > currentMax) {
            currentMax = currentAbsNumber;
            targetRow = rowIndex;
        }
    }

    return targetRow;
}

bool swapMatrixRows(const MatrixOutcome matrix,
                    const size_t A,
                    const size_t B) {
    bool isSwapped = false;

    const size_t rowSize = sizeof(MATRIX_TYPE) * matrix.matrix->cols;

    MATRIX_TYPE* buffer = malloc(rowSize);

    if (buffer != NULL &&
        A < matrix.matrix->rows &&
        B < matrix.matrix->rows) {
        if (A != B) {
            memcpy(buffer, matrix.matrix->data[B], rowSize);
            memcpy(matrix.matrix->data[B], matrix.matrix->data[A], rowSize);
            memcpy(matrix.matrix->data[A], buffer, rowSize);
        }

        isSwapped = true;
    }

    free(buffer);

    return isSwapped;
}

GaussMatrixResult transformMatrixByDirectAlgorithm(const MatrixOutcome matrix) {
    GaussMatrixResult result = {matrix, 0};

    size_t currentStep = 0;

    while (matrix.errorCode == NONE_ERROR &&
           currentStep < matrix.matrix->rows - 1) {
        swapMatrixRows(matrix,
                       currentStep,
                       findRowWithMaxAbsInColumn(matrix,
                                                    currentStep,
                                                    currentStep));

        const MATRIX_TYPE supportElement =
            matrix.matrix->data[currentStep][currentStep];

        for (size_t rowIter = currentStep + 1;
             rowIter < matrix.matrix->rows;
             rowIter++) {
            const MATRIX_TYPE factor =
                matrix.matrix->data[rowIter][currentStep] /
                supportElement;


            for (size_t colIter =  currentStep;
                 colIter < matrix.matrix->cols;
                 colIter++) {
                matrix.matrix->data[rowIter][colIter] -=
                    matrix.matrix->data[currentStep][colIter] * factor;
            }
        }

        currentStep++;
    }

    result.stepsCount = currentStep;

    return result;
}

GaussMatrixResult getGaussMatrix(const MatrixOutcome inputMatrix) {
    GaussMatrixResult result = {NULL, NONE_ERROR, 0};
    if (!inputMatrix.matrix) {
        result.matrix.errorCode = NULL_POINTER_ERROR;
    } else {
        const MatrixOutcome matrixCopy = getMatrixCopy(inputMatrix);
        result = transformMatrixByDirectAlgorithm(matrixCopy);
    }

    return result;
}

GaussSolutions getSolutionsByBackAlgorithm(const MatrixOutcome matrix) {
    GaussSolutions solutions = {NULL, 0};

    solutions.values = malloc((matrix.matrix->cols - 1) * sizeof(MATRIX_TYPE));

    if (solutions.values) {
        size_t currentRow = matrix.matrix->rows - 1;
        size_t currentColumn = matrix.matrix->cols - 2;

        const size_t columnsCount = matrix.matrix->cols;

        while ((int)currentRow >= 0) {
            for (size_t colIter = currentColumn + 1;
                 colIter < columnsCount - 1;
                 colIter++) {
                matrix.matrix->data[currentRow][columnsCount - 1] -=
                    matrix.matrix->data[currentRow][colIter];
            }

            const MATRIX_TYPE solution =
                matrix.matrix->data[currentRow][columnsCount - 1] /
                matrix.matrix->data[currentRow][currentColumn];

            solutions.values[currentColumn] = solution;
            ++solutions.count;

            for (size_t rowIter = currentRow - 1;
                 (int)rowIter >= 0;
                 rowIter--) {
                    matrix.matrix->data[rowIter][currentColumn] *= solution;
                 }

            currentColumn--;
            currentRow--;
        }
    }

    return solutions;
}

GaussSolutions getGaussSolutions(const GaussMatrixResult gaussMatrix) {
    GaussSolutions solutions = {NULL, 0};

    if (gaussMatrix.matrix.errorCode == NONE_ERROR && gaussMatrix.matrix.matrix) {
        const MatrixOutcome matrixCopy = getMatrixCopy(gaussMatrix.matrix);
        solutions = getSolutionsByBackAlgorithm(matrixCopy);
    }

    return solutions;
}
