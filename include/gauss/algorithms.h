#ifndef ALGORITHMS_H
#define ALGORITHMS_H

#include <matrix/matrix.h>

typedef struct {
    MatrixOutcome matrix;
    size_t stepsCount;
} GaussMatrixResult;

typedef struct {
    MATRIX_TYPE* values;
    size_t count;
} GaussSolutions;


#endif
