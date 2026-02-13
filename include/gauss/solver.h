#ifndef GAUSS_SOLVER_H
#define GAUSS_SOLVER_H

#include <gauss/algorithms.h>

GaussSolutions solveLinearSystemByGauss(MatrixOutcome coeffsMatrix);

MatrixDeterminant calcDeterminantByGauss(MatrixOutcome coeffsMatrix);

MatrixOutcome getInverseMatrixByGauss(MatrixOutcome coeffsMatrix);

#endif
