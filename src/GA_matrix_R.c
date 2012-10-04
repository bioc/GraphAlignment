/* ----------------------------------------------------------------------------
 * R package for graph alignment
 * ----------------------------------------------------------------------------
 *
 * Author: Joern P. Meier <mail@ionflux.org>
 * 
 * The package can be used freely for non-commercial purposes. If you use this 
 * package, the appropriate paper to cite is J. Berg and M. Laessig, 
 * "Cross-species analysis of biological networks by Bayesian alignment", 
 * PNAS 103 (29), 10967-10972 (2006)
 * 
 * This software is made available in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
 * or FITNESS FOR A PARTICULAR PURPOSE.
 * 
 * This software contains code for solving linear assignment problems which was 
 * written by Roy Jonker, MagicLogic Optimization Inc.. Please note that this 
 * code is copyrighted, (c) 2003 MagicLogic Systems Inc., Canada and may be 
 * used for non-commercial purposes only. See 
 * http://www.magiclogic.com/assignment.html for the latest version of the LAP 
 * code and details on licensing.
 *
 * ----------------------------------------------------------------------------
 * R utility functions for matrix types.
 * ----------------------------------------------------------------------------
 */

/** \file GA_matrix_R.c
 * \brief R utility functions for matrix types (implementation).
 */

#include "GA_alloc.h"
#include "GA_message.h"
#include "GA_matrix_R.h"

GAMatrixInt* GA_matrix_from_R_int(SEXP robj)
{
    PROTECT(robj);
    SEXPTYPE matrixType = TYPEOF(robj);
    if ((matrixType != INTSXP)
        && (matrixType != REALSXP))
    {
        char* message = GA_alloc(256, sizeof(char));
        snprintf(message, 256, 
            "[GA_matrix_from_R_int] Input is not a matrix of real or "
            "integer values (actual type: %i).", matrixType);
        GA_msg()(message, GA_MSG_ERROR);
        GA_free(message);
        UNPROTECT(1);
        return 0;
    }
    int numDims = LENGTH(GET_DIM(robj));
    if (numDims != 2)
    {
        GA_msg()("[GA_matrix_from_R_int] Input is not a "
            "two-dimensional matrix.", GA_MSG_ERROR);
        UNPROTECT(1);
        return 0;
    }
    int* dims = INTEGER(coerceVector(GET_DIM(robj), INTSXP));
    GAMatrixInt* matrix = GA_matrix_create_int(dims[0], dims[1]);
    if (matrix == 0)
    {
        UNPROTECT(1);
        return 0;
    }
    int* inputRaw = INTEGER(coerceVector(robj, INTSXP));
    int i;
    int j;
    for (i = 0; i < dims[0]; i++)
        for (j = 0; j < dims[1]; j++)
            matrix->elts[i][j] = inputRaw[j * dims[0] + i];
    UNPROTECT(1);
    return matrix;
}

SEXP GA_matrix_to_R_int(GAMatrixInt* matrix)
{
    SEXP result;
    PROTECT(result = allocMatrix(INTSXP, matrix->rows, matrix->cols));
    int* resultRaw = INTEGER(coerceVector(result, INTSXP));
    int i;
    int j;
    for (i = 0; i < matrix->rows; i++)
        for (j = 0; j < matrix->cols; j++)
            resultRaw[j * matrix->rows + i] = matrix->elts[i][j];
    UNPROTECT(1);
    return result;
}

GAMatrixReal* GA_matrix_from_R_real(SEXP robj)
{
    PROTECT(robj);
    SEXPTYPE matrixType = TYPEOF(robj);
    if ((matrixType != INTSXP)
        && (matrixType != REALSXP))
    {
        char* message = GA_alloc(256, sizeof(char));
        snprintf(message, 256, 
            "[GA_matrix_from_R_real] Input is not a matrix of real or "
            "integer values (actual type: %i).", matrixType);
        GA_msg()(message, GA_MSG_ERROR);
        GA_free(message);
        UNPROTECT(1);
        return 0;
    }
    int numDims = LENGTH(GET_DIM(robj));
    if (numDims != 2)
    {
        GA_msg()("[GA_matrix_from_R_real] Input is not a "
            "two-dimensional matrix.", GA_MSG_ERROR);
        UNPROTECT(1);
        return 0;
    }
    int* dims = INTEGER(coerceVector(GET_DIM(robj), INTSXP));
    GAMatrixReal* matrix = GA_matrix_create_real(dims[0], dims[1]);
    if (matrix == 0)
    {
        UNPROTECT(1);
        return 0;
    }
    double* inputRaw = REAL(coerceVector(robj, REALSXP));
    int i;
    int j;
    for (i = 0; i < dims[0]; i++)
        for (j = 0; j < dims[1]; j++)
            matrix->elts[i][j] = inputRaw[j * dims[0] + i];
    UNPROTECT(1);
    return matrix;
}

SEXP GA_matrix_to_R_real(GAMatrixReal* matrix)
{
    SEXP result;
    PROTECT(result = allocMatrix(REALSXP, matrix->rows, matrix->cols));
    double* resultRaw = REAL(coerceVector(result, REALSXP));
    int i;
    int j;
    for (i = 0; i < matrix->rows; i++)
        for (j = 0; j < matrix->cols; j++)
            resultRaw[j * matrix->rows + i] = matrix->elts[i][j];
    UNPROTECT(1);
    return result;
}
