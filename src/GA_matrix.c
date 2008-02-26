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
 * Matrix types.
 * ----------------------------------------------------------------------------
 */

/** \file GA_matrix.c
 * \brief Matrix types (implementation).
 */

#include <stdlib.h>
#include <stdio.h>
#include "GA_alloc.h"
#include "GA_message.h"
#include "GA_matrix.h"

GAMatrixInt* GA_matrix_create_int(int rows, int cols)
{
    GAMatrixInt* matrix = (GAMatrixInt*)GA_alloc(1, sizeof(GAMatrixInt));
    if (matrix == 0)
    {
        GA_msg()("[GA_matrix_create_int] "
            "Could not allocate matrix.", GA_MSG_ERROR);
        return 0;
    }
    matrix->rows = rows;
    matrix->cols = cols;
    matrix->refs = 1;
    matrix->elts = (int**)GA_alloc(rows, sizeof(int*));
    if (matrix->elts == 0)
    {
        GA_msg()("[GA_matrix_create_int] "
            "Could not allocate matrix rows.", GA_MSG_ERROR);
        GA_free((char*)matrix);
        return 0;
    }
    int i;
    for (i = 0; i < matrix->rows; i++)
    {
        matrix->elts[i] = (int*)GA_alloc(cols, sizeof(int));
        if (matrix->elts[i] == 0)
        {
            char* message = GA_alloc(256, sizeof(char));
            snprintf(message, 256, "[GA_matrix_create_int] "
                "Could not allocate matrix row %i.", i);
            GA_msg()(message, GA_MSG_ERROR);
            GA_free(message);
            int j;
            for (j = 0; j < i; j++)
                GA_free((char*)matrix->elts[j]);
            GA_free((char*)matrix->elts);
            GA_free((char*)matrix);
            return 0;
        }
    }
    return matrix;
}

GAMatrixInt* GA_matrix_create_square_int(int size)
{
    return GA_matrix_create_int(size, size);
}

GAMatrixInt* GA_matrix_ref_int(GAMatrixInt* matrix)
{
    matrix->refs++;
    return matrix;
}

void GA_matrix_destroy_int(GAMatrixInt* matrix)
{
    matrix->refs--;
    if (matrix->refs == 0)
    {
        if (matrix->elts != 0)
        {
            int i;
            for (i = 0; i < matrix->rows; i++)
                if (matrix->elts[i] != 0)
                {
                    GA_free((char*)matrix->elts[i]);
                    matrix->elts[i] = 0;
                }
            GA_free((char*)matrix->elts);
            matrix->elts = 0;
        }
        GA_free((char*)matrix);
    }
}

int* GA_matrix_get_elt_int(GAMatrixInt* matrix, int row, int col)
{
    if ((row < matrix->rows)
        && (col < matrix->cols))
        return matrix->elts[row] + col;
    if (row >= matrix->rows)
        GA_msg()("[GA_matrix_get_elt_int] "
            "Row index out of range.", GA_MSG_ERROR);
    if (col >= matrix->cols)
        GA_msg()("[GA_matrix_get_elt_int] "
            "Column index out of range.", GA_MSG_ERROR);
    return 0;
}

GAMatrixInt* GA_matrix_init_from_array_int(GAMatrixInt* matrix, int* source, 
    int srcSize)
{
    if ((matrix->rows * matrix->cols) != srcSize)
    {
        GA_msg()("[GA_matrix_init_from_array_int] "
            "Target matrix has wrong size.", GA_MSG_ERROR);
        return 0;
    }
    int i;
    int j;
    for (i = 0; i < matrix->rows; i++)
        for (j = 0; j < matrix->cols; j++)
            matrix->elts[i][j] = source[i * matrix->cols + j];
    return matrix;
}

GAMatrixInt* GA_matrix_create_from_array_int(int* source, int rows, int cols)
{
    GAMatrixInt* matrix = GA_matrix_create_int(rows, cols);
    if (matrix == 0)
        return 0;
    return GA_matrix_init_from_array_int(matrix, source, rows * cols);
}

GAMatrixInt* GA_matrix_init_zero_int(GAMatrixInt* matrix)
{
    int i;
    int j;
    for (i = 0; i < matrix->rows; i++)
        for (j = 0; j < matrix->cols; j++)
            matrix->elts[i][j] = 0;
    return matrix;
}

GAMatrixInt* GA_matrix_init_unit_int(GAMatrixInt* matrix)
{
    int i;
    int j;
    for (i = 0; i < matrix->rows; i++)
        for (j = 0; j < matrix->cols; j++)
            if (i != j)
                matrix->elts[i][j] = 0;
            else
                matrix->elts[i][j] = 1;
    return matrix;
}

GAMatrixInt* GA_matrix_print_int(GAMatrixInt* matrix)
{
    GA_msg()("[", GA_MSG_INFO);
    int j;
    for (j = 0; j < matrix->rows; j++)
    {
        GA_msg()("(", GA_MSG_INFO);
        int i;
        for (i = 0; i < matrix->cols; i++)
        {
            char* message = GA_alloc(100, sizeof(char));
            snprintf(message, 100, "%i", matrix->elts[j][i]);
            GA_msg()(message, GA_MSG_INFO);
            GA_free(message);
            if (i < (matrix->cols - 1))
                GA_msg()(", ", GA_MSG_INFO);
        }
        GA_msg()(")", GA_MSG_INFO);
        if (j < (matrix->rows - 1))
            GA_msg()(", ", GA_MSG_INFO);
    }
    GA_msg()("]", GA_MSG_INFO);
    return matrix;
}

GAMatrixReal* GA_matrix_create_real(int rows, int cols)
{
    GAMatrixReal* matrix = (GAMatrixReal*)GA_alloc(1, sizeof(GAMatrixReal));
    if (matrix == 0)
    {
        GA_msg()("[GA_matrix_create_real] "
            "Could not allocate matrix.", GA_MSG_ERROR);
        return 0;
    }
    matrix->rows = rows;
    matrix->cols = cols;
    matrix->refs = 1;
    matrix->elts = (double**)GA_alloc(rows, sizeof(double*));
    if (matrix->elts == 0)
    {
        GA_msg()("[GA_matrix_create_real] "
            "Could not allocate matrix rows.", GA_MSG_ERROR);
        GA_free((char*)matrix);
        return 0;
    }
    int i;
    for (i = 0; i < matrix->rows; i++)
    {
        matrix->elts[i] = (double*)GA_alloc(cols, sizeof(double));
        if (matrix->elts[i] == 0)
        {
            char* message = GA_alloc(256, sizeof(char));
            snprintf(message, 256, "[GA_matrix_create_real] "
                "Could not allocate matrix row %i.", i);
            GA_msg()(message, GA_MSG_ERROR);
            GA_free(message);
            int j;
            for (j = 0; j < i; j++)
                GA_free((char*)matrix->elts[j]);
            GA_free((char*)matrix->elts);
            GA_free((char*)matrix);
            return 0;
        }
    }
    return matrix;
}

GAMatrixReal* GA_matrix_create_square_real(int size)
{
    return GA_matrix_create_real(size, size);
}

GAMatrixReal* GA_matrix_ref_real(GAMatrixReal* matrix)
{
    matrix->refs++;
    return matrix;
}

void GA_matrix_destroy_real(GAMatrixReal* matrix)
{
    matrix->refs--;
    if (matrix->refs == 0)
    {
        if (matrix->elts != 0)
        {
            int i;
            for (i = 0; i < matrix->rows; i++)
                if (matrix->elts[i] != 0)
                {
                    GA_free((char*)matrix->elts[i]);
                    matrix->elts[i] = 0;
                }
            GA_free((char*)matrix->elts);
            matrix->elts = 0;
        }
        GA_free((char*)matrix);
    }
}

double* GA_matrix_get_elt_real(GAMatrixReal* matrix, int row, int col)
{
    if ((row < matrix->rows)
        && (col < matrix->cols))
        return matrix->elts[row] + col;
    if (row >= matrix->rows)
        GA_msg()("[GA_matrix_get_elt_real] "
            "Row index out of range.", GA_MSG_ERROR);
    if (col >= matrix->cols)
        GA_msg()("[GA_matrix_get_elt_real] "
            "Column index out of range.", GA_MSG_ERROR);
    return 0;
}

GAMatrixReal* GA_matrix_init_from_array_real(GAMatrixReal* matrix, 
    double* source, int srcSize)
{
    if ((matrix->rows * matrix->cols) != srcSize)
    {
        GA_msg()("[GA_matrix_init_from_array_real] "
            "Target matrix has wrong size.", GA_MSG_ERROR);
        return 0;
    }
    int i;
    int j;
    for (i = 0; i < matrix->rows; i++)
        for (j = 0; j < matrix->cols; j++)
            matrix->elts[i][j] = source[i * matrix->cols + j];
    return matrix;
}

GAMatrixReal* GA_matrix_create_from_array_real(double* source, int rows, 
    int cols)
{
    GAMatrixReal* matrix = GA_matrix_create_real(rows, cols);
    if (matrix == 0)
        return 0;
    return GA_matrix_init_from_array_real(matrix, source, rows * cols);
}

GAMatrixReal* GA_matrix_init_zero_real(GAMatrixReal* matrix)
{
    int i;
    int j;
    for (i = 0; i < matrix->rows; i++)
        for (j = 0; j < matrix->cols; j++)
            matrix->elts[i][j] = 0.0;
    return matrix;
}

GAMatrixReal* GA_matrix_init_unit_real(GAMatrixReal* matrix)
{
    int i;
    int j;
    for (i = 0; i < matrix->rows; i++)
        for (j = 0; j < matrix->cols; j++)
            if (i != j)
                matrix->elts[i][j] = 0.0;
            else
                matrix->elts[i][j] = 1.0;
    return matrix;
}

GAMatrixInt* GA_matrix_to_bin_real(GAMatrixReal* matrix, GAVectorReal* lookup, 
    GAClampMode clamp)
{
    GAMatrixInt* result = GA_matrix_create_int(matrix->rows, matrix->cols);
    if (result == 0)
        return 0;
    int i;
    int j;
    for (i = 0; i < matrix->rows; i++)
        for (j = 0; j < matrix->cols; j++)
            result->elts[i][j] = GA_get_bin_number(
                matrix->elts[i][j], lookup, clamp);
    /* ----- DEBUG ----- //
    GA_msg()("[GA_matrix_to_bin_real] Input matrix: ", GA_MSG_DEBUG);
    GA_matrix_print_real(matrix);
    GA_msg()("\n", GA_MSG_INFO);
    GA_msg()("[GA_matrix_to_bin_real] Lookup vector: ", GA_MSG_DEBUG);
    GA_vector_print_real(lookup);
    GA_msg()("\n", GA_MSG_INFO);
    GA_msg()("[GA_matrix_to_bin_int] result: ", GA_MSG_DEBUG);
    GA_matrix_print_int(result);
    GA_msg()("\n", GA_MSG_INFO);
    // ----- DEBUG ----- */
    return result;
}

GAMatrixReal* GA_matrix_print_real(GAMatrixReal* matrix)
{
    GA_msg()("[", GA_MSG_INFO);
    int j;
    for (j = 0; j < matrix->rows; j++)
    {
        GA_msg()("(", GA_MSG_INFO);
        int i;
        for (i = 0; i < matrix->cols; i++)
        {
            char* message = GA_alloc(100, sizeof(char));
            snprintf(message, 100, "%f", matrix->elts[j][i]);
            GA_msg()(message, GA_MSG_INFO);
            GA_free(message);
            if (i < (matrix->cols - 1))
                GA_msg()(", ", GA_MSG_INFO);
        }
        GA_msg()(")", GA_MSG_INFO);
        if (j < (matrix->rows - 1))
            GA_msg()(", ", GA_MSG_INFO);
    }
    GA_msg()("]", GA_MSG_INFO);
    return matrix;
}
