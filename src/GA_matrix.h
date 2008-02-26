#ifndef GA_MATRIX
#define GA_MATRIX
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

/** \file GA_matrix.h
 * \brief Matrix types.
 *
 * This module provides types which represent two-dimensional matrices of 
 * integers or real numbers. Memory management for the types is handled 
 * transparently by the API functions which are provided by the module. The 
 * type includes information about the size of the matrix, so no additional 
 * variables are required.
 */

#include "GA_vector.h"

#ifdef __cplusplus
extern "C"
{
#endif

/** A matrix of integers (implementation).
 *
 * This type holds the elements of a two-dimensional matrix of integers. 
 * Memory for this type is managed automatically if the interface functions 
 * from the API (as defined in GA_matrix.h) are used. To create a new matrix, 
 * use one of the GA_matrix_create* functions. To reference a matrix, use 
 * GA_matrix_ref_int(). To release a reference to a matrix, use 
 * GA_matrix_destroy_int(). The memory for the matrix instance and the matrix 
 * elements will be freed if the reference count reaches zero.
 *
 * Matrix elements may be accessed directly by dereferencing the elts pointer, 
 * for example
 *
 * <tt>m->elts[i][j]</tt>
 *
 * to access the element with index (i, j) of the matrix \c m. The size of 
 * the matrix is stored in the \c rows and \c cols members of the type. It is 
 * also possible to do a checked access using GA_matrix_get_elt_int(). If this 
 * function is used, an error message will be sent if an index is out of range.
 */
struct GAMatrixInt_Impl
{
    /** Elements of the matrix.
     */
    int** elts;
    /** Number of rows.
     */
    int rows;
    /** Number of columns.
     */
    int cols;
    /** Reference count.
     */
    int refs;
};

/** A matrix of integers.
 */
typedef struct GAMatrixInt_Impl GAMatrixInt;

/** Create matrix (int).
 *
 * Create a matrix of integers. The new matrix will be referenced and should 
 * be destroyed by using GA_matrix_destroy_int() when it is not needed 
 * anymore.
 *
 * \param rows Number of rows.
 * \param cols Number of columns.
 *
 * \return Pointer to a matrix, or 0 if an error occurs.
 *
 * \sa GA_matrix_destroy_int
 */
GAMatrixInt* GA_matrix_create_int(int rows, int cols);

/** Create square matrix (int).
 *
 * Create a square matrix of integers. The new matrix will be referenced 
 * and should be destroyed by using GA_matrix_destroy_int() when it is not 
 * needed anymore.
 *
 * \param size Size of the matrix.
 *
 * \return Pointer to a matrix, or 0 if an error occurs.
 *
 * \sa GA_matrix_destroy_int
 */
GAMatrixInt* GA_matrix_create_square_int(int size);

/** Add reference.
 *
 * Add a reference for a matrix. The user of this function is responsible 
 * for removing the reference using GA_matrix_destroy_int().
 *
 * \param matrix Matrix.
 *
 * \return The matrix.
 *
 * \sa GA_matrix_destroy_int
 */
GAMatrixInt* GA_matrix_ref_int(GAMatrixInt* matrix);

/** Destroy matrix.
 *
 * Remove a reference from a matrix. If the reference count drops to zero, 
 * all resources allocated for the matrix will be freed and the matrix itself 
 * will be destroyed.
 *
 * \param matrix Matrix.
 */
void GA_matrix_destroy_int(GAMatrixInt* matrix);

/** Get matrix element (int).
 *
 * Get the element of the matrix with the specified indices. An error will be 
 * reported if one of the indices is out of range.
 *
 * \param matrix matrix
 * \param row row index
 * \param col column index
 *
 * \return specified element, or 0 if the index is out of range
 */
int* GA_matrix_get_elt_int(GAMatrixInt* matrix, int row, int col);

/** Initialize matrix from array (int).
 *
 * Initialize a matrix of integers from an array of integers.
 *
 * \param matrix matrix
 * \param source source array
 * \param srcSize size of source array
 *
 * \return the matrix
 */
GAMatrixInt* GA_matrix_init_from_array_int(GAMatrixInt* matrix, int* source, 
    int srcSize);

/** Create matrix from array (int).
 *
 * Create a matrix of integers from an array of integers. The new matrix 
 * will be referenced and should be destroyed by using GA_matrix_destroy_int() 
 * when it is not needed anymore.
 *
 * \param source source array
 * \param rows number of rows of the new matrix
 * \param cols number of columns of the new matrix
 *
 * \return new matrix, or 0 if an error occurs
 */
GAMatrixInt* GA_matrix_create_from_array_int(int* source, int rows, int cols);

/** Initialize matrix to zero.
 *
 * Set all elements of a matrix to zero.
 *
 * \param matrix Matrix.
 *
 * \return The matrix.
 */
GAMatrixInt* GA_matrix_init_zero_int(GAMatrixInt* matrix);

/** Initialize matrix to unit.
 *
 * Set the elements of a matrix to the elements of the unit matrix of the 
 * appropriate size.
 *
 * \param matrix Matrix.
 *
 * \return The matrix.
 */
GAMatrixInt* GA_matrix_init_unit_int(GAMatrixInt* matrix);

/** Print matrix (int).
 *
 * Print a matrix of integers in a nice, readable way.
 *
 * \param matrix matrix
 *
 * \return the matrix
 */
GAMatrixInt* GA_matrix_print_int(GAMatrixInt* matrix);

/** A matrix of real numbers (implementation).
 *
 * This type holds the elements of a two-dimensional matrix of real numbers. 
 * Memory for this type is managed automatically if the interface functions 
 * from the API (as defined in GA_matrix.h) are used. To create a new matrix, 
 * use one of the GA_matrix_create* functions. To reference a matrix, use 
 * GA_matrix_ref_real(). To release a reference to a matrix, use 
 * GA_matrix_destroy_real(). The memory for the matrix instance and the matrix 
 * elements will be freed if the reference count reaches zero.
 *
 * Matrix elements may be accessed directly by dereferencing the elts pointer, 
 * for example
 *
 * <tt>m->elts[i][j]</tt>
 *
 * to access the element with index (i, j) of the matrix \c m. The size of 
 * the matrix is stored in the \c rows and \c cols members of the type. It is 
 * also possible to do a checked access using GA_matrix_get_elt_real(). If this 
 * function is used, an error message will be sent if an index is out of range.
 */
struct GAMatrixReal_Impl
{
    /** Elements of the matrix.
     */
    double** elts;
    /** Number of rows.
     */
    int rows;
    /** Number of columns.
     */
    int cols;
    /** Reference count.
     */
    int refs;
};

/** A matrix of real numbers.
 */
typedef struct GAMatrixReal_Impl GAMatrixReal;

/** Create matrix (real).
 *
 * Create a matrix of real numbers. The new matrix will be referenced and 
 * should be destroyed by using GA_matrix_destroy_real() when it is not needed 
 * anymore.
 *
 * \param rows Number of rows.
 * \param cols Number of columns.
 *
 * \return Pointer to a matrix, or 0 if an error occurs.
 *
 * \sa GA_matrix_destroy_real
 */
GAMatrixReal* GA_matrix_create_real(int rows, int cols);

/** Create square matrix (real).
 *
 * Create a square matrix of real numbers. The new matrix will be referenced 
 * and should be destroyed by using GA_matrix_destroy_real() when it is not 
 * needed anymore.
 *
 * \param size Size of the matrix.
 *
 * \return Pointer to a matrix, or 0 if an error occurs.
 *
 * \sa GA_matrix_destroy_real
 */
GAMatrixReal* GA_matrix_create_square_real(int size);

/** Add reference.
 *
 * Add a reference for a matrix. The user of this function is responsible 
 * for removing the reference using GA_matrix_destroy_real().
 *
 * \param matrix Matrix.
 *
 * \return The matrix.
 *
 * \sa GA_matrix_destroy_real
 */
GAMatrixReal* GA_matrix_ref_real(GAMatrixReal* matrix);

/** Destroy matrix.
 *
 * Remove a reference from a matrix. If the reference count drops to zero, 
 * all resources allocated for the matrix will be freed and the matrix itself 
 * will be destroyed.
 *
 * \param matrix Matrix.
 */
void GA_matrix_destroy_real(GAMatrixReal* matrix);

/** Get matrix element (real).
 *
 * Get the element of the matrix with the specified indices. An error will be 
 * reported if one of the indices is out of range.
 *
 * \param matrix matrix
 * \param row row index
 * \param col column index
 *
 * \return specified element, or 0 if the index is out of range
 */
double* GA_matrix_get_elt_real(GAMatrixReal* matrix, int row, int col);

/** Initialize matrix from array (real).
 *
 * Initialize a matrix of real numbers from an array of real numbers.
 *
 * \param matrix matrix
 * \param source source array
 * \param srcSize size of source array
 *
 * \return the matrix
 */
GAMatrixReal* GA_matrix_init_from_array_real(GAMatrixReal* matrix, 
    double* source, int srcSize);

/** Create matrix from array (real).
 *
 * Create a matrix of real numbers from an array of real numbers. The new 
 * matrix will be referenced and should be destroyed by using 
 * GA_matrix_destroy_real() when it is not needed anymore.
 *
 * \param source source array
 * \param rows number of rows of the new matrix
 * \param cols number of columns of the new matrix
 *
 * \return new matrix, or 0 if an error occurs
 */
GAMatrixReal* GA_matrix_create_from_array_real(double* source, int rows, 
    int cols);

/** Initialize matrix to zero.
 *
 * Set all elements of a matrix to zero.
 *
 * \param matrix Matrix.
 *
 * \return The matrix.
 */
GAMatrixReal* GA_matrix_init_zero_real(GAMatrixReal* matrix);

/** Initialize matrix to unit.
 *
 * Set the elements of a matrix to the elements of the unit matrix of the 
 * appropriate size.
 *
 * \param matrix Matrix.
 *
 * \return The matrix.
 */
GAMatrixReal* GA_matrix_init_unit_real(GAMatrixReal* matrix);

/** Convert to bin matrix (real).
 *
 * Convert a matrix of real numbers to a matrix of integer bin numbers 
 * according to the specified lookup vector. The new matrix will be 
 * referenced and should be destroyed by using GA_matrix_destroy_real() when 
 * it is not needed anymore.
 *
 * \param matrix Matrix
 * \param lookup Lookup vector.
 * \param clamp clamp values to the lookup range
 *
 * \return Matrix of bin numbers
 *
 * \sa GA_get_bin_number()
 */
GAMatrixInt* GA_matrix_to_bin_real(GAMatrixReal* matrix, GAVectorReal* lookup, 
    GAClampMode clamp);

/** Print matrix (real).
 *
 * Print a matrix of integers in a nice, readable way.
 *
 * \param matrix matrix
 *
 * \return the matrix
 */
GAMatrixReal* GA_matrix_print_real(GAMatrixReal* matrix);

#ifdef __cplusplus
}
#endif
#endif
