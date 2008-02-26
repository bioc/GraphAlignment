#ifndef GA_MATRIX_R
#define GA_MATRIX_R
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

/** \file GA_matrix_R.h
 * \brief R utility functions for matrix types.
 *
 * This module provides functions for handling conversions between the matrix 
 * types of the graph alignment package C implementation and the internal R 
 * objects (<tt>SEXP</tt>) which represent matrices.
 */

#include "R.h"
#include "Rinternals.h"
#include "Rdefines.h"
#include "GA_matrix.h"

#ifdef __cplusplus
extern "C"
{
#endif

/** Create matrix from R object (int).
 *
 * Create a matrix of integers from an R object. The new matrix will be 
 * referenced and should be destroyed by using GA_matrix_destroy_int() when 
 * it is not needed anymore.
 *
 * \param robj R object.
 *
 * \return Pointer to a matrix, or 0 if an error occurs.
 *
 * \sa GA_matrix_destroy_int
 */
GAMatrixInt* GA_matrix_from_R_int(SEXP robj);

/** Create R object from matrix (int).
 *
 * Create an R object from a matrix of integers.
 *
 * \param matrix Matrix.
 *
 * \return R object.
 */
SEXP GA_matrix_to_R_int(GAMatrixInt* matrix);

/** Create matrix from R object (real).
 *
 * Create a matrix of real numbers from an R object. The new matrix will be 
 * referenced and should be destroyed by using GA_matrix_destroy_real() when 
 * it is not needed anymore.
 *
 * \param robj R object.
 *
 * \return Pointer to a matrix, or 0 if an error occurs.
 *
 * \sa GA_matrix_destroy_real
 */
GAMatrixReal* GA_matrix_from_R_real(SEXP robj);

/** Create R object from matrix (real).
 *
 * Create an R object from a matrix of real numbers.
 *
 * \param matrix Matrix.
 *
 * \return R object.
 */
SEXP GA_matrix_to_R_real(GAMatrixReal* matrix);

#ifdef __cplusplus
}
#endif
#endif
