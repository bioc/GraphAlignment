#ifndef GA_VECTOR_R
#define GA_VECTOR_R
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
 * R utility functions for vector types.
 * ----------------------------------------------------------------------------
 */

/** \file GA_vector_R.h
 * \brief R utility functions for vector types.
 *
 * This module provides functions for handling conversions between the vector 
 * types of the graph alignment package C implementation and the internal R 
 * objects (<tt>SEXP</tt>) which represent vectors.
 */

#include "R.h"
#include "Rinternals.h"
#include "Rdefines.h"
#include "GA_vector.h"

#ifdef __cplusplus
extern "C"
{
#endif

/** Create vector from R object (int).
 *
 * Create a vector of integers from an R object. The new vector will be 
 * referenced and should be destroyed by using GA_vector_destroy_int() when 
 * it is not needed anymore.
 *
 * \param robj R object.
 *
 * \return Pointer to a vector, or 0 if an error occurs.
 *
 * \sa GA_vector_destroy_int
 */
GAVectorInt* GA_vector_from_R_int(SEXP robj);

/** Create R object from vector (int).
 *
 * Create an R object from a vector of integers.
 *
 * \param vec Vector.
 *
 * \return R object.
 */
SEXP GA_vector_to_R_int(GAVectorInt* vec);

/** Create vector from R object (real).
 *
 * Create a vector of real numbers from an R object. The new vector will be 
 * referenced and should be destroyed by using GA_vector_destroy_real() when 
 * it is not needed anymore.
 *
 * \param robj R object.
 *
 * \return Pointer to a vector, or 0 if an error occurs.
 *
 * \sa GA_vector_destroy_real
 */
GAVectorReal* GA_vector_from_R_real(SEXP robj);

/** Create R object from vector (real).
 *
 * Create an R object from a vector of real numbers.
 *
 * \param vec Vector.
 *
 * \return R object.
 */
SEXP GA_vector_to_R_real(GAVectorReal* vec);

/** Get clamp mode from R object (real).
 *
 * Get the clamp mode corresponding to the value of the specified R object.
 *
 * \param robj R object
 *
 * \return clamp mode
 */
GAClampMode GA_clamp_mode_from_R(SEXP robj);

#ifdef __cplusplus
}
#endif
#endif
