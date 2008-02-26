#ifndef GA_VECTOR
#define GA_VECTOR
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
 * Vector types.
 * ----------------------------------------------------------------------------
 */

/** \file GA_vector.h
 * \brief Vector types.
 *
 * This module provides types which represent vectors of integers or real 
 * numbers. Memory management for the types is handled transparently by the 
 * API functions which are provided by the module. The type includes 
 * information about the size of the vector, so no additional variables are 
 * required.
 */

#ifdef __cplusplus
extern "C"
{
#endif

/** A vector of integers (implementation).
 *
 * This type holds the elements of a vector of integers. Memory for this type 
 * is managed automatically if the interface functions from the API (as 
 * defined in GA_vector.h) are used. To create a new vector, use the 
 * GA_vector_create_int() function. To reference a vector, use 
 * GA_vector_ref_int(). To release a reference to a vector, use 
 * GA_vector_destroy_int(). The memory for the vector instance and the vector 
 * elements will be freed if the reference count reaches zero.
 *
 * Vector elements may be accessed directly by dereferencing the elts pointer, 
 * for example
 *
 * <tt>v->elts[i]</tt>
 *
 * to access the element with index i of the vector \c v. The size of a vector 
 * is stored in the \c size member of the type. It is also possible to do a 
 * checked access using GA_vector_get_elt_int(). If this function is used, an 
 * error message will be sent if the index is out of range.
 */
struct GAVectorInt_Impl
{
    /** Elements of the vector.
     */
    int* elts;
    /** Number of elements.
     */
    int size;
    /** Reference count.
     */
    int refs;
};

/** A vector of integers.
 */
typedef struct GAVectorInt_Impl GAVectorInt;

/** Create vector (int).
 *
 * Create a vector of integers. The new vector will be referenced and should 
 * be destroyed by using GA_vector_destroy_int() when it is not needed 
 * anymore.
 *
 * \param size Number of elements.
 *
 * \return Pointer to a vector, or 0 if an error occurs.
 *
 * \sa GA_vector_destroy_int
 */
GAVectorInt* GA_vector_create_int(int size);

/** Add reference (int).
 *
 * Add a reference for a vector. The user of this function is responsible 
 * for removing the reference using GA_vector_destroy_int().
 *
 * \param vec Vector.
 *
 * \return The vector.
 *
 * \sa GA_vector_destroy_int
 */
GAVectorInt* GA_vector_ref_int(GAVectorInt* vec);

/** Destroy vector (int).
 *
 * Remove a reference from a vector. If the reference count drops to zero, 
 * all resources allocated for the vector will be freed and the vector itself 
 * will be destroyed.
 *
 * \param vec Vector.
 */
void GA_vector_destroy_int(GAVectorInt* vec);

/** Get vector element (int).
 *
 * Get the element of the vector with the specified index. An error will be 
 * reported if the index is out of range.
 *
 * \param vec vector
 * \param index index
 *
 * \return specified element, or 0 if the index is out of range
 */
int* GA_vector_get_elt_int(GAVectorInt* vec, int index);

/** Initialize vector from array (int).
 *
 * Initialize a vector of integers from an array of integers.
 *
 * \param vec vector
 * \param source source array
 * \param srcSize size of source array
 *
 * \return the vector
 */
GAVectorInt* GA_vector_init_from_array_int(GAVectorInt* vec, int* source, 
    int srcSize);

/** Create vector from array (int).
 *
 * Create a vector of integers from an array of integers. The new vector 
 * will be referenced and should be destroyed by using GA_vector_destroy_int() 
 * when it is not needed anymore.
 *
 * \param source source array
 * \param srcSize size of source array
 *
 * \return new vector, or 0 if an error occurs
 */
GAVectorInt* GA_vector_create_from_array_int(int* source, int srcSize);

/** Initialize vector to zero (int).
 *
 * Set all elements of a vector to zero.
 *
 * \param vec vector
 *
 * \return the vector
 */
GAVectorInt* GA_vector_init_zero_int(GAVectorInt* vec);

/** Print vector (int).
 *
 * Print a vector of integers in a nice, readable way.
 *
 * \param vec vector
 *
 * \return the vector
 */
GAVectorInt* GA_vector_print_int(GAVectorInt* vec);

/** Invert permutation.
 *
 * Invert the permutation defined by the specified vector. The result will 
 * be referenced and should be destroyed by using GA_vector_destroy_int() 
 * when it is not needed anymore.
 *
 * \param vec vector
 *
 * \return inverse of the permutation
 */
GAVectorInt* GA_invert_permutation_int(GAVectorInt* vec);

/** A vector of real numbers (implementation).
 *
 * This type holds the elements of a vector of real numbers. Memory for this 
 * type is managed automatically if the interface functions from the API (as 
 * defined in GA_vector.h) are used. To create a new vector, use the 
 * GA_vector_create_real() function. To reference a vector, use 
 * GA_vector_ref_real(). To release a reference to a vector, use 
 * GA_vector_destroy_real(). The memory for the vector instance and the vector 
 * elements will be freed if the reference count reaches zero.
 *
 * Vector elements may be accessed directly by dereferencing the elts pointer, 
 * for example
 *
 * <tt>v->elts[i]</tt>
 *
 * to access the element with index i of the vector \c v. The size of a vector 
 * is stored in the \c size member of the type. It is also possible to do a 
 * checked access using GA_vector_get_elt_real(). If this function is used, an 
 * error message will be sent if the index is out of range.
 */
struct GAVectorReal_Impl
{
    /** Elements of the vector.
     */
    double* elts;
    /** Number of elements.
     */
    int size;
    /** Reference count.
     */
    int refs;
};

/** A vector of real numbers.
 */
typedef struct GAVectorReal_Impl GAVectorReal;

/** Create vector (real).
 *
 * Create a vector of real numbers. The new vector will be referenced and 
 * should be destroyed by using GA_vector_destroy_int() when it is not needed 
 * anymore.
 *
 * \param size Number of elements.
 *
 * \return Pointer to a vector, or 0 if an error occurs.
 *
 * \sa GA_vector_destroy_real
 */
GAVectorReal* GA_vector_create_real(int size);

/** Add reference (real).
 *
 * Add a reference for a vector. The user of this function is responsible 
 * for removing the reference using GA_vector_destroy_real().
 *
 * \param vec Vector.
 *
 * \return The vector.
 *
 * \sa GA_vector_destroy_real
 */
GAVectorReal* GA_vector_ref_real(GAVectorReal* vec);

/** Destroy vector (real).
 *
 * Remove a reference from a vector. If the reference count drops to zero, 
 * all resources allocated for the vector will be freed and the vector itself 
 * will be destroyed.
 *
 * \param vec vector.
 */
void GA_vector_destroy_real(GAVectorReal* vec);

/** Get vector element (real).
 *
 * Get the element of the vector with the specified index. An error will be 
 * reported if the index is out of range.
 *
 * \param vec vector
 * \param index index
 *
 * \return specified element, or 0 if the index is out of range
 */
double* GA_vector_get_elt_real(GAVectorReal* vec, int index);

/** Initialize vector from array (real).
 *
 * Initialize a vector of real numbers from an array of real numbers.
 *
 * \param vec vector
 * \param source source array
 * \param srcSize size of source array
 *
 * \return the vector
 */
GAVectorReal* GA_vector_init_from_array_real(GAVectorReal* vec, 
    double* source, int srcSize);

/** Create vector from array (real).
 *
 * Create a vector of real numbers from an array of real numbers. The new 
 * vector will be referenced and should be destroyed by using 
 * GA_vector_destroy_real() when it is not needed anymore.
 *
 * \param source source array
 * \param srcSize size of source array
 *
 * \return new vector, or 0 if an error occurs
 */
GAVectorReal* GA_vector_create_from_array_real(double* source, int srcSize);

/** Initialize vector to zero (real).
 *
 * Set all elements of a vector to zero.
 *
 * \param vec vector
 *
 * \return the vector
 */
GAVectorReal* GA_vector_init_zero_real(GAVectorReal* vec);

/** Clamp mode (implementation).
 *
 * The clamp mode specifies whether a value should be clamped to a range if 
 * it is outside of the boundaries of the range.
 */
enum GAClampMode_Impl
{
    /** Clamp mode: enabled.
     */
    GA_CLAMP_ENABLED = 1,
    /** Clamp mode: disabled.
     */
    GA_CLAMP_DISABLED = 0,
};

/** Clamp mode.
 */
typedef enum GAClampMode_Impl GAClampMode;

/** Get bin number.
 *
 * Get the bin number for the argument using the specified lookup vector. 
 * If \c clamp is set to GA_CLAMP_ENABLED, values which are outside of the 
 * boundaries of the lookup range will be clamped to the bins at the end of 
 * the range. If \c clamp is set to GA_CLAMP_DISABLED, out-of-range values 
 * will be reported as errors.
 *
 * \param x number which should be binned
 * \param lookup lookup vector
 * \param clamp clamp the result to the lookup range
 *
 * \return bin number, or -1 if an error occurs
 */
int GA_get_bin_number(double x, GAVectorReal* lookup, GAClampMode clamp);

/** Convert to bin vector (real).
 *
 * Convert a vector of real numbers to a vector of integer bin numbers 
 * according to the specified lookup vector. The new vector will be 
 * referenced and should be destroyed by using GA_vector_destroy_int() when 
 * it is not needed anymore.
 *
 * \param vec vector
 * \param lookup lookup vector
 * \param clamp clamp values to the lookup range
 *
 * \return vector of bin numbers
 *
 * \sa GA_get_bin_number()
 */
GAVectorInt* GA_vector_to_bin_real(GAVectorReal* vec, GAVectorReal* lookup, 
    GAClampMode clamp);

/** Print vector (real).
 *
 * Print a vector of real numbers in a nice, readable way.
 *
 * \param vec vector
 *
 * \return the vector
 */
GAVectorReal* GA_vector_print_real(GAVectorReal* vec);

#ifdef __cplusplus
}
#endif
#endif
