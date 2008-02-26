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

/** \file GA_vector.c
 * \brief Vector types implementation.
 */

#include <stdlib.h>
#include <stdio.h>
#include "GA_alloc.h"
#include "GA_message.h"
#include "GA_vector.h"

GAVectorInt* GA_vector_create_int(int size)
{
    GAVectorInt* vec = (GAVectorInt*)GA_alloc(1, sizeof(GAVectorInt));
    if (vec == 0)
    {
        GA_msg()("[GA_vector_create_int] "
            "Could not allocate vector.", GA_MSG_ERROR);
        return 0;
    }
    vec->size = size;
    vec->refs = 1;
    vec->elts = (int*)GA_alloc(size, sizeof(int));
    if (vec->elts == 0)
    {
        GA_msg()("[GA_vector_create_int] "
            "Could not allocate vector elements.", GA_MSG_ERROR);
        GA_free((char*)vec);
        return 0;
    }
    return vec;
}

GAVectorInt* GA_vector_ref_int(GAVectorInt* vec)
{
    vec->refs++;
    return vec;
}

void GA_vector_destroy_int(GAVectorInt* vec)
{
    vec->refs--;
    if (vec->refs == 0)
    {
        if (vec->elts != 0)
        {
            GA_free((char*)vec->elts);
            vec->elts = 0;
        }
        GA_free((char*)vec);
    }
}

int* GA_vector_get_elt_int(GAVectorInt* vec, int index)
{
    if (index < vec->size)
        return vec->elts + index;
    GA_msg()("[GA_vector_get_elt_int] "
        "Index out of range.", GA_MSG_ERROR);
    return 0;
}

GAVectorInt* GA_vector_init_from_array_int(GAVectorInt* vec, int* source, 
    int srcSize)
{
    if (vec->size != srcSize)
    {
        GA_msg()("[GA_vector_init_from_array_int] "
            "Target vector has wrong size.", GA_MSG_ERROR);
        return 0;
    }
    int i;
    for (i = 0; i < srcSize; i++)
        vec->elts[i] = source[i];
    return vec;
}

GAVectorInt* GA_vector_create_from_array_int(int* source, int srcSize)
{
    GAVectorInt* vec = GA_vector_create_int(srcSize);
    if (vec == 0)
        return 0;
    return GA_vector_init_from_array_int(vec, source, srcSize);
}

GAVectorInt* GA_vector_init_zero_int(GAVectorInt* vec)
{
    int i;
    for (i = 0; i < vec->size; i++)
        vec->elts[i] = 0;
    return vec;
}

GAVectorInt* GA_vector_print_int(GAVectorInt* vec)
{
    GA_msg()("(", GA_MSG_INFO);
    int i;
    for (i = 0; i < vec->size; i++)
    {
        char* message = GA_alloc(100, sizeof(char));
        snprintf(message, 100, "%i", vec->elts[i]);
        GA_msg()(message, GA_MSG_INFO);
        GA_free(message);
        if (i < (vec->size - 1))
            GA_msg()(", ", GA_MSG_INFO);
    }
    GA_msg()(")", GA_MSG_INFO);
    return vec;
}

GAVectorInt* GA_invert_permutation_int(GAVectorInt* vec)
{
    GAVectorInt* result = GA_vector_create_int(vec->size);
    int i;
    for (i = 0; i < vec->size; i++)
        *GA_vector_get_elt_int(result, vec->elts[i]) = i;
    return result;
}

GAVectorReal* GA_vector_create_real(int size)
{
    GAVectorReal* vec = (GAVectorReal*)GA_alloc(1, sizeof(GAVectorReal));
    if (vec == 0)
    {
        GA_msg()("[GA_vector_create_real] "
            "Could not allocate vector.", GA_MSG_ERROR);
        return 0;
    }
    vec->size = size;
    vec->refs = 1;
    vec->elts = (double*)GA_alloc(size, sizeof(double));
    if (vec->elts == 0)
    {
        GA_msg()("[GA_vector_create_real] "
            "Could not allocate vector elements.", GA_MSG_ERROR);
        GA_free((char*)vec);
        return 0;
    }
    return vec;
}

GAVectorReal* GA_vector_ref_real(GAVectorReal* vec)
{
    vec->refs++;
    return vec;
}

void GA_vector_destroy_real(GAVectorReal* vec)
{
    vec->refs--;
    if (vec->refs == 0)
    {
        if (vec->elts != 0)
        {
            GA_free((char*)vec->elts);
            vec->elts = 0;
        }
        GA_free((char*)vec);
    }
}

double* GA_vector_get_elt_real(GAVectorReal* vec, int index)
{
    if (index < vec->size)
        return vec->elts + index;
    GA_msg()("[GA_vector_get_elt_real] "
        "Index out of range.", GA_MSG_ERROR);
    return 0;
}

GAVectorReal* GA_vector_init_from_array_real(GAVectorReal* vec, 
    double* source, int srcSize)
{
    if (vec->size != srcSize)
    {
        GA_msg()("[GA_vector_init_from_array_int] "
            "Target vector has wrong size.", GA_MSG_ERROR);
        return 0;
    }
    int i;
    for (i = 0; i < srcSize; i++)
        vec->elts[i] = source[i];
    return vec;
}

GAVectorReal* GA_vector_create_from_array_real(double* source, int srcSize)
{
    GAVectorReal* vec = GA_vector_create_real(srcSize);
    if (vec == 0)
        return 0;
    return GA_vector_init_from_array_real(vec, source, srcSize);
}

GAVectorReal* GA_vector_init_zero_real(GAVectorReal* vec)
{
    int i;
    for (i = 0; i < vec->size; i++)
        vec->elts[i] = 0.0;
    return vec;
}

int GA_get_bin_number(double x, GAVectorReal* lookup, GAClampMode clamp)
{
    if (lookup->size == 0)
    {
        GA_msg()("[GA_get_bin_number] Lookup vector is empty.", 
            GA_MSG_ERROR);
        return -1;
    }
    if (lookup->size == 1)
    {
        /* There is no real lookup range. */
        if ((clamp == GA_CLAMP_DISABLED)
            && (x != lookup->elts[0]))
        {
            char* message = GA_alloc(256, sizeof(char));
            snprintf(message, 256, "[GA_get_bin_number] "
                "There is only a single lookup value and clamping is "
                "disabled, but the input value is not equal to the lookup "
                "value. Please make sure you have provided the correct "
                "lookup range and clamp mode (x = %f, lookup = %f).\n", 
                x, lookup->elts[0]);
            GA_msg()(message, GA_MSG_ERROR);
            GA_free(message);
            return -1;
        }
        /* Either the argument is clamped or it is the lookup value itself. */
        return 0;
    }
    if ((x < lookup->elts[0])
        || (x > lookup->elts[lookup->size - 1]))
    {
        /* There is a real lookup range, and the argument is outside of 
           the boundaries. */
        if (clamp == GA_CLAMP_DISABLED)
        {
            char* message = GA_alloc(256, sizeof(char));
            snprintf(message, 256, "[GA_get_bin_number] "
                "Argument is outside of lookup range and clamping is disabled. "
                "Please make sure you have provided the correct lookup range and "
                "clamp mode (x = %f, lower = %f, upper = %f).\n", 
                x, lookup->elts[0], lookup->elts[lookup->size - 1]);
            GA_msg()(message, GA_MSG_ERROR);
            GA_free(message);
            return -1;
        }
        if (x < lookup->elts[0])
            return 0;
        if (x > lookup->elts[lookup->size - 1])
            return lookup->size - 2;
    }
    /* Now, we have a real lookup range and the argument is within the 
       boundaries. */
    int result = 0;
    while (((result + 1) < (lookup->size - 1))
        && (x >= lookup->elts[result + 1]))
        result++;
    return result;
}

GAVectorInt* GA_vector_to_bin_real(GAVectorReal* vec, GAVectorReal* lookup, 
    GAClampMode clamp)
{
    GAVectorInt* result = GA_vector_create_int(vec->size);
    if (result == 0)
        return 0;
    int i;
    for (i = 0; i < vec->size; i++)
        result->elts[i] = GA_get_bin_number(vec->elts[i], lookup, clamp);
    /* ----- DEBUG ----- //
    GA_msg()("[GA_vector_to_bin_real] Input vector: ", GA_MSG_DEBUG);
    GA_vector_print_real(vec);
    GA_msg()("\n", GA_MSG_INFO);
    GA_msg()("[GA_vector_to_bin_real] Lookup vector: ", GA_MSG_DEBUG);
    GA_vector_print_real(lookup);
    GA_msg()("\n", GA_MSG_INFO);
    GA_msg()("[GA_vector_to_bin_int] result: ", GA_MSG_DEBUG);
    GA_vector_print_int(result);
    GA_msg()("\n", GA_MSG_INFO);
    // ----- DEBUG ----- */
    return result;
}

GAVectorReal* GA_vector_print_real(GAVectorReal* vec)
{
    GA_msg()("(", GA_MSG_INFO);
    int i;
    for (i = 0; i < vec->size; i++)
    {
        char* message = GA_alloc(100, sizeof(char));
        snprintf(message, 100, "%f", vec->elts[i]);
        GA_msg()(message, GA_MSG_INFO);
        GA_free(message);
        if (i < (vec->size - 1))
            GA_msg()(", ", GA_MSG_INFO);
    }
    GA_msg()(")", GA_MSG_INFO);
    return vec;
}
