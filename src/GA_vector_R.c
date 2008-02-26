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

/** \file GA_vector_R.c
 * \brief R utility functions for vector types (implementation).
 */

#include "GA_alloc.h"
#include "GA_message.h"
#include "GA_vector_R.h"

GAVectorInt* GA_vector_from_R_int(SEXP robj)
{
    PROTECT(robj);
    SEXPTYPE vecType = TYPEOF(robj);
    if ((vecType != INTSXP)
        && (vecType != REALSXP))
    {
        char* message = GA_alloc(256, sizeof(char));
        snprintf(message, 256, 
            "[GA_vector_from_R_int] Input is not a vector of real or "
            "integer values (actual type: %i).", vecType);
        GA_msg()(message, GA_MSG_ERROR);
        GA_free(message);
        UNPROTECT(1);
        return 0;
    }
    int size = LENGTH(robj);
    GAVectorInt* vec = GA_vector_create_int(size);
    if (vec == 0)
    {
        UNPROTECT(1);
        return 0;
    }
    int* inputRaw = INTEGER(coerceVector(robj, INTSXP));
    int i;
    for (i = 0; i < size; i++)
        vec->elts[i] = inputRaw[i];
    UNPROTECT(1);
    return vec;
}

SEXP GA_vector_to_R_int(GAVectorInt* vec)
{
    SEXP result;
    PROTECT(result = allocVector(INTSXP, vec->size));
    int* resultRaw = INTEGER(coerceVector(result, INTSXP));
    int i;
    for (i = 0; i < vec->size; i++)
        resultRaw[i] = vec->elts[i];
    UNPROTECT(1);
    return result;
}

GAVectorReal* GA_vector_from_R_real(SEXP robj)
{
    PROTECT(robj);
    SEXPTYPE vecType = TYPEOF(robj);
    if ((vecType != INTSXP)
        && (vecType != REALSXP))
    {
        char* message = GA_alloc(256, sizeof(char));
        snprintf(message, 256, 
            "[GA_vector_from_R_real] Input is not a vector of real or "
            "integer values (actual type: %i).", vecType);
        GA_msg()(message, GA_MSG_ERROR);
        GA_free(message);
        UNPROTECT(1);
        return 0;
    }
    int size = LENGTH(robj);
    GAVectorReal* vec = GA_vector_create_real(size);
    if (vec == 0)
    {
        UNPROTECT(1);
        return 0;
    }
    double* inputRaw = REAL(coerceVector(robj, REALSXP));
    int i;
    for (i = 0; i < size; i++)
        vec->elts[i] = inputRaw[i];
    UNPROTECT(1);
    return vec;
}

SEXP GA_vector_to_R_real(GAVectorReal* vec)
{
    SEXP result;
    PROTECT(result = allocVector(REALSXP, vec->size));
    double* resultRaw = REAL(coerceVector(result, REALSXP));
    int i;
    for (i = 0; i < vec->size; i++)
        resultRaw[i] = vec->elts[i];
    UNPROTECT(1);
    return result;
}

GAClampMode GA_clamp_mode_from_R(SEXP robj)
{
    PROTECT(robj);
    SEXPTYPE objType = TYPEOF(robj);
    if ((objType != LGLSXP)
        && (objType != INTSXP)
        && (objType != REALSXP))
    {
        char* message = GA_alloc(256, sizeof(char));
        snprintf(message, 256, 
            "[GA_clamp_mode_from_R] Input is not a logical, real or "
            "integer value (actual type: %i).", objType);
        GA_msg()(message, GA_MSG_ERROR);
        GA_free(message);
        UNPROTECT(1);
        return 0;
    }
    int size = LENGTH(robj);
    int* inputRaw = LOGICAL(coerceVector(robj, LGLSXP));
    if (inputRaw[0] == 0)
    {
        UNPROTECT(1);
        return GA_CLAMP_DISABLED;
    }
    UNPROTECT(1);
    return GA_CLAMP_ENABLED;
}
