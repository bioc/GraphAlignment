#ifndef GA_ALLOC
#define GA_ALLOC
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
 * Memory alllocation.
 * ----------------------------------------------------------------------------
 */

/** \file GA_alloc.h
 * \brief Memory allocation.
 *
 * This module provides memory allocation services to the package. The 
 * allocation functions to be used can be set at run-time.
 */

#include <stdlib.h>

#ifdef __cplusplus
extern "C"
{
#endif

/** Memory allocation function.
 *
 * This is the specification of the memory allocation function which can be set 
 * using GA_set_alloc_funcs.
 */
typedef char* (*GAAllocFunc)(size_t, int);

/** Memory freeing function.
 *
 * This is the specification of the memory freeing function which can be set 
 * using GA_set_alloc_funcs.
 */
typedef void (*GAFreeFunc)(char*);

/** Set allocation functions.
 *
 * Set the memory allocation functions to be used.
 *
 * \param allocFunc Allocation function.
 * \param freeFunc Freeing function.
 */
void GA_set_alloc_funcs(GAAllocFunc allocFunc, GAFreeFunc freeFunc);

/** Allocate memory.
 *
 * Allocate memory using the currently set allocation function.
 *
 * \param numElem Number of elements to allocate.
 * \param eltSize Element size.
 *
 * \return Result of allocation.
 */
char* GA_alloc(size_t numElem, int eltSize);

/** Free memory.
 *
 * Free memory using the currently set freeing function.
 *
 * \param memLoc Memory location to be freed.
 */
void GA_free(char* memLoc);

#ifdef __cplusplus
}
#endif
#endif
