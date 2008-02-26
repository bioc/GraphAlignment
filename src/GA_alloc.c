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
 * Memory allocation.
 * ----------------------------------------------------------------------------
 */

/** \file GA_alloc.c
 * \brief Memory allocation (implementation).
 */

#include <stdlib.h>
#include <stdio.h>
#include "GA_alloc.h"

/** Global allocation function.
 */
GAAllocFunc GA_ALLOC_FUNC = (GAAllocFunc)calloc;

/** Global freeing function.
 */
GAFreeFunc GA_FREE_FUNC = (GAFreeFunc)free;

void GA_set_alloc_funcs(GAAllocFunc allocFunc, GAFreeFunc freeFunc)
{
    GA_ALLOC_FUNC = allocFunc;
    GA_FREE_FUNC = freeFunc;
}

char* GA_alloc(size_t numElem, int eltSize)
{
    return GA_ALLOC_FUNC(numElem, eltSize);
}

void GA_free(char* memLoc)
{
    return GA_FREE_FUNC(memLoc);
}
