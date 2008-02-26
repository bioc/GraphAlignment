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
 * Message reporting.
 * ----------------------------------------------------------------------------
 */

/** \file GA_message.c
 * \brief Message reporting (implementation).
 */

#include <stdlib.h>
#include <stdio.h>
#include "GA_message.h"

void GA_msg_default(const char* text, GAMessageLevel level)
{
    if (level == GA_MSG_DEBUG)
        printf("DEBUG: %s", text);
    else
    if (level == GA_MSG_WARNING)
        printf("WARNING: %s", text);
    else
    if (level == GA_MSG_ERROR)
        printf("ERROR: %s", text);
    else
        printf("%s", text);
}

/// Global message function.
GAMessageFunc GA_MSG_FUNC = GA_msg_default;

GAMessageFunc GA_set_msg_func(GAMessageFunc msgfunc)
{
    GA_MSG_FUNC = msgfunc;
    return GA_MSG_FUNC;
}

GAMessageFunc GA_msg()
{
    return GA_MSG_FUNC;
}
