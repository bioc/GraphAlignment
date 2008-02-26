#ifndef LAP_GENERAL
#define LAP_GENERAL
/************************************************************************
*
*  gnrl.h
   version 1.0 - 21 june 1996
   author  Roy Jonker, MagicLogic Optimization Inc.
   
   header file for general stuff
*
**************************************************************************/

/** \file gnrl.h
 * \brief LAP header for general stuff.
 *
 * This is a header file required by the LAP implementation.
 */

/* Some changes by Joern P. Meier <mail@ionflux.org>:
   
   (2006-06-18) - Added doxygen documentation.
                - Added include guards.
 */

/*************** CONSTANTS *******************/

#if !defined TRUE
#define     TRUE        1
#endif
#if !defined FALSE
#define  FALSE        0
#endif

/*************** DATA TYPES *******************/

typedef int boolean;
#endif
