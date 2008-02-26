#ifndef LAP_MAIN
#define LAP_MAIN
/************************************************************************
*
*  lap.h
   version 1.0 - 21 june 1996
   author  Roy Jonker, MagicLogic Optimization Inc.
   
   header file for LAP
*
**************************************************************************/

/** \file lap.h
 * \brief Linear assignment problem solver.
 *
 * This module solves linear assignment problems according to:
 *
 * "A Shortest Augmenting Path Algorithm for Dense and Sparse Linear 
 * Assignment Problems," Computing 38, 325-340, 1987
 * 
 * by
 * 
 * R. Jonker and A. Volgenant, University of Amsterdam.
 *
 * (Author: Roy Jonker, MagicLogic Optimization Inc.)
 */

/* Some changes by Joern P. Meier <mail@ionflux.org>:
   
   (2006-06-18) - Added a namespace prefix to functions.
                - Added doxygen documentation.
                - Added include guards.
                - Added conditional extern "C" tags.
   (2006-07-12) - Changed BIG to something bigger to prevent segfaults 
                  with matrices containing large values.
 */

/*************** CONSTANTS  *******************/

  #define BIG 2147483647

/*************** TYPES      *******************/

  typedef int row;
  typedef int col;
  typedef int cost;

/*************** FUNCTIONS  *******************/

#ifdef __cplusplus
extern "C"
{
#endif

/** Solve linear assignment problem.
 *
 * Solve a linear assignment problem.
 *
 * \param dim problem size
 * \param assigncost cost matrix
 * \param rowsol column assigned to row in solution
 * \param colsol row assigned to column in solution
 * \param u dual variables, row reduction numbers
 * \param v dual variables, column reduction numbers
 */
int LAP_lap(int dim, int **assigncost, 
    int *rowsol, int *colsol, int *u, int *v);

/** Check linear assignment solution.
 *
 * Check a linear assignment solution (?).
 *
 * \param dim problem size
 * \param assigncost cost matrix
 * \param rowsol column assigned to row in solution
 * \param colsol row assigned to column in solution
 * \param u dual variables, row reduction numbers
 * \param v dual variables, column reduction numbers
 */
void LAP_checklap(int dim, int **assigncost,
                     int *rowsol, int *colsol, int *u, int *v);

#ifdef __cplusplus
}
#endif
#endif
