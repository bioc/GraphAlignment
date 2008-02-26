#ifndef GA_GRAPH_ALIGNMENT
#define GA_GRAPH_ALIGNMENT
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
 * Package API functions.
 * ----------------------------------------------------------------------------
 */

/** \file GraphAlignment.h
 * \brief Graph alignment module.
 *
 * This module provides various functions used by the graph alignment R 
 * package, such as the wrapper for the linear assignment problem solver and 
 * the function which computes the score matrix M.
 */

#include "R.h"
#include "Rinternals.h"
#include "Rdefines.h"
#include "R_ext/Rdynload.h"
#include "GA_alloc.h"
#include "GA_message.h"
#include "GA_vector.h"
#include "GA_matrix.h"

#ifdef __cplusplus
extern "C"
{
#endif

/** Module initialization hook.
 *
 * Hook function for module initialization.
 *
 * \param info dynamic link library info
 */
void R_init_GraphAlignment(DllInfo* info);

/** Module unloading hook.
 *
 * Hook function for module unloading.
 *
 * \param info dynamic link library info
 */
void R_unload_GraphAlignment(DllInfo* info);

/** Report a message.
 *
 * Report a message using the appropriate R functions.
 *
 * \param text text of the message
 * \param level message level
 */
void GA_msg_R(const char* text, GAMessageLevel level);

/** Solve linear assignment problem.
 *
 * Solve the linear assignment problem specified by the cost matrix. 
 * The solution returned is a vector of columns assigned to rows. The vector 
 * which is returned will be referenced and should be destroyed by using 
 * GA_vector_destroy_int() when it is not needed anymore.
 *
 * \param costMatrix cost matrix
 *
 * \return Solution of the linear assignment problem
 */
GAVectorInt* GA_linear_assignment_solve(GAMatrixInt* costMatrix);

/** Encode directed graph.
 *
 * Encode an adjacency matrix for a directed graph into a symmetric matrix.
 * 
 * This is done by setting entries (i, j) and (j, i) of the target matrix m' 
 * to 1 if entry m[i, j] = 1 and i > j and to -1 if m[i, j] = 1 and j > i. A 
 * permutation will be applied to i, j if specified as an argument.
 * 
 * The matrix which is returned will be referenced and should be destroyed by 
 * using GA_matrix_destroy_real() when it is not needed anymore.
 * 
 * \param matrix matrix
 * \param p permutation vector
 *
 * \return symmetric matrix
 */
GAMatrixReal* GA_encode_directed_graph(GAMatrixReal* matrix, GAVectorInt* p);

/** Encode directed graph (R).
 *
 * Encode an adjacency matrix for a directed graph into a symmetric matrix.
 * 
 * \param matrix matrix
 * \param p permutation vector
 *
 * \return symmetric matrix
 */
SEXP GA_encode_directed_graph_R(SEXP matrix, SEXP p);

/** Solve linear assignment problem (R).
 *
 * Solve the linear assignment problem specified by the cost matrix. 
 * The solution returned is a vector of rows assigned to columns.
 *
 * \param costMatrix cost matrix
 *
 * \return Solution of the linear assignment problem
 */
SEXP GA_linear_assignment_solve_R(SEXP costMatrix);

/** Directed mode (implementation).
 *
 * The directed mode specifies whether input matrices should be treated as 
 * adjacency matrices of directed graphs.
 */
enum GADirectedMode_Impl
{
    /** Directed mode: enabled.
     */
    GA_DIRECTED_ENABLED = 1,
    /** Directed mode: disabled.
     */
    GA_DIRECTED_DISABLED = 0,
};

/** Directed mode.
 */
typedef enum GADirectedMode_Impl GADirectedMode;

/** Get directed mode from R object (real).
 *
 * Get the directed mode corresponding to the value of the specified R object.
 *
 * \param robj R object
 *
 * \return directed mode
 */
GADirectedMode GA_directed_mode_from_R(SEXP robj);

/** Compute score matrix.
 *
 * Compute the complete score matrix M. The matrix which is returned will be
 * referenced and should be destroyed by using GA_matrix_destroy_real() when 
 * it is not needed anymore.
 *
 * \param na adjacency matrix for network A
 * \param nb adjacency matrix for network B
 * \param r node similarity matrix
 * \param p permutation vector
 * \param linkScore link score matrix
 * \param selfLinkScore self link score matrix
 * \param nodeScore1 node score matrix (1)
 * \param nodeScore2 node score matrix (2)
 * \param lookupLink link bin lookup table
 * \param lookupNode node bin lookup table
 * \param clamp clamp mode for bin lookups
 *
 * \return the score matrix M
 *
 * \sa GA_get_bin_number()
 */
GAMatrixReal* GA_compute_M(GAMatrixReal* na, GAMatrixReal* nb, 
    GAMatrixReal* r, GAVectorInt* p, GAMatrixReal* linkScore, 
    GAMatrixReal* selfLinkScore, GAVectorReal* nodeScore1, 
    GAVectorReal* nodeScore2, GAVectorReal* lookupLink, 
    GAVectorReal* lookupNode, GAClampMode clamp);

/** Compute score matrix (R).
 *
 * Compute the complete score matrix M.
 *
 * \param a adjacency matrix for network A
 * \param b adjacency matrix for network B
 * \param r node similarity matrix
 * \param p permutation vector
 * \param linkScore link score matrix
 * \param selfLinkScore self link score matrix
 * \param nodeScore1 node score matrix (1)
 * \param nodeScore2 node score matrix (2)
 * \param lookupLink link bin lookup table
 * \param lookupNode node bin lookup table
 * \param clamp clamp mode for bin lookups
 *
 * \return the score matrix M
 */
SEXP GA_compute_M_R(SEXP a, SEXP b, SEXP r, SEXP p, SEXP linkScore, 
    SEXP selfLinkScore, SEXP nodeScore1, SEXP nodeScore2, SEXP lookupLink, 
    SEXP lookupNode, SEXP clamp);

#ifdef __cplusplus
}
#endif
#endif
