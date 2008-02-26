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

/** \file GraphAlignment.c
 * \brief Graph alignment module (implementation).
 */

#include "GraphAlignment.h"
#include "GA_vector_R.h"
#include "GA_matrix_R.h"
#include "lap.h"

void GA_msg_R(const char* text, GAMessageLevel level)
{
    if (level == GA_MSG_DEBUG)
        Rprintf("DEBUG: %s", text);
    else
    if (level == GA_MSG_WARNING)
        Rprintf("WARNING: %s", text);
    else
    if (level == GA_MSG_ERROR)
        error("ERROR: %s", text);
    else
        Rprintf("%s", text);
}

GAVectorInt* GA_linear_assignment_solve(GAMatrixInt* costMatrix)
{
    if (costMatrix->rows != costMatrix->cols)
    {
        GA_msg()("[GA_linear_assignment_solve] "
            "Cost matrix is not a square matrix", GA_MSG_ERROR);
        return 0;
    }
    /* ----- DEBUG ----- //
    GA_msg()("[GA_linear_assignment_solve_R] "
        "Creating result vectors.\n", 
        GA_MSG_DEBUG);
    // ----- DEBUG ----- */
    GAVectorInt* resultRows = GA_vector_create_int(costMatrix->rows);
    if (resultRows == 0)
        return 0;
    GA_vector_init_zero_int(resultRows);
    GAVectorInt* resultCols = GA_vector_create_int(costMatrix->rows);
    if (resultCols == 0)
        return 0;
    GA_vector_init_zero_int(resultCols);
    GAVectorInt* u = GA_vector_create_int(costMatrix->rows);
    if (u == 0)
        return 0;
    GA_vector_init_zero_int(u);
    GAVectorInt* v = GA_vector_create_int(costMatrix->rows);
    if (v == 0)
        return 0;
    GA_vector_init_zero_int(v);
    /* ----- DEBUG ----- //
    GA_msg()("[GA_linear_assignment_solve_R] "
        "Calling LAP:\n", 
        GA_MSG_DEBUG);
    GA_msg()("[GA_linear_assignment_solve_R] cost matrix:\n", GA_MSG_DEBUG);
    GA_matrix_print_int(costMatrix);
    GA_msg()("\n", GA_MSG_INFO);
    GA_msg()("[GA_linear_assignment_solve_R] resultRows = ", GA_MSG_DEBUG);
    GA_vector_print_int(resultRows);
    GA_msg()("\n", GA_MSG_INFO);
    GA_msg()("[GA_linear_assignment_solve_R] resultCols = ", GA_MSG_DEBUG);
    GA_vector_print_int(resultCols);
    GA_msg()("\n", GA_MSG_INFO);
    GA_msg()("[GA_linear_assignment_solve_R] u = ", GA_MSG_DEBUG);
    GA_vector_print_int(u);
    GA_msg()("\n", GA_MSG_INFO);
    GA_msg()("[GA_linear_assignment_solve_R] v = ", GA_MSG_DEBUG);
    GA_vector_print_int(v);
    GA_msg()("\n", GA_MSG_INFO);
    char* message = GA_alloc(512, sizeof(char));
    snprintf(message, 512, "[GA_linear_assignment_solve_R] Input: "
        "costMatrix->rows = %i, costMatrix->elts = %p, "
        "resultRows->elts = %p, resultCols->elts = %p, "
        "u->elts = %p, v->elts = %p\n", costMatrix->rows, 
        costMatrix->elts, resultRows->elts, resultCols->elts, 
        u->elts, v->elts);
    GA_msg()(message, GA_MSG_DEBUG);
    GA_free(message);
    // ----- DEBUG ----- */
    LAP_lap(costMatrix->rows, costMatrix->elts, 
        resultRows->elts, resultCols->elts, u->elts, v->elts);
    /* ----- DEBUG ----- //
    GA_msg()("[GA_linear_assignment_solve_R] "
        "Destroying result vectors.\n", 
        GA_MSG_DEBUG);
    // ----- DEBUG ----- */
    GA_vector_destroy_int(resultRows);
    GA_vector_destroy_int(u);
    GA_vector_destroy_int(v);
    return resultCols;
}

GAMatrixReal* GA_encode_directed_graph(GAMatrixReal* matrix, GAVectorInt* p)
{
    if (matrix->rows != matrix->cols)
    {
        GA_msg()("[GA_encode_directed_graph] "
            "Input matrix is not a square matrix.", GA_MSG_ERROR);
        return 0;
    }
    GAMatrixReal* result = GA_matrix_create_square_real(matrix->rows);
    GA_matrix_init_zero_real(result);
    int i;
    int j;
    if (p == 0)
    {
        for (i = 0; i < result->rows; i++)
            for (j = 0; j < result->rows; j++)
                if (matrix->elts[i][j] == 1)
                    if (i <= j)
                    {
                        result->elts[i][j] = 1;
                        result->elts[j][i] = 1;
                    } else
                    {
                        result->elts[i][j] = -1;
                        result->elts[j][i] = -1;
                    }
    } else
    {
        if (result->rows > p->size)
        {
            GA_msg()("[GA_encode_directed_graph] "
                "Not enough elements in the permutation vector.", 
                    GA_MSG_ERROR);
            return 0;
        }
        for (i = 0; i < result->rows; i++)
            for (j = 0; j < result->rows; j++)
                if (matrix->elts[i][j] == 1)
                    if (p->elts[i] <= p->elts[j])
                    {
                        result->elts[i][j] = 1;
                        result->elts[j][i] = 1;
                    } else
                    {
                        result->elts[i][j] = -1;
                        result->elts[j][i] = -1;
                    }
    }
    return result;
}

SEXP GA_encode_directed_graph_R(SEXP matrix, SEXP p)
{
    PROTECT(matrix);
    PROTECT(p);
    GAMatrixReal* gaMatrix = GA_matrix_from_R_real(matrix);
    if (gaMatrix == 0)
    {
        UNPROTECT(2);
        return R_NilValue;
    }
    GAVectorInt* gaP = GA_vector_from_R_int(p);
    if (gaP == 0)
    {
        UNPROTECT(2);
        return R_NilValue;
    }
    GAMatrixReal* gaResult = GA_encode_directed_graph(gaMatrix, gaP);
    GA_matrix_destroy_real(gaMatrix);
    GA_vector_destroy_int(gaP);
    if (gaResult == 0)
    {
        UNPROTECT(2);
        return R_NilValue;
    }
    SEXP result;
    result = GA_matrix_to_R_real(gaResult);
    GA_matrix_destroy_real(gaResult);
    UNPROTECT(2);
    return result;
}

SEXP GA_linear_assignment_solve_R(SEXP costMatrix)
{
    PROTECT(costMatrix);
    /* ----- DEBUG ----- //
    GA_msg()("[GA_linear_assignment_solve_R] "
        "Creating cost matrix.\n", 
        GA_MSG_DEBUG);
    // ----- DEBUG ----- */
    GAMatrixInt* gaCostMatrix = GA_matrix_from_R_int(costMatrix);
    if (gaCostMatrix == 0)
    {
        UNPROTECT(1);
        return R_NilValue;
    }
    /* ----- DEBUG ----- //
    GA_msg()("[GA_linear_assignment_solve_R] "
        "Solving linear assignment problem.\n", 
        GA_MSG_DEBUG);
    // ----- DEBUG ----- */
    GAVectorInt* gaResult = GA_linear_assignment_solve(gaCostMatrix);
    /* ----- DEBUG ----- //
    GA_msg()("[GA_linear_assignment_solve_R] "
        "Destroying cost matrix.\n", 
        GA_MSG_DEBUG);
    // ----- DEBUG ----- */
    GA_matrix_destroy_int(gaCostMatrix);
    if (gaResult == 0)
    {
        UNPROTECT(1);
        return R_NilValue;
    }
    SEXP result;
    /* ----- DEBUG ----- //
    GA_msg()("[GA_linear_assignment_solve_R] "
        "Creating result vector.\n", 
        GA_MSG_DEBUG);
    // ----- DEBUG ----- */
    result = GA_vector_to_R_int(gaResult);
    GA_vector_destroy_int(gaResult);
    UNPROTECT(1);
    return result;
}

GADirectedMode GA_directed_mode_from_R(SEXP robj)
{
    PROTECT(robj);
    SEXPTYPE objType = TYPEOF(robj);
    if ((objType != LGLSXP)
        && (objType != INTSXP)
        && (objType != REALSXP))
    {
        char* message = GA_alloc(256, sizeof(char));
        snprintf(message, 256, 
            "[GA_directed_mode_from_R] Input is not a logical, real or "
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
        return GA_DIRECTED_DISABLED;
    }
    UNPROTECT(1);
    return GA_DIRECTED_ENABLED;
}

GAMatrixReal* GA_compute_M(GAMatrixReal* a, GAMatrixReal* b, 
    GAMatrixReal* r, GAVectorInt* p, GAMatrixReal* linkScore, 
    GAMatrixReal* selfLinkScore, GAVectorReal* nodeScore1, 
    GAVectorReal* nodeScore2, GAVectorReal* lookupLink, 
    GAVectorReal* lookupNode, GAClampMode clamp)
{
    GAVectorInt* pInv = GA_invert_permutation_int(p);
    if (pInv == 0)
    {
        GA_msg()("[GA_compute_M] "
            "Inverted permutation is null.", GA_MSG_ERROR);
        return 0;
    }
    /* Various sanity checks of input values. */
    if (a->rows != a->cols)
    {
        GA_msg()("[GA_compute_M] "
            "Adjacency matrix for network A is not a square matrix.", 
                GA_MSG_ERROR);
        return 0;
    }
    if (b->rows != b->cols)
    {
        GA_msg()("[GA_compute_M] "
            "Adjacency matrix for network B is not a square matrix.", 
                GA_MSG_ERROR);
        return 0;
    }
    /* More sanity checks. */
    if ((r->rows != a->rows) 
        || (r->cols != b->rows))
    {
        char* message = GA_alloc(256, sizeof(char));
        snprintf(message, 256, "[GA_compute_M] "
            "Node similarity matrix R has wrong dimensions (%i, %i) "
            "(expected (%i, %i)).", r->rows, r->cols, a->rows, b->rows);
        GA_msg()(message, GA_MSG_ERROR);
        GA_free(message);
        GA_vector_destroy_int(pInv);
        return 0;
    }
    if ((linkScore->rows < (lookupLink->size - 1))
        || (linkScore->cols < (lookupLink->size - 1)))
    {
        char* message = GA_alloc(256, sizeof(char));
        snprintf(message, 256, "[GA_compute_M] "
            "Link score matrix dimension does not match number of bins "
            "(dim(linkScore) = (%i, %i), length(lookupLink) = %i).", 
            linkScore->rows, linkScore->cols, lookupLink->size);
        GA_msg()(message, GA_MSG_ERROR);
        GA_free(message);
        GA_vector_destroy_int(pInv);
        return 0;
    }
    if ((selfLinkScore->rows < (lookupLink->size - 1))
        || (selfLinkScore->cols < (lookupLink->size - 1)))
    {
        char* message = GA_alloc(256, sizeof(char));
        snprintf(message, 256, "[GA_compute_M] "
            "Self link score matrix dimension does not match number of bins "
            "(dim(linkScore) = (%i, %i), length(lookupLink) = %i).", 
            linkScore->rows, linkScore->cols, lookupLink->size);
        GA_msg()(message, GA_MSG_ERROR);
        GA_free(message);
        GA_vector_destroy_int(pInv);
        return 0;
    }
    if (nodeScore1->size < (lookupNode->size - 1))
    {
        char* message = GA_alloc(256, sizeof(char));
        snprintf(message, 256, "[GA_compute_M] "
            "Length of node score vector (s1) does not match number of bins "
            "(length(nodeScore1) = %i, length(lookupNode) = %i).", 
            nodeScore1->size, lookupNode->size);
        GA_msg()(message, GA_MSG_ERROR);
        GA_free(message);
        GA_vector_destroy_int(pInv);
        return 0;
    }
    if (nodeScore2->size < (lookupNode->size - 1))
    {
        char* message = GA_alloc(256, sizeof(char));
        snprintf(message, 256, "[GA_compute_M] "
            "Length of node score vector (s2) does not match number of bins "
            "(length(nodeScore2) = %i, length(lookupNode) = %i).", 
            nodeScore2->size, lookupNode->size);
        GA_msg()(message, GA_MSG_ERROR);
        GA_free(message);
        GA_vector_destroy_int(pInv);
        return 0;
    }
    /* ----- DEBUG ----- //
    GA_msg()("[GA_compute_M] p:\n", GA_MSG_DEBUG);
    GA_vector_print_int(p);
    GA_msg()("\n", GA_MSG_INFO);
    GA_msg()("[GA_compute_M] Binning matrices...\n", GA_MSG_DEBUG);
    // ----- DEBUG ----- */
    GAMatrixInt* aBin = GA_matrix_to_bin_real(a, lookupLink, clamp);
    if (aBin == 0)
        return 0;
    GAMatrixInt* bBin = GA_matrix_to_bin_real(b, lookupLink, clamp);
    if (bBin == 0)
        return 0;
    GAMatrixInt* rBin = GA_matrix_to_bin_real(r, lookupNode, clamp);
    if (rBin == 0)
        return 0;
    /* ----- DEBUG ----- //
    GA_msg()("[GA_compute_M] Done.\n", GA_MSG_DEBUG);
    GA_msg()("[GA_compute_M] Calculating M...\n", GA_MSG_DEBUG);
    GAMatrixReal* linkScoreResult = GA_matrix_create_square_real(p->size);
    GAMatrixReal* selfLinkScoreResult = GA_matrix_create_square_real(p->size);
    GAMatrixReal* nodeScoreResult = GA_matrix_create_square_real(p->size);
    GA_matrix_init_zero_real(linkScoreResult);
    GA_matrix_init_zero_real(selfLinkScoreResult);
    GA_matrix_init_zero_real(nodeScoreResult);
    // ----- DEBUG ----- */
    GAMatrixReal* result = GA_matrix_create_square_real(p->size);
    GA_matrix_init_zero_real(result);
    int i;
    int j;
    for (i = 0; i < p->size; i++)
    {
        for (j = 0; j < p->size; j++)
        {
            int k;
            /* Sum up link scores. */
            /* ----- DEBUG ----- //
            char* message = GA_alloc(256, sizeof(char));
            snprintf(message, 256, "[GA_compute_M] "
                "Calculating link scores for M(%i, %i)...\n", i, j);
            GA_msg()(message, GA_MSG_DEBUG);
            // ----- DEBUG ----- */
            double linkScoreSum = 0.0;
            if ((i < b->rows)
                && (j < a->rows))
                for (k = 0; k < a->rows; k++)
                {
                    if ((k != j)
                        && (p->elts[k] != i)
                        && (p->elts[k] < b->rows))
                    {
                        /* ----- DEBUG ----- //
                        snprintf(message, 256, "[GA_compute_M] "
                            "(i, j, k, p[k], aBin->elts[j][k], "
                            "bBin->elts[i][p->elts[k]]) = "
                            "(%i, %i, %i, %i, %i, %i)\n", 
                            i, j, k, p->elts[k], aBin->elts[j][k], 
                            bBin->elts[i][p->elts[k]]);
                        GA_msg()(message, GA_MSG_DEBUG);
                        // ----- DEBUG ----- */
                        linkScoreSum += linkScore->elts[aBin->elts[j][k]]
                            [bBin->elts[i][p->elts[k]]];
                    }
                }
            /* Sum up self link scores. */
            /* ----- DEBUG ----- //
            snprintf(message, 256, "[GA_compute_M] "
                "Calculating self link scores for M(%i, %i)...\n", i, j);
            GA_msg()(message, GA_MSG_DEBUG);
            // ----- DEBUG ----- */
            double selfLinkScoreSum = 0.0;
            if ((i < b->rows)
                && (j < a->rows))
            {
                selfLinkScoreSum += selfLinkScore->elts[aBin->elts[j][j]]
                    [bBin->elts[i][i]];
            }
            /* Sum up node similarity scores. */
            /* ----- DEBUG ----- //
            snprintf(message, 256, "[GA_compute_M] "
                "Calculating node scores for M(%i, %i)...\n", i, j);
            GA_msg()(message, GA_MSG_DEBUG);
            // ----- DEBUG ----- */
            double nodeScoreSum = 0.0;
            if ((i < b->rows)
                && (j < a->rows))
            {
                nodeScoreSum += nodeScore1->elts[rBin->elts[j][i]];
                for (k = 0; k < a->rows; k++)
                {
                    if ((p->elts[k] >= b->rows)
                        && (k != j))
                        nodeScoreSum += nodeScore2->elts[rBin->elts[k][i]];
                }
                for (k = 0; k < b->rows; k++)
                {
                    if ((pInv->elts[k] >= a->rows)
                        && (k != i))
                        nodeScoreSum += nodeScore2->elts[rBin->elts[j][k]];
                }
            }
            /* Set the element of M. */
            result->elts[i][j] = linkScoreSum + selfLinkScoreSum 
                + nodeScoreSum;
            /* ----- DEBUG ----- //
            linkScoreResult->elts[i][j] = linkScoreSum;
            selfLinkScoreResult->elts[i][j] = selfLinkScoreSum;
            nodeScoreResult->elts[i][j] = nodeScoreSum;
            // ----- DEBUG ----- */
        }
    }
    /* ----- DEBUG ----- //
    GA_msg()("[GA_compute_M] Link score result:\n", GA_MSG_DEBUG);
    GA_matrix_print_real(linkScoreResult);
    GA_msg()("\n", GA_MSG_INFO);
    GA_msg()("[GA_compute_M] Self link score result:\n", GA_MSG_DEBUG);
    GA_matrix_print_real(selfLinkScoreResult);
    GA_msg()("\n", GA_MSG_INFO);
    GA_msg()("[GA_compute_M] Node score result:\n", GA_MSG_DEBUG);
    GA_matrix_print_real(nodeScoreResult);
    GA_msg()("\n", GA_MSG_INFO);
    GA_matrix_destroy_real(linkScoreResult);
    GA_matrix_destroy_real(selfLinkScoreResult);
    GA_matrix_destroy_real(nodeScoreResult);
    GA_msg()("[GA_compute_M] Done.\n", GA_MSG_DEBUG);
    // ----- DEBUG ----- */
    GA_vector_destroy_int(pInv);
    GA_matrix_destroy_int(aBin);
    GA_matrix_destroy_int(bBin);
    GA_matrix_destroy_int(rBin);
    return result;
}

SEXP GA_compute_M_R(SEXP a, SEXP b, SEXP r, SEXP p, SEXP linkScore, 
    SEXP selfLinkScore, SEXP nodeScore1, SEXP nodeScore2, SEXP lookupLink, 
    SEXP lookupNode, SEXP clamp)
{
    PROTECT(a);
    PROTECT(b);
    PROTECT(r);
    PROTECT(p);
    PROTECT(linkScore);
    PROTECT(selfLinkScore);
    PROTECT(nodeScore1);
    PROTECT(nodeScore2);
    PROTECT(lookupLink);
    PROTECT(lookupNode);
    PROTECT(clamp);
    static const int numArgs = 11;
    GAMatrixReal* gaA = GA_matrix_from_R_real(a);
    if (gaA == 0)
    {
        UNPROTECT(numArgs);
        return R_NilValue;
    }
    GAMatrixReal* gaB = GA_matrix_from_R_real(b);
    if (gaB == 0)
    {
        GA_matrix_destroy_real(gaA);
        UNPROTECT(numArgs);
        return R_NilValue;
    }
    GAMatrixReal* gaR = GA_matrix_from_R_real(r);
    if (gaR == 0)
    {
        GA_matrix_destroy_real(gaA);
        GA_matrix_destroy_real(gaB);
        UNPROTECT(numArgs);
        return R_NilValue;
    }
    GAVectorInt* gaP = GA_vector_from_R_int(p);
    if (gaP == 0)
    {
        GA_matrix_destroy_real(gaA);
        GA_matrix_destroy_real(gaB);
        GA_matrix_destroy_real(gaR);
        UNPROTECT(numArgs);
        return R_NilValue;
    }
    GAMatrixReal* gaLinkScore = GA_matrix_from_R_real(linkScore);
    if (gaLinkScore == 0)
    {
        GA_matrix_destroy_real(gaA);
        GA_matrix_destroy_real(gaB);
        GA_matrix_destroy_real(gaR);
        GA_vector_destroy_int(gaP);
        UNPROTECT(numArgs);
        return R_NilValue;
    }
    GAMatrixReal* gaSelfLinkScore = GA_matrix_from_R_real(selfLinkScore);
    if (gaSelfLinkScore == 0)
    {
        GA_matrix_destroy_real(gaA);
        GA_matrix_destroy_real(gaB);
        GA_matrix_destroy_real(gaR);
        GA_vector_destroy_int(gaP);
        GA_matrix_destroy_real(gaLinkScore);
        UNPROTECT(numArgs);
        return R_NilValue;
    }
    GAVectorReal* gaNodeScore1 = GA_vector_from_R_real(nodeScore1);
    if (gaNodeScore1 == 0)
    {
        GA_matrix_destroy_real(gaA);
        GA_matrix_destroy_real(gaB);
        GA_matrix_destroy_real(gaR);
        GA_vector_destroy_int(gaP);
        GA_matrix_destroy_real(gaLinkScore);
        GA_matrix_destroy_real(gaSelfLinkScore);
        UNPROTECT(numArgs);
        return R_NilValue;
    }
    GAVectorReal* gaNodeScore2 = GA_vector_from_R_real(nodeScore2);
    if (gaNodeScore2 == 0)
    {
        GA_matrix_destroy_real(gaA);
        GA_matrix_destroy_real(gaB);
        GA_matrix_destroy_real(gaR);
        GA_vector_destroy_int(gaP);
        GA_matrix_destroy_real(gaLinkScore);
        GA_matrix_destroy_real(gaSelfLinkScore);
        GA_vector_destroy_real(gaNodeScore1);
        UNPROTECT(numArgs);
        return R_NilValue;
    }
    GAVectorReal* gaLookupLink = GA_vector_from_R_real(lookupLink);
    if (gaLookupLink == 0)
    {
        GA_matrix_destroy_real(gaA);
        GA_matrix_destroy_real(gaB);
        GA_matrix_destroy_real(gaR);
        GA_vector_destroy_int(gaP);
        GA_matrix_destroy_real(gaLinkScore);
        GA_matrix_destroy_real(gaSelfLinkScore);
        GA_vector_destroy_real(gaNodeScore1);
        GA_vector_destroy_real(gaNodeScore2);
        UNPROTECT(numArgs);
        return R_NilValue;
    }
    GAVectorReal* gaLookupNode = GA_vector_from_R_real(lookupNode);
    if (gaLookupNode == 0)
    {
        GA_matrix_destroy_real(gaA);
        GA_matrix_destroy_real(gaB);
        GA_matrix_destroy_real(gaR);
        GA_vector_destroy_int(gaP);
        GA_matrix_destroy_real(gaLinkScore);
        GA_matrix_destroy_real(gaSelfLinkScore);
        GA_vector_destroy_real(gaNodeScore1);
        GA_vector_destroy_real(gaNodeScore2);
        GA_vector_destroy_real(gaLookupLink);
        UNPROTECT(numArgs);
        return R_NilValue;
    }
    GAClampMode gaClamp = GA_clamp_mode_from_R(clamp);
    GAMatrixReal* gaResult = GA_compute_M(gaA, gaB, gaR, gaP, gaLinkScore, 
        gaSelfLinkScore, gaNodeScore1, gaNodeScore2, gaLookupLink,
        gaLookupNode, gaClamp);
    SEXP result = GA_matrix_to_R_real(gaResult);
    GA_matrix_destroy_real(gaResult);
    GA_matrix_destroy_real(gaA);
    GA_matrix_destroy_real(gaB);
    GA_matrix_destroy_real(gaR);
    GA_vector_destroy_int(gaP);
    GA_matrix_destroy_real(gaLinkScore);
    GA_matrix_destroy_real(gaSelfLinkScore);
    GA_vector_destroy_real(gaNodeScore1);
    GA_vector_destroy_real(gaNodeScore2);
    GA_vector_destroy_real(gaLookupLink);
    GA_vector_destroy_real(gaLookupNode);
    UNPROTECT(numArgs);
    return result;
}

/** Methods used with the Call interface.
 */
R_CallMethodDef GA_callMethods[] = {
    {
        "GA_linear_assignment_solve_R",
        (DL_FUNC)&GA_linear_assignment_solve_R,
        1
    },
    {
        "GA_compute_M_R",
        (DL_FUNC)&GA_compute_M_R,
        11
    },
    {
        "GA_encode_directed_graph_R",
        (DL_FUNC)&GA_encode_directed_graph_R,
        2
    },
    {
        NULL,
        NULL,
        0
    }
};

/** Free function which does nothing.
 */
void GA_free_dummy(char* memLoc)
{
    /* Nothing */
}

void R_init_GraphAlignment(DllInfo* info)
{
    /* ----- DEBUG ----- //
    Rprintf("[R_init_GraphAlignment] DEBUG: "
        "Initialized GraphAlignment module.\n");
    // ----- DEBUG ----- */
    /* Use the R API for allocating memory. */
    GA_set_alloc_funcs(R_alloc, GA_free_dummy);
    /* Use the R API for printing messages. */
    GA_set_msg_func(GA_msg_R);
    R_registerRoutines(info, NULL, GA_callMethods, NULL, NULL);
}

void R_unload_GraphAlignment(DllInfo* info)
{
    /* ----- DEBUG ----- //
    Rprintf("[R_unload_GraphAlignment] DEBUG: "
        "Unloaded GraphAlignment module.\n");
    // ----- DEBUG ----- */
}
