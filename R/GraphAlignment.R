## ----------------------------------------------------------------------------
## R package for graph alignment
## ----------------------------------------------------------------------------
#
## Author: Joern P. Meier <mail@ionflux.org>
## 
## The package can be used freely for non-commercial purposes. If you use this 
## package, the appropriate paper to cite is J. Berg and M. Laessig, 
## "Cross-species analysis of biological networks by Bayesian alignment", 
## PNAS 103 (29), 10967-10972 (2006)
## 
## This software is made available in the hope that it will be useful, but 
## WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
## or FITNESS FOR A PARTICULAR PURPOSE.
## 
## This software contains code for solving linear assignment problems which was 
## written by Roy Jonker, MagicLogic Optimization Inc.. Please note that this 
## code is copyrighted, (c) 2003 MagicLogic Systems Inc., Canada and may be 
## used for non-commercial purposes only. See 
## http://www.magiclogic.com/assignment.html for the latest version of the LAP 
## code and details on licensing.
#
## ----------------------------------------------------------------------------
## R functions
## ----------------------------------------------------------------------------

.onLoad <- function(libname, pkgname)
{
}

.onAttach <- function(libname, pkgname)
{
}

.onUnload <- function(libpath)
{
}

.Last.lib <- function(libpath)
{
    library.dynam.unload("GraphAlignment", libpath)
}

LinearAssignment <- function(matrix)
{
    .Call("GA_linear_assignment_solve_R", matrix, PACKAGE="GraphAlignment") + 1
}

ComputeM <- function(A, B, R, P, linkScore, selfLinkScore, nodeScore1,
    nodeScore0, lookupLink, lookupNode, clamp=TRUE)
{
    .Call("GA_compute_M_R", A, B, R, P-1, linkScore, selfLinkScore, nodeScore1,
        nodeScore0, lookupLink, lookupNode, clamp, 
    PACKAGE="GraphAlignment")
}

AlignNetworks <- function (A, B, R, P, linkScore, selfLinkScore, nodeScore1,
  nodeScore0, lookupLink, lookupNode, bStart, bEnd, maxNumSteps=2, 
  clamp=TRUE, directed=FALSE)
{
  if (maxNumSteps <= 1)
    stop("[AlignNetworks] Maximum number of steps must be greater than 1.")
  bStep <- (bEnd - bStart)/(maxNumSteps - 1)
  bCur <- bStart

  if (directed)
  {
    ## generate 3x3 link scoring matrices for binary directed graphs
    linkScoreP <- matrix(0,3,3)
    linkScoreP[2:3,2:3] <- linkScore
    linkScoreP[1,1] <- linkScore[1,1]
    linkScoreP[1,2] <- linkScore[2,1]
    linkScoreP[1,3] <- linkScore[1,2]
    linkScoreP[2,1] <- linkScore[1,2]
    linkScoreP[3,1] <- linkScore[2,1]
    
    selfLinkScoreP <- matrix(0,3,3);
    selfLinkScoreP[2:3,2:3] <- selfLinkScore
  }
  
  for (i in 1:maxNumSteps)
  {
    if (directed)
      {
        DA <- EncodeDirectedGraph(A, P)
        DB <- EncodeDirectedGraph(B, P)
      M <- ComputeM(DA, DB, R, P, linkScoreP, selfLinkScoreP, nodeScore1,
        nodeScore0, c(-1.5,-.5,.5,1.5), lookupNode, clamp)
    }
    if (!directed)
      M <- ComputeM(A, B, R, P, linkScore, selfLinkScore, nodeScore1,
        nodeScore0, lookupLink, lookupNode, clamp)
    
    if (bStep != 0)
    {
      S <- matrix(rnorm(length(P)^2)/bCur, length(P), length(P))
      P <- LinearAssignment(round(-1000 * (M/max(abs(M)) + S)))
      bCur <- bCur + bStep
    } else
      P <- LinearAssignment(round(-1000 * (M/max(abs(M)))))
    PM <- M[P, ]
  }
  return(P)
}

InitialAlignment <- function(psize, r=NA, mode="random")
{
  asize <- dim(r)[1]
  bsize <- dim(r)[2]
  if (mode == "random")
  {
    if (!is.na(r) && (psize < min(dim(r))))
      stop("[InitialAlignment] Node similarity matrix R has been specified as well as psize, but psize is too small.")
    rank(rnorm(psize))
  } else
  if (mode == "reciprocal")
  {
    if (is.na(r[1]))
      stop("[InitialAlignment] Node similarity matrix R is required for mode 'reciprocal', but has not been specified.")
    p <- rep(NA, psize)
    pInv <- rep(NA, psize)
    ## Determine aligned nodes.
    for (k in 1:dim(r)[1])
    {
      rowMax <- which.max(r[k,])
      colMax <- which.max(r[,rowMax])
      ## Nodes are aligned if the row/column maxima are unique and 
      ## occur at the same element.
      if ((colMax == k)
          && (length(which(r[k,] == max(r[k,]))) == 1)
          && (length(which(r[,rowMax] == max(r[,rowMax]))) == 1))
      {
        p[k] <- rowMax
        pInv[rowMax] <- k
      }
    }
    ## Align remaining nodes in network A with dummy nodes from network B.
    dummyCountB <- 0
    for (k in 1:asize)
    {
      if (is.na(p[k]))
      {
        if ((bsize + dummyCountB + 1) > psize)
          stop("[InitialAlignment] Not enough dummy nodes in network B (try setting a higher psize).")
        p[k] <- bsize + dummyCountB + 1
        pInv[bsize + dummyCountB + 1] <- k
        dummyCountB <- dummyCountB + 1
      }
    }
    ## Align remaining nodes in network B with dummy nodes from network A.
    dummyCountA <- 0
    for (k in 1:bsize)
    {
      if (is.na(pInv[k]))
      {
        if ((asize + dummyCountA + 1) > psize)
          stop("[InitialAlignment] Not enough dummy nodes in network A (try setting a higher psize).")
        pInv[k] <- asize + dummyCountA + 1
        p[asize + dummyCountA + 1] <- k
        dummyCountA <- dummyCountA + 1
      }
    }
    ## Align remaining dummy nodes.
    l <- 1
    remainingA <- psize - (asize + dummyCountA)
    if (remainingA > 0)
      for (k in (asize + dummyCountA + 1):psize)
      {
        ## Find an unaligned node in B.
        while ((!is.na(pInv[l]))
          && (l < psize))
          l <- l + 1
        if (l > psize)
          stop("[InitialAlignment] Could not align all nodes (this should not happen. Please report a bug!).")
        p[k] <- l
        l <- l + 1
      }
    return(p)
  }
}

CreateScoreMatrix <- function(lookupX, lookupY)
{
  rows <- length(lookupX)
  cols <- length(lookupY)
  s <- matrix(0, rows, cols)
  for (i in 1:rows)
    for (j in 1:cols)
      s[i, j] = lookupX[i] * lookupY[j]
  s
}

InvertPermutation <- function(p)
{
  result <- vector()
  for (i in 1:length(p))
    result[p[i]] <- i
  result
}

Permute <- function(m, p, invertp=FALSE)
{
  if (invertp)
  {
    pInv <- InvertPermutation(p)
    m[pInv,pInv]
  } else
    m[p,p]
}

Trace <- function(m)
{
    sum(diag(m))
}

GetBinNumber <- function(x, lookup, clamp=TRUE)
{
    if (length(lookup) == 0)
        stop("[GetBinNumber] Lookup vector is empty.")
  if (length(lookup) == 1)
  {
    if (!clamp
      && (x != lookup[1]))
    {
      stop("[GetBinNumber] There is only a single lookup value and clamping is disabled, but the input value is not equal to the lookup value. Please make sure you have provided the correct lookup range and clamp mode")
    }
  }
  result = 1
    if ((x < lookup[1])
        || (x > lookup[length(lookup)]))
    {
    if (!clamp)
    {
          stop("[GetBinNumber] Argument is outside of lookup range and clamping is disabled. Please make sure you have provided the correct lookup range and clamp mode")
    }
    if (x < lookup[1])
      result = 1
    else
    if (x > lookup[length(lookup)])
      result = length(lookup) - 1
    } else
  {
      while (((result + 1) <= (length(lookup) - 1))
          && (x >= lookup[result + 1]))
              result <- result + 1
  }
    result
}

VectorToBin <- function(v, lookup, clamp=TRUE)
{
  result <- vector()
  for (i in 1:length(v))
        result[i] = GetBinNumber(v[i], lookup, clamp)
  result
}

MatrixToBin <- function(M, lookup, clamp=TRUE)
{
    result <- matrix(0, dim(M)[1], dim(M)[2])
    for (i  in 1:dim(M)[1])
        for (j in 1:dim(M)[2])
            result[i, j] = GetBinNumber(M[i, j], lookup, clamp)
  result
}

ComputeScores <- function(A, B, R, P, linkScore, selfLinkScore, nodeScore1,
  nodeScore0, lookupLink, lookupNode, symmetric=TRUE, clamp=TRUE)
{
  w <- function(i, j, p, pinv, da, db)
  {
      if ((p[i] <= db) && (pinv[j] <= da))
        0.5
      else
      if ((p[i] <= db) || (pinv[j] <= da))
        1
      else
        0
  }
  aBin<-MatrixToBin(A, lookupLink, clamp)
  bBin<-MatrixToBin(B, lookupLink, clamp)
  rBin<-MatrixToBin(R, lookupNode, clamp)
  dimA<-dim(A)[1]
  dimB<-dim(B)[1]
  pInv<-InvertPermutation(P)
  ## link scores and self link scores
  tsl <- 0
  c0 <- 1
  if (symmetric)
    c0<-0.5
  tsl     <- 0.0
  tslSelf <- 0.0
  ## pairs in alignment
  for (i in 1:dimA) if (P[i] <= dimB)
    for (j in 1:dimA) if ((P[j] <= dimB) && (i != j))
      tsl <- tsl + linkScore[aBin[i, j], bBin[P[i], P[j]]]
  ## nodes in alignment
  for (i in 1:dimA) if (P[i] <= dimB)
    tslSelf <- tslSelf + selfLinkScore[aBin[i, i], bBin[P[i], P[i]]]
  tsl <- c0 * tsl + tslSelf
  tsn <- 0
  for (i in 1:dimA)
    if (P[i] <= dimB)
      {
      ## aligned nodes from network A
        tsn <- tsn + nodeScore1[rBin[i, P[i]]]
      for (j in 1:dimB)
        ## nodes from network B which are not aligned to the current 
        ## node from network A
        if (j != P[i])
          tsn <- tsn + w(i, j, P, pInv, dimA, dimB) * nodeScore0[rBin[i, j]]
      }
  for (j in 1:dimB)
    if (pInv[j] <= dimA)
      {
      ## aligned nodes from network B
      for (i in 1:dimA)
        if (i != pInv[j])
          ## nodes from network A which are not aligned to the current 
          ## node from network B
          tsn <- tsn + w(i, j, P, pInv, dimA, dimB) * nodeScore0[rBin[i, j]]
      }
  list(sl=tsl, sn=tsn)
}

GenerateExample <- function(dimA, dimB, filling, covariance, 
  symmetric=FALSE, numOrths=0, correlated=NA)
{
    if(is.na(correlated[1]))
        correlated <- 1:min(dimA, dimB);

    x=matrix(rnorm(dimA * dimA), dimA, dimA);
    y=matrix(rnorm(dimB * dimB), dimB, dimB);

    ta <- x;
    tb <- y;
    ta[correlated,correlated] <- (sqrt(1 + covariance) 
        * x[correlated,correlated] + sqrt(1 - covariance) 
        * y[correlated,correlated])/sqrt(2);
    tb[correlated,correlated] <- (sqrt(1 + covariance) 
        * x[correlated,correlated] - sqrt(1 - covariance) 
        * y[correlated,correlated])/sqrt(2);

    maskA=matrix(runif(dimA * dimA), dimA, dimA);
    maskB=matrix(runif(dimA * dimA), dimA, dimA);
    maskB[correlated, correlated] <- maskA[correlated, correlated];
    ta[maskA > filling] <- 0;
    tb[maskB > filling] <- 0;

    #enforce symmetry
    if (symmetric) {
        ta[lower.tri(ta)] <- t(ta)[lower.tri(t(ta))]
        tb[lower.tri(tb)] <- t(tb)[lower.tri(t(tb))]
    }

    tr <- matrix(0, dimA, dimB)
    if (numOrths > 0)
    {
        if ((numOrths <= dimA) && (numOrths <= dimB))
            diag(tr)[1:numOrths] <- 1
        else
            diag(tr) <- 1
    }
    list(a=ta, b=tb, r=tr)
}

ComputeLinkParameters <- function(A, B, P, lookupLink, clamp=TRUE)
{
    aBin <- MatrixToBin(A, lookupLink, clamp)
    bBin <- MatrixToBin(B, lookupLink, clamp)
    dimA <- dim(A)[1]
    dimB <- dim(B)[1]
    numBins <- length(lookupLink) - 1
    ## frequency table for self links
    Q0 <- matrix(1, numBins, numBins)
    ## frequency table for other links
    Q1 <- matrix(1, numBins, numBins)
    ## calculate link pair frequencies
    for (i in 1:dimA)
      if (P[i] <= dimB)
      {
          for (j in 1:dimA)
              if (P[j] <= dimB)
              {
                  if (i == j)
                  {
                      ## self link
                      Q0[aBin[i, j], bBin[P[i], P[j]]] <- Q0[aBin[i, j], 
                          bBin[P[i], P[j]]] + 1;
                  } else
                  {
                      ## other link
                      Q1[aBin[i, j], bBin[P[i], P[j]]] <- Q1[aBin[i, j], 
                          bBin[P[i], P[j]]] + 1;
                  }
              }
      }
    Q0 <- Q0 / sum(Q0)
    Q1 <- Q1 / sum(Q1)
    ## calculate marginal distribution vectors
    ## self link score
    pa0 <- vector()
    pb0 <- vector()
    for (i in 1:numBins)
        {
        pa0[i] <- 0
        pb0[i] <- 0
            for (j in 1:numBins)
        {
                    pa0[i] <- pa0[i] + Q0[i, j]
            pb0[i] <- pb0[i] + Q0[j, i]
        }
    }
    ## link score
    pa1 <- vector()
    pb1 <- vector()
    for (i in 1:numBins)
        {
        pa1[i] <- 0
        pb1[i] <- 0
        for (j in 1:numBins)
        {
            pa1[i] <- pa1[i] + Q1[i, j]
            pb1[i] <- pb1[i] + Q1[j, i]
        }
    }
    ## compute link scores
    tSelfLinkScore <- matrix(0, numBins, numBins)
    tLinkScore <- matrix(0, numBins, numBins)
    for (i in 1:numBins)
        for (j in 1:numBins)
        {
            tSelfLinkScore[i, j] <- log(Q0[i, j] / (pa0[i] * pb0[j]))
            tLinkScore[i, j] <- log(Q1[i, j] / (pa1[i] * pb1[j]))
        }
    list(lsSelf=tSelfLinkScore, ls=tLinkScore)
}

ComputeNodeParameters <- function(dimA, dimB, R, P, lookupNode, clamp=TRUE)
{
    rBin <- MatrixToBin(R, lookupNode, clamp)
    q1 <- rep(1, length(lookupNode) - 1)
    q0 <- rep(1, length(lookupNode) - 1)
    p0 <- rep(1, length(lookupNode) - 1)
    pInv <- InvertPermutation(P)
    for (i in 1:dimA)
        for (j in 1:dimB)
        {
            if (P[i] == j)
            {
                ## the pair of nodes (i, j) is aligned
                q1[rBin[i, j]] <- q1[rBin[i, j]] + 1
                p0[rBin[i, j]] <- p0[rBin[i, j]] + 1
            } else
            if ((P[i] <= dimB) || (pInv[j] <= dimA))
            {
                ## at least one of the nodes in (i, j) is aligned, but the 
                ## pair (i, j) is not aligned
                q0[rBin[i, j]] <- q0[rBin[i, j]] + 1
                p0[rBin[i, j]] <- p0[rBin[i, j]] + 1
            }
        }
    ## normalize frequency vectors
    q1 <- q1 / sum(q1)
    q0 <- q0 / sum(q0)
    p0 <- p0 / sum(p0)
    ts1 <- log(q1 / p0)
    ts0 <- log(q0 / p0)
    list(s0=ts0, s1=ts1)
}

AnalyzeAlignment <- function(A, B, R, P, lookupNode, epsilon=0, clamp=TRUE)
{
    dimA <- dim(A)[1]
    dimB <- dim(B)[1]
    rBin<-MatrixToBin(R, lookupNode, clamp)
    tna <- 0
    tnb <- 0
    tnc <- 0
    for (i in 1:dimA)
    {
        j <- P[i]
        if (j <= dimB)
        {
            ## tna: number of aligned node pairs.
            tna <- tna + 1
            ## tnb: number of aligned node pairs (ia, ib), where no jb exists 
            ## such that R[ia, jb] > epsilon and no ja such that R[ja, ib] 
            ## > epsilon
            checkB <- TRUE
            for (k in 1:dimB)
                if (R[i, k] > epsilon)
                    checkB <- FALSE
            checkA <- TRUE
            for (k in 1:dimA)
                if (R[k, j] > epsilon)
                    checkA <- FALSE
            if (checkA && checkB)
                tnb <- tnb + 1
            ## tnc: number of aligned node pairs (ia, ib) with R[ia, ib] 
            ## < epsilon but jb or ja exists, such that R[ia, jb] > epsilon or 
            ## R[ja, ib] > epsilon
            if ((R[i, j] < epsilon)
                && (!checkA || !checkB))
                tnc <- tnc + 1
        }
    }
    list(na=tna, nb=tnb, nc=tnc)
}

AlignedPairs <- function(A, B, P)
{
    dimA <- dim(A)[1]
    dimB <- dim(B)[1]
    tna <- 0
    for (i in 1:dimA)
            if ((i <= length(P)) && (P[i] <= dimB))
            tna <- tna + 1
    alignedPairs <- matrix(0, tna, 2)
    pos <- 1
    for (i in 1:dimA)
            if ((i <= length(P)) && (P[i] <= dimB))
        {
            alignedPairs[pos, 1] <- i
            alignedPairs[pos, 2] <- P[i]
            pos <- pos + 1
        }
    alignedPairs
}

EncodeDirectedGraph <- function(matrix, P)
{
    .Call("GA_encode_directed_graph_R", matrix, P-1, PACKAGE="GraphAlignment")
}
