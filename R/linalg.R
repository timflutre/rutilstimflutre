## Contains functions useful for linear algebra.

##' LDL' decomposition
##'
##' This wrapper around \code{\link[base]{chol}} factors a square, symmetric, positive (semi-)definite matrix into the product of a lower triangular matrix (L), a diagonal matrix (D), and an upper triangular matrix (U) equal to the transpose of the lower triangular one (U=L').
##' It is also known as the square-root-free Cholesky decomposition.
##' This decomposition is notably used in quantitative genetics (Misztal and Perez-Enciso, 1993).
##' The determinant of the input matrix is then easily computed as the product of the diagonal entries of the decomposition.
##' @param x matrix
##' @param returnFormat a lower triangular matrix with the elements of D on the diagonal, or a list with two components, L as a matrix with 1 on the diagonal and the diagonal of D as a vector
##' @return matrix
##' @seealso \code{\link{getSparseInv}}
##' @examples
##' \dontrun{## example matrix from Misztal and Perez-Enciso (1993):
##' W <- matrix(NA, nrow=5, ncol=5)
##' W[1,] <- c(5, 0, 0, 3, 0)
##' W[2,-1] <- c(3, 0, 0, 1)
##' W[3,-c(1:2)] <- c(4, 1, 1)
##' W[4,-c(1:3)] <- c(4, 0)
##' W[5,-c(1:4)] <- 2
##' W
##' W[lower.tri(W)] <- t(W)[lower.tri(W)]
##' W
##'
##' ## decomposition/factorization:
##' (decompM <- LDLT(W, returnFormat="matrix"))
##' (decompL <- LDLT(W, returnFormat="list"))
##'
##' ## determinant:
##' prod(diag(decompM)) # ~= 162, as in the paper
##' prod(decompL$vecD)
##' }
##' @export
LDLT <- function(x, returnFormat="matrix"){
  stopifnot(nrow(x) == ncol(x),
            returnFormat %in% c("matrix","list"))

  ## chol() returns the upper triangular factor of the Cholesky decomposition
  ## i.e., R such that R^T R = x
  ## let L_c = R^T
  ## moreover, D' = diagonal matrix made from diag(L_c)
  ## we can also write L_c = L D'; thus L = L_c D'^-1
  ## we have x = L_c L_c^T = (L D') (L D')^T = L D'^2 L^T = L D L^T
  ## where D = D'^2
  L_c <- t(chol(x))
  diag_L_c <- diag(L_c)
  L <- L_c %*% diag(1 / diag_L_c)
  vecD <- diag_L_c^2

  if(returnFormat == "matrix"){
    out <- L
    diag(out) <- vecD
  } else if(returnFormat == "list")
    out <- list(L=L, vecD=vecD)

  return(out)
}

##' Entries of a sparse inverse
##'
##' Efficiently computes the entries of the sparse inverse (C) of a square, symmetric matrix (W) from its square-root-free Cholesky decomposition (W = L D L') as implemented in \code{\link{LDLT}}.
##' This is notably used in quantitative genetics as described initially by \href{http://doi.org/10.3168/jds.S0022-0302(93)77478-0}{Misztal and Perez-Enciso (1993)}.
##' See the section 15.7 "Computing entries of the inverse of a sparse matrix" from Duff et al (2017) for the correct equation.
##' @param vecD vector
##' @param L lower triangular matrix
##' @param reprSparse if not NULL, the output matrix will be a \code{\link[Matrix]{sparseMatrix}} of the given representation ("C", "R" or "T")
##' @param verbose verbosity level
##' @return matrix
##' @author Timothee Flutre
##' @seealso \code{\link{LDLT}}
##' @examples
##' \dontrun{## example matrix from Misztal and Perez-Enciso (1993):
##' W <- matrix(NA, nrow=5, ncol=5)
##' W[1,] <- c(5, 0, 0, 3, 0)
##' W[2,-1] <- c(3, 0, 0, 1)
##' W[3,-c(1:2)] <- c(4, 1, 1)
##' W[4,-c(1:3)] <- c(4, 0)
##' W[5,-c(1:4)] <- 2
##' W
##' W[lower.tri(W)] <- t(W)[lower.tri(W)]
##' W
##'
##' ## LDLT decomposition (a.k.a. square-root-free Cholesky):
##' decomp <- LDLT(W, returnFormat="list")
##'
##' ## C, inverse of W:
##' getSparseInv(decomp$vecD, decomp$L)
##' getSparseInv(decomp$vecD, decomp$L, "T")
##' }
##' @export
getSparseInv <- function(vecD, L, reprSparse=NULL, verbose=FALSE){
  stopifnot(is.vector(vecD),
            is.matrix(L))
  if(! is.null(reprSparse))
    stopifnot(reprSparse %in% c("C","R","T"))

  U <- t(L)
  n <- nrow(U)
  if(! is.null(reprSparse)){
    vec_i <- vector("integer")
    vec_j <- vector("integer")
    vec_x <- vector("numeric")
  } else
    C <- matrix(NA, nrow=n, ncol=n)

  for(j in n:1){
    if(verbose) cat(paste0("j=",j,"\n"))

    ## diagonal value:
    i <- j
    if(verbose) cat(paste0("i=",i, " -> C_ii\n"))
    tmp <- rep(0, n-(i+1)+1)
    if(i < n){
      for(k in (i+1):n){
        if(U[i,k] != 0){
          if(! is.null(reprSparse)){
            idx <- which(vec_i == i & vec_j == k)
            stopifnot(length(idx) == 1)
            C_ik <- vec_x[idx]
          } else
            C_ik <- C[i,k]
          tmp[k-i] <- U[i,k] * C_ik
        }
      }
    }
    if(verbose > 1) print(tmp)
    C_ii <- 1 / vecD[i] - sum(tmp)
    if(! is.null(reprSparse)){
      vec_i <- c(vec_i, i)
      vec_j <- c(vec_j, j)
      vec_x <- c(vec_x, C_ii)
    } else
      C[i,i] <- C_ii

    ## upper-triangular values:
    for(i in (j-1):1){
      if(i == 0)
        break
      if(verbose) cat(paste0("i=",i))
      if(i < j){
        if(U[i,j] == 0){
          if(verbose) cat(paste0(" -> U[",i,",",j,"]=0\n"))
        } else{
          if(verbose) cat(" -> C_ij\n")
          tmp <- rep(0, n-(i+1)+1)
          for(k in n:(i+1)){
            if(k <= j){
              if(U[k,j] != 0){
                if(! is.null(reprSparse)){
                  idx <- which(vec_i == k & vec_j == j)
                  stopifnot(length(idx) == 1)
                  C_kj <- vec_x[idx]
                } else
                  C_kj <- C[k,j]
                tmp[n-k+1] <- U[i,k] * C_kj # see Duff et al (2017), eq.15.7.5
              }
            } else{
              if(U[j,k] != 0){
                if(! is.null(reprSparse)){
                  idx <- which(vec_i == j & vec_j == k)
                  stopifnot(length(idx) == 1)
                  C_jk <- vec_x[idx]
                } else
                  C_jk <- C[j,k]
                tmp[n-k+1] <- U[i,k] * C_jk
              }
            }
          }
          if(verbose > 1) print(tmp)
          C_ij <- - sum(tmp)
          if(! is.null(reprSparse)){
            vec_i <- c(vec_i, i)
            vec_j <- c(vec_j, j)
            vec_x <- c(vec_x, C_ij)
          } else
            C[i,j] <- C_ij
        }
      }
    }
  }

  if(! is.null(reprSparse))
    C <- Matrix::sparseMatrix(i=vec_i, j=vec_j, x=vec_x, dims=c(n, n),
                              symmetric=TRUE, index1=TRUE, repr=reprSparse)

  return(C)
}
