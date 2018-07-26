## Contains functions useful for spatial statistics.

##' AR1xAR1 simulation
##'
##' Simulate random samples from a separable AR1xAR1 process as in \href{https://dx.doi.org/10.1016/0378-3758(95)00066-6}{Martin (1996)}.
##' @param n number of samples
##' @param R number of rows
##' @param C number of columns
##' @param rho.r correlation between rows
##' @param rho.c correlation between columns
##' @param sigma2 variance of the errors
##' @return array
##' @author Timothee Flutre
##' @export
simulAr1Ar1 <- function(n=1, R, C, sigma2, rho.r=0, rho.c=0){
  stopifnot(R > 1,
            C > 1,
            abs(rho.r) <= 1,
            abs(rho.c) <= 1,
            sigma2 > 0)

  epsilons <- array(data=stats::rnorm(n=R*C*n, mean=0, sd=sqrt(sigma2)),
                    dim=c(R, C, n))

  tmp <- lapply(1:n, function(i){
    Xi <- matrix(data=NA, nrow=R, ncol=C)
    Xi[1,1] <- sqrt(1/sigma2) * epsilons[1,1,i] # first cell
    for(c in 2:C) # first row
      Xi[1,c] <- rho.c * Xi[1,c-1] + sqrt((1 - rho.c^2) * 1 / sigma2) *
        epsilons[1,c,i]
    for(r in 2:R) # first column
      Xi[r,1] <- rho.r * Xi[r-1,1] + sqrt((1 - rho.r^2) * 1 / sigma2) *
        epsilons[r,1,i]
    for(r in 2:R)
      for(c in 2:C)
        Xi[r,c] <- rho.r * Xi[r-1,c] + rho.c * Xi[r,c-1] -
          rho.r * rho.c * Xi[r-1,c-1] + epsilons[r,c,i]
    return(Xi)
  })

  x <- array(data=do.call(c, tmp),
             dim=c(R, C, n))

  return(x)
}
