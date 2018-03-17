## Contains functions useful for growth curves.

##' Logistic growth
##'
##' Simulate phenotypes via the following model (\href{http://www.genetics.org/content/161/4/1751.abstract}{Ma et al (2002)}): y(t) = g(t) + epsilon(t) where g(t) = a / (1 + b exp(-r t))) and epsilon(t) ~ N(0, sigma^2).
##' @param t vector of time points
##' @param a asymptotic value of g when t tends to infinity
##' @param g.t0 initial value of g at time t = 0
##' @param r relative rate of growth
##' @param sigma2 variance of the errors
##' @return list
##' @author Timothee Flutre
##' @examples
##' \dontrun{set.seed(1859)
##' model <- simulLogistic(t=1:20, a=50, g.t0=1, r=0.6, sigma2=1)
##' plot(x=model$t, y=model$g.t, type="b", las=1, xlab="time (t)", ylab="g(t)",
##'      main="Logistic without noise")
##' plot(x=model$t, y=model$y.t, type="b", las=1, xlab="time (t)", ylab="y(t)",
##'      main="Logistic with noise")
##' }
##' @export
simulLogistic <- function(t=1:20, a=50, g.t0=1, r=0.6, sigma2=0){
  stopifnot(is.numeric(t),
            length(t) > 0,
            is.numeric(a),
            length(a) == 1,
            is.numeric(g.t0),
            length(g.t0) == 1,
            is.numeric(r),
            length(r) == 1,
            is.numeric(sigma2),
            length(sigma2) == 1,
            sigma2 >= 0)

  t <- t[! is.na(t)]

  b <- (a - g.t0) / g.t0

  g.t <- a / (1 + b * exp(- r * t))
  y.t <- g.t
  if(sigma2 > 0)
    y.t <- y.t + stats::rnorm(n=length(t), mean=0, sd=sqrt(sigma2))

  tI <- log(b) / r
  g.tI <- a / 2

  return(list(t=t, a=a, g.t0=g.t0, r=r, b=b,
              g.t=g.t, y.t=y.t, tI=tI, g.tI=g.tI))
}

##' Generalised logistic growth (Richards' curve)
##'
##' Simulate phenotypes via the following model (\href{https://doi.org/10.1093/jxb/10.2.290}{Richards (1959)}): y(t) = g(t) + epsilon(t) where g(t) = a + (k - a) / ((c + b exp(-r t))^(1/nu)) and epsilon(t) ~ N(0, sigma^2).
##' @param t vector of time points
##' @param a lower asymptote
##' @param k upper asymptote (if a=0, k is called the carrying capacity)
##' @param r growth rate
##' @param nu affects near which asymptote maximum growth occurs (> 0)
##' @param b related to the value at t=0
##' @param c typically takes a value of 1
##' @param sigma2 variance of the errors
##' @return list
##' @author Timothee Flutre
##' @examples
##' \dontrun{set.seed(1859)
##' model <- simulGeneralisedLogistic(t=seq(-1.5, 3.5, 0.1), a=0, k=1, r=3, nu=0.5, b=0.5, c=1, sigma2=0.01)
##' plot(x=model$t, y=model$g.t, type="b", las=1, xlab="time (t)", ylab="g(t)",
##'      main="Generalised logistic without noise")
##' plot(x=model$t, y=model$y.t, type="b", las=1, xlab="time (t)", ylab="y(t)",
##'      main="Generalised logistic with noise")
##' }
##' @export
simulGeneralisedLogistic <- function(t=seq(-1.5, 3.5, 0.5),
                                     a=0, k=1, r=3, nu=0.5, b=0.5, c=1,
                                     sigma2=0){
  stopifnot(is.numeric(t),
            length(t) > 0,
            is.numeric(a),
            length(a) == 1,
            is.numeric(k),
            length(k) == 1,
            is.numeric(r),
            length(r) == 1,
            is.numeric(nu),
            length(nu) == 1,
            nu > 0,
            is.numeric(b),
            length(b) == 1,
            is.numeric(c),
            length(c) == 1,
            is.numeric(sigma2),
            length(sigma2) == 1,
            sigma2 >= 0)

  t <- t[! is.na(t)]

  g.t <- a + (k - a) / ((c + b * exp(- r * t))^(1 / nu))
  y.t <- g.t
  if(sigma2 > 0)
    y.t <- y.t + stats::rnorm(n=length(t), mean=0, sd=sqrt(sigma2))

  return(list(t=t, a=a, k=k, r=r, nu=nu, b=b, c=c,
              g.t=g.t, y.t=y.t))
}
