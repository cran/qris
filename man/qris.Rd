% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qris.R
\name{qris}
\alias{qris}
\title{Estimate a quantile regression estimator of residual lifetime from survival data}
\usage{
qris(
  formula,
  data,
  t0 = 0,
  Q = 0.5,
  nB = 100,
  method = c("smooth", "iterative", "nonsmooth"),
  se = c("fmb", "pmb"),
  init = c("rq", "noeffect"),
  verbose = FALSE,
  control = qris.control()
)
}
\arguments{
\item{formula}{a formula expression, of the form \code{response ~ predictors}.
The \code{response} is a \code{Surv} object with right censoring.}

\item{data}{an optional data.frame in which to interpret the variables occurring in the \code{formula}.}

\item{t0}{is the followup time (or basetime of analysis). The default followup time is set to 0.}

\item{Q}{is the quantile. The default quantile is set to 0.5.}

\item{nB}{is number of multiplier bootstrapping for V matrix estimation. The default number of bootstrapping is set to 100.}

\item{method}{is an option for specifying the methods of parameters estimation.
("smooth" is default in which parameters estimates and their standard errors are obtained via induced smoothed estimating equations.
"nonsmooth" uses a L1-minimization method for non-smooth object functions in coefficient estimation.
"iterative" simultaneously estimates parameters and their standard errors based the iterative updates for parameter estimates.)}

\item{se}{is an option for specifying the methods of standard errors estimation
("pmb" is default in which parameters estimates' standard errors are obtained via partial multiplier bootstrapping. 
It is only for "smooth" and "iterative" options.
"fmb" uses a full multiplier bootstrapping in standard errors estimation. In "nonsmooth" method, "pmb" option is not available.}

\item{init}{is an option for specifying the initial values of the parameters estimates
("rq" is default in which the estimates from the non-smooth counterpart is specified,
User defined vector as an initial value)}

\item{verbose}{Show computation status.}

\item{control}{controls maximum number of iteration, tolerance of convergence and whether to display output for each iteration when method = "iterative".}
}
\value{
An object of class "\code{qris}" contains model fitting results.
The \code{qris} object is a list containing at least the following components:
\describe{
  \item{coefficient}{a vector of point estimates}
  \item{stderr}{a vector of standard error of point estiamtes}
  \item{vcov}{a matrix of the estimated variance-covariance matrix}
  \item{maxiter}{a number of itertation until convergence (only for iterative procedure)}
  }
}
\description{
Using three estimation methods
1. L1-minimization(non-smooth estimating equation)
2. Induced smoothing approach (smooth estimating equation)
3. Iterative procedure with induced smoothing approach (smooth estimating equation)
}
\examples{
## ######################################### 
## Simulated data
## #########################################

data.gen <- function(n) {
  r0 <- .2 * sqrt(log(2))
  r1 <- .1 * sqrt(log(2))
  dat <- data.frame(censoring = runif(n, 0, 24.35),
                    Time0 = sqrt(-log(1 - runif(n))),
                    X = rbinom(n, 1, .5))
  dat$Time0 <- ifelse(dat$X > 0, dat$Time0 / r1, dat$Time0 / r0)
  dat$Time <- pmin(dat$Time0, dat$censoring)
  dat$status <- 1 * (dat$Time0 < dat$censoring)
  subset(dat, select = c(Time, status, X))
}

set.seed(1)
dat <- data.gen(200)
fm <- Surv(Time, status) ~ X
fit1 <- qris(fm, data = dat, t0 = 1, Q = 0.5, nB = 100, "smooth", "pmb", c(1,1))
fit2 <- qris(fm, data = dat, t0 = 1, Q = 0.5, nB = 100, "nonsmooth", "fmb", "rq")
fit3 <- qris(fm, data = dat, t0 = 1, Q = 0.5, nB = 100, "iterative", "fmb", "rq",
             control = qris.control(maxit = 20, tol = 1e-3, trace = TRUE))

summary(fit1)
summary(fit2)
summary(fit3)

## #########################################
## Real data application
## #########################################
data(cancer, package = "survival")
lung2 <- subset(lung, select = c(time, status, age, sex))
## tidy up the data
lung2$status <- lung2$status - 1
lung2$sex <- lung2$sex - 1

fm <- Surv(time, status) ~ age + sex
fit1 <- qris(fm, data = lung2, t0 = 0, Q = 0.5, nB = 100, "iterative", "pmb", "rq")
fit2 <- qris(fm, data = lung2, t0 = 30, Q = 0.5, nB = 100, "nonsmooth", "fmb", c(1, 0, 1))
fit3 <- qris(fm, data = lung2, t0 = 100, Q = 0.5, nB = 100,"smooth", "pmb", "rq")

summary(fit1)
summary(fit2)
summary(fit3)

plot(fit2, Qs = 4:6 / 10)


}
