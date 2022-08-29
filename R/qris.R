#' Estimate a quantile regression estimator of residual lifetime from survival data
#'
#' Using three estimation methods
#' 1. L1-minimization(non-smooth estimating equation)
#' 2. Induced smoothing approach (smooth estimating equation)
#' 3. Iterative procedure with induced smoothing approach (smooth estimating equation)
#'
#' @param formula  a formula expression, of the form \code{response ~ predictors}.
#'     The \code{response} is a \code{Surv} object with right censoring.
#' @param data an optional data.frame in which to interpret the variables occurring in the \code{formula}.
#' @param t0 is the followup time (or basetime of analysis). The default followup time is set to 0.
#' @param Q is the quantile. The default quantile is set to 0.5.
#' @param nB is number of multiplier bootstrapping for V matrix estimation. The default number of bootstrapping is set to 100.
#' @param method is an option for specifying the methods of parameters estimation.
#'("smooth" is default in which parameters estimates and their standard errors are obtained via induced smoothed estimating equations.
#' "nonsmooth" uses a L1-minimization method for non-smooth object functions in coefficient estimation.
#' "iterative" simultaneously estimates parameters and their standard errors based the iterative updates for parameter estimates.)
#' @param se is an option for specifying the methods of standard errors estimation
#'("pmb" is default in which parameters estimates' standard errors are obtained via partial multiplier bootstrapping. 
#' It is only for "smooth" and "iterative" options.
#' "fmb" uses a full multiplier bootstrapping in standard errors estimation. In "nonsmooth" method, "pmb" option is not available.
#' @param init is an option for specifying the initial values of the parameters estimates
#' ("rq" is default in which the estimates from the non-smooth counterpart is specified,
#' User defined vector as an initial value)
#' @param verbose Show computation status.
#' @param control controls maximum number of iteration, tolerance of convergence and whether to display output for each iteration when method = "iterative".
#' 
#' @return An object of class "\code{qris}" contains model fitting results.
#' The \code{qris} object is a list containing at least the following components:
#' \describe{
#'   \item{coefficient}{a vector of point estimates}
#'   \item{stderr}{a vector of standard error of point estiamtes}
#'   \item{vcov}{a matrix of the estimated variance-covariance matrix}
#'   \item{maxiter}{a number of itertation until convergence (only for iterative procedure)}
#'   }
#'
#' @export
#' @importFrom survival Surv survfit is.Surv
#' @importFrom quantreg rq.wfit
#' @importFrom nleqslv nleqslv
#' @importFrom stats pnorm rnorm complete.cases
#' @importFrom stringr str_replace
#' @import Rcpp
#' @example inst/examples/ex_qris.R
qris <- function(formula, data, t0 = 0, Q = 0.5, nB = 100,
                 method = c("smooth", "iterative", "nonsmooth"),
                 se = c("fmb","pmb"),
                 init = c("rq", "noeffect"),
                 verbose = FALSE, control = qris.control()) {
  scall <- match.call()
  mnames <- c("", "formula", "data")
  cnames <- names(scall)
  cnames <- cnames[match(mnames, cnames, 0)]
  mcall <- scall[cnames]
  mcall[[1]] <- as.name("model.frame")
  m <- eval(mcall, parent.frame())
  mterms <- attr(m, "terms")
  obj <- unclass(m[,1])
  method <- match.arg(method)
  se <- match.arg(se)
  if (!is.Surv(m[[1]]) || ncol(obj) > 2)
    stop("qris only supports Surv object with right censoring.", call. = FALSE)
  formula[[2]] <- NULL
  ## Create data; the first 2 columns are from Surv(), e.g., time, status, x1, x2, ...
  if (formula == ~1) {
    stop("No covariates are detected.")
  } else {
    data <- cbind(obj, model.matrix(mterms, m))
  }
  data <- as.data.frame(data)
  X <- covariate <- as.matrix(data[, -(1:2), drop = FALSE])
  nc <- ncol(covariate)
  n <- nrow(covariate)
  ## Checks
  if(nc < 2) stop("Use at least one covariate")
  if(t0 < 0) stop("basetime must be 0 and positive number")
  if(length(Q) > 1) stop("Multiple taus not allowed in qris")
  if(Q <= 0 | Q >= 1) stop("Tau must be scalar number between 0 and 1")
  ## Suppress warning message
  logZ <- suppressWarnings(log(data[,1] - t0))
  I <- as.numeric(data[,1] >= t0)
  data <- cbind(time = data[, 1], logtime = logZ, I, data[, -1])
  data[is.na(data[, 2]), 2] <- -10
  data[data[,2] == -Inf, 2] <- -10
  data[n, 4] <- 1
  colnames(data)[1:4] <- c("Z", "log(Z-t0)", "I[Z>t0]","delta")
  data <- na.omit(data)
  n <- nrow(data)
  ## Rcpp IPCW with jump weight
  # Use survfit
  sv <- survfit(Surv(data[, 1], 1 - data[, 4]) ~ 1)
  if (t0 <= sv$time[1]) {ghatt0 <- 1
  } else {ghatt0 <- sv$surv[min(which(sv$time>t0))-1]}
  W <- data[,4] / sv$surv[findInterval(data[,1], sv$time)] * ghatt0
  W[is.na(W)] <- max(W, na.rm = TRUE)
  data[, ncol(data) + 1] <- W
  
  colnames(data)[ncol(data)] <- c("weight")
  logZ <- data[,2]
  I <- data[,3]
  H <- diag(1 / n, nc, nc)
  ## Define initial value
  if (is.character(init)) {
    init <- tryCatch(match.arg(init),
                     error = function(e) {
                       tmp <- conditionMessage(e)
                       tmp <- str_replace(tmp, "'arg'", "'init'")
                       stop(tmp, " or a numerical vector \n")
                       ## on.exit(options(show.error.messages = F))
                     })
    if (init == "rq") betastart <- as.vector(rq.wfit(X, data[,2], tau = Q, weights = W)$coef)
  } else {
    if (!is.numeric(init)) stop("User specified initial value must be a numerical vector")
    if (length(init) != nc) stop("User specified initial value must match the number of covariates")
    betastart <- as.vector(init)
  }
  ## collect all useful information
  info <- list(X = X, I = I, W = W, Q = Q, nB = nB, nc = nc, n = n, H = H, t0 = t0,
               logZ = logZ, data = data, betastart = betastart, se = se,
               verbose = verbose, control = control)
  ## pass to fit
  out <- qris.fit(info, method)
  out$call <- scall
  out$varNames <- colnames(covariate)
  out$para <- list(method = method, Q = Q, t0 = t0, nB = nB)
  if (info$control$trace) {
    rownames(out$trace.coefficient) <- rownames(out$trace.stderr) <- NULL
    colnames(out$trace.coefficient) <- colnames(out$trace.stderr) <- out$varNames
  }
  out <- out[order(names(out))]
  class(out) <- "qris"
  return(out)
}

#' Calculate the weighted Kaplan-Meier estimate
#'
#' @param Time is a vector of observed time, which is minimum of failure time and censored time
#' @param censor is a vector of censoring indicator (not censored = 1, censored = 0)
#' @param wgt is a vector of weight
#'
#' @return
#' A data frame containing the following components:
#' \describe{
#'   \item{deathtime}{the observed time}
#'   \item{ndeath}{a vector of number of subject who experienced event at deathtime}
#'   \item{nrisk}{a vector of number of subject who are possible to experience event at deathtime}
#'   \item{survp}{a vector of survival probability at deathtime}
#'   }
#' @export
ghat <- function(Time, censor, wgt = 1) {
  deathtime <- c(0, sort(unique(Time[censor > 0])))
  ndeath <- colSums(outer(Time, deathtime, "==") * censor * wgt)
  nrisk <- colSums(outer(Time, deathtime, ">=") * wgt)
  survp <- cumprod(1 - ndeath / nrisk)
  data.frame(deathtime, ndeath, nrisk, survp)
}

#' Auxiliary for Controlling qris
#'
#' Auxiliary function as user interface for \code{qris} fitting.
#'
#' When \code{trace} is TRUE, output for each iteration is printed to the screen.
#' 
#' @param maxiter max number of iteration.
#' @param tol tolerance of convergence
#' @param trace a binary variable, determine whether to save output for each iteration.
#'
#' @export
#' @return A list with the arguments as components.
#' @seealso \code{\link{qris}}
qris.control <- function(maxiter = 10, tol = 1e-3, trace = FALSE) {
  list(maxiter = maxiter, tol = tol, trace = trace)
}
