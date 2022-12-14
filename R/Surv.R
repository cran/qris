#' \code{Surv} function imported from \code{survival}
#'
#' This function is imported from the \code{survival} package. See
#' \code{\link[survival]{Surv}}.
#'
#' @importFrom survival Surv
#' @name export_Surv
#' @aliases Surv
#' @return An object of class \code{Surv}.  There are methods for \code{print}, \code{is.na}, and subscripting survival objects. 
#' \code{Surv} objects are implemented as a matrix of 2 or 3 columns that has further attributes. 
#' These include the type (left censored, right censored, counting process, etc.) and labels for the states for multi-state objects.
#' Any attributes of the input arguments are also preserved in \code{inputAttributes}.  
#' This may be useful for other packages that have attached further information to data items such as labels;
#' none of the routines in the survival package make use of these values, however. 
#' In the case of \code{is.Surv}, a logical value \code{TRUE} if \code{x} inherits from class \code{"Surv"}, otherwise an \code{FALSE}.
#' @export Surv
NULL
