#' Sequential Tests for the GEVr Model
#'
#' Sequentially performs the entropy difference (ED) test for the GEVr model.
#' @param data Data should be contain n rows, each a GEVr observation.
#' @param method string indicating the method. Only \code{ed} is supported.
#' @param ... currently ignored (unused arguments) for compatibility with former versions of the package
#' @examples
#' x <- rgevr(200, 5, loc = 0.5, scale = 1, shape = 0.25)
#' gevrSeqTests(x, method = "ed")
#' @return
#' Function returns a data frame containing the test statistics, estimates, and p-value results of the sequential tests.
#'
#' \item{r}{Value of r to be tested.}
#' \item{p.values}{Raw p-values from the individual tests at each value of r.}
#' \item{ForwardStop}{Transformed p-values according to the ForwardStop stopping rule.}
#' \item{StrongStop}{Transformed p-values according to the StrongStop stopping rule.}
#' \item{statistic}{Returned test statistics of each individual test.}
#' \item{est.loc}{Estimated location parameter for the given r.}
#' \item{est.scale}{Estimated scale parameter for the given r.}
#' \item{est.shape}{Estimated shape parameter for the given r.}
#' @details GEVr data (in matrix x) should be of the form \eqn{x[i,1] > x[i, 2] > \cdots > x[i, r]} for each observation \eqn{i = 1, \ldots, n}.
#' See function `pSeqStop' for details on transformed p-values.
#' @import parallel
#' @export
gevrSeqTests <- function(data, method = "ed", ...) {
  if(method !="ed"){stop("Invalid test statistic: method must be `ed`")}
  data <- as.matrix(data)
  R <- ncol(data)
  method <- match.arg(method)
  if(R == 1)
      stop("R must be at least two")
    result <- matrix(0, R-1, 8)
    for(i in 2:R) {
      result[i-1, 1] <- i
      fit <- gevrEd(data[, 1:i])
      result[i-1, 2] <- fit$p.value
      result[i-1, 5] <- fit$statistic
      result[i-1, 6:8] <- fit$theta
    }
  result[, 3] <- rev(pSeqStop(rev(result[, 2]))$ForwardStop)
  result[, 4] <- rev(pSeqStop(rev(result[, 2]))$StrongStop)
  colnames(result) <- c("r", "p.values", "ForwardStop", "StrongStop", "statistic", "est.loc", "est.scale", "est.shape")
  as.data.frame(result)
}
