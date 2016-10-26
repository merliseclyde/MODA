
#' Title
#'
#' @param formula an object of class \code{\link{formula}} (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which lm is called.
#' @param subset an optional vector specifying a subset of observations to be used in the fitting process.
#' @param weights  an optional vector of weights to be used in the fitting process. Should be NULL or a numeric vector. If non-NULL, weighted least squares is used with weights weights (that is, minimizing sum(w*e^2)); otherwise ordinary least squares is used.
#' @param na.action a function which indicates what should happen when the data contain NAs. The default is set by the na.action setting of options, and is na.fail if that is unset. The ‘factory-fresh’ default is na.omit. Another possible value is NULL, no action. Value na.exclude can be useful.
#' @param n.models  number of models to sample
#' @param betaprior  prior on coefficients
#' @param modelprior  prior on the model space
#' @param method  the method to be used
#' @param force.in variables that shoujld always be included
#' @param iterations  number of iterations
#'
#' @return None yet
#' @importFrom stats contrasts model.matrix model.response model.weights weighted.mean
#' @importFrom BAS beta.binomial uniform
#' @export
#' @references Clyde, MA, & Ghosh, J. (2012). Finite population estimators in stochastic search variable selection. Biometrika, 99 (4), 981-988. <10.1093/biomet/ass040>
#' @details  add more here
#' @examples
#' library(MASS)
#' data(UScrime)
#' oda.lm(y ~ ., data=UScrime)
#'
oda.lm = function(formula, data,  subset, weights,
                  na.action="na.omit",
                  n.models=NULL,  betaprior=GDP(),
                  modelprior=beta.binomial(1,1),
                  method="coordinate-ascent",
                  force.in,
                  iterations=NULL
                  )  {



  call = match.call()

  # from lm
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action",
               "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())

  #data = model.frame(formula, data, na.action=na.action, weights=weights)
  n.NA = length(attr(mf, 'na.action'))

  if (n.NA > 0) {
    warning(paste("dropping ", as.character(n.NA),
                  "rows due to missing data"))
  }

  Y = model.response(mf, "numeric")
  mt <- attr(mf, "terms")
  X = model.matrix(mt, mf, contrasts)
  #X = model.matrix(formula, mf)

  namesx = dimnames(X)[[2]]
  namesx[1] = "Intercept"
  n <- dim(X)[1]

  weights=as.vector(model.weights(mf))
  if (is.null(weights)) weights = rep(1, n)

  if (length(weights) != n) stop(simpleError(paste("weights are of length ", length(weights), "not of length ", n)))

  mean.x = apply(X[,-1, drop=F], 2, weighted.mean, w=weights)
  ones = X[,1]
  X = cbind(ones, sweep(X[, -1, drop=FALSE], 2, mean.x))
  p <-  dim(X)[2]  # with intercept



 # Modify below to call scala
    result = switch(method,
                  "coordinate-ascent" = rscala(),
                  "oda"= rscala(),
                  "enumerate"= rscala(),
                  "adapt"= rscala()
  )

  result$namesx=namesx
  result$n=n
  result$betaprior=betaprior
  result$modelprior=modelprior
  result$call=call

  class(result) = c("oda","bma")
    return(result)
}
