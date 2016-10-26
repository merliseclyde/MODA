oda.lm = function(formula, data,  subset, weights,
                  na.action="na.omit",
                  n.models=NULL,  prior=GDP(),
                  modelprior=beta.binomial(1,1),
                  method="coordinate-asscent", update=NULL,
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
                  "BAS" = .Call("sampleworep",
                                Yvec, X, sqrt(weights),
                                prob, modeldim,
                                incint=as.integer(int),
                                alpha= as.numeric(alpha),
                                method=as.integer(method.num), modelprior=modelprior,
                                update=as.integer(update),
                                Rbestmodel=as.integer(bestmodel),
                                plocal=as.numeric(prob.local),
                                PACKAGE="BAS"),
                  "MCMC+BAS"= .Call("mcmcbas",
                                    Yvec, X, sqrt(weights),
                                    prob, modeldim,
                                    incint=as.integer(int),
                                    alpha= as.numeric(alpha),
                                    method=as.integer(method.num), modelprior=modelprior,
                                    update=as.integer(update),
                                    Rbestmodel=as.integer(bestmodel),
                                    plocal=as.numeric(1.0 - prob.rw), as.integer(Burnin.iterations),
                                    as.integer(MCMC.iterations), as.numeric(lambda),as.numeric(delta),
                                    PACKAGE="BAS"),
                  "MCMC"= .Call("mcmc_new",
                                Yvec, X, sqrt(weights),
                                prob, modeldim,
                                incint=as.integer(int),
                                alpha= as.numeric(alpha),
                                method=as.integer(method.num), modelprior=modelprior,
                                update=as.integer(update),
                                Rbestmodel=as.integer(bestmodel),
                                plocal=as.numeric(1.0 - prob.rw), as.integer(Burnin.iterations),
                                as.integer(MCMC.iterations), as.numeric(lambda),as.numeric(delta),
                                as.integer(thin),
                                PACKAGE="BAS"),
                  "MCMC_old"= .Call("mcmc",
                                    Yvec, X, sqrt(weights),
                                    prob, modeldim,
                                    incint=as.integer(int),
                                    alpha= as.numeric(alpha),
                                    method=as.integer(method.num), modelprior=modelprior,
                                    update=as.integer(update),
                                    Rbestmodel=as.integer(bestmodel),
                                    plocal=as.numeric(1.0 - prob.rw), as.integer(Burnin.iterations),
                                    as.integer(MCMC.iterations), as.numeric(lambda),as.numeric(delta),
                                    as.integer(thin),
                                    PACKAGE="BAS"),
                  "AMCMC" = .Call("amcmc",
                                  Yvec, X, sqrt(weights),
                                  prob, modeldim,
                                  incint=as.integer(int),
                                  alpha= as.numeric(alpha),
                                  method=as.integer(method.num), modelprior=modelprior,
                                  update=as.integer(update),
                                  Rbestmodel=as.integer(bestmodel),
                                  plocal=as.numeric(1.0-prob.rw), as.integer(Burnin.iterations),
                                  as.integer(MCMC.iterations), as.numeric(lambda),as.numeric(delta),
                                  PACKAGE="BAS"),
                  #    "MAXeffect" = .Call("posisearch",
                  #     Yvec, X,
                  #     prob, modeldim,
                  #     incint=as.integer(int),
                  #     alpha= as.numeric(alpha),
                  #     method=as.integer(method.num), modelprior=modelprior,
                  #     update=as.integer(update),
                  #     Rbestmodel=as.integer(bestmodel),
                  #     Rbestmarg=as.numeric(bestmarg),
                  #     plocal=as.numeric(1.0-prob.rw), as.integer(Burnin.iterations),
                  #     as.integer(MCMC.iterations), as.numeric(lambda),as.numeric(delta),
                  #     PACKAGE="BAS"),
                  "deterministic" = .Call("deterministic",
                                          Yvec, X, sqrt(weights),
                                          prob, modeldim,
                                          incint=as.integer(int),
                                          alpha= as.numeric(alpha),
                                          method=as.integer(method.num),modelprior=modelprior,
                                          PACKAGE="BAS")
  )

  result$namesx=namesx
  result$n=length(Yvec)
  result$prior=prior
  result$modelprior=modelprior
  result$call=call

  class(result) = c("oda","bma")
    return(result)
}
