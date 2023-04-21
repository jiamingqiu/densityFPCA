## methods with given density family ==========================================

# some tools
set.knowFamEsti_control <- function(envir){
  # setting control list for knowFamEsti
  # args:
  # envir = the current within knowFamEsti.
  # returns:
  # a set control list.

  control <- envir$control
  control$return.scale <- match.arg(
    control$return.scale,
    c("parameter", "log", "origin", "optim")
  )
  if (control$return.scale %in% c("log", "origin")) {
    if (is.null(control$grid)) {
      stop('Need grid computing density.')
    }
    if (any(!is.finite(control$grid))) {
      stop('Non finite grid value.')
    }
  }
  return(control)
}

knowFamEsti <- function(mat.obsv, den.fam, esti.method, control = list()){
  # compute estimation with known family
  # args:
  # mat.obsv: a matrix or list where one row/element is one sample;
  # den.fam: a density family;
  # esti.method: a character array of the method name, say
  #              "MLE", "emTruncMixGauss";
  # control: a list of controling arguments:
  #   return.scale = what to return: "log", "origin", "parameter"(default), "optim";
  #   grid = grid to compute return when return.scale = "log"/"origin";
  #   init.par = initial value for optim/em,
  #             miss -> MLE use midpoint of den.fam$range, EM will use kmeans.
  #  calMLE only options:
  #   method = method to use for optimization, default = DIRECT-P;
  #   nl.control = list of options for nloptr, c.f. nloptr::nl.opts, default = list();
  #  emTruncMixGauss only options:, if missing, use emTruncMixGauss default.
  #   tol = tolenrence to stop iteration, default is 1e-10;
  #   maxit = maximum number of iteration, default 1000;
  #   min.var = minimum variance of component, if during iteration there
  #     would be a mixture reaching a smaller variance, its variance will
  #     be inflated to the mean of variance that are >= min.var, default 1e-3.
  #
  # returns:
  # a list of
  # - a matrix where estimated parameters / (log-)densities as rows,
  # - grid, NULL if not necessary;
  # - a data.frame for indexing obervations and etimating methods;
  # - list of error message, if any;
  # - time.

  # START LEGACY
  # grid: grid on which estimated density to be evaluated.
  # esti.option: list of additional argument for esti.method, currently
  #              mle.init = initPar for calMLE, if "auto", use means of denFam$gen_par;
  #              one can also specify directly by a matrix of nrow(new.obsv) X dim.par
  #              mle.multiStart = multiStart for calMLE, "ON" or "OFF";
  # TODO:BETTER? k = k, num of mixtures for MixGauss related methods;
  #              scale = 'origin'(default) or 'log'.
  # err: what to do if esti.method fails, default = fill in NA.
  #
  # returns:
  # a list of
  #   a matrix where estimated densities as rows;
  #   a data.frame for indexing the rows by idx.obsv, esti.method;
  #   a data.frame for timing;
  #   grid.
  # END LEGACY

  # input arguments
  control <- set.knowFamEsti_control(environment())
  esti.method <- match.arg(
    esti.method,
    c('MLE', 'emTruncMixGauss'),
    several.ok = TRUE
  )
  # if asked EM
  if (is.element("emTruncMixGauss", esti.method)) {
    num.mixture <- (ncol(den.fam$par.range) + 1) / 3
  }
  # # if asked KDE
  # if ("KDE" %in% esti.method & control$return.scale == 'parameter') {
  #   esti.method <- esti.method[esti.method != "KDE"]
  #   warning("KDE only supported when returning (log-)density, dropped here.")
  # }
  # if mat.obsv is not a matrix
  if (is.matrix(mat.obsv)) {
    mat.obsv <- mat2list(mat.obsv, by = 1)
  }

  # error handling function
  err.f <- function(err) err

  # handling init.par of calMLE, make it a matrix
  if (is.null(control$init.par)) {
    mle.init.par <- t(matrix(
      colMeans(den.fam$par.range),
      nrow = ncol(den.fam$par.range),
      ncol = length(mat.obsv)
    ))
  } else {
    if (!is.matrix(control$init.par)) {
      mle.init.par <- t(matrix(
        control$init.par,
        nrow = ncol(den.fam$par.range),
        ncol = length(mat.obsv)
      ))
    } else {
      if (ncol(control$init.par) != ncol(den.fam$par.range) |
          nrow(control$init.par) != length(mat.obsv)) {
        stop('Check dimension of init.par.')
      }
      mle.init.par <- control$init.par
    }
  }
  # list of functions used for estimation
  ls.f <- list()
  w.control <- control
  w.control$return.scale <- 'optim'
  for(i in seq_along(esti.method)){
    # get estimating function
    ls.f[[i]] <- switch (
      esti.method[i],
      # "KDE" = function(obsv){
      #   tryCatch(calKDE(obsv, control = control), error = err.f)
      # },
      "MLE" = function(obsv, init.par){
        w.control$init.par <- init.par
        tryCatch(calMLE(obsv, den.fam, control = w.control), error = err.f)
      },
      "emTruncMixGauss" = function(obsv, init.par){
        em.control <- w.control
        if(!is.null(control$init.par)){
          em.control$init.par <- tMixGauss.repar(work.par = init.par)
        }
        tryCatch(
          emTruncMixGauss(
            obsv, edge = den.fam$range,
            num.mixture = num.mixture, control = em.control),
          error = err.f
        )
      }
    )
  }

  # compute estimation
  mat.par <- matrix(
    0, ncol = ncol(den.fam$par.range),
    nrow = length(mat.obsv) * length(esti.method)
  )
  # list of timing.
  ls.t <- list()
  # list of error message.
  idx.message <- 1
  ls.message <- list()
  # get optim results first.
  idx.run <- 1
  ls.opt.res <- list()

  for (idx.esti in seq_along(ls.f)) {
    ls.t[[idx.esti]] <- system.time({
      for (i in seq_along(mat.obsv)) {
        tm.res <- (ls.f[[idx.esti]])(mat.obsv[[i]], mle.init.par[i, ])

        if (!inherits(tm.res, 'error')) {
          ls.opt.res[[idx.run]] <- tm.res
          mat.par[idx.run, ] <- tm.res$par
        } else {
          ls.opt.res[[idx.run]] <- NULL
          mat.par[idx.run, ] <- rep(NA, ncol(den.fam$par.range))
          ls.message[[idx.message]] <- list(
            idx.obsv = i,
            esti.method = esti.method[idx.esti],
            error = tm.res
          )
          idx.message <- idx.message + 1
        }
        idx.run <- idx.run + 1
      }
    })
  }

  # now mat.par in format of 1 row = 1 par, and in trunc of esti.method
  # indexing data.frame
  df.idx <- data.frame(
    idx.obsv = rep(seq_along(mat.obsv), times = length(esti.method)),
    esti.method = rep(esti.method, each = length(mat.obsv))
  )

  # get the time
  df.t <- do.call(rbind, ls.t)
  df.t <- rbind(df.t, colSums(df.t))
  df.t <- cbind(data.frame(step = c(esti.method, 'TOTAL')), df.t)

  # returning
  if (control$return.scale == 'optim') {
    return(list(
      res = ls.opt.res,
      idx = df.idx,
      grid = NULL,
      error = ls.message,
      time = df.t
    ))
  }
  if (control$return.scale == 'parameter') {
    return(list(
      res = mat.par,
      idx = df.idx,
      grid = NULL,
      error = ls.message,
      time = df.t
    ))
  }
  # the following need par2pdf
  mat.pdf <- par2pdf(den.fam, mat.par, grid = control$grid)
  if (control$return.scale == 'origin') {
    return(list(
      res = mat.pdf,
      idx = df.idx,
      grid = control$grid,
      error = ls.message,
      time = df.t
    ))
  }
  if (control$return.scale == 'log') {
    return(list(
      res = log(mat.pdf),
      idx = df.idx,
      grid = control$grid,
      error = ls.message,
      time = df.t
    ))
  }
}

# individual methods ==========================================================

# calculating MLE with known family
calMLE = function(newx, den.fam, control = list()){
  # calculating MLE with known density family.
  # args:
  # newx = array of new observations;
  # den.fam = density family;
  # control = a list of controling arguments:
  #   init.par = initial value for optim, missing -> midpoint of den.fam$range;
  #   method = method to use for optimization, default = DIRECT-P;
  #   nl.control = list of options for nloptr, c.f. nloptr::nl.opts, default = list();
  #   return.scale = what to return: "log", "origin", "parameter"(default), "optim";
  #   grid = grid to compute return when return.scale = "log"/"origin".
  # returns:
  # either the MLE or the corresponding log-density / density on grid,
  # or the optimization result, for debugging.

  # filling default and check input
  if (is.null(control$init.par)) {
    control$init.par <- as.numeric(colSums(den.fam$par.range) / 2)
  }
  if (is.null(control$method)) {
    control$method <- "DIRECT-P"
  }
  if (is.null(control$nl.control)) {
    control$nl.control <- list()
  }
  control$return.scale <- match.arg(
    control$return.scale,
    c("parameter", "log", "origin", "optim")
  )
  if (control$return.scale %in% c("log", "origin")) {
    if (is.null(control$grid)) {
      stop("Need grid to compute return.")
    }
  }

  # obtaining objective function and its gradient(if possible.)
  obj.f <- getNegLogll(newx, den.fam)
  if (!is.null(den.fam$logpdf.gr)) {
    obj.f.gr <- getNegLogll.gr(newx, den.fam)
  } else{
    obj.f.gr <- NULL
  }

  opt.res <- safeOptim(
    x0 = control$init.par, fn = obj.f, gr = obj.f.gr,
    lower = den.fam$par.range[1, ],
    upper = den.fam$par.range[2, ],
    method = control$method, backup.method = "DIRECT",
    nl.control = control$nl.control
  )

  if (control$return.scale == 'optim') {
    return(opt.res)
  }
  if (control$return.scale == 'log') {
    return(log(den.fam$pdf(grid, opt.res$par)))
  }
  if (control$return.scale == 'origin') {
    return(den.fam$pdf(grid, opt.res$par))
  }
  if (control$return.scale == 'parameter') {
    return(opt.res$par)
  }
}
# LEGACY
# calMLE = function(newx, denFam, initPar, grid, scale = 'log', multiStart = FALSE){
#   input.multiStart = multiStart
#   # just to translate into acceptable argument.
#   if(typeof(multiStart[1]) == 'logical'){
#     multiStart = c(1,0)
#   }
#   dimPar = length(initPar)
#
#   if(dimPar == 1){
#     opt_res = optim(initPar, function(par){-1*sum(log(denFam$pdf(newx, par)))}, method = 'Brent', lower = -1e+3, upper = 1e+3)
#   }else{
#     opt_res = smartOptim(
#       initPar, function(par){-1*sum(log(denFam$pdf(newx, par)))},
#       method = 'L-BFGS-B', multiStart = multiStart,
#       lower = denFam$par.range[1, ],
#       upper = denFam$par.range[2, ]
#     )
#   }
#   if(opt_res$convergence!=0){
#     warning(c('optim convergence code ', opt_res$convergence), ' while calculating MLE.')
#   }
#   if(is.na(opt_res$value)){
#     stop('NA as maximum.')
#   }
#   if (scale == 'log') {
#     return(log(denFam$pdf(grid, opt_res$par)))
#   }
#   if (scale == 'origin') {
#     return(denFam$pdf(grid, opt_res$par))
#   }
#   if (scale == 'parameter') {
#     return(c(opt_res$value, opt_res$par))
#   }
# }


# EM algorithm for Mixture Gaussian family (no truncation)
# Caution: no variance inflating step.
emMixGauss <- function(obsv, num.mixture, control = list()){
  # EM algorithm for Gaussian mixture model, at least k+1 observations.
  # args:
  # obsv = observations;
  # num.mixture = number of mixtures (k);
  # control = a list of controling arguments:
  #   return.scale = what to return: "log", "origin", "parameter"(default);
  #   grid = grid to compute return when return.scale = "log"/"origin";
  #   tol = tolenrence to stop iteration, default is 1e-10;
  #   maxit = maximum number of iteration, default 1000.
  # returns:
  # either the list of parameter estimates and number of iteration;
  # or the corresponding log-density / density on grid,
  # Detail:
  # If iteration determine there would be no observation in a mixture,
  # the number of mixture will be shrinked by 1, and algorithm restart.

  # handling control list
  control$return.scale <- match.arg(
    control$return.scale,
    c("parameter", "log", "origin")
  )
  if (control$return.scale %in% c("log", "origin")) {
    if (is.null(control$grid)) {
      stop('Need grid computing density.')
    }
    if (any(!is.finite(control$grid))) {
      stop('Non finite grid value.')
    }
  }
  if (is.null(control$tol)) {
    tol <- 1e-10
  } else {
    tol <- control$tol
  }
  if (is.null(control$maxit)) {
    maxit <- 1000
  } else {
    maxit <- control$maxit
  }

  # renaming
  x <- obsv
  k <- num.mixture
  n <- length(x)

  flg.done = FALSE
  while(!flg.done){
    # initial by k-means
    res.km = kmeans(x, centers = k)
    mixWeight = as.numeric(table(res.km$cluster) / n)
    mixMean = as.numeric(res.km$centers)
    mixVar = sapply(1:k, function(idx){var(x[res.km$cluster == idx])})

    matCoef = matrix(0, nrow = length(x), ncol = k)
    for(j in 1:k){
      matCoef[,j] = mixWeight[j] * dnorm(x, mean = mixMean[j], sd = sqrt(mixVar[j]))
    }
    matCoef = matCoef / rowSums(matCoef)

    count = 1
    flg.nextK = FALSE
    while(count <= maxit){
      newWeight = colSums(matCoef) / sum(matCoef)
      newMean = colSums(matCoef * x) / colSums(matCoef)
      # newMean = sort(newMean)
      newVar = colSums(matCoef * (x - rep(newMean, each = n))^2) / colSums(matCoef)
      if(anyNA(newVar)){
        flg.nextK = TRUE
        break
      }
      if(sum((newWeight - mixWeight)^2 + (newMean - mixMean)^2 + (newVar - mixVar)^2) < tol){
        break
      }
      mixWeight = newWeight
      mixMean = newMean
      mixVar = newVar
      for(j in 1:k){
        matCoef[,j] = mixWeight[j] * dnorm(x, mean = mixMean[j], sd = sqrt(mixVar[j]))
      }
      matCoef = matCoef / rowSums(matCoef)
      count = count + 1
    }
    if(flg.nextK){
      k = k - 1
    }else{
      flg.done = TRUE
    }
  }

  # returning
  if (control$return.scale == 'parameter') {
    return(list(
      weight = mixWeight[order(mixMean)],
      mean = mixMean[order(mixMean)],
      var = mixVar[order(mixMean)],
      num.mixture = k,
      iter = count
    ))
  }
  grid  <- control$grid
  tm.mat <- matrix(0, ncol = k, nrow = length(grid))
  for(j in 1:k){
    tm.mat[,j] <- mixWeight[j] * dnorm(grid, mean = mixMean[j], sd = sqrt(mixVar[j]))
  }
  if (control$return.scale == 'origin') {
    return(rowSums(tm.mat))
  }
  if (control$return.scale == 'log') {
    return(log(rowSums(tm.mat)))
  }
  # return(rowSums(tm.mat))
}

# a function calculating first and second moments of truncated normal r.v.
# higher order by pracma::fderiv
mTruncNrml = function(a, b, p = 1, mu = 0, sd = 1){

  # smpl = rnorm(1000, mu, sd)
  # smpl = smpl[smpl>=a & smpl<=b]
  # hist(smpl)
  # mean(smpl^p)

  w.a = (pmin(a, b)-mu)/sd
  w.b = (pmax(a, b)-mu)/sd

  # temp: 1st moment of N(0, 1) trunc. [w.a, w.b]
  tmM = -1 * (dnorm(w.b) - dnorm(w.a)) / (pnorm(w.b) - pnorm(w.a))
  if(p == 1){
    return(
      mu + sd * tmM
    )
  }
  if(p == 2){
    return(
      mu^2 + 2 * sd * mu * tmM + sd^2 - sd^2 * (w.b * dnorm(w.b) - w.a * dnorm(w.a)) / (pnorm(w.b) - pnorm(w.a))
    )
  }

  # the following only work if a & b are scales
  w.a = pmin(a,b) - mu
  w.b = pmax(a,b) - mu
  mgf = function(t){
    exp(mu*t + sd^2*t^2/2) *
      (pnorm(w.b/sd - sd*t) - pnorm(w.a/sd - sd*t)) /
      (pnorm(w.b/sd) - pnorm(w.a/sd))
  }

  return(
    pracma::fderiv(mgf, 0, n=p)
  )
}


emTruncMixGauss <- function(obsv, num.mixture, edge, control = list()){
  # EM algorithm for truncated mixture Gaussian. Truncation after mixing.
  # args:
  # obsv = observations;
  # num.mixture = number of mixtures (k);
  # edge = truncation boundary;
  # control = a list of controling arguments:
  #   init.par = a list of inital weight, mean and var, missing -> kmeans,
  #              the means should not be identical at the very least;
  #   return.scale = what to return: "log", "origin", "parameter"(default), "optim";
  #   grid = grid to compute return when return.scale = "log"/"origin";
  #   tol = tolenrence to stop iteration, default is 1e-10;
  #   maxit = maximum number of iteration, default 1000;
  #   min.var = minimum variance of component, if during iteration there
  #     would be a mixture reaching a smaller variance, its variance will
  #     be inflated to the mean of variance that are >= min.var, default 1e-3.
  # returns:
  # either the list of parameter estimates and number of iteration;
  # or the corresponding log-density / density on grid,
  # Detail:
  # If iteration determine there would be no observation in a mixture,
  # the number of mixture will be shrinked by 1, and algorithm restart,
  # if init.par provided, the mixture component with smallest mixture weight
  # will be dropped.
  # If return.scale = "optim", return a list similar to optim result, where
  # $par is parameters in the format of tMixGauss(), intended for internal use.
  # Legacy arguments:
  # (x, k, edge, returnDensity = FALSE, tol = 1e-10, maxit = 1000, minVar = 1e-3)

  # handling control list
  if (!is.null(control$init.par)) {
    ls.init.par <- control$init.par
    if (anyNA(ls.init.par, recursive = TRUE))
      stop('invalid init.par: contains NA')
    if (!all(is.element(c('mean', 'weight', 'var'), names(ls.init.par))))
      stop('invalid init.par: missing mean/weight/var')
    if (!all.equal(sum(ls.init.par$weight), 1))
      stop('invalid init.par: sum(weight) != 1')
    if (sd(ls.init.par$mean) < .Machine$double.eps ^ (1/2))
      warning('init.par: mixture mean too close, EM result may not be valid')
    # sort mixture components by decreasing weights
    ord.idx <- order(ls.init.par$weight, decreasing = TRUE)
    ls.init.par$weight <- ls.init.par$weight[ord.idx]
    ls.init.par$mean <- ls.init.par$mean[ord.idx]
    ls.init.par$var <- ls.init.par$var[ord.idx]
    # rm(ord.idx)
  }
  control$return.scale <- match.arg(
    control$return.scale,
    c("parameter", "log", "origin", "optim")
  )
  if (control$return.scale %in% c("log", "origin")) {
    if (is.null(control$grid)) {
      stop('Need grid computing density.')
    }
    if (any(!is.finite(control$grid))) {
      stop('Non finite grid value.')
    }
  }
  if (is.null(control$tol)) {
    tol <- 1e-10
  } else {
    tol <- control$tol
  }
  if (is.null(control$maxit)) {
    maxit <- 1000
  } else {
    maxit <- control$maxit
  }
  if (is.null(control$min.var)) {
    minVar <- 1e-3
  } else {
    minVar <- control$min.var
  }

  # handling truncation edge
  if(length(edge) == 1){
    w.a = -abs(edge)
    w.b = abs(edge)
  }else{
    w.a = min(edge)
    w.b = max(edge)
  }
  # renaming
  x <- obsv
  k <- num.mixture
  n <- length(x)

  flg.done = FALSE
  while(!flg.done){
    if(is.null(control$init.par)){
      # initial by k-means
      res.km = kmeans(x, centers = k)
      mixWeight = as.numeric(table(res.km$cluster) / n)
      mixMean = as.numeric(res.km$centers)
      mixVar = sapply(1:k, function(idx){var(x[res.km$cluster == idx])})
      # sort the components for comparison
      ord.idx <- order(mixMean)
      mixWeight <- mixWeight[ord.idx]
      mixMean <- mixMean[ord.idx]
      mixVar <- mixVar[ord.idx]
    }else{
      mixWeight <- ls.init.par$weight[seq_len(k)]
      mixWeight <- mixWeight / sum(mixWeight)
      mixMean <- ls.init.par$mean[seq_len(k)]
      mixVar <- ls.init.par$var[seq_len(k)]
      # sort the components for comparison
      ord.idx <- order(mixMean)
      mixWeight <- mixWeight[ord.idx]
      mixMean <- mixMean[ord.idx]
      mixVar <- mixVar[ord.idx]
    }

    matCoef = matrix(0, nrow = length(x), ncol = k)
    for(j in 1:k){
      matCoef[,j] =
        mixWeight[j] * dnorm(x, mean = mixMean[j], sd = sqrt(mixVar[j]))
    }
    matCoef = matCoef / rowSums(matCoef)

    count = 1
    flg.nextK = FALSE
    while(count <= maxit){
      newWeight = colSums(matCoef) / sum(matCoef)
      newMean = colSums(matCoef * x) / colSums(matCoef) -
        mTruncNrml(w.a-mixMean, w.b - mixMean, p = 1, mu = 0, sd = sqrt(mixVar))
      # newMean = sort(newMean)
      newVar =
        colSums(matCoef * (x - rep(newMean, each = n))^2) / colSums(matCoef) +
        mixVar -
        mTruncNrml(w.a-mixMean, w.b - mixMean, p = 2, mu = 0, sd = sqrt(mixVar))
      # if any mixture component has no points, restart with fewer k
      if(anyNA(newVar)){
        flg.nextK = TRUE
        break
      }
      # if any component collapses, inflat its variance
      if(any(newVar < minVar)){
        newVar[newVar < minVar] = mean(newVar[newVar>=minVar])
      }
      # sort the components before comparison
      ord.idx <- order(newMean)
      newWeight <- newWeight[ord.idx]
      newMean <- newMean[ord.idx]
      newVar <- newVar[ord.idx]
      # compare changes
      if(sum((newWeight - mixWeight)^2 + (newMean - mixMean)^2 + (newVar - mixVar)^2) < tol){
        break
      }
      mixWeight = newWeight
      mixMean = newMean
      mixVar = newVar
      for(j in 1:k){
        matCoef[,j] = mixWeight[j] * dnorm(x, mean = mixMean[j], sd = sqrt(mixVar[j]))
      }
      tmCoef = (
        pnorm(w.b, mean = mixMean, sd = sqrt(mixVar)) -
          pnorm(w.a, mean = mixMean, sd = sqrt(mixVar))
      )
      matCoef = matCoef / rep(tmCoef, each = n)
      matCoef = matCoef / rowSums(matCoef)
      count = count + 1
    }
    if(flg.nextK){
      k = k - 1
    }else{
      flg.done = TRUE
    }
  }

  # the following steps are needed since the EM calculation use
  # truncation before mixing (i.e., mixture of truncated Gaussian)
  # yet the return are parameterized as truncation after mixing.
  # compute normalizing constant
  tmCoef = pnorm(w.b, mean = mixMean, sd = sqrt(mixVar)) -
    pnorm(w.a, mean = mixMean, sd = sqrt(mixVar))
  # recovering the original weight
  mixWeight = mixWeight / tmCoef
  mixWeight = mixWeight / sum(mixWeight)

  # returning
  if (control$return.scale == 'parameter') {
    return(list(
      weight = mixWeight[order(mixMean)],
      mean = mixMean[order(mixMean)],
      var = mixVar[order(mixMean)],
      num.mixture = k,
      iter = count
    ))
  }
  if (control$return.scale == 'optim') {
    ls.par <- list(
      weight = rep(0, num.mixture),
      mean = rep(0, num.mixture),
      var = rep(1, num.mixture)
    )
    ord.idx <- order(mixMean)
    # work around: dropped component -> mixweight = 0
    ls.par$weight[seq_len(k)] <- mixWeight[ord.idx]
    ls.par$mean[seq_len(k)] <- mixMean[ord.idx]
    ls.par$var[seq_len(k)] <- mixVar[ord.idx]
    return(list(
      par = tMixGauss.repar(ls.par = ls.par),
      value = NA,
      convergence = 1,
      iter = count,
      message = paste(
        'this is a dummy optim list constructed from EM results,',
        'only the par and iter components are useful'
      )
    ))
  }
  # density related
  grid  <- control$grid
  tm.mat <- matrix(0, ncol = k, nrow = length(grid))
  for(j in 1:k){
    tm.mat[,j] <- mixWeight[j] * dnorm(grid, mean = mixMean[j], sd = sqrt(mixVar[j]))
  }
  # pdf
  if (control$return.scale == 'origin') {
    return(rowSums(tm.mat) / sum(tmCoef * mixWeight))
  }
  # log-pdf
  if (control$return.scale == 'log') {
    return(log(rowSums(tm.mat) / sum(tmCoef * mixWeight)))
  }
}


# LEGACY, not needed anymore, calMLE can handle this.
mleTruncMixGauss = function(x, k, edge, returnDensity = FALSE){
  # multistart MLE wrapper for truncated mixture
  # ONE SHALL NOT USE
  # the optimization here is not very careful

  if(length(edge) == 1){
    w.a = -abs(edge)
    w.b = abs(edge)
  }else{
    w.a = min(edge)
    w.b = max(edge)
  }

  objf = function(par){
    w.par = matrix(par, nrow = 3)
    w.par[1,] = w.par[1,]^2 / sum(w.par[1,]^2)

    nrml.const = pnorm(w.b, mean = w.par[2,], sd = exp(w.par[3,])) - pnorm(w.a, mean = w.par[2,], sd = exp(w.par[3,]))
    # adjusted weight to accommodate truncate after mixing
    w.par[1,] = w.par[1,] * nrml.const
    w.par[1,] = w.par[1,] / sum(w.par[1,])
    mat.Comp = sapply(x, function(obsv){
      w.par[1,] * dnorm(obsv, mean = w.par[2,], sd = exp(w.par[3,])) / nrml.const
    })
    return(
      -1*sum(log(colSums(mat.Comp)))
    )
  }
  n = length(x)

  res.km = kmeans(x, centers = k)
  mixWeight = as.numeric(table(res.km$cluster) / n)
  mixMean = as.numeric(res.km$centers)
  mixVar = sapply(1:k, function(idx){var(x[res.km$cluster == idx])})

  init_par = as.numeric(t(matrix(c(mixWeight, mixMean, log(mixVar)/2), ncol = 3)))
  # init_par[c(3,6,9)] = 0

  # invisible({
  #   tst = densityFPCA::smartOptim(init_par, objf, multiStart = c(50, 4))
  #   })
  multiStart = matrix(0, nrow = 3, ncol = k)
  multiStart[1,] = 2
  multiStart[2,] = (w.b-w.a)/2
  multiStart[3,] = 1.5 * IQR(log(mixVar)/2)
  multiStart = c(min(3*k, 10), as.numeric(multiStart))

  opt_res_multi = smartOptim(init_par, objf, multiStart = multiStart)
  # sorting
  res.tm = matrix(opt_res_multi$par, nrow = 3)
  res.tm = res.tm[,order(res.tm[2,])]

  mixWeight = res.tm[1,]^2 / sum(res.tm[1,]^2)
  mixMean = res.tm[2,]
  mixVar = exp(res.tm[3,])^2

  if(is.logical(returnDensity)){
    if(returnDensity) stop('Need to specify grid as returnDensity.')
    return(list(
      weight = mixWeight[order(mixMean)],
      mean = mixMean[order(mixMean)],
      var = mixVar[order(mixMean)]
    ))
  }else{
    # if one wants density
    grid = returnDensity
    if(anyNA(grid) | any(!is.finite(grid)) | any(is.null(grid))) stop('Non finite grid value.')
    tm.mat = matrix(0, ncol = k, nrow = length(grid))
    # compute normalizing constant
    tmCoef = pnorm(w.b, mean = mixMean, sd = sqrt(mixVar)) - pnorm(w.a, mean = mixMean, sd = sqrt(mixVar))
    for(j in 1:k){
      tm.mat[,j] = mixWeight[j] * dnorm(grid, mean = mixMean[j], sd = sqrt(mixVar[j]))
    }
    return(rowSums(tm.mat) / sum(tmCoef * mixWeight))
  }
}

# constant list of non-fpca based methods.
kListOtherMethod <- list(
  MLE = calMLE,
  emMixGauss = emMixGauss,
  emTruncMixGauss = emTruncMixGauss,
  mleTruncMixGauss = mleTruncMixGauss,
  KDE = stats::density
)
