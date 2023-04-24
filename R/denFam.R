## family declarations and some tools

# model based inference =======================================================

# TBD
# selectK.AIC <- function(obsv, fpca.res){
#   # select number of components for new observation with fitted family
#   # args:
#   #  obsv: a vector of new observations
#   #  fpca.res: fdapace::FPCA result
#   # returns:
#   #  a list of AIC and optimum k
#   fpca.k <- fpca.res$
#   getNegLogll
# }
# END TBD

# waldCI <- function(den.fam, obsv, mle, true.par){
#   # compute Wald Confidence interval of density estimation
#   # args:
#   # den.fam: density family
#   # obsv: observations;
#   # mle: MLE of the parameters;
#   # true.par(optional): true parameter value.
#   # returns:
#   # a function taking par & alpha = 0.95 computing conf. band.
#   # Detail:
#
#   if (missing(true.par)) {
#     true.par <- mle
#   }
#   # get observed information.
#   obsv.info <- empFisherInfo(obsv, den.fam)
#   mle + c(-1, 1) *
#
#   # THREAD OF CONSTRUCTION
#
#
# }

deltaMethodCI <- function(stat, theta.hat, theta.var, level){
  # compute confidence interval of stat with delta method
  # args:
  # stat: a scaler function;
  # theta.hat: argument of stat;
  # theta.var: variance or covariance matrix of theta.hat, should be small;
  # level(optional): confidence level.
  # returns:
  # 2 number or, if level missing, a function taking level and gives CI.

  h <- stat
  f.h.gr <- function(par) as.numeric(nloptr::nl.grad(par, h))
  f.h.gr2 <- function(par) nloptr::nl.jacobian(par, f.h.gr)
  h.gr <- f.h.gr(theta.hat)

  if (all(abs(h.gr) <= .Machine$double.eps^(1/3))) {
    h.gr2 <- f.h.gr2(theta.hat)

    smpl.gChisq <- rgChisq(
      50000, rep(0, length(theta.hat)),
      theta.var, h.gr2 / 2
    )

    res.f <- function(level){
      alpha <- 1 - level
      q.chi <- as.numeric(
        quantile(smpl.gChisq, probs = c(alpha / 2, 1 - alpha / 2))
      )
      return(
        stat(theta.hat) + q.chi
      )
    }
  } else {
    res.f <- function(level){
      alpha <- 1 - level
      w.var <- quadForm(theta.var, h.gr)
      q.norm <- qnorm(1 - alpha / 2)
      return(
        stat(theta.hat) + c(-1, 1) * rep((q.norm * sqrt(w.var)), 2)
      )
    }
  }

  if (missing(level)) {
    return(res.f)
  } else{
    return(res.f(level))
  }
}

fisherInfo <- function(den.fam){
  # compute the Fisher Information (theoretical).
  # args:
  # den.fam = density family;
  # returns:
  # a function computing the Fisher Information at par.

  res.f <- function(par){
    grid <- seq(den.fam$range[1], den.fam$range[2], length.out = 1000)
    vec.pdf <- den.fam$pdf(grid, par)
    # score function (1st derivative of log-pdf)
    if (is.null(den.fam$logpdf.gr)) {
      tm.f <- function(par) dotL2(log(den.fam$pdf(grid, par)), vec.pdf, grid = grid)
      score.f <- function(par) nloptr::nl.grad(par, tm.f)
    } else {
      score.f <- function(par){
        tm.mat <- den.fam$logpdf.gr(grid, par)
        res <- rep(0, ncol(tm.mat))
        for(i in seq_len(ncol(tm.mat))){
          res[i] <- dotL2(tm.mat[, i], vec.pdf, grid = grid)
        }
        return(res)
      }
    }
    # 2nd derivative of log-pdf
    return(-1 * nloptr::nl.jacobian(par, score.f))
  }
  return(res.f)
}

empFisherInfo <- function(obsv, den.fam){
  # compute the Fisher Information (empirical).
  # args:
  # den.fam = density family;
  # obsv = array of observation from den.fam given some parameter value.
  # returns:
  # function taking par and computes the empirical Fisher Information.

  # 1st derivative of negative log-likelihood
  if (is.null(den.fam$logpdf.gr)) {
    neglogll <- getNegLogll(obsv, den.fam)
    neglogll.gr <- function(par) nloptr::nl.grad(par, neglogll)
  } else {
    neglogll.gr <- getNegLogll.gr(obsv, den.fam)
  }

  neglogll.d3 <- function(par) nloptr::nl.jacobian(par, neglogll.gr)

  # NEED CORRECTION
  f.return <- function(par){
    return(
      neglogll.d3(par)  / length(obsv)
    )
  }
}

# fpca related ================================================================
fpcaEsti2pdf <- function(fpca.res, mat.par, num.k, grid){
  # wrapper for computing pdf on a grid w.r.t multiple estimates from fpca.res
  # args:
  # fpca.res: a result from fdapace::FPCA
  # mat.par: matrix where first few elements of a row is one set of parameters.
  # num.k: a nrow(mat.par) lenght array of numbers specifying numbers of
  #        eigenfuncs to use, if missing, use all, if length 1, recycle.
  # grid: on which to compute
  # returns:
  # a matrix whose row is one density curve

  # sanity
  if (!is.matrix(mat.par)) {
    mat.par <- matrix(mat.par, nrow = 1)
  }
  if (ncol(fpca.res$phi) < ncol(mat.par))
    stop('too much parameters for fpca.res, check ncol(mat.par)')
  if (missing(num.k))
    num.k <- min(ncol(fpca.res$phi), ncol(mat.par))
  if (length(num.k) == 1){
    den.fam <- fpca2DenFam(fpca.res, list(num.k = num.k))
    return(par2pdf(
      den.fam = den.fam,
      mat.par = mat.par[, seq_len(num.k), drop = FALSE],
      grid = grid
    ))
  }else{
    den.fam <- fpca2DenFam(fpca.res, list(num.k = max(num.k, ncol(mat.par))))
    w.mat.par <- matrix(0, nrow = nrow(mat.par), ncol = ncol(mat.par))
    for(i in seq_len(nrow(mat.par))){
      w.mat.par[i, seq_len(num.k[i])] <- mat.par[i, seq_len(num.k[i])]
    }
    return(par2pdf(
      den.fam = den.fam,
      mat.par = w.mat.par,
      grid = grid
    ))
  }

  # mat.pdf <- matrix(0, nrow = nrow(mat.par), ncol = length(grid))
  # for(i in seq_len(nrow(mat.par))){
  #   den.fam <- fpca2DenFam(fpca.res, list(num.k = num.k[i]))
  #   mat.pdf[i, ] <- den.fam$pdf(grid, mat.par[i, seq_len(num.k[i])])
  # }
  # return(mat.pdf)
}

#' Construct approximating family from FPCA result
#'
#' @param fpca.res `fdapace::FPCA` result.
#' @param control a list of controlling arguments:
#' - `num.k`: number of eigenfunc building the family, missing the use all;
#' - `extend.par.range`: `FALSE`(default), or positive numbers > 1 of length
#' one or `num.k`, determining how each dim of the `par.range` would be extended
#' from the `fpca.res$xiEst`;
#' - `check`: `FALSE`(default)/`TRUE`, whether to run `checkDenFamNumeric`;
#' - `check.sample`(optional): argument passed to `checkDenFamNumeric`, if
#' supplied, check will be set as `TRUE`;
#' - `get.prior`: `FALSE`(default)/`TRUE`, whether to compute prior information
#' from `fpca.res$xiEst`.
#'
#' @return a list of
#' - `pdf`: density;
#' - `logpdf.gr`: gradient of log density;
#' - `par.range`: a matrix for range of parameters (`xiEst`);
#' - `suff.f`: a function computes sufficient statistics.
#'
#' Note that no random generator included.
#' @export
#'
#' @examples
#' den.fam <- tNrml(2)
#' ls.obsv <- with(
#'   den.fam, apply(
#'     gen_par(100), 1, function(par) rpdf(rpois(1, 500), par),
#'     simplify = FALSE
#'   )
#' )
#' grid <- seq(min(den.fam$range), max(den.fam$range), length.out = 1024)
#' mat.pdf <- preSmooth.kde(ls.obsv, grid = grid, kde.opt = list(bw = 0.25))
#' ls.hilbert <- toHilbert(
#'   mat.pdf, grid,
#'   transFun = function(f, grid) orthLog(f, grid, against = 1 / 4),
#'   eps = .Machine$double.eps^(1/2)
#' )
#' fpca.res <- do.call(
#'   fdapace::FPCA, list(
#'     Ly = asplit(ls.hilbert$mat.curve, 1),
#'     Lt = replicate(nrow(ls.hilbert$mat.curve), expr = grid, simplify = FALSE),
#'     optns = list(lean = TRUE, methodSelectK = 2, useBinnedData = 'OFF')
#'   )
#' )
#' fpca.den.fam <- fpca2DenFam(fpca.res, control = list(num.k = 2))
#' matplot(grid, fpca.den.fam$suff.f(grid), type = 'l')
fpca2DenFam <- function(fpca.res, control = list()){
  # function constructing a density family list with FPCA result.
  # args:
  # fpca.res = fdapace::FPCA result
  # control = a list of controling arguments:
  #   num.k = number of eigenfunc building the family, missing -> use all;
  #   extend.par.range = FALSE(default) or one or num.k positive numbers > 1,
  #                      determining how each dim of the par.range would be
  #                      extended from the fpca.res$xiEst;
  #   check = FALSE(default)/TRUE, whether to run checkDenFamNumeric;
  #   check.sample = (optional) that of checkDenFamNumeric, if supplied, check
  #                  will be set as TRUE
  #   get.prior = FALSE(default)/TRUE, whether to compute prior information
  # return:
  # a list of pdf, logpdf.gr, par.range, the random generator parts ar omitted.
  # also a suff.f, a function: x -> sufficient statistics T(x).

  if (is.null(control$num.k)) {
    control$num.k <- ncol(fpca.res$phi)
  }
  if (control$num.k > ncol(fpca.res$phi)) {
    warning('Insufficient number of eigenfucntions, num.k reset.')
    control$num.k <- ncol(fpca.res$phi)
  }
  if (is.null(control$extend.par.range))
    control$extend.par.range <- FALSE
  if (!(length(control$extend.par.range) %in% c(1, control$num.k)))
    stop('check input length of extend.par.range')
  if (is.null(control$check))
    control$check <- FALSE
  if (is.null(control$check.sample)){
    control$check.sample <- NULL
  }else{
    control$check <- TRUE
  }
  if (is.null(control$get.prior))
    control$get.prior <- FALSE

  # range of parameters, just use range of xiEst for now
  par.range <- apply(
    fpca.res$xiEst[, seq_len(control$num.k), drop = FALSE], 2,
    range
    # function(x){
    #   c.x <- mean(range(x))
    #   span.x <- diff(range(x))
    #   r.x <- range(x)
    #   return(1.5 * span.x / 2 * (r.x - c.x) + c.x)
    # }
  )
  if(is.numeric(control$extend.par.range)){
    c.ex <- control$extend.par.range
    if(any(c.ex <= 0))
      stop('extend.par.range must be positive')
    mid.par.range <- colMeans(par.range)
    wid.par.range <- 0.5 * c.ex * (par.range[2, ] - par.range[1, ])
    par.range[1, ] <- mid.par.range - wid.par.range
    par.range[2, ] <- mid.par.range + wid.par.range
  }

  # prior information
  if(control$get.prior){
    # prepare quantities used in BLUP estimator to save computation
    ls.cov <- list()
    for(nk in seq_len(control$num.k)){
        ls.cov[[nk]] <- local({
          tm.cov <- condVarT(
            fpca.res = fpca.res, e.co = fpca.res$xiEst[, seq_len(nk), drop = FALSE],
            return.all = TRUE
          )
          E_mCo <- colMeans(tm.cov$m.co)
          cov_mCo <- cov(tm.cov$m.co)

          # E(Var(T|natural par))
          E_cond_cov_T <- matrix(0, nk, nk)
          E_cond_cov_T[lower.tri(E_cond_cov_T, diag = T)] <- # must be lower.tri, R index row first
            colMeans(tm.cov$cond.cov)
          E_cond_cov_T <- E_cond_cov_T + t(E_cond_cov_T)
          diag(E_cond_cov_T) <- diag(E_cond_cov_T) / 2
          list(
            m.co = tm.cov$m.co,
            E_mCo = E_mCo,
            cov_mCo = cov_mCo,
            E_cond_cov_T = E_cond_cov_T
          )
        })
      }
      names(ls.cov) <- seq_len(control$num.k)
  }else{
    ls.cov <- NULL
  }

  local.nrmlConst <- function(par) nrmlConst(par, fpca.res, control$num.k)
  suff.f <- function(x){
    # return: a length(x) X num.k matrix
    mat.x <- matrix(0, ncol = control$num.k, nrow = length(x))
    for(i in seq_len(control$num.k)){
      mat.x[, i] <- approx(
        x = fpca.res$workGrid,
        y = fpca.res$phi[, i],
        xout = x,
        rule = 2
      )$y
    }
    return(mat.x)
  }

  pdf <- function(x, par){
    # the density function
    # return a array where ith element = pdf(x[i], par)
    nrml.const <- local.nrmlConst(par)
    mat.x <- suff.f(x)
    # mat.x <- matrix(0, ncol = control$num.k, nrow = length(x))
    # for(i in seq_len(control$num.k)){
    #   mat.x[, i] <- approx(
    #     x = fpca.res$workGrid,
    #     y = fpca.res$phi[, i],
    #     xout = x,
    #     rule = 2
    #   )$y
    # }
    only.x <- approx(
      x = fpca.res$workGrid,
      y = fpca.res$mu,
      xout = x,
      rule = 2
    )$y

    return(
      as.numeric(exp(only.x + mat.x %*% par - nrml.const))
    )
  }
  logpdf.gr <- function(x, par){
    # the gradient of log-density
    # args:
    # x = observation;
    # par = parameters;
    # returns:
    # a matrix where ith row = gradient of log-pdf(x[i], par)
    nrml.const.gr <- nloptr::nl.grad(par, local.nrmlConst)
    mat.gr <- t(suff.f(x))
    mat.gr <- t(mat.gr - nrml.const.gr)
    # mat.gr <- matrix(0, ncol = control$num.k, nrow = length(x))
    # for(i in seq_len(control$num.k)){
    #   mat.gr[, i] <- approx(
    #     x = fpca.res$workGrid,
    #     y = fpca.res$phi[, i],
    #     xout = x,
    #     rule = 2
    #   )$y - nrml.const.gr[i]
    # }
    return(mat.gr)
  }
  if(control$check){
    num.check <- do.call(
      checkDenFamNumeric,
      list(
        den.fam = list(
          pdf = pdf,
          logpdf.gr = logpdf.gr,
          suff.f = suff.f,
          par.range = par.range,
          range = range(fpca.res$workGrid)
        ),
        check.sample = control$check.sample
      )
    )
    if(!is.null(num.check))
      warning('checkDenFamNumeric fails, optimization may not work')
  }

  moment <- function(r, par, grid){
    # compute moment EX^r with pdf
    # r: order
    # par: parameters, if shorter, fill with 0
    # grid: the grid to work with, if missing, use range(fpca.res$workGrid)
    if(length(par) < control$num.k)
      par <- c(par, rep(0, control$num.k - length(par)))
    if(base::missing(grid)){
      res <- dotL2(
        f1 = function(x) x^r,
        f2 = function(x) pdf(x, par)
        , range = range(fpca.res$workGrid)
      )
    }else{
      res <- dotL2(
        f1 = grid ^ r,
        f2 = pdf(grid, par)
        , grid = grid
      )
    }
    return(res)
  }
  fill.par <- function(par){
    # fill par to length = control$num.k with 0
    n.par <- control$num.k

    if(is.matrix(par)){
      stopifnot(ncol(par) <= n.par)
      res <- matrix(0, nrow = nrow(par), ncol = n.par)
      res[, seq(ncol(par))] <- par
    }else{
      stopifnot(length(par) <= n.par)
      res <-
        c(par, rep(0, n.par - length(par)))
    }
    return(res)
  }

  return(
    list(
      pdf = pdf,
      logpdf.gr = logpdf.gr,
      suff.f = suff.f,
      moment = moment,
      par.range = par.range,
      range = range(fpca.res$workGrid),
      prior = ls.cov,
      fill.par = fill.par
    )
  )
}

fpca2ObjFun <- function(obsv, fpca.res, control = list()){
  # obtaining the part of negLogll related to parameters, dropping irrelavents.
  # args:
  # obsv = array of obseration;
  # fpca.res = fdapace::FPCA result
  # control = a list of controling arguments:
  #   num.k = number of eigenfunc building the family, missing -> use all;
  #   ... TBD
  # return:
  # a list of neglogll(part of) and its gradient.

  if (is.null(control$num.k)) {
    control$num.k <- ncol(fpca.res$phi)
  }
  if (control$num.k > ncol(fpca.res$phi)) {
    warning('Insufficient number of eigenfucntions, num.k reset.')
    control$num.k <- ncol(fpca.res$phi)
  }

  local.nrmlConst <- function(par) nrmlConst(par, fpca.res, control$num.k)

  mat.x.t <- matrix(0, nrow = control$num.k, ncol = length(obsv))
  for(i in seq_len(control$num.k)){
    mat.x.t[i, ] <- approx(
      x = fpca.res$workGrid,
      y = fpca.res$phi[, i],
      xout = obsv,
      rule = 2
    )$y
  }
  s.x <- rowSums(mat.x.t)

  pNegLogll <- function(par){
    # the negLogll function
    # returns: a number.
    nrml.const <- local.nrmlConst(par)
    return(
      -1 * sum(s.x * par) + nrml.const * length(obsv)
    )
  }

  pNegLogll.gr <- function(par){
    # the gradient of negative loglikelihood
    # args:
    # par = parameters;
    # returns:
    # an array of gradient.
    nrml.const.gr <- nloptr::nl.grad(par, local.nrmlConst)
    # mat.gr <- matrix(0, ncol = control$num.k, nrow = length(obsv))
    # for(i in seq_len(control$num.k)){
    #   mat.gr[, i] <- approx(
    #     x = fpca.res$workGrid,
    #     y = fpca.res$phi[, i],
    #     xout = obsv,
    #     rule = 2
    #   )$y
    # }
    return(
      -1 * s.x + nrml.const.gr * length(obsv)
    )
  }

  return(
    list(
      pNegLogll = pNegLogll,
      pNegLogll.gr = pNegLogll.gr
    )
  )
}

nrmlConst <- function(par, fpca.res, k){
  # computing normalizing constant
  # args:
  # par = parameters;
  # fpca.res = result from fdapace::FPCA;
  # k = number of eigenfunc to use.
  # return:
  # a number.
  # Detail:
  # nrmlConst = \log( \int \exp( mu(x) + \sum_i par_i * phi_i(x) ) dx )
  return(
    log(dotL2(
      exp(fpca.res$mu),
      exp(fpca.res$phi[, seq_len(k), drop = FALSE] %*% par),
      grid = fpca.res$workGrid
    ))
  )
}


# Negative loglikelihood related ==============================================
getNegLogll <- function(obsv, den.fam){
  # obtaining negative loglikelihood function
  # args:
  # obsv = the array of observations;
  # den.fam = density family, a list with element pdf,
  # it must be vectorized w.r.t. observation.
  # returns:
  # a negative loglikelihood function taking parameters.

  neg.logll <- function(par){
    -1 * sum(log(
      den.fam$pdf(obsv, par)
    ))
  }

  return(neg.logll)
}

getNegLogll.gr <- function(obsv, den.fam){
  # obtaining gradient of negative loglikelihood function
  # args:
  # obsv = the array of observations;
  # den.fam = density family, a list with element logpdf.gr,
  # it must be vectorized w.r.t. observation.
  # returns:
  # the gradient of negative loglikelihood function taking parameters.
  if (is.null(den.fam$logpdf.gr)) {
    neg.logll.gr <- NULL
  } else{
    neg.logll.gr <- function(par){
      -1 * colSums(
        den.fam$logpdf.gr(obsv, par)
      )
    }
  }
  return(neg.logll.gr)
}

# for visually inspecting density families ====================================
sketchDenFam = function(denFam, ylim, n=10, scale = 'origin', denFam.name){
  # denFam.name used for plotting
  # scale = log or orgin

  if(missing(denFam.name)) denFam.name = NULL
  demo_par = denFam$gen_par(n)
  grid = seq(denFam$range[1], denFam$range[2], length.out = 1000)
  mat_pdf = matrix(0, nrow = n, ncol = length(grid))
  for(i in 1:n){
    mat_pdf[i, ] =  denFam$pdf(grid, par = demo_par[i,])
  }
  if(scale == 'log'){
    mat_pdf = log(mat_pdf)
    if(any(abs(mat_pdf) == Inf)) warning('Inf in log-density.')
  }
  if(missing(ylim)) ylim = range(mat_pdf, finite = TRUE)
  plot(x = grid, y = mat_pdf[1,], type = 'l', col = rgb(1,0,0,0.25),
       ylim = ylim,
       xlab = 'domain', ylab = 'density',
       main = denFam.name)
  for(i in 2:n){
    points(grid, mat_pdf[i,], type = 'l', col = rgb(1,0,0,0.25))
  }
}

checkSampler <- function(den.fam, size, n = 10, what = 'diff', kde.option = list()) {
  # plot diff between analytic pdf and KDE from generated sample.
  # args:
  # den.fam = density family;
  # size = size of individual sample;
  # n = number of samples.
  # what = what to plot,
  #   if "diff": difference between KDE and true
  #   if "kde":  KDE
  #   if "log-kde": log-KDE
  #   if "pdf": plot pdf
  # kde.option = list of stats::density options.
  # returns:
  # NULL but plot.

  # handle input
  what <- match.arg(what, c('diff', 'kde', 'log-kde', 'pdf'))
  # kde.option
  if (is.null(kde.option$bw)) kde.option$bw <- 'SJ'
  if (is.null(kde.option$kernel)) kde.option$kernel <- 'gaussian'

  mat.par <- den.fam$gen_par(n)
  grid <- seq(min(den.fam$range), max(den.fam$range), length.out = 500)
  mat.pdf <- matrix(0, ncol = length(grid), nrow = n)
  mat.sample <- matrix(0, ncol = size, nrow = n)
  mat.kde <- matrix(0, ncol = length(grid), nrow = n)
  mat.diff <- matrix(0, ncol = length(grid), nrow = n)
  for (i in seq_len(n)) {
    mat.sample[i, ] <-  den.fam$rpdf(ncol(mat.sample), mat.par[i, ])
    tm.kde <- do.call(
      stats::density,
      c(list(x = mat.sample[i, ]), kde.option)
    )
    mat.kde[i, ] <- approx(x = tm.kde$x, y = tm.kde$y, xout = grid, rule = 1)$y
    mat.pdf[i, ] <- den.fam$pdf(grid, mat.par[i, ])
    mat.diff[i, ] <- mat.pdf[i, ] - mat.kde[i, ]
  }
  if (what == 'diff') {
    matplot(x = grid, y = t(mat.diff), lty = 1, type = 'l', col = rgb(1,0,0,0.2), lwd = 1,
            ylab = 'True - KDE')
  }
  if (what == 'kde') {
    matplot(x = grid, y = t(mat.kde), lty = 1, type = 'l', col = rgb(1,0,0,0.2), lwd = 1,
            ylab = 'KDE')
  }
  if (what == 'log-kde') {
    log.mat.kde <- log(mat.kde)
    if (!all(is.finite(log.mat.kde))){
      num.inf <- sum(rowSums(!is.finite(log.mat.kde)) != 0)
      warning(paste(num.inf, 'inf after log.'))
    }
    matplot(x = grid, y = t(log.mat.kde), lty = 1, type = 'l', col = rgb(1,0,0,0.2), lwd = 1,
            ylab = 'logKDE')
  }
  if (what == 'pdf') {
    matplot(x = grid, y = t(mat.pdf), lty = 1, type = 'l', col = rgb(1,0,0,0.2), lwd = 1,
            ylab = 'PDF')
  }
}

# for numerically inspecting density families =================================

#' check whether negative log-likelihood and its gradient maintain
#' sensible values on vertices of parameter range.
#'
#' @param den.fam density family.
#' @param check.sample (optional) sample used to compute loglikelihood,
#' recommended if `den.fam` does not have random generator (`rpdf`).
#'
#' @return `NULL` or a list of sample used and a matrix of problematic points
#' and values. Cannot handle high dimensional parameters.
#' samples will be generated with randomly chosen parameter values.
#' @details If `den.fam` does not provide random generator `rpdf`, use
#' `Runuran` to draw sample based on density function `pdf`, which may be slow.
#' @export
#'
#' @examples
#' checkDenFamNumeric(tNrml(3))
#' den.fam <- tMixGauss(6)
#' den.fam$par.range[, 1] <- c(-10, 10)
#' checkDenFamNumeric(den.fam)
#' den.fam$par.range[, 1] <- c(-10.5, 11)
#' checkDenFamNumeric(den.fam)
checkDenFamNumeric <- function(den.fam, check.sample){
  # check whether negative log-likelihood and its gradient
  # maintain sensible values on vertices of par.range.
  # args:
  # den.fam = density family;
  # check.sample = (optional) sample used to compute loglikelihood
  #                recommanded if den.fam does not have rpdf
  # returns:
  # NULL or a list of sample used and a matrix of pts and problematic values.
  # Cannot handle more than 15 parameters
  # details:
  #   sample will be generated with a randomly chosen par value.

  # n.grid <- 100
  # par.dim <- ncol(den.fam$par.range)
  # cdd.pts <- matrix(0, ncol = n.grid, nrow = 2 * par.dim * n.grid)
  # cdd.edge <- matrix(0, ncol = par.dim, nrow = n.grid)
  # for(i in seq_len(ncol(cdd.edge))){
  #   cdd.edge[, i] <- seq(
  #     from = den.fam$par.range[1, i],
  #     to = den.fam$par.range[2, i],
  #     length.out = n.grid
  #   )
  # }

  par.dim <- ncol(den.fam$par.range)
  if(par.dim > 15) stop('Too many parameters.')

  # n.pts <- 2 ^ par.dim
  # cdd.pts <- matrix(0, ncol = par.dim, nrow = n.pts)
  # tm.arr.1 <- TRUE
  # tm.arr.2 <- FALSE
  # for(j in seq_len(ncol(cdd.pts))){
  #   tm.arr <- c(tm.arr.1, tm.arr.2)
  #   cdd.pts[, j] <- rep(tm.arr, length.out = n.pts)
  #   tm.arr.1 <- rep(tm.arr.1, 2)
  #   tm.arr.2 <- rep(tm.arr.2, 2)
  # }

  # ls.par <- as.list(as.data.frame(den.fam$par.range))
  # ls.par$stringsAsFactors = FALSE
  ls.par <- mat2list(den.fam$par.range, by = 2)
  cdd.pts <- as.matrix(do.call(expand.grid, ls.par))

  val.neglogll <- rep(0, nrow(cdd.pts))
  mat.gr <- matrix(0, ncol = ncol(den.fam$par.range), nrow = nrow(cdd.pts))
  # tst.par <- colMeans(den.fam$par.range)
  tst.par <- runif(
    ncol(den.fam$par.range),
    min = den.fam$par.range[1, ],
    max = den.fam$par.range[2, ]
  )
  if(missing(check.sample)){
    if(is.null(den.fam$rpdf)){
      warning('generating random sample with Runuran(slow)')
      gen <- Runuran::pinv.new(
        pdf = den.fam$pdf,
        lb = min(den.fam$range),
        ub = max(den.fam$range),
        center = mean(den.fam$range),
        par = tst.par
      )
      tst.obsv <- Runuran::ur(gen, 50)
    }else{
      tst.obsv <- den.fam$rpdf(50, tst.par)
    }
  }else{
    tst.obsv <- as.numeric(check.sample)
  }
  # get negative loglikelihood and evaluate
  neglogll <- getNegLogll(tst.obsv, den.fam)
  neglogll.gr <- ifelse(
    is.null(den.fam$logpdf.gr),
    function(par) nloptr::nl.grad(par, neglogll),
    getNegLogll.gr(tst.obsv, den.fam)
  )
  for(i in seq_len(nrow(mat.gr))){
    val.neglogll[i] <- neglogll(cdd.pts[i, ])
    mat.gr[i, ] <- neglogll.gr(cdd.pts[i, ])
  }
  # any nonfinite values in the neglogll or its derivatives?
  idx.err <- rep(0, length(val.neglogll))
  for(i in seq_len(length(val.neglogll))){
    idx.err[i] <- any(!is.finite(mat.gr[i, ]))
  }
  idx.err <- idx.err | !is.finite(val.neglogll)
  if (any(idx.err)) {
    mat.res <- cbind(
      cdd.pts[idx.err, , drop = FALSE],
      val.neglogll[idx.err],
      mat.gr[idx.err, , drop = FALSE]
    )
    colnames(mat.res) <- c(
      paste0("par.", seq_len(ncol(cdd.pts))),
      "neglogll",
      paste0("gr.", seq_len(ncol(cdd.pts)))
    )
    return(
      list(
        obsv = tst.obsv,
        pts = mat.res
      )
    )
  } else {
    return(NULL)
  }
}
# tst.den.fam <- tMixGauss(6)
# tst.den.fam$par.range[, 1] <- c(-10, 10)
# checkDenFamNumeric(tst.den.fam)
# tst.den.fam$par.range[, 1] <- c(-10.5, 11)
# checkDenFamNumeric(tst.den.fam)

#' Check if the gradient of loglikelihood is correct.
#'
#' @param den.fam density family;
#' @param tol tolerence of relative error;
#'
#' @return `NULL` if all good, or a list of problematic points and gradients.
#' @export
#'
#' @examples # See vignette.
checkDenFamGrident <- function(den.fam, tol = 1e-4){
  # for density family, check if the analytic gradient is correct, if provided.
  # args:
  # den.fam = density family;
  # tol = tolerence of relative error of diff / finite.diff.grad
  # returns:
  # nothing (NULL) if all good, or a list of points and gradients having problem.
  # remark:
  # it seems nloptr::check.derivatives is a bit problematic, uses homemade one.

  if (is.null(den.fam$logpdf.gr)) {
    message('Analytic gradient not provided, nothing to check.')
    return(NULL)
  }

  if(is.null(den.fam$gen_par)) {
    den.fam$gen_par <- function(n) {
      t(replicate(
        n, apply(den.fam$par.range, 2, function(x) runif(1, x[1], x[2]))
      ))
    }
  }
  true.par <- as.numeric(den.fam$gen_par(1))
  if(is.null(den.fam$rpdf)) {
    warning('generating random sample with Runuran(slow)')
    den.fam$rpdf <- function(n, par) {
      gen <- Runuran::pinv.new(
        pdf = den.fam$pdf,
        lb = min(den.fam$range),
        ub = max(den.fam$range),
        center = mean(den.fam$range),
        par = par
      )
      Runuran::ur(gen, n)
    }
  }
  obsv <- den.fam$rpdf(50, true.par)

  negLogll <- getNegLogll(obsv, den.fam)
  negLogll.gr <- getNegLogll.gr(obsv, den.fam)

  tst.par <- den.fam$gen_par(100)
  ls.res <- list()
  idx.err <- rep(FALSE, nrow(tst.par))
  for(i in seq_len(nrow(tst.par))){
    ls.res[[i]] <- checkDerivatives(
      tst.par[i, ],
      negLogll,
      negLogll.gr,
      tol
    )
    ls.res[[i]]$x <- tst.par[i, ]
    if (any(ls.res[[i]]$flag_derivative_warning)) {
      idx.err[i] <- TRUE
    }
  }

  if (any(idx.err)) {
    return(list(
      true.par = true.par,
      obsv = obsv,
      err.pts = ls.res[idx.err]
    ))
  } else {
    return(NULL)
  }
}

# utility functions ===========================================================

wrapDenFam = function(denFam){
  # put every thing to [0,1] (i.e., transform domain to [0,1])
  # TBD: add quantile function for MKE
  old.range = denFam$range
  mid.pt = mean(old.range)
  sclcoef = sum(range(old.range) * c(-0.5, 0.5))
  transMap = function(x){(x*2-1) * sclcoef + mid.pt}
  to01 = function(x) ((x - mid.pt)/sclcoef + 1)/2
  return(
    list(
      pdf = function(x, par){denFam$pdf(transMap(x), par) * sclcoef * 2},
      rpdf = function(n, par){to01(denFam$rpdf(n, par))},
      gen_par = denFam$gen_par,
      range = c(0,1)
    )
  )
}

r.truncate <- function(n, rpdf, edge, ...){
  # truncating r.v. generator
  # n = desired sample size to generate
  # edge = truncation window
  # ... = arguments for rpdf
  if (length(edge) == 1){
    w.a <- -1 * abs(edge)
    w.b <- abs(edge)
  } else {
    w.a <- min(edge)
    w.b <- max(edge)
  }
  cnt <- 0
  o_cnt <- 0
  acpt <- 1
  res <- rep(0, n)
  while(cnt < n){
    tm <- rpdf(n = floor((n-cnt) / acpt)+1, ...)
    pickIdx <- (tm >= w.a) & (tm <= w.b)
    acpt <- sum(pickIdx)
    if(acpt == 0){
      acpt <- 0.5
      next
    }
    cnt <- cnt + acpt
    res[(o_cnt + 1) : cnt] = tm[pickIdx]
    o_cnt <- cnt
    acpt <- acpt / length(pickIdx)
  }
  return(res[1:n])
}

#' wrapper for computing pdf on a grid w.r.t multiple sets of parameter
#'
#' @param den.fam density family.
#' @param mat.par matrix, one row is one set of parameters.
#' @param grid on which to compute.
#' @param n.rep number of repeats, duplicates the resulting density matrix.
#'
#' @return a matrix, one row is one density curve.
#' @export
#'
#' @examples
#' den.fam <- tNrml(3)
#' grid <- seq(-3, 3, 0.05)
#' matplot(grid, t(par2pdf(den.fam, den.fam$gen_par(3), grid)), type = 'l')
par2pdf <- function(den.fam, mat.par, grid, n.rep = 1){
  # wrapper for computing pdf on a grid w.r.t multiple parameters
  # args:
  # density.family: density family.
  # mat.par: matrix whose row is one set of parameters.
  # grid: on which to compute
  # n.rep: number of repeats, sorely for assErr.
  # returns:
  # a matrix whose row is one density curve

  if (!is.matrix(mat.par)) {
    mat.par <- matrix(mat.par, nrow = 1)
  }
  mat.pdf <- matrix(0, ncol = nrow(mat.par), nrow = length(grid))
  for(i in seq_len(nrow(mat.par))){
    mat.pdf[, i] <- den.fam$pdf(grid, mat.par[i, ])
  }
  # apply(mat.par, 1, function(par) den.fam$pdf(grid, par))
  mat.pdf <- rep(mat.pdf, times = n.rep)
  mat.pdf <- matrix(mat.pdf, nrow = length(grid))
  return(t(mat.pdf))
}


# static families ======================================================
# Normal family
nrml_family <- list(
  pdf = function(x, par){dnorm(x, par[1], (par[2]))},
  rpdf = function(n, par){rnorm(n, par[1], (par[2]))},
  gen_par = function(n){return(matrix(c(runif(n, -5,5),
                                       (runif(n,1,5))), ncol = 2))},
  range = c(-20,20),
  par.range = matrix(c(-5, 5, 1, 5), 2, 2)
)
# Truncated Normal family
tNrml_family = list(
  pdf = function(x, par){dnorm(x, par[1], (par[2]))/(
                                              pnorm(1, par[1], (par[2])) -
                                                pnorm(0, par[1], (par[2])) ) },
  rpdf = function(n, par){
   cnt = 0
   o_cnt = 0
   acpt = 1
   res = rep(0, n)
   while(cnt < n){
     tm = rnorm(floor((n-cnt) / acpt)+1, par[1], (par[2]))
     pickIdx = (tm >= 0) & (tm <= 1)
     acpt = sum(pickIdx)
     if(acpt == 0){
       acpt = 0.5
       next
     }
     cnt = cnt + acpt
     res[(o_cnt + 1) : cnt] = tm[pickIdx]
     o_cnt = cnt
     acpt = acpt / length(pickIdx)
   }
   return(res[1:n])
   },
  gen_par = function(n){return(matrix(c(runif(n, -5,5),
                                       (runif(n,1,5))), ncol = 2))},
  range = c(0, 1),
  par.range = matrix(c(-5, 5, 1, 5), 2, 2)
)

# Beta family
# slightly truncated to stay away from 0 and 1.
bRange = c(1e-3, 1-1e-3)
beta_family = list(
  pdf = function(x, par){
    dbeta(x, exp(par[1]), exp(par[2]))/
      (pbeta(bRange[2], exp(par[1]), exp(par[2])) -
         pbeta(bRange[1], exp(par[1]), exp(par[2])))
    },
  rpdf = function(n, par){
    cnt = 0
    o_cnt = 0
    acpt = 1
    res = rep(0, n)
    while(cnt < n){
      tm = rbeta(floor((n-cnt) / acpt)+1, exp(par[1]), exp(par[2]))
      pickIdx = (tm >= bRange[1]) & (tm <= bRange[2])
      acpt = sum(pickIdx)
      if(acpt == 0){
        acpt = 0.5
        next
      }
      cnt = cnt + acpt
      res[(o_cnt + 1) : cnt] = tm[pickIdx]
      o_cnt = cnt
      acpt = acpt / length(pickIdx)
    }
    return(res[1:n])
  },
  gen_par = function(n){
    return(log(matrix(
      c(runif(n, 0.5,5), runif(n,0.5,5)),
      ncol = 2)
      ))
    },
  range = bRange,
  par.range = log(matrix(c(0.5, 5, 0.5, 5), 2, 2))
)

# further truncation of beta
tbRange = c(0.25, 0.75)
tBeta_family = list(pdf = function(x, par){dbeta(x, exp(par[1]), exp(par[2]))/
    (pbeta(tbRange[2], exp(par[1]), exp(par[2])) -
       pbeta(tbRange[1], exp(par[1]), exp(par[2])))},
    rpdf = function(n, par){
      cnt = 0
      o_cnt = 0
      acpt = 1
      res = rep(0, n)
      while(cnt < n){
        tm = rbeta(floor((n-cnt) / acpt)+1, exp(par[1]), exp(par[2]))
        pickIdx = (tm >= tbRange[1]) & (tm <= tbRange[2])
        acpt = sum(pickIdx)
        if(acpt == 0){
          acpt = 0.5
          next
        }
        cnt = cnt + acpt
        res[(o_cnt + 1) : cnt] = tm[pickIdx]
        o_cnt = cnt
        acpt = acpt / length(pickIdx)
      }
      return(res[1:n])
    },
    gen_par = function(n){return(log(matrix(c(runif(n, 0.5,5),
                                              runif(n,0.5,5)), ncol = 2)))},
    range = tbRange,
    par.range = log(matrix(c(0.5, 5, 0.5, 5), 2, 2))
)
# sketchDenFam(wrapDenFam(tBeta_family), n = 100)

# assignInNamespace('beta_family', beta_family, ns = "densityFPCA")
# assignInNamespace('bRange', bRange, ns = "densityFPCA")

# Cauchy family
cauchy_family = list(
  pdf = function(x, par) dcauchy(x, location = par[1], scale = (par[2])),
  rpdf = function(n, par) rcauchy(n, location = par[1], scale = (par[2])),
  gen_par = function(n) matrix(c(runif(n,-5,5), (runif(n, 0.5,3))), ncol=2),
  range = c(-10,10),
  par.range = matrix(c(-5, 5, 0.5, 3), 2, 2)
  )
# sketchDenFam(wrapDenFam(cauchy_family), n=100)

# Truncated Cauchy on c(0, 1)
tCrange = c(0,1)
tCauchy_family = list(
  pdf = function(x, par){
    tm_pdf = function(x){dcauchy(x, location = par[1], scale = exp(par[2]))}
    nrml_const = pcauchy(tCrange[2], location = par[1], scale = exp(par[2])) - pcauchy(tCrange[1], location = par[1], scale = exp(par[2]))
    #nrml_const = integrate(tm_pdf, lower = tCrange[1], upper = tCrange[2])$value
    tm_pdf(x)/nrml_const
  },
  rpdf = function(n, par){
    cnt = 0
    o_cnt = 0
    acpt = 1
    res = rep(0, n)
    while(cnt < n){
      tm = rcauchy(floor((n-cnt) / acpt)+1, par[1], scale = exp(par[2]))
      pickIdx = tm >= tCrange[1] & tm <= tCrange[2]
      acpt = sum(pickIdx)
      if(acpt == 0){
        acpt = 0.5
        next
      }
      cnt = cnt + acpt
      res[(o_cnt + 1) : cnt] = tm[pickIdx]
      o_cnt = cnt
      acpt = acpt / length(pickIdx)
    }
    return(res[1:n])
  },
  gen_par = function(n){matrix(c(runif(n,0,1), log(runif(n, 0.5,1))), ncol=2)},
  range = tCrange,
  par.range = matrix(c(0, 1, log(0.5), log(1)), 2, 2)
)

#' Bimodal family
#'
#' @return a list of
#' - `pdf`: density function;
#' - `rpdf`: random generator;
#' - `gen_par`: generate parameters;
#' - `par.range`: matrix of parameter range, one column for one parameter
#' - `range`: domain of support.
#'
#' @details The density function is proportional to
#' \deqn{
#'   \exp( (4 + \theta) x - (26.5 + \theta) x^2 + 47 x^3 - 25 x^4 )
#' }
#' on \eqn{x \in (0, 1)}.
#'
#' @export
#'
biM_family <- function(){

  # biMrange <- c(0,1)
  edge <- c(0, 1)

  unPdf <- function(x, par){
    nGrid = length(x)
    exp(
      rep(4+par, nGrid) * x -
      rep(26.5 + par, nGrid) * x^2 +
      rep(47, nGrid) * x^3
      - rep(25, nGrid) * x^4
    )
  }
  list(
    pdf = function(x, par){
      nrml_const = integrate(
        unPdf, lower = edge[1], upper = edge[2],
        par = par
      )$value
      return(unPdf(x, par) / nrml_const)
    },
    rpdf = function(n, par){
      # if(tolVar(function(x) unPdf(x, par), edge) < .Machine$double.eps^(1/3)){
      #   return(
      #     runif(n, min = min(edge), max = max(edge))
      #   )
      # }else{
      #   gen <- Runuran::pinv.new(
      #     pdf = unPdf,
      #     lb = min(edge), ub = max(edge),
      #     center = mean(edge),
      #     par = par
      #   )
      #   return(Runuran::ur(gen, n))
      # }
      # SimDesign::rejectionSampling(
      #   max(100, 2*n),
      #   df = function(x){
      #     nGrid = length(x)
      #     exp(rep(4+par, nGrid) * x - rep(26.5 + par, nGrid) * x^2 + rep(47, nGrid) * x^3 - rep(25, nGrid) * x^4)
      #   },
      #   dg = function(x) dunif(x, min = edge[1], max = edge[2]),
      #   rg = function(n) runif(n, min = edge[1], max = edge[2])
      # )[1:n]
      rejSampling(
        n, df = function(x){
          exp((4 + par) * x - (26.5 + par) * x^2 + 47 * x^3 - 25 * x^4)
        },
        dg = function(x) dunif(x, min = edge[1], max = edge[2]),
        rg = function(x) runif(x, min = edge[1], max = edge[2]),
        support = edge
      )
    },
    gen_par = function(n){matrix(runif(n, -5, 5), ncol = 1)},
    range = edge,
    par.range = matrix(c(-5, 5), ncol = 1)
  )
}
# den.fam <- biM_family()
# checkSampler(den.fam, size = 500, n = 100, what = 'diff')
# checkSampler(den.fam, size = 5000, n = 100, what = 'kde')
# checkSampler(den.fam, size = 50, n = 100, what = 'log')
# checkSampler(den.fam, size = 50, n = 100, what = 'pdf')
# checkDenFamGrident(den.fam)
# checkDenFamNumeric(den.fam)

xdBiM <- function(){
  # identical to biM_family, but different parameter space,
  # so as to be consistent with the writeup.
  res <- biM_family()
  res$gen_par <- function(n){matrix(runif(n, 0, 10), ncol = 1)}
  res$par.range <- matrix(c(0, 10), ncol = 1)
  return(res)
}
# den.fam <- xdBiM()
# checkSampler(den.fam, size = 500, n = 100, what = 'diff')
# checkSampler(den.fam, size = 500, n = 100, what = 'kde')
# checkSampler(den.fam, size = 50, n = 100, what = 'log')
# checkSampler(den.fam, size = 50, n = 100, what = 'pdf')
# checkDenFamGrident(den.fam)
# checkDenFamNumeric(den.fam)


## exp(polynomial)
expPoly <- function(edge = c(0, 1), deg = 3){
  # edge <- c(0, 1)
  unPdf <- function(x, par){
    mat.x <- sapply(seq_len(deg), function(i) x^i, simplify = FALSE)
    mat.x <- matrix(unlist(mat.x), ncol = deg)
    exp(
      colSums(par * t(mat.x))
    )
  }
  nrmlConst <- function(par){
    grid <- seq(from = min(edge), to = max(edge), length.out = 512)
    nGrid = length(grid)
    return(
      dotL2(
        rep(1, nGrid),
        unPdf(grid, par),
        grid = grid
      )
    )
  }
  return(list(
    pdf = function(x, par){
      return(unPdf(x, par) / nrmlConst(par))
    },
    rpdf = function(n, par){
      # it seems Runuran has some issue with flat pdf
      if(tolVar(function(x) unPdf(x, par), edge) < .Machine$double.eps^(1/3)){
        return(
          runif(n, min = min(edge), max = max(edge))
        )
      }else{
        gen <- Runuran::pinv.new(
          pdf = unPdf,
          lb = min(edge), ub = max(edge),
          center = mean(edge),
          par = par
        )
        return(Runuran::ur(gen, n))
      }
    },
    gen_par = function(n){
      matrix(runif(deg * n, -1, 1), ncol = deg)
    },
    range = edge,
    par.range = matrix(c(-1, 1), ncol = deg, nrow = 2)
  ))
}

## exp(sin)
expSin <- function(){
  edge <- c(0, 2 * pi)
  unPdf <- function(x, par) exp(sin(par[1] * x + par[2]))
  nrmlConst <- function(par){
    grid <- seq(from = min(edge), to = max(edge), length.out = 1000)
    n.grid <- length(grid)
    return(
      integrate(
        unPdf,
        lower = min(edge), upper = max(edge),
        par = par
      )$value
    )
  }
  return(list(
    pdf = function(x, par) unPdf(x, par) / nrmlConst(par),
    rpdf = function(n, par){
      # it seems Runuran has some issue with flat pdf
      if(tolVar(function(x) unPdf(x, par), edge) < .Machine$double.eps^(1/3)){
        return(
          runif(n, min = min(edge), max = max(edge))
        )
      }else{
        gen <- Runuran::pinv.new(
          pdf = unPdf,
          lb = min(edge), ub = max(edge),
          center = mean(edge),
          par = par
        )
        return(Runuran::ur(gen, n))
      }
    },
    gen_par = function(n){
      matrix(runif(2 * n, -3, 3), ncol = 2)
    },
    par.range = matrix(c(-3, 3), ncol = 2, nrow = 2),
    range = c(0, 2 * pi)
  ))
}


## mixture Gaussian
mixGauss_family = list(
  pdf = function(x, par){
    dnorm(x, par[1], exp(par[4]))/3 +
    dnorm(x, par[2], exp(par[5]))/3 +
    dnorm(x, par[3], exp(par[6]))/3},
  rpdf = function(n, par){
    idx = replicate(n, sample.int(3,1))
    rnorm(n, par[idx], exp(par[3+idx]))},
  gen_par = function(n){return(t(
    replicate(n,
              c(runif(3, -5,5), log(runif(3,1,5)))
              )
    ))},
  range = c(-20,20),
  par.range = matrix(c(
    -5, 5, -5, 5, -5, 5,
    log(1), log(5), log(1), log(5), log(1), log(5)
  ), nrow = 2)
)

## mixture Gaussian with truncation
tMixGauss_family = list(
  pdf = function(x, par){
    (dnorm(x, par[1], exp(par[4]))/3 +
      dnorm(x, par[2], exp(par[5]))/3 +
      dnorm(x, par[3], exp(par[6]))/3) /
      (pnorm(7, par[1], exp(par[4]))/3 +
         pnorm(7, par[2], exp(par[5]))/3 +
         pnorm(7, par[3], exp(par[6]))/3 -
       pnorm(-7, par[1], exp(par[4]))/3 -
         pnorm(-7, par[2], exp(par[5]))/3 -
         pnorm(-7, par[3], exp(par[6]))/3)},
  rpdf = function(n, par){
    cnt = 0
    o_cnt = 0
    acpt = 1
    res = rep(0, n)
    while(cnt < n){
      idx = replicate(floor((n-cnt) / acpt)+1, sample.int(3,1))
      tm = rnorm(floor((n-cnt) / acpt)+1, par[idx], exp(par[3+idx]))
      pickIdx = (tm >= -7) & (tm <= 7)
      acpt = sum(pickIdx)
      if(acpt == 0){
        acpt = 0.5
        next
      }
      cnt = cnt + acpt
      res[(o_cnt + 1) : cnt] = tm[pickIdx]
      o_cnt = cnt
      acpt = acpt / length(pickIdx)
    }
    return(res[1:n])
    },

  gen_par = function(n){return(t(
    replicate(n,
              c(runif(3, -5,5), log(runif(3,1,5)))
    )
  ))},
  range = c(-7,7),
  par.range = matrix(c(
    -5, 5, -5, 5, -5, 5,
    log(1), log(5), log(1), log(5), log(1), log(5)
  ), nrow = 2)
)

# functions returning truncated family ================================
# with truncation specified by argument

# truncated location-scale t distribution
tStudent <- function(edge){
  if(length(edge) == 1){
    edge = c(-1, 1) * abs(edge)
  }
  edge = sort(edge)

  nrmlizer <- function(par){
    extraDistr::plst(edge[2], df = par[1], mu = par[2], sigma = par[3]) -
      extraDistr::plst(edge[1], df = par[1], mu = par[2], sigma = par[3])
  }
  return(
    list(
      pdf = function(x, par){
        extraDistr::dlst(x, df = par[1], mu = par[2], sigma = par[3]) / nrmlizer(par)
      },
      rpdf = function(n, par) r.truncate(
        n = n, rpdf = extraDistr::rlst, edge = edge,
        df = par[1], mu = par[2], sigma = par[3]
      ),
      gen_par = function(n){
        matrix(c(
          runif(n, 0.1, 3),
          # rnorm(n, 0, 3),
          runif(n, -7, 7),
          runif(n, 0.5, 5)),
          ncol = 3)
      },
      range = edge,
      par.range = matrix(
        c(0.1, 3,  -7, 7,    0.5, 5), ncol = 3
      )
    )
  )
}
# # LEGACY
# tStudent <- function(edge){
#   if(length(edge) == 1){
#     edge = c(-1, 1) * abs(edge)
#   }
#   edge = sort(edge)
#
#   nrmlizer <- function(par){
#     pt(edge[2], df = par[1], ncp = par[2]) -
#       pt(edge[1], df = par[1], ncp = par[2])
#   }
#   return(
#     list(
#       pdf = function(x, par){
#         dt(x, df = par[1], ncp = par[2]) / nrmlizer(par)
#       },
#       rpdf = function(n, par) r.truncate(
#         n = n, rpdf = rt, edge = edge,
#         df = par[1], ncp = par[2]
#       ),
#       gen_par = function(n){
#         matrix(c(
#           runif(n, 1, 5),
#           runif(n, -5, 5)),
#         ncol = 2)
#       },
#       range = edge,
#       par.range = matrix(
#         c(1, 5,    -5, 5), ncol = 2
#       )
#     )
#   )
# }
# # END LEGACY

# truncated gamma distribution, left side set as 0.1 if not specified.
tGamma <- function(edge){
  if(length(edge) == 1){
    edge <- c(0.1, abs(edge))
  }
  edge <- sort(edge)
  if (edge[1] < 0) stop('Check input edge.')

  nrmlizer <- function(par){
    pgamma(edge[2], shape = par[1], rate = par[2]) -
      pgamma(edge[1], shape = par[1], rate = par[2])
  }
  return(
    list(
      pdf = function(x, par){
        dgamma(x, shape = par[1], rate = par[2]) / nrmlizer(par)
      },
      rpdf = function(n, par) r.truncate(
        n = n, rpdf = rgamma, edge = edge,
        shape = par[1], rate = par[2]
      ),
      gen_par = function(n) matrix(runif(2*n, 0.1, 5), ncol = 2),
      range = edge,
      par.range = matrix(
        c(0.1, 5, 0.1, 5), 2, 2
      )
    )
  )
}


#' Truncated Normal family
#'
#' @param edge domain of support, single value taken as `c(-1, 1) * abs(edge)`;
#' @param easy whether the parameters would be generated with fewer variation.
#'
#' @return a list of
#' - `pdf`: density function;
#' - `rpdf`: random generator;
#' - `gen_par`: generate parameters;
#' - `par.range`: matrix of parameter range, one column for one parameter
#' - `range`: domain of support.
#'
#' @export
#'
tNrml = function(edge, easy = FALSE){
  if(length(edge) == 1){
    edge = c(-1, 1) * abs(edge)
  }
  edge = sort(edge)
  if(!easy){
    return(
      list(
        pdf = function(x, par){
          dnorm(x, par[1], (par[2]))/(
          pnorm(edge[2], par[1], (par[2])) -
            pnorm(edge[1], par[1], (par[2])) ) },

        rpdf = function(n, par){
          cnt = 0
          o_cnt = 0
          acpt = 1
          res = rep(0, n)
          while(cnt < n){
            tm = rnorm(floor((n-cnt) / acpt)+1, par[1], (par[2]))
            pickIdx = (tm >= edge[1]) & (tm <= edge[2])
            acpt = sum(pickIdx)
            if(acpt == 0){
              acpt = 0.5
              next
            }
            cnt = cnt + acpt
            res[(o_cnt + 1) : cnt] = tm[pickIdx]
            o_cnt = cnt
            acpt = acpt / length(pickIdx)
          }
          return(res[1:n])
        },

        gen_par = function(n){return(matrix(c(runif(n, -5, 5),
                                              (runif(n, 1, 5))), ncol = 2))},
        range = edge,
        par.range = matrix(c(
          -5,     5,
          (1), (5)
        ), nrow = 2)
      )
    )
  }else{
    return(
      list(
        pdf = function(x, par){
          dnorm(x, par[1], (par[2]))/(
            pnorm(edge[2], par[1], (par[2])) -
              pnorm(edge[1], par[1], (par[2])) ) },

        rpdf = function(n, par){
          cnt = 0
          o_cnt = 0
          acpt = 1
          res = rep(0, n)
          while(cnt < n){
            tm = rnorm(floor((n-cnt) / acpt)+1, par[1], (par[2]))
            pickIdx = (tm >= edge[1]) & (tm <= edge[2])
            acpt = sum(pickIdx)
            if(acpt == 0){
              acpt = 0.5
              next
            }
            cnt = cnt + acpt
            res[(o_cnt + 1) : cnt] = tm[pickIdx]
            o_cnt = cnt
            acpt = acpt / length(pickIdx)
          }
          return(res[1:n])
        },

        gen_par = function(n){return(matrix(c(runif(n, -2, 2),
                                              (runif(n, 2, 4))), ncol = 2))},
        range = edge,
        par.range = matrix(c(
          -2,   2,
           2,   4
        ), nrow = 2)
      )
    )
  }
}

# sketchDenFam(tNrml(20), n = 100)
# sketchDenFam(tNrml(20), n = 100, scale = 'log')
#
# sketchDenFam(tNrml(10), n = 100)
# sketchDenFam(tNrml(10), n = 100, scale = 'log')
#
# sketchDenFam(tNrml(5), n = 100)
# sketchDenFam(tNrml(5), n = 100, scale = 'log')
#
# sketchDenFam(tNrml(c(0,1)), n = 100)
# sketchDenFam(tNrml(c(0,1)), n = 100, scale = 'log')

tCauchy = function(edge){
  if(length(edge) == 1){
    edge = c(-1, 1) * abs(edge)
  }
  edge = sort(edge)

  return(
    list(
      pdf = function(x, par){
        dcauchy(x, location = par[1], scale = (par[2])) /
          (pcauchy(edge[2], location = par[1], scale = (par[2])) -
             pcauchy(edge[1], location = par[1], scale = (par[2])))
      },

      rpdf = function(n, par){
        cnt = 0
        o_cnt = 0
        acpt = 1
        res = rep(0, n)
        while(cnt < n){
          tm = rcauchy(floor((n-cnt) / acpt)+1, par[1], scale = (par[2]))
          pickIdx = tm >= edge[1] & tm <= edge[2]
          acpt = sum(pickIdx)
          if(acpt == 0){
            acpt = 0.5
            next
          }
          cnt = cnt + acpt
          res[(o_cnt + 1) : cnt] = tm[pickIdx]
          o_cnt = cnt
          acpt = acpt / length(pickIdx)
        }
        return(res[1:n])
      },

      gen_par = function(n){matrix(c(runif(n,-5,5), (runif(n, 0.5, 2))), ncol=2)},
      range = edge,
      par.range = matrix(c(
        -5,  5,
        0.5, 2
      ), nrow = 2)
    )
  )
}

# sketchDenFam(tCauchy(100), n = 100)
# sketchDenFam(tCauchy(100), n = 100, scale = 'log')
# sketchDenFam(tCauchy(5), n = 100)
# sketchDenFam(tCauchy(5), n = 100, scale = 'log')
# sketchDenFam(tCauchy(2), n = 100)
# sketchDenFam(tCauchy(2), n = 100, scale = 'log')
# sketchDenFam(tCauchy(7), n = 100)
# sketchDenFam(tCauchy(7), n = 100, scale = 'log')


xdCauchy <- function(){
  # return the Cauchy family as that in the writeup
  # i.e., x in (0, 1), density f(x) \prop 1 / (1 + ((x - l) / s)^2)
  # and l ~ Unif(0, 1), s ~ Unif(0.1, 0.5).
  edge <- c(0, 1)
  return(
    list(
      pdf = function(x, par){
        dcauchy(x, location = par[1], scale = (par[2])) /
          (pcauchy(edge[2], location = par[1], scale = (par[2])) -
             pcauchy(edge[1], location = par[1], scale = (par[2])))
      },

      rpdf = function(n, par){
        r.truncate(
          n, rpdf = rcauchy, edge = edge,
          location = par[1], scale = par[2]
        )
      },

      gen_par = function(n){
        matrix(c(
          runif(n, 0, 1),
          runif(n, 0.1, 0.5)
        ), ncol=2)
      },
      range = edge,
      par.range = matrix(c(
        0,  1,
        0.1, 0.5
      ), nrow = 2)
    )
  )
}

# den.fam <- xdCauchy()
# checkSampler(den.fam, size = 500, n = 100, what = 'diff')
# checkSampler(den.fam, size = 50, n = 100, what = 'log')
# checkSampler(den.fam, size = 50, n = 100, what = 'pdf')
# checkDenFamGrident(den.fam)
# checkDenFamNumeric(den.fam)

#' Truncated Gaussian mixture family.
#'
#' @param edge domain of support, single value taken as `c(-1, 1) * abs(edge)`.
#' @param num.mixture number of mixture.
#'
#' @return a list of
#' - `pdf`: density function;
#' - `logpdf.gr`: gradient of loglikelihood;
#' - `rpdf`: random generator;
#' - `gen_par`: generate parameters;
#' - `par.range`: matrix of parameter range, one column for one parameter
#' - `range`: domain of support;
#' - `identify_par`: function to sort parameter (per mixture) for
#' identifiability.
#'
#' @export
#'
tMixGauss <- function(edge, num.mixture = 3){
  # construct a truncated Gaussian mixture family.

  # gaussian mixture N(mu[i], sigma[i]^2)
  # mixing weight alpha[i]
  # truncated according to feed edge
  # mu, sigma, theta(polar) are parameters, listed in par as
  # 1st k: mix.mean = mu
  # 2nd k: mix.sd = sigma
  # 3rd k(drop 1): theta
  # where theta s.t.
  #   lambda <- pol2cart(c(1, theta))
  #   alpha <- lambda / sum(lambda)

  if (num.mixture <= 0) {
    stop('check input num.mixture.')
  }
  if (length(edge) == 1) {
    edge = c(-1, 1) * abs(edge)
  }
  edge = sort(edge)

  if (num.mixture == 1) {
    return(
      tNrml(edge)
    )
  }

  pdf <- function(x, par){
    # vectorized density function

    # translating
    # num.mixture <- (length(par) + 1) / 3
    mix.mean <- par[seq_len(num.mixture)]
    mix.sd <- par[num.mixture + seq_len(num.mixture)]
    mix.polar <- par[2 * num.mixture + seq_len(num.mixture - 1)]
    mix.weight <- as.numeric(pol2cart(c(1, mix.polar)))
    mix.weight <- mix.weight / sum(mix.weight)

    mat.tm <- matrix(0, nrow = length(x), ncol = num.mixture)
    for (i in seq_len(num.mixture)) {
      mat.tm[, i] <- mix.weight[i] * dnorm(x, mix.mean[i], mix.sd[i])
    }

    nrml.const <-
      sum(mix.weight * pnorm(edge[2], mix.mean, mix.sd)) -
      sum(mix.weight * pnorm(edge[1], mix.mean, mix.sd))

    return(
      rowSums(mat.tm) / nrml.const
    )
  }

  logpdf.gr <- function(x, par){
    # vectorized gradient function of log-density

    # translating
    mix.mean <- par[seq_len(num.mixture)]
    mix.sd <- par[num.mixture + seq_len(num.mixture)]
    mix.polar <- par[2 * num.mixture + seq_len(num.mixture - 1)]
    mix.lambda <- as.numeric(pol2cart(c(1, mix.polar)))
    mix.weight <- mix.lambda / sum(mix.lambda)
    a <- min(edge)
    b <- max(edge)

    # derivatives of the part related to cdf (noted as nPhi)
    # where nPhi is a matrix, 1 row correspond to 1 observation in x
    nPhi <- pnorm(b, mix.mean, mix.sd) - pnorm(a, mix.mean, mix.sd)
    nPhi <- t(matrix(nPhi, nrow = num.mixture, ncol = length(x)))
    tm.a <- dnorm(a, mix.mean, mix.sd)
    tm.b <- dnorm(b, mix.mean, mix.sd)
    d.nPhi.mu <- (tm.a - tm.b)
    d.nPhi.mu <- t(matrix(d.nPhi.mu, nrow = num.mixture, ncol = length(x)))
    d.nPhi.sigma <- ((a - mix.mean) * tm.a - (b - mix.mean) * tm.b) / mix.sd
    d.nPhi.sigma <-
      t(matrix(d.nPhi.sigma, nrow = num.mixture, ncol = length(x)))

    # derivatives of the part related to pdf (noted as phi)
    # where phi is a matrix, 1 row correspond to 1 observation in x
    phi <- matrix(0, ncol = num.mixture, nrow = length(x))
    d.phi.mu <- matrix(0, ncol = num.mixture, nrow = length(x))
    d.phi.sigma <- matrix(0, ncol = num.mixture, nrow = length(x))
    for(i in seq_len(num.mixture)){
      phi[, i] <- dnorm(x, mix.mean[i], mix.sd[i])
      d.phi.mu[, i] <- (x - mix.mean[i]) / mix.sd[i]^2 * phi[, i]
      d.phi.sigma[, i] <-
        ((x - mix.mean[i])^2 - mix.sd[i]^2) / mix.sd[i]^3 * phi[, i]
    }

    # some temp vals
    tm.n1 <- colSums(mix.weight * t(phi))
    tm.n2 <- colSums(mix.weight * t(nPhi))
    mat.mix.weight <-
      t(matrix(mix.weight, nrow = num.mixture, ncol = length(x)))

    # gradient w.r.t. mix.weight
    d.l.weight <- phi / tm.n1 - nPhi / tm.n2

    # adjust due to reparameterization ==
    # jacobian of d lambda / d theta
    D.polar2cart <- jacob.pol2cart(c(1, mix.polar))
    D.polar2cart <- D.polar2cart[, -1, drop = FALSE]
    # jacobian of d mix.weight / d lambda
    D.lambda2weight <- matrix(-1 * mix.lambda, num.mixture, num.mixture)
    diag(D.lambda2weight) <- sum(mix.lambda) - mix.lambda
    D.lambda2weight <- D.lambda2weight / sum(mix.lambda) ^ 2
    # adjustment matrix
    D.adj <- D.lambda2weight %*% D.polar2cart
    # gradient w.r.t. actual parameters cooresponding to polar coord
    d.l.polar <- d.l.weight %*% D.adj

    # gradient w.r.t. mix.mean
    d.l.mu <- mat.mix.weight * d.phi.mu / tm.n1 -
      mat.mix.weight * d.nPhi.mu / tm.n2

    # gradient w.r.t. mix.sd
    d.l.sigma <- mat.mix.weight * d.phi.sigma / tm.n1 -
      mat.mix.weight * d.nPhi.sigma / tm.n2

    return(cbind(
      d.l.mu,
      d.l.sigma,
      d.l.polar
    ))

  }

  rpdf <- function(n, par){
    # random generating function

    # translating
    # num.mixture <- (length(par) + 1) / 3
    mix.mean <- par[seq_len(num.mixture)]
    mix.sd <- par[num.mixture + seq_len(num.mixture)]
    mix.polar <- par[2 * num.mixture + seq_len(num.mixture - 1)]
    mix.weight <- as.numeric(pol2cart(c(1, mix.polar)))
    mix.weight <- mix.weight / sum(mix.weight)

    cnt = 0
    o_cnt = 0
    acpt = 1
    res = rep(0, n)
    while(cnt < n){
      idx <- sample.int(
        num.mixture, size = floor((n-cnt) / acpt) + 1,
        prob = mix.weight, replace = TRUE
      )
      tm <- rnorm(floor((n-cnt) / acpt) + 1, mix.mean[idx], mix.sd[idx])
      pickIdx <- (tm >= edge[1]) & (tm <= edge[2])
      acpt <- sum(pickIdx)
      if(acpt == 0){
        acpt <- 0.5
        next
      }
      cnt <- cnt + acpt
      res[(o_cnt + 1) : cnt] <- tm[pickIdx]
      o_cnt <- cnt
      acpt <- acpt / length(pickIdx)
    }
    return(res[1:n])
  }
  # range of parameters
  mix.mean.range <- c(-5, 5)
  mix.sd.range <- c(0.5, 5)
  # the range matrix
  par.range = matrix(
    c(
    rep(mix.mean.range, num.mixture),
    rep(mix.sd.range, num.mixture),
    rep(c(0, pi/2), num.mixture)
    ),
    nrow = 2
  )
  par.range <- par.range[, -ncol(par.range), drop = FALSE]

  gen_par = function(n){
    # parameter generation

    # generate mixture probability from Dirichlet distribution
    mat.weight <- extraDistr::rdirichlet(n, rep(1, times = num.mixture))
    # transform into polar coord
    mat.polar <- (cart2pol(mat.weight))[, -1]

    mat.par = matrix(0, ncol = 3 * num.mixture - 1, nrow = n)
    # the mixture mean
    mat.par[, seq_len(num.mixture)] <-
      runif(num.mixture * n, min(mix.mean.range), max(mix.mean.range))
    # the mixture sd
    mat.par[, num.mixture + seq_len(num.mixture)] <-
      runif(num.mixture * n, min(mix.sd.range), max(mix.sd.range))
    # the mixture weight (will not be auto exclude if num.mixture = 1)
    mat.par[, 2 * num.mixture + seq_len(num.mixture - 1)] <- mat.polar

    return(mat.par)
  }

  identify_par <- function(par){
    # make parameter identifiable, namely sorting the mixture by their means.
    # par could be a matrix where one row is one parameter

    if (num.mixture == 1){
      return(par)
    }

    if (!is.matrix(par)){
      par <- matrix(par, nrow = 1)
    }

    # translating
    # num.mixture <- (length(par) + 1) / 3
    mix.mean <- par[, seq_len(num.mixture), drop = FALSE]
    mix.sd <- par[, num.mixture + seq_len(num.mixture), drop = FALSE]

    mix.polar <- matrix(1, ncol = num.mixture, nrow = nrow(par))
    mix.polar[, -1] <-
      par[, 2 * num.mixture + seq_len(num.mixture - 1), drop = FALSE]
    mix.weight <- pol2cart(mix.polar)
    mix.weight <- mix.weight / rowSums(mix.weight)

    # sorting
    for(i in seq_len(nrow(mix.mean))){
      idx <- order(mix.mean[i, ])
      mix.mean[i, ] <- mix.mean[i, idx]
      mix.sd[i, ] <- mix.sd[i, idx]
      mix.weight[i, ] <- mix.weight[i, idx]
    }

    mix.polar <- (cart2pol(mix.weight))[, -1, drop = FALSE]

    # returning
    mat.par = matrix(0, ncol = 3 * num.mixture - 1, nrow = nrow(par))
    # the mixture mean
    mat.par[, seq_len(num.mixture)] <- mix.mean
    # the mixture sd
    mat.par[, num.mixture + seq_len(num.mixture)] <- mix.sd
    # the mixture weight (will not be auto exclude if num.mixture = 1)
    mat.par[, 2 * num.mixture + seq_len(num.mixture - 1)] <- mix.polar

    return(mat.par)
  }

  return(list(
    pdf = pdf,
    logpdf.gr = logpdf.gr,
    rpdf = rpdf,
    gen_par = gen_par,
    par.range = par.range,
    range = edge,
    identify_par = identify_par,
    convex = FALSE
  ))
}

#' Translate a list of weight, mean, var to/from working parameterization
#'
#' @param ls.par  a named list with `weight`, `mean`, and `var` for Gaussian
#' mixture family.
#' @param work.par an array of working parameters, as used in `tMixGauss()`.
#'
#' @return a array of working parameter if supplying `ls.par`, or a list of
#' `weight`, `mean`, and `var` if `work.par` supplied. Supplying both lead to
#' error.
#' @export
#'
tMixGauss.repar <- function(ls.par, work.par){
  # transplate a list of weigth, mean, var into/from working parameterization
  # args:
  # ls.par = a list with weight, mean, var
  # work.par =  an array of working parameters, formated as tMixGauss() family
  # returns:
  # a array of working parameter or a list of weight, mean, var
  # details:
  # by default (if only one argument) will do ls.par -> w.par,
  # say work.par <- tMixGauss.repar(whatever_a_list)
  # or ls.par <- tMixGauss.repar(work.par), translating to list.
  # will raise error if both supplied.
  if(missing(work.par)){
    polar <- cart2pol(ls.par$weight)
    num.k <- length(ls.par$weight)
    work.par <- rep(0, 3 * num.k - 1)
    work.par[seq_len(num.k)] <- ls.par$mean
    work.par[seq_len(num.k) + num.k] <- sqrt(ls.par$var)
    work.par[seq_len(num.k - 1) + 2 * num.k] <- polar[-1]
    return(work.par)
  }
  if(missing(ls.par)){
    num.k <- (length(work.par) + 1) / 3
    ls.par <- list()
    ls.par$mean <- work.par[seq_len(num.k)]
    ls.par$var <- (work.par[seq_len(num.k) + num.k]) ^ 2
    mix.polar <- work.par[2 * num.k + seq_len(num.k - 1)]
    mix.weight <- as.numeric(pol2cart(c(1, mix.polar)))
    mix.weight <- mix.weight / sum(mix.weight)
    ls.par$weight <- mix.weight
    return(ls.par)
  }
  stop('too much input')
}

# LEGACY
# tMixGauss <- function(edge, num.mixture = 3){
#   # gaussian mixture N(mu[i], sigma[i]^2)
#   # mixing weight p[i]
#   # truncated according to feed edge
#   # mu, sigma, p are parameters, listed as
#   # p[i] = par[(i-1) * k + 1]
#   # mu[i] = par[(i-1) * k + 2]
#   # sigma[i] = exp(par[(i-1) * k + 3])
#   # where k = num.mixture
#
#   if (num.mixture <= 0) {
#     stop('check input num.mixture.')
#   }
#   if (length(edge) == 1) {
#     edge = c(-1, 1) * abs(edge)
#   }
#   edge = sort(edge)
#
#   return(list(
#     pdf = function(x, par){
#       # make sure positive weight sum = 1
#       mix.weight <-
#         par[seq_len(length(par)) %% 3 == 1]^2
#       mix.weight <- mix.weight / sum(mix.weight)
#
#       mix.mean <- par[seq_len(length(par)) %% 3 == 2]
#       mix.sd <- exp(par[seq_len(length(par)) %% 3 == 0])
#
#       mat.tm <- matrix(0, ncol = length(x), nrow = num.mixture)
#       for (i in seq_len(num.mixture)) {
#         mat.tm[i, ] <- mix.weight[i] * dnorm(x, mix.mean[i], mix.sd[i])
#       }
#
#       nrml.const <-
#         sum(mix.weight * pnorm(edge[2], mix.mean, mix.sd)) -
#         sum(mix.weight * pnorm(edge[1], mix.mean, mix.sd))
#
#       return(
#         colSums(mat.tm) / nrml.const
#       )
#
#       # # make sure positive weight sum = 1
#       # par[c(1,4,7)] = par[c(1,4,7)]^2 / sum(par[c(1,4,7)]^2)
#       # (
#       #   dnorm(x, par[2], exp(par[3]))*par[1] +
#       #     dnorm(x, par[5], exp(par[6]))*par[4] +
#       #     dnorm(x, par[8], exp(par[9]))*par[7]
#       # )/(
#       #   pnorm(edge[2], par[2], exp(par[3]))*par[1] +
#       #     pnorm(edge[2], par[5], exp(par[6]))*par[4] +
#       #     pnorm(edge[2], par[8], exp(par[9]))*par[7] -
#       #     pnorm(edge[1], par[2], exp(par[3]))*par[1] -
#       #     pnorm(edge[1], par[5], exp(par[6]))*par[4] -
#       #     pnorm(edge[1], par[8], exp(par[9]))*par[7]
#       # )
#     },
#
#     rpdf = function(n, par){
#       # make sure positive weight sum = 1
#       mix.weight <-
#         par[seq_len(length(par)) %% 3 == 1]^2
#       mix.weight <- mix.weight / sum(mix.weight)
#
#       mix.mean <- par[seq_len(length(par)) %% 3 == 2]
#       mix.sd <- exp(par[seq_len(length(par)) %% 3 == 0])
#
#       cnt = 0
#       o_cnt = 0
#       acpt = 1
#       res = rep(0, n)
#       while(cnt < n){
#         idx <- sample.int(
#           num.mixture, size = floor((n-cnt) / acpt) + 1,
#           prob = mix.weight, replace = TRUE
#         )
#         tm <- rnorm(floor((n-cnt) / acpt) + 1, mix.mean[idx], mix.sd[idx])
#         pickIdx <- (tm >= edge[1]) & (tm <= edge[2])
#         acpt <- sum(pickIdx)
#         if(acpt == 0){
#           acpt <- 0.5
#           next
#         }
#         cnt <- cnt + acpt
#         res[(o_cnt + 1) : cnt] <- tm[pickIdx]
#         o_cnt <- cnt
#         acpt <- acpt / length(pickIdx)
#       }
#       return(res[1:n])
#     },
#
#     gen_par = function(n){
#       # generate mixture probability from Dirichlet distribution
#       mat.weight <- extraDistr::rdirichlet(n, rep(1, times = num.mixture))
#
#       # tm_weight = matrix(runif(2*n, 0, 1), ncol=2)
#       # mat_weight = matrix(0, ncol = 3, nrow = n)
#       # mat_weight[,2] = pmax(tm_weight[,1], tm_weight[,2])
#       # mat_weight[,3] = 1 - mat_weight[,2]
#       # mat_weight[,1] = pmin(tm_weight[,1], tm_weight[,2])
#       # mat_weight[,2] = mat_weight[,2] - mat_weight[,1]
#       # mat_weight = sqrt(mat_weight)
#
#       mat.par = matrix(0, ncol = 3 * num.mixture, nrow = n)
#       mat.par[, seq_len(ncol(mat.par)) %% 3 == 1] <-
#         sqrt(mat.weight)
#       mat.par[, seq_len(ncol(mat.par)) %% 3 == 2] <-
#         runif(num.mixture * n, -5,5)
#       mat.par[, seq_len(ncol(mat.par)) %% 3 == 0] <-
#         log(runif(num.mixture * n, 0.5, 5))
#       return(mat.par)
#     },
#     range = edge,
#     par.range = matrix(c(
#       0,        1,
#       -5,       5,
#       log(0.5), log(5),
#       0,        1,
#       -5,       5,
#       log(0.5), log(5),
#       0,        1,
#       -5,       5,
#       log(0.5), log(5)
#     ), nrow = 2)
#   ))
# }
#
# sketchDenFam(tMixGauss(100), n = 200)
# sketchDenFam(tMixGauss(100), n = 200, scale = 'log')
# sketchDenFam(tMixGauss(5), n = 200)
# sketchDenFam(tMixGauss(5), n = 200, scale = 'log')
# sketchDenFam(tMixGauss(2), n = 200)
# sketchDenFam(tMixGauss(2), n = 200, scale = 'log')
# sketchDenFam(tMixGauss(7), n = 200)
# sketchDenFam(tMixGauss(7), n = 200, scale = 'log')


# TBD: hierarchical model
# hierarc = list(
#   # a hierarchical model
#
#   # N ~ Poisson(par[1]) + 1
#   # mu[i] ~ Unif(-5,5)
#   # exp(sigma[i]) ~ Unif(1, 5)
#
#   # weight | N ~ Dirichlet(1/N, ..., 1/N)
#   # mixIdx | weight ~ multinomial(weight)
#   # X | N, mixIdx = i, mu, sigma ~(i.i.d.)~ N(mu[i], sigma^2[i])
#
#   pdf = function(x, par){}
#   rpdf = function(n, par){}
#   gen_par = function(n){
#
#   }
#
# )

# sketchDenFam(wrapDenFam(nrml_family), n=100)
# sketchDenFam(wrapDenFam(tNrml_family), n=100)
#
# sketchDenFam(wrapDenFam(beta_family), n=100, ylim = c(0,5))
#
# sketchDenFam(wrapDenFam(cauchy_family), n=100)
#
# sketchDenFam(wrapDenFam(tCauchy_family), n=100)
#
# sketchDenFam(wrapDenFam(biM_family), n=100)
#
# sketchDenFam(wrapDenFam(mixGauss_family), n=250)
# sketchDenFam(wrapDenFam(tMixGauss_family), n=250)
