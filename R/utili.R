# some useful =================================================================

rejSampling <- function(n, df, dg, rg, support){
  # rejection sampling
  # n: number of samples to generate
  # df: pdf (function) of target, vectorized
  # dg: pdf (function) of generator, vectorized
  # rg: generator
  # support: common support

  logll <- function(x) log(df(x)) - log(dg(x))
  # maximum
  log.bound.m <- optimize(logll, interval = support, maximum = TRUE)$maximum

  res <- numeric(0)
  n.batch <- max(100, 2 * n)
  while(length(res) < n){
    idx.unif <- runif(n.batch, 0, 1)
    propose <- rg(n.batch)
    res <- c(res, propose[log(idx.unif) < logll(propose) - log.bound.m])
  }
  return(res[seq(n)])
}

orthLog <- function(f, grid, against = 1){
  # compute log(f) - < logf, against >_{L^2},
  # args:
  #  f: a vector for f values.
  #  grid: the grid on which f is evaluated.
  #  against: a vector of same size, or recycled/truncated to the same size.
  # return:
  #  a vector of same length as f, as evaluated on grid.

  # sanity check.
  if(any(!is.finite(f)) | any(!is.finite(grid)) | any(!is.finite(against)))
    stop('nonfinite input.')
  if(any(!is.finite(log(f))))
    stop('nonfinite logf.')
  n.grid <- length(grid)
  if(length(f) != n.grid)
    stop('check input length.')
  if(length(against) != n.grid)
    against <- array(against, dim = n.grid)

  return(log(f) - cpp_dotL2(log(f), against, grid))
}

mat2list <- function(mat, by = 1){
  # convert either rows (by = 1) or columns (by = 2) of matrix to list of arraies.
  # args:
  # mat: matrix to convert,
  #      numeric array will be convert to matrix with as.matrix.
  # by: 1 then row, 2 then column.
  # return:
  # a list of arraies.
  mat <- as.matrix(mat)
  if(by == 2){
    return(
      lapply(seq_len(ncol(mat)), function(i) mat[, i])
    )
  }else{
    return(
      lapply(seq_len(nrow(mat)), function(i) mat[i, ])
    )
  }
}
# microbenchmark::microbenchmark(
#   list = alist(
#     l.d = as.list(as.data.frame(den.fam$par.range)),
#     m2l = mat2list(den.fam$par.range, by = 2)
#   )
# )


df2list <- function(df, y.name, by.name){
  # convert data.frame to list of arraies.
  # args:
  # df: data frame to convert,
  # y.name: which variable to use as response
  # by.name: which variable to use as grouping variable
  #          note: NA in by.name will be dropped!
  # return:
  # a list of arraies.
  #DEBUG
  # df <- data.frame(a = c(1,3,4,2,6, 8, 9), idx = c(1,2,1,1,2, NA, NaN))
  # y.name <- 'a'
  # by.name <- 'idx'

  if(!all(is.element(c(y.name, by.name), names(df)))){
    stop('check input: y.name and by.name')
  }
  idx <- factor(df[[by.name]])
  return(split(df[[y.name]], idx))
}

normalize <- function(mat.pdf, grid){
  # normalize a matrix of pdfs, one row treated as one curve.
  # grid: on which the densities are evaluated.
  grid.size <- length(grid)
  if(grid.size != ncol(mat.pdf)) stop('Check input dimension.')

  mat.res <- mat.pdf
  for(i in seq_len(nrow(mat.res))){
    mat.res[i, ] <- mat.res[i, ] /
      dotL2(mat.pdf[i, ], rep(1, grid.size), grid = grid)
  }
  return(mat.res)
}

relErr <- function(true, esti, eps = .Machine$double.eps^(1/3)){
  # compute relative error, true and esti need to be in same dimension.
  err <- true - esti
  # abs.err <- abs(err)
  rel.err <- err / abs(true)
  true[is.na(true)] <- 0
  rel.err[abs(true) <= eps] <- err[abs(true) <= eps]
  return(rel.err)
}

qgChisq <- function(p, mu, mat.var, mat.A){
  # computing the quantile of generalized Chi-square ditribution
  # args:
  # p: probability;
  # mu, mat.var: mean and covariance matrix;
  # mat.A: matrix of quadratic form.
  # returns:
  # quantiles.
  # details:
  # generalized Chi-square = X^T A X, with X \sim N(mu, mat.var)

  smpl <- rgChisq(10000, mu, mat.var, mat.A)
  return(
    quantile(smpl, probs = p)
  )
}
rgChisq <- function(n, mu, mat.var, mat.A){
  # computing the quantile of generalized Chi-square ditribution
  # args:
  # n: number of observations to draw;
  # mu, mat.var: mean and covariance matrix;
  # mat.A: matrix of quadratic form.
  # returns:
  # sample.
  # details:
  # generalized Chi-square = X^T A X, with X \sim N(mu, mat.var)

  d <- length(mu)
  x <- matrix(rnorm(n * d), nrow = d)
  mat.sd <- pracma::sqrtm(mat.var)$B
  x <- mat.sd %*% x + mu
  smpl <- as.numeric(colSums((t(mat.A) %*% x) * x))

  return(smpl)
}

ksTest <- function(obsv, density, grid, ...){
  # wrapper for Kolmogorov-Smirnov test.
  # args:
  # obsv = an array of observation;
  # density = either a pdf function or its value on grid;
  # grid = necessary when density is an array;
  # ...: additional arguments passed to ks.test.
  # returns:
  # what ks.test return.
  if (class(density) == 'function') {
    w.cdf <- function(x){
      integrate(
        density,
        lower = -Inf,
        upper = x
      )
    }
  } else {
    if (missing(grid)) stop('Need grid.')
    if (!is.numeric(density)) stop('Check input density.')
    probs <- cumsum(density[-length(density)] * diff(grid))
    probs <- c(0, probs) / max(probs)
    w.cdf = approxfun(x = grid, y = probs, yleft = 0, yright = 1)
  }
  return(
    ks.test(
      x = obsv,
      y = "w.cdf",
      ...
    )
  )
}

rowProd <- function(mat){
  # row product
  return(
    exp(rowSums(log(mat)))
  )
}

quadForm <- function(M, x){
  # compute t(x) %*% M %*% x
  # args:
  # M = a symmetric matrix
  # x = an array
  # returns:
  # a number.
  # Remark: c.f. emulator::quad.form().
  return(
    crossprod(crossprod(M, x), x)
  )
}

checkDerivatives <- function(x, fn, gr, check_derivatives_tol = 1e-04,
   ...){
  # does what nloptr::check.derivatives does, quieter and correct.
  # args: same as nloptr::check.derivatives.
  # returns: same as nloptr::check.derivatives.
  # different: better precision and better handling finite.diff.grad = 0.
    analytic.grad <- gr(x, ...)
    if (length(fn(x)) == 1) {
      finite.diff.grad <- nloptr::nl.grad(x, fn, ...)
    } else {
      finite.diff.grad <- nloptr::nl.jacobian(x, fn, ...)
    }
    relative.error <- ifelse(
      abs(finite.diff.grad) <= .Machine$double.eps^(1/3),
      analytic.grad,
      abs((analytic.grad - finite.diff.grad)/finite.diff.grad)
    )
    flag.derivative.warning <- relative.error > check_derivatives_tol
    return(
      list(
        analytic = analytic.grad,
        finite_difference = finite.diff.grad,
        relative_error = relative.error,
        flag_derivative_warning = flag.derivative.warning
      )
    )
}

# polar <-> cartesian coord ===================================================

pol2cart <- function(polar){
  # transform polar coord -> cartesian in R^n, n >= 2.
  # args:
  # polar = a matrix whose row is (r, theta1, ...)
  # returns:
  # a matrix whose row is (x1, x2, ...)

  if (!is.matrix(polar)) {
    polar <- matrix(polar, nrow = 1)
  }

  r <- polar[, 1]
  theta <- polar[, -1, drop = FALSE]
  n.dim <- ncol(theta) + 1
  n.pts <- length(r)
  res <- matrix(0, ncol = n.dim, nrow = n.pts)
  tm <- r
  for (i in seq_len(n.dim - 1)) {
    res[, i] <- tm * cos(theta[, i])
    tm <- tm * sin(theta[, i])
  }
  res[, n.dim] <- tm
  return(res)
}

jacob.pol2cart <- function(polar){
  # the Jacobian matrix of coordinate transform: polar -> cartesian
  # args:
  #  polar: a array of (r, theta_1, ..., theta_{n-1})
  # returns:
  #  a n X n matrix where the row = i, col = j element being \partial cartesian_i / \partial polar_j
  n <- length(polar)
  theta <- polar[-1, drop = FALSE]
  r <- polar[1]

  # matrix whose row products are cartesian coord, each col -> polar coord
  mat.coord <- matrix(1, n, n)
  # matrix of terms d(cartesian coord)/d(polar coord)
  mat.d <- matrix(0, n, n)
  mat.coord[, 1] <- r
  mat.d[, 1] <- 1
  for (j in seq_len(n - 1)) {
    mat.coord[j, j + 1] <- cos(theta[j])
    mat.coord[(j + 1) : (n), j + 1] <- sin(theta[j])

    mat.d[j, j + 1] <- -1 * sin(theta[j])
    mat.d[(j + 1) : (n), j + 1] <- cos(theta[j])
  }

  res <- matrix(0, ncol = n, nrow = n)
  for (j in seq_len(n)) {
    res[, j] <- mat.d[, j, drop = FALSE] * rowProd(mat.coord[, -j, drop = FALSE])
  }

  return(
    res
  )

}

cart2pol <- function(cart){
  # inverse of pol2cart
  # args:
  # cart = a matrix of coords, one row <-> one point. at least in R^2
  # returns:
  # a matrix of polar coords.

  # implementing by noting:
  # \sum(x_i^2) = r^2
  # \sum(x_{-1}^2) = r^2 * sin^2 theta_1
  # x_1^2 = r^2 cos^2 theta_1
  if (!is.matrix(cart)) {
    cart <- matrix(cart, nrow = 1)
  }

  r <- sqrt(rowSums(cart ^ 2))

  n.dim <- ncol(cart)
  n.pts <- nrow(cart)
  theta <- matrix(0, ncol = n.dim - 1, nrow = n.pts)

  theta[, 1] <- atan( sqrt( rowSums(cart[, -1, drop = FALSE]^2) / cart[, 1]^2) )
  tm <- r * sin(theta[, 1])
  for(i in seq_len(n.dim - 2) + 1){
    theta[, i] <- acos( cart[, i] / tm )
    tm <- tm * sin(theta[, i])
  }
  return(cbind(
    r,
    theta
  ))
}

# functions computing error between density estimates =========================

calErr = function(f, hatf, grid, grid2, method = c('ISE', 'KL', 'rev.KL', 'FisherRao', 'Wass..5', 'Wass.1', 'Wass.2')){
# grid = grid of f and hatf unles
# grid2 is not missing, in that case
# grid = grid of f, grid2 = grid of hatf

  # same my time if any NA in f or hatf or grid or grid2
  if(anyNA(f) | !all(is.finite(f)) | any(is.null(f)) |
     anyNA(hatf) | !all(is.finite(hatf)) | any(is.null(hatf)) |
     anyNA(grid) | !all(is.finite(grid)) | any(is.null(grid))){
    return(rep(NA, length(method)))
  }
  if(!missing(grid2)){
    if(anyNA(grid2) | !all(is.finite(grid2)) | any(is.null(grid2))){
      return(rep(NA, length(method)))
  }}

# too lazy to read through base::substitute and all the expression stuff
# just use the most intuitive way

  ls.err <- list()

  if(!missing(grid2)){
    if(any('ISE' == method)){
      ls.err[['ISE']] = tryCatch(
        calISE(f, hatf, grid = grid, grid2 = grid2),
        error = function(e){return(NA)}
      )
    }
    if(any('KL' == method)){
      ls.err$KL = tryCatch(
        calKL(f, hatf, grid = grid, grid2 = grid2),
        error = function(e){return(NA)}
      )
    }
    if(any('rev.KL' == method)){
      ls.err$rev.KL = tryCatch(
        calKL(hatf, f, grid = grid2, grid2 = grid),
        error = function(e){return(NA)}
      )
    }
    if(any('FisherRao' == method)){
      ls.err$FisherRao = tryCatch(
        calFR(f, hatf, grid = grid, grid2 = grid2),
        error = function(e){return(NA)}
      )
    }
    if(any('Wass..5' == method)){
      ls.err$Wass..5 = tryCatch(
        distWass(f, hatf, p = 0.5, grid = grid, grid2 = grid2),
        error = function(e){return(NA)}
      )
    }
    if(any('Wass.1' == method)){
      ls.err$Wass.1 = tryCatch(
        distWass(f, hatf, p = 1, grid = grid, grid2 = grid2),
        error = function(e){return(NA)}
      )
    }
    if(any('Wass.2' == method)){
      ls.err$Wass.2 = tryCatch(
        distWass(f, hatf, p = 2, grid = grid, grid2 = grid2),
        error = function(e){return(NA)}
      )
    }
  }else{
    # well, if no grid2
    if(any('ISE' == method)){
      ls.err$ISE = tryCatch(
        calISE(f, hatf, grid = grid),
        error = function(e){return(NA)}
      )
    }
    if(any('KL' == method)){
      ls.err$KL = tryCatch(
        calKL(f, hatf, grid = grid),
        error = function(e){return(NA)}
      )
    }
    if(any('rev.KL' == method)){
      ls.err$rev.KL = tryCatch(
        calKL(hatf, f, grid = grid),
        error = function(e){return(NA)}
      )
    }
    if(any('FisherRao' == method)){
      ls.err$FisherRao = tryCatch(
        calFR(f, hatf, grid = grid),
        error = function(e){return(NA)}
      )
    }
    if(any('Wass..5' == method)){
      ls.err$Wass..5 = tryCatch(
        distWass(f, hatf, p = 0.5, grid = grid),
        error = function(e){return(NA)}
      )
    }
    if(any('Wass.1' == method)){
      ls.err$Wass.1 = tryCatch(
        distWass(f, hatf, p = 1, grid = grid),
        error = function(e){return(NA)}
      )
    }
    if(any('Wass.2' == method)){
      ls.err$Wass.2 = tryCatch(
        distWass(f, hatf, p = 2, grid = grid),
        error = function(e){return(NA)}
      )
    }
  }
  return(setNames(as.numeric(unlist(ls.err[method])), method))
}

calFR = function(f1, f2, range, grid, grid2){
  # calculate Fisher-Rao metric defined as
  # acos( < sqrt(f1), sqrt(f2) >_{L2} )

  if(missing(grid) & missing(grid2)){
    w.f1 = function(x){sqrt(f1(x))}
    w.f2 = function(x){sqrt(f2(x))}
    return(acos(
      dotL2(w.f1, w.f2, range = range) /
        normL2(w.f1, range = range) / normL2(w.f2, range = range)
      ))
  }

  if(!missing(grid2)){
    cmGrid = sort(unique(c(grid, grid2)))
    w.f1 = sqrt(approx(x = grid, y = f1, xout = cmGrid, yleft = 0, yright = 0)$y)
    w.f2 = sqrt(approx(x = grid2, y = f2, xout = cmGrid, yleft = 0, yright = 0)$y)
    # browser()
  }else{
    cmGrid = grid
    w.f1 = sqrt(f1)
    w.f2 = sqrt(f2)
  }
  #browser()
  if(anyNA(w.f1) | !all(is.finite(w.f1)) | any(is.null(w.f1)) |
     anyNA(w.f2) | !all(is.finite(w.f2)) | any(is.null(w.f2)) ){
    stop('non-finite function value')
  }
  return(acos(cpp_dotL2(w.f1, w.f2, cmGrid) /
                sqrt(cpp_dotL2(w.f1, w.f1, cmGrid) * cpp_dotL2(w.f2, w.f2, cmGrid))))
}

calISE = function(f1, f2, range, grid, grid2){
  # calculating ISE of f1 and f2
  # either feeding functions and range
  # or feeding array of func values and the grids
  # if only grid input, will be used for both array.
  # if 2 grids, linear interpolate f1 & f2 on their union.
  if(missing(grid) & missing(grid2)){
    return(normL2(function(x){f1(x) - f2(x)}, range = range)^2)
  }
  if(!missing(grid2)){
    cmGrid = sort(unique(c(grid, grid2)))
    f1 = approx(x = grid, y = f1, xout = cmGrid, yleft = 0, yright = 0)$y
    f2 = approx(x = grid2, y = f2, xout = cmGrid, yleft = 0, yright = 0)$y
    # err = w.f1 - w.f2
    # browser()
  }else{
    cmGrid = grid
    # err = f1 - f2
  }

  # normalize f1 & f2
  # w.f1 = f1 / dotL2(f1, rep(1, length(cmGrid)), grid = cmGrid)
  # w.f2 = f2 / dotL2(f2, rep(1, length(cmGrid)), grid = cmGrid)

  # error
  err = f1 - f2

  #browser()
  if(anyNA(err) | !all(is.finite(err)) | any(is.null(err))){
    stop('non-finite function value')
  }
  return(cpp_dotL2(err,err, cmGrid))
}


calKL = function(f1, f2, range, grid, grid2){
  # calculate D_{KL}(f1 || f2) = E_{X~f1} log(f1(X)/f2(X))
  # range = range on which f1 & f2 sit
  # grid, grid2 = grid points on which f1 &(or) f2 are evaluated,
  # if f1 & f2 are array instead of functions.
  # NOTE: if passing through functions, the user must make sure
  # f1 << f2, otherwise may run without error but
  # give a wrong answer.

  if(missing(grid) & missing(grid2)){
    lg = function(x){log(f1(x)/f2(x))}
    return(dotL2(f1, lg, range = range))
  }

  if(!missing(grid2)){
    cmGrid = sort(unique(c(grid, grid2)))
    f1 = approx(x = grid, y = f1, xout = cmGrid, yleft = 0, yright = 0)$y
    f2 = approx(x = grid2, y = f2, xout = cmGrid, yleft = 0, yright = 0)$y
    # browser()
  }else{
    cmGrid = grid
  }

  # normalize f1 & f2 to avoid negative D_KL
  nrml_const1 = dotL2(f1, rep(1, length(cmGrid)), grid = cmGrid)
  nrml_const2 = dotL2(f2, rep(1, length(cmGrid)), grid = cmGrid)
  if(nrml_const1 != 0) f1 = f1 / nrml_const1
  if(nrml_const2 != 0) f2 = f2 / nrml_const2

  # check for absolute continuity, i.e. f1 << f2
  # if there is any f2 = 0 but f1 != 0, return NA
  # and some times due to machine precision
  # dropping off some, say 1e-5 / 1e-315 gives Inf
  lg = f1 / f2
  if(any(f2 == 0 & f1 != 0) | any(!is.finite(lg[f1!=0]))){
    warning('No absolute continuity.')
    return(NA)
  }
  lg = log(lg)
  lg[f1 == 0] = 0

  #browser()
  if(anyNA(lg) | !all(is.finite(lg)) | any(is.null(lg))){
    warning('non-finite logarithm value')
    return(NaN)
  }
  return(cpp_dotL2(f1,lg, cmGrid))
}

# optimization with nloptr ====================================================

safeOptim <- function(
  x0,
  fn, gr = NULL,
  lower, upper,
  method = "MLSL",
  backup.method = "DIRECT",
  nl.control = list(), ...
){
  # a safe wrapper for NLOptim.
  # args:
  # x0: initial point, missing -> midpoint of the box, omitted for some methods;
  # fn, gr: objective function and its gradient (optional);
  # lower, upper: bounds;
  # method: method to use, say DIRECT, DIRECT-P, LBFGS, MLSL, Nelder-Mead, SubPlex;
  # backup.method: what to use if method fails;
  # nl.control: control list to pass to nloptr functions, c.f. nloptr::nl.opts;
  # ...: additional arguments passed to fn and gr.
  # returns:
  # what nloptr returns, i.e., a list of par, value, iter, e.t.c;
  # will raise stop if convergence < 0, i.e. error, or value = -Inf.

  res <- tryCatch(
    NLOptim(
      x0,
      fn, gr,
      lower, upper,
      method, nl.control, ...
    ),
    error = function(e){
      list(value = NA, convergence = -10, message = as.character(e))
    }
  )
  if (res$convergence > 0 && is.finite(res$value)) {
    return(res)
  } else {
    b.res <- tryCatch(
      NLOptim(
        x0,
        fn, gr,
        lower, upper,
        backup.method, nl.control, ...
      ),
      error = function(e){
        list(value = NA, convergence = -10, message = as.character(e))
      }
    )
    if (b.res$convergence > 0 && is.finite(b.res$value)) {
      return(b.res)
    }
  }

  # if unfortunately end up here, either both fails (code < 0) or
  # both being -Inf, in such case code = 2, stopval = -Inf by default.
  # if set larger stopval, will not see this.
  stop(paste(
    method, 'fails, code', res$convergence, ';\n', res$message, '\n',
    backup.method, 'fails, code', b.res$convergence, ';\n', b.res$message
  )
  )
}

NLOptim <- function(
  x0,
  fn, gr = NULL,
  lower, upper,
  method = "DIRECT-P", nl.control = list(), ...){
  # nonlinear optimization (minimization) with box constraints using nloptr
  # args:
  # x0: initial point, missing -> midpoint of the box, omitted for some methods;
  # fn, gr: objective function and its gradient (optional);
  # lower, upper: bounds;
  # method: method to use, say DIRECT, DIRECT-P, LBFGS, MLSL, Nelder-Mead, SubPlex;
  # nl.control: control list to pass to nloptr functions, c.f. nloptr::nl.opts;
  # ...: additional arguments passed to fn and gr.
  # returns:
  # what nloptr returns, i.e., a list of par, value, iter, convergence and message,

  method <- match.arg(
    arg = method,
    choices = c("DIRECT", "DIRECT-P", "LBFGS", "MLSL", "Nelder-Mead", "SubPlex")
  )

  if (!(method %in% c("DIRECT", "DIRECT-P")) && missing(x0)) {
    x0 = (lower + upper) / 2
  }

  res <- switch(
    method,
    "DIRECT" = {
      nloptr::direct(
        fn = fn,
        lower = lower, upper = upper,
        control = nl.control, ...
      )
    },
    "DIRECT-P" = {
      directPolish(
        fn = fn, gr = gr,
        lower = lower, upper = upper,
        control = nl.control, ...
      )
    },
    "LBFGS" = {
      nloptr::lbfgs(
        x0 = x0,
        fn = fn, gr = gr,
        lower = lower, upper = upper,
        control = nl.control, ...
      )
    },
    "MLSL" = {
      nloptr::mlsl(
        x0 = x0,
        fn = fn, gr = gr,
        lower = lower, upper = upper,
        # low.discrepancy = TRUE,
        control = nl.control, ...
      )
    },
    "Nelder-Mead" = {
      nloptr::neldermead(
        x0 = x0,
        fn = fn,
        lower = lower, upper = upper,
        control = nl.control, ...
      )
    },
    "SubPlex" = {
      nloptr::sbplx(
        x0 = x0,
        fn = fn,
        lower = lower, upper = upper,
        control = nl.control, ...
      )
    }
  )

  return(res)
}

directPolish <- function(fn, gr = NULL, lower, upper, control = list(), ...){
  # DIRECT-P for optimization
  # args:
  # fn, gr: function to be minimized and its gradient;
  # lower, upper: bounds;
  # control: control list to be passed to nloptr functions;
  # ...: additional arguments for fn and gr.
  # returns:
  # what nloptr returns, i.e., a list of par, value, iter, convergence and message.
  # details:
  # The method DIRECT-P will return the final result generated by L-BFGS
  # initialized with result from a piloting DIRECT method, unless it fail. In such case,
  # switch to Subplex method and try again with same initialization. Note If also fail,
  # result from DIRECT method is return with a warning. (polishing fails.)

  direct.res <- nloptr::direct(
    fn = fn,
    lower = lower, upper = upper,
    control = control, ...
  )

  res <- nloptr::lbfgs(
    x0 = direct.res$par,
    fn = fn,
    gr = gr,
    lower = lower,
    upper = upper,
    control = control,
    ...
  )
  if (res$convergence > 0 & res$value <= direct.res$value) {
    return(res)
  } else {
    res <- nloptr::sbplx(
      x0 = direct.res$par,
      fn = fn,
      lower = lower,
      upper = upper,
      control = control,
      ...
    )
    if (res$convergence > 0) {
      return(res)
    } else {
      warning('polishing fails.')
      return(direct.res)
    }
  }
}



#LEGACY Wrapper for built-in optim ============================================

smartOptim = function(par, fn, method = 'Nelder-Mead', multiStart = c(1,0), ...){
  # a wrapper for optim in R so as to handle errors
  # only for length(par) >= 2
  # par = initial value
  # fn = objective function to be minimized
  # multiStart = number of initial values to be generated, and variance of it.
  #              for example, c(3, 2, 3, 4) means in addition to par whose length = 3, also initial at
  #              replicate(2, par + runif(length(par), min = -c(2,3,4), max = c(2,3,4)))
  #              if 1st number = 1, 2nd number ignored.
  #              in short, length(multiStart) = 2 or length(par) + 1
  #              if lenth(multiStart) = 2, the last element is recycled
  # ... : additional control argument passing to optim
  if(multiStart[1] == 1){
    par = list(par)
  }else{
    if(length(multiStart) != 2 & length(multiStart) - 1 != length(par)) stop('check length of par and multiStart.')
    par = c(
      list(par),
      replicate(multiStart[1] - 1,
                par + runif(length(par), min = -multiStart[-1], max = multiStart[-1]),
                simplify = F)
    )
  }

  # idx.init = 1

  opt_res = list()
  # -1 for not touched
  # 0 for works but may be cov.code
  # 1 for totally broken
  # though make thing overly slow, not implemented.
  broken_opt = rep(-1, length(par))
  idx.optim = 1

  # get a list of methods with the specified one at first
  lst_method = c("Nelder-Mead", "BFGS", "CG", "SANN")
  lst_method = c(method, lst_method[!(lst_method == method)])
  idx = 1
  for(idx.init in 1:length(par)){
    flg.nextMethod = TRUE
    while(flg.nextMethod){
      opt_res[[idx.init]] = tryCatch(
        {optim(par = par[[idx.init]], fn = fn, method = lst_method[idx], ...)},
        error = function(e){list(err = e, is.err = TRUE)}
      )
      # if optim works, gooooood
      if(is.null(opt_res[[idx.init]]$is.err)){
        flg.nextMethod = FALSE
        # This will make things too slow, impractical
        # if(opt_res[[idx.init]]$convergence == 0){
        #   flg.nextMethod = FALSE
        #   broken_opt[idx.init] = 0
        # }

        # label the minimum result, and note down its method
        if(opt_res[[idx.init]]$value <= opt_res[[idx.optim]]$value){
          idx.optim = idx.init
          opt_res[[idx.optim]]$method = lst_method[idx]
        }
      }else{
        # if all methods do not work
        if(idx == length(lst_method)){
          broken_opt[idx.init] = 1
          # if(broken_opt[idx.init] == -1) broken_opt[idx.init] = 1
          flg.nextMethod = FALSE
        }else{
          # otherwise try next
          idx = idx + 1
        }
      }
    }
  }
  # if all methods fail for all initial values
  if(all(broken_opt == 1)) stop('Optim broken.')

  # now just keep the best opt_res
  opt_res = opt_res[[idx.optim]]
  msg = NULL
  if(opt_res$method != method) msg = paste0(msg, 'Optim using ', opt_res$method, '. ')
  if(opt_res$convergence != 0) msg = paste0(msg, 'Convergence code ', opt_res$convergence, '.')
  if(!is.null(msg)) warning(paste0(msg, '\n'))
  return(opt_res)
}

