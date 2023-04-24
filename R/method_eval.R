# Into modular ============================

#' Prepare simulation data
#'
#' @param density.family density family.
#' @param n.train number of training samples.
#' @param train.size number of observations per training sample. Set to `Inf`
#' to obtain complete density curves.
#' @param n.test number of training samples.
#' @param test.size number of observations per testing sample.
#' @param grid.size when `train.size = Inf`, the number of grid points.
#' @param control omitted.
#'
#' @return A list of matrices of training samples, testing samples and true
#' parameters, one row for one sample. Grid will be NA if `train.size != Inf`,
#' with `grid.size` ignored.
#' @export
#'
#' @examples
#' dat <- genSimData(tNrml(3), 10, 100, 3, 5)
#' str(dat)
genSimData <- function(
  density.family, n.train, train.size, n.test, test.size,
  grid.size = 128, control = list()) {
  # generate simulation data given density family
  # args:
  #   density.family: a density family.
  #   n.train: number of training samples to generate.
  #   train.size: number of observations per training sample.
  #               return complete density if = Inf(default)
  #   n.test: number of testing samples.
  #   test.size: number of observations per testing sample.
  #   grid.size: when train.size = Inf, the number of grid points.
  #TBD
  #   control: a list of controling arguments:
  #     gen.par = "random"(default)/"lattice" how parameters are generated
  #     gen.train.par = same as gen.par, only affects training samples
  #     gen.test.par = same as gen.par, only affects test samples
  #TBD
  # returns:
  #   A list consists of training samples, testing samples and true parameters,
  #   one row for one sample.
  #   grid will be NA if train.size != Inf, and grid.size ignored.

  # DEBUG
  # density.family <- tNrml(6)
  # n.train <- 100
  # train.size <- 10
  # n.test <- 20
  # test.size <- 5
  # grid.size <- 128
  # set.seed(100)

  # sanity check
  if (any(!is.finite(test.size)))
    stop('infinite test.size.')
  if (any(!is.finite(train.size)) & length(train.size) > 1)
    stop('train.size need to be either all finite or all infinite.')
  if (!all(c('pdf', 'rpdf', 'gen_par', 'range') %in% names(density.family)))
    stop('density.family missing element.')

  # some renaming and some default
  df.pdf <- density.family$pdf
  df.rpdf <- density.family$rpdf
  df.genPar <- density.family$gen_par
  df.range <- density.family$range
  grid <- NA

  # generating parameters
  train.par <- df.genPar(n.train)
  test.par <- df.genPar(n.test)

  # recycle train.size and/or test.size
  train.size <- rep(train.size, length.out = n.train)
  test.size <- rep(test.size, length.out = n.test)

  # getting training samples
  if (all(train.size == Inf)) {
    grid <- seq(from = min(df.range), to = max(df.range), length.out = grid.size)
    train.sample <- apply(train.par, 1, function(par){df.pdf(x = grid, par = par)})
    train.sample <- t(train.sample) # so that 1 row ~ 1 sample
  } else {
    train.sample <- mapply(
      function(n, par){df.rpdf(n = n, par = par)},
      n = train.size, par = mat2list(train.par, by = 1)
    )
    # if train.size all equal, make it a matrix
    if (all(train.size == train.size[1]))
      train.sample <- t(as.matrix(train.sample))
  }

  # getting testing samples
  test.sample <- mapply(
    function(n, par){df.rpdf(n = n, par = par)},
    n = test.size, par = mat2list(test.par, by = 1)
  )
  # if test.size all equal, make it a matrix
  if (all(test.size == test.size[1]))
    test.sample <- t(as.matrix(test.sample))

  return(
    list(
      train.par = train.par,
      test.par = test.par,
      train.sample = train.sample,
      test.sample = test.sample,
      grid = grid
    )
  )
}


#' Presmoothing with KDE
#'
#' @param obsv a matrix of observations, one row is one sample; or a list of
#' samples.
#' @param grid grid points to evaluate estimated density.
#' @param kde.opt list of options for `stats::density`.
#'
#' @return KDE of the samples in the form of a `nrow(obsv) X length(grid)`
#' matrix.
#' @export
#'
#' @examples
#' den.fam <- biM_family()
#' ls.obsv <- with(
#'   den.fam, apply(
#'     gen_par(5), 1, function(par) rpdf(rpois(1, 50), par),
#'     simplify = F
#'   )
#' )
#' grid <- seq(min(den.fam$range), max(den.fam$range), length.out = 1024)
#' mat.pdf <- preSmooth.kde(ls.obsv, grid = grid, kde.opt = list(bw = 'sj'))
#' matplot(x = grid, y = t(mat.pdf), type = 'l')
preSmooth.kde <- function(obsv, grid, kde.opt = list()){
  # Presmoothing with stats::density, return values interpolated with stats::approx.
  # args:
  # obsv: a matrix of observations, one row is one sample, or
  #       a list of arrays where each one is one sample.
  # grid: grid points to evaluate estimated density
  # kde.opt: list of options for stats::density
  # return:
  # a nrow(obsv) X length(grid) matrix.
  # Note: the returned values are not necessarily normalized.

  # set.seed(1)
  # obsv <- matrix(runif(250), ncol = 50)
  # breaks <- seq(0, 1, length.out = 10)
  if(!is.list(obsv) & !is.matrix(obsv)){
    stop('check input type, must be list or matrix.')
  }
  if(is.matrix(obsv)){
    obsv <- mat2list(obsv, by = 1)
  }
  n.sample <- length(obsv)
  mat.pdf <- matrix(0, ncol = length(grid), nrow = n.sample)

  for(i in seq_len(nrow(mat.pdf))){
    tm.res <- do.call(
      stats::density,
      c(list(x = obsv[[i]]), kde.opt)
    )
    mat.pdf[i, ] <- approx(
      x = tm.res$x, y = tm.res$y, xout = grid,
      rule = 1
    )$y
  }

  return(mat.pdf)
}

#' Presmoothing with KDE
#'
#' @param obsv a matrix of observations, one row is one sample; or a list of
#' samples.
#' @param grid grid points to evaluate estimated density.
#' @param lsp.opt a list of options for `logspline::logspline`.
#'
#' @return A `nrow(obsv) X length(grid)` matrix, one row for one sample.
#' Note that the returned densities are not necessarily normalized (integrate
#' to one).
#'
#' @export
#'
#' @examples
#' den.fam <- biM_family()
#' ls.obsv <- with(
#'   den.fam, apply(
#'     gen_par(5), 1, function(par) rpdf(rpois(1, 50), par),
#'     simplify = F
#'   )
#' )
#' grid <- seq(min(den.fam$range), max(den.fam$range), length.out = 1024)
#' mat.pdf <- preSmooth.logspline(ls.obsv, grid = grid, lsp.opt = list(
#'   lbound = min(den.fam$range), ubound = max(den.fam$range)
#' ))
#' matplot(x = grid, y = t(mat.pdf), type = 'l')
preSmooth.logspline <- function(obsv, grid, lsp.opt = list()){
  # Presmoothing with logspline::logspline
  # args:
  # obsv: a matrix of observations, one row is one sample, or
  #       a list of arrays where each one is one sample.
  # grid: grid points to evaluate estimated density
  # lsp.opt: list of options for logspline::logspline,
  #          for example, lbound and ubound.
  # return:
  # a n.sample X length(grid) matrix.
  # Note: the returned values are not necessarily normalized.

  # set.seed(1)
  # obsv <- matrix(runif(250), ncol = 50)
  # breaks <- seq(0, 1, length.out = 10)

  if(!is.list(obsv) & !is.matrix(obsv)){
    stop('check input type, must be list or matrix.')
  }
  if(is.matrix(obsv)){
    obsv <- mat2list(obsv, by = 1)
  }
  n.sample <- length(obsv)
  mat.pdf <- matrix(0, ncol = length(grid), nrow = n.sample)
  for(i in seq_len(nrow(mat.pdf))){
    tm.res <- do.call(
      logspline::logspline,
      c(list(x = obsv[[i]]), lsp.opt)
    )
    mat.pdf[i, ] <- logspline::dlogspline(q = grid, fit = tm.res)
  }
  return(mat.pdf)
}

preSmooth.locfit <- function(obsv, grid, lp.opt = list()){
  # Presmoothing with locfit::locfit(~lp(obsv))
  # args:
  # obsv: a matrix of observations, one row is one sample, or
  #       a list of arrays where each one is one sample.
  # grid: grid points to evaluate estimated density
  # lp.opt: list of options for locfit::lp,
  #          for example, h and deg.
  # return:
  # a nrow(obsv) X length(grid) matrix.
  # Note: the returned values are not necessarily normalized.

  # set.seed(1)
  # obsv <- matrix(runif(250) + 2, ncol = 50)
  # lp.opt <- list(deg = 1)
  # grid <- seq(2, 3, length.out = 128)

  if(!is.list(obsv) & !is.matrix(obsv)){
    stop('check input type, must be list or matrix.')
  }
  if(is.matrix(obsv)){
    obsv <- mat2list(obsv, by = 1)
  }
  n.sample <- length(obsv)
  mat.pdf <- matrix(0, ncol = length(grid), nrow = n.sample)
  for(i in seq_len(nrow(mat.pdf))){
    tm.lp <- do.call(locfit::lp, c(list(obsv[[i]]), lp.opt))
    tm.res <- locfit::locfit(
      ~tm.lp,
      ev = locfit::lfgrid(mg = 512)
    )
    mat.pdf[i, ] <- locfit:::predict.locfit(tm.res, grid)
  }

  # plot(tm.res)
  # lines(grid, mat.pdf[i, ] / ( sum(mat.pdf[i, ]) * mean(diff(grid))), col = 3)
  # lines(density(obsv[i, ], bw = 'sj'), col = 2)

  return(mat.pdf)
}


#' Translate density function into Hilbert space.
#'
#' @param mat.curve a matrix of density functions. One row for one curve.
#' @param grid on which the curves are evaluated.
#' @param transFun what transform function to use, default 'log'.
#' @param eps value smaller than esp will be treated as `NA`.
#'
#' @return a list of
#'   - `mat.curve`: the trajectories after transformation;
#'   - `grid`: grid;
#'   - `idx.drop`: an array of `TRUE`/`FALSE` indicating whether the original
#'   curve was dropped;
#'   - `ref.pt`: the index of the grid point used as reference. `NULL` unless
#'   `transFun = 'log'`.
#' @export
#'
#' @details The default `transFun = 'log'` is only applying log-transformation,
#' for centered log-transformation (centered log-ratio), use `orthLog`.
#' To use custom transformation function, supply `transFun` with
#' `function(f, grid)`, where `f` and `grid` are arguments for the density and
#' grid, respectively.
#'
#' @examples
#' den.fam <- tNrml(1)
#' ls.obsv <- with(
#'   den.fam, apply(
#'     gen_par(5), 1, function(par) rpdf(rpois(1, 500), par),
#'     simplify = F
#'   )
#' )
#' grid <- seq(min(den.fam$range), max(den.fam$range), length.out = 1024)
#' mat.pdf <- preSmooth.kde(ls.obsv, grid = grid, kde.opt = list(bw = 'sj'))
#' ls.hilbert <- toHilbert(
#'   mat.pdf, grid,
#'   transFun = function(f, grid) orthLog(f, grid, against = 1),
#'   eps = .Machine$double.eps^(1/2)
#' )
#' matplot(x = ls.hilbert$grid, y = t(ls.hilbert$mat.curve), type = 'l')
toHilbert <- function(mat.curve, grid, transFun = 'log', eps = .Machine$double.eps^(1/3)){
  # transform the curves to a (hopefully) Hilbert space.
  # args:
  # mat.curve: matrix of curves, one row is one curve
  # grid: on which the curves are evaluated
  # transFun: what transform function to use, default 'log'.
  # (for log transform only)
  #   eps: value smaller than esp will be treated as NA.
  # All curves with nonfinite values will be dropped.
  # return:
  # a list of
  #   mat.curve: the trajectories after transformation;
  #   grid: grid;
  #   idx.drop: an array of TRUE/FALSE indicating whether the original curve was dropped;
  #   ref.pt: the index of the grid point used as reference.

  # sanity check
  if (ncol(mat.curve) != length(grid)) {
    stop('check input dimension: ncol(mat.curve) != length(grid).')
  }

  if (class(transFun) != 'function'){
    transFun <- match.arg(transFun, c('log'))

    if (transFun == 'log'){
      idx.drop <- (rowSums(abs(mat.curve) < eps) != 0)
      if (any(idx.drop)) {
        if (sum(idx.drop) == nrow(mat.curve))
          stop('All curves are dropped, consider change eps.')
        warning(paste(
          'Dropping logInf samples,',
          100 * sum(idx.drop) / length(idx.drop),
          '% dropped.\n'
        ))
        mat.curve <- mat.curve[!idx.drop, ]
      }

      ref.pt <- ceiling(ncol(mat.curve) / 2)
      mat.after <- log(mat.curve) - log(mat.curve[, ref.pt])
    }

  }else{
    ### TBD
    w.transFun <- function(f, grid){
      tryCatch(
        transFun(f, grid = grid),
        error = function(e) rep(NA, length(grid))
      )
    }
    mat.after <- t(apply(
      mat.curve, 1,
      w.transFun, grid = grid
    ))
    idx.drop <- (rowSums(!is.finite(mat.after)) != 0)
    if (any(idx.drop)) {
      if (sum(idx.drop) == nrow(mat.curve))
        stop('All curves are dropped, consider change transformation.')
      warning(paste(
        'Dropping curves with nonfinite values,',
        100 * sum(idx.drop) / length(idx.drop),
        '% dropped.\n'
      ))
      mat.after <- mat.after[!idx.drop, ]
    }
    ref.pt <- NULL
  }

  return(list(
    mat.curve = mat.after,
    grid = grid,
    idx.drop = idx.drop,
    ref.pt = ref.pt
  ))

  # ls.res <- list(
  #   Ly = lapply(seq_len(nrow(mat.pdf)), function(i) mat.pdf[i, ]),
  #   Lt = replicate(n = nrow(mat.pdf), expr = grid, simplify = FALSE)
  # )
}

assErr <- function(actual, esti, grid, err.method, idx.df){
  # assessing error with err.method
  # args:
  # actual: a matrix of actual curves, one row = one curve.
  # esti: a matrix of estimated curves, one row = one curve.
  # grid: grid on which all curves sits, length must = ncol(actual) & ncol(esti).
  # err.method: a character array of names of error methods passed to calErr.
  # idx.df(optional): a data.frame(or matrix) to index the results,
  #                   if missing, filled as 1:n.
  # returns:
  # a list of
  # a tbl_df in long format, one row = one type of error of one pair of actual vs esti.
  # and timing.

  # DEBUG
  # dummy <- genSimData(tNrml(6), n.train = 40, train.size = Inf, n.test = 20, test.size = 10)
  # actual <- dummy$train.sample[1:20, ]
  # esti <- dummy$train.sample[21:40, ]
  # grid <- dummy$grid
  # err.method <- c('ISE', 'KL', 'rev.KL', 'FisherRao', 'Wass..5', 'Wass.1', 'Wass.2')
  # idx.df <- dplyr::as_tibble(expand.grid(idx.test = 1:5, idx.train = 1:4))

  # sanity check
  if (length(grid) != ncol(actual) | length(grid) != ncol(esti)) {
    stop('check input dimension: ncol != length(grid).')
  }
  if (nrow(actual) != nrow(esti)) {
    stop('check input dimension: nrow of actual and esti.')
  }

  if (missing(idx.df)) {
    idx.df <- data.frame(idx = seq_len(nrow(actual)))
  }

  # fill in default
  if (missing(err.method)) {
    err.method <- c('ISE', 'KL', 'rev.KL', 'FisherRao', 'Wass..5', 'Wass.1', 'Wass.2')
  }
  err.method <- sort(err.method)

  # w.actual <- as.list(data.frame(t(actual)))
  # w.esti <- as.list(data.frame(t(esti)))
  # tm.err <- purrr::map2(
  #   w.actual, w.esti,
  #   .f = calErr, grid = grid, method = err.method
  #   )

  err.t <- system.time({
    tm.err <- matrix(0, ncol = length(err.method), nrow = nrow(actual))
    for (i in seq_len(nrow(actual))) {
      tm.err[i, ] <- calErr(f = actual[i, ], hatf = esti[i, ],
                            grid = grid, method = err.method)
    }
    colnames(tm.err) <- err.method
    # tm.df <- dplyr::bind_cols(idx.df, dplyr::as_tibble(tm.err))
    # res.df <- tidyr::gather(tm.df, key = "type.err", value = "err",
    #                         err.method, factor_key = FALSE)

    # tm.df <- dplyr::bind_cols(idx.df, as.data.frame(tm.err))
    ls.df <- list()
    for(i in seq_along(err.method)){
      ls.df[[i]] <- data.frame(
        type.err = rep(err.method[i], nrow(actual)),
        err = tm.err[, i],
        stringsAsFactors = FALSE
      )
      ls.df[[i]] <- dplyr::bind_cols(idx.df, ls.df[[i]])
    }
    res.df <- dplyr::bind_rows(ls.df)
  })
  return(
    list(
      err = res.df,
      time = err.t
    )
  )
}

assErr.cv <- function(obsv, esti.f, grid, err.method, idx.df, control = list()){
  # assessing error with err.method
  # args:
  # obsv: a list/matrix of samples, one element/row = one sample.
  # esti.f: a function taking samples giving density estimates as a function
  #         i.e., esti.f(obsv[[i]]) is a function; or a list where one of its
  #         elements is such function, other elements could be used for
  #         auxiliary results, here one need to specify the key to the function
  #         in control.
  # grid: an array for grid on which to evaluate functions
  # err.method: a character array, 'cv.ISE', 'cv.KL', 'cross.entropy'.
  # idx.df(optional): a data.frame(or matrix) to index the samples,
  #                   if missing, filled as 1:length(obsv).
  # control: a list of control arguments
  #     cv.cntrl: a list for args of kFoldIdx k & n.rep, missing --> LOOcv
  #     return.esti: TRUE/FALSE(default) or a chracter vector, whether or
  #                  which of the estimates would be returned. If TRUE, the
  #                  estimates will be returned in a list corresponding to the
  #                  order of obsv, note this is esti.f(obsv[[i]]), not those
  #                  used in LOOCV or CV.
  #     key.fhat: (optional) a character string, if esti.f(obsv[[i]]) returns
  #               a list, this will be used to index the function inside.
  #               If missing, then esti.f(obsv[[i]]) should return a function.
  # returns:
  # a data.frame in long format, one row = one type of error of one sample,
  # unless return.esti is not FALSE, then the data.frame and a list of estimates
  # named after obsv.

  # DEBUG
  # browser()
  # err.method <- c('cv.KL', 'cv.ISE', 'cross.entropy')
  # den.fam <- tNrml(6)
  # range <- den.fam$range
  # set.seed(1)
  # dat <- genSimData(
  #   den.fam, n.train = 40, train.size = 20, n.test = 1, test.size = 1
  # )
  # obsv <- dat$train.sample
  # esti.f <- fpcaEsti.f(
  #   learnFamily(tNrml(6), trainSize = 50, maxK = 20, gridSize = 512),
  #   c('FPCA_MLE', 'FPCA_MAP', 'FPCA_BLUP'),
  #   list(
  #     num.k = 'AIC', max.k = 10,
  #     method = 'LBFGS', return.scale = 'function',
  #     extend.par.range = 3
  #   )
  # )
  # esti.f <- esti.f$FPCA_BLUP
  # # esti.f <- function(x){
  # #   tm.pdf <- density(
  # #     x = x, bw = 'sj', from = den.fam$range[1], to = den.fam$range[2]
  # #   )
  # #   approxfun(x = tm.pdf$x, y = tm.pdf$y, yleft = 0, yright = 0)
  # # }
  # # esti.f <- function(x){
  # #   esti.par <- calMLE(x, den.fam, list(method = 'LBFGS'))
  # #   function(x) den.fam$pdf(x, esti.par)
  # # }
  # grid <- seq(den.fam$range[1], den.fam$range[2], length.out = 512)
  # control <- list(
  #   cv.cntrl = list(k = 5, n.rep = 0)
  # )
  # END DEBUG

  # sanity check
  err.method <- match.arg(
    err.method,
    c('cv.ISE', 'cv.KL', 'cross.entropy'),
    several.ok = TRUE
  )

  if(!is.list(obsv)){
    if(is.matrix(obsv)){
      obsv <- mat2list(obsv, by = 1)
    }else{
      stop('check input type: obsv')
    }
  }
  len.obsv <- sapply(obsv, length)
  if(any(len.obsv <= 1))
    stop('check input sample size: need at least 2 for cross validation')
  if(missing(idx.df))
    idx.df <- data.frame(idx.obsv = seq_along(obsv))
  if(is.null(control$return.esti))
    control$return.esti <- FALSE
  # stopifnot(is.logical(control$return.esti) & length(control$return.esti) == 1)

  # construct CV data
  loo.obsv <- list()
  evl.obsv <- list()
  if(is.null(control$cv.cntrl)){
    for(idx in seq_along(obsv)){
      tm.obsv <- obsv[[idx]]
      # LOOCV
      loo.obsv[[idx]] <- lapply(
        seq_len(length(tm.obsv)),
        function(i) tm.obsv[-i]
      )
      evl.obsv[[idx]] <- as.list(tm.obsv)
    }
  }else{
    for(idx in seq_along(obsv)){
      tm.obsv <- obsv[[idx]]
      # k-fold
      tm.idx <- do.call(
        kFoldIdx,
        c(
          list(len.obsv[idx]),
          control$cv.cntrl$k, control$cv.cntrl$n.rep, TRUE
        )
      )
      # some adjustment for diff kFoldIdx output structure
      # if(control$cv.cntrl$n.rep != 0)
      #   tm.idx <- unlist(tm.idx, recursive = FALSE)

      loo.obsv[[idx]] <- lapply(
        tm.idx,
        function(i) tm.obsv[i[[1]]]
      )
      evl.obsv[[idx]] <- lapply(
        tm.idx,
        function(i) tm.obsv[i[[2]]]
      )
    }
  }

  # get esti.f at each loo.obsv
  ls.cv.fhat <- list()
  ls.itgl <- list()
  ls.eval <- list()
  ls.fhat <- list()
  if(is.null(control$key.fhat)){ # if esti.f(obsv[[i]]) is function
    for(idx in seq_along(loo.obsv)){
      ls.cv.fhat[[idx]] <- lapply(
        loo.obsv[[idx]], esti.f
      )
      ls.eval[[idx]] <- mapply(
        function(f, arg) do.call(f, list(arg)),
        ls.cv.fhat[[idx]],
        evl.obsv[[idx]],
        SIMPLIFY = FALSE
      )
      nodrop.fhat <- esti.f(obsv[[idx]])
      ls.fhat[[idx]] <- nodrop.fhat  # recording
      nodrop.fhat <- nodrop.fhat(grid)
      # R(f hat) = \int fhat^2
      ls.itgl[[idx]] <- dotL2(
        nodrop.fhat,
        nodrop.fhat,
        grid = grid
      )
    }
  }else{ # if esti.f(obsv[[i]]) is list
    w.esti.f <- function(obsv) (esti.f(obsv))[[control$key.fhat]]
    for(idx in seq_along(loo.obsv)){
      ls.cv.fhat[[idx]] <- lapply(
        loo.obsv[[idx]], w.esti.f
      )
      ls.eval[[idx]] <- mapply(
        function(f, arg) do.call(f, list(arg)),
        ls.cv.fhat[[idx]],
        evl.obsv[[idx]],
        SIMPLIFY = FALSE
      )
      nodrop.fhat <- esti.f(obsv[[idx]])
      ls.fhat[[idx]] <- nodrop.fhat  # recording
      nodrop.fhat <- nodrop.fhat[[control$key.fhat]](grid)
      # R(f hat) = \int fhat^2
      ls.itgl[[idx]] <- dotL2(
        nodrop.fhat,
        nodrop.fhat,
        grid = grid
      )
    }
  }
  # cv.error
  ls.df <- list()
  for(w.err in err.method){
    if(w.err == 'cv.KL'){  ## cv.KL
      ls.df[[w.err]] <- data.frame(
        idx.df, type.err = w.err,
        err = -1 *
          sapply(ls.eval, function(x) mean(log(unlist(x)), na.rm = FALSE)),
        stringsAsFactors = FALSE
      )
    }
    if(w.err == 'cv.ISE'){  ## cv.ISE (ucv)
      ls.df[[w.err]] <- data.frame(
        idx.df, type.err = w.err,
        err = unlist(ls.itgl) -
          2 * sapply(ls.eval, function(x) mean(unlist(x), na.rm = FALSE)),
        stringsAsFactors = FALSE
      )
    }
    if(w.err == 'cross.entropy'){  ## cross-entropy
      tm.eval <- mapply(
        function(fhat, x) fhat[[control$key.fhat]](x),
        fhat = ls.fhat, x = obsv, SIMPLIFY = FALSE
      )
      ls.df[[w.err]] <- data.frame(
        idx.df, type.err = w.err,
        err = -1 * sapply(tm.eval, function(x) mean(log(unlist(x)), na.rm = FALSE)),
        stringsAsFactors = FALSE
      )
    }
  }
  names(ls.fhat) <- names(obsv)
  if(control$return.esti[1] == TRUE){
    return(list(
      err = dplyr::bind_rows(ls.df),
      esti = ls.fhat
    ))
  }
  if(control$return.esti[1] == FALSE){
    return(dplyr::bind_rows(ls.df))
  }
  # otherwise
  for(i in seq_along(ls.fhat)){
    ls.fhat[[i]] <- ls.fhat[[i]][as.character(control$return.esti)]
  }
  return(list(
    err = dplyr::bind_rows(ls.df),
    esti = ls.fhat
  ))
}

# for convenience only
learnFamily <- function(
  density_family, trainSize, maxK = 5, gridSize = 128,
  refPoint = 0, simple_result = TRUE
){
  # return a FPCA object learning generated trajectories from density_family
  # unless simple_result = FALSE, then also returning generated parameters.

  pdf = density_family$pdf
  rpdf = density_family$rpdf
  gen_par = density_family$gen_par
  range = density_family$range

  grid = seq(from = min(range), to = max(range), length.out = gridSize)
  Lt = replicate(n = trainSize, expr = grid, simplify = F)
  train_par = gen_par(trainSize)
  Ly = apply(train_par, 1, function(par){log(pdf(grid, par)) - log(pdf(refPoint, par))})
  Ly = as.list(data.frame(Ly))
  names(Ly) = rep('', trainSize)
  fpca_res = fdapace::FPCA(
    Ly, Lt,
    list(
      plot = F, error = F, useBinnedData = 'OFF', maxK = maxK,
      FVEthreshold = 1
    )
  )
  # limit maximum number of components to avoid overfitting
  # TBD: FVE cutoff instead of maxK
  # if(length(fpca_res$lambda) > maxK){
  #   fpca_res$lambda = (fpca_res$lambda)[1:maxK]
  #   fpca_res$phi = (fpca_res$phi)[, 1:maxK, drop = FALSE]
  #   fpca_res$xiEst = (fpca_res$xiEst)[, 1:maxK, drop = FALSE]
  # }
  if(simple_result){
    return(
      fpca_res
    )
  }else{
    return(
      list(
        fpca_res = fpca_res,
        parameters = train_par
      )
    )
  }
}
