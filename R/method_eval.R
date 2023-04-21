# Into modular ============================

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

# seemingly no longer necessary
# fitData <- function(new.obsv, train.sample, grid, esti.method, esti.option, fpca.option){
#   # non-parametric density estimation of new.obsv based on FPCA upon train.sample.
#   #
#   # args:
#   # new.obsv: a list of new observations, one element = one sample.
#   # train.sample: a matrix of completely observed training curves.
#   # grid: grid on which train.sample and estimated densities lies.
#   # esti.method: a character array of names.
#   # esti.option: list of options for estimating pass to fpcaEsti.
#   # fpca.option: list of options for fdapace::FPCA
#   #
#   # returns:
#   # a list of
#   #   a matrix where one row is one esti pdf of newobsv[idx.new, ] with one method;
#   #   a data.frame for indexing the esti matrix;
#   #   a FPCA result;
#   #   a data.frame for timing.
#
#   # DEBUG
#   # set.seed(1)
#   # dummy <- genSimData(tNrml(6), 50, Inf, 5, 5)
#   # train.sample <- dummy$train.sample
#   # new.obsv <- dummy$test.sample
#   # grid <- dummy$grid
#   # esti.method <- c('FPCA_MLE', 'FPCA_BLUP')
#
#   # fill in default esti.method
#   if (missing(esti.method)) {
#     esti.method <- names(kListFpcaMethod)
#   }
#   # fill in default esti.option
#   if (missing(esti.option)) {
#     esti.option <- list()
#   }
#   default.esti.option <- list(
#     multiStart = FALSE,
#     fpca.k = "all",
#     scale = 'origin'
#   )
#   for (nm in names(default.esti.option)) {
#     if (is.null(esti.option[[nm]])) {
#       esti.option[[nm]] <- default.esti.option[[nm]]
#     }
#   }
#   if (!all(esti.method %in% names(kListFpcaMethod))) {
#     stop('unknown esti.method.')
#   }
#
#   # train
#   fpca.res <- trainModel(train.sample, grid, fpca.option)
#
#   # before estimating
#   # some preparing, first, set multiStart
#   if(!esti.option$multiStart){
#     w.multiStart <-  c(1,0)
#   }else{
#     # determine initial values
#     tm.Range <- matrix(apply(fpca.res$xiEst, 2, range), nrow = 2)
#     # set as so since mean(xiEst) = 0, choose the hyper-rectangle.
#     w.multiStart <- c(
#       min(9, length(fpca.res$lambda)*2 - 1),
#       pmin(abs(tm.Range[2,]), abs(tm.Range[1,])) / 2
#     )
#   }
#   # fill in working esti.option
#   w.esti.option <- esti.option
#   w.esti.option$multiStart = w.multiStart
#   w.esti.option$scale = 'origin'
#
#   # estimating
#   ls.res <- list()
#   for (i in seq_len(length(esti.method))) {
#     # getting estimating functions
#     esti.fun <- fpcaEstiFun(esti.method[i], w.esti.option)
#     ls.res[[i]] <- fpcaEsti(fpca.res, new.obsv, esti.method = esti.fun)
#     # note down the name, for easier debugging
#     ls.res[[i]]$esti.method <- esti.method[i]
#   }
#
#   # DEBUG
#   # tst <- lapply(ls.res, function(x) t(x$mat.esti) )
#   # all(tst[[1]] == tst[[2]])
#   # p2 <- fpca_BLUP(new.obsv[1, ], fpca.res, scale = 'origin')
#   # p1 <- fpca_MLE(new.obsv[1, ], fpca.res, scale = 'origin')
#   # all(p1 == p2)
#   # all(tst[[1]][,1] == p1)
#   # all(tst[[2]][,1] == p2)
#
#   mat.esti <- do.call(
#     c,
#     lapply(ls.res, function(x) t(x$mat.esti) ) # a list of mat with col = 1 pdf.
#   )
#   mat.esti <- matrix(
#     mat.esti, nrow = length(grid),
#     ncol = length(esti.method) * nrow(new.obsv)
#   )
#   # now in format of 1 row = 1 pdf, and in trunc of esti.method
#   mat.esti <- t(mat.esti)
#   # indexing data.frame
#   df.idx <- data.frame(
#     idx.obsv = rep(seq_len(nrow(new.obsv)), times = length(esti.method)),
#     esti.method = rep(esti.method, each = nrow(new.obsv))
#   )
#
#   # get the time
#   df.t <- do.call(rbind, lapply(ls.res, '[[', 'time'))
#   df.t <- rbind(df.t, colSums(df.t))
#   df.t <- cbind(data.frame(step = c(esti.method, 'TOTAL')), df.t)
#
#   # return
#   return(
#     list(
#       mat.esti = mat.esti,
#       idx = df.idx,
#       time = df.t
#     )
#   )
# }

# old ways 09/26/19============================================================
# replaced by preSmooth.kde + toHilbert.
preTrain <- function(train.sample, grid, transFun, sample.size = Inf){
  # replaced by preSmooth.kde + toHilbert.
  # preprocessing the training samples and sanity check, format into FPCA Lt, Ly style.
  #
  # args:
  # train.sample: a matrix where one row is one sample.
  # grid: a vector on which the output train.sample curves are evaluated.
  #       if sample.size == Inf, grid should be that for train.sample.
  # transFun: the function to apply on each sample,
  #           2 args: a sample (a array) and grid, return a numeric array.
  # sample.size: Inf or anything else, means train.sample are curves or not.
  # returns:
  # a list of processed train.sample (now they are curves) in Ly and grid in Lt.

  # TODO: presmoothing, maybe put as stand along function.
  # if (is.finte(sample.size)) {
  #    ... to presmoothing
  # }

  # DEBUG
  # gen.dat <- genSimData(tNrml(6), n.train = 50, train.size = Inf, n.test = 50, test.size = 5)
  # train.sample <- gen.dat$train.sample
  # grid <- gen.dat$grid

  if (missing(transFun)) {
    transFun <- function(x, grid){
      idx.mid <- ceiling(length(grid)/2)
      return(
        log(x) - log(x[idx.mid])
      )
    }
  }

  after.proc <- apply(train.sample, 1, transFun, grid = grid)
  # sanity check, we do want to know if anything goes wrong
  if (any(!is.finite(after.proc)) | any(is.null(after.proc))) {
    stop('non finite value after preprocessing.')
  }

  fpca.ly <- lapply(seq_len(ncol(after.proc)), function(i) after.proc[, i])
  fpca.lt <- replicate(n = nrow(train.sample), expr = grid, simplify = FALSE)

  return(
    list(
      Ly = fpca.ly,
      Lt = fpca.lt
    )
  )
}

# Use preSmooth.kde + toHilbert + mat2list instead.
trainModel <- function(train.sample, grid, fpca.option){
  # Use preSmooth.kde + toHilbert + mat2list instead.
  # fill in default FPCA options
  if (missing(fpca.option)) {
    fpca.option <- list(
      error = FALSE,
      lean = TRUE,
      # maxK = 20,
      plot = FALSE,
      useBinnedData = 'OFF'
    )
  }
  # sanity check
  if (ncol(train.sample) != length(grid)) {
    stop('check input dimension: train.sample, grid.')
  }

  # pre-processing
  proc.train <- preTrain(train.sample = train.sample, grid = grid)
  proc.train$optns <- fpca.option
  # training
  fpca.res <- do.call(fdapace::FPCA, proc.train)
  return(fpca.res)
}


# foo <- function(n){
#   mat <- matrix(runif(100), nrow = n,10)
#   res <- rep(0,n)
#   for(i in 1:n){
#     res <- mean(mat[i,])
#   }
#   return(res)
# }
# app <- function(n){
#   mat <- matrix(runif(100), nrow = n,10)
#   res <- rep(0,n)
#   return(apply(mat, 1, mean))
# }
# lapp <- function(n){
#   mat <- matrix(runif(100), nrow = n,10)
#   res <- rep(0,n)
#   return(lapply(as.list(data.frame(t(mat))), mean))
# }
# microbenchmark::microbenchmark(
#   foo(100),
#   app(100),
#   lapp(100),
#   times = 1000
# )


# old ways ====================================================================
method_eval = function(density_family, nTrain, trainSize, nTest, testSize,
                       maxK = 5, gridSize = 128, refPoint = 0,
                       multiStart = FALSE, silent = FALSE, returnFPCA = TRUE, to01 = FALSE,
                       esti.method = c(
                         'FPCA_MLE',
                         'FPCA_MAP',
                         'FPCA_BLUP',
                         'FPCA_MKE',
                         'MLE',
                         'KDE'
                       ),
                       err.method = c(
                         'ISE',
                         'KL',
                         'rev.KL',
                         'FisherRao',
                         'Wass..5',
                         'Wass.1',
                         'Wass.2'
                       )){
  # density_family = list(
  #   pdf = function(x, par), value of density function at x given parameter, able to take vector for x
  #   rpdf = function(sample_size, par), generate iid random sample given parameter
  #   gen_par = function(num), generate num sets of parameters, return a matrix of num X dim(par)
  #   range = c(., .), a array of length 2, onwhich the family lies.
  # )
  # nTrain = number of times training performed
  # trainSize = number of training trajectories feed to FPCA
  # nTest = number of testing cases
  # testSize = sample size of a new observation sample for testing
  # gridSize = size of working grid
  # refPoint = reference point
  # multiStart = whether one should use multiple initial value during optim
  # to01 = whether domain of distribution should be transport to [0,1]
  # esti.method = list of estimating methods to compare, note: essentially no one should use emMixGauss, always truncate.
  # err.method = list of errors to use

  # DEBUG
  # set.seed(100)
  # density_family = beta_family
  # nTrain = 1
  # trainSize = 10
  # nTest = 1
  # testSize = 10
  # gridSize = 512
  # refPoint = 0.5
  #
  # density_family = nrml_family
  # nTrain = 1
  # trainSize = 50
  # nTest = 1
  # testSize = 2
  # gridSize = 512
  # refPoint = 0
  #
  # silent = FALSE
  # multiStart = FALSE
  # to01 = TRUE

  t_ttl = proc.time()
  if(
    all(is.element(err.method, set = c('ISE', 'KL', 'rev.KL', 'FisherRao', 'Wass..5', 'Wass.1', 'Wass.2'))) &
    all(is.element(esti.method, set = c('FPCA_MLE', 'FPCA_MAP', 'FPCA_BLUP', 'FPCA_MKE',
                                        'MLE', 'mleTruncMixGauss', 'emMixGauss', 'emTruncMixGauss', 'KDE')))
  ){
    err.name = err.method
    esti.name = esti.method[esti.method != 'KDE']
    # KDE will be handled individually to accommodate potentially different grids
    flg.kde = any(esti.method == 'KDE')
  }else{
    stop('Invalid methods.')
  }

  if(missing(nTrain)) nTrain = 1
  if(to01) density_family = wrapDenFam(density_family)

  pdf = density_family$pdf
  rpdf = density_family$rpdf
  gen_par = density_family$gen_par
  range = density_family$range
  grid = seq(from = min(range), to = max(range), length.out = gridSize)

  if(!silent) print('Training...')

  t_train = system.time({
    pca_res = replicate(nTrain,
                        learnFamily(density_family, trainSize = trainSize, maxK = maxK,
                                    gridSize = gridSize, refPoint = refPoint, simple_result = FALSE),
                        simplify = FALSE
    )
  })

  train_par = lapply(pca_res, "[[", "parameters")
  pca_res = lapply(pca_res, "[[", "fpca_res")

  if(!silent) print('Estimating...')
  # testing
  test_par = gen_par(nTrain * nTest)
  test_obs = t(apply(test_par, 1, function(par){rpdf(testSize, par)}))
  test_true = apply(test_par, 1, function(par){pdf(grid, par)})

  # note down input for later use
  input.multiStart = multiStart
  # just to translate into acceptable argument.
  if(!input.multiStart){
    multiStart = c(1,0)
  }else{
    # determine initial values
    tm.Range = matrix(apply(pca_res[[1]]$xiEst, 2, range), nrow = 2)
    # set as so since mean(xiEst) = 0, choose the hyper-rectangle so as no to go outside.
    multiStart = c(min(9, length(pca_res[[1]]$lambda)*2 - 1), pmin(abs(tm.Range[2,]), abs(tm.Range[1,])) / 2)
  }

  # list storing estimated densities and timing
  ls.res.esti = list()
  ls.t.esti = list()

  # a giant, ugly
  for(idx.esti in 1:length(esti.name)){

    # the estimated density is stored as a matrix
    # one column correspond to one estimated density (values at grid)
    # and the (nTest*(i-1) + 1) : (nTest*i) columns are fitted with
    # the ith training result.

    if(esti.name[idx.esti] == 'FPCA_MLE'){
      ls.t.esti[[idx.esti]] = system.time({
        test_fpcaMLE = sapply(1:nTrain,
                              function(i){
                                apply(test_obs[(nTest*(i-1)+1) : (nTest*i), , drop = FALSE], 1,
                                      function(newx){
                                        tryCatch(
                                          fpca_MLE(newx,
                                                   pca_res[[i]],
                                                   multiStart = multiStart,
                                                   scale = 'origin'),
                                          error = function(e){rep(NA, gridSize)}
                                        )
                                      })
                              },
                              simplify = FALSE)
        ls.res.esti[[idx.esti]] = do.call(cbind, test_fpcaMLE)
      })
    }

    if(esti.name[idx.esti] == 'FPCA_MAP'){
      ls.t.esti[[idx.esti]] = system.time({
        test_fpcaMAP = sapply(1:nTrain,
                              function(i){
                                apply(test_obs[(nTest*(i-1)+1) : (nTest*i), , drop = FALSE], 1,
                                      function(newx){
                                        tryCatch(
                                          fpca_MAP(newx,
                                                   pca_res[[i]],
                                                   multiStart = multiStart,
                                                   scale = 'origin'),
                                          error = function(e){rep(NA, gridSize)}
                                        )
                                      })
                              },
                              simplify = FALSE)
        ls.res.esti[[idx.esti]] = do.call(cbind, test_fpcaMAP)
      })
    }

    if(esti.name[idx.esti] == 'FPCA_BLUP'){
      ls.t.esti[[idx.esti]] = system.time({
        test_fpcaBLUP = sapply(1:nTrain,
                               function(i){
                                 apply(test_obs[(nTest*(i-1)+1) : (nTest*i), , drop = FALSE], 1,
                                       function(newx){
                                         tryCatch(
                                           fpca_BLUP(newx,
                                                     pca_res[[i]],
                                                     multiStart = multiStart,
                                                     scale = 'origin'),
                                           error = function(e){rep(NA, gridSize)}
                                         )
                                       })
                               },
                               simplify = FALSE)
        ls.res.esti[[idx.esti]] = do.call(cbind, test_fpcaBLUP)
      })
    }

    if(esti.name[idx.esti] == 'FPCA_MKE'){
      ls.t.esti[[idx.esti]] = system.time({
        # not exactly sure about the convexity of objective func, multistart anyway.
        if(!input.multiStart){
          # in this case, others are not using multiStart, need to calculate its own
          # determine initial values
          tm.Range = matrix(apply(pca_res[[1]]$xiEst, 2, range), nrow = 2)
          # set as so since mean(xiEst) = 0, choose the hyper-rectangle so as no to go outside.
          MKE.multiStart = c(min(9, length(pca_res[[1]]$lambda)*2 - 1), pmin(abs(tm.Range[2,]), abs(tm.Range[1,])) / 2)
        }else{
          # if all multiStart, steal
          MKE.multiStart = multiStart
        }
        test_fpcaMKE = sapply(1:nTrain,
                              function(i){
                                apply(test_obs[(nTest*(i-1)+1) : (nTest*i), , drop = FALSE], 1,
                                      function(newx){
                                        tryCatch(
                                          fpca_MKE(newx,
                                                   pca_res[[i]],
                                                   multiStart = MKE.multiStart,
                                                   scale = 'origin'),
                                          error = function(e){rep(NA, gridSize)}
                                        )
                                      })
                              },
                              simplify = FALSE)
        ls.res.esti[[idx.esti]] = do.call(cbind, test_fpcaMKE)
      })
    }

    if(esti.name[idx.esti] == 'MLE'){
      # calculating MLE with known family
      dimPar = ncol(train_par[[1]])
      # now if multiStart , its range should be determined by the actual parameter
      if(!input.multiStart){
        MLE.multiStart = c(1,0)
        init_par = colMeans(matrix(train_par[[1]], ncol = dimPar))
      }else{
        # determine initial values
        tm.Range = matrix(apply(train_par[[1]], 2, range), ncol = dimPar)
        MLE.multiStart = c(min(9, dimPar * 2 - 1), (tm.Range[2,] - tm.Range[1,])/3)
        init_par = as.numeric(colMeans(tm.Range))
      }
      ls.t.esti[[idx.esti]] = system.time({
        ls.res.esti[[idx.esti]] = apply(test_obs, 1,
                                        function(newx){
                                          tryCatch(
                                            calMLE(newx, denFam = density_family, initPar = init_par,
                                                   grid = grid, scale = 'origin',
                                                   multiStart = MLE.multiStart),
                                            error = function(e){rep(NA, gridSize)}
                                          )
                                          # tryCatch(
                                          #   tmpFunc_MLE(newx),
                                          #   error = function(e){rep(NA, gridSize)}
                                          # )
                                        })
      })
      # tmpFunc_MLE = function(newx){
      #   opt_res = smartOptim(init_par, function(par){-1*sum(log(pdf(newx, par)))},
      #                        method = 'Nelder-Mead', multiStart = multiStart)
      #   # if(opt_res$convergence!=0){
      #   #   warning(c('optim convergence code ', opt_res$convergence), ' while calculating MLE.')
      #   # }
      #   if(is.na(opt_res$value)){
      #     stop('NA as maximum.')
      #   }
      #   return(pdf(grid, opt_res$par))
      # }
      # if(dimPar == 1){
      #   tmpFunc_MLE = function(newx){
      #     opt_res = optim(rep(1, dimPar), function(par){-1*sum(log(pdf(newx, par)))}, method = 'Brent', lower = -1e+3, upper = 1e+3)
      #     if(opt_res$convergence!=0){
      #       warning(c('optim convergence code ', opt_res$convergence), ' while calculating MLE.')
      #     }
      #     if(is.na(opt_res$value)){
      #       stop('NA as maximum.')
      #     }
      #     return(pdf(grid, opt_res$par))
      #   }
      # }
    }

    if(esti.name[idx.esti] == 'emMixGauss'){
      ls.t.esti[[idx.esti]] = system.time({
        test_em = sapply(1:nTrain,
                         function(i){
                           apply(test_obs[(nTest*(i-1)+1) : (nTest*i), , drop = FALSE], 1,
                                 function(newx){
                                   tryCatch(
                                     emMixGauss(x = newx, k = 3, returnDensity = grid),
                                     error = function(e){rep(NA, gridSize)}
                                   )
                                 })
                         },
                         simplify = FALSE)
        ls.res.esti[[idx.esti]] = do.call(cbind, test_em)
      })
    }

    if(esti.name[idx.esti] == 'emTruncMixGauss'){
      ls.t.esti[[idx.esti]] = system.time({
        test_emT = sapply(1:nTrain,
                          function(i){
                            apply(test_obs[(nTest*(i-1)+1) : (nTest*i), , drop = FALSE], 1,
                                  function(newx){
                                    tryCatch(
                                      emTruncMixGauss(x = newx, k = 3, edge = density_family$range, returnDensity = grid),
                                      error = function(e){rep(NA, gridSize)}
                                    )
                                  })
                          },
                          simplify = FALSE)
        ls.res.esti[[idx.esti]] = do.call(cbind, test_emT)
      })
    }

    if(esti.name[idx.esti] == 'mleTruncMixGauss'){
      ls.t.esti[[idx.esti]] = system.time({
        test_mleTMG = sapply(1:nTrain,
                             function(i){
                               apply(test_obs[(nTest*(i-1)+1) : (nTest*i), , drop = FALSE], 1,
                                     function(newx){
                                       tryCatch(
                                         mleTruncMixGauss(x = newx, k = 3, edge = density_family$range, returnDensity = grid),
                                         error = function(e){rep(NA, gridSize)}
                                       )
                                     })
                             },
                             simplify = FALSE)
        ls.res.esti[[idx.esti]] = do.call(cbind, test_mleTMG)
      })
    }



  }

  # compute KDE in asked
  if(flg.kde){
    # In the following way errors could be inflated, especially when newx are highly concentrated.
    # say Beta(0.5, 100) at x = 0. But it seems ok, I guess. Most of the time the following is
    # doing better than the method I thought would be better.
    ls.t.esti[[idx.esti +1]] = system.time({
      test_kde = apply(test_obs, 1, function(newx){density(newx, bw = 'SJ',
                                                           from = min(range), to = max(range), n = gridSize)})
    })
    # I thought this would be better, but most the time not. Only better for highly concentrated data
    # t_kde = system.time({
    #   test_kde = apply(test_obs, 1, stats::density, bw = 'SJ')
    # })
  }


  ##### hopefully this is human readable
  if(!silent) print('Calculating error...')
  # err.name = c('ISE', 'KL', 'rev.KL', 'FisherRao', 'Wass..5', 'Wass.1', 'Wass.2')
  # # order of esti.name should be consistant with the test_fpca etc. in esti_pdf
  # esti.name = c('FPCA_MLE', 'FPCA_MAP', 'FPCA_BLUP', 'FPCA_MKE', 'MLE', 'KDE')

  # esti_pdf is a data.frame where each column is one estimated density,
  # columns ordered as esti.method --> trainCase --> testCase
  esti_pdf = data.frame(
    do.call(cbind, ls.res.esti)
  )
  # esti_pdf = data.frame(
  #                       matrix(c(test_fpcaMLE, test_fpcaMAP, test_fpcaBLUP, test_fpcaMKE, test_MLE), nrow = gridSize)
  #             )

  t_err = system.time({
    # tmErr will be a list of length = #of esti density where
    # one element in the list correspond to errors of the estimated density w.r.t. different types of error
    # the elements of the list are in the order of esti.method --> trainCase --> testCase
    # i.e., the first nTest elements are errors fitted by fpcaMLE with 1st training result.
    tmErr = ( purrr::map2(as.list(esti_pdf), rep(as.list(data.frame(test_true)), length(esti.name)),
                          function(esti,actual){
                            calErr(actual, esti, grid = grid, method = err.name)
                          }
    ))
  })

  # calculate error for KDE
  if(flg.kde){
    t_err = t_err +
      system.time({
        tmErr = c(tmErr,
                  purrr::map2(test_kde, as.list(data.frame(test_true)),
                              function(esti, actual){
                                calErr(actual, esti$y, grid = esti$x, grid2 = grid, method = err.name)
                              })
        )
      })
    # and put names back
    esti.name = c(esti.name, 'KDE')
  }

  # translate into long format
  Err = data.frame(
    # id of train result
    idx.train = rep(rep(1:nTrain, each = nTest), times = length(esti.name)),
    # id of test case
    idx.test = rep(rep(1:nTest, each = length(err.name)), times = nTrain * length(esti.name)),
    # id of estimation method
    method.esti = rep(esti.name, each = length(err.name) * nTest * nTrain),
    # id of err measurement
    type.err = rep(err.name, times = nTest * nTrain * length(esti.name)),
    # actual errors
    err = as.numeric(unlist(tmErr))
  )

  # Err = list()
  # lst_mtd_name = c('ISE', 'KL', 'rev.KL', 'FisherRao')
  # for(i in 1:4){
  #   tm = as.numeric(tmErr[i,])
  #   m.tm = matrix(colMeans(matrix(tm, nrow = nTest)), nrow = nTrain)
  #   tm = matrix(tm, nrow = nTest * nTrain)
  #   colnames(m.tm) = c('FPCA_MLE', 'FPCA_MAP', 'FPCA_BLUP', 'MLE', 'KDE')
  #   colnames(tm) = c('FPCA_MLE', 'FPCA_MAP', 'FPCA_BLUP', 'MLE', 'KDE')
  #   m.tm = data.frame(m.tm)
  #   tm = data.frame(tm)
  #   Err[[ lst_mtd_name[i] ]] = tm
  #   Err[[ paste0('M', lst_mtd_name[i]) ]] = m.tm
  # }

  # t_err = system.time({
  #   ISE = as.numeric( purrr::map2(esti_pdf, rep(as.list(data.frame(test_true)), 5),
  #      function(x,y){
  #        calISE(x, y, grid = grid)
  #      }
  #   ))
  # })
  # # formating output
  # MISE = matrix(colMeans(matrix(ISE, nrow = nTest)), nrow = nTrain)
  # ISE = matrix(ISE, nrow = nTest * nTrain)
  # colnames(MISE) = c('FPCA_MLE', 'FPCA_MAP', 'FPCA_BLUP', 'MLE', 'KDE')
  # colnames(ISE) = c('FPCA_MLE', 'FPCA_MAP', 'FPCA_BLUP', 'MLE', 'KDE')
  # MISE = data.frame(MISE)
  # ISE = data.frame(ISE)

  # Another way doing this, same spead.
  # {  esti_pdf = matrix(c(test_fpcaMLE, test_fpcaMAP, test_fpcaBLUP, test_MLE, test_kde), nrow = gridSize)
  #   t_err = system.time({
  #     ISE2 =
  #       sapply(1:5, function(i){
  #         apply(test_true - esti_pdf[,((i-1)*(nTrain*nTest)+1):(i*(nTrain*nTest)), drop = FALSE], 2,
  #                                   function(x){
  #                                     dotL2(x, x, grid = grid)
  #                                   }
  #         )
  #       })
  #   })



  # Print timing
  duration = do.call(rbind, ls.t.esti)
  duration = rbind(t_train, duration, t_err)
  row.names(duration) = c('Training', esti.name, 'Err_Cal')
  if(!silent) print(duration)

  if(returnFPCA){
    return(list(
      FPCA = pca_res, Err = Err,
      train_par = train_par,
      test_obs = test_obs,
      test_par = test_par,
      esti_method = esti.name,
      # note that KDE esti. is not returned
      esti_pdf = esti_pdf,
      duration = duration,
      elapsed = proc.time() - t_ttl
    ))
  }else{
    simpleFPCA = list()
    for(i in 1:length(pca_res)){
      simpleFPCA[[i]] = list(
        lambda = pca_res[[i]]$lambda,
        phi = pca_res[[i]]$phi,
        cumFVE = pca_res[[i]]$cumFVE
      )
    }
    return(list(
      simpleFPCA = simpleFPCA,
      Err = Err,
      train_par = train_par,
      test_obs = test_obs,
      test_par = test_par,
      esti_method = esti.name,
      # note that KDE esti. is not returned
      esti_pdf = esti_pdf,
      grid = grid,
      duration = duration,
      elapsed = proc.time() - t_ttl
    ))
  }
}

# Some markup during composing:
# tst = matrix(rnorm(4),2,2)
# split(tst, row(tst)) %>% microbenchmark(times = 10000)
# split(tst, c(row(tst))) %>% microbenchmark(times = 10000)
# as.list(data.frame(t(tst))) %>% microbenchmark(times = 10000)
#
# tst = matrix(rnorm(6), 2,3)
# apply(tst, 1, function(x){x+1})
# tst
# tst1 = matrix(rnorm(6),2,3)
# tst2 = matrix(rnorm(8),2,4)
# map2(as.list(data.frame(t(tst1))), as.list(data.frame(t(tst2))), function(x,y){range(x)+min(y)})
#
# map(as.list(data.frame(t(tst1))), max) %>% microbenchmark
# apply(tst1, 1, max) %>% microbenchmark


# Examples
# nrml_family = list(pdf = function(x, par){dnorm(x, par[1], par[2])},
#                    rpdf = function(n, par){rnorm(n, par[1], par[2])},
#                    gen_par = function(n){matrix(c(runif(n, -5,5),
#                                                   runif(n,1,5)), ncol = 2)}
#                    )
# func_ISE = method_eval(nrml_family, trainSize = 1000, nTest = 1000, testSize = 10, range = c(-20,20), gridSize = 512)
# par(mfrow = c(1,2))
# boxplot(func_ISE)
#

# beta_family = list(pdf = function(x, par){dbeta(x, exp(par[1]), exp(par[2]))},
#                    rpdf = function(n, par){rbeta(n,exp(par[1]), exp(par[2]))},
#                    gen_par = function(n){return(log(matrix(c(runif(n, 0.5,5),
#                                                   runif(n,0.5,5)), ncol = 2)))}
# )
# set.seed(100)
# beta_res = method_eval(beta_family, trainSize = 1000, nTest = 10000, testSize = 10, range = c(1e-6,1-1e-6), gridSize = 512, refPoint = 0.5)
# beepr::beep(4)
# boxplot(beta_res$Err$ISE, ylim = c(0, 2))
# summary(beta_res$Err$ISE)
# fpca_MLE(rep(1e-7, 10), beta_res$FPCA)
#
#
# beta

#
# gma_family = list(pdf = function(x, par){dgamma(x, par[1], par[2])},
#                    rpdf = function(n, par){rgamma(n, par[1], par[2])},
#                    gen_par = function(n){matrix(c(runif(n, 0.5,5),
#                                                   runif(n,0.5,5)), ncol = 2)}
# )
# gma_ISE = method_eval(gma_family, trainSize = 1000, nTest = 1000, testSize = 10, range = c(0.01,10), gridSize = 512, refPoint = 3)
# boxplot(gma_ISE, ylim = c(0, 2))
# summary(gma_ISE)

# problem:
# 1. convergence code





learnFamily = function(density_family, trainSize, maxK = 5, gridSize = 128, refPoint = 0, simple_result = TRUE){
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

