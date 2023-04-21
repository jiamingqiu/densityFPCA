## functions for fpca related estimates

fpcaEsti <- function(mat.obsv, fpca.res, esti.method, control = list()){
  # compute fpca based estimation of new.obsv given fpca.res with esti.method
  # args:
  # mat.obsv = a matrix or list where one row is one sample;
  # fpca.res = a result from fdapacd::FPCA;
  # esti.method = a character array of the method names, say
  #              "FPCA_BLUP", "FPCA_MLE", "FPCA_MAP", "FPCA_MKE";
  # control = a list of controling arguments:
  #   num.k = number of eigenfunc to use, missing -> use all, can be "AIC" as well;
  #   max.k = maximum number of eigenfunc to use, neglected unless num.k = 'AIC';
  #   init.par = initial value for optim, missing -> 0, neglected if num.k not number;
  #   method = method to use for optimization, missing -> DIRECT-P;
  #   nl.control = list of options for nloptr, c.f. nloptr::nl.opts, missing -> list();
  #   return.scale = what to return: "log", "origin", "parameter"(default), "optim";
  #   grid = grid to compute return when return.scale = "log"/"origin",
  #          missing -> fpca.res$workGrid;
  #   wass.p = (FPCA_MKE only) power of Wasserstein metric.
  #   extend.par.range, check, check.sample = passed to fpca2DenFam, use default
  #   if missing, set.fpcaEstFun_control may alter some default values.
  #
  # returns:
  # a list of
  # - a matrix where estimated parameters / (log-)densities as rows,
  #   or list of optim results;
  # - grid, NULL if not necessary;
  # - a data.frame for indexing obervations and etimating methods;
  # - list of error message, if any;
  # - time.
  # Note:
  #   num.k = "AIC" implemented for all method, but only justified for FPCA_MLE.

  # DEBUG
  # err = 'NA'
  # dummy <- genSimData(tNrml(6), 50, Inf, 5, 5)
  # train.sample <- dummy$train.sample
  # new.obsv <- dummy$test.sample
  # grid <- dummy$grid
  # esti.method <- c('FPCA_MLE', 'FPCA_MKE')
  # fpca.option <- list(
  #   error = FALSE,
  #   lean = TRUE,
  #   maxK = 5,
  #   plot = FALSE,
  #   useBinnedData = 'OFF'
  # )
  # proc.train <- preTrain(train.sample = train.sample, grid = grid)
  # proc.train$optns <- fpca.option
  # fpca.res <- do.call(fdapace::FPCA, proc.train)
  # esti.method <- fpcaEstiFun(esti.method, list(multiStart = c(1, 0), scale = 'origin'))
  # esti.method <- esti.method[[1]]

  # if (class(fpca.res) != 'FPCA') {
  #   stop('check input fpca.res.')
  # }

  # match estimation method names
  esti.method <- match.arg(
    esti.method,
    c("FPCA_BLUP", "FPCA_MAP", "FPCA_MKE", "FPCA_MLE", "FPCA_MOM"),
    several.ok = TRUE
  )
  # filling default and check input control
  control <- set.fpcaEstFun_control(environment())
  # if mat.obsv is not a list
  if(!is.list(mat.obsv) & !is.matrix(mat.obsv)){
    stop('check input type, must be list or matrix.')
  }
  if(is.matrix(mat.obsv)){
    mat.obsv <- mat2list(mat.obsv, by = 1)
  }

  # LEGACY if mat.obsv is not a matrix
  # if (!is.matrix(mat.obsv)) {
  #   mat.obsv <- matrix(mat.obsv, nrow = 1)
  # } END LEGACY

  # some preparation:
  return.grid.size <- length(control$grid)
  w.control <- control
  w.control$return.scale <- "optim"
  # error handling function
  err.f <- function(err) err
  # list of controls
  n.esti.method <- length(esti.method)
  w.control <- replicate(n.esti.method, w.control, simplify = FALSE)
  if('FPCA_BLUP' %in% esti.method){
    # prepare quantities used in BLUP estimator to save computation
    ls.cov <- list()
    if(control$num.k == 'AIC'){                  # if AIC
      for(nk in seq_len(control$max.k)){
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
            E_mCo = E_mCo,
            cov_mCo = cov_mCo,
            E_cond_cov_T = E_cond_cov_T
          )
        })
      }
      names(ls.cov) <- seq_len(control$max.k)
    }else{                                     # if specify num.k
      ls.cov[[1]] <- local({
        tm.cov <- condVarT(
          fpca.res = fpca.res,
          e.co = fpca.res$xiEst[, seq_len(control$num.k), drop = FALSE],
          return.all = TRUE
        )
        E_mCo <- colMeans(tm.cov$m.co)
        cov_mCo <- cov(tm.cov$m.co)

        # E(Var(T|natural par))
        E_cond_cov_T <- matrix(0, control$num.k, control$num.k)
        E_cond_cov_T[lower.tri(E_cond_cov_T, diag = T)] <- # must be lower.tri, R index row first
          colMeans(tm.cov$cond.cov)
        E_cond_cov_T <- E_cond_cov_T + t(E_cond_cov_T)
        diag(E_cond_cov_T) <- diag(E_cond_cov_T) / 2
        list(
          E_mCo = E_mCo,
          cov_mCo = cov_mCo,
          E_cond_cov_T = E_cond_cov_T
        )
      })
      names(ls.cov) <- control$num.k
    }
    w.control[[which(esti.method == 'FPCA_BLUP')]]$ls.cov <- ls.cov
  }
  # list of estimating functions
  ls.f <- list()
  for(i in seq_len(n.esti.method)){
    ls.f[[i]] <- get_fpca_esti_fun(
      esti.method = esti.method[i],
      fpca.res = fpca.res,
      control = w.control[[i]],
      err.f = err.f
    )
    # LEGACY
    # ls.f[[i]] <- switch(
    #   esti.method[i],
    #   "FPCA_BLUP" = function(obsv){
    #     tryCatch(fpca_BLUP(obsv, fpca.res, w.control), error = err.f)
    #   },
    #   "FPCA_MAP" = function(obsv){
    #     tryCatch(fpca_MAP(obsv, fpca.res, w.control), error = err.f)
    #   },
    #   "FPCA_MKE" = function(obsv){
    #     tryCatch(fpca_MKE(obsv, fpca.res, w.control), error = err.f)
    #   },
    #   "FPCA_MLE" = function(obsv){
    #     tryCatch(fpca_MLE(obsv, fpca.res, w.control), error = err.f)
    #   }
    # )
  }


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
        tm.res <- (ls.f[[idx.esti]])(mat.obsv[[i]])

        if (!inherits(tm.res, 'error')) {
          ls.opt.res[[idx.run]] <- tm.res
        } else {
          ls.opt.res[[idx.run]] <- NULL
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

  # now in format of 1 element in ls.opt.res = 1 opt.res, and in trunc of esti.method
  # indexing data.frame
  if(is.numeric(control$num.k)){
    lbl.num.k <- control$num.k
    lbl.aic <- NULL
  }else{
    lbl.num.k <- sapply(
      ls.opt.res,
      function(x){
        if(is.null(x))
          return(NA)
        return(x$num.k)
      }
    )
    lbl.aic <- sapply(
      ls.opt.res,
      function(x){
        if(is.null(x))
          return(NA)
        return(x$AIC)
      }
    )
  }
  df.idx <- data.frame(
    idx.obsv = rep(seq_along(mat.obsv), times = length(esti.method)),
    esti.method = rep(esti.method, each = length(mat.obsv))
  )
  df.idx$method.num.k <- as.character(control$num.k)
  df.idx$num.k <- lbl.num.k
  df.idx$AIC <- lbl.aic

  # get the time
  df.t <- do.call(rbind, ls.t)
  df.t <- rbind(df.t, colSums(df.t))
  df.t <- cbind(data.frame(step = c(esti.method, 'TOTAL')), df.t)

  # return
  if (control$return.scale == 'optim') {
    return(list(
      res = ls.opt.res,
      idx = df.idx,
      grid = NULL,
      error = ls.message,
      time = df.t
    ))
  }
  # the following return.scale need parameters
  mat.par <- matrix(
    NA,
    nrow = length(mat.obsv) * n.esti.method,
    ncol = control$max.k
  )
  for(i in seq_len(nrow(mat.par))){
    if (is.null(ls.opt.res[[i]])) {
      mat.par[i, ] <- rep(NA, control$max.k)
    } else {
      mat.par[i, seq_len(df.idx$num.k[i])] <- ls.opt.res[[i]]$par
    }
  }
  # trim unnecessary columns when num.k is specified.
  if (is.numeric(control$num.k)){
    mat.par <- mat.par[, seq_len(control$num.k), drop = FALSE]
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
  mat.pdf <- fpcaEsti2pdf(
    fpca.res, mat.par,
    num.k = df.idx$num.k,
    grid = control$grid
  )
  # LEGACY, equal num.k case
  # fpca.den.fam <- fpca2DenFam(
  #   fpca.res,
  #   list(num.k = control$num.k)
  # )
  # mat.pdf <- par2pdf(fpca.den.fam, mat.par, grid = control$grid)
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

fpcaEsti.f <- function(fpca.res, esti.method, control = list()){
  # get functions for fpca based density estimation given fpca.res with esti.method
  # args:
  # fpca.res = a result from fdapacd::FPCA;
  # esti.method = a character array of the method names, say
  #              "FPCA_BLUP", "FPCA_MLE", "FPCA_MAP", "FPCA_MKE";
  # control = a list of controling arguments:
  #   num.k = number of eigenfunc to use, missing -> use all, can be "AIC" as well;
  #   max.k = maximum number of eigenfunc to use, neglected unless num.k = 'AIC';
  #   init.par = initial value for optim, missing -> 0, neglected if num.k not number;
  #   method = method to use for optimization, missing -> DIRECT-P;
  #   nl.control = list of options for nloptr, c.f. nloptr::nl.opts, missing -> list();
  #   return.scale = what to return: "parameter"(default), "optim", "function";
  #   wass.p = (FPCA_MKE only) power of Wasserstein metric.
  #   extend.par.range, check, check.sample = passed to fpca2DenFam, use default
  #   if missing, set.fpcaEstFun_control may alter some default values.
  #   prior: (optional, internal use only) a list of condVarT results to save time.
  # returns:
  # a named list of estimating functions, taking one sample giving whatever
  # desired as specified in return.scale.
  #
  # Note:
  #   num.k = "AIC" implemented for all method, but only justified for FPCA_MLE.
  #   acceptable return.scale is slightly different to fpcaEsti.
  #   for example, when return.scale = "function", fpcaEsti.f returns f s.t.
  #   f(obsv) gives h, where h is a function s.t. h(x) is estimated density at x.
  #   and can handle multiple return.scale, in such case the f(obsv) will be a
  #   list.

  # DEBUG
  # err = 'NA'
  # dummy <- genSimData(tNrml(6), 50, Inf, 5, 5)
  # train.sample <- dummy$train.sample
  # new.obsv <- dummy$test.sample
  # grid <- dummy$grid
  # esti.method <- c('FPCA_MLE', 'FPCA_MKE')
  # fpca.option <- list(
  #   error = FALSE,
  #   lean = TRUE,
  #   maxK = 5,
  #   plot = FALSE,
  #   useBinnedData = 'OFF'
  # )
  # proc.train <- preTrain(train.sample = train.sample, grid = grid)
  # proc.train$optns <- fpca.option
  # fpca.res <- do.call(fdapace::FPCA, proc.train)
  # esti.method <- fpcaEstiFun(esti.method, list(multiStart = c(1, 0), scale = 'origin'))
  # esti.method <- esti.method[[1]]

  # if (class(fpca.res) != 'FPCA') {
  #   stop('check input fpca.res.')
  # }

  # match estimation method names
  esti.method <- match.arg(
    esti.method,
    c("FPCA_BLUP", "FPCA_MAP", "FPCA_MKE", "FPCA_MLE", "FPCA_MOM"),
    several.ok = TRUE
  )
  # filling default and check input control, some workaround
  tm.return.scale <- control$return.scale
  control$return.scale <- NULL
  control <- set.fpcaEstFun_control(environment())
  control$return.scale <- tm.return.scale
  control$return.scale <- match.arg(
    control$return.scale,
    c("parameter", "optim", "function"),
    several.ok = TRUE
  )
  # some preparation:
  w.control <- control
  w.control$return.scale <- "optim"
  # error handling function
  err.f <- function(err) err

  # list of controls
  n.esti.method <- length(esti.method)
  w.control <- replicate(n.esti.method, w.control, simplify = FALSE)

  if('FPCA_BLUP' %in% esti.method){
    if(!is.null(control$prior)){
      ls.cov <- control$prior
    }else{
      # prepare quantities used in BLUP estimator to save computation
      ls.cov <- list()
      if(control$num.k == 'AIC'){                  # if AIC
        for(nk in seq_len(control$max.k)){
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
              E_mCo = E_mCo,
              cov_mCo = cov_mCo,
              E_cond_cov_T = E_cond_cov_T
            )
          })
        }
        names(ls.cov) <- seq_len(control$max.k)
      }else{                                       # if specify num.k
        ls.cov[[1]] <- local({
          tm.cov <- condVarT(
            fpca.res = fpca.res,
            e.co = fpca.res$xiEst[, seq_len(control$num.k), drop = FALSE],
            return.all = TRUE
          )
          E_mCo <- colMeans(tm.cov$m.co)
          cov_mCo <- cov(tm.cov$m.co)

          # E(Var(T|natural par))
          E_cond_cov_T <- matrix(0, control$num.k, control$num.k)
          E_cond_cov_T[lower.tri(E_cond_cov_T, diag = T)] <- # must be lower.tri, R index row first
            colMeans(tm.cov$cond.cov)
          E_cond_cov_T <- E_cond_cov_T + t(E_cond_cov_T)
          diag(E_cond_cov_T) <- diag(E_cond_cov_T) / 2
          list(
            E_mCo = E_mCo,
            cov_mCo = cov_mCo,
            E_cond_cov_T = E_cond_cov_T
          )
        })
        names(ls.cov) <- control$num.k
      }
      ls.cov <- ls.cov
    }
    w.control[[which(esti.method == 'FPCA_BLUP')]]$ls.cov <- ls.cov
  }

  # list of estimating functions.
  ls.f <- list()
  for(i in seq_along(esti.method)){
    ls.f[[i]] <- get_fpca_esti_fun(
      esti.method = esti.method[i],
      fpca.res = fpca.res,
      control = w.control[[i]],
      err.f = err.f
    )
  }
  # browser()
  res <- list()
  for(i in seq_along(esti.method)){
    res[[i]] <- local({
      # viva la lazy eval,
      # c.f. https://stackoverflow.com/questions/31556088/r-defining-functions-within-a-loop
      i <- i
      function(obsv){
        opt.res <- ls.f[[i]](obsv)
        den.fam <- fpca2DenFam(fpca.res, list(num.k = length(opt.res$par)))
        res <- list()
        for(nm.return in as.character(control$return.scale)){
          res[[nm.return]] <- switch(
            nm.return,
            "parameter" = opt.res$par,
            "optim" = opt.res,
            "function" = function(x){
                den.fam$pdf(x, opt.res$par)
            }
          )
        }
        if(length(res) == 1)
          return(res[[1]])
        else
          return(res)
      }
      # switch(
      #   control$return.scale,
      #   "parameter" = function(obsv){
      #     opt.res <- ls.f[[i]](obsv)
      #     return(opt.res$par)
      #   },
      #   "optim" = ls.f[[i]],
      #   "function" = function(obsv){
      #     opt.res <- ls.f[[i]](obsv)
      #     den.fam <- fpca2DenFam(fpca.res, list(num.k = length(opt.res$par)))
      #     return(function(x){
      #       den.fam$pdf(x, opt.res$par)
      #     })
      #   }
      # )
    })
  }
  return(setNames(res, esti.method))
}

# no longer needed after unifying call syntax for individual methods.
# fpcaEstiFun <- function(esti.method, esti.option){
#   # get the actual fpca based estimating function from its name
#   # args:
#   # esti.method: the name (character) of estimating method.
#   # esti.option: a list of options to pass to actual estimating function.
#   # return:
#   # the desired estimating function.
#
#   ls.esti.fun <- kListFpcaMethod
#
#   if (!all(esti.method %in% names(ls.esti.fun))) {
#     stop('unknown method.')
#   }
#
#   if (length(esti.method) == 1) {
#     res.fun <- function(newx, fpca.res){
#       (ls.esti.fun[[esti.method]])(
#         newx = newx, fpca_res = fpca.res,
#         scale = esti.option$scale, multiStart = esti.option$multiStart
#       )
#     }
#   } else {
#     # TODO: this part seems to be broken, idk why.
#     res.fun <- list()
#     for (i in seq_len(length(esti.method))) {
#       tm.name <- esti.method[i]
#       res.fun[[ tm.name ]] <- function(newx, fpca.res){
#         #	  cat(paste0(tm.name, '\n'))
#         (ls.esti.fun[[ tm.name ]])(
#           newx = newx, fpca_res = fpca.res,
#           scale = esti.option$scale, multiStart = esti.option$multiStart
#         )
#       }
#     }
#   }
#   return(res.fun)
#
# }


### individual methods ========================================================

# some tools ==================================================================

set.fpcaEstFun_control <- function(envir){
  # setting control list for individual methods below
  # args:
  # envir = the current within function environment with control, fpca.res.
  # return:
  # a set control list.

  control <- envir$control

  ##### num.k related
  if (is.null(control$num.k)) {
    control$num.k <- ncol(envir$fpca.res$phi)
  }
  if (!is.numeric(control$num.k)) {
    control$num.k <- match.arg(control$num.k, 'AIC')
    if (is.null(control$max.k)) {
      control$max.k <- min(
        which(envir$fpca.res$cumFVE >= 99.999),
        ncol(envir$fpca.res$phi)
      )
    }
    if (control$max.k > ncol(envir$fpca.res$phi)) {
      warning('Insufficient number of eigenfucntions, max.k reset.')
      control$max.k <- ncol(envir$fpca.res$phi)
    }
  } else {
    if (control$num.k > ncol(envir$fpca.res$phi)) {
      warning('Insufficient number of eigenfucntions, num.k reset.')
      control$num.k <- ncol(envir$fpca.res$phi)
    }
    if (is.null(control$init.par)) {
      control$init.par <- rep(0, control$num.k)
    } else {
      if (length(control$init.par) != control$num.k) {
        stop('Check length of init.par.')
      }
    }
    # max.k will be omitted when specifying num.k, for integrity only.
    control$max.k <- ncol(envir$fpca.res$phi)
  }

  ##### optimization related
  if (is.null(control$method)) {
    control$method <- "DIRECT-P"
  }
  if (is.null(control$nl.control)) {
    control$nl.control <- list()
  }

  ##### return related
  control$return.scale <- match.arg(
    control$return.scale,
    c("parameter", "log", "origin", "optim")
  )
  if (control$return.scale %in% c("log", "origin")) {
    if (is.null(control$grid)) {
      control$grid <- envir$fpca.res$workGrid
    }
    if (any(!is.finite(control$grid))) {
      stop('Non finite grid value.')
    }
  }

  ##### used only in fpca_MKE
  if (is.null(control$wass.p)) {
    control$wass.p <- 2
  } else {
    if (control$wass.p < 1) {
      stop(paste('check wass.p, input is', control$wass.p))
    }
  }

  ##### used only for fpca2DenFam
  # args: check, check.sample  (TBD)
  if (is.null(control$extend.par.range))
    control$extend.par.range <- 2

  return(control)

  # LEGACY
  # cran check will complain about using global variable here. CAREFUL.
  # so use other way
  # w.f <- function(){
  #   if (is.null(control$num.k)) {
  #     control$num.k <- length(fpca.res$lambda)
  #   }
  #   if (control$num.k > length(fpca.res$lambda)) {
  #     warning('Insufficient number of eigenfucntions, num.k reset.')
  #     control$num.k <- length(fpca.res$lambda)
  #   }
  #   if (is.null(control$init.par)) {
  #     control$init.par <- rep(0, control$num.k)
  #   }
  #   if (is.null(control$method)) {
  #     control$method <- "DIRECT-P"
  #   }
  #   if (is.null(control$nl.control)) {
  #     control$nl.control <- list()
  #   }
  #   control$return.scale <- match.arg(
  #     control$return.scale,
  #     c("parameter", "log", "origin", "optim")
  #   )
  #   if (control$return.scale %in% c("log", "origin")) {
  #     if (is.null(control$grid)) {
  #       control$grid <- fpca.res$workGrid
  #     }
  #   }
  #   return(control)
  # }
  # environment(w.f) <- envir
  # return(w.f())
}

get_fpca_esti_fun <- function(esti.method, fpca.res, control, err.f){
  # get the actual fpca based estimating function from its name and controls,
  # used internally in fpcaEsti.
  # args:
  # esti.method: the name (character) of estimating method
  # control: a list of options to pass to actual estimating function,
  #          return.scale will be overwritten to 'optim'
  # err.f: a function handling errors, default(missing) is just returning the err
  # return:
  # the desired estimating function f(obsv)

  # DEBUG
  # esti.method <- 'FPCA_MLE'
  # control <- list(num.k = 'AIC', max.k = 4)
  # err.f <- function(err) err
  # wgrid <- seq(-3, 3, length.out = 512)
  # set.seed(1)
  # fpca.res <- list(
  #   mu = rep(0, 512),
  #   phi = matrix(
  #     c(wgrid, wgrid^2, rnorm(512, sd = 0.5), rnorm(512, sd = 0.1)),
  #     ncol = 4
  #     ),
  #   xiEst = matrix(rnorm(512 * 4), ncol = 4),
  #   workGrid = wgrid
  # )
  # obsv <- rnorm(25)

  if(missing(err.f))
    err.f <- function(err) err

  esti.method <- match.arg(
    esti.method,
    c("FPCA_BLUP", "FPCA_MAP", "FPCA_MKE", "FPCA_MLE", "FPCA_MOM"),
    several.ok = FALSE
  )
  w.control <- control
  w.control$return.scale <- 'optim'
  # if specifying num.k
  if(is.numeric(control$num.k)){
    res.f <- switch(
      esti.method,
      "FPCA_BLUP" = function(obsv){
        tryCatch(fpca_BLUP(obsv, fpca.res, w.control), error = err.f)
      },
      "FPCA_MAP" = function(obsv){
        tryCatch(fpca_MAP(obsv, fpca.res, w.control), error = err.f)
      },
      "FPCA_MKE" = function(obsv){
        tryCatch(fpca_MKE(obsv, fpca.res, w.control), error = err.f)
      },
      "FPCA_MLE" = function(obsv){
        tryCatch(fpca_MLE(obsv, fpca.res, w.control), error = err.f)
      },
      "FPCA_MOM" = function(obsv){
        tryCatch(fpca_MOM(obsv, fpca.res, w.control), error = err.f)
      }
    )
    return(res.f)
  }

  # if num.k = 'AIC'
  control$num.k <- match.arg(control$num.k, 'AIC')
  if(is.null(control$max.k))
    control$max.k <- ncol(fpca.res$phi)
  w.control <- control
  w.control$return.scale <- 'optim'
  tm.f <- switch(
    esti.method,
    "FPCA_BLUP" = function(obsv, num.k){
      w.control$num.k <- num.k
      tryCatch(fpca_BLUP(obsv, fpca.res, w.control), error = err.f)
    },
    "FPCA_MAP" = function(obsv, num.k){
      w.control$num.k <- num.k
      tryCatch(fpca_MAP(obsv, fpca.res, w.control), error = err.f)
    },
    "FPCA_MKE" = function(obsv, num.k){
      w.control$num.k <- num.k
      tryCatch(fpca_MKE(obsv, fpca.res, w.control), error = err.f)
    },
    "FPCA_MLE" = function(obsv, num.k){
      w.control$num.k <- num.k
      tryCatch(fpca_MLE(obsv, fpca.res, w.control), error = err.f)
    },
    "FPCA_MOM" = function(obsv, num.k){
      w.control$num.k <- num.k
      tryCatch(fpca_MOM(obsv, fpca.res, w.control), error = err.f)
    }
  )
  w.tm.f <- function(obsv, num.k){
    tm.res <- tm.f(obsv, num.k)
    if(inherits(tm.res, 'error'))
      return(tm.res)
    tm.res$AIC <- 2 * num.k + 2 *
      do.call(
        getNegLogll(
          obsv, fpca2DenFam(fpca.res, list(num.k = num.k))
        ),
        list(par = as.numeric(tm.res$par))
      )
    return(tm.res)
  }
  res.f <- downHill.f(
    fun = w.tm.f, over = list(num.k = seq_len(control$max.k)),
    what = 'AIC'
  )
  # res.f(obsv)
  return(res.f)
}

fpca_MLE = function(newx, fpca.res, control = list()){
  # MLE of newx given FPCA result.
  # args:
  # newx = array of new observations;
  # fpca.res = what fdapace::FPCA returns
  # control = a list of controling arguments:
  #   num.k = number of eigenfunc to use, missing -> use all;
  #   init.par = initial value for optim, missing -> 0;
  #   method = method to use for optimization, missing -> DIRECT-P;
  #   nl.control = list of options for nloptr, c.f. nloptr::nl.opts, missing -> list();
  #   return.scale = what to return: "log", "origin", "parameter"(default), "optim";
  #   grid = grid to compute return when return.scale = "log"/"origin",
  #          missing -> fpca.res$workGrid.
  #   extend.par.range, check, check.sample = passed to fpca2DenFam, use default
  #   if missing, set.fpcaEstFun_control may alter some default values.
  # returns:
  # either the MLE or the corresponding log-density / density on grid,
  # or the optimization result, for debugging.

  # filling default and check input control
  control <- set.fpcaEstFun_control(environment())
  # if (is.null(control$num.k)) {
  #   control$num.k <- length(fpca.res$lambda)
  # }
  # if (control$num.k > length(fpca.res$lambda)) {
  #   warning('Insufficient number of eigenfucntions, num.k reset.')
  #   control$num.k <- length(fpca.res$lambda)
  # }
  # if (is.null(control$init.par)) {
  #   control$init.par <- rep(0, control$num.k)
  # }
  # if (is.null(control$method)) {
  #   control$method <- "DIRECT-P"
  # }
  # if (is.null(control$nl.control)) {
  #   control$nl.control <- list()
  # }
  # control$return.scale <- match.arg(
  #   control$return.scale,
  #   c("parameter", "log", "origin", "optim")
  # )
  # if (control$return.scale %in% c("log", "origin")) {
  #   if (is.null(control$grid)) {
  #     control$grid <- fpca.res$workGrid
  #   }
  # }

  # getting a dummy density family.
  den.fam <- fpca2DenFam(fpca.res, control = control)

  ls.objf <- fpca2ObjFun(newx, fpca.res, control = list(num.k = control$num.k))
  opt.res <- safeOptim(
    x0 = control$init.par,
    fn = ls.objf$pNegLogll,
    gr = ls.objf$pNegLogll.gr,
    lower = den.fam$par.range[1, ],
    upper = den.fam$par.range[2, ],
    method = control$method, backup.method = "DIRECT",
    nl.control = control$nl.control
  )

  if (control$return.scale == 'optim') {
    return(opt.res)
  }
  if (control$return.scale == 'log') {
    return(log(den.fam$pdf(control$grid, opt.res$par)))
  }
  if (control$return.scale == 'origin') {
    return(den.fam$pdf(control$grid, opt.res$par))
  }
  if (control$return.scale == 'parameter') {
    return(opt.res$par)
  }

}

fpca_MAP = function(newx, fpca.res, control = list()){
  # MAP of newx given FPCA result.
  # args:
  # newx = array of new observations;
  # fpca.res = what fdapace::FPCA returns
  # control = a list of controling arguments:
  #   num.k = number of eigenfunc to use, missing -> use all;
  #   init.par = initial value for optim, missing -> 0;
  #   method = method to use for optimization, missing -> DIRECT-P;
  #   nl.control = list of options for nloptr, c.f. nloptr::nl.opts, missing -> list();
  #   return.scale = what to return: "log", "origin", "parameter"(default), "optim";
  #   grid = grid to compute return when return.scale = "log"/"origin",
  #          missing -> fpca.res$workGrid.
  #   extend.par.range, check, check.sample = passed to fpca2DenFam, use default
  #   if missing, set.fpcaEstFun_control may alter some default values.
  # returns:
  # either the MAP or the corresponding log-density / density on grid,
  # or the optimization result, for debugging.

  # filling default and check input control
  control <- set.fpcaEstFun_control(environment())

  # getting a dummy density family.
  den.fam <- fpca2DenFam(fpca.res, control)
  # getting partial neglogll.
  ls.pNegLogll <- fpca2ObjFun(newx, fpca.res, list(num.k = control$num.k))

  # parametric normal prior (its log and gradient of its log)
  cov.prior <- cov(fpca.res$xiEst[, seq_len(control$num.k), drop = FALSE])
  inv.cov <- MASS::ginv(cov.prior)
  logPrior <- function(x) -0.5 * quadForm(inv.cov, x)
  logPrior.gr <- function(x) -1 * crossprod(inv.cov, x)

  ls.objf <- list(
    objf = function(par) ls.pNegLogll$pNegLogll(par) - logPrior(par),
    objf.gr = function(par) ls.pNegLogll$pNegLogll.gr(par) - logPrior.gr(par)
  )

  opt.res <- safeOptim(
    x0 = control$init.par,
    fn = ls.objf$objf,
    gr = ls.objf$objf.gr,
    lower = den.fam$par.range[1, ],
    upper = den.fam$par.range[2, ],
    method = control$method, backup.method = "DIRECT",
    nl.control = control$nl.control
  )

  if (control$return.scale == 'optim') {
    return(opt.res)
  }
  if (control$return.scale == 'log') {
    return(log(den.fam$pdf(control$grid, opt.res$par)))
  }
  if (control$return.scale == 'origin') {
    return(den.fam$pdf(control$grid, opt.res$par))
  }
  if (control$return.scale == 'parameter') {
    return(opt.res$par)
  }
}

# getLogll_optimOnly = function(fpca_res){
# # build a loglikelihood function given FPCA results,
# # but ONLY the part that is proportional to parameter.
#   return(
#     function(x, par){
#     # x:  input data, length = n, from n iid observation
#     # par: parameters
#     n = length(x)
#     nrml_const = log( dotL2(exp(fpca_res$mu), exp(fpca_res$phi %*% par), grid = fpca_res$workGrid) )
#     # in approx, rule = 2 to handle out of range new observation
#     mat_tm = apply(fpca_res$phi, 2, function(input){approx(x = fpca_res$workGrid, y = input, xout = x, rule = 2)$y})
#     return(sum(mat_tm %*% par) - n*nrml_const)
#   })
# }

m2eCoord = function(m.co, fpca.res, control = list()){
  # m-coordinate -> nature coordinate (e-coord)
  # args:
  # mco = m-coord, an array;
  # fpca.res = fdapace::FPCA result
  # control = a list of controling arguments:
  #   num.k = number of eigenfunc to use, missing -> use all;
  #   init.par = initial value for optim, missing -> 0;
  #   method = method to use for optimization, missing -> DIRECT-P;
  #   nl.control = list of options for nloptr, c.f. nloptr::nl.opts, missing -> list();
  #   return.scale = what to return: "log", "origin", "parameter"(default), "optim";
  #   grid = grid to compute return when return.scale = "log"/"origin",
  #          missing -> fpca.res$workGrid;
  #   extend.par.range, check, check.sample = passed to fpca2DenFam, use default
  #   if missing, set.fpcaEstFun_control may alter some default values.
  # returns:
  # corresponding e-coord, or optim result, or (log-)density on grid.
  # details:
  # by Legendre transformation, namely
  # maximizing m.co %*% par - nrml_const(par) for par to get e.co,
  # default initial value of optim set to be the column means of fpca_res$xiEst
  # which is of course 0.

  if (!is.matrix(m.co)) {
    m.co <- matrix(m.co, nrow = 1)
  }
  # filling default and check input control
  control <- set.fpcaEstFun_control(environment())
  # getting a dummy density family.
  den.fam <- fpca2DenFam(fpca.res, control)

  # the normalizing const
  local.nrmlConst <- function(e.co) nrmlConst(e.co, fpca.res, control$num.k)
  local.nrmlConst.gr <- function(e.co) nloptr::nl.grad(e.co, local.nrmlConst)

  objf <- function(e.co){
    return(-1 * sum(e.co * m.co) + local.nrmlConst(e.co))
  }
  objf.gr <- function(e.co){
    return(-1 * m.co + local.nrmlConst.gr(e.co))
  }

  opt.res <- safeOptim(
    x0 = control$init.par,
    fn = objf,
    gr = objf.gr,
    lower = den.fam$par.range[1, ],
    upper = den.fam$par.range[2, ],
    method = control$method, backup.method = "DIRECT",
    nl.control = control$nl.control
  )

  if (control$return.scale == 'optim') {
    return(opt.res)
  }
  if (control$return.scale == 'log') {
    return(log(den.fam$pdf(control$grid, opt.res$par)))
  }
  if (control$return.scale == 'origin') {
    return(den.fam$pdf(control$grid, opt.res$par))
  }
  if (control$return.scale == 'parameter') {
    return(opt.res$par)
  }
}


e2mCoord <- function(e.co, fpca.res, nrml.const, control = list()){
  # m-coordinate <- nature coordinate (e-coord)
  # args:
  # e.co = an array of natrue coordinate;
  # fpca.res = fdapace::FPCA result;
  # nrml.const = normalizing constant, provide s.t. no repeat computation;
  # control = a list of controling arguments:
  #   num.k = number of eigenfunc to use, missing -> use all.
  # returns:
  # corresponding m-coord.
  # details:
  # just integrate, since
  # m-coord = E(sufficient T(X)) with X ~ f(.|e-coord)

  if (is.null(control$num.k)) {
    control$num.k <- ncol(fpca.res$phi)
  }

  if(missing(nrml.const)){
    nrml.const = nrmlConst(e.co, fpca.res, control$num.k)
  }

  # get pdf on grid
  vec.pdf <- exp(
    fpca.res$mu +
      fpca.res$phi[, seq_len(control$num.k), drop = FALSE] %*% e.co -
      nrml.const
  )

  # now doing integral
  m.co <- rep(0, control$num.k)
  for(i in seq_len(length(m.co))){
    m.co[i] <- dotL2(
      fpca.res$phi[, i],
      vec.pdf,
      grid = fpca.res$workGrid
    )
  }
  return(m.co)
}

condVarT <- function(fpca.res, e.co, m.co, nrml.const, return.all = FALSE){
  # Compute conditional covariance matrix of T(X)|natural_parameter(e.co)
  # where pdf of X is defined by the fpca.res using length(e.co) many eigenfunc.
  # args:
  #   fpca.res = fdapace::FPCA result;
  #   e.co = natural parameter
  #   m.co & nrml.const = to save computation, see detail;
  #   return.all = should nrml.const and m.co be returned.
  # returns:
  #   a matrix where one row is one conditional cov matrix;
  #   or a list of m.co, nrml.const, cond.cov
  # Details:
  #   namely integrate(T_1 * T_2  * pdf)
  #   T here is phi, eigen functions.
  #   returns a matrix of (1+selectK)*selectK/2 \times number of trajectories in fpca_res input
  #   each row is the conditional covariance matrix of T,
  #   corresponding to location (1,1), (1,2), (1,3),... in the cov matrix
  #   if m.co (i.e., ET) and nrml.const of the density function are provided,
  #   then will not recalculated again.
  #   if m.co and nrml.const are both missing, will compute at e.co = fpca.res$xiEst.
  #   e.co or m.co must be a matrix with ncol = length(fpca.res$lambda), otherwise cut.
  #   Length(nrml.const) = nrow(e.co) or nrow(m.co)


  if (missing(e.co)) {
    e.co = fpca.res$xiEst
  }
  if (!is.matrix(e.co)) {
    e.co <- matrix(e.co, nrow = 1)
  }
  nk = ncol(e.co)


  # (willbe) matrix of densities on grid, one row correspond to one e.co
  log.pdf <- matrix(0, nrow = nrow(e.co), ncol = length(fpca.res$workGrid))
  for(i in seq_len(nrow(log.pdf))){
    log.pdf[i, ] <-
      fpca.res$phi[, seq_len(nk), drop = FALSE] %*% e.co[i, ]
  }

  # calculate normalizing constant for all pdf corresponding to e.co
  if (missing(nrml.const)) {
    nrml.const <- rep(0, nrow(e.co))
    for(i in seq_len(length(nrml.const))){
      nrml.const[i] <-
        log(dotL2(
          exp(fpca.res$mu),
          exp(log.pdf[i, ]),
          grid = fpca.res$workGrid
        ))
    }
    # apply(e.co, 1, function(par) nrmlConst(par, fpca.res, nk))
  }

  # now log.pdf is what it claims to be
  for(i in seq_len(nrow(log.pdf))){
    log.pdf[i, ] <- exp(fpca.res$mu + log.pdf[i, ] - nrml.const[i])
  }

  # computing m-coord if not provided.
  if (missing(m.co)) {
    m.co = matrix(0, nrow = nrow(e.co), ncol = nk)
    for(i in seq_len(nrow(m.co))){
      for(j in seq_len(ncol(m.co))){
        m.co[i, j] <- dotL2(
          fpca.res$phi[, j],
          log.pdf[i, ],
          grid = fpca.res$workGrid
        )
      }
      # m.co[i, ] <- e2mCoord(e.co, fpca.res, nrml.const[i], list(num.k = nk))
    }
  }

  # get index for combination, the indices are in lexical order
  idx.mat <- arrangements::combinations(nk, 2, replace = T)
  # matrix whose columns are T_i * T_j
  prod.mat <- matrix(
    fpca.res$phi[, idx.mat[,1]] * fpca.res$phi[, idx.mat[,2]],
    ncol = nk * (nk+1) / 2
  )
  # matrix whose columns are E(T_i | eCo_k) * E(T_j | eCo_k) at kth row
  mean.mat <- matrix(
    m.co[, idx.mat[,1]] * m.co[, idx.mat[,2]],
    ncol = nk * (nk+1) / 2
  )
  # matrix whose columns are E(T_i * T_j | eCo_k) at kth row
  cond.cov <- matrix(0, nrow = nrow(e.co), ncol = nk * (nk + 1) / 2)
  for(i in seq_len(nrow(cond.cov))){
    for(j in seq_len(ncol(cond.cov))){
      cond.cov[i, j] <-
      dotL2(
        prod.mat[, j],
        log.pdf[i, ],
        grid = fpca.res$workGrid
      )
    }
  }
  # now cond.cov is a matrix with nrow = nrow(eCo), ncol = nk * (nk+1) / 2, to get condVar
  # substract crossprod of m.co
  cond.cov = cond.cov - mean.mat
  if (!return.all) {
    return(cond.cov)
  } else {
    return(list(
      m.co = m.co,
      nrml.const = nrml.const,
      cond.cov = cond.cov
    ))
  }

}

fpca_BLUP = function(newx, fpca.res, control = list()){
  # FPCA-BLUP method of newx given FPCA result.
  # args:
  # newx = array of new observations;
  # fpca.res = what fdapace::FPCA returns
  # control = a list of controling arguments:
  #   num.k = number of eigenfunc to use, missing -> use all;
  #   init.par = initial value for optim, missing -> 0;
  #   method = method to use for optimization, missing -> DIRECT-P;
  #   nl.control = list of options for nloptr, c.f. nloptr::nl.opts, missing -> list();
  #   return.scale = what to return: "log", "origin", "parameter"(default), "optim";
  #   grid = grid to compute return when return.scale = "log"/"origin",
  #          missing -> fpca.res$workGrid.
  #   extend.par.range, check, check.sample = passed to fpca2DenFam, use default
  #   if missing, set.fpcaEstFun_control may alter some default values.
  #   ls.cov: (optional, internal use only) a list of condVarT results to save time.
  # returns:
  # either the BLUP e-coord or the corresponding log-density / density on grid,
  # or the optimization result from m2eCoord.
  # details:
  #   Use m <--> e tranlation, get BLUP of m-coord with new observation
  #   then translate m --> e-coord for ease of calculation of density function.
  #   ls.cov should be a list of lists labeled by num.k, each sub-list should contain
  #   E_mCo, cov_mCo and E_cond_cov_T used for BLUP with num.k eigenfunctions.

  # filling default and check input control
  control <- set.fpcaEstFun_control(environment())
  nk = control$num.k
  # getting a dummy density family.
  den.fam <- fpca2DenFam(fpca.res, control)

  if(is.null(control$ls.cov[[as.character(nk)]])){
    # let condVarT compute m.co, nrml.const and conditional covariance.
    ls.cov <- condVarT(
      fpca.res = fpca.res, e.co = fpca.res$xiEst[, seq_len(nk), drop = FALSE],
      return.all = TRUE
    )
    E_mCo <- colMeans(ls.cov$m.co)
    cov_mCo <- cov(ls.cov$m.co)

    # E(Var(T|natural par))
    E_cond_cov_T <- matrix(0, nk, nk)
    E_cond_cov_T[lower.tri(E_cond_cov_T, diag = T)] <- # must be lower.tri, R index row first
      colMeans(ls.cov$cond.cov)
    E_cond_cov_T <- E_cond_cov_T + t(E_cond_cov_T)
    diag(E_cond_cov_T) <- diag(E_cond_cov_T) / 2
  }else{
    ls.cov <- control$ls.cov[[as.character(nk)]]
    E_mCo <- ls.cov$E_mCo
    cov_mCo <- ls.cov$cov_mCo
    # E(Var(T|natural par))
    E_cond_cov_T <- ls.cov$E_cond_cov_T
  }

  # var T, more precisely, var of mean of T
  varT <- E_cond_cov_T / length(newx) + cov_mCo
  # new \bar T from newx
  Tnew <- den.fam$suff.f(newx)
  Tnew <- colMeans(Tnew)
  # BLUP of new mCo
  mCoEst <- cov_mCo %*% MASS::ginv(varT) %*% (Tnew - E_mCo) + E_mCo

  # translate back to e-coord
  m2e.control <- control
  m2e.control$return.scale <- 'optim'
  opt.res <- m2eCoord(mCoEst, fpca.res, control = m2e.control)
  eCoEst <- opt.res$par

  # return
  if (control$return.scale == 'optim') {
    return(opt.res)
  }
  if (control$return.scale == 'log') {
    return(log(den.fam$pdf(control$grid, opt.res$par)))
  }
  if (control$return.scale == 'origin') {
    return(den.fam$pdf(control$grid, opt.res$par))
  }
  if (control$return.scale == 'parameter') {
    return(opt.res$par)
  }
}

# set.seed(100)
# plot(y = fpca_BLUP(rnorm(10, 0,10), fpca_res, scale = '1'), x = fpca_res$workGrid, type = 'l')
# points(x = fpca_res$workGrid, y = dnorm(fpca_res$workGrid, 0, 10), type = 'l', col = 'red')

fpca_MKE = function(newx, fpca.res, control = list(wass.p = 2)){
  # computing MKE, namely minimizing Wasserstein distance
  # args:
  # newx = array of new observations;
  # fpca.res = what fdapace::FPCA returns;
  # wass.p = power of Wasserstein metric;
  # control = a list of controling arguments:
  #   num.k = number of eigenfunc to use, missing -> use all;
  #   init.par = initial value for optim, missing -> 0;
  #   method = method to use for optimization, missing -> DIRECT-P;
  #   nl.control = list of options for nloptr, c.f. nloptr::nl.opts, missing -> list();
  #   return.scale = what to return: "log", "origin", "parameter"(default), "optim";
  #   grid = grid to compute return when return.scale = "log"/"origin",
  #          missing -> fpca.res$workGrid;
  #   wass.p = pth Wasserstein distance;
  #   extend.par.range, check, check.sample = passed to fpca2DenFam, use default
  #   if missing, set.fpcaEstFun_control may alter some default values.
  # returns:
  # either the MAP or the corresponding log-density / density on grid,
  # or the optimization result, for debugging.

  # filling default and check input control
  control <- set.fpcaEstFun_control(environment())
  # getting a dummy density family.
  den.fam <- fpca2DenFam(fpca.res, control)
  # matrix of eigenfunctions, one column is one function.
  mat.eigenfunc <- fpca.res$phi[, seq_len(control$num.k), drop = FALSE]
  # objective function.
  objf <- function(par){
    # obtain values of density on workGrid, given par
    w.pdf <- exp(fpca.res$mu + mat.eigenfunc %*% par)
    # compute semi-discrete Wasserstein distance and return
    return(
      semiDscWass(
        obsv = newx, target = w.pdf,
        p = control$wass.p, grid = fpca.res$workGrid
    ))
  }

  # optimizing
  opt.res <- safeOptim(
    x0 = control$init.par,
    fn = objf, gr = NULL,
    lower = den.fam$par.range[1, ],
    upper = den.fam$par.range[2, ],
    method = control$method, backup.method = "DIRECT",
    nl.control = control$nl.control
  )

  # returning
  if (control$return.scale == 'optim') {
    return(opt.res)
  }
  if (control$return.scale == 'log') {
    return(log(den.fam$pdf(control$grid, opt.res$par)))
  }
  if (control$return.scale == 'origin') {
    return(den.fam$pdf(control$grid, opt.res$par))
  }
  if (control$return.scale == 'parameter') {
    return(opt.res$par)
  }
}

fpca_MOM = function(newx, fpca.res, control = list()){
  # Method of moment estimates of newx given FPCA result.
  # args:
  # newx = array of new observations;
  # fpca.res = what fdapace::FPCA returns
  # control = a list of controling arguments:
  #   num.k = number of eigenfunc to use, missing -> use all;
  #   init.par = initial value for optim, missing -> 0;
  #   method = method to use for optimization, missing -> DIRECT-P;
  #   nl.control = list of options for nloptr, c.f. nloptr::nl.opts, missing -> list();
  #   return.scale = what to return: "log", "origin", "parameter"(default), "optim";
  #   grid = grid to compute return when return.scale = "log"/"origin",
  #          missing -> fpca.res$workGrid.
  #   extend.par.range, check, check.sample = passed to fpca2DenFam, use default
  #   if missing, set.fpcaEstFun_control may alter some default values.
  # returns:
  # either the MOM estimate or the corresponding log-density / density on grid,
  # or the optimization result, for debugging.

  # browser()

  # filling default and check input control
  control <- set.fpcaEstFun_control(environment())

  # getting sample moment
  sample.moment <- rep(0, control$num.k)
  for(k in seq(control$num.k)){
    sample.moment[k] <- mean(newx ^ k)
  }
  # getting a dummy density family.
  den.fam <- fpca2DenFam(fpca.res, control)
  # function for moments
  quard.mom <- function(par){
    par.moment <- sapply(
      seq(control$num.k),
      den.fam$moment,
      par = par, grid = fpca.res$workGrid
    )
    return(sum((sample.moment - par.moment)^2))
  }

  opt.res <- safeOptim(
    x0 = control$init.par,
    fn = quard.mom,
    gr = NULL,
    lower = den.fam$par.range[1, ],
    upper = den.fam$par.range[2, ],
    method = control$method, backup.method = "DIRECT",
    nl.control = control$nl.control
  )

  if (control$return.scale == 'optim') {
    return(opt.res)
  }
  if (control$return.scale == 'log') {
    return(log(den.fam$pdf(control$grid, opt.res$par)))
  }
  if (control$return.scale == 'origin') {
    return(den.fam$pdf(control$grid, opt.res$par))
  }
  if (control$return.scale == 'parameter') {
    return(opt.res$par)
  }
}

# # constant list of fpca based methods.
# kListFpcaMethod <- list(
#   FPCA_MLE = fpca_MLE,
#   FPCA_MAP = fpca_MAP,
#   FPCA_BLUP = fpca_BLUP,
#   FPCA_MKE = fpca_MKE
# )
