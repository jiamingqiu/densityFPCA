# some draft on constructing a general k-fold CV routine.
# so that can be used in selecting presmoothing bw and padding (preSmooth.kde)
# or number of knots and degree (preSmooth.logspline)

kFoldCV.f <- function(dat, f, k, control = list(), ...){
  # k-fold cross-validation for tuning f with dat.
  # args:
  #  dat: a matrix where one row is one sample, or a list
  #  f: the objective function with argument train, test, tune, ...
  #     must be compatible with the format of input dat, i.e.
  #     list for list, matrix for matrix
  #  k: number of folds
  #  control: a list of
  #    n.rep: a non-negative integer for repeated k-fold CV,
  #    0(default/missing) means split data deterministically(no shuffle)
  #    1 means shuffle and split;
  #    aggr: function to combine evaluations of objective function f,
  #    mean(default/missing) takes mean.
  #  ...: other arguments for f.
  # returns:
  #  a function g with argument tune, being
  #  g(tune, ...) = mean( aggr(f_tune over folds) ), mean over n.rep.
  # note: the data splition would not change once g is created, therefore one
  #  should NOT run such returned function multiple times for repeated CV.

  # #DEBUG
  # f <- function(train, test, tune, deg){
  #   fit <- mean(train) + tune
  #   loss <- mean(abs(test - fit) ^ deg)
  #   return(loss)
  # }
  # set.seed(1)
  # dat <- matrix(rnorm(500 * 2), ncol = 2)
  # k <- 500
  # control <- list(n.rep = 10, aggr = function(x) prod(abs(x)))

  # sanity check
  if(is.matrix(dat)){
    dat <- mat2list(dat, by = 1)
  }
  if(!is.list(dat))
    stop('check input type, must be list or matrix.')
  n.sample <- length(dat)
  if(k > n.sample)
    stop('k larger than number of samples.')
  if(is.null(control$n.rep))
    control$n.rep <- 0
  if(is.null(control$aggr))
    control$aggr <- mean

  ls.split.idx <- kFoldIdx(n.sample, k, control$n.rep, TRUE)
  if(control$n.rep %in% c(0, 1))
    ls.split.idx <- list(ls.split.idx)

  ls.split.dat <- list()
  for(idx.rep in seq_along(ls.split.idx)){
    ls.split.dat[[idx.rep]] <- lapply(
      ls.split.idx[[idx.rep]],
      function(idx){
        list(
          train = dat[idx[[1]]],
          test = dat[idx[[2]]]
        )
      }
    )
  }

  cv.f <- function(tune, ...){
    arr.f <- rep(0, length(ls.split.dat))
    for(idx.rep in seq_along(ls.split.dat)){
      arr.f[idx.rep] <- control$aggr(
        sapply(
          ls.split.dat[[idx.rep]],
          function(dat){
            f(train = dat$train, test = dat$test, tune = tune, ...)
          }
        )
      )
    }
    return(mean(arr.f))
  }

  # #DEBUG
  # cv.f(0, deg = 2)
  # optimize(cv.f, c(-1, 1), deg = 2)

  return(cv.f)
}


kFoldIdx <- function(n, k, n.rep, traintest = FALSE){
  # generating index for spliting data for k-fold CV
  # args:
  #  n: number of samples
  #  k: k
  #  n.rep: a non-negative integer for repeated k-fold CV,
  #  0(default/missing) means split data deterministically(no shuffle)
  #  1 means shuffle and split;
  #  traintest: return indices for partitioning data into k groups (when FALSE) or
  #  further into labels for training and testing.
  # returns:
  #  a list of n.rep lists where each sublist contains k lists of
  #  indices for spliting data, if n.rep = 0, 1, then just list of k lists.
  #  if traintest = TRUE, each of those k lists will be, insteda for spliting data,
  #  yet another layer of list of 2 arrays, where the 1st is indices for training,
  #  and 2nd for testing.
  #
  if(n < k)
    stop('n < k.')
  if(k <= 1)
    stop('k must be at least 2.')
  if(missing(n.rep))
    n.rep <- 0

  sf.idx <- seq_len(n)
  lbl <- seq_len(k)

  cut.idx <- cut(sf.idx, breaks = k, labels = lbl)

  if(n.rep == 0){
    ls.res <- split(x = sf.idx, f = cut.idx)
    names(ls.res) <- NULL
    ls.res <- list(ls.res)
  }else{
    # shuffle and draw
    ls.sf.idx <- replicate(
      n = n.rep,
      expr = sample.int(n = n, size = n),
      simplify = FALSE
    )
    ls.res <- lapply(
      ls.sf.idx,
      function(sf.idx){
        ls.res <- split(x = sf.idx, f = cut.idx)
        names(ls.res) <- NULL
        return(ls.res)
      }
    )
  }
  if(!traintest){
    if(n.rep %in% c(0, 1))
      ls.res <- ls.res[[1]]
    return(ls.res)
  }

  # if further into training/testing
  for(idx.rep in seq_along(ls.res)){
    ls.res[[idx.rep]] <- lapply(
      ls.res[[idx.rep]],
      function(x){
        list(
          sf.idx[!(sf.idx %in% x)],
          x
        )
      }
    )
  }
  if(n.rep %in% c(0, 1))
    ls.res <- ls.res[[1]]
  return(ls.res)
}

# some function for searching

downHill.f <- function(fun, over, what){
  # construct a function for downhill search
  # args:
  # fun: objective function, one of whose arguments to be searched
  # over: a list of all possible values, name correspond to the arg searched
  # what: default(missing) then fun return a number, otherwise this is the name
  #       of what fun returns to be used as loss
  # returns:
  # a function identical to fun but "minizied" w.r.t. the arg in over, whose return
  # is a list of the original returns(if missing what, key is "value") and the
  # value of arg minimizing it.
  # details:
  # start from the first value in over, stop when the loss starts to increase.
  # values of over leading to errors will be dropped.
  # NOTE: if over and the list return of fun have some identical names,
  # names in return could be duplicated.
  # IDK how to check return of fun is list or not, be careful not to drop what.

  if(length(over) != 1)
    stop('check over, should be a list of length 1.')
  if(!all(names(over) %in% names(formals(fun))))
    stop('arg in over is not args of fun.')
  grid.cdd <- over[[1]]
  arg.nm <- names(over)
  # wrap fun if necessary
  if(missing(what)){
    what <- 'value'
    w.fun <- function(...){
      tm.res <- fun(...)
      if(!is.numeric(tm.res))
        return(tm.res)
      return(list(value = tm.res))
    }
  }else{
    w.fun <- fun
  }
  # construct return
  res.f <- function(...){
    ls.res <- list()
    ls.opt <- list()
    tm.opt <- list()
    arr.what <- numeric(length(grid.cdd))
    flg.cdd <- 0
    for(idx.cdd in seq_along(grid.cdd)){
      # set temp options
      tm.opt[[arg.nm]] <- grid.cdd[idx.cdd]
      ls.opt[[idx.cdd]] <- tm.opt
      # evaluate fun at tm.opt with ...
      ls.res[[idx.cdd]] <- do.call(
        w.fun,
        c(
          list(...),
          tm.opt
        )
      )
      # sanity, if not complying, stop search
      if(is.null(ls.res[[idx.cdd]][[what]])){
        warning(paste(
          'fun does not return list contain what, search stopped at',
          paste(names(tm.opt), collapse = ', '), '=',
          paste(unlist(tm.opt, use.names = FALSE), collapse = ', ')
        ))
        ls.res[[idx.cdd]][[what]] <- NA
        flg.cdd <- idx.cdd - 1
        break
      }else{
        arr.what[idx.cdd] <- ls.res[[idx.cdd]][[what]]
      }
      # stop if loss start climbing
      if(arr.what[idx.cdd] > arr.what[max(idx.cdd - 1, 1)]){
        flg.cdd <- idx.cdd - 1
        break
      }
      # if loss keeps descending
      flg.cdd <- idx.cdd
    }
    idx.return <- flg.cdd
    if(idx.return == 0)
      stop('non finite fun returns for initial value in over.')
    # browser()
    return(c(
      ls.res[[idx.return]],
      ls.opt[[idx.return]]
    ))
  }
  return(res.f)
}

gridSearch.f <- function(fun, over, what){
  # construct a function for grid search
  # args:
  # fun: objective function, some of whose arguments to be searched
  # over: a list of all possible values, names correspond to the args searched
  # what: default(missing) then fun return a number, otherwise this is the name
  #       of what fun returns to be used as loss
  # returns:
  # a function identical to fun but minizied w.r.t. the args in over, whose return
  # is a list of the original returns(if missing what, key is "value") and the
  # args minimizing it.
  # details:
  # values of over leading to errors will be dropped.
  # NOTE: if over and the list return of fun have some identical names,
  # names in return could be duplicated.
  # IDK how to check return of fun is list or not, be careful not to drop what.

  if(!all(names(over) %in% names(formals(fun))))
    stop('not all args in over are args of fun.')
  grid.cdd <- do.call(expand.grid, over)
  # wrap fun if necessary
  if(missing(what)){
    what <- 'value'
    w.fun <- function(...){
      tm.res <- fun(...)
      if(!is.numeric(tm.res))
        return(tm.res)
      return(list(value = tm.res))
    }
  }else{
    w.fun <- fun
  }
  # construct return
  res.f <- function(...){
    ls.res <- list()
    ls.opt <- list()
    tm.opt <- list()
    for(idx.cdd in seq_len(nrow(grid.cdd))){
      # set temp options
      for(nm in names(grid.cdd)){
        tm.opt[[nm]] <- grid.cdd[[nm]][idx.cdd]
      }
      ls.opt[[idx.cdd]] <- tm.opt
      # evaluate fun at tm.opt with ...
      ls.res[[idx.cdd]] <- do.call(
        w.fun,
        c(
          list(...),
          tm.opt
        )
      )
      # sanity
      if(is.null(ls.res[[idx.cdd]][[what]])){
        warning(paste(
          'fun does not return list contain what, candidate dropped at',
          paste(names(tm.opt), collapse = ', '), '=',
          paste(unlist(tm.opt, use.names = FALSE), collapse = ', ')
        ))
        ls.res[[idx.cdd]][[what]] <- NA
      }
    }
    # get the smallest one
    arr.what <- sapply(ls.res, function(x) x[[what]])
    idx.return <- which.min(arr.what)
    if(length(idx.return) == 0)
      stop('non finite fun returns for all values in over.')
    # browser()
    return(c(
      ls.res[[idx.return]],
      ls.opt[[idx.return]]
    ))
  }
  return(res.f)
}

