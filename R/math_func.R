dotL2 = function(f1, f2, range, grid = NA){
  # f1 & f2 are functions, range are range, say c(0,5)
  # if grid != NA, then
  #   f2 & f2 are vectors, i.e. f1(grid), per say.
  #   and in this case grid should be a vector
  #   of the same size.

  if(is.na(grid[1])){
    return(integrate(function(x){f1(x)*f2(x)}, min(range), max(range))$value)
  }else{
    # error handling
    nGrid = length(grid)
    if(length(f1) != nGrid || length(f2) != nGrid){
      stop('Check Input vector length.')
    }
    # do not join them together, will slow down X 100
    #tm = c(f1, f2, grid)
    if(anyNA(f1) | !all(is.finite(f1)) | any(is.null(f1)) |
       anyNA(f2) | !all(is.finite(f2)) | any(is.null(f2)) |
       anyNA(grid) | !all(is.finite(grid)) | any(is.null(grid))){
      stop('non-finite function value')
    }
    # # special case so that we can use all.equal out later
    # if(nGrid == 2){
    #   return(sum(f1*f2)/2 * diff(grid))
    # }

    # # equally spaced grid, native R
    # tm_diff = diff(grid)
    # if(TRUE == do.call(all.equal, as.list(tm_diff))){
    #   return( (sum(f1*f2) - 0.5 * f1[1]*f2[1] - 0.5 * f1[nGrid] * f2[nGrid]) * tm_diff[1] )
    # }
    return(cpp_dotL2(f1, f2, grid))

    # Legacy  of dotL2=============
    # work.f1 = f1
    # work.f2 = f2
    # tm = diff(grid)
    # n = length(tm)
    # tm = rep(tm, c(1, rep(2, n-2), 1))
    # tm = colSums(matrix(tm, nrow = 2))/2
    # res = sum( work.f1 * work.f2 * c( (grid[2]-grid[1])/2, tm, (grid[n+1]-grid[n])/2 ) )
    # return(res)

    # it seems this, while more accurate, lead to some Inf problem
    # return(
    #   integrate(
    #   splinefun(
    #     x = grid, y =  f1 * f2, method = 'natural'),
    #   min(grid), max(grid))$value
    #   )

    # return(
    #   integrate(
    #   approxfun(
    #     x = grid, y =  f1 * f2, method = 'linear', rule = 2),
    #   min(grid), max(grid))$value
    #   )
  }
}
# assignInNamespace('dotL2', dotL2, "densityFPCA")
# test with integrate(splinefun)
# pretty much the same
# library(microbenchmark)
# integrate(splinefun(x = nrml_trial$workGrid, y = nrml_trial$phi[,1], method = 'natural'), min(nrml_trial$workGrid), max(nrml_trial$workGrid)) %>% microbenchmark(times=10000)
# dotL2(nrml_trial$phi[,1], rep(1, length(nrml_trial$workGrid)), grid = nrml_trial$workGrid)%>% microbenchmark(times=10000)

normL2 = function(f, range, grid = NA){
  return(sqrt(dotL2(f,f,range,grid)))
}

distHausdorff = function(space1, space2, range, grid = NA){
  dim1 = length(space1)
  dim2 = length(space2)
  if(dim1 != dim2){
    return(1)
  }else{
    p = dim1
    # orthogonize first
    # TBC

    mat_g = matrix(0, p,p)
    for(i in 1:p){
      for(j in 1:p){
        mat_g[i,j] = dotL2(space1[[i]], space2[[j]], range = range, grid = grid)
      }
    }

  }
  return(max(1 - (svd(mat_g)$d)^2))
}

matGram = function(orth1, orth2, range, grid = NA){
  # compute "Gram" matrix between orth1 and orth2
  # orth1 & orth2: list of functions or arrays.
  # range: the domain if provided functions.
  # grid: grid point if provided arrays.
  # return:
  # a length(orth1) X length(orth2) matrix whose (i, j) element is
  # < orth1[[i]], orth2[[j]]>_{L^2}

  p = length(orth1)
  q = length(orth2)

  mat_g = matrix(0, p, q)
  for(i in 1:p){
    for(j in 1:q){
      mat_g[i,j] = dotL2(orth1[[i]], orth2[[j]], range, grid)
    }
  }
  return(mat_g)
}

orthogonalize = function(base, range, grid=NA){
  p = length(base)
  res = list()

  res[[1]] = function(x){
    base[[1]](x)/sqrt(dotL2(base[[1]],base[[1]], range, grid))
  }

  for(i in 2:p){
    proj = function(x){0}
    for(j in 1:i){
      proj = function(x){ proj(x) +
        dotL2(base[[i]], res[[j]], range, grid) * res[[j]](x)
      }
    }
    tm = function(x){base[[j]](x) - proj(x)}
    res[[i]] = function(x){tm(x)/sqrt(dotL2(tm,tm, range, grid))}
  }
  return(res)
}


semiDscWass = function(obsv, target, p, grid){
  # compute semi discrete Wasserstein_p distance
  # obsv = array of observations
  # target = quantile function or array of pdf values evaluated at grid
  #          or values proportional to the actual pdf. (i.e., normalizing not needed)
  # check input TBD

  if(typeof(target) == 'closure'){
    if(!missing(grid)) warnings('Unused arguement: grid.')
    # compute ecdf
    cdf = ecdf(obsv)
    # getting its inverse defined as F^{-1}(t) = inf{x | F(x) >= t}
    # and there will be no ties in values of ecdf
    # ecdf knots and values are distinct, last value = 1
    # remark: the documentation for stepfun is bloody broken
    # it seems the writer was lost about what it means to be left continuous.
    w.qfunc = stepfun(x = environment(cdf)$y,
                      y = (environment(cdf)$x)[c(1, 1:length(environment(cdf)$x))],
                      f = 1,
                      right = FALSE)
    # a better documented function is approxfun
    # w.qfunc = approxfun(x = c(0,environment(cdf)$y),
    #                          y = c(environment(cdf)$x, max(grid)),
    #                          method = 'constant',
    #                          ties = min, yleft = min(grid), yright = max(grid))

    # integration
    diff.func = function(x){abs(w.qfunc(x) - target(x))^p}
    return(integrate(diff.func, 0, 1)$value)
  }

  # if target is density
  if(!missing(grid) & typeof(target) == 'double'){
    nGrid = length(grid)
    if(nGrid != length(target)) stop('Check length of grid and target.')
    # calculate cdf from target
    w.cdf = c(0, cumsum(0.5*diff(grid)*(target[2:nGrid] + target[(2:nGrid) - 1])) )
    # re-unify
    w.cdf = w.cdf/max(w.cdf)

    # construct common grid
    cmGrid = sort(unique(c(obsv, grid)))
    # pad cmGrid so that both w.cdf and obsv.cdf starts with 0
    cmGrid = c(min(cmGrid) - 1, cmGrid)
    nGrid = length(cmGrid)

    # interpolate w.cdf onto cmGrid
    w.cdf = approx(x = grid, y = w.cdf, xout = cmGrid,
                   yleft = 0, yright = 1)$y

    # calculate ecdf from obsv
    obsv.cdf = (ecdf(obsv))(cmGrid)

    if(anyNA(w.cdf) | !all(is.finite(w.cdf)) | any(is.null(w.cdf)) |
       anyNA(obsv.cdf) | !all(is.finite(obsv.cdf)) | any(is.null(obsv.cdf))){
      stop('Non finite values of cdf.')
    }
    return(cpp_distWass(obsv.cdf, w.cdf, cmGrid, as.double(p)))
  }
}
# Legacy of semiDscWass ===================
#   w.grid = seq(0.001, 0.999, length.out =500)
#   w.obsv = quantile(obsv, probs = w.grid)
#
#   if(typeof(target) == 'closure'){
#     if(!missing(grid)) warnings('Unused arguement: grid.')
#     w.qfunc = target(w.grid)
#   }
#   if(!missing(grid) & typeof(target) == 'double'){
#     nGrid = length(grid)
#     if(nGrid != length(target)) stop('Check length of grid and target.')
#     # if(anyNA(target) | !all(is.finite(target)) | any(is.null(target))) stop('Non finite values of target.')
#     # w.cdf = cpp_calCdf(pdf = target, grid = grid)
#     w.cdf = c(0, cumsum(0.5*diff(grid)*(target[2:nGrid] + target[(2:nGrid) - 1])) )
#     # re-unify
#     w.cdf = w.cdf/max(w.cdf)
#     #if(anyNA(w.cdf)) stop('NA in working cdf.')
#     # ties = min, so that if cdf flat, min(grid where flat) is taken
#     # since quantile func(t) := inf(x| F(x) >= t)
#     w.qfunc = approx(x = w.cdf, y = grid, xout = w.grid,
#                      rule = 2, ties = min)$y
#   }
#   w.diff = abs(w.obsv - w.qfunc)^p
#   if(anyNA(w.diff) | !all(is.finite(w.diff)) | any(is.null(w.diff))){
#     stop('Non finite values of |* - *|^p.')
#   }
#   return((cpp_dotL2(w.diff, rep(1, length(w.grid)), grid = w.grid))^(1/p))
#   # return((integrate(approxfun(x = w.grid, y = w.diff, rule = 2), lower = min(w.grid), upper = max(w.grid)))$value^(1/p))



distWass = function(f1, f2, p, grid, grid2){
  # compute Wasserstein_p distance
  # f1 & f2 = quantile function or array of pdf values evaluated at grid
  #           or separately at grid and grid2,
  #           note that they suffice to be proportional to pdf.
  #           and if quantile function, must not return Inf
  # grid will be used as grid of the 1st of f1 & f2 that is not a function
  # grid2 will be used as grid of the 2nd of f1 & f2 that is not a function
  # if one wish to compute wasserstein distant for two samples, use for example
  # distWass(f1 = as.numeric(table(sample1))/length(sample1),
  #          f2 = as.numeric(table(sample2))/length(sample2),
  #          p = ...,
  #          grid = sort(unique(sample1)),
  #          grid2 = sort(unique(sample2)))

  # check input TBD

  if(typeof(f1) == 'closure' & typeof(f2) == 'closure'){
    return(
      integrate(function(x){abs(f1(x) - f2(x))^p}, 0, 1)$value
    )
  }

  # swap order if necessary, s.t. function comes the first
  if(typeof(f1) == 'double' & typeof(f2) == 'closure'){
    tm = f2
    f2 = f1
    f1 = tm
    rm(tm)
  }

  if(typeof(f1) == 'double'){ # so that both f1 & f2 are array of density
    if(missing(grid) & missing(grid2)) stop('check input type and arguements.')
    if(!missing(grid2)){
      if(length(f1) != length(grid) | length(f2) != length(grid2)) stop('Check input length.')
      cmGrid = sort(unique(c(grid, grid2)))
      f1 = approx(x = grid, y = f1, xout = cmGrid, yleft = 0, yright = 0)$y
      f2 = approx(x = grid2, y = f2, xout = cmGrid, yleft = 0, yright = 0)$y
      # browser()
    }else{
      if(length(f1) != length(grid) | length(f2) != length(grid)) stop('Check input length.')
      cmGrid = grid
    }
    nGrid = length(cmGrid)
    w.cdf.1 = c(0, cumsum(0.5*diff(cmGrid)*(f1[2:nGrid] + f1[(2:nGrid) - 1])) )
    w.cdf.2 = c(0, cumsum(0.5*diff(cmGrid)*(f2[2:nGrid] + f2[(2:nGrid) - 1])) )
    # re-unify
    w.cdf.1 = w.cdf.1/max(w.cdf.1)
    w.cdf.2 = w.cdf.2/max(w.cdf.2)
    # sanity check
    if(anyNA(w.cdf.1) | !all(is.finite(w.cdf.1)) | any(is.null(w.cdf.1)) |
       anyNA(w.cdf.2) | !all(is.finite(w.cdf.2)) | any(is.null(w.cdf.2))){
      stop('Non finite values of cdf.')
    }
    return(
      cpp_distWass(w.cdf.1, w.cdf.2, cmGrid, p)
    )

    # some old code
    #if(anyNA(w.cdf)) stop('NA in working cdf.')
    # ties = min, so that if cdf flat, min(grid where flat) is taken
    # since quantile func(t) := inf(x| F(x) >= t)
    # w.qfunc.1 = approx(x = w.cdf.1, y = cmGrid, xout = w.grid,
    #                    rule = 2, ties = min)$y
  }

  # reach here only it f1 is quantile function and f2 is density
  nGrid = length(grid)
  w.grid = seq(0.001, 0.999, length.out =1e+3)
  w.qfunc.1 = f1(w.grid)
  if(typeof(f2) == 'double'){
    if(length(f2) != length(grid)) stop('Check input length.')
    w.cdf.2 = c(0, cumsum(0.5*diff(grid)*(f2[2:nGrid] + f2[(2:nGrid) - 1])) )
    # re-unify
    w.cdf.2 = w.cdf.2/max(w.cdf.2)
    #if(anyNA(w.cdf)) stop('NA in working cdf.')
    w.qfunc.2 = approx(x = w.cdf.2, y = grid, xout = w.grid,
                       rule = 2, ties = min)$y
  }

  w.diff = abs(w.qfunc.1 - w.qfunc.2)^p
  if(anyNA(w.diff) | !all(is.finite(w.diff)) | any(is.null(w.diff))){
    stop('Non finite values of |* - *|^p.')
  }
  return((cpp_dotL2(w.diff, rep(1, length(w.grid)), grid = w.grid))^(1/p))
  # return((integrate(approxfun(x = w.grid, y = w.diff, rule = 2), lower = min(w.grid), upper = max(w.grid)))$value^(1/p))
}


tolVar <- function(f, range){
  # total variation of function f
  # args:
  # f = an array of values of function on grid, or a function;
  # range(optional) = needed only if f is function;
  # returns:
  # a number.
  # Detail:
  # numeric approx with f linearly interpolated.

  if (!missing(range) & class(f) == 'function') {
    abs.deriv.f <- function(x) abs(nloptr::nl.grad(x, f))
    w.f <- Vectorize(abs.deriv.f)
    return(integrate(
      w.f,
      min(range),
      max(range)
    )$value)
  }
  if (typeof(f) == 'double' && all(is.finite(f))) {
    return(sum(abs(diff(f))))
  }
  stop('Check input.')
}
