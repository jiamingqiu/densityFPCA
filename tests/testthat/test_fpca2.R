# testing fpca to other objects
context('fpca2...')

den.fam <- tNrml(20)
# den.fam <- tMixGauss(5)
# den.fam <- xdBiM()
set.seed(42)
dat <- genSimData(
  den.fam,
  n.train = 20, train.size = Inf,
  n.test = 1, test.size = 100,
  grid.size = 512
)
dat$train.curve <- toHilbert(
  dat$train.sample,
  dat$grid,
  transFun = function(f, grid){
    orthLog(f, grid, 1 / diff(den.fam$range))
  }
)$mat.curve
# matplot(t(dat$train.curve), type = 'l')

fpca.res <- do.call(
  fdapace::FPCA,
  list(
    Ly = mat2list(dat$train.curve, by = 1),
    Lt = replicate(n = nrow(dat$train.curve), expr = dat$grid, simplify = FALSE),
    optns = list(
      error = FALSE,
      lean = TRUE,
      methodSelectK = 'FVE',
      FVEthreshold = 0.9999,
      plot = FALSE,
      useBinnedData = 'OFF'
    )
  )
)

test_that('fpca2DenFam',{
  # fpca.den.fam <- fpca2DenFam(
  #   fpca.res,
  #   list(
  #     extend.par.range = 500,
  #     check = FALSE
  #   )
  # )
  # Testing checks, should throw an error
  # set.seed(42)
  # checkDenFamNumeric(fpca.den.fam, dat$test.sample[1, ])
  # calMLE(
  #   dat$test.sample[1, ], fpca.den.fam, list(
  #     method = 'DIRECT-P',
  #     return.scale = 'optim'
  #   )
  # )
  expect_warning(
    fpca2DenFam(
      fpca.res,
      list(
        extend.par.range = 6,
        check = TRUE,
        check.sample = dat$test.sample[1, ]
      )
    )
  )
  # a dummy fpca.res
  fpca.den.fam <- fpca2DenFam(
    list(
      workGrid = dat$grid,
      mu = rep(0, length(dat$grid)),
      phi = matrix(
        c(dat$grid, dat$grid^2),
        ncol = 2
      ),
      xiEst = fpca.res$xiEst[, seq(2)],
      lambda = fpca.res$lambda[seq(2)]
    ),
    control = list(
      extend.par.range = 1.5,
      check = FALSE,
      get.prior = FALSE
    )
  )
  # if supplied par is too long
  expect_error(
    fpca.den.fam$fill.par(seq(ncol(fpca.res$xiEst) + 2))
  )
  # 1st and 2nd moment of N(0, 1) and N(1, 1)
  # note using canonical par: (mu / sigma^2, -1 / 2 / sigma^2)
  expect_equal(
    c(
      fpca.den.fam$moment(1, c(0, -1/2), grid = dat$grid)
      , fpca.den.fam$moment(2, c(0, -1/2), grid = dat$grid)
      , fpca.den.fam$moment(1, c(1, -1/2), grid = dat$grid)
      , fpca.den.fam$moment(2, c(1, -1/2), grid = dat$grid)
    ),
    c(0, 1, 1, 2)
  )

})


# # testing assErr.cv
# set.seed(42)
# test.dat <- genSimData(
#   den.fam, n.train = 0, train.size = Inf,
#   # n.test = 5000, test.size = rep(seq(50) * 20, each = 5)
#   n.test = 50, test.size = seq(50) * 50
# )
# cv.res <- assErr.cv(
#   obsv = test.dat$test.sample,
#   esti.f = function(obsv){
#     function(x) dnorm(x, mean = mean(obsv), sd = sd(obsv))
#   },
#   grid = fpca.res$workGrid,
#   err.method = c('cv.ISE', 'cv.KL'),
#   idx.df = data.frame(n = sapply(test.dat$test.sample, length))
#   # , control = list(cv.cntrl = list(k = 5, n.rep = 2))
# )
# cv.res %>%
#   plyr::ddply(c('n', 'type.err'), function(df){
#     data.frame(
#       n = df$n[1], type.err = df$type.err[1],
#       err = mean(df$err)
#     )
#   }) %>%
#   dplyr::filter(type.err == 'cv.ISE') %>%
#   with(plot(n, err))
# assErr(
#   actual = par2pdf(
#     den.fam, test.dat$test.par, grid = fpca.res$workGrid
#   ),
#   esti = t(sapply(test.dat$test.sample, function(obsv){
#     dnorm(fpca.res$workGrid, mean(obsv), sd(obsv))
#   })),
#   grid = fpca.res$workGrid,
#   err.method = c('ISE', 'KL'),
#   idx.df = data.frame(n = sapply(test.dat$test.sample, length))
# ) %>%
#   .$err %>%
#   plyr::ddply(c('n', 'type.err'), function(df){
#     data.frame(
#       n = df$n[1], type.err = df$type.err[1],
#       err = mean(df$err)
#     )
#   }) %>%
#   dplyr::filter(type.err == 'KL') %>%
#   with(plot(n, log10(err)))


