if (FALSE)
  skip("Too much time.")

context('Testing estimating functions.')

den.fam <- tNrml(20)
set.seed(42)
dat <- genSimData(
  den.fam,
  n.train = 100, train.size = 500,
  n.test = 1, test.size = 1
)

# computing MLE with known family
mle.res <- knowFamEsti(
  dat$train.sample,
  den.fam,
  'MLE',
  list(
    method = 'LBFGS'
  )
)
# for checking sensitivity on initial values
set.seed(42)
rand.init <- runif(
  nrow(dat$train.par) * ncol(den.fam$par.range),
  min = den.fam$par.range[1, ],
  max = den.fam$par.range[2, ]
)
ni.mle.res <- knowFamEsti(
  dat$train.sample,
  den.fam,
  'MLE',
  list(
    method = 'LBFGS',
    init.par = rand.init
  )
)

test_that(
  "MLE, init val", {
  expect_equal(
    mle.res$res, ni.mle.res$res, tolerance = 1e-6,
    info = 'MLE estimates sensitive to initial values'
  )
})

theo.mle <- matrix(0, ncol = 2, nrow = nrow(dat$train.sample))
theo.mle[, 1] <- rowMeans(dat$train.sample)
theo.mle[, 2] <- apply(
  dat$train.sample, 1,
  function(x){
    sqrt(sum((x - mean(x))^2) / length(x))
  }
)
test_that(
  'MLE with known family',
  expect_equal(
    mle.res$res, theo.mle, tolerance = 1e-2
  )
)

### testing fpca related methods

# fpca_MLE
# provided with a dummy fpca.res, fpca_MLE should provide estimates
# close to natural parameters.
grid <- seq(-20, 20, length.out = 512)
canon.par <- dat$train.par
canon.par[, 2] <- -0.5 / canon.par[, 2]^2
canon.par[, 1] <- dat$train.par[, 1] / dat$train.par[, 2] ^ 2
dummy.fpca.res <- list(
  workGrid = grid,
  mu = rep(0, length(grid)),
  phi = matrix(
    c(grid, grid^2, rnorm(512, sd = 0.5), rnorm(512, sd = 0.1)),
    ncol = 4
  ),
  xiEst = cbind(
    canon.par,
    matrix(rnorm(nrow(dat$train.sample) * 2, sd = 0.25), ncol = 2)
  ),
  lambda = c(apply(canon.par, 2, var), 0.25, 0.01)
)
# testing fpca2DenFam
test_that(
  'fpca2DenFam', {
    d.den.fam <- fpca2DenFam(dummy.fpca.res, list(num.k = 2))
    expect_equal(
      d.den.fam$suff.f(grid),
      matrix(c(grid, grid^2), ncol = 2),
      info = 'fpca2DenFam: wrong suff.f.'
    )
    expect_equal(
      checkDenFamNumeric(den.fam),
      NULL,
      info = 'fpca2DenFam: potentially numeric instable'
    )
    expect_equal(
      den.fam$pdf(grid, c(0, 1)),
      d.den.fam$pdf(grid, c(0, -1/2)),
      info = 'fpca2DenFam: wrong pdf'
    )
  }
)
# dummy.fpca.res$phi %>% matplot(type = 'l')
n2.res <- fpcaEsti(
  dat$train.sample, dummy.fpca.res, 'FPCA_MLE',
  control = list(num.k = 2, method = 'LBFGS')
)
test_that(
  'FPCA_MLE',
  expect_equal(
    n2.res$res, canon.par,
    tolerance = 1e-1
  )
)

# # testing FPCA_MOM
# dummy.fpca.den.fam <- fpca2DenFam(dummy.fpca.res, list(num.k = 2))
# idx <- 2
# tst.par <- fpca_MOM(
#   dat$train.sample[idx, ], dummy.fpca.res,
#   list(
#     num.k = 2, method = 'LBFGS',
#     return.scale = 'optim',
#     init.par = n2.res$res[idx, ]
#   )
# )
# tst.par
# tst.par <- tst.par$par
# dummy.fpca.den.fam$moment(1, tst.par)
# fpca2DenFam(dummy.fpca.res)$moment(2, tst.par) - fpca2DenFam(dummy.fpca.res)$moment(1, tst.par)^2
# dat$train.sample[idx, ] %>% mean
# dat$train.sample[idx, ] %>% var
# plot(
#   x = grid,
#   y = par2pdf(dummy.fpca.den.fam, dummy.fpca.den.fam$fill.par(tst.par), grid = grid) %>% as.numeric
#   ,type = 'l'
# )
# lines(
#   x = grid,
#   y = par2pdf(
#     den.fam,
#     dat$train.par[idx, ],
#     grid = grid
#   )
#   , col = 'red'
# )
#
# mom.res <- fpcaEsti(
#   dat$train.sample, dummy.fpca.res, 'FPCA_MOM',
#   control = list(num.k = 2, method = 'LBFGS')
# )
# mat.moment <- cbind(dat$train.par, dat$train.par)
# mat.moment[, 1] <- rowMeans(dat$train.sample)
# mat.moment[, 2] <- rowMeans(dat$train.sample ^ 2)
# mat.moment[, 3] <- apply(mom.res$res, 1, function(par) dummy.fpca.den.fam$moment(1, par))
# mat.moment[, 4] <- apply(mom.res$res, 1, function(par) dummy.fpca.den.fam$moment(2, par))
# plot(mat.moment[, seq(2)] - mat.moment[, seq(2) + 2])
# summary(mat.moment[, seq(2)] - mat.moment[, seq(2) + 2])

# dummy.fpca.den.fam$moment(1, mom.res$res[2, ])
# dummy.fpca.den.fam$moment(2, mom.res$res[2, ])
# mat.moment[2, ]

# TBD: how do we test for fpca_BLUP and fpca_MAP?

# testing AIC implementation
set.seed(42)
dat <- genSimData(
  den.fam,
  n.train = 100, train.size = 25,
  n.test = 1, test.size = 1
)

aic.res <- fpcaEsti(
  dat$train.sample, dummy.fpca.res, 'FPCA_MLE',
  control = list(num.k = 'AIC', method = 'LBFGS', max.k = 4)
)
n1.res <- fpcaEsti(
  dat$train.sample, dummy.fpca.res, 'FPCA_MLE',
  list(num.k = 1, method = 'LBFGS')
)
n2.res <- fpcaEsti(
  dat$train.sample, dummy.fpca.res, 'FPCA_MLE',
  control = list(num.k = 2, method = 'LBFGS')
)
n3.res <- fpcaEsti(
  dat$train.sample, dummy.fpca.res, 'FPCA_MLE',
  list(num.k = 3, method = 'LBFGS')
)
n3.aic <- 2 * sapply(
  seq_len(nrow(dat$train.sample)),
  function(i){
    getNegLogll(
      obsv = dat$train.sample[i, ],
      den.fam = fpca2DenFam(dummy.fpca.res, list(num.k = 3))
    )(n3.res$res[i, ])
  }
) + 2 * 3
n2.aic <- 2 * sapply(
  seq_len(nrow(dat$train.sample)),
  function(i){
    getNegLogll(
      obsv = dat$train.sample[i, ],
      den.fam = fpca2DenFam(dummy.fpca.res, list(num.k = 2))
    )(n2.res$res[i, ])
  }
) + 2 * 2
n1.aic <- 2 * sapply(
  seq_len(nrow(dat$train.sample)),
  function(i){
    getNegLogll(
      obsv = dat$train.sample[i, ],
      den.fam = fpca2DenFam(dummy.fpca.res, list(num.k = 1))
    )(n1.res$res[i, ])
  }
) + 2 * 1
mat.aic <- matrix(c(n1.aic, n2.aic, n3.aic), ncol = 3)
test_that(
  'minize AIC?', {
  expect_true(
    all(apply(mat.aic, 1, which.min) == aic.res$idx$num.k),
    'fpcaEst fail to get minimized AIC.'
  )
  expect_true(
    all(aic.res$idx$num.k[aic.res$idx$AIC == mat.aic[, 2]] == 2) &
    all(aic.res$idx$num.k[aic.res$idx$AIC == mat.aic[, 1]] == 1) &
    all(aic.res$idx$num.k[aic.res$idx$AIC == mat.aic[, 3]] == 3),
    "Error in calculating AIC, values mismatch for num.k = 'AIC' and num.k = some number."
  )
})

# testing wrapper fpcaEsti
esti.res <- fpcaEsti(
  dat$train.sample, dummy.fpca.res,
  esti.method = c('FPCA_MLE', 'FPCA_MAP', 'FPCA_BLUP', 'FPCA_MKE'), #, 'FPCA_MOM'),
  list(
    num.k = 2, method = 'LBFGS'
  )
)

# check with functional wrapper
ls.esti.f <- fpcaEsti.f(
  dummy.fpca.res,
  c('FPCA_MLE', 'FPCA_MAP', 'FPCA_BLUP', 'FPCA_MKE'), #, 'FPCA_MOM'),
  control = list(
    num.k = 2, method = 'LBFGS',
    return.scale = c('par', 'optim', 'func')
  )
)

esti.res.f <- list()
for(i in seq_along(ls.esti.f)){
  esti.res.f[[i]] <- apply(dat$train.sample, 1, ls.esti.f[[i]])
}
names(esti.res.f) <- names(ls.esti.f)

test_that('fpcaEsti funcitonal wrapper', {
  expect_equal(
    esti.res$res,
    do.call(
      rbind,
      lapply(esti.res.f, function(x){
        t(sapply(x, '[[', 'parameter'))
      })
    ),
    info = 'result mismatch between fpcaEsti and fpcaEsti.f'
  )
  expect_equal(
    fpca2DenFam(dummy.fpca.res, list(num.k = 2))$pdf(
      grid,
      esti.res.f[[1]][[1]][['parameter']]
    ),
    esti.res.f[[1]][[1]][['function']](grid),
    info = 'fpcaEsti.f, par results and func results inconsistent.'
  )
})

# check for mismatch
fpcamle.res <- t(apply(
  dat$train.sample, 1,
  fpca_MLE,
  fpca.res = dummy.fpca.res,
  control = list(num.k = 2, method = 'LBFGS', return.scale = 'par')
))
fpcamap.res <- t(apply(
  dat$train.sample, 1,
  fpca_MAP,
  fpca.res = dummy.fpca.res,
  control = list(num.k = 2, method = 'LBFGS', return.scale = 'par')
))
fpcablup.res <- t(apply(
  dat$train.sample, 1,
  fpca_BLUP,
  fpca.res = dummy.fpca.res,
  control = list(num.k = 2, method = 'LBFGS', return.scale = 'par')
))
fpcamke.res <- t(apply(
  dat$train.sample, 1,
  fpca_MKE,
  fpca.res = dummy.fpca.res,
  control = list(num.k = 2, method = 'LBFGS', return.scale = 'par')
))

test_that('fpcaEsti wrapper', {
  expect_equal(
    fpcamle.res,
    esti.res$res[esti.res$idx$esti.method == 'FPCA_MLE', ],
    info = 'result mismatch for method FPCA_MLE.'
  )
  expect_equal(
    fpcamap.res,
    esti.res$res[esti.res$idx$esti.method == 'FPCA_MAP', ],
    info = 'result mismatch for method FPCA_MAP.'
  )
  expect_equal(
    fpcablup.res,
    esti.res$res[esti.res$idx$esti.method == 'FPCA_BLUP', ],
    info = 'result mismatch for method FPCA_BLUP.'
  )
  expect_equal(
    fpcamke.res,
    esti.res$res[esti.res$idx$esti.method == 'FPCA_MKE', ],
    info = 'result mismatch for method FPCA_MKE.'
  )
})

# testing EM/MLE for Gaussian mixture family
rm(list = ls())
test_that('EM/MLE for tMixGauss', {
  den.fam <- tMixGauss(5, num.mixture = 3)
  dat <- list()
  dat$train.par <- t(matrix(c(
    tMixGauss.repar(
      ls.par = list(
        mean = c(-2, 0, 2),
        var = rep(0.5, 3),
        weight = c(0.3, 0.4, 0.3)
      )
    ),
    tMixGauss.repar(
      ls.par = list(
        mean = c(-3, 0, 2),
        var = rep(0.5, 3),
        weight = c(0.4, 0.3, 0.3)
      )
    )), ncol = 2))
  # den.fam$pdf(seq(-5, 5, 0.1), dat$train.par[2, ]) %>% plot(type = 'l')
  dat$train.sample <- matrix(0, nrow = nrow(dat$train.par), ncol = 2000)
  set.seed(1)
  dat$train.sample[1, ] <- den.fam$rpdf(2000, dat$train.par[1, ])
  dat$train.sample[2, ] <- den.fam$rpdf(2000, dat$train.par[2, ])


  # EM algorithm
  em.res <- knowFamEsti(
    dat$train.sample,
    den.fam,
    'em',
    list(
      return.scale = 'par',
      init.par = dat$train.par
    )
  )

  expect_equal(
    em.res$res, dat$train.par, tolerance = 1e-1,
    info = 'large error in EM estimates'
  )
  # em.res$res - dat$train.par

  mle.res <- knowFamEsti(
    dat$train.sample,
    den.fam,
    'MLE',
    list(method = 'LBFGS', init.par = dat$train.par)
  )
  mle.res$res <- den.fam$identify_par(mle.res$res)
  expect_equal(
    dat$train.par, mle.res$res, tolerance = 1e-1,
    info = 'large error in MLE estimates'
  )
  # dat$train.par - mle.res$res

  expect_equal(
    mle.res$res, em.res$res, tolerance = 1e-3,
    info = 'large difference between EM and MLE'
  )
  # mle.res$res - em.res$res

  # double check wrapper indexing
  b.res <- knowFamEsti(
    dat$train.sample,
    den.fam,
    c('MLE', 'em'),
    list(method = 'LBFGS', init.par = dat$train.par)
  )
  expect_identical(
    den.fam$identify_par(b.res$res[b.res$idx$esti.method == 'MLE', ]),
    mle.res$res,
    info = 'knowFamEsti result mismatch for MLE'
  )
  expect_identical(
    b.res$res[b.res$idx$esti.method == 'emTruncMixGauss', ],
    em.res$res,
    info = 'knowFamEsti result mismatch for EM'
  )

  # # check sensitivity on initial values
  # set.seed(42)
  # rand.init <- runif(
  #   2 * ncol(den.fam$par.range),
  #   min = den.fam$par.range[1, ],
  #   max = den.fam$par.range[2, ]
  # )
  # ni.em.res <- knowFamEsti(
  #   dat$train.sample,
  #   den.fam,
  #   'em',
  #   list(
  #     return.scale = 'par',
  #     init.par = rand.init
  #   )
  # )
  # expect_equal(
  #   em.res$res, ni.em.res$res, tolerance = 1e-3,
  #   info = 'EM estimates sensitive to initial values'
  # )
  # though, it seems if initial mixture means are close
  # EM fail to produce valid estimates.
  # emTruncMixGauss(
  #   obsv = dat$train.sample[1, ], num.mixture = 3,
  #   edge = 5,
  #   list(
  #     init.par = tMixGauss.repar(work.par = colMeans(den.fam$par.range)),
  #     tol = 1e-12
  #   )
  # )

  # check what if removing analytic grad
  den.fam$logpdf.gr <- NULL
  suppressWarnings(
    ng.mle.res <- knowFamEsti(
      dat$train.sample,
      den.fam,
      'MLE',
      list(method = 'LBFGS', init.par = dat$train.par)
    )
  )
  ng.mle.res$res <- den.fam$identify_par(ng.mle.res$res)
  expect_equal(
    ng.mle.res$res, mle.res$res,
    info = 'removing analytic gradient affects MLE computation'
  )
  # ng.mle.res$res - mle.res$res

  # check will number of components get shrinked.
  set.seed(1)
  obsv <- tNrml(3)$rpdf(150, c(0, 1))
  fin.em.res <- emTruncMixGauss(
    obsv,
    num.mixture = 50,
    edge = 3
  )
  expect_true(
    fin.em.res$num.mixture < 50,
    info = 'emTruncMixGauss: does not drop mixtures'
  )

})
