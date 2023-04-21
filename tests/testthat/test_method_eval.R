context('Testing method_eval.R')
local.tol = 1e-8

n.train <- 200
train.size <- Inf
n.test <- 100
test.size <- 200

ls.den = list(
  nrml_family = nrml_family,
  # tNrml_family = tNrml_family,
  beta_family = beta_family,
  # tBeta_family = tBeta_family,
  # cauchy_family = cauchy_family,
  # tCauchy_family = tCauchy_family,
  xdCauchy_family = xdCauchy(),
  biM_family = biM_family(),
  xdBiM_family = xdBiM(),
  mixGauss_family = mixGauss_family,
  # tMixGauss_family = tMixGauss_family,
  expSin = expSin(),
  expPoly = expPoly()
)

ls.tDen <- list(
  tNrml = tNrml,
  tCauchy = tCauchy,
  tMixGauss = tMixGauss,
  tGamma = tGamma,
  tStudent = tStudent
)
for (idx.tf in seq_along(ls.tDen)) {
  for (edge in c(1,2,5,7,10,20)) {
    trunc.fam <- ls.tDen[[idx.tf]]
    nm <- paste(names(ls.tDen)[idx.tf], 'at', edge)
    ls.den[[nm]] <- trunc.fam(edge)
  }
}

for(i in seq_along(ls.den)){
  name.fam <- names(ls.den)[i]
  den.fam <- ls.den[[i]]
  set.seed(i)
  sim.dat <- genSimData(
    density.family = den.fam,
    n.train = n.train, train.size = train.size,
    n.test = n.test, test.size = test.size,
    grid.size = 256
  )

  grid <- sim.dat$grid

  # checking sampler
  arr.pval <- rep(0, nrow(sim.dat$test.sample))
  for(j in seq_len(length(arr.pval))){
    tm.pdf <- den.fam$pdf(grid, sim.dat$test.par[j, ])
    tm.pdf <- tm.pdf / max(tm.pdf)
    arr.pval[j] <- ksTest(sim.dat$test.sample[j, ], tm.pdf, grid)$p.value
  }
  q.pval <- quantile(arr.pval, c(0.25, 0.75))

  # hist(sim.dat$test.sample[1, ])
  tol.p <- c(0.1, 0.9)
  if(name.fam == 'xdBiM_family') tol.p[1] <- 0.02
  test_that(
    paste(name.fam, 'random generator'), {
      expect_true(!anyNA(sim.dat$test.sample), info = 'NA in sample!')
      expect_true(
        q.pval[1] > tol.p[1] & q.pval[2] < tol.p[2],
        info = sprintf(
          '%s p-val Q1 & Q3 as %f %f', name.fam, q.pval[1], q.pval[2]
        )
      )
    }
  )

  #
  # fpcaEsti.f(
  #   fpca.res,
  #   esti.method = c('FPCA_MLE', 'FPCA_MAP', 'FPCA_BLUP'),
  #   control = list(
  #     num.k = pick.num.k, max.k = 10,
  #     method = 'LBFGS', return.scale = c('function', 'optim'),
  #     extend.par.range = 5
  #   )
  # )
  #
  # assErr.cv(
  #   obsv = mat2list(sim.dat$test.sample, by = 1)[seq(10)],
  #   esti.f = list(
  #
  #   )
  # )

}


#
# res.fit <- fitData(
#   new.obsv = simDat$test.sample,
#   train.sample = simDat$train.sample,
#   grid = simDat$grid
# )
# mat.true <- matDensity(tNrml(6), simDat$test.par, grid = simDat$grid, n.rep = 4)
# res.err <- assErr(
#   actual = amt.true,
#   esti = res.fit$mat.esti,
#   grid = simDat$grid,
#   idx.df = res.fit$idx
# )
#
# str(res.err)
# res.err$time
# res.err$err %>% filter(idx.obsv == 1, type.err == 'KL')
#
# hand.fpca <- preTrain(simDat$train.sample, simDat$grid)
# hand.fpca$optns <- list(
#   error = FALSE,
#   lean = TRUE,
#   maxK = 5,
#   plot = FALSE,
#   useBinnedData = 'OFF'
# )
#
# fpca.res <- do.call(fdapace::FPCA, hand.fpca)
# idx.pick = 3
# hand.esti <- fpca_BLUP(simDat$test.sample[idx.pick, ], fpca.res, scale = 'origin')
#
# res.err$err %>% filter(idx.obsv == idx.pick, esti.method == 'FPCA_BLUP') %>%
#   select(err) %>% as.matrix ==
#   calErr(mat.true[idx.pick, ], hand.esti, grid = simDat$grid)
