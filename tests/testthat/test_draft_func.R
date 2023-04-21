if (TRUE)
  skip("LEGACY")

lst_den = list(nrml_family,
               tNrml_family,
               beta_family,
               tBeta_family,
               cauchy_family,
               tCauchy_family,
               xdCauchy(),
               biM_family(),
               xdBiM(),
               mixGauss_family,
               tMixGauss_family)

lst_name = c(
  "nrml_family",
  "tNrml_family",
  "beta_family",
  "tBeta_family",
  "cauchy_family",
  "tCauchy_family",
  "xdCauchy",
  "biM_family",
  "xdBiM",
  "mixGauss_family",
  "tMixGauss_family"
)

lst_method = list(
  fpca_BLUP,
  fpca_MAP,
  fpca_MKE,
  fpca_MLE
)

context('Testing m-e transform.')

for(idx.denFam in 1:length(lst_den)){

  denFam = lst_den[[idx.denFam]]
  set.seed(10)
  fpca_res = learnFamily(denFam, trainSize = 20, maxK = 10,
              gridSize = 64, refPoint = mean(denFam$range), simple_result = TRUE)
  # fpca_res = learnFamily(denFam, trainSize = 20, gridSize = 512)

  # determine initial values
  tm.Range = matrix(apply(fpca_res$xiEst, 2, range), nrow = 2)
  # set as so since mean(xiEst) = 0, choose the hyper-rectangle so as no to go outside.
  multiStart = c(min(9, length(fpca_res$lambda)*2 - 1), pmin(abs(tm.Range[2,]), abs(tm.Range[1,])) / 2)

  test_that(lst_name[[idx.denFam]],{

    tm_mco = matrix(t(apply(fpca_res$xiEst, 1, e2mCoord, fpca.res = fpca_res)), ncol = ncol(fpca_res$xiEst))
    recal_eco = matrix(
      t(apply(tm_mco, 1, m2eCoord,
            fpca.res = fpca_res, control = list(method = 'LBFGS'))),
      ncol = ncol(fpca_res$xiEst))
    err = abs(recal_eco - fpca_res$xiEst)
    relative.err <- ifelse(
      abs(fpca_res$xiEst) <= .Machine$double.eps^(1/3),
      err,
      err / abs(fpca_res$xiEst)
    )
    expect_true(
      max(relative.err) < 1e-3,
      info = paste('large m<->e translate error, max', max(relative.err))
    )
  })

  # test_that(paste(lst_name[idx.denFam], 'estimates'),{
  #
  #   set.seed(100)
  #   para = denFam$gen_par(1)
  #   newx = denFam$rpdf(7, para)
  #
  #   expect_silent(fpca_MLE(newx, fpca_res, multiStart = multiStart))
  #   expect_silent(fpca_MAP(newx, fpca_res, multiStart = multiStart))
  #   expect_silent(fpca_BLUP(newx, fpca_res, multiStart = multiStart))
  #   expect_silent(fpca_MKE(newx, fpca_res, multiStart = multiStart))
  #
  #   # tst = rnorm(50)
  #   # microbenchmark::microbenchmark(
  #   #   fpca_MLE(tst, nrml_res),
  #   #   fpca_MKE(tst, nrml_res, 1),
  #   #   fpca_MAP(tst, nrml_res),
  #   #   fpca_BLUP(tst, nrml_res),
  #   #   times = 10
  #   # )
  #   #
  #   # tst = rnorm(5)
  #   # plot(y = fpca_MKE(tst, nrml_res, 1, scale = 'origin'), x = nrml_res$workGrid, type = 'l', xlim = c(-3,3), col='blue')
  #   # points(x = nrml_res$workGrid, y = dnorm(nrml_res$workGrid), type = 'l', col = 'red', lwd=2)
  #   # points(y = fpca_MLE(tst, nrml_res, scale = 'origin'), x = nrml_res$workGrid, type = 'l', col = 'green')
  #
  # })
}





#   set.seed(100)
#   para = nrml_family$gen_par(1)
#   obsv = nrml_family$rpdf(5, para)
#
#
#   esti.blup = fpca_BLUP(newx = obsv,
#             fpca_res = nrml_res,
#             scale = 'origin',
#             multiStart = multiStart)
#   grid = nrml_res$workGrid
#   calErr(nrml_family$pdf(grid, para), esti.blup, grid = grid)
#   plot(x = grid, y = esti.blup, type = 'l')



# truncated families ========================
lst_den = list(tNrml, tCauchy, tMixGauss, tGamma, tStudent)
lst_name = c('truncatedNormal', 'truncatedCauchy', 'truncatedMixGauss', 'truncatedGamma', 'truncatedT')
w.edge = 7

for(idx.denFam in 1:length(lst_den)){
  denFam = lst_den[[idx.denFam]](w.edge)
  set.seed(100)
  fpca_res = learnFamily(denFam, trainSize = 50, maxK = 20,
                         gridSize = 64, refPoint = mean(denFam$range), simple_result = TRUE)

  # determine initial values
  tm.Range = matrix(apply(fpca_res$xiEst, 2, range), nrow = 2)
  # set as so since mean(xiEst) = 0, choose the hyper-rectangle so as no to go outside.
  multiStart = c(min(9, length(fpca_res$lambda)*2 - 1), pmin(abs(tm.Range[2,]), abs(tm.Range[1,])) / 2)

  test_that(lst_name[[idx.denFam]], {

    tm_mco = matrix(t(apply(fpca_res$xiEst, 1, e2mCoord, fpca.res = fpca_res)), ncol = ncol(fpca_res$xiEst))
    recal_eco = matrix(t(
      apply(
        tm_mco, 1, m2eCoord,
        fpca.res = fpca_res, control = list(method = 'LBFGS')
        )),
      ncol = ncol(fpca_res$xiEst))
    err = abs(recal_eco - fpca_res$xiEst)
    relative.err <- ifelse(
      abs(fpca_res$xiEst) <= .Machine$double.eps^(1/3),
      err,
      err / abs(fpca_res$xiEst)
    )
    expect_true(
      max(relative.err) < 1e-3,
      info = paste('large m<->e translate error, max', max(relative.err))
    )
  })

  # test_that(paste(lst_name[idx.denFam], 'estimates'),{
  #
  #   set.seed(100)
  #   para = denFam$gen_par(1)
  #   newx = denFam$rpdf(7, para)
  #
  #   expect_silent(fpca_MLE(newx, fpca_res, multiStart = multiStart))
  #   expect_silent(fpca_MAP(newx, fpca_res, multiStart = multiStart))
  #   expect_silent(fpca_BLUP(newx, fpca_res, multiStart = multiStart))
  #   expect_silent(fpca_MKE(newx, fpca_res, multiStart = multiStart))
  #
  # })
}

