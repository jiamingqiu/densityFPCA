context('Testing pre-specified density families')

tol = 0.1
num.rep.smpl <- 50
size.smpl <- 50

lst_den = list(
  nrml_family = nrml_family,
  tNrml_family = tNrml_family,
  beta_family = beta_family,
  tBeta_family = tBeta_family,
  cauchy_family = cauchy_family,
  tCauchy_family = tCauchy_family,
  biM_family = biM_family(),
  xdBiM_family = xdBiM(),
  xdCauchy = xdCauchy(),
  mixGauss_family = mixGauss_family,
  tMixGauss_family = tMixGauss_family,
  expPoly = expPoly(),
  expSin = expSin()
)

lst_tDen <- list(
  tNrml = tNrml,
  tCauchy = tCauchy,
  tMixGauss = tMixGauss,
  tGamma = tGamma,
  tStudent = tStudent
)
for (idx.tf in seq_along(lst_tDen)) {
  for (edge in c(1,2,5,7,10,20)) {
    trunc.fam <- lst_tDen[[idx.tf]]
    nm <- paste(names(lst_tDen)[idx.tf], 'at', edge)
    lst_den[[nm]] <- trunc.fam(edge)
  }
}

epsh <- 0.0005

# for(i in seq_along(lst_den)){
for(i in 39:43){
  name.fam <- names(lst_den)[i]
  denFam = lst_den[[i]]
  grid = seq(denFam$range[1], denFam$range[2], by = epsh)

  # checking random genenrator by ks.test()
  ks.p <- rep(0, num.rep.smpl)
  err.cdf <- rep(0, num.rep.smpl)
  set.seed(555)
  para = denFam$gen_par(num.rep.smpl)
  # testing parameter range
  w.par.range <- apply(para, 2, range)
  flg.par.range <- any(
    w.par.range[1, ] < denFam$par.range[1, ] |
    w.par.range[2, ] > denFam$par.range[2, ]
  )
  # observation range flag
  flg.smpl.range <- FALSE
  actual.smpl.range <- denFam$range

  # generating samples and test
  for (idx in seq_len(num.rep.smpl)) {
    set.seed(idx)
    smpl = denFam$rpdf(size.smpl, para[idx, ])

    smpl.range <- range(smpl)
    flg.smpl.range <- flg.smpl.range || (
      smpl.range[1] < denFam$range[1] ||
      smpl.range[2] > denFam$range[2]
    )
    actual.smpl.range[1] <- min(actual.smpl.range[1], smpl.range[1])
    actual.smpl.range[2] <- max(actual.smpl.range[2], smpl.range[2])

    probs <- cumsum(c(0, denFam$pdf(grid[-length(grid)], para[idx, ])) * epsh)
    err.cdf[idx] <- abs(max(probs) - 1)

    probs = probs / max(probs)
    cdf = approxfun(x = grid, y = probs, yleft = 0, yright = 1)
    ks.res = ks.test(x = smpl, y = "cdf")
    ks.p[idx] <- ks.res$p.value
  }

  test_that(name.fam,{
    # testing parameter range
    expect_true(
      !flg.par.range,
      info = paste(
        name.fam,
        'parameter out of range.'
      )
    )
    # testing for cdf range 0 -> 1
    if(mean(err.cdf) >= 1e-2){
      warning(paste0(name.fam, ', large numerical error in cdf, mean error is ', mean(err.cdf)))
    }
    # testing for sample range
    expect_true(
      !flg.smpl.range,
      info = paste(
        name.fam,
        'sample out of range, specified as',
        toString(denFam$range),
        'actual as',
        toString(actual.smpl.range)
      )
    )
    # testing for random generator.
    expect_true(mean(ks.p)>tol, info = paste(name.fam, 'mean p-val as', mean(ks.p)))
    # checking analytic derivatives if provided.
    expect_true(
      is.null(checkDenFamGrident(denFam, tol = 1e-3)),
      info = paste(
        name.fam,
        'error in analytic gradient.'
      )
    )
    # checking numeric issues (if NaN or Inf in neglogll and its gradient)
    numeric.iss <- checkDenFamNumeric(denFam)
    if(!is.null(numeric.iss)){
      warning(sprintf(
        '%s, numeric issues, NaN/Inf in neglogll and/or gradient.',
        name.fam
      ))
      numeric.iss
    }
    expect_true(
      is.null(numeric.iss)
    )
  })
}

test_that('MixGauss only', {
  den.fam <- tMixGauss(5)
  dat <- genSimData(
    den.fam,
    n.train = 10, train.size = 1,
    n.test = 1, test.size = 1
  )
  # checking par repar
  re.par <- list()
  re.re.par <- array(0, dim = dim(dat$train.par))
  for(i in seq_len(nrow(dat$train.par))){
    re.par[[i]] <- tMixGauss.repar(work.par = dat$train.par[i, ])
    re.re.par[i, ] <- tMixGauss.repar(ls.par = re.par[[i]])
  }
  expect_equal(
    dat$train.par,
    re.re.par,
    info = 'tMixGauss.repar: large translation error'
  )
})

#
# # test for wrapped family
# for(i in 1:length(lst_den)){
# test_that(paste0('wrapped ', lst_name[i]),{
#   denFam = wrapDenFam(lst_den[[i]])
#   set.seed(100)
#   para = denFam$gen_par(1)
#   grid = seq(denFam$range[1], denFam$range[2], by = 0.0005)
#   smpl = denFam$rpdf(50, para)
#   # empirical = density(smpl,
#   #                     from = (denFam$range)[1], to = denFam$range[2],
#   #                     bw = 'SJ', kernel = "triangular")
#   probs = cumsum(c(0, denFam$pdf(grid[-length(grid)], para))*0.0005)
#   # testing for cdf range 0 -> 1
#   if(abs(max(probs) - 1) >= 1e-2){
#     warning(paste0('wrapped ', lst_name[i], ', large numerical error in cdf, max is ', max(probs)))
#   }
#
#   probs = probs / max(probs)
#   cdf = approxfun(x = grid, y = probs, yleft = 0, yright = 1)
#   # cdf = splinefun(x = grid, y = cumsum(c(0, denFam$pdf(grid[-length(grid)], para))*0.0005), method = 'natural')
#
#   # plot(cdf)
#
#   # ks.test(x = smpl, y = "cdf")
#
#   # theoCdf = function(x){pcauchy(x, location = para[1], scale = exp(para[2]))/(sum(c(1,-1)*pcauchy(c(10,-10), location = para[1], scale = exp(para[2]))))}
#   # plot(x = grid, y = cumsum(c(0, denFam$pdf(grid[-length(grid)], para))*0.05), type = 'l')
#   # points(x = grid, y = pcauchy(grid, location = para[1], scale = exp(para[2])), type = 'l', col = 'red')
#
#   # ks.test(x = smpl, y = "cdf")$p.value
#   # ks.test(x = smpl, y = 'theoCdf')
#   # ks.test(x = smpl, y = pcauchy)
#   # ks.test(x = rcauchy(10000, 0,1), y = pcauchy)
#   #
#   # qqplot(smpl, rcauchy(500000, location = para[1], scale = exp(para[2])))
#   #
#   # theoretical = denFam$pdf(empirical$x, para)
#   # plot(empirical)
#   # points(y = theoretical, x = empirical$x, type = 'l', col = 'red')
#   # max(abs(empirical$y - theoretical))
#
#   ks.res = ks.test(x = smpl, y = "cdf")
#   # cat(c(lst_name[i], "\n"))
#   # print(range(environment(cdf)$y))
#   # print(ks.res)
#   expect_true(ks.res$p.value>tol, info = paste0('wrapped ', lst_name[i], ' p-val = ', ks.res$p.value))
#   })
# }
#
# context('Testing truncated family.')
# tol = 0.05
# lst_edge = c(1,2,5,7,10,20)
# lst_den = list(tNrml, tCauchy, tMixGauss)
# lst_name = c('truncatedNormal', 'truncatedCauchy', 'truncatedMixGauss')
# for(i in 1:length(lst_den))
#   for(edge in lst_edge){{
#     test_that(paste0(lst_name[i], ' truncated at ', edge),{
#         denFam = lst_den[[i]](edge)
#         for(r.sd in 1:5){
#           set.seed(r.sd)
#           para = denFam$gen_par(1)
#           grid = seq(denFam$range[1], denFam$range[2], by = 0.0005)
#           smpl = denFam$rpdf(10000, para)
#
#
#           probs = cumsum(c(0, denFam$pdf(grid[-length(grid)], para))*0.0005)
#           # testing for cdf range 0 -> 1
#           # test_that(paste0(lst_name[i], ' CDF -> 1.'),
#           #           expect_true(abs(max(probs) - 1) < 5e-2, info = paste0(lst_name[i], ', large numerical error in cdf, max is ', max(probs)))
#           # )
#           if(abs(max(probs) - 1) >= 1e-2){
#             warning(paste0(lst_name[i], ', large numerical error in cdf, max is ', max(probs)))
#           }
#
#           probs = probs / max(probs)
#           cdf = approxfun(x = grid, y = probs, yleft = 0, yright = 1)
#
#           ks.res = ks.test(x = smpl, y = "cdf")
#           expect_true(ks.res$p.value>tol, info = paste0(lst_name[i], ' ', edge, ' p-val = ', ks.res$p.value))
#         }
#     })
# }}
#
# # test for wrapped family
# for(i in 1:length(lst_den))
#   for(edge in lst_edge){{
# test_that(paste0('wrapped ', lst_name[i], ' truncated at ', edge),{
#     denFam = wrapDenFam(lst_den[[i]](edge))
#     set.seed(100)
#     para = denFam$gen_par(1)
#     grid = seq(denFam$range[1], denFam$range[2], by = 0.0005)
#     smpl = denFam$rpdf(50, para)
#
#
#     probs = cumsum(c(0, denFam$pdf(grid[-length(grid)], para))*0.0005)
#     # testing for cdf range 0 -> 1
#     # test_that(paste0(lst_name[i], ' CDF -> 1.'),
#     #           expect_true(abs(max(probs) - 1) < 5e-2, info = paste0(lst_name[i], ', large numerical error in cdf, max is ', max(probs)))
#     # )
#     if(abs(max(probs) - 1) >= 1e-2){
#       warning(paste0('wrapped ', lst_name[i], ', large numerical error in cdf, max is ', max(probs)))
#     }
#
#     probs = probs / max(probs)
#     cdf = approxfun(x = grid, y = probs, yleft = 0, yright = 1)
#
#
#     ks.res = ks.test(x = smpl, y = "cdf")
#     expect_true(ks.res$p.value>tol, info = paste0(lst_name[i], ' ', edge, ' p-val = ', ks.res$p.value))
#     })
# }}


###############################################################################
context('Testing m-e transform.')
lst_name <- names(lst_den)
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
      max(relative.err) < 5e-3,
      info = paste('large m<->e translate error, max', max(relative.err))
    )
  })
}

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
      max(relative.err) < 5e-3,
      info = paste('large m<->e translate error, max', max(relative.err))
    )
  })

}
