context('Testing HausdorffDist.R')

test_that('L2 inner product',{
  n = 5000
  grid = seq(0,1, length.out = n)
  grid2 = seq(0, 1, length.out = n/2)
  expect_equal(dotL2(cos, sin, range = c(-pi/2, pi/2)), 0)
  expect_equal(dotL2(grid^2, rep(1, n), grid = grid), 1/3)
  expect_equal(dotL2(grid^4, rep(1, n), grid = grid^2), 1/3)
  expect_equal(dotL2(function(x){x^2}, function(x){1}, range = c(0,1)), 1/3)
  expect_equal(dotL2(c(1,1), c(1,1), grid = c(0,1)), 1)
})


test_that('Hausdorff Distance in L2', {
  n = 25000
  # tol = pi^3 / n^2 / 12 * 10
  tol = sqrt(.Machine$double.eps)
  tst_grid = seq(-pi/2, pi/2, length.out = n)
  base1 = list(
    function(x){x*sqrt(12/pi^3)},
    function(x){x^2*sqrt(80/pi^5)}
  )
  base1.g = list(
    tst_grid*sqrt(12/pi^3),
    tst_grid^2*sqrt(80/pi^5)
  )

  base2 = list(
    function(x){sqrt(2/pi)*cos(x)}, function(x){sqrt(2/pi)*sin(x)}
  )

  base2.g = list(
    sqrt(2/pi)*cos(tst_grid),
    sqrt(2/pi)*sin(tst_grid)
  )

  expect_equal(matGram(base1, base1, c(-pi/2, pi/2)),         matrix(c(1,0,0,1),2,2), tolerance = tol)
  expect_equal(matGram(base1.g, base1.g, grid = tst_grid),    matrix(c(1,0,0,1),2,2), tolerance = tol)

  expect_equal(matGram(base2, base2, c(-pi/2, pi/2)),         matrix(c(1,0,0,1),2,2), tolerance = tol)
  expect_equal(matGram(base2.g, base2.g, grid = tst_grid),    matrix(c(1,0,0,1),2,2), tolerance = tol)

  expect_equal(distHausdorff(base1, base1, c(-pi/2, pi/2)),      0, tolerance = tol)
  expect_equal(distHausdorff(base1.g, base1.g, grid = tst_grid), 0, tolerance = tol)

})

#
# distHausdorff(base1, base2, c(-pi/2, pi/2))
# distHausdorff(base1.g, base2.g, grid = tst_grid)

#
# tm = function(x){x^3}
# dotL2(tm, tm, range = c(-pi/2, pi/2))
# tm %>% rm
#
# base3 = list(
#   function(x){x^4*sqrt(2304/pi^9)},
#   function(x){x^5*sqrt(11264/pi^11)}
# )
# base3.1 = list(
#   function(x){x^4*sqrt(2304/pi^9)},
#   function(x){-x^5*sqrt(11264/pi^11)}
# )
#
# distHausdorff(base1, base3, range = c(-pi/2, pi/2))
# distHausdorff(base1, base3.1, range = c(-pi/2, pi/2))
# matGram(base1, base3, range = c(-pi/2, pi/2))


test_that('Semi discrete Wasserstein distance',{

  skip('developing')

  set.seed(100)
  tst = rnorm(5000)
  grid = seq(-15, 15, by = 0.001)
  ts1 = semiDscWass(tst, qnorm, p = 2)
  ts2 = semiDscWass(tst, dnorm(grid), grid = grid, p = 2)
  ts3 = semiDscWass(tst, dnorm(grid)*5, grid = grid, p = 2)
  expect_true(ts1 < 0.05, info = 'semiDscWass to qnorm')
  expect_true(ts2 < 0.05, info = 'semiDscWass to dnorm')
  expect_true(ts3 < 0.05, info = 'semiDscWass to dnorm')
  expect_equal(ts1 - ts2, 0, tol = 1e-5)
  expect_equal(ts2 - ts3, 0, info = 'semiDscWass should include re-unify step.')

  # set.seed(100)
  # tst = rcauchy(5000)
  # grid = seq(min(tst), max(tst), length.out = 1e+7)
  # ts1 = semiDscWass(tst, qcauchy, p = 2)
  # ts2 = semiDscWass(tst, dcauchy(grid), grid = grid, p = 2)
  # expect_true(ts1 < 0.05, info = 'semiDscWass to qcauchy')
  # expect_true(ts2 < 0.05, info = 'semiDscWass to dcauchy')
  # expect_equal(ts1 - ts2, 0)

  set.seed(100)
  tst = runif(500)
  grid = seq(0, 1, by = 0.001)
  ts1 = semiDscWass(tst, qunif, p = 2)
  ts2 = semiDscWass(tst, dunif(grid), grid = grid, p = 2)
  expect_true(ts1 < 0.05, info = 'semiDscWass to qunif')
  expect_true(ts2 < 0.05, info = 'semiDscWass to dunif')
  expect_equal(ts1 - ts2, 0, tol = 1e-4)

  # heavy tail distribution need to be truncated to work properly
  set.seed(100)
  tst = rt(5000, 1)
  tst = tst[tst>=-10 & tst<=10]
  # grid = seq(min(tst), max(tst), length.out = 1e+6)
  # grid = sort(unique(c(tst, seq(-1e+3, 1e+3, length.out = 1e+4))))
  grid = seq(-10, 10, length.out = 1000)
  # ts1 = semiDscWass(tst, function(x) qt(x*pt(10, 1) + (1-x) * pt(-10,1), 1), p = 1)
  # qfunc = function(x) qt(x*pt(10, 1) + (1-x) * pt(-10,1), 1)
  # e.cdf = ecdf(tst)
  # e.qfunc = stepfun(x = environment(e.cdf)$y,
  #                   y = (environment(e.cdf)$x)[c(1, 1:length(environment(e.cdf)$x))],
  #                   f = 1,
  #                   right = FALSE)
  # itpl.e.qfunc = approxfun(x = c(0,environment(e.cdf)$y),
  #                          y = c(environment(e.cdf)$x, 10),
  #                          method = 'constant',
  #                          ties = min, yleft = -10, yright = 10)
  #
  # probs = seq(0,1, 1e-6)
  # plot(x = probs, y = e.qfunc(probs), type = 'l')
  # points(x = probs, y = qfunc(probs), type = 'l', col = 'green')
  # plot(x = probs, y = e.qfunc(probs) - qfunc(probs), type = 'l')
  # points(x = probs, y = itpl.e.qfunc(probs) - qfunc(probs), type = 'l', col = 'blue')
  #
  # integrate(function(x){itpl.e.qfunc(x) - qfunc(x)}, 0, 1, stop.on.error = F)
  # 1e-6 * sum(abs(itpl.e.qfunc(probs) - qfunc(probs)))
  #
  # integrate(function(x){e.qfunc(x) - qfunc(x)}, 0.1,0.9, stop.on.error = F)

  ts2 = semiDscWass(tst, as.numeric(dt(grid, 1) / (pt(10,1) - pt(-10,1))), grid = grid, p = 1)
  expect_true(ts2 < 0.05, info = paste0('semiDscWass to dt, ', ts2, ' is not small enough.'))

  # expect_true(ts1 < 0.5, info = 'semiDscWass to qt')
  # expect_true(ts1 - ts2 < 1e-3, info = 'to qt vs to dt')

  set.seed(100)
  tst = rbeta(1000, 2,2)
  grid = seq(min(tst), max(tst), by = 0.001)
  ts1 = semiDscWass(tst, function(x) qbeta(x, 2,2), p = 2)
  ts2 = semiDscWass(tst, dbeta(grid, 2,2), grid = grid, p = 2)
  expect_true(ts1 < 0.01, info = 'semiDscWass to qbeta')
  expect_true(ts2 < 0.01, info = 'semiDscWass to dbeta')
  expect_true(ts1 - ts2 < 1e-4, info = 'to qbeta vs to dbeta')
})


test_that('Wasserstein distance',{

  skip('developing')

  expect_equal(distWass(qnorm, qnorm, p=1), 0)
  expect_equal(distWass(qcauchy, qcauchy, p=2), 0)

  grid = seq(-7,7, length.out = 1e+5)
  grid2 = seq(-7,7, length.out = 5e+4)

  expect_equal(distWass(qnorm, dnorm(grid), p=2, grid = grid), 0)
  expect_equal(distWass(dnorm(grid), qnorm, p=2, grid = grid), 0)
  expect_equal(
    distWass(dnorm(grid), dnorm(grid2), p=2, grid = grid, grid2 = grid2),
    0,
    tolerance = 1e-5
  )
  expect_error(
    distWass(dnorm(grid), dnorm(grid2), p=2, grid = grid2, grid2 = grid),
    regexp = 'Check input length.'
  )

  # tail is very problematic, TBD
  # grid = unique(c(seq(-5000, -10, length.out = 10000), seq(-10,10, length.out = 1e+7), seq(10, 5000, length.out = 10000)))
  # grid2 = seq(-50,50, length.out = 5e+4)
  # distWass(qcauchy, dcauchy(grid), p=2, grid = grid)
  # distWass(dcauchy(grid), dcauchy(grid2), p=2, grid = grid, grid2 = grid2)

  grid = seq(-7,7,length.out = 1e+5)
  grid2 = seq(-7, 7, length.out = 1e+6)
  tq1 = function(x){qnorm(x*pnorm(7, 1, 1) + (1-x) * pnorm(-7, 1,1), 1,1)}
  tq2 = function(x){qnorm(x*pnorm(7, -1, 1) + (1-x) * pnorm(-7, -1,1),-1,1)}
  trueWass = integrate(function(x){abs(tq1(x) - tq2(x))}, 0,1)
  expect_equal(
    distWass(dnorm(grid, 1,1), dnorm(grid, -1,1), grid = grid, p=1),
    trueWass$value,
    tolerance = 1e-5
  )

  expect_equal(
    distWass(
      function(x){qcauchy(x*pcauchy(7) + (1-x)*pcauchy(-7))},
      dcauchy(grid), p=2, grid = grid
    ), 0, tolerance = 1e-6)
  expect_equal(
    distWass(dcauchy(grid), dcauchy(grid2), p=2, grid = grid, grid2 = grid2),
    0, tolerance = 1e-6
  )

  # grid = seq(-20,20, length.out = 1e+5)
  # distWass(dnorm(grid, 1,1), dnorm(grid, -1,1), grid = grid, p=1)
  # tq1 = function(x){qnorm(x*pnorm(20, 1, 1) + (1-x) * pnorm(-20, 1,1), 1,1)}
  # tq2 = function(x){qnorm(x*pnorm(20, 0, 1) + (1-x) * pnorm(-20, 0,1),0,1)}
  # trueWass = integrate(function(x){abs(tq1(x) - tq2(x))}, 0, 1)
  # trueWass
  # cpp_distWass(c(0,pnorm(grid, 1,1)/pnorm(20,1,1)), c(0, pnorm(grid, 0,1)/pnorm(20,0,1)), grid = c(-8, grid), p=1)
  # cpp_distWass(c(0,pnorm(grid, 1,1)), c(0, pnorm(grid, 0,1)), grid = c(-21, grid), p=1)
  #
  #
  #
  # spl1 = c(1,2,3)
  # spl2 = c(1,2,5)
  # cmGrid = sort(unique(c(spl1, spl2)))
  # pdf1 = rep(0, length(cmGrid))
  # pdf2 = rep(0, length(cmGrid))
  # pdf1[is.element(cmGrid, spl1)] = 1/length(spl1)
  # pdf2[is.element(cmGrid, spl2)] = 1/length(spl2)
  # cdf1 = cumsum(pdf1)
  # cdf2 = cumsum(pdf2)
  # cpp_distWass(cdf1, cdf2, cmGrid, 1)
  # sum(abs(sort(spl1) - sort(spl2))  ) / length(spl1)
  #
  # spl1 = rnorm(10)
  # spl2 = rnorm(5)
  # cmGrid = sort(unique(c(spl1, spl2)))
  # pdf1 = rep(0, length(cmGrid))
  # pdf2 = rep(0, length(cmGrid))
  # pdf1[is.element(cmGrid, spl1)] = 1/length(spl1)
  # pdf2[is.element(cmGrid, spl2)] = 1/length(spl2)
  # cdf1 = cumsum(pdf1)
  # cdf2 = cumsum(pdf2)
  # cpp_distWass(c(0,cdf1), c(0,cdf2), c(cmGrid[1] - 1, cmGrid), 1)
  # sum(abs(cdf1 - cdf2)[-length(cdf1)] * diff(cmGrid))
  #
  # # so as to demonstrate the algorithm in cpp_distWass is ok for point mass
  # # yet one need to rewrite it so that it take cdf, in which case let R
  # # handle pdf --> cdf to enable more flexibility such as handling semiDsc or
  # # both point mass or etc etc.
  # set.seed(1)
  # spl1 = rnorm(50)
  # spl2 = rnorm(100)
  # cmGrid = sort(unique(c(spl1, spl2)))
  # pdf1 = rep(0, length(cmGrid))
  # pdf2 = rep(0, length(cmGrid))
  # pdf1[is.element(cmGrid, spl1)] = 1/length(spl1)
  # pdf2[is.element(cmGrid, spl2)] = 1/length(spl2)
  # sum(abs(sort(spl1) - sort(spl2))  ) / length(spl1)
  #
  # # t.cdf[i] means on [cmGrid[i], cmGrid[i+1]), cdf = t.cdf[i]
  # # cmGrid[maxlength+1] = Inf
  # t.cdf1 = cumsum(pdf1)
  # t.cdf2 = cumsum(pdf2)
  # # padding so that the feed cdf always starts with 0 to its left
  # cdf1 = c(0, t.cdf1)
  # cdf2 = c(0, t.cdf2)
  # grid = c(min(cmGrid) - 1, cmGrid)
  #
  # from = 0
  # to = 0
  # idx1 = 1
  # idx2 = 1
  # res = 0
  # p = 1
  #
  # while(from < 1.0){
  #   # cdf1 jump next?
  #   if(cdf1[idx1 + 1] < cdf2[idx2 + 1]){
  #     to = cdf1[idx1 + 1]
  #
  #     res = res + (to - from) * abs(grid[idx1+1] - grid[idx2+1])^p
  #     from = to
  #     idx1 = idx1 + 1
  #   }
  #   # cdf2 jump next?
  #   if(cdf1[idx1 + 1] > cdf2[idx2 + 1]){
  #     to = cdf2[idx2 + 1]
  #
  #     res = res + (to - from) * (abs(grid[idx1+1] - grid[idx2+1]))^p
  #     from = to
  #     idx2 = idx2 + 1
  #   }
  #   if(cdf1[idx1 + 1] == cdf2[idx2 + 1]){
  #     to = cdf2[idx2 + 1]
  #
  #     res = res + (to - from) * (abs(grid[idx1+1] - grid[idx2+1]))^p
  #     from = to;
  #     idx1 = idx1 + 1
  #     idx2 = idx2 + 1
  #   }
  # }
  #
  # res
  # sum(abs(t.cdf1 - t.cdf2)[-length(t.cdf1)] * diff(cmGrid))
  # cpp_distWass(cdf1, cdf2, grid, 1)

})


