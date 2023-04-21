context('Testing utili.R')

test_that('calISE',{
  n = 10000
  grid = seq(0,1, length.out = n)
  grid2 = seq(0, 1, length.out = n/2)
  f1 = function(x){x^2}
  f2 = function(x){0}
  lst.f1 = grid^2
  lst.f2 = rep(0, n)
  lst.f1.2 = grid2^2
  true = 1/5
  expect_equal(calISE(f1, f2, range = c(0,1)), true)
  expect_equal(calISE(lst.f1, lst.f2, grid = grid),      true)
  expect_equal(calISE(lst.f1.2, lst.f2, grid = grid2, grid2 = grid), true)
})

test_that('calKL',{
  n = 10000
  grid = seq(-5,2, length.out = n)
  grid2 = seq(-5, 2, length.out = n/2)
  lst.std = dnorm(grid)
  lst.std.g2 = dnorm(grid2)
  lst.sps = dnorm(grid, sd = 10)
  lst.sps.g2 = dnorm(grid2, sd = 10)
  t2 = integrate(function(x){
    dnorm(x) / (pnorm(2) - pnorm(-5)) * log(dnorm(x)/(pnorm(2) - pnorm(-5)) /
                     dnorm(x, sd = 10) * (pnorm(2, sd = 10) - pnorm(-5, sd = 10))
    )
  }, -5, 2)$value

  # this is the un-normalized value, cal. by Mathematica
  # 1.0*Integrate[
  #   PDF[NormalDistribution[0, 1], x]*
  #     Log[PDF[NormalDistribution[0, 1], x]/
  #           PDF[NormalDistribution[0, 10], x]], {x, -5, 2}]
  u.t2 = 1.81991651195262

  expect_equal(calKL(dnorm, dnorm, range = c(-5,2)), 0)
  expect_equal(calKL(dnorm, function(x){dnorm(x, sd = 10)}, range = c(-5,2)), u.t2)

  expect_equal(calKL(lst.std, lst.std, grid = grid),0)
  expect_equal(calKL(lst.std, lst.std.g2, grid = grid, grid2 = grid2),0)
  expect_equal(calKL(lst.std.g2, lst.std, grid = grid2, grid2 = grid),0)

  expect_equal(calKL(lst.sps, lst.sps, grid = grid),0)
  expect_equal(calKL(lst.std, lst.sps, grid = grid), t2)
  expect_equal(calKL(lst.std, lst.sps.g2,    grid = grid,  grid2 = grid2),t2)
  #expect_equal(calKL(lst.std.g2, lst.sps,    grid = grid2, grid2 = grid ),t2)
  expect_equal(calKL(lst.std.g2, lst.sps.g2, grid = grid2, grid2 = grid2),t2, tol = 1e-7)

  expect_warning(calKL(rep(1,1000), rep(1,1000),grid = seq(0,1,length.out = 1000), grid2 = seq(0.5, 1.5, length.out = 1000)),
                 regexp = 'No absolute continuity.')

  # this test does not make sense mathematically, just in case
  expect_equal(calKL(rep(0,1000), rep(1, 1000), grid = seq(0,1, length.out = 1000)), 0)
})

test_that('calFR',{

  n = 30000
  g1 = seq(0, pi/2, length.out = n)
  g2 = seq(0, pi/2, length.out = n/2)
  f1 = function(x){sin(x)^2}
  f2 = function(x){cos(x)^2}
  expect_equal(calFR(f1, f2, range = c(0, pi/2)), acos(2/pi))
  expect_equal(calFR(f1(g1), f2(g1), grid = g1), acos(2/pi))
  expect_equal(calFR(f1(g1), f2(g2), grid = g1, grid2 = g2), acos(2/pi))

})


test_that('assErr', {
  grid.size <- 512
  grid <- seq(5, 6, length.out = grid.size)
  actual <- t(matrix(c(grid, grid^2, grid^3), ncol = 3))
  set.seed(1)
  esti <- actual + rnorm(length(actual), sd = 0.5)
  ass.res <- assErr(
    actual, esti, grid = grid,
    err.method = c('FisherRao', 'ISE', 'KL', 'rev.KL')
  )
  rec.res <- rep(0, nrow(ass.res$err))
  for(i in seq_along(rec.res)){
    idx.err <- ass.res$err$type.err[i]
    idx.row <- ass.res$err$idx[i]
    true.pdf <- actual[idx.row, ]
    esti.pdf <- esti[idx.row, ]
    if(idx.err == 'FisherRao'){
      rec.res[i] <- calFR(
        f1 = true.pdf,
        f2 = esti.pdf,
        grid = grid
      )
    }
    if(idx.err == 'ISE'){
      rec.res[i] <- calISE(
        f1 = true.pdf,
        f2 = esti.pdf,
        grid = grid
      )
    }
    if(idx.err == 'KL'){
      rec.res[i] <- calKL(
        f1 = true.pdf,
        f2 = esti.pdf,
        grid = grid
      )
    }
    if(idx.err == 'rev.KL'){
      rec.res[i] <- calKL(
        f1 = esti.pdf,
        f2 = true.pdf,
        grid = grid
      )
    }
  }
  expect_true(
    all(rec.res == ass.res$err$err),
    'label mismatched in wrapper assErr.'
  )
})
