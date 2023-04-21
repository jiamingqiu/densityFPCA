context('Testing cv_routine.R')

test_that('gridSearch.f', {
  # testing scale returns
  tst.f <- function(a, b) (a - b)^2
  gs.f <- gridSearch.f(tst.f, over = list(b = seq_len(10)))
  expect_equal(
    gs.f(2),
    list(
      value = 0, b = 2
    )
  )
  tst.f <- function(a, b){
    if(a != b){
      (a - b)^2
    }else{
      simpleError('a = b')
    }
  }
  gs.f <- gridSearch.f(tst.f, over = list(b = seq_len(10)))
  expect_warning(
    expect_equal(
      gs.f(2),
      list(
        value = 1, b = 1
      )
    )
    # 'fun does not return list contain what, candidate dropped at b = 2'
  )

  # testing list returns.
  tst.f <- function(a, b){
    list(
      la = (a - b) ^ 2
    )
  }
  gs.f <- gridSearch.f(tst.f, over = list(b = seq_len(10)), what = 'la')
  expect_equal(
    gs.f(2),
    list(la = 0, b = 2)
  )
  # testing errors
  tst.f <- function(a, b){
    simpleError('well')
  }
  gs.f <- gridSearch.f(tst.f, over = list(b = seq_len(10)), what = 'la')
  expect_warning(expect_error(gs.f(2)))

  gs.f <- gridSearch.f(tst.f, over = list(b = seq_len(10)))
  expect_warning(expect_error(gs.f(2)))

  # tst.f <- function(a, b){
  #   list(
  #     value = b,
  #     object = (a - b) ^ 2
  #   )
  # }
  # gs.f <- gridSearch.f(tst.f, over = list(b = seq_len(10)))
  # gs.f(2) # work but not as intended
})

test_that('downHill.f', {
  # testing scale returns
  tst.f <- function(a, b) (a - b)^2
  dh.f <- downHill.f(tst.f, over = list(b = seq_len(10)))
  expect_equal(
    dh.f(2),
    list(
      value = 0, b = 2
    )
  )
  expect_equal(
    dh.f(-1),
    list(
      value = 4, b = 1
    )
  )
  tst.f <- function(a, b){
    if(a != b){
      (a - b)^2
    }else{
      simpleError('a = b')
    }
  }
  dh.f <- downHill.f(tst.f, over = list(b = seq_len(10)))
  expect_warning(
    expect_equal(
      dh.f(2),
      list(
        value = 1, b = 1
      )
    )
    # 'fun does not return list contain what, candidate dropped at b = 2'
  )

  # testing list returns.
  tst.f <- function(a, b){
    list(
      la = (a - b) ^ 2
    )
  }
  dh.f <- downHill.f(tst.f, over = list(b = seq_len(10)), what = 'la')
  expect_equal(
    dh.f(2),
    list(la = 0, b = 2)
  )
  expect_equal(
    dh.f(0),
    list(la = 1, b = 1)
  )

  # testing errors
  tst.f <- function(a, b){
    simpleError('well')
  }
  dh.f <- downHill.f(tst.f, over = list(b = seq_len(10)), what = 'la')
  expect_warning(expect_error(dh.f(2)))

  dh.f <- downHill.f(tst.f, over = list(b = seq_len(10)))
  expect_warning(expect_error(dh.f(2)))

  # tst.f <- function(a, b){
  #   list(
  #     value = b,
  #     object = (a - b) ^ 2
  #   )
  # }
  # dh.f <- downHill.f(tst.f, over = list(b = seq_len(10)))
  # dh.f(2) # work but not as intended
})

test_that('k-fold CV', {
  f <- function(train, test, tune, deg){
    train <- t(matrix(unlist(train), nrow = 2))
    test <- t(matrix(unlist(test), nrow = 2))

    fit <- colMeans(train) + tune
    loss <- mean(abs(t(test) - fit) ^ deg)
    return(loss)
  }
  set.seed(1)
  dat <- matrix(rnorm(500 * 2), ncol = 2)
  k <- 50
  control <- list(n.rep = 10)
  cv.f <- kFoldCV.f(dat, f, k, control)
  expect_equal(
    0,
    optimize(cv.f, c(-1, 1), deg = 2)$minimum,
    info = 'MSE should be minimized at mean.'
  )
})




