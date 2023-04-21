# Inspecting results ==========================================================

alignCurve <- function(mat, grid, by = 2){
  # aligning the matrix of curves using the 1st row/col as reference s.t.
  # curves dotL2 ref > 0, by mltiplying -1 (or not)
  # args:
  # mat: matrix of curves,
  # grid: the grid,
  # by: treat 1 for row, 2 for column as one curve.
  # return:
  # a matrix of same dimension.
  if(by == 1) mat <- t(mat)

  res <- mat
  for(i in seq_len(ncol(mat))){
    dot.res <- dotL2(mat[, 1], mat[, i], grid = grid)
    if(dot.res < 0)
      res[, i] <- -1 * res[, i]
  }
  if(by == 1)
    return(t(res))

  return(res)
}


pltPrior <- function(fpca.res, mat.esti, idx.df, num.k){
  # visual inspection of prior and potential estimation. (2-d only)
  # args:
  # fpca.res: result from fdapace::FPCA
  # mat.esti and idx.df: optional, a matrix and indexing data frame of the
  #                      estimated parameters. idx.df$esti.method is used for
  #                      colour, idx.df$idx.obsv used for label.
  # num.k: number of components to consider, default = ncol(fpca.res$xiEst)
  # return: a ggplot2::ggplot object

  # ### DEBUG
  #   fpca.res <- list(
  #     xiEst = matrix(rnorm(500, sd = 2), ncol = 4)
  #   )
  #   mat.esti <- matrix(rnorm(100), ncol = 2)
  #   idx.df <- data.frame(
  #     esti.method = rep(c('m.a', 'm.b'), each = 25),
  #     idx.obsv = rep(seq_len(25), 2)
  #     )
  #   num.k <- 3
  # ### END DEBUG

  if(!missing(mat.esti) & !missing(idx.df))
    if(nrow(mat.esti) != nrow(idx.df))
      stop('check input rows: mat.esti and idx.df.')
  if(missing(num.k))
    num.k <- ncol(fpca.res$xiEst)

  fpca.df <- data.frame(
    xiEst = fpca.res$xiEst[, seq_len(num.k), drop = FALSE],
    esti.method = 'prior_xiEst'
  )
  names(fpca.df) <- c(
    paste0('xiEst.', seq_len(num.k)),
    'esti.method'
  )

  if(missing(idx.df))
    idx.df <- data.frame(esti.method = NULL)
  if(missing(mat.esti)){
    esti.df <- data.frame(
      xiEst = matrix(0, ncol = num.k, nrow = 0)
    )
  }else{
     # if(num.k > ncol(mat.esti))
     #   stop('check input col: fpca.res$xiEst and mat.esti.')
    colnames(mat.esti) <- NULL
    esti.df <- data.frame(
      xiEst = mat.esti[, seq_len(min(ncol(mat.esti), num.k)), drop = FALSE]
    )
  }

  if(is.null(idx.df$esti.method)){
    esti.df$esti.method <- 'estimates'
  }else{
    esti.df$esti.method <- idx.df$esti.method
  }
  if(is.null(idx.df$idx.obsv)){
    esti.df$idx.obsv <- NA
    flg.text <- FALSE
  }else{
    esti.df$idx.obsv <- idx.df$idx.obsv
    flg.text <- TRUE
  }

  plt.df <- suppressWarnings(dplyr::bind_rows(
    fpca.df,
    esti.df
  ))

  # some local temp functions
  noFillDensity <- function(data, mapping){
    w.mapping <- mapping
    w.mapping$shape <- NULL
    w.mapping$fill <- w.mapping$colour
    ggplot2::ggplot(
      data = data,
      mapping = w.mapping
    ) +
      ggplot2::stat_density(
        alpha = 0.4, position='identity',
        bw = 'SJ'
        )
  }

  tmpLower.text <- function(data, mapping){
    alt.mapping <- mapping
    alt.mapping$label <- NULL
    res <- ggplot2::ggplot(
      data = dplyr::filter(data, esti.method != 'prior_xiEst'),
      mapping = mapping
    ) +
    ggplot2::geom_text(size = 5, alpha = 0.5) +
    ggplot2::geom_density2d(
      data = dplyr::filter(data, esti.method == 'prior_xiEst'),
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = dplyr::filter(data, esti.method == 'prior_xiEst'),
      mapping = alt.mapping,
      show.legend = FALSE,
      size = 0.75
    )
    res
    # plotly::ggplotly(res)
  }

  tmpLower.pts <- function(data, mapping){
    res <- ggplot2::ggplot(
      data = dplyr::filter(data, esti.method != 'prior_xiEst'),
      mapping = mapping
    ) +
      ggplot2::geom_point(size = 1.75) +
      ggplot2::geom_density2d(
        data = dplyr::filter(data, esti.method == 'prior_xiEst'),
        show.legend = FALSE
      ) +
      ggplot2::geom_point(
        data = dplyr::filter(data, esti.method == 'prior_xiEst'),
        mapping = mapping,
        show.legend = FALSE,
        size = 0.75
      )
    res
    # plotly::ggplotly(res)
  }
  if(flg.text){
    res <- GGally::ggpairs(
      plt.df, title = 'xiEst plot',
      columns = seq_len(num.k),
      mapping = ggplot2::aes(shape = esti.method, color = esti.method, label = idx.obsv),
      lower = list(continuous = tmpLower.text),
      diag = list(
        continuous = noFillDensity
      )
    )
  }else{
    res <- GGally::ggpairs(
      plt.df, title = 'xiEst plot',
      columns = seq_len(num.k),
      mapping = ggplot2::aes(shape = esti.method, color = esti.method),
      lower = list(continuous = tmpLower.pts),
      diag = list(
        continuous = noFillDensity
      )
    )
  }
  return(res)

  # res.plt <- ggplot2::ggplot(
  #   plt.df,
  #   ggplot2::aes(x = x, y = y, shape = esti.method, col = esti.method)
  # ) +
  #   ggplot2::geom_point(size = 2) +
  #   ggplot2::geom_density2d(
  #     data = fpca.df,
  #     aes(x = x, y = y),
  #     show.legend = FALSE
  #   )
  # return(res.plt)
}


# Merging simulation results from multiple batches ============================

mergeSimRes = function(res_path, simple = T){
  # dir is the dir where all folders of individual batch of simulation results are

  ls_batch = list.dirs(path = res_path, recursive = F)
  if(length(ls_batch) == 0) stop('No subfolder.')

  ls.res = list()
  ls.seed = list()
  ls.simSettings = list()

  ls.err = list()
  ls.duration = list()

  idx.batch = 1
  for(path_batch in ls_batch){
    rm(simRes)
    invisible(gc())
    if(memory.size() > 0.975 * memory.limit()){
      stop('Insufficient RAM.')
    }

    cat(paste0('Reading batch ', idx.batch, '.\n'))
    load(paste0(path_batch, '/simRes.RData'))
    load(paste0(path_batch, '/seed.RData'))

    if(length(simRes) != nrow(simSettings)) stop('result list length differ from settings.')
    row.names(simSettings) = NULL

    cat(paste0('Processing batch ', idx.batch, '.\n'))
    Err = lapply(simRes, "[[", "Err")
    duration = lapply(simRes, "[[", "duration")
    for(i in 1:length(Err)){
      cat(paste0('Setting: ', i, ';\n'))
      tmp_simCase = simSettings[rep(i, nrow(Err[[i]])),]
      Err[[i]] = cbind(tmp_simCase, Err[[i]])
      Err[[i]]$idx.batch = idx.batch
      # Err[[i]]$idx.train = paste0(Err[[i]]$idx.train, '.', )
      Err[[i]]$seed = lst_seed[i]
      duration[[i]] = as.data.frame(duration[[i]])
      duration[[i]]$step = row.names(duration[[i]])
      duration[[i]]$idx.batch = idx.batch
      duration[[i]] = cbind(simSettings[rep(i, nrow(duration[[i]])),], duration[[i]])
    }
    cat('combining...\n')
    ls.err[[idx.batch]] = dplyr::bind_rows(Err)
    ls.duration[[idx.batch]] = dplyr::bind_rows(duration)
    # Err = dplyr::select(Err, c(colnames(simSettings), 'method.esti', 'type.err', 'err'))
    row.names(ls.err[[idx.batch]]) = NULL
    row.names(ls.duration[[idx.batch]]) = NULL
    if(!simple){
      if(memory.size() > 0.85 * memory.limit()){
        warning('Result too large for fully loaded into RAM, force simple = T.')
        simple = T
        idx.batch = idx.batch + 1
        rm(ls.res)
        gc()
        next
      }
      ls.res[[idx.batch]] = simRes
      ls.simSettings[[idx.batch]] = simSettings
      ls.seed[[idx.batch]] = lst_seed
    }
    idx.batch = idx.batch + 1
  }
  if(simple){
    return(list(Err = do.call(rbind, ls.err),
                duration = do.call(rbind, ls.duration)))
  }else{
    return(list(Err = do.call(rbind, ls.err),
                duration = do.call(rbind, ls.duration),
                simRes = ls.res,
                simSettings = ls.simSettings,
                seed = ls.seed))
  }
}

# res_path = '../results/sim25fin'
# mergeRes = mergeSimRes(res_path = res_path)
# save(mergeRes, file = paste0(res_path, '/merged.RData'))
# load(paste0(res_path,'/merged.RData'))

### summerising results ===============================================

# TBD: need update to accommodate new Err output
interpretSimRes = function(mergeRes, printSummary = FALSE, plotISE = FALSE){
  # producing readable result
  # mergRes should be produced by mergeSimRes
  # For now, this function does not produce return value.
  err.dt = dplyr::group_by(mergeRes$Err, den_fam, method.esti, type.err)
  t.dt = dplyr::group_by(mergeRes$duration, den_fam, step)

  if(printSummary){
    cat('Error ========================\n')
    print(
      dplyr::summarise(err.dt,
                     min = min(err, na.rm = T),
                     Q1 = quantile(err, 0.25, na.rm = T),
                     mean = mean(err, na.rm = T),
                     median = median(err, na.rm = T),
                     Q3 = quantile(err, 0.75, na.rm = T),
                     max = max(err, na.rm = T),
                     eff.n = sum(!is.na(err))
    ), n=Inf)
    cat('Duration =====================\n')
    print(
      dplyr::summarise(t.dt,
                     min = min(elapsed),
                     Q1 = quantile(elapsed, 0.25),
                     mean = mean(elapsed),
                     median = median(elapsed),
                     Q3 = quantile(elapsed, 0.75),
                     max = max(elapsed)
    ), n=Inf)
  }

  if(plotISE){
    ls_type.err = as.character(unique(mergeRes$Err$type.err))
    ls_bx_plt = list()
    i=1
    for(t.err in ls_type.err){
      plt_df = dplyr::filter(mergeRes$Err, type.err == t.err)
      plt_df = dplyr::summarise(
        dplyr::group_by(plt_df,
                 den_fam, trainSize, testSize, gridSize, method.esti, idx.train, idx.batch
                 ),
        mean.err = mean(err, na.rm = T)
      )
      plt.up = quantile(plt_df$mean.err, probs = 0.95, na.rm = T)
      ls_bx_plt[[i]] = ggplot2::ggplot(plt_df, ggplot2::aes(x = den_fam, y = mean.err, fill = method.esti)) +
        ggplot2::geom_boxplot() +
        ggplot2::coord_cartesian(ylim = c(0, plt.up)) +
        ggplot2::ggtitle(paste0('Type.err:  ', t.err))
      i=i+1
    }
    ls_bx_plt[['ncol']] = 2
    plt = do.call(cowplot::plot_grid, ls_bx_plt)
    print(plt)
  }

  ### TDS
  # error handling in simRes

  ### some testing
  # do.call is better
  # microbenchmark::microbenchmark(
  #   do.call(rbind, ISE),
  #   t(matrix(unlist(ISE), nrow = ncol(simRes[[1]]$Err$ISE)))
  # )

}

# res_path = '../results/test_merge'
# interpretSimRes(mergeSimRes(res_path), T, T)

# inspect = function(simRes, res_path, ...){
# # simRes should be produced by mergeSimRes
# # if res_path not missing, mergeSimRes will be called.
#
#   if(!missing(simRes)){
#     dat = simRes
#   }else{
#   if(!missing(res_path)){
#     dat = mergeSimRes(res_patch, simple = T)
#   }else{
#     stop('Check input argument.')
#   }}
#
# }
