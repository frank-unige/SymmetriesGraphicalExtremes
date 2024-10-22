## Plotting function from the code to the paper "Graphical models for multivariate extremes" of Sebastian Engelke, 
## Manuel Hentschel, Michaël Lalancette and Frank Röttger, see https://arxiv.org/abs/2402.02187

MAX_N_EDGES_TICKS <- 7

plot_fitted_params <- function(G0, G1, xlab = 'Empirical', ylab = 'Fitted', colors = NULL, isFirst = TRUE){
  x <- G0[upper.tri(G0)]
  y <- G1[upper.tri(G1)]
  if(is.matrix(colors)){
    stopifnot(identical(dim(colors), dim(G0)))
    colors <- colors[upper.tri(colors)]
    # Order so that colored points are plotted "on top"
    ord <- order(colors, na.last = FALSE)
    colors <- colors[ord]
    x <- x[ord]
    y <- y[ord]
  } else if(!is.null(colors)){
    stop('Argument `colors` needs to be a matrix or NULL.')
  }
  rng <- range(x, y)
  ggp <- (
    ggplot()
    + geom_point(aes(
      x = x,
      y = y,
      col = colors
    ))
    + scale_color_discrete(na.value = 'black')
    # + scale_colour_discrete()
    + geom_abline(slope = 1, intercept = 0)
    + xlab(xlab)
    + ylab(ylab)
    + xlim(rng)
    + ylim(rng)
    + theme(legend.position = "none")
  )
  if(!isFirst){
    ggp <- ggp + theme(
      # axis.ticks.y = element_blank(),
      axis.title.y = element_blank()
      # axis.text.y = element_blank()
    )
  }
  return(ggp)
}

makeColorMat <- function(graph, parti, char = TRUE){
  d <- igraph::vcount(graph)
  edges <- igraph::ends(graph, igraph::E(graph))
  ret <- matrix(NA, d, d)
  for(i in seq_along(parti)){
    for(j in parti[[i]]){
      xy <- edges[j,]
      ret[xy[1], xy[2]] <- i
      ret[xy[2], xy[1]] <- i
    }
  }
  if(char){
    ret[] <- as.character(ret[])
  }
  return(ret)
}