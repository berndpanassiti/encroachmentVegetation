# Function for plotting metaMDS objects from package vegan
# x                     metaMDS object
# dims                  which axes to be plotted
# var                   logical: whether variables are shown
# varlab                logical: whether labels should be plotted for
#                       variables
# obslab                logical: whether labels should be plotted for
#                       observations
# obscolor              char colour or vector with colours for 
#                       observations
# obssize               size of obs points or vectors with values that
#                       are used for scaling the points
# color_legend_title    title for the color legend
# size_legend_title     title for the size legend
# surf
# theme
# ...                   passed to sub-functions

autoplot.metaMDS <- function(x,
                             dims = c(1, 2),
                             var = FALSE,
                             varlab = FALSE,
                             obslab = FALSE,
                             obscolor = NULL,
                             obssize = NULL,
                             color_legend_title = "Groups",
                             size_legend_title = NULL,
                             surf = NULL, # beta
                             theme = theme_minimal(),
                             ...) {
  require(vegan)
  require(ggplot2)
  require(ggrepel)

  # Observation scores
  mdsobs <- as.data.frame(vegan::scores(x, display = "sites"))[, dims]
  
  # Variable scores (in case NMDS was not based on a distance object)
  if ("matrix" %in% class(x$species)) {
    mdsvar <- as.data.frame(vegan::scores(x, display = "species"))[, dims]
  } else {
    var <- FALSE
    varlab <- FALSE
    mdsvar <- NULL
  }

  namx <- colnames(mdsobs)[1]
  namy <- colnames(mdsobs)[2]

  # Build a color vector
  # Function to check if all elements are valid colors
  is_valid_color <- function(col) {
    if (is.numeric(col)) {
      return(FALSE)
    }
    tryCatch({
      grDevices::col2rgb(col)
      TRUE
    }, error = function(e) {
      FALSE
    })
  }
  # Add color to mdsobj
  if (is.null(obscolor)) {
    mdsobs$obscolor <- "#333333"
  } else if (is_valid_color(obscolor)) {
    mdsobs$obscolor <- obscolor
  } else { 
    mdsobs$obscolor <- as.factor(obscolor)
  }
  
  collegend <- ifelse(length(unique(mdsobs$obscolor)) == 1, FALSE, TRUE)

  # Add point size to mdsobs if provided
  if (!is.null(obssize)) {
    mdsobs$obssize <- obssize
    sizlegend <- ifelse(length(unique(mdsobs$obssize)) == 1, FALSE, TRUE)
  } else {
    sizlegend <- FALSE
  }

  # Legend titles
  if (is.null(color_legend_title)) {
    color_legend_title <- deparse(substitute(obscolor))
  }

  g <- ggplot(mdsobs, aes(x = get(namx), y = get(namy)))


  if (collegend && sizlegend) {
    g <- g + geom_point(aes(colour = obscolor, size = obssize), alpha = 0.7)
  } else if (collegend) {
    g <- g + geom_point(aes(colour = obscolor), alpha = 0.7)
  } else if (sizlegend) {
    g <- g + geom_point(aes(size = obssize), alpha = 0.7)
  } else {
    g <- g + geom_point(colour = mdsobs$obscolor, alpha = 0.7)
  }

  g <- g + labs(x = namx, y = namy, colour = color_legend_title, size = size_legend_title) + coord_fixed(ratio = 1)

  if (var) {
    g <- g + geom_point(data = mdsvar,
                        aes(x = mdsvar[, 1], y = mdsvar[, 2]),
                        colour = "darkgrey", shape = 3, inherit.aes = FALSE)
  }

  if (varlab) {
    g <- g + ggrepel::geom_text_repel(data = mdsvar,
                                      aes(x = mdsvar[, 1], y = mdsvar[, 2], label = rownames(mdsvar)),
                                      colour = "red", show.legend = FALSE, inherit.aes = FALSE, ...)
  }

  if (obslab) {
    g <- g + ggrepel::geom_text_repel(data = mdsobs,
                                      aes(x = get(namx), y = get(namy), label = rownames(mdsobs)),
                                      show.legend = FALSE, inherit.aes = FALSE, ...)
  }
  
  # Set axis limits to match total ranges
  xlim_range <- range(c(mdsvar[[namx]], mdsobs[[namx]]), na.rm = TRUE)
  ylim_range <- range(c(mdsvar[[namy]], mdsobs[[namy]]), na.rm = TRUE)
  g <- g + xlim(xlim_range) + ylim(ylim_range)

  # Add isolines from vegan::ordisurf
  if (!is.null(surf)) {
    grid_vals <- expand.grid(x = surf$grid$x, y = surf$grid$y)
    grid_vals$z <- as.vector(surf$grid$z)
    g <- g + stat_contour(data = grid_vals, aes(x = x, y = y, z = z, colour = ..level..)) +
      scale_colour_gradient(low = "green4", high = "red")
  }
  
  # Apply the specified theme
  g <- g + theme
  
  print(g)
}






# Function for plotting metaMDS objects from package vegan, S. Schmidtlein
# x                     metaMDS object
# dims                  which axes to be plotted
# var                   logical: whether variables are shown
# varlab                logical: whether labels should be plotted for
#                       variables
# obslab                logical: whether labels should be plotted for
#                       observations
# obscolor              char colour or vector with colours for 
#                       observations
# obssize               size of obs points or vectors with values that
#                       are used for scaling the points
# color_legend_title    title for the color legend
# size_legend_title     title for the size legend
# surf
# theme
# ...                   passed to sub-functions

autoplot.metaMDS <- function(x,
                             dims = c(1, 2),
                             var = FALSE,
                             varlab = FALSE,
                             obslab = FALSE,
                             obscolor = NULL,
                             obssize = NULL,
                             color_legend_title = "Groups",
                             size_legend_title = NULL,
                             surf = NULL, # beta
                             theme = theme_minimal(),
                             ...) {
  require(vegan)
  require(ggplot2)
  require(ggrepel)
  
  # Observation scores
  mdsobs <- as.data.frame(vegan::scores(x, display = "sites"))[, dims]
  
  # Variable scores (in case NMDS was not based on a distance object)
  if ("matrix" %in% class(x$species)) {
    mdsvar <- as.data.frame(vegan::scores(x, display = "species"))[, dims]
  } else {
    var <- FALSE
    varlab <- FALSE
    mdsvar <- NULL
  }
  
  namx <- colnames(mdsobs)[1]
  namy <- colnames(mdsobs)[2]
  
  # Build a color vector
  # Function to check if all elements are valid colors
  is_valid_color <- function(col) {
    if (is.numeric(col)) {
      return(FALSE)
    }
    tryCatch({
      grDevices::col2rgb(col)
      TRUE
    }, error = function(e) {
      FALSE
    })
  }
  # Add color to mdsobj
  if (is.null(obscolor)) {
    mdsobs$obscolor <- "#333333"
  } else if (is_valid_color(obscolor)) {
    mdsobs$obscolor <- obscolor
  } else { 
    mdsobs$obscolor <- as.factor(obscolor)
  }
  
  collegend <- ifelse(length(unique(mdsobs$obscolor)) == 1, FALSE, TRUE)
  
  # Add point size to mdsobs if provided
  if (!is.null(obssize)) {
    mdsobs$obssize <- obssize
    sizlegend <- ifelse(length(unique(mdsobs$obssize)) == 1, FALSE, TRUE)
  } else {
    sizlegend <- FALSE
  }
  
  # Legend titles
  if (is.null(color_legend_title)) {
    color_legend_title <- deparse(substitute(obscolor))
  }
  
  g <- ggplot(mdsobs, aes(x = get(namx), y = get(namy)))
  
  
  if (collegend && sizlegend) {
    g <- g + geom_point(aes(colour = obscolor, size = obssize), alpha = 0.7)
  } else if (collegend) {
    g <- g + geom_point(aes(colour = obscolor), alpha = 0.7)
  } else if (sizlegend) {
    g <- g + geom_point(aes(size = obssize), alpha = 0.7)
  } else {
    g <- g + geom_point(colour = mdsobs$obscolor, alpha = 0.7)
  }
  
  g <- g + labs(x = namx, y = namy, colour = color_legend_title, size = size_legend_title) + coord_fixed(ratio = 1)
  
  if (var) {
    g <- g + geom_point(data = mdsvar,
                        aes(x = mdsvar[, 1], y = mdsvar[, 2]),
                        colour = "darkgrey", shape = 3, inherit.aes = FALSE)
  }
  
  if (varlab) {
    g <- g + ggrepel::geom_text_repel(data = mdsvar,
                                      aes(x = mdsvar[, 1], y = mdsvar[, 2], label = rownames(mdsvar)),
                                      colour = "red", show.legend = FALSE, inherit.aes = FALSE, ...)
  }
  
  if (obslab) {
    g <- g + ggrepel::geom_text_repel(data = mdsobs,
                                      aes(x = get(namx), y = get(namy), label = rownames(mdsobs)),
                                      show.legend = FALSE, inherit.aes = FALSE, ...)
  }
  
  # Set axis limits to match total ranges
  xlim_range <- range(c(mdsvar[[namx]], mdsobs[[namx]]), na.rm = TRUE)
  ylim_range <- range(c(mdsvar[[namy]], mdsobs[[namy]]), na.rm = TRUE)
  g <- g + xlim(xlim_range) + ylim(ylim_range)
  
  # Add isolines from vegan::ordisurf
  if (!is.null(surf)) {
    grid_vals <- expand.grid(x = surf$grid$x, y = surf$grid$y)
    grid_vals$z <- as.vector(surf$grid$z)
    g <- g + stat_contour(data = grid_vals, aes(x = x, y = y, z = z, colour = ..level..)) +
      scale_colour_gradient(low = "green4", high = "red")
  }
  
  # Apply the specified theme
  g <- g + theme
  # Add stress value as annotation
  g <- g + annotate("text", x = Inf, y = -Inf, 
                    label = paste0("Stress = ", round(x$stress, 3)), 
                    hjust = 1.1, vjust = -0.5, size = 4)
  
  print(g)
}

