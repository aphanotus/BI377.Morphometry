# BI347 January 2026

# Template


# Copy the data file at the URL below to your labtop
# https://github.com/aphanotus/BI377.Morphometry/raw/refs/heads/main/wing.silhouettes.zip
# Then upload the zip file to your user folder on Rstudio. It will automatically unzip
# to create a folder called `wing.silhouettes`
image.folder <- "wing.silhouettes"

# Set your working directory to one level down from this folder
# If you've never messed with working directories this should be `~`

?setwd

# Load packages
library(stringr)
library(magrittr)
library(Momocs)
library(RRPP)


# Defining some useful functions.
{
  review.outlines <- function(x) {
    number.of.screens <- ceiling(length(x)/4)
    par(mfrow = c(2, 2))
    if (is(x)[1] == "list") {
      coo.count <- unique(sort(unlist(lapply(x, dim))))[-1]
      message(paste0("Coordinate number min: ",min(coo.count), "  median: ",median(coo.count), "  max: ",max(coo.count)))
      for (i in 1:number.of.screens) {
        j <- i*4 -3
        plot(x[[j]], main = paste0(j,": ",names(x)[j]), asp=c(1,1), xlab = "x", ylab="y")
        j <- i*4 -2
        if (j <= length(x)) { plot(x[[j]], main = paste0(j,": ",names(x)[j]), asp=c(1,1), xlab = "x", ylab="y") }
        j <- i*4 -1
        if (j <= length(x)) { plot(x[[j]], main = paste0(j,": ",names(x)[j]), asp=c(1,1), xlab = "x", ylab="y") }
        j <- i*4 -0
        if (j <= length(x)) { plot(x[[j]], main = paste0(j,": ",names(x)[j]), asp=c(1,1), xlab = "x", ylab="y") }
        invisible(readline(prompt=paste("Screen",i,"of",number.of.screens," - Press [ENTER] to continue ([ESC] to quit)")))
      }
    } else {
      if (is(x)[1] == "Out") {
        coo.count <- unique(sort(unlist(lapply(x$coo, function(x) {dim(as.data.frame(x))}))))[-1]
        message(paste0("Coordinate number min: ",min(coo.count), "  median: ",median(coo.count), "  max: ",max(coo.count)))
        for (i in 1:number.of.screens) {
          j <- i*4 -3
          coo_plot(x[j], main = paste0(j,": ",names(x)[j]))
          j <- i*4 -2
          if (j <= length(x)) { coo_plot(x[j], main = paste0(j,": ",names(x)[j])) }
          j <- i*4 -1
          if (j <= length(x)) { coo_plot(x[j], main = paste0(j,": ",names(x)[j])) }
          j <- i*4 -0
          if (j <= length(x)) { coo_plot(x[j], main = paste0(j,": ",names(x)[j])) }
          invisible(readline(prompt=paste("Screen",i,"of",number.of.screens," - Press [ENTER] to continue ([ESC] to quit)")))
        }
      } else {
        warning("Unrecognised data structure.")
      }
    }
    par(mfrow = c(1, 1))
  } # End of function

  images.to.Out <- function(x, sampling.depth = NULL) {
    if (is(x)[1] != "list") {
      stop("Unrecognised file format")
    }
    if (is.null(sampling.depth)) {
      sampling.depth <- min(unique(sort(unlist(lapply(x, dim))))[-1])
      message(paste("Sampling depth set to:", sampling.depth))
    }
    if (sampling.depth > min(unique(sort(unlist(lapply(x, dim))))[-1])) {
      sampling.depth <- min(unique(sort(unlist(lapply(x, dim))))[-1])
      warning(paste("Sampling depth was set too high! Using the minimum coordinate number:", sampling.depth))
    }
    coo.list <- list()
    for (i in 1:length(x)) {
      sampled.points <- sample(1:(dim(x[[i]])[1]), sampling.depth, replace = FALSE)
      # This step finds the contour points from a binary image
      coo.list[[i]] <- coo_extract(x[[i]], ids = sort(sampled.points))
    }
    names(coo.list) <- names(x)
    return(Out(coo.list))
  } # End of function

  check.harmonic.power <- function (x, max.harmonic = 15, tolerance = 0.0001, show.plot = TRUE) {
    # Error checking: `x` should be a 2-column matrix
    if (is.null(x) | ((dim(x)[2]) != 2) | (!is.matrix(x))) {
      stop("Argument `x` should be a 2-column matrix or Outline data type.")
    }
    if(is.null(tolerance)) {
      tolerance <- 1/(max.harmonic)
    }

    # Calculate cumulative harmonic power
    cumulative.harmonic.power <- x %>%
      efourier(nb.h = max.harmonic) %>%
      harm_pow() %>%
      cumsum()

    # Find the median value
    median.harmonic.power <- median(cumulative.harmonic.power)[1]

    # Find the relative deviation from the median
    deviance <- (median.harmonic.power - abs(median.harmonic.power - cumulative.harmonic.power)) / median.harmonic.power

    # Check for convergence
    if ((cumulative.harmonic.power[max.harmonic]/median.harmonic.power) - (cumulative.harmonic.power[max.harmonic-round(max.harmonic/3)]/median.harmonic.power) > tolerance) {
      if (show.plot) {
        plot(
          cumulative.harmonic.power,
          type='o',
          xlab='Harmonic rank',
          ylab='Cumulative harmonic power'
        )
        abline(h=median.harmonic.power, col = "darkred")
      }
      stop(paste0("Process may not have reached stability at `max.harmonic = ",max.harmonic,"`. Consider increasing this value."))
    }

    # Flag the first harmonic within tolerance of the median as optimal
    optimal.harmonic <- which(deviance >=  (1 - tolerance))[1]

    # Plots
    if (show.plot) {
      # Two plots panels
      par(mfrow = c(1,2))

      plot(
        cumulative.harmonic.power,
        type='o',
        xlab='Harmonic rank',
        ylab='Cumulative harmonic power'
      )
      abline(h=median.harmonic.power, col = "darkred")
      abline(v=optimal.harmonic, col = "darkred")

      x %T>%
        coo_plot() %>%
        # coo_slide(ldk = 2) %>%
        efourier(
          nb.h = optimal.harmonic,
          norm = FALSE
        ) %>%
        efourier_i() %T>%
        coo_draw(
          coo = .,
          border="darkred",
          lwd = 2
        )
      par(mfrow = c(1,1))
    }
    return(optimal.harmonic)
  } # End of function

}

# Get the file names
(image.files <- list.files(image.folder, full.names=TRUE))

# Check one of the images, just to confirm
?img_plot


# Import all outlines based on these images
# Note that the `threshold` argument can be adjusted based on the images' contrast
wing.outlines <- import_jpg(image.files, threshold = 0.5)

# Review the images, checking that each one is (approximately) the expected shape.
review.outlines(   )

# Convert to the outline (`Out`) format
# A very large number of coordinate points can cause problems for later steps.
# The `sampling.depth` argument can be used to reduced the number of points.
wing.outlines <- images.to.Out(    )

# Add metadata
# For this exercise metadata is embedded in file names

?names

?str_split_fixed



wing.outlines$fac <- tibble::tibble(

)

names(wing.outlines)




# Review the images again for quality control




# Also possible using tools from `Momocs`






# Check the balance of sampling






############################################################
# Generalized Procrustes Analysis (GPA)
############################################################

# Note that the function does not calculate centroid size
# Instead this workflow saves `Csize` as a separate object
?fgProcrustes




# Test harmonic power

check.harmonic.power(   )




# EFA
efa.wing.outlines <- gpa.wing.outlines %>%
  efourier(
    nb.h =   ,
    norm = FALSE
  )



# The matrix of coefficients




# Plot the values of the harmonic coefficients
efa.wing.outlines %>%
  rm_harm() %>%
  boxplot()




# Harmonic contribution to shape
?hcontrib





############################################################
# Principal Component Analysis (PCA)
############################################################


?PCA


?plot_PCA


plot_PCA(
  x = ,
  f = ~ ,
  chullfilled = ,
  labelgroups = ,
  morphospace_position = ,
  center_origin =
)




# Separate FW and HW

?which











############################################################
# Modeling (using multivariate ANOVA a.k.a MANOVA)
############################################################

# There are many possible ways to implement MANOVA in R.
# While the `Momocs` package actually has a function `MANOVA`, I don't recommend it.
# It is difficult to use and there's no reliable implementation of a pairwise
# post hoc test.
# Instead, the code here uses multivariate linear modeling with randomized
# residual permutation, as implemented in the function `lm.rrpp` from the
# `RRPP` package.

# The function expects a list that contains a matrix of the
# multivariate data and other elements that contain the predictive factors.
list.wing.outlines <- list(
  coe = efa.wing.outlines$coe, # the harmonic coefficients
  wing = efa.wing.outlines$fac$wing,
  treatment = efa.wing.outlines$fac$treatment
)

# Build the model
model.coe.by.wing <- lm.rrpp(
  ~ ,
  data = list.wing.outlines,
  iter = (1e4)-1
)

anova(model.coe.by.wing)


# Post hoc pairwise tests (for model factors with more than 2 levels)

pw.coe.by.wing <- pairwise(
  fit = model.coe.by.wing,
  # fit.null = ?,  # If you don't specify a null model, the function will guess using `reveal.model.designs`
  groups = list.wing.outlines$wing
)

summary(pw.coe.by.wing)


# More complex models



# Just forewings



# Just hindwings





############################################################
# Linear Discriminant Analysis (LDA)
############################################################

?LDA

lda.wing.outlines <- LDA(
  x = efa.wing.outlines,
  fac = ~ treatment
)

plot_LDA(
  x = lda.wing.outlines,
  chullfilled = TRUE,
  labelgroups = TRUE,
  morphospace_position = "range",
  center_origin = FALSE
)

plot_CV(
  x = lda.wing.outlines
)

plot_CV2(
  x = lda.wing.outlines
)



# The End
