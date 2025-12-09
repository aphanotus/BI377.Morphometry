# BI377 Fall 2025
# In-Class Code
# Week 13

# REFERENCE

############################################################
############################################################
# Outline-based Elliptical Fourier Analysis (EFA) of shape
############################################################
############################################################

# Load the library `Momocs`
library(Momocs)

# Set your working directory to a place one level down from 
# the folder containing the images for analysis 
# setwd("~/Documents/5.teaching/openEd/BI377.Morphometry/class4.September29/")

# Defining some useful functions (To be added to `borealis` after this semester)
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
# This command assumes high-contrast image files are in a folder named `wing.outlines`
(image.files <- list.files("wing.outlines", full.names=TRUE))

# Check one of the images, just to confirm
img_plot(image.files[1])


# Note regarding image quality
# The `Momocs` function for importing images `import_jpg` attempts to detect an 
# object’s edges using a random seed. So it may work, of it may outline a small
# feature that is not the intended object. It may works sometimes and not others. 
# A solution is to select the object of interest in each image using the "magic" 
# selection tool in Photoshop (which detects edges more reliably), then delete it 
# and replace it with black (as the foreground color), then invert the selection
# and delete it in favor of the background (white). This reliably results in a 
# successful import. However, it can be tedious for a large number of image files.


# Import all outlines based on these images
# Note that the `threshold` argument can be adjusted based on the images' contrast
wing.outlines <- import_jpg(image.files, threshold = 0.5)

# Review the images, checking that each one is (approximately) the expected shape.
review.outlines(wing.outlines)

# Convert to the outline (`Out`) format
# A very large number of coordinate points can cause problems for later steps.
# The `sampling.depth` argument can be used to reduced the number of points.
wing.outlines <- images.to.Out(wing.outlines, sampling.depth = 2000)

# Add metadata
wing.outlines$fac <- tibble::tibble(
  species = factor(c("fervidus","borealis","borealis","borealis","borealis","borealis","borealis","borealis","fervidus","fervidus")),
  Csize = coo_centsize(wing.outlines)
)


# Review the images again for quality control
review.outlines(wing.outlines)

# Also possible using tools from `Momocs`
panel(wing.outlines, fac = "species", names=TRUE)

stack(wing.outlines, fac = "species")

############################################################
# Generalized Procrustes Analysis (GPA)
############################################################

# Note that the function does not calculate centroid size
# Instead this workflow saves `Csize` as a separate object
gpa.wing.outlines <- fgProcrustes(wing.outlines)

stack(gpa.wing.outlines, fac = "species")


# Test harmonic power
# I suggest examining harmonic power for a few representative outlines,
# but not necessarily every one in a large dataset. The optimal number 
# may differ if you look at different individuals. However, for downstream
# analyses that use harmonic coefficients, it will be important that all
# specimens in the dataset undergo EFA with the same number of harmonics.
check.harmonic.power(gpa.wing.outlines[1])

# EFA with 7 harmonics
efa.wing.outlines <- gpa.wing.outlines %>%
  # coo_slide(ldk = 2) %>%
  efourier(
    nb.h = 7,
    norm = FALSE
  )
# The `norm` argument can results in bad EFA fits for outlines
# that are circular or strongly bilaterally symmetrical.

############################################################
# Principal Component Analysis (PCA)
############################################################

pca.wing.outlines <- PCA(efa.wing.outlines)
  
?plot_PCA

plot_PCA(
  x = pca.wing.outlines, 
  f = ~ species, 
  chullfilled = TRUE,
  labelgroups = TRUE,
  morphospace_position = "range_axes", # options: "range", "full", "circle", "xy", "range_axes", "full_axes"
  center_origin = FALSE
)

############################################################
# Linear Discriminant Analysis (LDA)
############################################################

lda.wing.outlines <- LDA(
  x = efa.wing.outlines,
  fac = ~ species
)

plot_LDA(
  x = lda.wing.outlines, 
  chullfilled = TRUE,
  labelgroups = TRUE,
  morphospace_position = "range_axes", # options: "range", "full", "circle", "xy", "range_axes", "full_axes"
  center_origin = FALSE
)

############################################################
# Modeling (using multivariate ANOVA a.k.a MANOVA)
############################################################

# There are many possible ways implement MANOVA in R.
# While the `Momocs` package actually has a function `MANOVA`, I don't recommend it.
# It is difficult to use and there's no reliable implementation of a pairwise 
# post hoc test.
# Instead, the code here uses multivariate linear modeling with randomized 
# residual permutation, as implemented in the function `lm.rrpp` from the 
# `RRPP` package, which is one of the packages under the hood of `geomorph`)

# The function expects a list that contains a matrix of the 
# multivariate data and other elements that contain the predictive factors.
list.wing.outlines <- list(
  coe = efa.wing.outlines$coe, # the harmonic coefficients 
  species = efa.wing.outlines$species,
  Csize = wing.outlines$fac$Csize
)

# Build the model just like we did using `procD.lm`
model.wings.outlines.by.size <- lm.rrpp(
  coe ~ Csize, 
  data = list.wing.outlines, 
  iter = (1e4)-1
)
anova(model.wings.outlines.by.size)

model.wings.outlines.by.species <- lm.rrpp(
  coe ~ species, 
  data = list.wing.outlines, 
  iter = (1e4)-1
)
anova(model.wings.outlines.by.species)

# The results can be interpreted just like the ANOVA tables 
# generated by `geomorph::procD.lm`

# Post hoc pairwise tests are implements the same way as in our 
# landmark-based workflow

pw.wings.outlines.by.species <- pairwise(
  fit = model.wings.outlines.by.species,
  # fit.null = ?,  # If you don't specify a null model, the function will guess using `reveal.model.designs`
  groups = list.wing.outlines$species
)

summary(pw.wings.outlines.by.species)


############################################################
# Disparity analysis
############################################################

library(geomorph)

# The function `morphol.disparity` also expects a list that contains a matrix of the 
# multivariate data and other elements that contain predictive factors.

morphol.disparity(
  coe ~ 1, 
  groups = ~ species, 
  data = list.wing.outlines, 
  iter = 999
)



