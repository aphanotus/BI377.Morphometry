# BI377 Fall 2025
# In-Class Code
# Week 13

# REFERENCE

############################################################
# Outline-based Elliptical Fourier Analysis (EFA) of shape
############################################################

# Load the library `Momocs`
library(Momocs)

# setwd("~/Documents/5.teaching/openEd/BI377.Morphometry/class4.September29/")

# Defining some additional useful functions (soon to be added to `borealis`!)
images.to.Out <- function(x, sampling.depth = NULL) {
  if (is.null(sampling.depth)) {
    sampling.depth <- min(unique(sort(unlist(lapply(fish, dim))))[-1])
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

# Get the file names
# This command assumes high-contrast image files are in a folder named `wing.outlines`
(image.files <- list.files("wing.outlines", full.names=TRUE))

# Check one of the images, just to confirm
img_plot(image.files[1])

# Import all outlines based on these images
# Note that the `threshold` argument can be adjusted based on the images' contrast
wing.outlines <- import_jpg(image.files, threshold = 0.5)

# Convert to the outline (`Out`) format
# The `sampling.depth` argument can be used to reduced the number of points in large images 
wing.outlines <- images.to.Out(wing.outlines, sampling.depth = 2000)

# Add metadata
wing.outlines$fac <- tibble::tibble(
  species = factor(c("fervidus","borealis","borealis","borealis","borealis","borealis","borealis","borealis","fervidus","fervidus")),
  Csize = coo_centsize(wing.outlines)
)


# Check out the outlines
coo_plot(wing.outlines[1])

panel(wing.outlines, fac = "species", names=TRUE)

stack(wing.outlines, fac = "species")


# Generalized Procrustes Analysis
# Note that the function does not calculate centroid size
# Instead this workflow saves `Csize` as a seperate object
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


# Principal Component Analysis (PCA)
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

# Linear Discriminant Analysis (LDA)
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

# Modeling (using multivariate ANOVA a.k.a MANOVA)

# The argument `retain` indicates the number of harmonics to include in the analysis.
# This code uses all of them and draws that number from the EFA object directly.

# The `test` argument determines the metric for effect size. The options, in order of
# robustness, are "Roy" (Roy's largest root), "Hotelling-Lawley", (Hotelling's trace),
# "Wilks" (Wilks' Lambda -- the recommendation) and "Pillai" (Pillai's trace).

MANOVA(
  x = efa.wing.outlines,
  fac = ~ species,
  test = "Wilks", 
  retain = (ncol(efa.wing.outlines$coe) / 4 )
)

# There is no reliable implementation of a pairwise post hoc test for Momocs::MANOVA
# However, you could perform pairwise MANOVAs, subsetting the dataset to pairs each time.
# Then apply a multiple-test correction to interpret the p-values.


# Disparity analysis
library(geomorph)

gdf.wing.outlines <- geomorph.data.frame(
  coords = efa.wing.outlines$coe, 
  species = efa.wing.outlines$species
)

# To group the "shape space" variance of groups
morphol.disparity(
  coords ~ 1, 
  groups = ~ species, 
  data = gdf.wing.outlines, 
  iter = 999
)



