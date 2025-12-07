# BI377 Fall 2025
# In-Class Code
# Week 13

# REFERENCE

############################################################
############################################################
# Analysis of multivariate linear data as "shape"
############################################################
############################################################

# First some general adivce:
# As we discussed early in the semester, multiple linear measurements can be used
# as an approximation of shape. In some cases this way be the only method possible,
# for example with soft and flexible subject, where landmarks or outlines may not be
# consistent, but where measurements along or across a structure might be consistent.
# However, there are other situations in which multiple linear measurements might 
# provide predictive (independent) variables. In most examples this semester we've 
# dealt with categorical predictive factors (like species or sex), but we have also
# used continuous predictive variables like centroid size. It's also possible to use
# multiple continuous variables as predictors. Typically this is done using PC1 
# (and maybe PC2) following PCA on those variables, then modeling the dependent 
# "shape" variables.

library(tidyverse)
library(borealis)
library(viridis)

# An example dataset: `mtcars`
# Data from the 1974 Motor Trend US magazine, and comprises fuel consumption 
# and 10 aspects of automobile design and performance for 32 automobiles 
# (1973--74 models).
head(mtcars)
dim(mtcars)

# You don't need to know anything about cars to follow this example,
# but the meaning of each column is given below.
# To my point above, some of these variables can be thought of as 
# "independent" or intrinsic aspects of each car's "anatomy" (or shape).
# While other variables might be considered "dependent" or 
# fitness/performace measurements.

# column  description                               
# mpg	    Miles/(US) gallon                         dependent
# cyl	    Number of cylinders                       independent
# disp	  Displacement (cu.in.)                     independent
# hp	    Gross horsepower                          dependent
# drat	  Rear axle ratio                           independent
# wt	    Weight (1000 lbs)                         independent
# qsec	  1/4 mile time                             dependent
# vs	    Engine (0 = V-shaped, 1 = straight)       independent
# am	    Transmission (0 = automatic, 1 = manual)  independent
# gear	  Number of forward gears                   independent
# carb    Number of carburetors                     independent

# I'll start by simple moving the car's name into a column, rather 
# the the rownames where they are now. 
cars <- mtcars %>% 
  rownames_to_column(var = "model.name") 

# It will be useful to add a factor for somer analysis. So I'll add
# a column with the nationality of the car makers. 
cars$nationality <- c(
  "Japan", "Japan", "Japan", "United States", "United States", "United States", "United States",
  "Germany", "Germany", "Germany", "Germany", "Germany", "Germany", "Germany",
  "United States", "United States", "United States", "Italy", "Japan", "Japan", "Japan",
  "United States", "United States", "United States", "United States", "Italy", "Germany", "Northern Europe",
  "United States", "Italy", "Italy", "Northern Europe"
)
# Lotus was based in the UK, and Volvo is famously Swedish, but I'm 
# combining them here.

############################################################
# Principal Component Analysis (PCA)
############################################################

# Be sure to filter out specimens with NA values (alternatively exclude 
# columns with NA values) 
# The function `dplyr` function `drop_na` is very useful for this!

# Indicate the data columns
independent.vars <- which(
  colnames(cars) %in%
    c("cyl","disp","drat","wt","vs","am","gear","crab")
)

# It's often the case that different measurements have very different magnitudes.
# If so, it's necessary to add the argument `scale = TRUE` so that variables with
# the largest magnitude don't completely dominant the results.
pca.cars <- prcomp(cars[,independent.vars], scale = TRUE)

pcvar(pca.cars)
scree.plot(pca.cars)

shape.space(
  x= pca.cars,
  group = cars$nationality,
  convex.hulls = TRUE
)

# You can examine the association (or "loading") of each variable to each PC axis.
# Larger values indicate a stronger association of that trait with the PC. 
pca.cars$rotation

biplot(pca.cars, scale = 0)

# We can also examine variation in the independent variables using PCA

dependent.vars <- which(
  colnames(cars) %in% c("mpg","hp","qsec")
)

pca.car.performance <- prcomp(cars[,dependent.vars], scale = TRUE)

pcvar(pca.car.performance)

biplot(pca.car.performance, scale = 0)

shape.space(
  x= pca.car.performance,
  group = cars$nationality,
  convex.hulls = TRUE
)

# We might choose to store PC1 as a combined metric of all 
# these performance measurements
cars$performance.PC1 <- pca.car.performance$x[,1]


############################################################
# Linear Discriminant Analysis
############################################################

# Load the MASS package, which contains the lda() function
library(MASS)

# Fit an LDA model. 
# The formula `nationality ~ .` indicates that 'nationality' is the grouping 
# variable and all other variables are used as predictors.
# `prior` sets the prior probabilities for each group (equal in this case).
# `CV = TRUE` performs leave-one-out cross-validation for classification.

number.of.levels <- length(unique(cars$nationality))

lda.cars <- lda(
  formula = cars$nationality ~ ., 
  data = cars[,independent.vars], 
  prior = rep(1,number.of.levels)/number.of.levels,
  CV = TRUE
)

# View the posterior probabilities of group membership for each observation
# That is, for each car, what's the probability of it's nationality based
# on the LDA analysis?
lda.cars$posterior %>% 
  tibble() %>% 
  mutate(model.name = cars$model.name)

# You can also fit an LDA model without cross-validation to get the model details
lda.full.cars <- lda(
  formula = cars$nationality ~ ., 
  data = cars[,independent.vars]
)

# View the details of the full LDA model, including group means and 
# discriminant coefficients
print(lda.full.cars)

# Make predictions on new data (or a subset of the original data)
# For example, let's predict on the first 10 observations

test_data <- cars[1:10, independent.vars]

predictions <- predict(
  object = lda.full.cars,
  newdata = test_data
)

# View the predicted classes and posterior probabilities
predictions$posterior %>% 
  tibble() %>% 
  mutate(model.name = cars$model.name[1:10])

# Note that the Datsun actually gets mis-classified!
# But we are doing LDA with a relatively small dataset here.

# Add the LD scores and model names to the original test data for plotting
test_data$LD1 <- predictions$x[, 1]
test_data$LD2 <- predictions$x[, 2] # If applicable (more than 2 groups)
test_data$model.name <- cars$model.name[1:10]
test_data$nationality <- cars$nationality[1:10]

# Plot the first two linear discriminants
test_data %>% 
  ggplot(aes(x = LD1, y = LD2, color = nationality, label = model.name)) +
  theme_bw() +
  geom_point(size = 3) +
  geom_text(vjust=2, show.legend = FALSE) + 

############################################################
# Modeling (using multivariate ANOVA a.k.a MANOVA)
############################################################

# We'll perform multivariate linear modeling with randomized residual 
# permutation, as implemented in the function `lm.rrpp` from the `rrpp` 
# package, which is one of the packages supporting `geomorph`.

# The function we'll use expects a list that contains a matrix of the 
# multivariate data and other elements that contain the predictive factors.
car.list <- list(
  m = as.matrix(cars[,independent.vars]),
  nationality = cars$nationality,
  mpg = cars$mpg,
  hp = cars$hp,
  qsec = cars$qsec,
  performance.PC1 = cars$performance.PC1
)

# Build the model just like we did using `procD.lm`
model.cars.by.nationality <- lm.rrpp(
  m ~ nationality, 
  data = car.list, 
  iter = (1e4)-1
)

anova(model.cars.by.nationality)
#             Df     SS    MS     Rsq      F     Z Pr(>F)    
# nationality  4 339257 84814 0.71219 16.703 5.443  1e-04 ***
# Residuals   27 137098  5078 0.28781                        
# Total       31 476355       
 
# The results can be interpreted just like the ANOVA tables 
# generated by `geomorph::procD.lm`

# Post hoc pairwise tests are implements the same way as in our 
# landmark-based workflow

pw.cars.by.nationality <- pairwise(
  fit = model.cars.by.nationality,
  # fit.null = ?,  # If you don't specify a null model, the function will guess using `reveal.model.designs`
  groups = car.list$nationality
)

summary(pw.cars.by.nationality)
# Pairwise distances between means, plus statistics
#                                        d UCL (95%)          Z Pr > d
# Germany:Italy                  45.398975  147.0998 -0.1318661 0.5647
# Germany:Japan                  80.506184  130.4233  0.7859513 0.2295
# Germany:Northern Europe        88.290564  184.3426  0.3819662 0.3827
# Germany:United States         161.878339  109.3973  2.3688049 0.0027
# Italy:Japan                    35.127808  156.3528 -0.4561829 0.6721
# Italy:Northern Europe          42.906885  204.9362 -0.5572536 0.7039
# Italy:United States           207.264896  138.8105  2.4099939 0.0019
# Japan:Northern Europe           7.835148  191.1275 -1.6428428 0.9384
# Japan:United States           242.378562  119.9530  3.0774552 0.0001
# Northern Europe:United States 250.157852  177.1503  2.2314011 0.0034


############################################################
# Disparity analysis
############################################################

library(geomorph)

# For disparity analysis using linear measurements, think carefully about the
# scale of the measurements. Those with larger magnitude will dominate.
# If you have many measurements in the same units and their means are within an
# order of magnitude, then it may be fine to proceed with unscaled values.

morphol.disparity(
  m ~ 1, 
  groups = ~ nationality, 
  data = car.list, 
  iter = 999
)

# However, in this example `disp` is much larger than the other. 
# head(car.list$m)

# Therefore it makes more sense to based a disparity analysis on 
# scaled PC ("shape space") values.

pca.cars$nationality <- cars$nationality

morphol.disparity(
  x ~ 1, 
  groups = ~ nationality, 
  data = pca.cars,
  iter = 999
)

