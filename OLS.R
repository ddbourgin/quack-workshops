#------------------------------------------------------------
#       Deriving the OLS Estimator Using Linear Algebra
#                 QuACK Workshop, 10.13.16
#------------------------------------------------------------

# install and load various packages
installed <- installed.packages()

if (!"scatterplot3d" %in% installed) {install.packages("scatterplot3d")}
if (!"MASS" %in% installed) {install.packages("MASS")}
if (!"car" %in% installed) {install.packages("car")}

library("scatterplot3d")
library("MASS")
library("car")

# seed the rng for reproducibility
set.seed(9608)

k <- 2   # number of IVs
n <- 150 # number of datapoints

# sample mean and variances from a uniform
# distribution over [0, 1]
means <- runif(k, min=0, max=1)

# create a k-by-k symmetric covariance matrix
covariance <- matrix(rnorm(k^2, 0, 0.5), c(k, k))
covariance[lower.tri(covariance)] <- t(covariance)[lower.tri(covariance)]
covariance <- t(covariance) %*% covariance

# generate a fake dataset by sampling n observations from a
# multivariate normal distribution
df <- mvrnorm(n, means, covariance)

# make df into a design matrix by adding a column of 1's
df <- cbind(rep.int(1, n), df)

# generate a vector of normally distributed noise from a
# normal distribution with mean 0 and variance 1
error <- rnorm(n, 0, 1)

# sample the "true" coefficient values for our
# regression model from a uniform distribution
# on [-1, 1]
beta_true <- runif(k + 1, min=-1, max=1)

# generate noisy observations of our DV
y <- df %*% beta_true + error

# convert df into a dataframe for plotting
df <- cbind(df, y)
df <- as.data.frame(df)
colnames(df) <- c("Intercept", "X1", "X2", "y")

# get some descriptive stats
summary(df)

# scatterplot matrix for pairwise comparisons
scatterplotMatrix(~ y + X1 + X2,
                  data = df,
                  ellipse = TRUE,
                  diagonal = "density",
                  smoother = FALSE,
                  reg.line = FALSE,
                  main = "df scatterplot matrix")


# Create the appropriate matrices from our dataframe
X <- cbind(df$Intercept, df$X1, df$X2)
y <- as.matrix(df$y)
X.T <- t(X)

# solve for beta_hat
# beta_hat = (X'X)^(-1) X'y
beta_hat <- solve(X.T %*% X) %*% X.T %*% y

layout(matrix(c(1, 2, 1, 3), 2, 2, byrow = TRUE),
       widths = c(1.5, 1, 1.5, 1),
       heights = c(2, 2, 1, 1))

# plot the regression fit associated with beta_hat
plt.fit <- scatterplot3d(x = df$X1,
                         y = df$y,
                         z = df$X2,
                         xlab = "X1",
                         ylab = "X2",
                         zlab = "y",
                         pch = 19,
                         scale.y = 0.5,
                         color = 'blue',
                         main = "Regression fit (hand)")
plt.fit$plane3d(beta_hat, draw_polygon = TRUE)

# create a scatter plot for X1 vs. y with the regression
# fit in red and true function in green
plot(df$X1, df$y, xlab = 'X1', ylab = 'y', main = 'X1 vs. y')
abline(beta_hat[1], beta_hat[2], col = 'red')
abline(beta_true[1], beta_true[2], col = 'green')

# create a scatter plot for X2 vs. y with the regression
# fit in red and true function in green
plot(df$X2, df$y, xlab = 'X2', ylab = 'y', main = 'X2 vs. y')
abline(beta_hat[1], beta_hat[3], col = 'red')
abline(beta_true[1], beta_true[3], col = 'green')

# generate regression fits using the lm command to verify
# we get the same answer
my.lm <- lm(y ~ X1 + X2, data=df)

layout(matrix(c(1,2,1,3), 2, 2, byrow = TRUE),
       widths=c(1.5,1,1.5,1),
       heights=c(2, 2, 1, 1))

plt.fit2 <- scatterplot3d(x = df$X1,
                          y = df$y,
                          z = df$X2,
                          xlab = "X1",
                          ylab = "X2",
                          zlab = "y",
                          pch = 19,
                          scale.y = 0.5,
                          color = 'blue',
                          main = "Regression fit (lm)")
plt.fit2$plane3d(my.lm, draw_polygon = TRUE)

# extract coefficients from lm call
lmcoefs <- my.lm$coefficients

# create a scatter plot for X1 vs. y with the regression
# fit in red and true function in green
plot(df$X1, df$y, xlab = 'X1', ylab = 'y', main='X1 vs. y')
abline(lmcoefs[1], lmcoefs[2], col = 'red')
abline(beta_true[1], beta_true[2], col = 'green')

# create a scatter plot for X2 vs. y with the regression
# fit in red and true function in green
plot(df$X2, df$y, xlab = 'X2', ylab = 'y', main='X2 vs. y')
abline(lmcoefs[1], lmcoefs[3], col = 'red')
abline(beta_true[1], beta_true[3], col='green')

# verify that we arrived at the same values
cat('Coefficients found using R\'s lm command: \n', lmcoefs)
cat('Coefficients found "by hand:" \n', beta_hat)
cat('True coefficients:\n', beta_true)