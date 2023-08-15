library(BITS)

# Generate toy data
set.seed(123)
maf <- 0.25
n <- 3000
p <- 50
X <- matrix(sample(0:2, p * n, replace = TRUE,
                   prob = c((1-maf)^2, 1-(1-maf)^2-maf^2, maf^2)), ncol = p)
truth <- X[,1] + X[,2] * (X[,3] > 0) + 2 * X[,2] * (X[,4] < 2) * (2-X[,5])
y <- truth + rnorm(n, 0, sd(truth))

train.ind <- 1:1000; val.ind <- 1001:2000; test.ind <- 2001:3000

# Identify and apply modes of inheritance
X.moi <- applyMOI(X, MOI(X[train.ind,], y[train.ind]))

# Fit BITS models for a grid of gamma and lambda values
model <- BITS.complete(X[train.ind,], y[train.ind], max.iter = function(p) 1000,
                       gsteps = 10)
# Get ideal hyperparameters using validation data
ideal.mod <- get.ideal.model(model, X[val.ind,], y[val.ind])
# Plot validation data performance
plot(ideal.mod$val.res)

# Fit the final BITS model using the optimal hyperparameters
model <- BITS(X[c(train.ind, val.ind),], y[c(train.ind, val.ind)],
              gamma = ideal.mod$best.g, max.iter = function(p) 1000)
model$s <- ideal.mod$best.s

# Normalized mean squared error
calcNMSE(predict(model, X[test.ind,]), y[test.ind])
# Terms included in the final model
get.included.vars(model)
