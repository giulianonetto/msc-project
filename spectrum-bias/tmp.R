v <- 5
S <- diag(2, 3, 3)
S[upper.tri(S)] <- c(-0.2, 0, 0) * sqrt(v)^2  # set correlations
S[lower.tri(S)] <- t(S)[lower.tri(S)]
X <- cbind(
  1,
  MASS::mvrnorm(
    n_pop,
    mu = c(0, 0, 0),
    Sigma = S
  )
) %>% as.matrix()

{
  pp <-plogis(X %*% c(-1.5, log(1.2), log(.3), log(3)))
  print(mean(pp))
  hist(pp)
}
