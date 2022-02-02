# Define auxiliary functions ---------------------------------------------------
get_X <- function(n, v = 5, corr = c(-0.2, 0, 0), mu = c(0, 0, 0), .seed = 123) {
  S <- diag(v, length(corr), length(corr))
  S[upper.tri(S)] <- corr * sqrt(v)^2  # set correlations
  S[lower.tri(S)] <- t(S)[lower.tri(S)]
  set.seed(.seed)
  x <- MASS::mvrnorm(
    n_pop,
    mu = mu,
    Sigma = S
  )
  X <- as.matrix(cbind(1, x))
  
  return(X)
}

odds <- function(p) p/(1-p)
inv_odds <- function(o) o/(1 + o)
logit <- function(p) log(odds(p))
inv_logit <- function(p) plogis(p)
add_noise <- function(p, sd = 0.5) inv_logit(
  logit(p) + rnorm(length(p), mean = 0, sd = sd)
)
get_auc <- function(y, x, .digits = 3) round(
  pROC::auc(y, x), .digits
)

get_auc_from_fit <- function(obj) {
  0.5 + validate(obj, B = 200)[1,5]/2
}

get_prior_modeling_offset <- function(p_hat, sampling_ratio) {
  log(sampling_ratio *(1-p_hat)/p_hat)
}

sample_case_control <- function(df, n = 200,
                                noise = NULL, cutoffs = NULL) {

  "
  instead of outcome risk, use predictor value for sampling probability?
  "
  
  if (!is.null(cutoffs)) {
    stopifnot(length(cutoffs) == 2)
    high_risk_cases_ix <- which(df$y == 1 & df$p > cutoffs[1])
    low_risk_controls_ix <- which(df$y == 0 & df$p < cutoffs[2])
    df <- df[union(high_risk_cases_ix, low_risk_controls_ix), ]
    p_cases <- p_controls <- NULL
    
  } else {
    
    p_cases <- df$p[df$y == 1]
    p_controls <- 1 - df$p[df$y == 0]
    
    if (!is.null(noise)) {
      p_cases <- add_noise(p = p_cases, sd = noise)
      p_controls <- add_noise(p = p_controls, sd = noise)
    }
    
  }
  
  cases_ix <- sample(
    which(df$y == 1),
    round(n/2),
    prob = p_cases
  )
  controls_ix <- sample(
    which(df$y == 0),
    round(n/2),
    prob = p_controls
  )

  stopifnot(  # make sure samples are not mixed up
    length(intersect(cases_ix, controls_ix)) == 0
  )
  
  
  return(df[c(cases_ix, controls_ix),])
}

sample_cross_sectional <- function(df, n) {
  cs_sample_ix <- sample(1:nrow(df), n)
  return(df[cs_sample_ix, ])
}

group_differences <- \(d) {
  d %>% 
    group_by(y) %>% 
    summarise_if(is.numeric, mean) %>% 
    select(y, contains('x')) %>% 
    pivot_wider(names_from = y,
                values_from = contains('x')) %>% 
    mutate(x1_diff = x1_1 - x1_0,
           x2_diff = x2_1 - x2_0,
           x3_diff = x3_1 - x3_0) %>% 
    select(contains('diff'))
}

plot_disease_prob_stacked <- function(df) {
  df %>% 
    ggplot(aes(x = p, fill = factor(y))) +
    geom_histogram(bins = 30) +
    labs(
      x = "Disease probability",
      fill = "Disease label"
    ) +
    scale_fill_manual(values = c('0' = "skyblue", '1' = "magenta")) +
    scale_x_continuous(breaks = scales::pretty_breaks(5)) 
}

plot_disease_prob_nonstacked <- function(df) {
  ggarrange(
    df %>% 
      filter(y == 0) %>% 
      ggplot(aes(x = p, fill = factor(y))) +
      geom_histogram(fill = "skyblue", bins = 30) +
      labs(x = NULL) +
      scale_x_continuous(breaks = scales::pretty_breaks(5)) +
      coord_cartesian(xlim = c(0, 1)) +
      theme_minimal(),
    df %>% 
      filter(y == 1) %>% 
      ggplot(aes(x = p, fill = factor(y))) +
      geom_histogram(fill = "magenta", bins = 30) +
      scale_x_continuous(breaks = scales::pretty_breaks(5)) +
      coord_cartesian(xlim = c(0, 1)) +
      theme_minimal() +
      labs(x = "Disease probability"),
    nrow = 2
  )
}

plot_disease_prob <- function(df, title = NULL) {
  p <- ggarrange(
    plot_disease_prob_stacked(df),
    plot_disease_prob_nonstacked(df),
    ncol = 2
  )
  if (!is.null(title)) {
    p <- annotate_figure(p, top = title)
  }
  p
}

format2 <- function(x, .digits = 3) {
  round(unname(x), .digits)
}
