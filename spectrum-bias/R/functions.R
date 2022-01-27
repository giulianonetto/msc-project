# Define auxiliary functions ---------------------------------------------------
sample_case_control <- function(df, n_sample = 200) {
  logit <- \(p) log(p/(1-p))
  sampling_prob <- \(p, sd=0.5) plogis(logit(p) + rnorm(length(p), sd=sd))
  cases_ix <- sample(
    which(df$y == 1),
    round(n_sample/2),
    prob = sampling_prob(df$p[df$y == 1], sd = .5)
  )
  controls_ix <- sample(
    which(df$y == 0),
    round(n_sample/2),
    prob = sampling_prob(1-df$p[df$y == 0], sd = .5)
  )
  
  # cases_ix <- sample(
  #   which(df_pop$y == 1 & df_pop$p > 0.2),
  #   round(n_sample/2)
  # )
  # controls_ix <- sample(
  #   which(df_pop$y == 0 & df_pop$p < 0.6),
  #   round(n_sample/2)
  # )
  
  stopifnot(  # make sure samples are not mixed up
    length(intersect(cases_ix, controls_ix)) == 0
  )
  
  
  return(df[c(cases_ix, controls_ix),])
}

diff_groups <- \(d) {
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
