# Load packages ----------------------------------------------------------------
library(tidyverse)
library(ggpubr)


# Source auxiliary R files -----------------------------------------------------
source("R/functions.R")


# Analysis ---------------------------------------------------------------------

## Define random seed ----

SEED <- 237429

## Define population size ----

n_pop <- 1e6

## Simulate Predictor (eg biomarker or age) ----
set.seed(SEED)
X <- cbind(
  1, 
  x1 = rnorm(n = n_pop, mean = 0, sd = 5),
  x2 = rnorm(n = n_pop, mean = 0, sd = 5),
  x3 = rnorm(n = n_pop, mean = 0, sd = 5)
)

## Define true coefficients ----

# intercept and betas such that average prob. is ~ 25%
beta <- c(-1.5, log(1.2), log(.8), log(1))

## Simulate disease probabilities ----
p <- plogis(X %*% beta)

## Simulate disease labels ----
y <- rbinom(n = n_pop, size = 1, prob = p)

## Define population data set ----
df_pop <- data.frame(
  patient_id = paste0("patient-", 1:n_pop),
  y = y,
  p = p, 
  x1 = X[, 2],
  x2 = X[, 3],
  x3 = X[, 4]
)
rm(X, p, y)

## Plot disease probabilities distribution ----
theme_set(theme_bw())
set.seed(SEED)
ix <- sample(1:n_pop, 1e4)  # subset 1000 random samples for plotting
p <- df_pop[ix, ] %>% 
  ggplot(aes(x = p, fill = factor(y))) +
  geom_histogram() +
  labs(
    x = "Disease probability",
    fill = "Disease label"
  )
ggsave("output/population_disease_probability.png",
       p, width = 8, height = 4.5)
## Simulate case-control experiment ----


df_sample <- sample_case_control(df_pop, n_sample = 400)
diff_groups(df_pop)
diff_groups(df_sample)
### Plot disease probabilities from case-control sample ----

p <- df_sample %>% 
  ggplot(aes(x = p, fill = factor(y))) +
  geom_histogram() +
  labs(
    x = "Disease probability",
    fill = "Disease label",
    title = "Case-control sample"
  )
ggsave("output/case-control-sample_disease_probability.png",
       p, width = 8, height = 4.5)

### Fit logistic regression model ----

#### Calculate offset for debiasing intercept ----
prev_hat <- mean(df_pop$p[sample(1:n_pop, 500)])  # mimic external info
df_sample$.offset <- log(0.5/0.5 *(1-prev_hat)/prev_hat)
fit <- glm(
  formula = y ~ x1 + x2 + x3 + offset(.offset), 
  data = df_sample, 
  family = binomial
)

#### validate in cross-sectional data ----

df_val <- df_pop %>% 
  filter(!(patient_id %in% df_sample$patient_id)) %>% 
  sample_n(1e3) %>% 
  mutate(.offset = df_sample$.offset[1])

df_val$lp <- predict(fit, newdata = df_val) - df_sample$.offset[1]
df_val$p_hat <- plogis(df_val$lp)

plot(df_val$p_hat, df_val$p); abline(0,1,col='red')
glm(y ~ lp, family = binomial, data = df_val) %>% coef()

dca <- rmda::decision_curve(
  y ~ p_hat, 
  data = df_val %>% mutate(y = as.numeric(y == 1)), 
  fitted.risk = T
)
rmda::plot_decision_curve(dca, standardize = F)
pROC::auc(df_val[, "y"], df_val[, "p_hat"])

### validate in case-control data ----

df_val2 <- sample_case_control(
  df_pop %>% 
    filter(!(patient_id %in% df_sample$patient_id)),
  n_sample = 1e3
) %>% 
  mutate(.offset = df_sample$.offset[1])

df_val2$lp <- predict(fit, newdata = df_val2) - df_sample$.offset[1]
df_val2$p_hat <- plogis(df_val2$lp)

plot(df_val2$p_hat, df_val2$p); abline(0,1,col='red')
glm(y ~ lp, family = binomial, data = df_val2) %>% coef()

dca <- rmda::decision_curve(
  y ~ p_hat, 
  data = df_val2 %>% mutate(y = as.numeric(y == 1)), 
  fitted.risk = T
)
rmda::plot_decision_curve(dca, standardize = F)
pROC::auc(df_val2[, "y"], df_val2[, "p_hat"])

"
basically, we want to show two problems:
(1) problem with training with case-control data
(2) problem with validating with case-control data

(1) may look at bad performance: disc, cal, nb
(2) may look at bias in performance estimates: disc, cal, nb (external at first, but also bias in CV estimates)
"

