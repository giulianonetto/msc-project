# Load packages ----------------------------------------------------------------
library(tidyverse)
library(ggpubr)
library(patchwork)
library(rms)
library(rmda)
library(pROC)
#library(caret) # only createFolds function used
theme_set(theme_bw(base_size = 14))
source("R/functions.R") 
# create output dir for figures and stuff
dir.create("output/manuscript", showWarnings = FALSE)

# Simulate pop ----

SEED <- 123456
n_pop <- 1e6
X <- get_X(
  n = n_pop,
  v = 2,                  # common variance
  corr = c(0, 0, 0),   # correlations x1-x2, x1-x3, x2-x3
  mu = c(0, 0, 0),        # all centered
  .seed = SEED
)

beta <- c(-1.6, log(1.2), log(.3), log(3))
p <- plogis(X %*% beta)
set.seed(SEED)
y <- rbinom(n = n_pop, size = 1, prob = p)
df_pop <- data.frame(
  patient_id = paste0("patient-", 1:n_pop),
  y = y,
  p = p, 
  x1 = X[, 2],
  x2 = X[, 3],
  x3 = X[, 4]
)
head(df_pop)


# Max AUC in pop ----

.title <- paste0(
  "max AUC ", get_auc(df_pop$y, df_pop$p)
)
plot_disease_prob(df_pop, title = .title)

round(mean(between(df_pop[df_pop$y == 0, "p"], 0.8, 1))*100,2)
round(mean(between(df_pop[df_pop$y == 1, "p"], 0, 0.3))*100,2) # 20% of cases have risk btw 0 and 30%

# Samples ----

## cross-sectional ----
n_sample <- 2e4
df_sample_cs <- sample_cross_sectional(df = df_pop, n = n_sample, .seed = SEED)
.title <- paste0(
  "best possible AUC with this sample ", get_auc(df_sample_cs$y, df_sample_cs$p)
)
plot_disease_prob(df = df_sample_cs)

#### Figure 1B ----
.p <- df_sample_cs %>% 
  mutate(
    y = ifelse(y==1, "Cases", "Controls") %>% 
      factor(levels = c("Controls", "Cases"))
  ) %>% 
  plot_disease_prob_stacked()  +
  labs(fill = NULL, y = "Frequency", x = "True disease probability") +
  theme(legend.position = c(.8, .8))
ggsave(
  "output/manuscript/disease-probabilities-distribution-cs-sample.png",
  .p, width = 6.5, height = 4, bg = "white", dpi = 600
)

round(mean(between(df_sample_cs[df_sample_cs$y == 0, "p"], 0.8, 1))*100, 2)
round(mean(between(df_sample_cs[df_sample_cs$y == 1, "p"], 0, 0.3))*100,2)
## Case-control ----

### sample ----
cc_cutoffs <- c(0.3, 0.7)
df_sample_cc <- sample_case_control(df = df_pop, n = n_sample, 
                                    cutoffs = cc_cutoffs,
                                    .seed = SEED)
.title <- paste0(
  "best possible AUC with this sample ", get_auc(df_sample_cc$y, df_sample_cc$p)
)
plot_disease_prob(df = df_sample_cc, title = .title)

#### Figure 1C ----
.p <- df_sample_cc %>% 
  mutate(
    y = ifelse(y==1, "Cases", "Controls") %>% 
      factor(levels = c("Controls", "Cases"))
  ) %>% 
  plot_disease_prob_stacked()  +
  labs(fill = NULL, y = "Frequency", x = "True disease probability") +
  theme(legend.position = c(.8, .8))
ggsave(
  "output/manuscript/disease-probabilities-distribution-cc-sample.png",
  .p, width = 6.5, height = 4, bg = "white", dpi = 600
)


# Fit models ----

## CS sample ----
fit_cs <- lrm(y ~ x1 + x2 + x3, data = df_sample_cs, x = TRUE, y = TRUE)
coef(fit_cs) %>% print_coef()


## CC sample ----
fit_cc <- lrm(y ~ x1 + x2 + x3, data = df_sample_cc, x = TRUE, y = TRUE)
coef(fit_cc) %>% print_coef()

## AUCs ----
auc_cs <- get_auc_from_fit(fit_cs)
cat("\nEstimated coeffs: ", exp(coef(fit_cs)[-1]),
    "\nEstimated AUC (internal validation): ", auc_cs)
auc_cc <- get_auc_from_fit(fit_cc)
cat("\nEstimated coeffs: ", exp(coef(fit_cc)[-1]),
    "\nEstimated AUC  (internal validation): ", auc_cc)

## Cross-validation ----

### CS model ----
folds_cs <- caret::createFolds(df_sample_cs$y, k = 10)
cv_cs_model <- map_df(folds_cs, ~{
  fold_testing <- df_sample_cs[.x, ]
  fold_training <- df_sample_cs[!(df_sample_cs$patient_id %in% fold_testing$patient_id), ]
  fold_fit <- lrm(y ~ x1 + x2 + x3, data = fold_training, x = TRUE, y = TRUE)
  fold_p_hat <- predict(fold_fit, newdata = fold_testing, type = "fitted")
  .roc <- pROC::roc(fold_testing$y, fold_p_hat)
  tibble(
    auc = .roc$auc[1],
    sens = .roc$sensitivities,
    spec = .roc$specificities
  )
}, .id = 'fold_id')

### CC model ----
folds_cc <- caret::createFolds(df_sample_cc$y, k = 10)
cv_cc_model <- map_df(folds_cc, ~{
  fold_testing <- df_sample_cc[.x, ]
  fold_training <- df_sample_cc[!(df_sample_cc$patient_id %in% fold_testing$patient_id), ]
  fold_fit <- lrm(y ~ x1 + x2 + x3, data = fold_training, x = TRUE, y = TRUE)
  fold_p_hat <- predict(fold_fit, newdata = fold_testing, type = "fitted")
  .roc <- pROC::roc(fold_testing$y, fold_p_hat)
  tibble(
    auc = .roc$auc[1],
    sens = .roc$sensitivities,
    spec = .roc$specificities
  )
}, .id = 'fold_id')

bind_rows(
  cv_cc_model %>% mutate(model = "CC model"),
  cv_cs_model %>% mutate(model = "CS model")
) %>% 
  ggplot(aes(1-spec, sens, group = fold_id)) +
  geom_line(alpha = 0.4) +
  geom_text(
    data = . %>% 
      group_by(model) %>% 
      summarise(
        m = mean(auc)
      ) %>% 
      mutate(
        .text = paste0(
          "AUC ", round(m,3)
        )
      ),
    aes(label = .text, x = .5, y=.3), inherit.aes = F
  ) +
  facet_wrap(~model) +
  labs(x = "1 - Specificity",
       y = "Sensitivity")

## External validation ----

development_patients <- c(
  df_sample_cs$patient_id,
  df_sample_cc$patient_id
)
n_val <- 20000
df_val_cs <- sample_cross_sectional(
  df = df_pop %>% 
    filter(!(patient_id %in% development_patients)),
  n = n_val,
  .seed = SEED
)
df_val_cc <- sample_case_control(
  df = df_pop %>% 
    filter(!(patient_id %in% development_patients)),
  n = n_val,
  cutoffs = cc_cutoffs,
  .seed = SEED
)

### CC external data ----

#### AUC/ROC ----

p_hat_cs <- predict(fit_cs, newdata = df_val_cc, type = "fitted")
p_hat_cc <- predict(fit_cc, newdata = df_val_cc, type = "fitted")
(auc_cs <- get_auc(df_val_cc$y, p_hat_cs))
(auc_cc <- get_auc(df_val_cc$y, p_hat_cc))

.roc_cc <- pROC::roc(df_val_cc$y, p_hat_cc)
.roc_cs <- pROC::roc(df_val_cc$y, p_hat_cs)
tibble(
  sens = c(
    .roc_cc$sensitivities,
    .roc_cs$sensitivities
  ),
  spec = c(
    .roc_cc$specificities,
    .roc_cs$specificities
  ),
  model = rep(c("CC model", "CS model"),
              each = length(.roc_cc$specificities))
) %>% 
  ggplot(aes(1-spec, sens, group = model)) +
  geom_line() +
  geom_text(
    data = tibble(
      model = c("CC model", "CS model"),
      .auc = c(.roc_cc$auc[1], .roc_cs$auc[1]) %>% 
        round(3),
      .text = paste0("AUC ", .auc)
    ),
    aes(label = .text, x = .5, y=.3), inherit.aes = F
  ) +
  theme_bw() +
  facet_wrap(~model) +
  labs(x = "1 - Specificity",
       y = "Sensitivity",
       title = "Case-control test set")


#### Calibration ----

png("output/manuscript/calibration-cc-test-cc-training.png",
    res = 300, width = 3600/2, height = 1600)
cal_cc <- val.prob(y = df_val_cc$y, p = p_hat_cc,
                   ylab = "Observed proportions",
                   smooth = F, statloc = F)
title(main = 'Model trained on case-control data')
dev.off()

png("output/manuscript/calibration-cc-test-cs-training.png",
    res = 300, width = 3600/2, height = 1600)
cal_cs <- val.prob(y = df_val_cc$y, p = p_hat_cs,
                   ylab = "Observed proportions",
                   smooth = F, statloc = F)
title(main = 'Model trained on cross-sectional data')
dev.off()

#### True probs ----


par(mfrow = c(1, 2))
plot(p_hat_cc, df_val_cc$p,
     xlab = "Predicted probabilities",
     ylab = "True disease probabilities",
     main = 'Model trained on case-control data')
abline(0, 1, lty = 2, col = "red", lwd = 2)
plot(p_hat_cs, df_val_cc$p,
     xlab = "Predicted probabilities",
     ylab = "True disease probabilities",
     main = 'Model trained on cross-sectional data')
abline(0, 1, lty = 2, col = "red", lwd = 2)
dev.off()

### CS external data ----

#### AUC/ROC ----

p_hat_cs <- predict(fit_cs, newdata = df_val_cs, type = "fitted")
p_hat_cc <- predict(fit_cc, newdata = df_val_cs, type = "fitted")
(auc_cs <- get_auc(df_val_cs$y, p_hat_cs))
(auc_cc <- get_auc(df_val_cs$y, p_hat_cc))

.roc_cc <- pROC::roc(df_val_cs$y, p_hat_cc)
.roc_cs <- pROC::roc(df_val_cs$y, p_hat_cs)
tibble(
  sens = c(
    .roc_cc$sensitivities,
    .roc_cs$sensitivities
  ),
  spec = c(
    .roc_cc$specificities,
    .roc_cs$specificities
  ),
  model = rep(c("CC model", "CS model"),
              each = length(.roc_cc$specificities))
) %>% 
  ggplot(aes(1-spec, sens, group = model)) +
  geom_line() +
  geom_text(
    data = tibble(
      model = c("CC model", "CS model"),
      .auc = c(.roc_cc$auc[1], .roc_cs$auc[1]) %>% 
        round(3),
      .text = paste0("AUC ", .auc)
    ),
    aes(label = .text, x = .5, y=.3), inherit.aes = F
  ) +
  theme_bw() +
  facet_wrap(~model) +
  labs(x = "1 - Specificity",
       y = "Sensitivity",
       title = "Cross-sectional test set")

#### Calibration ----


png("output/manuscript/calibration-cs-test-cc-training.png",
    res = 300, width = 3600/2, height = 1600)
cal_cc <- val.prob(y = df_val_cs$y, p = p_hat_cc,
                   ylab = "Observed proportions",
                   smooth = F, statloc = F)
title(main = 'Model trained on case-control data')
dev.off()

png("output/manuscript/calibration-cs-test-cs-training.png",
    res = 300, width = 3600/2, height = 1600)
cal_cs <- val.prob(y = df_val_cs$y, p = p_hat_cs,
                   ylab = "Observed proportions",
                   smooth = F, statloc = F)
title(main = 'Model trained on cross-sectional data')
dev.off()

#### True probs ----

## True probs

par(mfrow = c(1, 2))
plot(p_hat_cc, df_val_cs$p,
     xlab = "Predicted probabilities",
     ylab = "True probabilities",
     main = "CC model")
abline(0, 1, lty = 2, col = "red", lwd = 2)
plot(p_hat_cs, df_val_cs$p,
     xlab = "Predicted probabilities",
     ylab = "True probabilities",
     main = "CS model")
abline(0, 1, lty = 2, col = "red", lwd = 2)
par(mfrow = c(1, 1))


# Tyranny of classifiers ----

x0 <- p_hat_cc[near(df_val_cs$p, .2, 1e-4)][1]
plot(p_hat_cc, df_val_cs$p,
     xlab = "Predicted probabilities",
     ylab = "True probabilities",
     main = 'Model trained on case-control data')
abline(0, 1, lty = 2, col = "red", lwd = 2)
segments(
  x0 = x0,
  x1 = x0,
  y0 = 0,
  y1 = 0.2,
  lwd = 2, lty = 2,
  col = "blue"
)
segments(
  x0 = 0,
  x1 = x0,
  y0 = 0.2,
  y1 = 0.2,
  lwd = 2, lty = 2,
  col = "blue"
)
text(x0, 0.315, str_glue("({round(x0, 3)}, 0.20)"))


## How many people like John Doe ----
df_val_cs$p_cc <- p_hat_cc
df_val_cs$p_cs <- p_hat_cs

npv <- round(100*mean(df_val_cs$y[df_val_cs$p_cc < 0.1] == 0))
cat(
  paste0("NPV for model trained in case-control data: ", npv, "% (cutoff 10%)")
)
should_have_proportion <- round(100*mean(df_val_cs$p[df_val_cs$p_cc < 0.1] > 0.1))
df_val_cs %>% 
  filter(p_hat_cc < 0.1) %>% 
  ggplot(aes(x = p)) +
  geom_segment(
    aes(
      x = .105, xend = 0.21,
      y = 600,
      yend = 600
    ),
    color = "red"
  ) +
  annotate(
    'text',
    x = 0.16, y = 700,
    label = paste0(
      "Erroneously untreated\n",
      "(", should_have_proportion, "%)"
    ),
    color = "red",
    fontface = "bold",
    size = 4
  ) +
  geom_histogram() +
  geom_segment(
    aes(
      x = 0.1, xend = 0.1,
      y = 0, yend = 1000
    ),
    color = "black", lwd = 1.4, lty = 2
  ) +
  labs(x = latex2exp::TeX("True disease probabilities"),
       y = "Count")


# Decision consequences ----

thresholds <- seq(0.01, 0.99, 0.01)
consequences_cc <- get_decision_consequences(
  df_val_cs,
  phat = "p_cc",
  .thresholds = thresholds
)
consequences_cs <- get_decision_consequences(
  df_val_cs,
  phat = "p_cs",
  .thresholds = thresholds
) 


df_consequences <- bind_rows(consequences_cc,
                             consequences_cs) %>% 
  select(-contains('pv')) %>% 
  pivot_longer(
    cols = -c(threshold, estimator)
  ) %>% mutate(
    estimator = case_when(
      estimator == "p_cc" ~ "CC model",
      estimator == "p_cs" ~ "CS model",
      T ~ NA_character_
    )
  )

### undertreatment ----


df_undertreat <- df_consequences %>% 
  filter(name == "undertreat")

undertreat_at_t10 <- df_undertreat %>% 
  pivot_wider(
    names_from = estimator,
    values_from = value
  ) %>% 
  filter(threshold == 0.1)
undertreat_at_t10

.p_under <- df_undertreat %>% 
  ggplot(aes(threshold, value, color = estimator)) +
  geom_line(lwd = 2, lty = 1) +
  scale_y_continuous(labels = scales::percent,
                     limits = c(0,1)) +
  scale_x_continuous(labels = scales::percent) +
  labs(
    x = "Decision threshold", 
    y = "Undertreatment  probability",
    color = NULL
  ) +
  theme(
    legend.position = c(0.7, 0.8)
  ) +
  scale_color_discrete(
    label = list(
      "CC model" = "Model trained on case-control data",
      "CS model" = "Model trained on cross-sectional data"
    )
  )
ggsave("output/manuscript/undertreatment-probability.png", 
       .p_under, width = 7.5, height = 4, dpi = 600)

### Overtreatment ----


df_overtreat <- df_consequences %>% 
  filter(name == "overtreat")

overtreat_at_t50 <- df_overtreat %>% 
  pivot_wider(
    names_from = estimator,
    values_from = value
  ) %>% 
  filter(threshold == 0.5)
overtreat_at_t50

.p_over <- df_overtreat %>% 
  ggplot(aes(threshold, value, color = estimator)) +
  geom_line(lwd = 2, lty = 1) +
  scale_y_continuous(labels = scales::percent,
                     limits = c(0,1)) +
  scale_x_continuous(labels = scales::percent) +
  labs(
    y = "Overtreatment  probability",
    x = "Decision threshold", 
    y = NULL,
    color = NULL
  ) +
  theme(
    legend.position = "none"
  ) +
  scale_color_discrete(
    label = list(
      "CC model" = "Model trained on case-control data",
      "CS model" = "Model trained on cross-sectional data"
    )
  )
ggsave("output/manuscript/overtreatment-probaility.png", 
       .p_over, width = 7.5, height = 4, dpi = 600)

# Recalibration ----

df_val_cs$lp_cc <- logit(df_val_cs$p_cc) - coef(fit_cc)[1]
set.seed(SEED)
ix <- sample(1:nrow(df_val_cs), 1000)
fit_recal <- glm(y ~ lp_cc, data = df_val_cs[ix,], family = binomial())

fit_cc_updated <- fit_cc
fit_cc_updated$coefficients[1] <- coef(fit_recal)[1]
fit_cc_updated$coefficients[-1] <- coef(fit_recal)[2]*coef(fit_cc)[-1]
df_val_cs$p_cc_updated <-  predict(fit_cc_updated, newdata = df_val_cs, 
                                   type = "fitted")
# AUC is unchanged
pROC::auc(df_val_cs$y[-ix], df_val_cs$p_cc_updated[-ix])

## True probs ----

png("output/manuscript/recalibration-cc-training-true-probs.png",
    res = 300, width = 3600/2, height = 1600)
plot(df_val_cs$p_cc_updated[-ix], df_val_cs$p[-ix],
     xlab = "Predicted probabilities",
     ylab = "True disease probabilities",
     main = "CC model updated")
abline(0, 1, lty = 2, col = "red", lwd = 2)
dev.off()

## Undertreatment ----
consequences_cc_updated <- get_decision_consequences(
  df_val_cs[-ix,],
  phat = "p_cc_updated",
  .thresholds = thresholds
)

consequences_cc_updated %>% 
  filter(threshold == 0.1) %>% 
  select(contains('treat'))


