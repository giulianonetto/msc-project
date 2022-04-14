# Load packages ----------------------------------------------------------------
library(tidyverse)
library(ggpubr)
library(patchwork)
library(rms)
library(rmda)
library(pROC)
#library(caret) # only createFolds function used
install.packages('predtools'); library(predtools)
theme_set(theme_bw())
source("R/functions.R") 
# create output dir for figures and stuff
dir.create("output/report", showWarnings = FALSE)

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
.p <- plot_disease_prob(df_pop, title = .title)
ggsave("output/report/pop-risk-distribution.png", 
       .p, width = 10, height = 4, bg = 'white')

round(mean(between(df_pop[df_pop$y == 0, "p"], 0.8, 1))*100,2)
round(mean(between(df_pop[df_pop$y == 1, "p"], 0, 0.3))*100,2) # 20% of cases have risk btw 0 and 30%

# Samples ----

## cross-sectional ----
n_sample <- 2e4
df_sample_cs <- sample_cross_sectional(df = df_pop, n = n_sample, .seed = SEED)
.title <- paste0(
  "max AUC ", get_auc(df_sample_cs$y, df_sample_cs$p)
)
.p <- plot_disease_prob(df = df_sample_cs, title = .title)
ggsave("output/report/sample-cs-risk-distribution.png", 
       .p, width = 10, height = 4, bg = 'white')

round(mean(between(df_sample_cs[df_sample_cs$y == 0, "p"], 0.8, 1))*100, 2)
round(mean(between(df_sample_cs[df_sample_cs$y == 1, "p"], 0, 0.3))*100,2)
## Case-control ----

### mechanism ----

vlines <- list(
  geom_vline(xintercept = c(.3, .7), lty = 2, lwd = 1.5, alpha = 1,
             color = c("magenta", "skyblue"))
)
.p <- ggarrange(
  plot_disease_prob_stacked(df_sample_cs) + vlines,
  plot_disease_prob_nonstacked(df_sample_cs, .geom = vlines),
  ncol = 2
) + bgcolor("white")
ggsave("output/report/sample-cc-mechanism-risk-distribution.png", 
       .p, width = 10, height = 4)

### sample ----
cc_cutoffs <- c(0.3, 0.7)
df_sample_cc <- sample_case_control(df = df_pop, n = n_sample, 
                                    cutoffs = cc_cutoffs,
                                    .seed = SEED)
.title <- paste0(
  "max AUC ", get_auc(df_sample_cc$y, df_sample_cc$p)
)
.p <- plot_disease_prob(df = df_sample_cc, title = .title)
ggsave("output/report/sample-cc-risk-distribution.png", 
       .p, width = 10, height = 4)

## Group differences ----

get_avg_diff <- function(.data, conf.int = TRUE) {
  .data %>% 
    pivot_longer(contains('x')) %>% 
    group_by(name = str_replace(name, "x", "Predictor ")) %>% 
    group_modify(~broom::tidy(lm(value ~ y, data = .x), conf.int=T)) %>% 
    filter(term == "y") %>% 
    select(name, estimate, conf.low, conf.high) %>% 
    mutate(
      label = ifelse(
        isTRUE(conf.int),
        paste0(
          "Avg. diff.\n",
          round(estimate, 2),
          " [",
          round(conf.low, 2),
          " \u2012 ",
          round(conf.high, 2),
          "]"
        ),
        paste0(
          "Avg. diff.\n",
          round(estimate, 2)
        )
      )
    )
}

gdiff_pop <- df_pop %>% 
  pivot_longer(cols = contains('x')) %>% 
  mutate(name = str_replace(name, "x", "Predictor "),
         y = paste0("Y = ", y)) %>% 
  ggplot(aes(x = value, group=y, fill = y)) +
  geom_density(alpha = .7) +
  facet_wrap(~name, ncol=3) +
  geom_text(
    data = get_avg_diff(df_pop, conf.int = FALSE),
    aes(x = 5.5, y = 0.3, label = label),
    inherit.aes = F
  ) +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(size = 18, face = 'bold'),
        plot.title.position = "plot",
        plot.title = element_text(size = 14, face = 'bold')) +
  scale_fill_brewer(type = 'qual', palette = 2) +
  labs(x = NULL, fill = NULL, y = "density",
       title = "Population") +
  guides(fill = 'none') +
  coord_cartesian(xlim = c(-9, 9))

gdiff_cs <- df_sample_cs %>% 
  pivot_longer(cols = contains('x')) %>% 
  mutate(name = str_replace(name, "x", "Predictor "),
         y = paste0("Y = ", y)) %>% 
  ggplot(aes(x = value, group=y, fill = y)) +
  geom_density(alpha = .7) +
  geom_text(
    data = get_avg_diff(df_sample_cs),
    aes(x = 5.5, y = 0.3, label = label),
    inherit.aes = F
  ) +
  facet_wrap(~name, ncol=3) +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_blank(),
        plot.title.position = "plot",
        plot.title = element_text(size = 14, face = 'bold')) +
  scale_fill_brewer(type = 'qual', palette = 2) +
  labs(x = NULL, fill = NULL, y = "density",
       title = 'Cross-sectional sample') +
  coord_cartesian(xlim = c(-9, 9))

gdiff_cc <- df_sample_cc %>% 
  pivot_longer(cols = contains('x')) %>% 
  mutate(name = str_replace(name, "x", "Predictor "),
         y = paste0("Y = ", y)) %>% 
  ggplot(aes(x = value, group=y, fill = y)) +
  geom_density(alpha = .7) +
  geom_text(
    data = get_avg_diff(df_sample_cc),
    aes(x = 5.5, y = 0.3, label = label),
    inherit.aes = F
  ) +
  facet_wrap(~name, ncol=3) +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_blank(),
        plot.title.position = "plot",
        plot.title = element_text(size = 14, face = 'bold')) +
  scale_fill_brewer(type = 'qual', palette = 2) +
  labs(x = NULL, fill = NULL, y = "density",
       title = 'Case-control sample') +
  guides(fill = 'none') +
  coord_cartesian(xlim = c(-9, 9))


.p <- gdiff_pop/gdiff_cs/gdiff_cc 

ggsave("output/report/group-differences.png", 
       .p, width = 14, height = 10)

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

cv_roc <- bind_rows(
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
ggsave("output/report/cv_roc.png",
       cv_roc, width = 10, height = 4)
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
.p <- tibble(
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
ggsave("output/report/test_set_roc_cc.png",
       .p, width = 10, height = 4)

#### Calibration ----

png("output/report/calibration-cc-data.png",
    res = 300, width = 3600, height = 1600)
par(mfrow = c(1, 2))
cal_cc <- val.prob(y = df_val_cc$y, p = p_hat_cc,
                   ylab = "Observed probabilities",
                   smooth = F, statloc = F)
title(main = 'CC model')
cal_cs <- val.prob(y = df_val_cc$y, p = p_hat_cs,
                   ylab = "Observed probabilities",
                   smooth = F, statloc = F)
title(main = 'CS model')
dev.off()

#### True probs ----

png("output/report/calibration-cc-data-true-probs.png",
    res = 300, width = 3600, height = 1600)
par(mfrow = c(1, 2))
plot(p_hat_cc, df_val_cc$p,
     xlab = "Predicted probabilities",
     ylab = "True underlying probabilities",
     main = "CC model")
abline(0, 1, lty = 2, col = "red", lwd = 2)
plot(p_hat_cs, df_val_cc$p,
     xlab = "Predicted probabilities",
     ylab = "True underlying probabilities",
     main = "CS model")
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
.p <- tibble(
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
ggsave("output/report/test_set_roc_cs.png",
       .p, width = 10, height = 4)

#### Calibration ----

png("output/report/calibration-cs-data.png",
    res = 300, width = 3600, height = 1600)
par(mfrow = c(1, 2))
cal_cc <- val.prob(y = df_val_cs$y, p = p_hat_cc,
                   ylab = "Observed probabilities",
                   smooth = F, statloc = F)
title(main = 'CC model')
cal_cs <- val.prob(y = df_val_cs$y, p = p_hat_cs,
                   ylab = "Observed probabilities",
                   smooth = F, statloc = F)
title(main = 'CS model')
dev.off()

#### True probs ----

## True probs

png("output/report/calibration-cs-data-true-probs.png",
    res = 300, width = 3600, height = 1600)
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
dev.off()

#### mROC ----

{
  png("output/report/mROC.png",
      res = 300, width = 3600, height = 1600)
  par(mfrow = c(1, 2))
  r <- roc(df_val_cs$y, p_hat_cs)
  plot(r, main = "CS model")
  r2 <- roc(df_sample_cs$y, predict(fit_cs, type = 'fitted'))
  plot(r2, add = TRUE, col = 'blue')
  lines(mROC(p_hat_cs), col = "red", lty = 2, lwd = 3)
  legend("bottomright", legend = c("ROC", "ROC training", "mROC"), 
         col = c('black', 'blue', 'red'), 
         lty = c(1, 1, 2), cex = 0.7)
  
  r <- roc(df_val_cs$y, p_hat_cc)
  plot(r, main = "CC model")
  r2 <- roc(df_sample_cc$y, predict(fit_cc, type = 'fitted'))
  plot(r2, add = TRUE, col = 'blue')
  lines(mROC(p_hat_cc), col = "red", lty = 2, lwd = 3)
  legend("bottomright", legend = c("ROC", "ROC training", "mROC"), 
         col = c('black', 'blue', 'red'), 
         lty = c(1, 1, 2), cex = 0.7)
  dev.off()
}

#### DCA ----
library(dcurves)
decision_curve_analisys <- dca(
  y ~ CS_model + CC_model,
  data = data.frame(
    y = df_val_cs$y,
    CS_model = p_hat_cs, 
    CC_model = p_hat_cc
  ), 
  label = list(CS_model = "CS model", CC_model = "CC model")
)
.p <- decision_curve_analisys %>% 
  plot(smooth = F) +
  geom_vline(xintercept = .1, lty = 2, alpha = 0.75) +
  geom_line(size = 1.5) +
  labs(x = "Decision threshold")

ggsave("output/report/dca-cs-data.png",
       .p, width = 8.5, height = 4)

decision_curve_analisys$dca %>% 
  filter(threshold == 0.1)
## Distribution of predictions

.p <- data.frame(
  model = rep(c("CC model", "CS model"),
              each = n_val),
  p = c(p_hat_cc, p_hat_cs),
  y = factor(c(df_val_cs$y, df_val_cs$y))
) %>% 
  ggplot(aes(x = p)) +
  geom_histogram(aes(fill = y),
                 position = "stack") +
  facet_wrap(~ model) +
  labs(x = "Predicted probabilities",
       fill = "Has clotting issues") +
  theme(legend.position = "top")

ggsave("output/report/predicted-probs-distribution.png",
       .p, width = 8.5, height = 4)

## ROC curves

png("output/report/roc-curves-cs-data.png",
    res = 300, width = 3600, height = 1600)
par(mfrow = c(1,2))
pROC::plot.roc(df_val_cs$y, p_hat_cs, main = "CS model")
pROC::plot.roc(df_val_cs$y, p_hat_cc, main = "CC model")
dev.off()

# Tyranny of classifiers ----

x0 <- p_hat_cc[near(df_val_cs$p, .2, 1e-4)][1]
png("output/report/tyranny-of-classifiers.png",
    res = 300, width = 2600, height = 1600)
par(mfrow = c(1,1))
plot(p_hat_cc, df_val_cs$p,
     xlab = "Predicted probabilities",
     ylab = "True probabilities",
     main = "CC model")
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
dev.off()

## How many people like John Doe ----
df_val_cs$p_cc <- p_hat_cc
df_val_cs$p_cs <- p_hat_cs

npv <- round(100*mean(df_val_cs$y[df_val_cs$p_cc < 0.1] == 0))
should_have_proportion <- round(100*mean(df_val_cs$p[df_val_cs$p_cc < 0.1] > 0.1))
.p <- df_val_cs %>% 
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
  labs(x = latex2exp::TeX("True underlying probabilities"),
       y = "Count")

ggsave("output/report/how-many-people-like-john-doe.png", 
       .p, width = 7, height = 4, dpi = 600)

## Decision consequences ----


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

.p_under <- df_undertreat %>% 
  ggplot(aes(threshold, value, color = estimator)) +
  geom_line(lwd = 2, lty = 1) +
  geom_segment(
    aes(
      x = .1, xend = .1,
      y = 0,
      yend = undertreat_at_t10[['CC model']]
    ), 
    color = "red", lwd = 1.1, lty = 2
  ) +
  geom_segment(
    aes(
      y = undertreat_at_t10[['CC model']], 
      yend = undertreat_at_t10[['CC model']],
      x = 0,
      xend = .1
    ), 
    color = "red", lwd = 1.1, lty = 2
  ) +
  annotate(
    'text', x = .6, y = .45,
    label = paste0(
      "Every 100 patients not treated using threshold of 10%,\nabout ",
      round(undertreat_at_t10[['CC model']]*100),
      " actually should have been treated"
    ),
    color = "red",
    fontface = "bold",
    size = 3.5
  ) +
  scale_y_continuous(labels = scales::percent,
                     limits = c(0,1)) +
  scale_x_continuous(labels = scales::percent) +
  labs(
    # title = "Fail-to-treat probability",
    title = "Undertreatment probability",
    x = "Decision threshold", 
    y = NULL,
    color = NULL
  )
ggsave("output/report/how-many-john-does.png", 
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

.p_over <- df_overtreat %>% 
  ggplot(aes(threshold, value, color = estimator)) +
  geom_line(lwd = 2, lty = 1) +
  geom_segment(
    aes(
      x = .5, xend = .5,
      y = 0,
      yend = overtreat_at_t50[['CC model']]
    ), 
    color = "red", lwd = 1.1, lty = 2
  ) +
  geom_segment(
    aes(
      y = overtreat_at_t50[['CC model']], 
      yend = overtreat_at_t50[['CC model']],
      x = 0,
      xend = .5
    ), 
    color = "red", lwd = 1.1, lty = 2
  ) +
  annotate(
    'text', x = .4, y = .65,
    label = paste0(
      "Every 100 patients treated using threshold of 50%,\nabout ",
      round(overtreat_at_t50[['CC model']]*100),
      " actually should not have been treated"
    ),
    color = "red",
    fontface = "bold",
    size = 3.5
  ) +
  scale_y_continuous(labels = scales::percent,
                     limits = c(0,1)) +
  scale_x_continuous(labels = scales::percent) +
  labs(
    # title = "Fail-to-treat probability",
    title = "Overtreatment probability",
    x = "Decision threshold", 
    y = NULL,
    color = NULL
  )
ggsave("output/report/overtreatment-probaility.png", 
       .p_over, width = 7.5, height = 4, dpi = 600)


.p <- (.p_under/.p_over) + plot_layout(guides = 'collect')
ggsave("output/report/under-and-over-treatment-probaility.png", 
       .p, width = 7.5, height = 8, dpi = 600)

### undertreat and sick ----

df_undertreat_sick <- df_consequences %>% 
  filter(name == "undertreat_sick")

undertreat_sick_at_t10 <- df_undertreat_sick %>% 
  pivot_wider(
    names_from = estimator,
    values_from = value
  ) %>% 
  filter(threshold == 0.1)

.p <- df_undertreat_sick %>% 
  ggplot(aes(threshold, value, color = estimator)) +
  geom_segment(
    aes(
      x = .1, xend = .1,
      y = 0,
      yend = undertreat_sick_at_t10[['CC model']]
    ), 
    color = "red", lwd = 1.1, lty = 3
  ) +
  geom_segment(
    aes(
      y = undertreat_sick_at_t10[['CC model']], 
      yend = undertreat_sick_at_t10[['CC model']],
      x = 0,
      xend = .1
    ), 
    color = "red", lwd = 1.1, lty = 3
  ) +
  annotate(
    'text', x = .3, y = .04,
    label = paste0(
      "Every 100 patients not treated using threshold of 10%,\nabout ",
      round(undertreat_sick_at_t10[['CC model']]*100),
      " actually should have been treated\n",
      "(all of which can sue us)"
    ),
    color = "red",
    fontface = "bold",
    size = 3
  ) +
  geom_line(lwd = 2, lty = 2) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(labels = scales::percent) +
  labs(
    # title = "Fail-to-treat probability",
    title = "Percentage of people not treated who *should* have been treated",
    subtitle = "P(true probability > t | prediction < t)",
    x = "Threshold (t)", 
    y = NULL,
    color = NULL
  )
ggsave("output/report/how-many-john-does-sick.png", 
       .p, width = 7, height = 4, dpi = 600)