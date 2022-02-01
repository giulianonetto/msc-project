library(tidyverse)
library(sjPlot)
library(rms)
theme_set(theme_bw())

df <- read_csv("data/df_giu.csv") %>% 
  filter(
    plasma_ptau181 < quantile(plasma_ptau181, probs = .975)
  )
df$y <- as.numeric(df$DX == "Dementia")
dd <- datadist(df)
options(datadist = "dd")

pROC::roc(df$y, df$plasma_ptau181) %>% 
  pROC::ggroc()
pROC::auc(df$y, df$plasma_ptau181)
cutoff <- 20
(ppv <- mean(df[df$plasma_ptau181 > cutoff, "DX"] == "Dementia"))


fit <- lrm(y ~ rcs(plasma_ptau181, 3),
           data = df, x = T, y = T)

sjPlot::plot_model(fit, "pred", terms = "plasma_ptau181 [all]") +
  geom_vline(xintercept = 20) +
  geom_hline(yintercept = ppv)
df %>% 
  ggplot(aes(x = AV45)) +
  geom_histogram() +
  facet_wrap(~y)
ecdf(df$plasma_ptau181[df$y==0])(30) - ecdf(df$plasma_ptau181[df$y==0])(10)
df$p_hat <- predict(fit, type = 'fitted')
df$outcomes <- df$y
threshperf_plot(df, outcome = "outcomes",
                prediction = "p_hat")
predict(fit, newdata = data.frame(plasma_ptau181 = cutoff), type='fitted')
