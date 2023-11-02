####################
## Description:
##  - In this file we compare characteristics and outcomes of semaglutide vs dulaglutide/liraglutide
####################


## Load libraries
library(tidyverse)
library(ggplot2)
library(PSweight)
library(patchwork)

## Set up directory path to save files (stagered to ensure folders are created)

dir.create("Samples")
dir.create("Samples/Sema_vs_Lira_Dula")

output_path <- "Samples/Sema_vs_Lira_Dula/Aurum"
dir.create(output_path)

## make directory for outputs
dir.create("Plots")

## male directory for outputs
dir.create(paste0(output_path, "/semaglutide"))


###############################################################################
###############################################################################
############################### Read Data In ##################################
###############################################################################
###############################################################################

## Load functions required

source("01.slade_aurum_set_data.R")

###############################################################################
###############################################################################
########################## General variables ##################################
###############################################################################
###############################################################################


## Load dataset
dataset <- set_up_data_sglt2_glp1(dataset.type="full.cohort") %>%
  filter(drugsubstances %in% c("Dulaglutide", "Liraglutide", "Canagliflozin", "Dapagliflozin", "Empagliflozin", "Semaglutide")) %>%
  filter(!is.na(stopdrug_6m_3mFU)) %>%
  mutate(drugclass = ifelse(drugclass == "SGLT2", "SGLT2i", "GLP1-RA"),
         drugclass = factor(drugclass, levels = c("SGLT2i", "GLP1-RA")),
         ncurrtx = factor(ncurrtx, levels = c("1", "2", "3", "4", "5+"), labels = c("1", "2", "3", "4", "4"))) %>%
  filter(yrdrugstart > 2018)



###############################################################################
###############################################################################
########################## Match individuals ##################################
###############################################################################
###############################################################################


###    HbA1c outcomes
#:----------------------------------------------

dataset_formatted_hba1c <- dataset %>%
  filter(!is.na(posthba1cfinal)) %>%
  filter(!is.na(prehba1c)) %>%
  mutate(drugclass_study = ifelse(drugsubstances %in% c("Canagliflozin", "Dapagliflozin", "Empagliflozin"), "SGLT2i", 
                                  ifelse(drugsubstances %in% c("Dulaglutide", "Liraglutide"), "GLP1-RA", drugsubstances)),
         drugclass_study = factor(drugclass_study))



h.pscores <- SumStat(ps.formula = drugclass_study ~ agetx + sex + t2dmduration + drugline + ncurrtx + prehba1c,
                     data = dataset_formatted_hba1c,
                     weight = c("IPW", "overlap"))


lm.hba1c <- lm(formula = posthba1cfinal ~ drugclass_study + prehba1c + hba1cmonth,
               data = dataset_formatted_hba1c,
               weights = h.pscores$pw.weights$overlap)

patient.prediction <- data.frame(
  prehba1c = 75,
  drugclass_study = unique(dataset_formatted_hba1c$drugclass_study),
  hba1cmonth = 6
)

hba1c_estimates <- as.matrix(predict(lm.hba1c, patient.prediction, interval = "confidence") - 75) %>%
  as.data.frame() %>%
  cbind(drugclass_study = patient.prediction$drugclass_study) %>%
  left_join(
    dataset_formatted_hba1c %>%
      select(drugclass_study, drugclass) %>%
      group_by(drugclass_study) %>%
      mutate(n.patients = n()) %>%
      ungroup() %>%
      unique(),
    by = "drugclass_study"
  ) %>%
  mutate(drugclass = factor(drugclass, levels = c("SGLT2i", "GLP1-RA")))


plot_hba1c_estimates <- hba1c_estimates %>%
  ggplot(aes(x = drugclass_study, y = fit)) +
  geom_point() +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.5) +
  facet_grid(~drugclass, scales = "free_x") +
  theme_bw() +
  labs(
    y = "Adjusted Response (%)"
  ) +
  theme(axis.title.x = element_blank())




###    Weight outcomes
#:----------------------------------------------


dataset_formatted_weight <- dataset %>%
  filter(!is.na(postweight6m)) %>%
  filter(!is.na(preweight)) %>%
  mutate(drugclass_study = ifelse(drugsubstances %in% c("Canagliflozin", "Dapagliflozin", "Empagliflozin"), "SGLT2i", 
                                  ifelse(drugsubstances %in% c("Dulaglutide", "Liraglutide"), "GLP1-RA", drugsubstances)),
         drugclass_study = factor(drugclass_study))



w.pscores <- SumStat(ps.formula = drugclass_study ~ agetx + sex + t2dmduration + drugline + ncurrtx + preweight,
                     data = dataset_formatted_weight,
                     weight = c("IPW", "overlap"))


lm.weight <- lm(formula = postweight6m ~ drugclass_study + preweight,
                data = dataset_formatted_weight,
                weights = w.pscores$pw.weights$overlap)


patient.prediction <- data.frame(
  preweight = 95,
  drugclass_study = unique(dataset_formatted_weight$drugclass_study)
)


weight_estimates <- as.matrix(predict(lm.weight, patient.prediction, interval = "confidence")-95) %>%
  as.data.frame() %>%
  cbind(drugclass_study = patient.prediction$drugclass_study) %>%
  left_join(
    dataset_formatted_weight %>%
      select(drugclass_study, drugclass) %>%
      group_by(drugclass_study) %>%
      mutate(n.patients = n()) %>%
      ungroup() %>%
      unique(),
    by = "drugclass_study"
  ) %>%
  mutate(drugclass = factor(drugclass, levels = c("SGLT2i", "GLP1-RA")))


plot_weight_estimates <- weight_estimates %>%
  ggplot(aes(x = drugclass_study, y = fit)) +
  geom_point() +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.5) +
  facet_grid(~drugclass, scales = "free_x") +
  theme_bw() +
  labs(
    y = "Adjusted Response (kg)"
  ) +
  theme(axis.title.x = element_blank())



lm.weight.sex <- lm(formula = postweight6m ~ drugclass_study + preweight + drugclass_study*sex,
                    data = dataset_formatted_weight,
                    weights = h.pscores$pw.weights$overlap)


patient.prediction <- expand.grid(
  preweight = 95,
  drugclass_study = unique(dataset_formatted_weight$drugclass_study),
  sex = c("Male", "Female")
)


weight_estimates_sex <- as.matrix(predict(lm.weight.sex, patient.prediction, interval = "confidence")-95) %>%
  as.data.frame() %>%
  cbind(drugclass_study = patient.prediction$drugclass_study,
        sex = patient.prediction$sex) %>%
  left_join(
    dataset_formatted_weight %>%
      select(drugclass_study, drugclass, sex) %>%
      group_by(drugclass_study, sex) %>%
      mutate(n.patients = n()) %>%
      ungroup() %>%
      unique(),
    by = c("drugclass_study", "sex")
  ) %>%
  mutate(drugclass = factor(drugclass, levels = c("SGLT2i", "GLP1-RA")))


###    Discontinuation outcomes
#:----------------------------------------------


dataset_formatted_discontinuation <- dataset %>%
  filter(!is.na(stopdrug_6m_3mFU)) %>%
  mutate(drugclass_study = ifelse(drugsubstances %in% c("Canagliflozin", "Dapagliflozin", "Empagliflozin"), "SGLT2i", 
                                  ifelse(drugsubstances %in% c("Dulaglutide", "Liraglutide"), "GLP1-RA", drugsubstances)),
         drugclass_study = factor(drugclass_study))



d.pscores <- SumStat(ps.formula = drugclass_study ~ agetx + sex + t2dmduration + drugline + ncurrtx,
                     data = dataset_formatted_discontinuation,
                     weight = c("IPW", "overlap"))


lm.discontinuation <- lm(formula = stopdrug_6m_3mFU ~ drugclass_study,
                         data = dataset_formatted_discontinuation,
                         weights = w.pscores$pw.weights$overlap)


patient.prediction <- data.frame(
  drugclass_study = unique(dataset_formatted_discontinuation$drugclass_study)
)


discontinuation_estimates <- predict(lm.discontinuation, patient.prediction, interval = "confidence") %>%
  as.data.frame() %>%
  cbind(drugclass_study = patient.prediction$drugclass_study) %>%
  left_join(
    dataset_formatted_discontinuation %>%
      select(drugclass_study, drugclass) %>%
      group_by(drugclass_study) %>%
      mutate(n.patients = n()) %>%
      ungroup() %>%
      unique(),
    by = "drugclass_study"
  ) %>%
  mutate(drugclass = factor(drugclass, levels = c("SGLT2i", "GLP1-RA")))


plot_discontinuation_estimates <- discontinuation_estimates %>%
  ggplot(aes(x = drugclass_study, y = fit)) +
  geom_point() +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.5) +
  facet_grid(~drugclass, scales = "free_x") +
  theme_bw() +
  scale_y_continuous(labels=scales::percent) +
  labs(
    y = "Adjusted Response (%)"
  ) +
  theme(axis.title.x = element_blank())































m.hba1c <- matchit(formula = drugclass_study ~ agetx + sex + t2dmduration + prehospitalisation + drugline + ncurrtx + prehba1c , # removed due to missingness: ethnicity
                   data = dataset_formatted_hba1c)

## check standardised mean differences
plot(summary(m.hba1c), var.order = "unmatched")

dataset_formatted_hba1c <- dataset_formatted_hba1c %>%
  cbind(
    p.scores = m.hba1c$distance,
    matched = m.hba1c$weights
  )


## distribution of propensity scores
plot_pre_pscores <- dataset_formatted_hba1c %>%
  ggplot() +
  geom_density(aes(x = p.scores, fill = drugclass_study), alpha = 0.7) +
  labs(
    title = "Pre-match propensity scores",
    fill = "Subtype"
  ) +
  theme(
    legend.position = "bottom"
  )

plot_post_pscores <- dataset_formatted_hba1c %>%
  filter(matched == 1) %>%
  ggplot() +
  geom_density(aes(x = p.scores, fill = drugclass_study), alpha = 0.7) +
  labs(
    title = "Post-match propensity scores",
    fill = "Subtype"
  ) +
  theme(
    legend.position = "bottom"
  )


plot_pscores_hba1c <- plot_pre_pscores + plot_post_pscores +
  plot_layout(guides = 'collect') &
  theme(
    legend.position = "bottom"
  )


### regression
lm.hba1c <- lm(formula = posthba1cfinal ~ drugclass_study + prehba1c + hba1cmonth,
               data = dataset_formatted_hba1c %>%
                 filter(matched == 1))


###    Weight outcomes
#:----------------------------------------------

dataset_formatted_weight <- dataset_formatted %>%
  filter(!is.na(postweight6m)) %>%
  filter(!is.na(preweight)) %>%
  mutate(drugclass_study = ifelse(drugsubstances == "Semaglutide", "Semaglutide", "Liraglutide/Dulaglutide"),
         drugclass_study = factor(drugclass_study))


m.weight <- matchit(formula = drugclass_study ~ agetx + sex + t2dmduration + prehospitalisation + drugline + ncurrtx + preweight , # removed due to missingness: ethnicity
                   data = dataset_formatted_weight)

## check standardised mean differences
plot(summary(m.weight), var.order = "unmatched")

dataset_formatted_weight <- dataset_formatted_weight %>%
  cbind(
    p.scores = m.weight$distance,
    matched = m.weight$weights
  )


plot_pre_pscores <- dataset_formatted_weight %>%
  ggplot() +
  geom_density(aes(x = p.scores, fill = drugclass_study), alpha = 0.7) +
  labs(
    title = "Pre-match propensity scores",
    fill = "Subtype"
  ) +
  theme(
    legend.position = "bottom"
  )

plot_post_pscores <- dataset_formatted_weight %>%
  filter(matched == 1) %>%
  ggplot() +
  geom_density(aes(x = p.scores, fill = drugclass_study), alpha = 0.7) +
  labs(
    title = "Post-match propensity scores",
    fill = "Subtype"
  ) +
  theme(
    legend.position = "bottom"
  )


plot_pscores_weight <- plot_pre_pscores + plot_post_pscores +
  plot_layout(guides = 'collect') &
  theme(
    legend.position = "bottom"
  )



### regression
lm.weight <- lm(formula = postweight6m ~ drugclass_study + preweight,
               data = dataset_formatted_weight %>%
                 filter(matched == 1))




###    eGFR outcomes
#:----------------------------------------------

dataset_formatted_egfr <- dataset_formatted %>%
  filter(!is.na(preegfr)) %>%
  filter(!is.na(postegfr6m)) %>%
  mutate(drugclass_study = ifelse(drugsubstances == "Semaglutide", "Semaglutide", "Liraglutide/Dulaglutide"),
         drugclass_study = factor(drugclass_study))


m.egfr <- matchit(formula = drugclass_study ~ agetx + sex + t2dmduration + prehospitalisation + drugline + ncurrtx + preegfr , # removed due to missingness: ethnicity
                    data = dataset_formatted_egfr)

## check standardised mean differences
plot(summary(m.egfr), var.order = "unmatched")

dataset_formatted_egfr <- dataset_formatted_egfr %>%
  cbind(
    p.scores = m.egfr$distance,
    matched = m.egfr$weights
  )


plot_pre_pscores <- dataset_formatted_egfr %>%
  ggplot() +
  geom_density(aes(x = p.scores, fill = drugclass_study), alpha = 0.7) +
  labs(
    title = "Pre-match propensity scores",
    fill = "Subtype"
  ) +
  theme(
    legend.position = "bottom"
  )

plot_post_pscores <- dataset_formatted_egfr %>%
  filter(matched == 1) %>%
  ggplot() +
  geom_density(aes(x = p.scores, fill = drugclass_study), alpha = 0.7) +
  labs(
    title = "Post-match propensity scores",
    fill = "Subtype"
  ) +
  theme(
    legend.position = "bottom"
  )


plot_pscores_egfr <- plot_pre_pscores + plot_post_pscores +
  plot_layout(guides = 'collect') &
  theme(
    legend.position = "bottom"
  )



### regression
lm.egfr <- lm(formula = postegfr6m ~ drugclass_study + preegfr,
                data = dataset_formatted_egfr %>%
                  filter(matched == 1))



###    Discontinuation outcomes
#:----------------------------------------------

dataset_formatted_discontinuation <- dataset_formatted %>%
  filter(!is.na(stopdrug_6m_3mFU)) %>%
  filter(!is.na(prehba1c)) %>%
  filter(!is.na(preweight)) %>%
  mutate(drugclass_study = ifelse(drugsubstances == "Semaglutide", "Semaglutide", "Liraglutide/Dulaglutide"),
         drugclass_study = factor(drugclass_study))


m.discontinuation <- matchit(formula = drugclass_study ~ agetx + sex + t2dmduration + prehospitalisation + drugline + ncurrtx + prehba1c + preweight, # removed due to missingness: ethnicity
                             data = dataset_formatted_discontinuation)

## check standardised mean differences
plot(summary(m.discontinuation), var.order = "unmatched")

dataset_formatted_discontinuation <- dataset_formatted_discontinuation %>%
  cbind(
    p.scores = m.discontinuation$distance,
    matched = m.discontinuation$weights
  )


plot_pre_pscores <- dataset_formatted_discontinuation %>%
  ggplot() +
  geom_density(aes(x = p.scores, fill = drugclass_study), alpha = 0.7) +
  labs(
    title = "Pre-match propensity scores",
    fill = "Subtype"
  ) +
  theme(
    legend.position = "bottom"
  )

plot_post_pscores <- dataset_formatted_discontinuation %>%
  filter(matched == 1) %>%
  ggplot() +
  geom_density(aes(x = p.scores, fill = drugclass_study), alpha = 0.7) +
  labs(
    title = "Post-match propensity scores",
    fill = "Subtype"
  ) +
  theme(
    legend.position = "bottom"
  )


plot_pscores_discontinuation <- plot_pre_pscores + plot_post_pscores +
  plot_layout(guides = 'collect') &
  theme(
    legend.position = "bottom"
  )



### regression
lm.discontinuation <- lm(formula = stopdrug_6m_3mFU ~ drugclass_study,
              data = dataset_formatted_discontinuation %>%
                filter(matched == 1))



















