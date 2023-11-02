####################
## Description:
##  - In this file we compare outcomes between dulaglutide and liraglutide and SGLT2 drugs
####################


## Load libraries
library(tidyverse)
library(tableone)
library(PSweight)
library(patchwork)
library(scales)

## Set up directory path to save files (stagered to ensure folders are created)

dir.create("Samples")
dir.create("Samples/Sema_vs_Lira_Dula")

output_path <- "Samples/Sema_vs_Lira_Dula/Aurum"
dir.create(output_path)

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
dataset <- set_up_data_sglt2_glp1(dataset.type="full.cohort", diagnosis = FALSE) %>%
  filter(drugsubstances %in% c("Dulaglutide", "Liraglutide", "Canagliflozin", "Dapagliflozin", "Liraglutide", "Empagliflozin", "Semaglutide")) %>%
  filter(!is.na(stopdrug_6m_3mFU)) %>%
  mutate(drugclass = ifelse(drugclass == "SGLT2", "SGLT2i", "GLP1-RA"),
         drugclass = factor(drugclass, levels = c("SGLT2i", "GLP1-RA")),
         ncurrtx = factor(ncurrtx, levels = c("1", "2", "3", "4", "5+"), labels = c("1", "2", "3", "4", "4")))



###############################################################################
###############################################################################
############################ Check variables ##################################
###############################################################################
###############################################################################

###    Yr of prescription
#:----------------------------------------------

# table_yr_prescriptions <- table(dataset$yrdrugstart, dataset$drugsubstances, useNA = "ifany")

plot_year_prescription <- dataset %>%
  ggplot() +
  geom_bar(aes(x = yrdrugstart, fill = drugsubstances), position = "dodge")


dataset_formatted <- dataset %>%
  filter(yrdrugstart > 2018)


## Table of characteristics
#:----------------------------------------------

###    Categorical - ncurrtx, drugline, sex, smoke, prehospitalisation, preangina, precld, 
###         prediabeticnephropathy, preheartfailure, prehypertension, preihd, premyocardialinfarction, 
###         preneuropathy, prepad, preretinopathy, prerevasc, prestroke, pretia, preaf, preckd

###    Continuous - agetx, t2dmduration, prehba1c, prebmi, preegfr, preacr, prealbuminblood, prealt,
###         preast, prebilirubin, prefastingglucose, prehaematocrit, prehaemoglobin, prehdl, premap,
###         pretotalcholesterol, pretriglyceride, preweight

table_characteristics <- CreateTableOne(data = dataset_formatted,
                                        strata = "drugsubstances",
                                        vars = c("ncurrtx", "drugline", "sex", "ethnicity", "agetx",
                                                 "t2dmduration", "prehba1c", "prebmi", "preegfr", "preweight", "posthba1cfinal",
                                                 "prehypertension", "preihd", "premyocardialinfarction", "prehospitalisation", 
                                                 "preneuropathy", "prepad", "preretinopathy", "prerevasc", 
                                                 "prestroke", "pretia", "preaf", "preckd", "deprivation", "smoke",
                                                 "preacr", "preangina", "precld", "prediabeticnephropathy", "preheartfailure",
                                                 "prealbuminblood", "prealt", "preast", "prebilirubin", "prefastingglucose", 
                                                 "prehaematocrit", "prehaemoglobin", "prehdl", "premap", "pretotalcholesterol", 
                                                 "pretriglyceride"))

table_characteristics_print <- print(table_characteristics, nonnormal = c("agetx", "t2dmduration", "prehba1c", "prebmi", "preegfr", "preacr", 
                                                                          "prealbuminblood", "prealt", "preast", "prebilirubin", "prefastingglucose", 
                                                                          "prehaematocrit", "prehaemoglobin", "prehdl", "premap", "pretotalcholesterol", 
                                                                          "pretriglyceride", "preweight", "posthba1cfinal"), formatOptions = list(big.mark = ","), test = FALSE,
                                     quote = FALSE, noSpaces = TRUE, printToggle = FALSE)

write.csv(table_characteristics_print, file = paste0(output_path, "/sema_vs_glp1_sglt2.csv"))



###############################################################################
###############################################################################
########################## Match individuals ##################################
###############################################################################
###############################################################################

### HbA1c outcomes
#:----------------------------------------------

dataset_formatted_hba1c <- dataset_formatted %>%
  filter(!is.na(posthba1cfinal)) %>%
  filter(!is.na(prehba1c))


h.pscores <- SumStat(ps.formula = drugsubstances ~ agetx + sex + t2dmduration + drugline + ncurrtx + prehba1c,
                     data = dataset_formatted_hba1c,
                     weight = c("IPW", "overlap"))


lm.hba1c <- lm(formula = posthba1cfinal ~ drugsubstances + prehba1c + hba1cmonth,
               data = dataset_formatted_hba1c,
               weights = h.pscores$pw.weights$overlap)

patient.prediction <- data.frame(
  prehba1c = 75,
  drugsubstances = unique(dataset_formatted_hba1c$drugsubstances),
  hba1cmonth = 6
)

hba1c_estimates <- as.matrix(predict(lm.hba1c, patient.prediction, interval = "confidence") - 75) %>%
  as.data.frame() %>%
  cbind(drugsubstances = patient.prediction$drugsubstances) %>%
  left_join(
    dataset_formatted_hba1c %>%
      select(drugsubstances, drugclass) %>%
      group_by(drugsubstances) %>%
      mutate(n.patients = n()) %>%
      ungroup() %>%
      unique(),
    by = "drugsubstances"
  ) %>%
  mutate(drugclass = factor(drugclass, levels = c("SGLT2i", "GLP1-RA")))
  
  
plot_hba1c_estimates <- hba1c_estimates %>%
  ggplot(aes(x = drugsubstances, y = fit)) +
  geom_point() +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.5) +
  facet_grid(~drugclass, scales = "free_x") +
  theme_bw() +
  labs(
    y = "Adjusted Response (mmol/mol)"
  ) +
  theme(axis.title.x = element_blank())
  

lm.hba1c.sex <- lm(formula = posthba1cfinal ~ drugsubstances + prehba1c + hba1cmonth + sex + drugsubstances*sex,
                   data = dataset_formatted_hba1c,
                   weights = h.pscores$pw.weights$overlap)

patient.prediction <- expand.grid(
  prehba1c = 75,
  drugsubstances = unique(dataset_formatted_hba1c$drugsubstances),
  hba1cmonth = 6,
  sex = c("Male", "Female")
)
  

hba1c_estimates_sex <- as.matrix(predict(lm.hba1c.sex, patient.prediction, interval = "confidence") - 75) %>%
  as.data.frame() %>%
  cbind(drugsubstances = patient.prediction$drugsubstances,
        sex = patient.prediction$sex) %>%
  left_join(
    dataset_formatted_hba1c %>%
      select(drugsubstances, drugclass, sex) %>%
      group_by(drugsubstances, sex) %>%
      mutate(n.patients = n()) %>%
      ungroup() %>%
      unique(),
    by = c("drugsubstances", "sex")
  ) %>%
  mutate(drugclass = factor(drugclass, levels = c("SGLT2i", "GLP1-RA")))


plot_hba1c_estimates_sex <- hba1c_estimates_sex %>%
  ggplot(aes(x = drugsubstances, y = fit, colour = sex)) +
  geom_point(position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.5, position=position_dodge(width=0.5)) +
  facet_grid(~drugclass, scales = "free_x") +
  theme_bw() +
  labs(
    y = "Adjusted Response (mmol/mol)",
    colour = "Sex"
  ) +
  theme(axis.title.x = element_blank(),
        legend.position = "bottom")


###    Weight outcomes
#:----------------------------------------------

dataset_formatted_weight <- dataset_formatted %>%
  filter(!is.na(postweight6m)) %>%
  filter(!is.na(preweight))


w.pscores <- SumStat(ps.formula = drugsubstances ~ agetx + sex + t2dmduration + drugline + ncurrtx + preweight,
                     data = dataset_formatted_weight,
                     weight = c("IPW", "overlap"))


lm.weight <- lm(formula = postweight6m ~ drugsubstances + preweight,
               data = dataset_formatted_weight,
               weights = w.pscores$pw.weights$overlap)


patient.prediction <- data.frame(
  preweight = 95,
  drugsubstances = unique(dataset_formatted_weight$drugsubstances)
)


weight_estimates <- as.matrix(predict(lm.weight, patient.prediction, interval = "confidence")-95) %>%
  as.data.frame() %>%
  cbind(drugsubstances = patient.prediction$drugsubstances) %>%
  left_join(
    dataset_formatted_weight %>%
      select(drugsubstances, drugclass) %>%
      group_by(drugsubstances) %>%
      mutate(n.patients = n()) %>%
      ungroup() %>%
      unique(),
    by = "drugsubstances"
  ) %>%
  mutate(drugclass = factor(drugclass, levels = c("SGLT2i", "GLP1-RA")))


plot_weight_estimates <- weight_estimates %>%
  ggplot(aes(x = drugsubstances, y = fit)) +
  geom_point() +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.5) +
  facet_grid(~drugclass, scales = "free_x") +
  theme_bw() +
  labs(
    y = "Adjusted Response (kg)"
  ) +
  theme(axis.title.x = element_blank())




lm.weight.sex <- lm(formula = postweight6m ~ drugsubstances + preweight + drugsubstances*sex,
                   data = dataset_formatted_weight,
                   weights = h.pscores$pw.weights$overlap)


patient.prediction <- expand.grid(
  preweight = 95,
  drugsubstances = unique(dataset_formatted_weight$drugsubstances),
  sex = c("Male", "Female")
)


weight_estimates_sex <- as.matrix(predict(lm.weight.sex, patient.prediction, interval = "confidence")-95) %>%
  as.data.frame() %>%
  cbind(drugsubstances = patient.prediction$drugsubstances,
        sex = patient.prediction$sex) %>%
  left_join(
    dataset_formatted_weight %>%
      select(drugsubstances, drugclass, sex) %>%
      group_by(drugsubstances, sex) %>%
      mutate(n.patients = n()) %>%
      ungroup() %>%
      unique(),
    by = c("drugsubstances", "sex")
  ) %>%
  mutate(drugclass = factor(drugclass, levels = c("SGLT2i", "GLP1-RA")))


plot_weight_estimates_sex <- weight_estimates_sex %>%
  ggplot(aes(x = drugsubstances, y = fit, colour = sex)) +
  geom_point(position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.5, position=position_dodge(width=0.5)) +
  facet_grid(~drugclass, scales = "free_x") +
  theme_bw() +
  labs(
    y = "Adjusted Response (kg)",
    colour = "Sex"
  ) +
  theme(axis.title.x = element_blank(),
        legend.position = "bottom")




# ###    eGFR outcomes
# #:----------------------------------------------
# 
# dataset_formatted_egfr <- dataset_formatted %>%
#   filter(!is.na(preegfr)) %>%
#   filter(!is.na(postegfr6m))
# 
# 
# e.pscores <- SumStat(ps.formula = drugsubstances ~ agetx + sex + t2dmduration + drugline + ncurrtx + preegfr,
#                      data = dataset_formatted_egfr,
#                      weight = c("IPW", "overlap"))
# 
# 
# lm.egfr <- lm(formula = postegfr6m ~ drugsubstances + preegfr,
#                 data = dataset_formatted_egfr,
#                 weights = w.pscores$pw.weights$overlap)
# 
# 
# patient.prediction <- data.frame(
#   preegfr = 95,
#   drugsubstances = unique(dataset_formatted_egfr$drugsubstances)
# )
# 
# 
# egfr_estimates <- as.matrix(predict(lm.egfr, patient.prediction, interval = "confidence")-95) %>%
#   as.data.frame() %>%
#   cbind(drugsubstances = patient.prediction$drugsubstances) %>%
#   left_join(
#     dataset_formatted_egfr %>%
#       select(drugsubstances, drugclass) %>%
#       group_by(drugsubstances) %>%
#       mutate(n.patients = n()) %>%
#       ungroup() %>%
#       unique(),
#     by = "drugsubstances"
#   ) %>%
#   mutate(drugclass = factor(drugclass, levels = c("SGLT2i", "GLP1-RA")))
# 
# 
# plot_egfr_estimates <- egfr_estimates %>%
#   ggplot(aes(x = drugsubstances, y = fit)) +
#   geom_point() +
#   geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.5) +
#   facet_grid(~drugclass, scales = "free_x") +
#   theme_bw() +
#   labs(
#     y = "Adjusted Response (egfr units)"
#   ) +
#   theme(axis.title.x = element_blank())
# 
# 
# 
# lm.egfr.sex <- lm(formula = postegfr6m ~ drugsubstances + preegfr + drugsubstances*sex,
#                     data = dataset_formatted_egfr,
#                     weights = h.pscores$pw.weights$overlap)
# 
# 
# patient.prediction <- expand.grid(
#   preegfr = 95,
#   drugsubstances = unique(dataset_formatted_egfr$drugsubstances),
#   sex = c("Male", "Female")
# )
# 
# 
# egfr_estimates_sex <- as.matrix(predict(lm.egfr.sex, patient.prediction, interval = "confidence")-95) %>%
#   as.data.frame() %>%
#   cbind(drugsubstances = patient.prediction$drugsubstances,
#         sex = patient.prediction$sex) %>%
#   left_join(
#     dataset_formatted_egfr %>%
#       select(drugsubstances, drugclass, sex) %>%
#       group_by(drugsubstances, sex) %>%
#       mutate(n.patients = n()) %>%
#       ungroup() %>%
#       unique(),
#     by = c("drugsubstances", "sex")
#   ) %>%
#   mutate(drugclass = factor(drugclass, levels = c("SGLT2i", "GLP1-RA")))
# 
# 
# plot_egfr_estimates_sex <- egfr_estimates_sex %>%
#   ggplot(aes(x = drugsubstances, y = fit, colour = sex)) +
#   geom_point(position=position_dodge(width=0.5)) +
#   geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.5, position=position_dodge(width=0.5)) +
#   facet_grid(~drugclass, scales = "free_x") +
#   theme_bw() +
#   labs(
#     y = "Adjusted Response (egfr units)",
#     colour = "Sex"
#   ) +
#   theme(axis.title.x = element_blank(),
#         legend.position = "bottom")





###    Discontinuation outcomes
#:----------------------------------------------

dataset_formatted_discontinuation <- dataset_formatted %>%
  filter(!is.na(stopdrug_6m_3mFU))


d.pscores <- SumStat(ps.formula = drugsubstances ~ agetx + sex + t2dmduration + drugline + ncurrtx,
                     data = dataset_formatted_discontinuation,
                     weight = c("IPW", "overlap"))


lm.discontinuation <- lm(formula = stopdrug_6m_3mFU ~ drugsubstances,
              data = dataset_formatted_discontinuation,
              weights = w.pscores$pw.weights$overlap)


patient.prediction <- data.frame(
  drugsubstances = unique(dataset_formatted_discontinuation$drugsubstances)
)


discontinuation_estimates <- predict(lm.discontinuation, patient.prediction, interval = "confidence") %>%
  as.data.frame() %>%
  cbind(drugsubstances = patient.prediction$drugsubstances) %>%
  left_join(
    dataset_formatted_discontinuation %>%
      select(drugsubstances, drugclass) %>%
      group_by(drugsubstances) %>%
      mutate(n.patients = n()) %>%
      ungroup() %>%
      unique(),
    by = "drugsubstances"
  ) %>%
  mutate(drugclass = factor(drugclass, levels = c("SGLT2i", "GLP1-RA")))


plot_discontinuation_estimates <- discontinuation_estimates %>%
  ggplot(aes(x = drugsubstances, y = fit)) +
  geom_point() +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.5) +
  facet_grid(~drugclass, scales = "free_x") +
  theme_bw() +
  scale_y_continuous(labels=percent) +
  labs(
    y = "Adjusted Response (%)"
  ) +
  theme(axis.title.x = element_blank())



lm.discontinuation.sex <- lm(formula = stopdrug_6m_3mFU ~ drugsubstances + drugsubstances*sex,
                  data = dataset_formatted_discontinuation,
                  weights = h.pscores$pw.weights$overlap)


patient.prediction <- expand.grid(
  drugsubstances = unique(dataset_formatted_discontinuation$drugsubstances),
  sex = c("Male", "Female")
)


discontinuation_estimates_sex <- predict(lm.discontinuation.sex, patient.prediction, interval = "confidence") %>%
  as.data.frame() %>%
  cbind(drugsubstances = patient.prediction$drugsubstances,
        sex = patient.prediction$sex) %>%
  left_join(
    dataset_formatted_discontinuation %>%
      select(drugsubstances, drugclass, sex) %>%
      group_by(drugsubstances, sex) %>%
      mutate(n.patients = n()) %>%
      ungroup() %>%
      unique(),
    by = c("drugsubstances", "sex")
  ) %>%
  mutate(drugclass = factor(drugclass, levels = c("SGLT2i", "GLP1-RA")))


plot_discontinuation_estimates_sex <- discontinuation_estimates_sex %>%
  ggplot(aes(x = drugsubstances, y = fit, colour = sex)) +
  geom_point(position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.5, position=position_dodge(width=0.5)) +
  facet_grid(~drugclass, scales = "free_x") +
  theme_bw() +
  scale_y_continuous(labels=percent) +
  labs(
    y = "Adjusted Response (%)",
    colour = "Sex"
  ) +
  theme(axis.title.x = element_blank(),
        legend.position = "bottom")




###############################################################################
###############################################################################
########################## Combine Results ####################################
###############################################################################
###############################################################################

plot_estimates <- patchwork::wrap_plots(
  list(
    plot_hba1c_estimates +
      ggtitle("Average HbA1c response (baseline 75)"), 
    plot_weight_estimates +
      ggtitle("Average weight response (baseline 95)"), 
    plot_discontinuation_estimates +
      ggtitle("Average discontinuation")
  )) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

plot_estimates_sex <- patchwork::wrap_plots(
  list(
    plot_hba1c_estimates_sex +
      ggtitle("Average HbA1c response (baseline 75)"), 
    plot_weight_estimates_sex +
      ggtitle("Average weight response (baseline 95)"), 
    plot_discontinuation_estimates_sex +
      ggtitle("Average discontinuation")
  )) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")




pdf(height = 6, width = 19, paste0(output_path, "/comparison_post_2018.pdf"))
plot_estimates
plot_estimates_sex
dev.off()






# just hba1c

patchwork::wrap_plots(
  list(
    plot_hba1c_estimates +
      ggtitle("Average HbA1c response (baseline 75)") +
      ylim(-20, -4), 
    plot_hba1c_estimates_sex +
      ggtitle("Average HbA1c response (baseline 75)") +
      ylim(-20, -4)
  )) &
  theme(legend.position = "bottom")

# just weight

patchwork::wrap_plots(
  list(
    plot_weight_estimates +
      ggtitle("Average weight response (baseline 95)") +
      ylim(-7, -1), 
    plot_weight_estimates_sex +
      ggtitle("Average weight response (baseline 95)") +
      ylim(-7, -1)
  )) &
  theme(legend.position = "bottom")


# just discontinuation

patchwork::wrap_plots(
  list(
    plot_discontinuation_estimates +
      ggtitle("Average weight response (baseline 95)") +
      scale_y_continuous(labels=percent, limits = c(0.13, 0.28)), 
    plot_discontinuation_estimates_sex +
      ggtitle("Average weight response (baseline 95)") +
      scale_y_continuous(labels=percent, limits = c(0.13, 0.28))
  )) &
  theme(legend.position = "bottom")








