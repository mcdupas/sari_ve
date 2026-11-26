# =============================================================================
# 03_vaccination_effectiveness_analysis.R
# VEBIS training project - influenza A(H1N1)pdm09 (synthetic data)
# =============================================================================

# vaccine effectiveness analysis for influenza A(H1N1)pdm09 using synthetic data
# VE = (1 - aOR) * 100, where aOR is the adjusted odds ratio
#                       of vaccination among cases vs controls
# aOR: cofounder-adjusted odds ratio of vaccination among cases vs controls
#      from logistic regression model.
#      adjuster for age, sex, at least one condition, onset date and country


# 1. Load required packages ----------------------------------------------------
library(dplyr)
library(broom)
library(ggplot2)

# 2. Load data -----------------------------------------------------------------
h1n1_cd1 <- readRDS("data/processed/H1N1_CD1.rds")
h1n1_cd2 <- readRDS("data/processed/H1N1_CD2.rds")


# 3. Prepare data for logistic regression -------------------------------------

# replace "Negative" with 0 and "Positive" with 1 in the flu_result column
h1n1_cd1 <- h1n1_cd1 %>%
  mutate(test_positive = ifelse(flu_result == "Negative", 0, 1))

h1n1_cd2 <- h1n1_cd2 %>%
  mutate(test_positive = ifelse(flu_result == "Negative", 0, 1))

# add a column comorbidity if the patient has at least one condition:
# "diabetes", "heart_disease", "lung_disease",
# "asthma", "immunodeficiency"
h1n1_cd1 <- h1n1_cd1 %>%
  mutate(comorbidity = ifelse(diabetes == 1 |
                                heart_disease == 1 |
                                lung_disease == 1 |
                                asthma == 1 |
                                immunodeficiency == 1, 1, 0))

h1n1_cd2 <- h1n1_cd2 %>%
  mutate(comorbidity = ifelse(diabetes == 1 |
                                heart_disease == 1 |
                                lung_disease == 1 |
                                asthma == 1 |
                                immunodeficiency == 1, 1, 0))

# Verify the contingency table
table(h1n1_cd1$flu_vaccinated_2023_24, h1n1_cd1$test_positive)


# 4. Logistic regression analysis ----------------------------------------------
# helper function to calculate vaccine effectiveness and
# 95% CI from logistic regression model
vaccine_effectiveness <- function(model, variable) {
  beta <- coef(model)[variable]
  SE <- summary(model)$coefficients[variable, "Std. Error"]

  # Odds ratio and 95% CI
  OR <- exp(beta)
  L <- exp(beta - 1.96 * SE)
  U <- exp(beta + 1.96 * SE)

  # Vaccine effectiveness and 95% CI
  VE <- (1 - OR) * 100
  VE_lower <- (1 - U) * 100
  VE_upper <- (1 - L) * 100

  list(VE = VE, VE_lower = VE_lower, VE_upper = VE_upper)
}

# Adjusted logistic regression
model_non_adj_cd1 <- glm(
  test_positive ~ flu_vaccinated_2023_24,
  data = h1n1_cd1,
  family = binomial(link = "logit")
)

summary(model_non_adj_cd1)


# Regression coefficient and standard error
ve <- vaccine_effectiveness(model_non_adj_cd1, "flu_vaccinated_2023_24")
# Display results
sprintf(
  "VE = %.1f%% (95%% CI: %.1f%% to %.1f%%)",
  ve$VE, ve$VE_lower, ve$VE_upper
)

######### Adjusted logistic regression with confounders
model_adj_cd1 <- glm(test_positive ~ flu_vaccinated_2023_24 +
                       age + sex + comorbidity + onset_date + country,
                     data = h1n1_cd1, family = binomial(link = "logit"))

summary(model_adj_cd1)

# Vaccine effectiveness and 95% CI
ve_adj <- vaccine_effectiveness(model_adj_cd1, "flu_vaccinated_2023_24")

# Display results
sprintf(
  "Adjusted VE = %.1f%% (95%% CI: %.1f%% to %.1f%%)",
  ve_adj$VE, ve_adj$VE_lower, ve_adj$VE_upper
)
