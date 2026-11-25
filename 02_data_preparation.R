# =============================================================================
# 02_data_preparation_H1N1_dplyr.R
# VEBIS training project - influenza A(H1N1)pdm09 (synthetic data)
# ==============================================================================

library(dplyr)

# 1. Import data ---------------------------------------------------------------
patients <- readRDS("data/synthetic_vebis_patients.rds")
calendar <- read.csv("data/synthetic_country_calendar.csv") %>%
  mutate(flu_campaign_start = as.Date(flu_campaign_start))

dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)

# 2. Helper: filter patients and record exclusions -----------------------------
filter_step <- function(sample, step, ...) {
  n_before <- nrow(sample$data)
  sample$data <- filter(sample$data, ...)
  sample$exclusions <- bind_rows(sample$exclusions, tibble(
    step = step,
    excluded = n_before - nrow(sample$data),
    remaining = nrow(sample$data)
  ))
  sample
}

# 3. Define SARI-CD1 and SARI-CD2 ----------------------------------------------
patients <- patients %>%
  mutate(
    SARI_CD1 = hospital_stay_hours >= 24 &
      (fever == 1 |
         cough == 1 |
         dyspnea == 1 |
         anosmia == 1 |
         ageusia == 1 |
         dysgeusia == 1),

    SARI_CD2 = hospital_stay_hours >= 24 &
      fever == 1 &
      cough == 1
  )

# check that SARI-CD2 is a subset of SARI-CD1
stopifnot(all(patients$SARI_CD2 <= patients$SARI_CD1))


# 4. Apply general influenza eligibility criteria -----------------------------
flu <- list(data = patients, exclusions = tibble())

flu <- filter_step(flu, "SARI_CD1", flu$data$SARI_CD1)
flu <- filter_step(flu, "Influenza test available", !is.na(flu$data$flu_result))

chronic <- c("diabetes", "heart_disease", "lung_disease",
             "asthma", "immunodeficiency")
required <- c(chronic, "age", "sex", "onset_date",
              "admission_date", "swab_date", "flu_vaccinated_2023_24")

flu <- filter_step(
  flu, "Missing key variables",
  if_all(all_of(required), ~ !is.na(.x)),
  flu$data$flu_vaccinated_2023_24 == 0 | !is.na(flu$data$flu_vaccination_date)
)

flu <- filter_step(
  flu, "Swab within 7 days of onset",
  between(as.numeric(flu$data$swab_date - flu$data$onset_date), 0, 7)
)

flu <- filter_step(
  flu, "Admission within 14 days of onset",
  as.numeric(admission_date - onset_date) <= 14
)

flu <- filter_step(
  flu, "Eligible population",
  age >= 65, long_term_care == 0,
  healthcare_worker == 0, vaccine_contraindication == 0
)

flu$data <- flu$data %>%
  mutate(
    any_chronic = as.integer(if_any(all_of(chronic), ~ .x == 1)),
    days_since_vaccination = as.numeric(onset_date - flu_vaccination_date)
  )

flu <- filter_step(
  flu, "Vaccination 1-13 days before onset",
  flu_vaccinated_2023_24 == 0 |
    days_since_vaccination <= 0 |
    days_since_vaccination >= 14
)

flu$data <- flu$data %>%
  mutate(
    vaccinated = as.integer(flu_vaccinated_2023_24 == 1 &
                              !is.na(days_since_vaccination) &
                              days_since_vaccination >= 14),
    age_group = if_else(age < 80, "65-79", "80+")
  )

# 5. H1N1 cases and controls --------------------------------------------------
# Keep H1N1-positive cases and influenza-negative controls
h1n1 <- flu

h1n1 <- filter_step(
  h1n1,
  "H1N1 cases and controls",
  h1n1$data$flu_result %in% c("A(H1N1)pdm09", "Negative")
)

# Create case status: 1 = case, 0 = control
h1n1$data <- h1n1$data %>%
  mutate(
    case = as.integer(flu_result == "A(H1N1)pdm09")
  )

# Restrict to the end of the study period
h1n1 <- filter_step(
  h1n1,
  "Study period",
  h1n1$data$onset_date <= as.Date("2024-05-19")
)

# Country-specific analysis periods were not applied; data are synthetic




# 6. CD1 and CD2: keep countries with >=5 cases and controls -------------------
valid_countries <- function(df) {
  df %>%
    group_by(country) %>%
    summarise(
      cases = sum(case == 1),
      controls = sum(case == 0),
      .groups = "drop"
    ) %>%
    filter(cases >= 5, controls >= 5) %>%
    pull(country)
}

h1n1_cd1 <- filter_step(
  h1n1, "Countries with >=5 cases and controls",
  country %in% valid_countries(h1n1$data)
)

h1n1_cd2 <- filter_step(h1n1, "SARI-CD2", SARI_CD2)
h1n1_cd2 <- filter_step(
  h1n1_cd2, "Countries with >=5 cases and controls",
  country %in% valid_countries(h1n1_cd2$data)
)

# 7. Export --------------------------------------------------------------------
saveRDS(h1n1_cd1$data, "data/processed/H1N1_CD1.rds")
saveRDS(h1n1_cd2$data, "data/processed/H1N1_CD2.rds")

exclusion_summary <- bind_rows(
  mutate(h1n1_cd1$flow, definition = "CD1"),
  mutate(h1n1_cd2$flow, definition = "CD2")
)
write.csv(exclusion_summary, "data/processed/H1N1_exclusion_summary.csv",
          row.names = FALSE)

cat("H1N1 CD1:", nrow(h1n1_cd1$data), "patients\n")
cat("H1N1 CD2:", nrow(h1n1_cd2$data), "patients\n")
print(table(h1n1_cd1$data$case))
print(table(h1n1_cd2$data$case))
