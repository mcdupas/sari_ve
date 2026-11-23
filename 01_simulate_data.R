# ==============================================================================
# 01_simulate_data.R
# Projet d'entrainement R - VEBIS / test-negative design (TND)
# ==============================================================================
# OBJECTIF
#   Creer une cohorte SYNTHETIQUE de patients hospitalises pour une affection
#   respiratoire pendant la saison 2023/24. Les donnees serviront a pratiquer :
#   - la reconstruction des definitions SARI-CD1 et SARI-CD2 ;
#   - l'application de criteres d'eligibilite ;
#   - les analyses d'efficacite vaccinale influenza / SARS-CoV-2 en TND.
#
# SOURCE METHODOLOGIQUE (PAS SOURCE DES DONNEES)
#   Marques et al., Vaccine (2026), doi:10.1016/j.vaccine.2026.129234,
#   article et supplement (notamment figure 1, tableaux S1-S7).
#
# IMPORTANT
#   100 % des observations sont fictives. Les probabilites, effets des vaccins,
#   dates de campagnes, dates de predominance JN.1 et effets pays ci-dessous
#   sont des choix pedagogiques et NE SONT PAS des estimations VEBIS.
#   Les effectifs et VE publies ne sont pas vises et ne doivent pas etre cites
#   comme des resultats retrouves avec ces donnees.
#
# PREREQUIS : R >= 4.0, aucun package externe.
# EXECUTION  : Rscript 01_simulate_data.R
# SORTIES    : data/synthetic_vebis_patients.csv
#              data/synthetic_vebis_patients.rds
#              data/synthetic_country_calendar.csv
# ===============================================================================

set.seed(20261123)
n_patients <- 16000L
output_dir <- "data"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- Fonctions utilitaires ----------------------------------------------------
bernoulli <- function(p) {
  # Les probabilites sont bornees pour eviter les valeurs numeriques invalides.
  p <- pmin(1, pmax(0, p))
  as.integer(stats::runif(length(p)) < p)
}

inject_na <- function(x, proportion) {
  if (!is.numeric(proportion) || length(proportion) != 1L ||
      proportion < 0 || proportion > 1) {
    stop("proportion doit etre comprise entre 0 et 1")
  }
  if (length(x) == 0L || proportion == 0) return(x)
  positions <- sample.int(length(x), size = round(length(x) * proportion))
  x[positions] <- NA
  x
}

# --- 1. Calendriers et parametres pays : valeurs FICTIVES ----------------------
# Les huit pays correspondent au schema initial de la figure 1 du papier.
country_calendar <- data.frame(
  country = c("CZ", "ES", "HR", "HU", "IE", "LT", "PT", "RO"),
  flu_campaign_start = as.Date(c(
    "2023-10-05", "2023-10-02", "2023-10-09", "2023-10-02",
    "2023-10-03", "2023-10-09", "2023-10-02", "2023-10-10"
  )),
  covid_campaign_start = as.Date(c(
    "2023-10-02", "2023-09-25", "2023-10-02", "2023-10-02",
    "2023-09-25", "2023-10-02", "2023-09-25", "2023-10-09"
  )),
  jn1_predominance_start = as.Date(c(
    "2023-12-11", "2023-12-04", "2023-12-11", "2023-12-11",
    "2023-12-04", "2023-12-18", "2023-12-04", "2023-12-18"
  )),
  # Seuils d'age simplifiés pour l'exercice, inspires du tableau S1 ;
  # le script 02 pourra les utiliser pour generer ses propres indicateurs.
  flu_priority_age = c(50L, 60L, 65L, 60L, 65L, 65L, 60L, 65L),
  covid_priority_age = c(50L, 60L, 65L, 60L, 50L, 65L, 60L, 18L),
  stringsAsFactors = FALSE
)

country_weight <- c(CZ = 0.05, ES = 0.29, HR = 0.06, HU = 0.15,
                    IE = 0.13, LT = 0.11, PT = 0.14, RO = 0.07)
country_effect <- c(CZ = -0.15, ES = 0.15, HR = -0.05, HU = 0.02,
                    IE = 0.18, LT = -0.10, PT = 0.12, RO = -0.22)
country <- sample(names(country_weight), size = n_patients,
                  replace = TRUE, prob = country_weight)
country_index <- match(country, country_calendar$country)

# --- 2. Demographie et comorbidites -------------------------------------------
patient_id <- sprintf("SYN%06d", seq_len(n_patients))

age_group_generated <- sample(
  c("40-64", "65-79", "80-100"), n_patients, replace = TRUE,
  prob = c(0.25, 0.34, 0.41)
)
age <- integer(n_patients)
age[age_group_generated == "40-64"] <- sample(
  40:64, sum(age_group_generated == "40-64"), replace = TRUE)
age[age_group_generated == "65-79"] <- sample(
  65:79, sum(age_group_generated == "65-79"), replace = TRUE)
age[age_group_generated == "80-100"] <- sample(
  80:100, sum(age_group_generated == "80-100"), replace = TRUE)

sex <- sample(c("Female", "Male"), n_patients, replace = TRUE,
              prob = c(0.49, 0.51))

diabetes <- bernoulli(plogis(-1.10 + 0.018 * (age - 65)))
heart_disease <- bernoulli(plogis(-0.45 + 0.047 * (age - 65)))
lung_disease <- bernoulli(plogis(-1.25 + 0.017 * (age - 65)))
asthma <- bernoulli(plogis(-2.05 + 0.005 * (age - 65)))
immunodeficiency <- bernoulli(plogis(-2.55 + 0.010 * (age - 65)))

# Indicateur interne, a recalculer dans le script 02 avec les colonnes brutes.
any_chronic <- as.integer(diabetes == 1L | heart_disease == 1L |
                          lung_disease == 1L | asthma == 1L |
                          immunodeficiency == 1L)

long_term_care <- bernoulli(plogis(-4.0 + 0.045 * (age - 65)))
healthcare_worker <- bernoulli(ifelse(age < 65, 0.035, 0.004))
vaccine_contraindication <- bernoulli(rep(0.002, n_patients))

# --- 3. Dates de debut des symptomes et hospitalisation ------------------------
# Pic hivernal fictif, sans tentative de reproduire exactement la figure 2.
possible_onsets <- seq(as.Date("2023-10-01"), as.Date("2024-05-19"), by = "day")
relative_days <- as.numeric(possible_onsets - as.Date("2023-10-01"))
onset_weights <- (0.35 +
                  2.8 * exp(-0.5 * ((relative_days - 87) / 32)^2) +
                  1.0 * exp(-0.5 * ((relative_days - 145) / 38)^2))
onset_date <- sample(possible_onsets, size = n_patients,
                     replace = TRUE, prob = onset_weights)

# La plupart des hospitalisations surviennent rapidement ; quelques observations
# ont des dates incompatibles ou tardives pour exercer les filtres d'exclusion.
onset_to_admission_days <- sample(
  -2:19, n_patients, replace = TRUE,
  prob = c(1, 1, 20, 21, 18, 13, 10, 8, 6, 5, 4, 3,
           2, 2, 2, 2, 1, 1, 1, 1, 1, 1)
)
admission_date <- onset_date + onset_to_admission_days

# Le prelevement a lieu le plus souvent autour de l'admission.
swab_offset_from_admission <- sample(
  -3:4, n_patients, replace = TRUE,
  prob = c(1, 4, 12, 40, 23, 11, 6, 3)
)
swab_date <- admission_date + swab_offset_from_admission

# Quelques hospitalisations <24 h et des issues cliniques simulées.
hospital_stay_hours <- pmax(6L, round(24 * (1 + stats::rpois(n_patients, 3))))
short_stay <- bernoulli(rep(0.04, n_patients)) == 1L
hospital_stay_hours[short_stay] <- sample(6:23, sum(short_stay), replace = TRUE)

# --- 4. Vaccination : comportements correles influenza / COVID -----------------
# La propension latente cree une correlation de comportement vaccinal.
healthcare_propensity <- stats::rnorm(n_patients)

p_flu_uptake <- plogis(
  0.05 + 0.029 * (age - 65) + 0.30 * any_chronic +
    0.58 * healthcare_propensity + unname(country_effect[country])
)
p_covid_uptake <- plogis(
  -0.18 + 0.026 * (age - 65) + 0.34 * any_chronic +
    0.67 * healthcare_propensity + unname(country_effect[country])
)

will_get_flu_vaccine <- bernoulli(p_flu_uptake) == 1L
will_get_covid_vaccine <- bernoulli(p_covid_uptake) == 1L

flu_campaign_for_patient <- country_calendar$flu_campaign_start[country_index]
covid_campaign_for_patient <- country_calendar$covid_campaign_start[country_index]

flu_planned_date <- flu_campaign_for_patient +
  round(pmax(0, pmin(135, stats::rnorm(n_patients, mean = 27, sd = 24))))
covid_planned_date <- covid_campaign_for_patient +
  round(pmax(0, pmin(145, stats::rnorm(n_patients, mean = 24, sd = 24))))

# Aucune date de vaccination future apres l'hospitalisation n'est inscrite.
flu_vaccination_date <- flu_planned_date
flu_vaccination_date[!will_get_flu_vaccine |
                       flu_planned_date > admission_date] <- as.Date(NA_character_)
covid_vaccination_date <- covid_planned_date
covid_vaccination_date[!will_get_covid_vaccine |
                         covid_planned_date > admission_date] <- as.Date(NA_character_)

flu_vaccinated_2023_24 <- as.integer(!is.na(flu_vaccination_date))
covid_vaccinated_2023_24 <- as.integer(!is.na(covid_vaccination_date))

# Variables internes non exportees : au moins 14 jours AVANT le debut des symptomes.
flu_days_since_vax <- as.numeric(onset_date - flu_vaccination_date)
covid_days_since_vax <- as.numeric(onset_date - covid_vaccination_date)
flu_protected <- !is.na(flu_days_since_vax) & flu_days_since_vax >= 14
covid_protected <- !is.na(covid_days_since_vax) & covid_days_since_vax >= 14

# Categorie de produit vaccinal (optionnelle pour analyses descriptives).
flu_vaccine_product <- rep(NA_character_, n_patients)
flu_vaccine_product[flu_vaccinated_2023_24 == 1L] <- sample(
  c("Standard", "High-dose", "Adjuvanted", "Cell-based"),
  sum(flu_vaccinated_2023_24 == 1L), replace = TRUE,
  prob = c(0.20, 0.31, 0.36, 0.13)
)
covid_vaccine_product <- rep(NA_character_, n_patients)
covid_vaccine_product[covid_vaccinated_2023_24 == 1L] <- sample(
  c("Pfizer", "Moderna", "Other"),
  sum(covid_vaccinated_2023_24 == 1L), replace = TRUE,
  prob = c(0.84, 0.01, 0.15)
)

# Deux indicateurs fictifs additionnels permettent de s'exercer aux exclusions
# specifiques COVID-19 presentees dans la figure 1.
covid_prior_dose_within_180d <- bernoulli(rep(0.015, n_patients))
covid_vaccine_authorised <- rep(1L, n_patients)
not_authorised <- covid_vaccinated_2023_24 == 1L &
  bernoulli(rep(0.006, n_patients)) == 1L
covid_vaccine_authorised[not_authorised] <- 0L

# --- 5. Etiologie virologique dependante du temps et de la vaccination ---------
# Hypothese generatrice : efficacite protectrice moderee, possiblement moindre
# a >=120 jours. Elle est entierement SYNTHETIQUE et n'est pas celle du papier.
day_number <- as.numeric(onset_date - as.Date("2023-10-01"))
wave_h1 <- exp(-0.5 * ((day_number - 86) / 36)^2)
wave_h3 <- exp(-0.5 * ((day_number - 112) / 34)^2)
wave_covid <- exp(-0.5 * ((day_number - 94) / 43)^2)

flu_risk_multiplier <- ifelse(
  flu_protected, ifelse(flu_days_since_vax >= 120, 0.79, 0.60), 1
)
covid_risk_multiplier <- ifelse(
  covid_protected, ifelse(covid_days_since_vax >= 120, 0.87, 0.66), 1
)

country_infection_factor <- exp(unname(country_effect[country]) * 0.20)
p_h1 <- (0.028 + 0.120 * wave_h1) * flu_risk_multiplier *
  country_infection_factor
p_h3 <- (0.024 + 0.075 * wave_h3) * flu_risk_multiplier *
  country_infection_factor
p_flu_b <- rep(0.009, n_patients)
p_flu_unsubtyped <- rep(0.020, n_patients)

# Tirage exclusif des sous-types influenza (une seule categorie par patient).
flu_draw <- stats::runif(n_patients)
flu_result <- rep("Negative", n_patients)
flu_result[flu_draw < p_h1] <- "A(H1N1)pdm09"
flu_result[flu_draw >= p_h1 & flu_draw < p_h1 + p_h3] <- "A(H3N2)"
flu_result[flu_draw >= p_h1 + p_h3 &
           flu_draw < p_h1 + p_h3 + p_flu_b] <- "B"
flu_result[flu_draw >= p_h1 + p_h3 + p_flu_b &
           flu_draw < p_h1 + p_h3 + p_flu_b + p_flu_unsubtyped] <- "A_unsubtyped"

p_covid <- pmin(0.90, (0.035 + 0.130 * wave_covid) *
                  covid_risk_multiplier * country_infection_factor)
covid_result <- ifelse(bernoulli(p_covid) == 1L, "Positive", "Negative")

# --- 6. Symptomes : dependance a l'infection ET a la vaccination --------------
# Conditionner les symptomes a l'infection permet d'explorer la selection
# potentielle produite par SARI-CD2. Un effet d'attenuation symptomatique
# chez les vaccines infectes est ajoute UNIQUEMENT pour la simulation.
influenza_positive <- flu_result != "Negative"
respiratory_positive <- influenza_positive | covid_result == "Positive"
attenuation <- as.numeric((influenza_positive & flu_protected) |
                          (covid_result == "Positive" & covid_protected))

p_fever <- plogis(-0.15 + 0.92 * influenza_positive +
                    0.55 * (covid_result == "Positive") - 0.45 * attenuation)
p_cough <- plogis(1.05 + 0.57 * influenza_positive +
                    0.38 * (covid_result == "Positive") - 0.20 * attenuation)
p_dyspnea <- plogis(-0.95 + 0.025 * (age - 65) +
                      0.25 * any_chronic + 0.37 * respiratory_positive)
p_anosmia <- plogis(-4.5 + 1.6 * (covid_result == "Positive"))
p_ageusia <- plogis(-4.7 + 1.5 * (covid_result == "Positive"))
p_dysgeusia <- plogis(-4.8 + 1.4 * (covid_result == "Positive"))

fever <- bernoulli(p_fever)
cough <- bernoulli(p_cough)
dyspnea <- bernoulli(p_dyspnea)
anosmia <- bernoulli(p_anosmia)
ageusia <- bernoulli(p_ageusia)
dysgeusia <- bernoulli(p_dysgeusia)

# Ne pas exporter SARI-CD1/CD2 : le script 02 sera l'occasion de les creer.
cd1_for_qc <- fever == 1L | cough == 1L | dyspnea == 1L |
  anosmia == 1L | ageusia == 1L | dysgeusia == 1L
cd2_for_qc <- fever == 1L & cough == 1L
stopifnot(all(!cd2_for_qc | cd1_for_qc))

# --- 7. Issues cliniques (descriptif, pas outcomes TND) -------------------------
p_icu <- plogis(-3.55 + 0.35 * respiratory_positive +
                  0.015 * (age - 65) + 0.15 * any_chronic)
icu_admission <- bernoulli(p_icu)
p_death <- plogis(-3.20 + 0.025 * (age - 65) +
                    0.50 * icu_admission + 0.13 * any_chronic)
in_hospital_death <- bernoulli(p_death)

# --- 8. Introduire un peu de donnees manquantes -------------------------------
# Appliquer les NA APRES la generation des relations causales permet de
# conserver des simulations coherentes tout en testant le nettoyage.
age <- inject_na(age, 0.003)
sex <- inject_na(sex, 0.003)
onset_date <- inject_na(onset_date, 0.008)
swab_date <- inject_na(swab_date, 0.006)
flu_result <- inject_na(flu_result, 0.045)
covid_result <- inject_na(covid_result, 0.075)
diabetes <- inject_na(diabetes, 0.015)
heart_disease <- inject_na(heart_disease, 0.012)
lung_disease <- inject_na(lung_disease, 0.012)
asthma <- inject_na(asthma, 0.010)
immunodeficiency <- inject_na(immunodeficiency, 0.010)
flu_vaccinated_2023_24 <- inject_na(flu_vaccinated_2023_24, 0.012)
covid_vaccinated_2023_24 <- inject_na(covid_vaccinated_2023_24, 0.014)

# Certaines personnes vaccines avec statut connu ont une date manquante.
flu_known_vaccinated <- which(!is.na(flu_vaccinated_2023_24) &
                               flu_vaccinated_2023_24 == 1L)
covid_known_vaccinated <- which(!is.na(covid_vaccinated_2023_24) &
                                 covid_vaccinated_2023_24 == 1L)
flu_vaccination_date[sample(flu_known_vaccinated,
  size = round(length(flu_known_vaccinated) * 0.012))] <- as.Date(NA_character_)
covid_vaccination_date[sample(covid_known_vaccinated,
  size = round(length(covid_known_vaccinated) * 0.015))] <- as.Date(NA_character_)

# --- 9. Export des donnees brutes synthetiques --------------------------------
synthetic_patients <- data.frame(
  patient_id = patient_id,
  country = country,
  age = age,
  sex = sex,
  diabetes = diabetes,
  heart_disease = heart_disease,
  lung_disease = lung_disease,
  asthma = asthma,
  immunodeficiency = immunodeficiency,
  long_term_care = long_term_care,
  healthcare_worker = healthcare_worker,
  vaccine_contraindication = vaccine_contraindication,
  onset_date = onset_date,
  admission_date = admission_date,
  swab_date = swab_date,
  hospital_stay_hours = hospital_stay_hours,
  fever = fever,
  cough = cough,
  dyspnea = dyspnea,
  anosmia = anosmia,
  ageusia = ageusia,
  dysgeusia = dysgeusia,
  flu_result = flu_result,
  covid_result = covid_result,
  flu_vaccinated_2023_24 = flu_vaccinated_2023_24,
  flu_vaccination_date = flu_vaccination_date,
  flu_vaccine_product = flu_vaccine_product,
  covid_vaccinated_2023_24 = covid_vaccinated_2023_24,
  covid_vaccination_date = covid_vaccination_date,
  covid_vaccine_product = covid_vaccine_product,
  covid_prior_dose_within_180d = covid_prior_dose_within_180d,
  covid_vaccine_authorised = covid_vaccine_authorised,
  icu_admission = icu_admission,
  in_hospital_death = in_hospital_death,
  stringsAsFactors = FALSE
)

stopifnot(nrow(synthetic_patients) == n_patients)
stopifnot(!anyDuplicated(synthetic_patients$patient_id))
stopifnot(all(is.na(synthetic_patients$flu_vaccination_date) |
              synthetic_patients$flu_vaccination_date <=
                synthetic_patients$admission_date))
stopifnot(all(is.na(synthetic_patients$covid_vaccination_date) |
              synthetic_patients$covid_vaccination_date <=
                synthetic_patients$admission_date))

utils::write.csv(synthetic_patients,
  file.path(output_dir, "synthetic_vebis_patients.csv"),
  row.names = FALSE, na = "")
saveRDS(synthetic_patients,
  file.path(output_dir, "synthetic_vebis_patients.rds"))
utils::write.csv(country_calendar,
  file.path(output_dir, "synthetic_country_calendar.csv"),
  row.names = FALSE, na = "")

# --- 10. Controle rapide : uniquement sur donnees fictives --------------------
cat("\n=== DONNEES SYNTHETIQUES : CONTROLES ===\n")
cat("Patients generes : ", nrow(synthetic_patients), "\n", sep = "")
cat("Pays : ", paste(sort(unique(country)), collapse = ", "), "\n", sep = "")
cat("SARI-CD1 (avant exclusions) : ", sum(cd1_for_qc), "\n", sep = "")
cat("SARI-CD2 (avant exclusions) : ", sum(cd2_for_qc), "\n", sep = "")
cat("Sous-ensemble CD2 inclus dans CD1 : OUI\n")
cat("\nResultats influenza (NA inclus) :\n")
print(table(synthetic_patients$flu_result, useNA = "ifany"))
cat("\nResultats SARS-CoV-2 (NA inclus) :\n")
print(table(synthetic_patients$covid_result, useNA = "ifany"))
cat("\nFichiers crees dans : ", normalizePath(output_dir), "\n", sep = "")
cat("RAPPEL : simulation pedagogique, pas donnees ni VE reelles VEBIS.\n")
