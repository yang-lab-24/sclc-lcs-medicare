# 02_sclc_study.R
#
# Sections:
#   I.   Setup and variable definitions
#   II.  Survival information
#   III. Unadjusted analyses
#   IV.  Propensity score matched analysis
#   V.   Covariate balance
#   VI.  Delivered dataset
#
# Last edited: 2026-09-01

library(arrow)
library(dplyr)
library(tidyr)
library(stringr)
library(survival)
library(MatchIt)
library(cobalt)
library(tableone)
library(readr)

source("utilities/CharlsonDeyo.R")


# =============================================================================
# I. Setup and variable definitions
# =============================================================================

DATA_PATH  <- "output/sclc/sclc_cohort.parquet"
# mirrors 01_sclc_extraction.R so a stale cohort file cannot carry the patient through
DEATH_DATE_EXCLUDE_IDS <- c("lnK2022U0887795")
OUTPUT_DIR <- "output/sclc/analysis"
DX_PATH_CANDIDATES <- c(
  "output/extracted/lung_diagnosis.parquet",
  "output/sclc/extracted/lung_diagnosis.parquet"
)

IL_FIPS <- "17"

RISK_ADJUSTMENT_COVARS <- c("charlson_bin", "nci_index")
RISK_TABLE_VARS <- c("charlson_score", "charlson_bin", "nci_index")

# duplicates dropped from the delivered dataset; screened is recoverable as scan_type == "Screened"
DELIVERY_DROP_VARS <- c(
  "ajcc_staging", "missing_ajcc", "sex_n", "no_cod_info",
  "fup_months", "tumor_size_clean", "screened"
)

IMAGING_TABLE_VARS <- c("brain_mri_pre12post3", "brain_mri_pre12post3_prior",
                        "pet_pre12post3")

TREATMENT_TABLE_VARS <- c("chemo_post12", "immuno_post12",
                          "radiation_post12", "resection_post12")

# gcc is modality-level; the component variables are descriptive and none enters gcc
GCC_VARS           <- c("gcc", "gcc_180d")
GCC_COMPONENT_VARS <- c("platinum_post12", "etoposide_post12",
                        "irinotecan_post12", "platinum_etoposide_post12",
                        "platinum_doublet_es_post12", "io_firstline_post12",
                        "resection_anatomic_post12",
                        "chemort_any_post12", "chemort_concurrent_post12",
                        "chemort_concurrent_tight_post12",
                        "chemort_sequential_post12",
                        "no_treatment_post12")
GCC_STRATA_VARS    <- c("gcc_rule", "gcc_stage_basis", "sclc_stage_group",
                        "nccn_stage", "gcc_category", "gcc_guideline")

FUP_RESTRICT_VAR   <- "tx_fup_analytic"
FUP_COMPONENT_VARS <- c("img_fup_complete", "img_fup_observable", "img_fup_analytic",
                        "tx_fup_complete",  "tx_fup_observable",  "tx_fup_analytic")
IMAGING_EXTRA_COLS <- c("brain_mri_pre12post3_count", "pet_pre12post3_count",
                        "brain_mri_n_pre",            "pet_n_pre",
                        "brain_mri_n_post",           "pet_n_post",
                        "brain_mri_days_from_dx",     "pet_days_from_dx",
                        "brain_mri_first_date",       "pet_first_date",
                        "brain_mri_last_date",        "pet_last_date")

COMORB_PREFIXES <- if (exists("CHARLSON_DEYO_RULES")) {
  names(CHARLSON_DEYO_RULES)
} else {
  c(
    "acute_mi", "history_mi", "chf", "pvd", "cvd", "copd", "dementia",
    "paralysis", "diabetes", "diabetes_comp", "renal_disease",
    "mild_liver_disease", "liver_disease", "ulcers", "rheum_disease", "aids"
  )
}
COMORB_COUNT_COLS  <- c(paste0(COMORB_PREFIXES, "_pre_count"), paste0(COMORB_PREFIXES, "_pre_month_count"))
COMORB_SOURCE_COLS <- c(paste0(COMORB_PREFIXES, "_medpar"), paste0(COMORB_PREFIXES, "_otherfile"))
COMORB_STATUS_COLS <- paste0(COMORB_PREFIXES, "_status")
COMORB_BINARY_COLS <- paste0(COMORB_PREFIXES, "_binary")

# character covariates enter the propensity model as factors; age and yr as continuous terms
PSM_COVARS <- c(
  "age", "sex_n", "race_n", "yr", "urban_rural_n",
  "YOST_ACS_2013_2017_QUINTILE", "MEDIANHOUSEHOLDINCOMEINFLATION", "STATE",
  RISK_ADJUSTMENT_COVARS
)


# fills a column missing from an older cohort file so the script still runs
guard_col <- function(df, col, fill) {
  if (!col %in% names(df)) df[[col]] <- fill
  df
}

km_landmarks <- function(fit, times = c(12, 24, 36)) {
  s <- summary(fit, times = times)
  cat(sprintf("  %-20s  %6s  %6s  %13s\n", "Group", "Time", "Surv", "95% CI"))
  cat("  ", strrep("-", 52), "\n", sep = "")
  for (i in seq_along(s$time)) {
    grp <- if (!is.null(s$strata)) as.character(s$strata[i]) else ""
    cat(sprintf("  %-20s  %6.0f  %6.3f  (%.3f, %.3f)\n",
                grp, s$time[i], s$surv[i], s$lower[i], s$upper[i]))
  }
}

km_plot <- function(fit, title,
                    xlab  = "Time (months)",
                    ylab  = "Survival Probability",
                    xlim  = c(0, 36),
                    labs  = c("Unscreened", "Screened")) {
  plot(fit, col = c("steelblue", "firebrick"), lty = 1, lwd = 2,
       xlab = xlab, ylab = ylab, main = title,
       xlim = xlim, conf.int = FALSE, mark.time = FALSE)
  legend("topright", legend = labs,
         col = c("steelblue", "firebrick"), lty = 1, lwd = 2, bty = "n")
}


cat("Loading cohort...\n")
cohort <- read_parquet(DATA_PATH)
n_death_excluded <- sum(cohort$PATIENT_ID %in% DEATH_DATE_EXCLUDE_IDS)
if (n_death_excluded > 0) {
  cohort <- cohort %>% filter(!PATIENT_ID %in% DEATH_DATE_EXCLUDE_IDS)
  cat(sprintf("  Dropped %d patient(s) with an unrepairable MBSF death date (re-run 01_sclc_extraction.R)\n",
              n_death_excluded))
}
cat(sprintf("  N = %d\n", nrow(cohort)))

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)


cat("Defining variables...\n")

cohort <- cohort %>%
  guard_col("SEER_REGISTRY",                  NA_character_)  %>%
  guard_col("LATERALITY",                      NA_integer_)  %>%
  guard_col("SEERCAUSESPECIFICDEATHCLASSIFIC", NA_integer_)  %>%
  guard_col("SEEROTHERCAUSEOFDEATHCLASSIFICA", NA_integer_)  %>%
  guard_col("YOST_ACS_2013_2017_QUINTILE",     NA_character_) %>%
  guard_col("MEDIANHOUSEHOLDINCOMEINFLATION",  NA_character_) %>%
  guard_col("STATE",                           NA_character_)

for (col in c(IMAGING_TABLE_VARS, TREATMENT_TABLE_VARS,
              GCC_VARS, GCC_COMPONENT_VARS,
              "gcc_eligible", "early_death_30d", "early_death_60d")) {
  cohort <- guard_col(cohort, col, NA_integer_)
}
for (col in GCC_STRATA_VARS) {
  cohort <- guard_col(cohort, col, NA_character_)
}
cohort <- guard_col(cohort, FUP_RESTRICT_VAR, NA_integer_)

if (all(is.na(cohort$gcc))) {
  cat("  Warning: gcc is entirely missing -- re-run 01_sclc_extraction.R to derive it.\n")
}

# 395 days back and 30 forward from dx_date_end with the rule-out filter; 0 if unavailable
dx_path <- DX_PATH_CANDIDATES[file.exists(DX_PATH_CANDIDATES)][1]
risk_scores <- NULL

if (is.na(dx_path)) {
  cat("  Warning: diagnosis claims parquet not found; Charlson-Deyo scores set to 0.\n")
} else {
  cat(sprintf("Computing Charlson-Deyo risk adjustment from: %s\n", dx_path))
  risk_scores <- tryCatch(
    {
      dx_claims <- read_parquet(dx_path) %>%
        filter(PATIENT_ID %in% cohort$PATIENT_ID) %>%
        select(PATIENT_ID, CLM_THRU_DT, any_of("CLM_FROM_DT"), source, matches("^ICD_DGNS|^PRNCPAL_DGNS|^FST_DGNS"))

      calculate_charlson_deyo(
        dx_claims = dx_claims,
        seer_df = cohort %>% transmute(PATIENT_ID, dx.date.end = as.Date(dx_date_end)),
        window_before_days = 395L,
        window_after_days = 30L,
        ruleout = TRUE
      ) %>%
        rename(
          charlson_score = Charlson,
          charlson_bin   = charlson.bin,
          nci_index      = NCI_index
        )
    },
    error = function(e) {
      cat(sprintf("  Warning: Charlson-Deyo scoring failed (%s); scores set to 0.\n", e$message))
      NULL
    }
  )
}

# extraction comorbidity columns sharing a name with the Charlson output get an extraction_ prefix
if (!is.null(risk_scores)) {
  clash <- setdiff(intersect(names(cohort), names(risk_scores)), "PATIENT_ID")
  if (length(clash) > 0) {
    cat(sprintf("  Renaming %d extraction column(s) to extraction_*: %s\n",
                length(clash), paste(clash, collapse = ", ")))
    cohort <- cohort %>% rename_with(~ paste0("extraction_", .x), all_of(clash))
  }
  cohort <- cohort %>% left_join(risk_scores, by = "PATIENT_ID")
}

cohort <- cohort %>%
  guard_col("charlson_score", 0) %>%
  guard_col("charlson_bin",   0L) %>%
  guard_col("nci_index",      0)

for (col in COMORB_COUNT_COLS) {
  cohort <- guard_col(cohort, col, 0L)
}
for (col in COMORB_SOURCE_COLS) {
  cohort <- guard_col(cohort, col, FALSE)
}
for (col in COMORB_STATUS_COLS) {
  cohort <- guard_col(cohort, col, "DNE")
}
for (col in COMORB_BINARY_COLS) {
  cohort <- guard_col(cohort, col, 0L)
}

cohort <- cohort %>%
  mutate(
    across(any_of(COMORB_COUNT_COLS),  ~ as.integer(coalesce(as.integer(.), 0L))),
    across(any_of(COMORB_SOURCE_COLS), ~ as.logical(coalesce(., FALSE))),
    across(any_of(COMORB_STATUS_COLS), ~ coalesce(as.character(.), "DNE")),
    across(any_of(COMORB_BINARY_COLS), ~ as.integer(coalesce(as.integer(.), 0L))),
    charlson_score = coalesce(charlson_score, 0),
    charlson_bin   = as.integer(coalesce(charlson_bin, 0L)),
    nci_index      = coalesce(nci_index, 0)
  )

cohort <- cohort %>%
  mutate(
    # missing_ajcc: registries 61-66; no_cod_info: 63 and 65; missing YOST quintile becomes level 9
    ajcc_staging = case_when(
      str_detect(coalesce(stage_simple, ""), "Missing Data|Data available") ~ "Missing/Unknown",
      is.na(stage_simple)                                                   ~ "Missing/Unknown",
      TRUE                                                                   ~ stage_simple
    ),
    missing_ajcc = as.integer(!is.na(SEER_REGISTRY) & between(as.integer(SEER_REGISTRY), 61L, 66L)),

    no_cod_info  = as.integer(!is.na(SEER_REGISTRY) & as.integer(SEER_REGISTRY) %in% c(63L, 65L)),

    tumor_size_clean = as.numeric(size),

    tumor_location = case_when(
      PRIMARY_SITE == "C341" & as.integer(LATERALITY) == 1 ~ "Right upper lobe",
      PRIMARY_SITE == "C342"                               ~ "Right middle lobe",
      PRIMARY_SITE == "C343" & as.integer(LATERALITY) == 1 ~ "Right lower lobe",
      PRIMARY_SITE == "C341" & as.integer(LATERALITY) == 2 ~ "Left upper lobe",
      PRIMARY_SITE == "C343" & as.integer(LATERALITY) == 2 ~ "Left lower lobe",
      as.integer(LATERALITY) == 4                          ~ "Bilateral",
      TRUE                                                  ~ NA_character_
    ),

    sex_n        = factor(sex,        levels = c("Male", "Female", "Unknown")),
    race_n       = factor(race,       levels = c("White", "Black", "Other or unknown")),
    urban_rural_n = factor(rural_urban, levels = c("Urban", "Rural", "Unknown")),

    YOST_ACS_2013_2017_QUINTILE = if_else(
      is.na(YOST_ACS_2013_2017_QUINTILE), "9", as.character(YOST_ACS_2013_2017_QUINTILE)
    )
  )

cat("Charlson-Deyo risk summary:\n")
print(cohort %>% count(charlson_bin) %>% mutate(pct = round(100 * n / sum(n), 1)))
cat(sprintf("Mean NCI index: %.3f\n", mean(cohort$nci_index, na.rm = TRUE)))


# =============================================================================
# II. Survival information
# =============================================================================

cohort <- cohort %>%
  mutate(
    fup_months     = tt_months_trunc,

    vital_status   = os_status,

    # cause-specific class 1 = lung cancer death; 8 and blank = NA; all else, including alive, = 0
    cod_class = na_if(str_trim(as.character(SEERCAUSESPECIFICDEATHCLASSIFIC)), ""),
    lung_cancer_death = case_when(
      is.na(cod_class) | cod_class == "8" ~ NA_integer_,
      cod_class == "1"                    ~ 1L,
      TRUE                                ~ 0L
    ),
    other_death = case_when(
      is.na(lung_cancer_death)                        ~ NA_integer_,
      vital_status == 1L & lung_cancer_death == 0L    ~ 1L,
      TRUE                                            ~ 0L
    )
  ) %>%
  select(-cod_class)


# =============================================================================
# III. Unadjusted analyses
# =============================================================================

cat("\n--- Section III: Unadjusted Analyses ---\n")

cat("\nScreening prevalence:\n")
print(cohort %>% count(screened) %>% mutate(pct = round(100 * n / sum(n), 1)))

cat("\nScreening by diagnosis year:\n")
print(
  cohort %>%
    group_by(yr) %>%
    summarize(
      n          = n(),
      n_screened = sum(screened == 1L),
      pct        = round(100 * n_screened / n, 1),
      .groups    = "drop"
    )
)

cat("\nSEER summary stage by screening status:\n")
stage_tab <- cohort %>%
  mutate(stage_label = case_when(
    staging == 1 ~ "Local",
    staging == 2 ~ "Regional",
    staging == 3 ~ "Distant",
    TRUE         ~ "Unknown"
  )) %>%
  group_by(screened, stage_label) %>%
  summarize(n = n(), .groups = "drop") %>%
  group_by(screened) %>%
  mutate(pct = round(100 * n / sum(n), 1))
print(stage_tab)

cat("\nAJCC stage by screening status (registries with AJCC data):\n")
ajcc_tab <- cohort %>%
  filter(missing_ajcc != 1) %>%
  group_by(screened, ajcc_staging) %>%
  summarize(n = n(), .groups = "drop") %>%
  group_by(screened) %>%
  mutate(pct = round(100 * n / sum(n), 1))
print(ajcc_tab)

cat("\nImaging and systemic therapy by screening status:\n")
summarize_measures <- function(d, label) {
  d %>%
    group_by(screened) %>%
    summarize(
      n = n(),
      across(any_of(c(IMAGING_TABLE_VARS, TREATMENT_TABLE_VARS)),
             ~ round(100 * mean(.x == 1L, na.rm = TRUE), 1),
             .names = "pct_{.col}"),
      .groups = "drop"
    ) %>%
    mutate(denominator = label, .before = 1)
}

imaging_by_screening <- bind_rows(
  summarize_measures(cohort, "Full cohort"),
  if (FUP_RESTRICT_VAR %in% names(cohort)) {
    summarize_measures(cohort %>% filter(.data[[FUP_RESTRICT_VAR]] == 1L),
                       "Complete 12-month post-dx coverage")
  }
)
print(imaging_by_screening)
cat("  (imaging window: 12 months before to 3 months after diagnosis;\n")
cat("   therapy window: diagnosis month through 12 months after)\n")
write_csv(imaging_by_screening, file.path(OUTPUT_DIR, "imaging_by_screening.csv"))

cat("\n--- NCCN guideline-concordant treatment ---\n")
if (all(is.na(cohort$gcc))) {
  cat("  Skipped: gcc not present in the cohort file.\n")
} else {
  cat(sprintf("  Scorable: %d of %d (%.1f%%); %d unstageable by both AJCC TNM and SEER Summary Stage\n",
              sum(!is.na(cohort$gcc)), nrow(cohort),
              100 * mean(!is.na(cohort$gcc)), sum(is.na(cohort$gcc))))

  if (!all(is.na(cohort$gcc_guideline))) {
    cat("  Guideline applied:\n")
    print(cohort %>% count(gcc_guideline, yr) %>%
            group_by(gcc_guideline) %>%
            summarize(dx_years = paste(range(yr), collapse = "-"),
                      n = sum(n), .groups = "drop"))
  }

  gcc_by_rule <- cohort %>%
    filter(!is.na(gcc_rule)) %>%
    group_by(gcc_rule, gcc_stage_basis) %>%
    summarize(n              = n(),
              pct_gcc        = round(100 * mean(gcc == 1L), 1),
              pct_gcc_180d   = round(100 * mean(gcc_180d == 1L), 1),
              pct_untreated  = round(100 * mean(no_treatment_post12 == 1L), 1),
              .groups = "drop") %>%
    arrange(gcc_rule, gcc_stage_basis)
  print(gcc_by_rule, n = Inf)
  write_csv(gcc_by_rule, file.path(OUTPUT_DIR, "gcc_by_stage_rule.csv"))

  tier_gap <- gcc_by_rule %>%
    group_by(gcc_rule) %>%
    summarize(gap = if (n() == 2) abs(diff(pct_gcc)) else NA_real_, .groups = "drop")
  cat(sprintf("  AJCC-vs-summary-stage gap within rule (pp): %s\n",
              paste(sprintf("%s %.1f", tier_gap$gcc_rule, tier_gap$gap), collapse = " | ")))

  cat("\n  By screening status:\n")
  # nested denominators: all scorable; complete coverage; complete coverage minus 60-day deaths
  gcc_denominators <- list(
    "All scorable"                          = rep(TRUE, nrow(cohort)),
    "Complete 12-mo coverage"               = cohort$gcc_eligible == 1L,
    "Complete coverage, excl. 60-day death" = cohort$gcc_eligible == 1L &
                                              cohort$early_death_60d == 0L
  )
  gcc_by_screening <- bind_rows(lapply(names(gcc_denominators), function(lbl) {
    keep <- (gcc_denominators[[lbl]] %in% TRUE) & !is.na(cohort$gcc)
    if (!any(keep)) return(NULL)
    cohort[keep, ] %>%
      group_by(screened) %>%
      summarize(n            = n(),
                pct_gcc      = round(100 * mean(gcc == 1L), 1),
                pct_gcc_180d = round(100 * mean(gcc_180d == 1L), 1),
                .groups = "drop") %>%
      mutate(denominator = lbl, .before = 1)
  }))
  print(gcc_by_screening, n = Inf)
  write_csv(gcc_by_screening, file.path(OUTPUT_DIR, "gcc_by_screening.csv"))

  cat("\n  Three-level category by screening status:\n")
  gcc_cat <- cohort %>%
    filter(!is.na(gcc_category)) %>%
    count(screened, gcc_category) %>%
    group_by(screened) %>%
    mutate(pct = round(100 * n / sum(n), 1)) %>%
    ungroup()
  print(gcc_cat, n = Inf)
  write_csv(gcc_cat, file.path(OUTPUT_DIR, "gcc_category_by_screening.csv"))

  cat("\n  By diagnosis year:\n")
  print(cohort %>%
    filter(!is.na(gcc)) %>%
    group_by(yr) %>%
    summarize(n            = n(),
              pct_gcc      = round(100 * mean(gcc == 1L), 1),
              pct_gcc_180d = round(100 * mean(gcc_180d == 1L), 1),
              .groups = "drop"))

  cat("\n  Treatment detail (descriptive; none of these enters gcc):\n")
  print(cohort %>%
    summarize(across(all_of(GCC_COMPONENT_VARS),
                     ~ round(100 * mean(.x == 1L, na.rm = TRUE), 1))) %>%
    tidyr::pivot_longer(everything(), names_to = "component", values_to = "pct"),
    n = Inf)

  cat("\n  Risk-adjusted association with screening (logistic, scorable patients):\n")
  # logistic model of gcc on screening, unknown sex dropped; screened is reported as an odds ratio
  gcc_covars <- intersect(c("age", "sex_n", "race_n", "yr", "urban_rural_n",
                            RISK_ADJUSTMENT_COVARS, "gcc_rule"),
                          names(cohort))
  gcc_model <- tryCatch(
    glm(reformulate(c("screened", gcc_covars), response = "gcc"),
        data = cohort %>% filter(!is.na(gcc), sex != "Unknown"),
        family = binomial()),
    error = function(e) { cat(sprintf("  Model failed: %s\n", e$message)); NULL }
  )
  if (!is.null(gcc_model)) {
    gcc_or <- summary(gcc_model)$coefficients
    keep <- rownames(gcc_or) == "screened"
    cat(sprintf("  Screened vs not: OR %.2f (95%% CI %.2f-%.2f), p = %.4g\n",
                exp(gcc_or[keep, 1]),
                exp(gcc_or[keep, 1] - 1.96 * gcc_or[keep, 2]),
                exp(gcc_or[keep, 1] + 1.96 * gcc_or[keep, 2]),
                gcc_or[keep, 4]))
    write_csv(
      tibble::rownames_to_column(as.data.frame(gcc_or), "term"),
      file.path(OUTPUT_DIR, "gcc_logistic_model.csv")
    )
  }
}

cat("\nStage and first-course treatment by screening status:\n")
cat_vars <- intersect(
  c("seer_t_group", "seer_n_group", "seer_stage_group", "seer_stage_source",
    "radiation_yes_no_unk", "chemotherapy",
    "chemo_post12", "chemo_post12_enumerated",
    "chemo_post12_admin", "chemo_post12_drug",
    "chemo_post12_icd_dx", "chemo_post12_icd_px",
    "immuno_post12", "immuno_post12_icd_px", "immuno_post12_any",
    "radiation_post12", "radiation_post12_hcpcs",
    "radiation_post12_icd_px", "radiation_post12_icd_dx",
    "radiation_post12_any",
    "brain_mri_inpatient_icd_px", "pet_inpatient_icd_px",
    "resection_post12", "resection_post12_hcpcs", "resection_post12_icd_px",
    "screened_primary_only",
    GCC_VARS, GCC_STRATA_VARS,
    FUP_COMPONENT_VARS),
  names(cohort)
)
stage_tx_by_screening <- bind_rows(lapply(cat_vars, function(v) {
  cohort %>%
    count(screened, level = as.character(.data[[v]])) %>%
    mutate(variable = v)
})) %>%
  group_by(variable, screened) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  ungroup() %>%
  select(variable, level, screened, n, pct) %>%
  arrange(variable, level, screened)
print(stage_tx_by_screening, n = Inf)
write_csv(stage_tx_by_screening, file.path(OUTPUT_DIR, "stage_treatment_by_screening.csv"))

cat("\nOverall Survival — log-rank test:\n")
km_os <- survfit(Surv(fup_months, vital_status) ~ screened, data = cohort)
print(survdiff(Surv(fup_months, vital_status) ~ screened, data = cohort))
cat("Overall Survival at 12, 24, 36 months:\n")
km_landmarks(km_os)
km_plot(km_os, title = "SCLC Overall Survival by Screening Status")

cat("\nOverall Survival by Year Group (2015-2019 vs 2020-2021):\n")
cohort <- cohort %>% mutate(year_group = case_when(
  yr >= 2015 & yr <= 2019 ~ "2015-2019",
  yr >= 2020 & yr <= 2021 ~ "2020-2021",
  TRUE ~ NA_character_
))
for (yg in c("2015-2019", "2020-2021")) {
  cat(paste0("\nYear group: ", yg, "\n"))
  sub <- cohort %>% filter(year_group == yg)
  if (nrow(sub) > 0) {
    fit <- survfit(Surv(fup_months, vital_status) ~ screened, data = sub)
    print(survdiff(Surv(fup_months, vital_status) ~ screened, data = sub))
    cat("Overall Survival at 12, 24, 36 months:\n")
    km_landmarks(fit)
    km_plot(fit, title = paste("SCLC OS by Screening (", yg, ")"))
  } else {
    cat("  No data for this year group.\n")
  }
}

cat("\nOverall Survival by Illinois vs Other States:\n")
cohort <- cohort %>%
  # STATE is a FIPS code, zero-padded to two digits
  mutate(is_IL = if_else(str_pad(as.character(STATE), 2, pad = "0") == IL_FIPS,
                         "IL", "Other"))
for (il_grp in c("IL", "Other")) {
  cat(paste0("\nState group: ", il_grp, "\n"))
  sub <- cohort %>% filter(is_IL == il_grp)
  if (nrow(sub) > 0) {
    fit <- survfit(Surv(fup_months, vital_status) ~ screened, data = sub)
    print(survdiff(Surv(fup_months, vital_status) ~ screened, data = sub))
    cat("Overall Survival at 12, 24, 36 months:\n")
    km_landmarks(fit)
    km_plot(fit, title = paste("SCLC OS by Screening (", il_grp, ")"))
  } else {
    cat("  No data for this state group.\n")
  }
}

# cause-specific survival excludes registries without cause of death
has_cod <- any(!is.na(cohort$SEERCAUSESPECIFICDEATHCLASSIFIC))
if (has_cod) {
  cohort_cod <- cohort %>% filter(no_cod_info != 1)

  cat("\nLung Cancer-Specific Survival — log-rank test:\n")
  km_lc <- survfit(Surv(fup_months, lung_cancer_death) ~ screened, data = cohort_cod)
  print(survdiff(Surv(fup_months, lung_cancer_death) ~ screened, data = cohort_cod))
  cat("Lung cancer death at 12, 24, 36 months:\n")
  km_landmarks(km_lc)
  km_plot(km_lc, title = "SCLC Lung Cancer-Specific Survival by Screening Status")

  cat("\nOther-Cause Survival — log-rank test:\n")
  km_oth <- survfit(Surv(fup_months, other_death) ~ screened, data = cohort_cod)
  print(survdiff(Surv(fup_months, other_death) ~ screened, data = cohort_cod))
  cat("Other-cause death at 12, 24, 36 months:\n")
  km_landmarks(km_oth)
  km_plot(km_oth, title = "SCLC Other-Cause Survival by Screening Status")
}

if ("time_diff" %in% names(cohort)) {
  cat("\nTime between last LDCT and diagnosis (screened patients):\n")
  hist(cohort$time_diff[cohort$screened == 1], breaks = 40,
       xlab = "Days Between LDCT Screening Exam and Diagnosis Date",
       ylab = "Proportion", freq = FALSE,
       main = "LDCT-to-Diagnosis Interval")
}

cat("\nBaseline characteristics (unmatched):\n")
tab_vars <- intersect(
  c("age", "sex_n", "race_n", "yr", "urban_rural_n",
    "YOST_ACS_2013_2017_QUINTILE", "MEDIANHOUSEHOLDINCOMEINFLATION",
    RISK_TABLE_VARS, IMAGING_TABLE_VARS, TREATMENT_TABLE_VARS, GCC_VARS),
  names(cohort)
)
tbl1 <- CreateTableOne(
  vars       = tab_vars,
  factorVars = intersect(c(IMAGING_TABLE_VARS, TREATMENT_TABLE_VARS, GCC_VARS), tab_vars),
  strata     = "screened",
  data       = cohort %>% filter(sex != "Unknown"),
  addOverall = TRUE
)
print(tbl1, showAllLevels = TRUE, quote = FALSE, noSpaces = TRUE)
write_csv(
  tibble::rownames_to_column(
    as.data.frame(print(tbl1, showAllLevels = TRUE, printToggle = FALSE,
                        quote = FALSE, noSpaces = TRUE)),
    "variable"
  ),
  file.path(OUTPUT_DIR, "table1_unmatched.csv")
)


# =============================================================================
# IV. Propensity score matched analysis
# =============================================================================

cat("\n--- Section IV: Propensity Score Matched Analysis ---\n")

# 1:1 nearest neighbour without replacement, caliper 0.01 on the probability scale, fixed seed
available_psm <- intersect(PSM_COVARS, names(cohort))
psm_data <- cohort %>%
  filter(sex != "Unknown") %>%
  select(PATIENT_ID, screened, fup_months, vital_status,
         lung_cancer_death, other_death, no_cod_info,
         staging, ajcc_staging, missing_ajcc,
         all_of(available_psm),
         any_of(c(RISK_TABLE_VARS, IMAGING_TABLE_VARS,
                  TREATMENT_TABLE_VARS, GCC_VARS, GCC_STRATA_VARS,
                  "gcc_eligible", "early_death_60d", FUP_RESTRICT_VAR))) %>%
  tidyr::drop_na(screened, age)

psm_formula <- reformulate(available_psm, response = "screened")

set.seed(123456)
psm_fit <- tryCatch(
  matchit(
    formula     = psm_formula,
    data        = psm_data,
    method      = "nearest",
    distance    = "glm",
    link        = "logit",
    caliper     = 0.01,
    std.caliper = FALSE,
    replace     = FALSE,
    ratio       = 1
  ),
  error = function(e) { cat(sprintf("PSM failed: %s\n", e$message)); NULL }
)

if (!is.null(psm_fit)) {
  matched <- match.data(psm_fit)
  cat(sprintf("  Matched cohort: %d patients (%d screened / %d unscreened pairs)\n",
              nrow(matched),
              sum(matched$screened == 1),
              sum(matched$screened == 0)))

  cat("\nSEER summary stage by screening status (matched):\n")
  stage_tab_m <- matched %>%
    mutate(stage_label = case_when(
      staging == 1 ~ "Local",
      staging == 2 ~ "Regional",
      staging == 3 ~ "Distant",
      TRUE         ~ "Unknown"
    )) %>%
    group_by(screened, stage_label) %>%
    summarize(n = n(), .groups = "drop") %>%
    group_by(screened) %>%
    mutate(pct = round(100 * n / sum(n), 1))
  print(stage_tab_m)

  cat("\nAJCC stage by screening status (matched, registries with AJCC data):\n")
  ajcc_tab_m <- matched %>%
    filter(missing_ajcc != 1) %>%
    group_by(screened, ajcc_staging) %>%
    summarize(n = n(), .groups = "drop") %>%
    group_by(screened) %>%
    mutate(pct = round(100 * n / sum(n), 1))
  print(ajcc_tab_m)

  cat("\nMatched Overall Survival — stratified log-rank test:\n")
  # log-rank tests stratified on the matched pair (subclass)
  km_os_m <- survfit(Surv(fup_months, vital_status) ~ screened, data = matched)
  print(survdiff(Surv(fup_months, vital_status) ~ screened + strata(subclass), data = matched))
  cat("Matched Overall Survival at 12, 24, 36 months:\n")
  km_landmarks(km_os_m)
  km_plot(km_os_m, title = "SCLC Overall Survival by Screening Status (Matched)",
          labs = c("Not Screened", "Screened"))

  if (has_cod) {
    matched_cod <- matched %>% filter(no_cod_info != 1)

    cat("\nMatched Lung Cancer-Specific Survival — stratified log-rank test:\n")
    km_lc_m <- survfit(Surv(fup_months, lung_cancer_death) ~ screened, data = matched_cod)
    print(survdiff(Surv(fup_months, lung_cancer_death) ~ screened + strata(subclass),
                   data = matched_cod))
    cat("Matched lung cancer death at 12, 24, 36 months:\n")
    km_landmarks(km_lc_m)
    km_plot(km_lc_m, title = "SCLC Lung Cancer-Specific Survival (Matched)",
            labs = c("Not Screened", "Screened"))

    cat("\nMatched Other-Cause Survival — stratified log-rank test:\n")
    km_oth_m <- survfit(Surv(fup_months, other_death) ~ screened, data = matched_cod)
    print(survdiff(Surv(fup_months, other_death) ~ screened + strata(subclass),
                   data = matched_cod))
    cat("Matched other-cause death at 12, 24, 36 months:\n")
    km_landmarks(km_oth_m)
    km_plot(km_oth_m, title = "SCLC Other-Cause Survival (Matched)",
            labs = c("Not Screened", "Screened"))
  }

  cat("\nBaseline characteristics (matched):\n")
  tab_vars_m <- intersect(tab_vars, names(matched))
  tbl1_m <- CreateTableOne(
    vars       = tab_vars_m,
    factorVars = intersect(c(IMAGING_TABLE_VARS, TREATMENT_TABLE_VARS, GCC_VARS), tab_vars_m),
    strata     = "screened",
    data       = matched,
    addOverall = TRUE
  )
  print(tbl1_m, showAllLevels = TRUE, quote = FALSE, noSpaces = TRUE)
  write_csv(
    tibble::rownames_to_column(
      as.data.frame(print(tbl1_m, showAllLevels = TRUE, printToggle = FALSE,
                          quote = FALSE, noSpaces = TRUE)),
      "variable"
    ),
    file.path(OUTPUT_DIR, "table1_matched.csv")
  )
}


# =============================================================================
# V. Covariate balance
# =============================================================================

if (!is.null(psm_fit)) {
  cat("\n--- Section V: Covariate Balance ---\n")

  ps_all <- data.frame(
    pscore  = psm_fit$distance,
    treated = as.integer(psm_data$screened)
  )
  plot(
    density(ps_all$pscore[ps_all$treated == 1]),
    col = "firebrick", lwd = 2, lty = 1,
    main = "Propensity Score Distribution Before Matching",
    xlab = "Propensity Score"
  )
  lines(density(ps_all$pscore[ps_all$treated == 0]),
        col = "steelblue", lwd = 2, lty = 2)
  legend("topright", legend = c("Screened", "Unscreened"),
         col = c("firebrick", "steelblue"), lwd = 2, lty = c(1, 2), bty = "n")

  ps_matched <- matched %>% select(distance, screened)
  plot(
    density(ps_matched$distance[ps_matched$screened == 1]),
    col = "firebrick", lwd = 2, lty = 1,
    main = "Propensity Score Distribution After Matching",
    xlab = "Propensity Score"
  )
  lines(density(ps_matched$distance[ps_matched$screened == 0]),
        col = "steelblue", lwd = 2, lty = 2)
  legend("topright", legend = c("Screened", "Unscreened"),
         col = c("firebrick", "steelblue"), lwd = 2, lty = c(1, 2), bty = "n")

  cat("\nStandardized mean differences (before and after matching):\n")
  bal <- bal.tab(psm_fit, un = TRUE, binary = "std", thresholds = c(m = 0.1))
  print(bal)

  love.plot(
    psm_fit,
    binary    = "std",
    thresholds = c(m = 0.1),
    title     = "Covariate Balance Before and After PSM",
    colors    = c("steelblue", "firebrick")
  )

  write_csv(
    tibble::rownames_to_column(as.data.frame(bal$Balance), "covariate"),
    file.path(OUTPUT_DIR, "psm_balance_smd.csv")
  )
}


# =============================================================================
# VI. Delivered dataset
# =============================================================================

cohort_out <- cohort %>% select(-any_of(DELIVERY_DROP_VARS))
dropped <- intersect(DELIVERY_DROP_VARS, names(cohort))
cat(sprintf("Dropping %d duplicate variable(s) from the delivered dataset: %s\n",
            length(dropped), paste(dropped, collapse = ", ")))

cohort_out_parquet <- file.path(OUTPUT_DIR, "final_analysis_cohort.parquet")
write_parquet(cohort_out, cohort_out_parquet)
write_csv(cohort_out, file.path(OUTPUT_DIR, "final_analysis_cohort_info.csv"))

charlson_cols <- intersect(
  c("charlson_score", "charlson_bin", "nci_index",
    COMORB_COUNT_COLS, COMORB_SOURCE_COLS, COMORB_STATUS_COLS, COMORB_BINARY_COLS),
  names(cohort_out)
)
cat(sprintf("Saved enriched cohort (%d rows, %d columns, %d Charlson-Deyo) to:\n  %s\n  %s\n",
            nrow(cohort_out), ncol(cohort_out), length(charlson_cols),
            cohort_out_parquet,
            file.path(OUTPUT_DIR, "final_analysis_cohort_info.csv")))

cat("\nAnalysis complete.\n")
