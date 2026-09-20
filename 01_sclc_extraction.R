# 01_sclc_extraction.R
#
# Sections:
#   I.    Configuration and code lists
#   II.   Helpers
#   III.  SEER: cohort selection
#   IV.   MBSF: enrollment, death and follow-up windows
#   V.    SEER variable recoding
#   VI.   Claims: screening, peri-diagnosis imaging and treatment
#   VII.  NCCN guideline-concordant treatment
#   VIII. Comorbidity flags
#   IX.   Outputs
#
# Last edited: 2026-09-01

library(dplyr)
library(tidyr)
library(lubridate)
library(stringr)
library(arrow)
library(readr)

source("code_definition_notes.R")


# =============================================================================
# I. Configuration and code lists
# =============================================================================

EXTRACTED  <- "output/extracted"
OUTPUT_DIR <- "output/sclc"
START_YEAR <- 2015L
END_YEAR   <- 2021L
LAST_FUP   <- ymd("2022-12-31")

# peri-diagnosis imaging window in days around dx_date_end
IMG_PRE_DAYS  <- 365L
IMG_POST_DAYS <- 90L

# claims dated on or after this carry ICD-10-CM diagnosis codes
ICD10_START <- ymd("2015-10-01")

# MBSF death date precedes later treatment claims and cannot be repaired; the patient is dropped
DEATH_DATE_EXCLUDE_IDS <- c("lnK2022U0887795")

SCLC_HIST <- c(8002L, 8041L:8045L)

# S8032 is a legacy code; screened_primary_only repeats the classification without it
LCS_CPTS_PRIMARY <- c("G0297", "71271")
LCS_CPTS_LEGACY   <- c("S8032")
LCS_CPTS  <- c(LCS_CPTS_PRIMARY, LCS_CPTS_LEGACY)

BRAIN_MRI_CPTS       <- c("70551", "70552", "70553")
# two-code subset used by the earlier brain MRI study
BRAIN_MRI_CPTS_PRIOR <- c("70551", "70553")

# inpatient (MEDPAR) imaging codes feed component columns only, never the primary flags
BRAIN_MRI_PX_ICD10_PREFIX <- "B030"
BRAIN_MRI_PX_ICD9         <- "8891"

PET_CPTS <- c("78811", "78812", "78813", "78814", "78815", "78816")

PET_PX_ICD10_PATTERN <- "^C.3"

# enumerated lists define chemo_post12_enumerated; the primary flag is range-based (CHEMO_*_RANGE)
CHEMO_ADMIN_CPTS <- c(
  "96401", "96402", "96405", "96406", "96409", "96411", "96413", "96415",
  "96416", "96417", "96420", "96422", "96423", "96425", "96440",
  "96446", "96450", "96521", "96522", "96523", "96542", "96549"
)

CHEMO_DRUG_HCPCS <- c(
  "J9060",
  "J9045",
  "J9181", "J8560",
  "J9351",
  "J9223",
  "J9206",
  "J9070",
  "J9000",
  "J9370", "J9371",
  "J9264", "J9265", "J9267",
  "J9171",
  "J9198", "J9201",
  "J9390",
  "J9033", "J9034", "J9036"
)

IMMUNO_HCPCS_SCLC <- c("J9022", "J9173", "J9299", "J9271")

# wider list, removed from the chemotherapy J-code range; defines no flag
IMMUNO_HCPCS_ALL <- c("J9271", "J9272", "J9299", "J9022", "J9023", "J9119",
                      "J9173", "J9228", "J9286", "J9355", "C9399")

# MEDPAR table 3E0, substance 0 with qualifier M/2/3; sensitivity columns only
IMMUNO_PX_ICD10_PATTERN <- "^3E0..0[M23]$"

# delivery codes define radiation_post12; supportive codes count only in radiation_post12_any
RADIATION_DELIVERY_CPT <- c(
  "77401", "77402", "77407", "77412",
  "G6003", "G6004", "G6005", "G6006", "G6007", "G6008", "G6009", "G6010",
  "G6011", "G6012", "G6013", "G6014",
  "77385", "77386", "G6015", "G6016",
  "77371", "77372", "77373",
  "77520", "77522", "77523", "77525",
  "77750",
  "77761", "77762", "77763", "77776", "77777", "77778", "77789", "77790",
  "77785", "77786", "77787",
  "77770", "77771", "77772",
  "77767", "77768"
)

RADIATION_SUPPORTIVE_CPT <- c(
  "77280", "77285", "77290", "77295",
  "77301", "77321",
  "77331", "77332", "77333", "77334", "77336", "77370",
  "77417",
  "77427", "77431", "77432", "77435", "77470",
  "77499"
)

RADIATION_CPT <- c(RADIATION_DELIVERY_CPT, RADIATION_SUPPORTIVE_CPT)

# ICD-10-PCS section D, excluding modality Y with qualifier 8 or K (thermal ablation)
RADIATION_PX_ICD9_PREFIX <- c("922", "923")

RADIATION_PX_ICD10_SECTION  <- "^D"
RADIATION_PX_ICD10_EXCL_MQ  <- c("8", "K")

RADIATION_DX_ICD10 <- c("Z510")
RADIATION_DX_ICD9  <- c("V580")

RESECTION_CPT <- c(
  "32440",
  "32442",
  "32445",
  "32488",
  "32671",
  "32480",
  "32482",
  "32486",
  "32663",
  "32670",
  "32484",
  "32669",
  "32505", "32506", "32507",
  "32666", "32667", "32668"
)

RESECTION_PX_ICD9_PREFIX <- c("3220", "3229", "323", "324", "325")

# root operation B (excision) and T (resection) of the lung, open and thoracoscopic approaches
RESECTION_PX_ICD10 <- c(
  "0BBC4ZX", "0BBC4ZZ", "0BBC0ZX", "0BBC0ZZ", "0BBD4ZX", "0BBD4ZZ", "0BBD0ZX", "0BBD0ZZ",
  "0BBF4ZX", "0BBF4ZZ", "0BBF0ZX", "0BBF0ZZ", "0BBG4ZX", "0BBG4ZZ", "0BBG0ZX", "0BBG0ZZ",
  "0BBH4ZX", "0BBH4ZZ", "0BBH0ZX", "0BBH0ZZ", "0BBJ4ZX", "0BBJ4ZZ", "0BBJ0ZX", "0BBJ0ZZ",
  "0BBK4ZX", "0BBK4ZZ", "0BBK0ZX", "0BBK0ZZ", "0BBL4ZX", "0BBL4ZZ", "0BBL0ZX", "0BBL0ZZ",
  "0BTC4ZZ", "0BTC0ZZ", "0BTD4ZZ", "0BTD0ZZ", "0BTF4ZZ", "0BTF0ZZ", "0BTG4ZZ", "0BTG0ZZ",
  "0BTH4ZZ", "0BTH0ZZ", "0BTJ4ZZ", "0BTJ0ZZ", "0BBM4ZX", "0BBM4ZZ", "0BBM0ZX", "0BBM0ZZ",
  "0BTM4ZZ", "0BTM0ZZ", "0BTK4ZZ", "0BTK0ZZ", "0BTL4ZZ", "0BTL0ZZ"
)

CHEMO_CPTS <- c(CHEMO_ADMIN_CPTS, CHEMO_DRUG_HCPCS)

# primary chemotherapy definition: 96400-96549, Q0083-Q0085 and J9000-J9999 minus IMMUNO_HCPCS_ALL
CHEMO_ADMIN_RANGE <- c(sprintf("%d", 96400:96549), "Q0083", "Q0084", "Q0085")

CHEMO_DRUG_RANGE <- setdiff(sprintf("J9%03d", 0:999), IMMUNO_HCPCS_ALL)

CHEMO_DX_ICD10 <- c("Z5111")
CHEMO_DX_ICD9  <- c("V581", "V5811")

# MEDPAR 3E0 with substance 0 and qualifier 0/5; XW codes are an allow-list, never a prefix match
CHEMO_PX_ICD9  <- c("9925", "0010", "1770", "9928")

CHEMO_PX_ICD10_PATTERN <- "^3E0..0[05]$"

CHEMO_PX_ICD10_NEWTECH <- c("XW03387", "XW04387", "XW033B3", "XW043B3")


# treatment-detail subsets are descriptive; none enters gcc
PLATINUM_HCPCS <- c("J9060", "J9062", "J9045")

ETOPOSIDE_HCPCS <- c("J9181", "J9182", "J8560")

IRINOTECAN_HCPCS <- c("J9206")

IO_FIRSTLINE_ES_HCPCS <- c("J9022", "J9173")

RESECTION_ANATOMIC_CPT <- c(
  "32440", "32442", "32445", "32488", "32671",
  "32480", "32482", "32486", "32663", "32670",
  "32484", "32669",
  "32507", "32668"
)

RESECTION_ANATOMIC_PX_ICD9_PREFIX <- c("323", "324", "325")

# root operation T subset of the lung list
RESECTION_ANATOMIC_PX_ICD10 <- RESECTION_PX_ICD10[
  str_starts(RESECTION_PX_ICD10, "0BT")
]

# load-time subset checks; a subset code outside the modality lists would give an all-zero flag
local({
  claimable <- c(CHEMO_CPTS, CHEMO_ADMIN_RANGE, CHEMO_DRUG_RANGE,
                 IMMUNO_HCPCS_ALL, RADIATION_CPT, RESECTION_CPT)
  for (nm in c("PLATINUM_HCPCS", "ETOPOSIDE_HCPCS", "IRINOTECAN_HCPCS",
               "IO_FIRSTLINE_ES_HCPCS", "RESECTION_ANATOMIC_CPT")) {
    missing <- setdiff(get(nm), claimable)
    if (length(missing) > 0) {
      stop(sprintf(paste0("%s contains code(s) not pulled from claims: %s. ",
                          "Add them to the modality list they belong to, or ",
                          "the flag will be silently all-zero."),
                   nm, paste(missing, collapse = ", ")))
    }
  }
  if (!all(RESECTION_ANATOMIC_CPT %in% RESECTION_CPT)) {
    stop("RESECTION_ANATOMIC_CPT must be a subset of RESECTION_CPT")
  }
  if (!all(IO_FIRSTLINE_ES_HCPCS %in% IMMUNO_HCPCS_SCLC)) {
    stop("IO_FIRSTLINE_ES_HCPCS must be a subset of IMMUNO_HCPCS_SCLC")
  }
  if (!all(RESECTION_ANATOMIC_PX_ICD10 %in% RESECTION_PX_ICD10)) {
    stop("RESECTION_ANATOMIC_PX_ICD10 must be a subset of RESECTION_PX_ICD10")
  }
  if (length(RESECTION_ANATOMIC_PX_ICD10) == 0L) {
    stop("RESECTION_ANATOMIC_PX_ICD10 is empty -- the 0BT filter matched nothing")
  }
  if (!all(RESECTION_ANATOMIC_PX_ICD9_PREFIX %in% RESECTION_PX_ICD9_PREFIX)) {
    stop("RESECTION_ANATOMIC_PX_ICD9_PREFIX must be a subset of RESECTION_PX_ICD9_PREFIX")
  }
})


PROC_CPTS <- unique(c(LCS_CPTS, BRAIN_MRI_CPTS, PET_CPTS,
                      CHEMO_CPTS, CHEMO_ADMIN_RANGE, CHEMO_DRUG_RANGE,
                      IMMUNO_HCPCS_ALL, RADIATION_CPT, RESECTION_CPT))

# descriptor for every enumerated code (observed_codes.csv); the script stops if one is missing
CODE_LABELS <- c(
  # Lung cancer screening LDCT
  "71271" = "CT, thorax, low dose for lung cancer screening, without contrast (2021+)",
  "G0297" = "Low dose CT scan (LDCT) for lung cancer screening (deleted 12/31/2020, replaced by 71271)",
  "S8032" = "Low-dose CT for lung cancer screening (added 10/2014, deleted 9/30/2016)",

  # Brain MRI
  "70551" = "MRI brain (including brain stem), without contrast",
  "70552" = "MRI brain (including brain stem), with contrast",
  "70553" = "MRI brain (including brain stem), without contrast followed by with contrast",

  # PET / PET-CT
  "78811" = "PET imaging, limited area (chest, head/neck)",
  "78812" = "PET imaging, skull base to mid-thigh",
  "78813" = "PET imaging, whole body",
  "78814" = "PET with concurrently acquired CT, limited area (chest, head/neck)",
  "78815" = "PET with concurrently acquired CT, skull base to mid-thigh",
  "78816" = "PET with concurrently acquired CT, whole body",

  # Chemotherapy administration
  "96401" = "Chemotherapy administration, SC/IM, non-hormonal anti-neoplastic",
  "96402" = "Chemotherapy administration, SC/IM, hormonal anti-neoplastic",
  "96405" = "Chemotherapy administration, intralesional, up to 7 lesions",
  "96406" = "Chemotherapy administration, intralesional, more than 7 lesions",
  "96409" = "Chemotherapy administration, IV push, single or initial substance",
  "96411" = "Chemotherapy administration, IV push, each additional substance",
  "96413" = "Chemotherapy administration, IV infusion, up to 1 hour, single/initial substance",
  "96415" = "Chemotherapy administration, IV infusion, each additional hour",
  "96416" = "Chemotherapy administration, IV infusion, initiation of prolonged infusion requiring pump",
  "96417" = "Chemotherapy administration, IV infusion, each additional sequential infusion",
  "96420" = "Chemotherapy administration, intra-arterial, push",
  "96422" = "Chemotherapy administration, intra-arterial infusion, up to 1 hour",
  "96423" = "Chemotherapy administration, intra-arterial infusion, each additional hour",
  "96425" = "Chemotherapy administration, intra-arterial infusion, prolonged, requiring pump",
  "96440" = "Chemotherapy administration into pleural cavity, requiring thoracentesis",
  "96446" = "Chemotherapy administration into peritoneal cavity via indwelling port/catheter",
  "96450" = "Chemotherapy administration into CNS, requiring lumbar puncture",
  "96521" = "Refilling and maintenance of portable pump",
  "96522" = "Refilling and maintenance of implantable pump or reservoir",
  "96523" = "Irrigation of implanted venous access device",
  "96542" = "Chemotherapy injection, subarachnoid/intraventricular via reservoir",
  "96549" = "Unlisted chemotherapy procedure",

  # Chemotherapy drugs (CMS HCPCS descriptors)
  "J9000" = "Injection, doxorubicin hydrochloride, 10 mg",
  "J9025" = "Injection, azacitidine, 1 mg",
  "J9045" = "Injection, carboplatin, 50 mg",
  "J9046" = "Injection, bortezomib (Dr. Reddy's), 0.1 mg (code created 2023)",
  "J9060" = "Injection, cisplatin, powder or solution, 10 mg",
  "J9062" = "Injection, cisplatin, powder or solution, 50 mg (retired; 0 claims in this extract)",
  "J9070" = "Cyclophosphamide, 100 mg (deleted 3/31/2024)",
  "J9171" = "Injection, docetaxel, 1 mg",
  "J9177" = "Injection, enfortumab vedotin-ejfv, 0.25 mg",
  "J9181" = "Injection, etoposide, 10 mg",
  "J9182" = "Injection, etoposide, 100 mg (retired; 0 claims in this extract)",
  "J8560" = "Etoposide, oral, 50 mg",
  "J9190" = "Injection, fluorouracil, 500 mg",
  "J9198" = "Injection, gemcitabine hydrochloride (Infugem), 100 mg (2020+)",
  "J9200" = "Injection, floxuridine, 500 mg",
  "J9201" = "Injection, gemcitabine hydrochloride, not otherwise specified, 200 mg",
  "J9205" = "Injection, irinotecan liposome, 1 mg",
  "J9206" = "Injection, irinotecan, 20 mg",
  "J9207" = "Injection, ixabepilone, 1 mg",
  "J9223" = "Injection, lurbinectedin, 0.1 mg (2021+)",
  "J9250" = "Methotrexate sodium, 5 mg (deleted 3/31/2024)",
  "J9264" = "Injection, paclitaxel protein-bound particles, 1 mg (nab-paclitaxel)",
  "J9265" = "Injection, paclitaxel, 30 mg (deleted 12/31/2014)",
  "J9267" = "Injection, paclitaxel, 1 mg",
  "J9280" = "Injection, mitomycin, 5 mg",
  "J9305" = "Injection, pemetrexed, not otherwise specified, 10 mg",
  "J9351" = "Injection, topotecan, 0.1 mg",
  "J9360" = "Injection, vinblastine sulfate, 1 mg",
  "J9370" = "Vincristine sulfate, 1 mg",
  "J9371" = "Injection, vincristine sulfate liposome, 1 mg (deleted 6/30/2024)",
  "J9390" = "Injection, vinorelbine tartrate, 10 mg",
  "J9033" = "Injection, bendamustine hydrochloride, 1 mg",
  "J9034" = "Injection, bendamustine hydrochloride (Bendeka), 1 mg (2017+)",
  "J9036" = "Injection, bendamustine hydrochloride (Belrapzo), 1 mg (2019+)",

  # Immunotherapy
  "J9022" = "Injection, atezolizumab, 10 mg",
  "J9023" = "Injection, avelumab, 10 mg",
  "J9119" = "Injection, cemiplimab-rwlc, 1 mg",
  "J9173" = "Injection, durvalumab, 10 mg",
  "J9228" = "Injection, ipilimumab, 1 mg",
  "J9271" = "Injection, pembrolizumab, 1 mg",
  "J9272" = "Injection, dostarlimab-gxly, 10 mg (code created 2022)",
  "J9286" = "Injection, glofitamab-gxbm, 2.5 mg (code created 2024; not tremelimumab)",
  "J9299" = "Injection, nivolumab, 1 mg",
  "J9355" = "Injection, trastuzumab, excludes biosimilar, 10 mg (not trastuzumab deruxtecan)",
  "C9399" = "Unclassified drugs or biologicals",

  # Radiation therapy, supportive
  "77280" = "Therapeutic radiology simulation-aided field setting, simple",
  "77285" = "Therapeutic radiology simulation-aided field setting, intermediate",
  "77290" = "Therapeutic radiology simulation-aided field setting, complex",
  "77295" = "3-dimensional radiotherapy plan, including dose-volume histograms",
  "77301" = "IMRT treatment planning",
  "77321" = "Special teletherapy port plan",
  "77331" = "Special dosimetry (e.g. abutting/overlapping fields, small-field dose confirmation)",
  "77332" = "Treatment devices, simple (e.g. simple block, simple bolus)",
  "77333" = "Treatment devices, intermediate (e.g. multiple blocks, stents, bite blocks)",
  "77334" = "Treatment devices, complex (e.g. irregular blocks, compensators, molds)",
  "77336" = "Continuing medical physics consultation, weekly",
  "77370" = "Special medical radiation physics consultation",
  "77417" = "Therapeutic radiology port image(s)",
  "77427" = "Radiation treatment management, five treatments",
  "77431" = "Radiation therapy management, one or two fractions only",
  "77432" = "Stereotactic radiation treatment management, cranial",
  "77435" = "Stereotactic body radiation therapy (SBRT), treatment management, up to 5 fractions",
  "77470" = "Special treatment procedure (unusual time/effort, e.g. TBI, hemibody)",
  "77499" = "Unlisted procedure, therapeutic radiology treatment management",

  # Radiation therapy, delivery
  "77401" = "Radiation treatment delivery, superficial/orthovoltage",
  "77402" = "Radiation treatment delivery, >1 MeV, simple/single area (revised 1/1/2015; OPPS -- MPFS uses G6003-G6006)",
  "77407" = "Radiation treatment delivery, >1 MeV, intermediate/2 areas (revised 1/1/2015; OPPS -- MPFS uses G6007-G6010)",
  "77412" = "Radiation treatment delivery, >1 MeV, complex/3+ areas (revised 1/1/2015; OPPS -- MPFS uses G6011-G6014)",
  "77385" = "IMRT delivery, simple",
  "77386" = "IMRT delivery, complex",
  "77371" = "Stereotactic radiosurgery, complete course, cranial, 1 session (replaces G0173, dead since 2015)",
  "77372" = "Stereotactic radiosurgery, linear accelerator based, delivery including collimator changes and custom plugging, fractionated treatment, all lesions, per course of treatment",
  "77373" = "Stereotactic body radiation therapy (SBRT), delivery, per fraction (replaces G0340, dead since 2015)",
  "G6003" = "Radiation treatment delivery, single area, single port, up to 5 MeV (MPFS)",
  "G6004" = "Radiation treatment delivery, single area, single port, 6-10 MeV (MPFS)",
  "G6005" = "Radiation treatment delivery, single area, single port, 11-19 MeV (MPFS)",
  "G6006" = "Radiation treatment delivery, single area, single/parallel-opposed ports, >=20 MeV (MPFS)",
  "G6007" = "Radiation treatment delivery, 2 separate treatment areas, up to 5 MeV (MPFS)",
  "G6008" = "Radiation treatment delivery, 2 separate treatment areas, 6-10 MeV (MPFS)",
  "G6009" = "Radiation treatment delivery, 2 separate treatment areas, 11-19 MeV (MPFS)",
  "G6010" = "Radiation treatment delivery, 2 separate treatment areas, >=20 MeV (MPFS)",
  "G6011" = "Radiation treatment delivery, 3+ separate treatment areas, up to 5 MeV (MPFS)",
  "G6012" = "Radiation treatment delivery, 3+ separate treatment areas, 6-10 MeV (MPFS)",
  "G6013" = "Radiation treatment delivery, 3+ separate treatment areas, 11-19 MeV (MPFS)",
  "G6014" = "Radiation treatment delivery, 3+ separate treatment areas, >=20 MeV (MPFS)",
  "G6015" = "IMRT delivery, dynamic MLC, hardware-based (MPFS equivalent of CPT 77385/77386)",
  "G6016" = "IMRT delivery, compensator-based modulated beam (MPFS equivalent of CPT 77385/77386)",
  "77520" = "Proton treatment delivery, simple, without compensation",
  "77522" = "Proton treatment delivery, simple, with compensation",
  "77523" = "Proton treatment delivery, intermediate",
  "77525" = "Proton treatment delivery, complex",
  "77750" = "Infusion or instillation of radioelement solution",
  "77761" = "Intracavitary radiation source application, simple",
  "77762" = "Intracavitary radiation source application, intermediate",
  "77763" = "Intracavitary radiation source application, complex",
  "77776" = "Interstitial radiation source application, simple",
  "77777" = "Interstitial radiation source application, intermediate",
  "77778" = "Interstitial radiation source application, complex",
  "77789" = "Surface application of radiation source",
  "77790" = "Supervision, handling, loading of radiation source",
  "77785" = "HDR brachytherapy delivery, 1 channel (through 2015; replaced by 77770 in 2016)",
  "77786" = "HDR brachytherapy delivery, 2-12 channels (through 2015; replaced by 77771 in 2016)",
  "77787" = "HDR brachytherapy delivery, >12 channels (through 2015; replaced by 77772 in 2016)",
  "77770" = "HDR brachytherapy delivery, 1 channel (2016+; replaces 77785)",
  "77771" = "HDR brachytherapy delivery, 2-12 channels (2016+; replaces 77786)",
  "77772" = "HDR brachytherapy delivery, >12 channels (2016+; replaces 77787)",
  "77767" = "Skin-surface HDR brachytherapy delivery, <=1 cm skin distance (2016+)",
  "77768" = "Skin-surface HDR brachytherapy delivery, >1 cm skin distance (2016+)",

  # Lung resection
  "32440" = "Pneumonectomy",
  "32442" = "Sleeve pneumonectomy",
  "32445" = "Extrapleural pneumonectomy",
  "32488" = "Completion pneumonectomy",
  "32671" = "Thoracoscopic (VATS) pneumonectomy",
  "32480" = "Lobectomy",
  "32482" = "Bilobectomy",
  "32486" = "Sleeve lobectomy",
  "32663" = "Thoracoscopic (VATS) lobectomy",
  "32670" = "Thoracoscopic (VATS) bilobectomy",
  "32484" = "Segmentectomy",
  "32669" = "Thoracoscopic (VATS) segmentectomy",
  "32505" = "Thoracotomy, with therapeutic wedge resection (e.g. mass, nodule), initial",
  "32506" = "Thoracotomy, with therapeutic wedge resection, each additional resection",
  "32507" = "Thoracotomy, with diagnostic wedge resection followed by anatomic lung resection",
  "32666" = "Thoracoscopic (VATS) therapeutic wedge resection, initial unilateral",
  "32667" = "Thoracoscopic (VATS) therapeutic wedge resection, each additional resection",
  "32668" = "Thoracoscopic (VATS) diagnostic wedge resection followed by anatomic lung resection"
)


# =============================================================================
# II. Helpers
# =============================================================================

claim_date_of <- function(x) as_date(x, tz = "UTC")

# MBSF wide to patient-month; ffs_month = non-HMO (indicator 0/4) with Part A and B buy-in (3/C)
mbsf_to_long <- function(d) {
  hmo <- d %>%
    select(PATIENT_ID, BENE_ENROLLMT_REF_YR, starts_with("HMO_IND_")) %>%
    pivot_longer(starts_with("HMO_IND_"),
                 names_to = "col", values_to = "hmo_ind") %>%
    mutate(month = as.integer(str_extract(col, "\\d+$"))) %>%
    select(-col)

  buyin <- d %>%
    select(PATIENT_ID, BENE_ENROLLMT_REF_YR, starts_with("MDCR_ENTLMT_BUYIN_IND_")) %>%
    pivot_longer(starts_with("MDCR_ENTLMT_BUYIN_IND_"),
                 names_to = "col", values_to = "buyin") %>%
    mutate(month = as.integer(str_extract(col, "\\d+$"))) %>%
    select(-col)

  hmo %>%
    inner_join(buyin, by = c("PATIENT_ID", "BENE_ENROLLMT_REF_YR", "month")) %>%
    mutate(
      ref_yr    = as.integer(BENE_ENROLLMT_REF_YR),
      date      = ceiling_date(
                    ymd(sprintf("%04d-%02d-01", ref_yr, month)), "month"
                  ) - days(1),
      ffs_month = as.integer(hmo_ind) %in% c(0L, 4L) & buyin %in% c("3", "C")
    ) %>%
    select(-ref_yr)
}


# =============================================================================
# III. SEER: cohort selection
# =============================================================================
cat("Loading SEER...\n")

# dx_date_end anchors enrollment, screening, imaging, comorbidities; dx_date_first, treatment
seer_all <- read_parquet(file.path(EXTRACTED, "seer_lung.parquet")) %>%
  mutate(
    yr  = as.integer(YEAR_OF_DIAGNOSIS),
    mo  = as.integer(MONTH_OF_DIAGNOSIS),
    age = as.numeric(AGERECODEWITHSINGLE_AGES_AND_90),
    dx_date       = ymd(if_else(
      !is.na(yr) & !is.na(mo),
      sprintf("%04d%02d15", yr, mo),
      NA_character_
    )),
    dx_date_end   = ceiling_date(
                      ymd(sprintf("%04d-%02d-01", yr, mo)), "month"
                    ) - days(1),
    dx_date_first = floor_date(dx_date, "month")
  )

seer <- seer_all %>%
  filter(
    yr >= START_YEAR,
    yr <= END_YEAR,
    as.integer(HISTOLOGIC_TYPE_ICD_O_3) %in% SCLC_HIST,
    !is.na(dx_date),
    age >= 65,
    age <= 78
  )

# attrition counts for the flow diagram
cat(sprintf("  SCLC dx %d-%d (any age, known/unknown dx date): %d\n",
  START_YEAR, END_YEAR,
  seer_all %>%
    filter(yr >= START_YEAR, yr <= END_YEAR,
           as.integer(HISTOLOGIC_TYPE_ICD_O_3) %in% SCLC_HIST) %>%
    nrow()))
cat(sprintf("  Missing diagnosis date: %d\n",
  seer_all %>%
    filter(yr >= START_YEAR, yr <= END_YEAR,
           as.integer(HISTOLOGIC_TYPE_ICD_O_3) %in% SCLC_HIST,
           is.na(dx_date)) %>%
    nrow()))
cat(sprintf("  Outside age 65-78: %d\n",
  seer_all %>%
    filter(yr >= START_YEAR, yr <= END_YEAR,
           as.integer(HISTOLOGIC_TYPE_ICD_O_3) %in% SCLC_HIST,
           !is.na(dx_date), (age < 65 | age > 78)) %>%
    nrow()))
cat(sprintf("  After inclusion criteria: %d\n", nrow(seer)))

pids <- unique(seer$PATIENT_ID)


# =============================================================================
# IV. MBSF: enrollment, death and follow-up windows
# =============================================================================
cat("Loading MBSF...\n")

mbsf_raw <- read_parquet(file.path(EXTRACTED, "mbsf.parquet")) %>%
  filter(PATIENT_ID %in% pids)

# latest MBSF death date wins; an erroneous earlier-year date is not overwritten when corrected
mbsf_death <- mbsf_raw %>%
  mutate(BENE_DEATH_DT = ymd(as.character(BENE_DEATH_DT))) %>%
  group_by(PATIENT_ID) %>%
  summarize(
    n_death_dates = n_distinct(BENE_DEATH_DT[!is.na(BENE_DEATH_DT)]),
    death_date    = if (n_death_dates == 0L) as.Date(NA) else max(BENE_DEATH_DT, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(last_fup_date = if_else(!is.na(death_date), death_date, LAST_FUP))

cat(sprintf("  Conflicting BENE_DEATH_DT across MBSF years (latest used): %d patient(s)\n",
            sum(mbsf_death$n_death_dates > 1L)))
mbsf_death <- mbsf_death %>% select(-n_death_dates)

mbsf_long <- mbsf_to_long(mbsf_raw)


cat("Checking enrollment...\n")

# fee-for-service months in the 365 days before and after dx_date_end
seer_dx <- seer %>% select(PATIENT_ID, dx_date_end)

mbsf_seer <- mbsf_long %>%
  right_join(seer_dx, by = "PATIENT_ID") %>%
  mutate(
    days_before = as.numeric(dx_date_end - date),
    days_after  = as.numeric(date - dx_date_end),
    prior_month = days_before >= 0 & days_before < 365,
    post_month  = days_after  >= 0 & days_after  < 365
  )

enrollment <- mbsf_seer %>%
  group_by(PATIENT_ID) %>%
  summarize(
    pre_months  = sum(ffs_month & prior_month, na.rm = TRUE),
    post_months = sum(ffs_month & post_month,  na.rm = TRUE),
    .groups = "drop"
  )

cat("Pre-dx FFS month distribution:\n")
print(enrollment %>% count(pre_months))

mbsf_final <- enrollment %>% inner_join(mbsf_death, by = "PATIENT_ID")

# the TRUE ~ 9L branch of os_status is unreachable and is kept as written
cohort_pre <- seer %>%
  inner_join(mbsf_final, by = "PATIENT_ID") %>%
  mutate(
    os_status = case_when(
      !is.na(death_date) ~ 1L,
      is.na(death_date)  ~ 0L,
      TRUE               ~ 9L
    ),
    tt_days         = as.numeric(difftime(last_fup_date, dx_date_end, units = "days")),
    tt_months       = tt_days / 30
  )

cohort <- cohort_pre %>% filter(pre_months == 12)

cat(sprintf("  Excluded (<12 pre-dx FFS months): %d\n", nrow(cohort_pre) - nrow(cohort)))

n_death_excluded <- sum(cohort$PATIENT_ID %in% DEATH_DATE_EXCLUDE_IDS)
cohort <- cohort %>% filter(!PATIENT_ID %in% DEATH_DATE_EXCLUDE_IDS)
cat(sprintf("  Excluded (unrepairable MBSF death date): %d\n", n_death_excluded))
cat(sprintf("  Cohort after enrollment check:    %d\n", nrow(cohort)))


cat("Deriving complete follow-up windows...\n")

# FFS coverage over months k_from..k_to from the diagnosis month, capped at death and LAST_FUP
coverage_complete <- function(cohort_df, mbsf_months, k_from, k_to, prefix) {
  cohort_df %>%
    select(PATIENT_ID, dx_date_end, death_date) %>%
    tidyr::crossing(k = k_from:k_to) %>%
    mutate(
      month_start = add_with_rollback(floor_date(dx_date_end, "month"), months(k)),
      month_end   = ceiling_date(month_start, "month") - days(1)
    ) %>%
    filter(is.na(death_date) | month_start <= death_date, month_end <= LAST_FUP) %>%
    left_join(mbsf_months %>% select(PATIENT_ID, date, ffs_month),
              by = c("PATIENT_ID", "month_end" = "date")) %>%
    mutate(ffs_month = coalesce(ffs_month, FALSE)) %>%
    group_by(PATIENT_ID) %>%
    summarize(exp_m = n(), obs_m = sum(ffs_month), .groups = "drop") %>%
    transmute(
      PATIENT_ID,
      "{prefix}_expected_months" := as.integer(exp_m),
      "{prefix}_ffs_months"      := as.integer(obs_m),
      "{prefix}_complete"        := as.integer(obs_m == exp_m),
      "{prefix}_observable"      := as.integer(exp_m > 0)
    )
}

cohort <- cohort %>%
  # img_fup: 11 months before to 3 after the diagnosis month; tx_fup: diagnosis month and 11 after
  left_join(coverage_complete(cohort, mbsf_long, -11L,  3L, "img_fup"), by = "PATIENT_ID") %>%
  left_join(coverage_complete(cohort, mbsf_long,   0L, 11L, "tx_fup"),  by = "PATIENT_ID") %>%
  mutate(across(c(img_fup_expected_months, img_fup_ffs_months,
                  img_fup_complete, img_fup_observable,
                  tx_fup_expected_months,  tx_fup_ffs_months,
                  tx_fup_complete,  tx_fup_observable),
                ~ as.integer(replace_na(.x, 0L))))

cohort <- cohort %>%
  mutate(
    img_fup_analytic = as.integer(img_fup_complete == 1L & img_fup_observable == 1L),
    tx_fup_analytic  = as.integer(tx_fup_complete  == 1L & tx_fup_observable  == 1L)
  )

cat(sprintf("  Complete coverage, imaging window (12 mo pre through 3 mo post, capped at death): %d of %d (%.1f%%)\n",
            sum(cohort$img_fup_analytic), nrow(cohort),
            100 * mean(cohort$img_fup_analytic)))
cat(sprintf("  Complete coverage, therapy window (dx month + 11, capped at death):               %d of %d (%.1f%%)\n",
            sum(cohort$tx_fup_analytic), nrow(cohort),
            100 * mean(cohort$tx_fup_analytic)))
cat(sprintf("  For contrast, the old all-12-calendar-months rule:                                %d of %d (%.1f%%)\n",
            sum(cohort$post_months == 12), nrow(cohort),
            100 * mean(cohort$post_months == 12)))
cat(sprintf("  Patients with no observable month in the therapy window (excluded):               %d\n",
            sum(cohort$tx_fup_observable == 0L)))
cat(sprintf("  Death date recorded before diagnosis (data error, excluded above):                %d\n",
            sum(!is.na(cohort$death_date) & cohort$death_date < cohort$dx_date_first)))


# =============================================================================
# V. SEER variable recoding
# =============================================================================
cat("Cleaning SEER variables...\n")

cohort <- cohort %>%
  rename(
    sex            = SEX,
    race           = RACE_RECODE_WHITE_BLACK_OTHER,
    marital_status = MARITAL_STATUS_AT_DIAGNOSIS
  ) %>%
  mutate(
    sex = case_when(
      as.integer(sex) == 1 ~ "Male",
      as.integer(sex) == 2 ~ "Female",
      TRUE                 ~ "Unknown"
    ),
    race = case_when(
      as.integer(race) == 1 ~ "White",
      as.integer(race) == 2 ~ "Black",
      TRUE                  ~ "Other or unknown"
    ),
    marital_status = case_when(
      marital_status == "1"                  ~ "Never married",
      marital_status == "2"                  ~ "Married",
      marital_status %in% c("3","4","5","6") ~ "Other",
      TRUE                                   ~ "Unknown"
    ),
    microscopically_confirmed = as.integer(DIAGNOSTIC_CONFIRMATION) %in% c(1, 2, 3, 4),
    histology_code = sprintf("%s/%s", HISTOLOGIC_TYPE_ICD_O_3, BEHAVIOR_CODE_ICD_O_3),

    ajcc_t7 = suppressWarnings(as.numeric(DERIVED_AJCC_T_7TH_ED_2010_2015)),
    ajcc_n7 = suppressWarnings(as.numeric(DERIVED_AJCC_N_7TH_ED_2010_2015)),
    ajcc_m7 = suppressWarnings(as.numeric(DERIVED_AJCC_M_7TH_ED_2010_2015)),

    # T/N/M from the era field: EOD 2018 strings, 2016-2017 c/p strings, or AJCC 7th numeric codes
    tnm_t = case_when(
      str_detect(coalesce(DERIVED_EOD_2018_T_2018, ""),          "^T1") |
        str_detect(coalesce(DERIVEDSEERCOMBINED_T_2016_2017, ""), "^[cp]1") |
        between(ajcc_t7, 100, 199) | between(ajcc_t7, 800, 810)    ~ "1",
      str_detect(coalesce(DERIVED_EOD_2018_T_2018, ""),          "^T2") |
        str_detect(coalesce(DERIVEDSEERCOMBINED_T_2016_2017, ""), "^[cp]2") |
        between(ajcc_t7, 200, 299)                                 ~ "2",
      str_detect(coalesce(DERIVED_EOD_2018_T_2018, ""),          "^T3") |
        str_detect(coalesce(DERIVEDSEERCOMBINED_T_2016_2017, ""), "^[cp]3") |
        between(ajcc_t7, 300, 399)                                 ~ "3",
      str_detect(coalesce(DERIVED_EOD_2018_T_2018, ""),          "^T4") |
        str_detect(coalesce(DERIVEDSEERCOMBINED_T_2016_2017, ""), "^[cp]4") |
        between(ajcc_t7, 400, 499)                                 ~ "4",
      str_detect(coalesce(DERIVED_EOD_2018_T_2018, ""),          "^Tis") |
        str_detect(coalesce(DERIVEDSEERCOMBINED_T_2016_2017, ""), "^[cp]IS") |
        ajcc_t7 == 50                                              ~ "is",
      str_detect(coalesce(DERIVED_EOD_2018_T_2018, ""),          "^T0") |
        str_detect(coalesce(DERIVEDSEERCOMBINED_T_2016_2017, ""), "^[cp]0") |
        ajcc_t7 == 0                                               ~ "0",
      str_detect(coalesce(DERIVED_EOD_2018_T_2018, ""),          "^TX") |
        str_detect(coalesce(DERIVEDSEERCOMBINED_T_2016_2017, ""), "^[cp]X") |
        ajcc_t7 == 999                                             ~ "X",
      TRUE ~ NA_character_
    ),
    tnm_n = case_when(
      str_detect(coalesce(DERIVED_EOD_2018_N_2018, ""),          "^N0") |
        str_detect(coalesce(DERIVEDSEERCOMBINED_N_2016_2017, ""), "^[cp]0") |
        between(ajcc_n7, 0, 40)                                    ~ "0",
      str_detect(coalesce(DERIVED_EOD_2018_N_2018, ""),          "^N1") |
        str_detect(coalesce(DERIVEDSEERCOMBINED_N_2016_2017, ""), "^[cp]1") |
        between(ajcc_n7, 100, 199)                                 ~ "1",
      str_detect(coalesce(DERIVED_EOD_2018_N_2018, ""),          "^N2") |
        str_detect(coalesce(DERIVEDSEERCOMBINED_N_2016_2017, ""), "^[cp]2") |
        between(ajcc_n7, 200, 299)                                 ~ "2",
      str_detect(coalesce(DERIVED_EOD_2018_N_2018, ""),          "^N3") |
        str_detect(coalesce(DERIVEDSEERCOMBINED_N_2016_2017, ""), "^[cp]3") |
        between(ajcc_n7, 300, 399)                                 ~ "3",
      str_detect(coalesce(DERIVED_EOD_2018_N_2018, ""),          "^NX") |
        str_detect(coalesce(DERIVEDSEERCOMBINED_N_2016_2017, ""), "^[cp]X") |
        ajcc_n7 == 999                                             ~ "X",
      TRUE ~ NA_character_
    ),
    tnm_m = case_when(
      str_detect(coalesce(DERIVED_EOD_2018_M_2018, ""),          "^M0") |
        str_detect(coalesce(DERIVEDSEERCOMBINED_M_2016_2017, ""), "^[cp]0") |
        between(ajcc_m7, 0, 10)                                    ~ "0",
      str_detect(coalesce(DERIVED_EOD_2018_M_2018, ""),          "^M1") |
        str_detect(coalesce(DERIVEDSEERCOMBINED_M_2016_2017, ""), "^[cp]1") |
        between(ajcc_m7, 100, 199)                                 ~ "1",
      str_detect(coalesce(DERIVED_EOD_2018_M_2018, ""),          "^MX") |
        str_detect(coalesce(DERIVEDSEERCOMBINED_M_2016_2017, ""), "^[cp]X") |
        ajcc_m7 == 999                                             ~ "X",
      TRUE ~ NA_character_
    ),

    # CS size codes 991-995 mean under 1-5 cm; 2016+ code 990 is a microscopic focus
    cs_size_num = suppressWarnings(as.numeric(CS_TUMOR_SIZE_2004_2015)),
    size_lt2015 = case_when(
      between(cs_size_num, 1, 989) ~ cs_size_num / 10,
      cs_size_num == 0             ~ 0,
      cs_size_num == 991           ~ 0.99,
      cs_size_num == 992           ~ 1.99,
      cs_size_num == 993           ~ 2.99,
      cs_size_num == 994           ~ 3.99,
      cs_size_num == 995           ~ 4.99,
      TRUE                         ~ NA_real_
    ),
    ts_size_num = suppressWarnings(as.numeric(TUMOR_SIZE_SUMMARY_2016)),
    size_gt2015 = case_when(
      ts_size_num > 0 & ts_size_num <= 989 ~ ts_size_num / 10,
      ts_size_num == 990                   ~ 0,
      TRUE                                 ~ NA_real_
    ),
    size = if_else(!is.na(size_lt2015), size_lt2015, size_gt2015),

    # AJCC 8th T from size and the era T code; the first matching branch wins
    t_stage_8 = case_when(
      size > 0   & size < 1 & tnm_t == "1"                ~ "T1a",
      size >= 1  & size < 2 & tnm_t == "1"                ~ "T1b",
      size >= 2  & size < 3 & tnm_t == "1"                ~ "T1c",
      (size >= 3 & size < 4) | (tnm_t == "2" & size < 4)  ~ "T2a",
      size >= 4  & size < 5                                ~ "T2b",
      (size >= 5 & size < 7) | (tnm_t == "3" & size < 7)  ~ "T3",
      (size >= 7 & size < 100) | tnm_t == "4"             ~ "T4",
      TRUE                                                 ~ NA_character_
    ),

    stage_simple = case_when(
      t_stage_8 %in% c("T1a","T1b","T1c") & tnm_n == "0" & tnm_m == "0" ~ "Stage IA",
      t_stage_8 == "T2a"                  & tnm_n == "0" & tnm_m == "0" ~ "Stage IB",
      t_stage_8 == "T2b"                  & tnm_n == "0" & tnm_m == "0" ~ "Stage IIA",
      (t_stage_8 %in% c("T1a","T1b","T1c","T2a","T2b") & tnm_n == "1" & tnm_m == "0") |
        (t_stage_8 == "T3" & tnm_n == "0" & tnm_m == "0")               ~ "Stage IIB",
      (t_stage_8 %in% c("T1a","T1b","T1c","T2a","T2b") & tnm_n == "2" & tnm_m == "0") |
        (t_stage_8 %in% c("T3","T4") & tnm_n == "1" & tnm_m == "0") |
        (t_stage_8 == "T4" & tnm_n == "0" & tnm_m == "0")               ~ "Stage IIIA",
      (t_stage_8 %in% c("T1a","T1b","T1c","T2a","T2b") & tnm_n == "3" & tnm_m == "0") |
        (t_stage_8 %in% c("T3","T4") & tnm_n == "2" & tnm_m == "0")     ~ "Stage IIIB",
      t_stage_8 %in% c("T3","T4") & tnm_n == "3" & tnm_m == "0"         ~ "Stage IIIC",
      tnm_m == "1"                                                        ~ "Stage IV",
      !is.na(t_stage_8) | !is.na(tnm_n) | !is.na(tnm_m)                 ~ "Data available but not classified",
      TRUE                                                                ~ "Missing Data"
    ),

    # summary stage 1/2/7 recoded 1/2/3, all else 9; fallback stage tier for the concordance rules
    staging = case_when(
      as.integer(COMBINED_SUMMARY_STAGE_2004) == 1 ~ 1L,
      as.integer(COMBINED_SUMMARY_STAGE_2004) == 2 ~ 2L,
      as.integer(COMBINED_SUMMARY_STAGE_2004) == 7 ~ 3L,
      TRUE                                         ~ 9L
    ),

    rural_urban = case_when(
      RURAL_URBAN_CONTINUUM_CODE %in% c("01","02","03") ~ "Urban",
      RURAL_URBAN_CONTINUUM_CODE %in% c("11","12")      ~ "Rural",
      TRUE                                              ~ "Unknown"
    ),

    # non-positive follow-up floored at 0.01 months
    tt_months_trunc = case_when(
      tt_months >  0 ~ tt_months,
      tt_months <= 0 ~ 0.01,
      TRUE           ~ NA_real_
    )
  ) %>%
  select(-cs_size_num, -ts_size_num, -ajcc_t7, -ajcc_n7, -ajcc_m7)


# registry-reported first-course treatment, kept beside the claims-based flags
cat("Recoding SEER treatment variables...\n")

cohort <- cohort %>%
  mutate(
    chemo_code = str_trim(as.character(CHEMOTHERAPY_RECODE_YES_NO_UNK)),
    chemotherapy = case_when(
      chemo_code == "1" ~ "Yes",
      chemo_code == "0" ~ "No or unknown",
      TRUE              ~ "No or unknown"
    ),
    chemotherapy_yes = as.integer(chemotherapy == "Yes"),

    rad_code = str_trim(as.character(RADIATION_RECODE)),
    radiation = case_when(
      rad_code == "0" ~ "None/unknown",
      rad_code == "1" ~ "Beam radiation",
      rad_code == "2" ~ "Radioactive implants",
      rad_code == "3" ~ "Radioisotopes",
      rad_code == "4" ~ "Combination of beam with implants/isotopes",
      rad_code == "5" ~ "Radiation, NOS",
      rad_code == "6" ~ "Other radiation",
      rad_code == "7" ~ "Refused",
      rad_code == "8" ~ "Recommended, unknown if administered",
      rad_code == "9" ~ "Unknown",
      TRUE            ~ "Not reported by registry"
    ),
    radiation_reported = as.integer(radiation != "Not reported by registry"),
    radiation_yes_no_unk = case_when(
      rad_code %in% c("1","2","3","4","5","6")  ~ "Yes",
      rad_code %in% c("0","7","8","9")          ~ "No or unknown",
      TRUE                                       ~ "Not reported by registry"
    ),
    radiation_yes = case_when(
      radiation_yes_no_unk == "Yes"           ~ 1L,
      radiation_yes_no_unk == "No or unknown" ~ 0L,
      TRUE                                    ~ NA_integer_
    )
  ) %>%
  select(-chemo_code, -rad_code)


cat("Deriving SEER stage at diagnosis...\n")

# strips the c/p prefix and the T/N letter; 88/888 (not applicable) and blanks become Unknown
seer_tn_group <- function(x, letter) {
  v <- str_to_upper(str_trim(as.character(x)))
  v <- str_remove(v, "^[CP]")
  v <- str_remove(v, paste0("^", letter))
  n <- suppressWarnings(as.numeric(v))
  ajcc7 <- str_detect(v, "^[0-9]{3}$")

  grp <- case_when(
    is.na(v) | v == "" | v %in% c("88", "888")               ~ NA_character_,
    ajcc7 & n == 999                                         ~ "X",
    ajcc7 & n == 0                                           ~ "0",
    ajcc7 & n == 50                                          ~ "is",
    ajcc7 & (between(n, 100, 199) | between(n, 800, 810))    ~ "1",
    ajcc7 & between(n, 200, 299)                             ~ "2",
    ajcc7 & between(n, 300, 399)                             ~ "3",
    ajcc7 & between(n, 400, 499)                             ~ "4",
    ajcc7                                                    ~ NA_character_,
    str_starts(v, "X")                                       ~ "X",
    str_starts(v, "IS")                                      ~ "is",
    str_starts(v, "0")                                       ~ "0",
    str_starts(v, "1")                                       ~ "1",
    str_starts(v, "2")                                       ~ "2",
    str_starts(v, "3")                                       ~ "3",
    str_starts(v, "4")                                       ~ "4",
    TRUE                                                     ~ NA_character_
  )
  if_else(is.na(grp), "Unknown/not available", paste0(letter, grp))
}

era_value <- function(yr, v2015, v1617, v2018) {
  raw <- case_when(
    yr >= 2018            ~ as.character(v2018),
    yr %in% c(2016, 2017) ~ as.character(v1617),
    TRUE                  ~ as.character(v2015)
  )
  na_if(str_trim(raw), "")
}

cohort <- cohort %>%
  mutate(
    seer_t = era_value(yr, DERIVED_AJCC_T_7TH_ED_2010_2015,
                       DERIVEDSEERCOMBINED_T_2016_2017, DERIVED_EOD_2018_T_2018),
    seer_n = era_value(yr, DERIVED_AJCC_N_7TH_ED_2010_2015,
                       DERIVEDSEERCOMBINED_N_2016_2017, DERIVED_EOD_2018_N_2018),
    seer_m = era_value(yr, DERIVED_AJCC_M_7TH_ED_2010_2015,
                       DERIVEDSEERCOMBINED_M_2016_2017, DERIVED_EOD_2018_M_2018),
    seer_stage_group = era_value(yr, DERIVEDAJCCSTAGEGROUP7THED20102,
                                 DERIVEDSEERCMBSTG_GRP_2016_2017,
                                 DERIVEDEOD2018_STAGE_GROUP_2018),

    seer_stage_source = case_when(
      yr >= 2018            ~ "Derived EOD 2018 (2018+)",
      yr %in% c(2016, 2017) ~ "Derived SEER Combined (2016-2017)",
      TRUE                  ~ "Derived AJCC 7th ed (2015)"
    ),

    seer_t_group = seer_tn_group(seer_t, "T"),
    seer_n_group = seer_tn_group(seer_n, "N"),

    # registries with no TNM fields in this delivery; 02 uses 61-66 for missing_ajcc
    tnm_reported_registry = as.integer(
      !SEER_REGISTRY %in% c("61", "62", "63", "65", "66")
    )
  )


# =============================================================================
# VI. Claims: screening, peri-diagnosis imaging and treatment
# =============================================================================

cat("Identifying LCS screening...\n")

# every HCPCS line in any code list, read once
proc_claims <- open_dataset(file.path(EXTRACTED, "lung_procedures.parquet")) %>%
  select(PATIENT_ID, HCPCS_CD, CLM_THRU_DT) %>%
  filter(HCPCS_CD %in% PROC_CPTS) %>%
  collect() %>%
  filter(PATIENT_ID %in% pids) %>%
  mutate(proc_date = claim_date_of(CLM_THRU_DT)) %>%
  filter(!is.na(proc_date)) %>%
  select(PATIENT_ID, HCPCS_CD, proc_date)

lcs_claims <- proc_claims %>%
  filter(HCPCS_CD %in% LCS_CPTS) %>%
  select(PATIENT_ID, lcs_date = proc_date)

lcs_pre <- cohort %>%
  select(PATIENT_ID, dx_date_end) %>%
  left_join(lcs_claims, by = "PATIENT_ID") %>%
  filter(!is.na(lcs_date), as.numeric(dx_date_end - lcs_date) >= 0) %>%
  group_by(PATIENT_ID) %>%
  arrange(lcs_date, .by_group = TRUE) %>%
  mutate(
    scan_number = row_number(),
    last_scan   = last(lcs_date),
    total_scans = n()
  ) %>%
  ungroup()

lcs_summary <- lcs_pre %>%
  distinct(PATIENT_ID, .keep_all = TRUE) %>%
  select(PATIENT_ID, last_scan, total_scans)

cohort <- cohort %>%
  left_join(lcs_summary, by = "PATIENT_ID") %>%
  mutate(
    # Screened = last pre-diagnosis LDCT within 365 days of dx_date_end
    time_diff = as.numeric(dx_date_end - last_scan),
    scan_type = case_when(
      time_diff <  365                     ~ "Screened",
      time_diff >= 365 & !is.na(time_diff) ~ "Missed Scan",
      is.na(time_diff)                     ~ "Never Screened",
      TRUE                                 ~ "Other"
    ),
    screened = as.integer(scan_type == "Screened")
  )

cat("Screening classification:\n")
print(cohort %>% count(scan_type))

lcs_claims_primary <- proc_claims %>%
  filter(HCPCS_CD %in% LCS_CPTS_PRIMARY) %>%
  select(PATIENT_ID, lcs_date = proc_date)

lcs_summary_primary <- cohort %>%
  select(PATIENT_ID, dx_date_end) %>%
  left_join(lcs_claims_primary, by = "PATIENT_ID") %>%
  filter(!is.na(lcs_date), as.numeric(dx_date_end - lcs_date) >= 0) %>%
  group_by(PATIENT_ID) %>%
  summarize(last_scan_primary = max(lcs_date), .groups = "drop")

cohort <- cohort %>%
  left_join(lcs_summary_primary, by = "PATIENT_ID") %>%
  mutate(
    screened_primary_only = as.integer(
      !is.na(last_scan_primary) &
        as.numeric(dx_date_end - last_scan_primary) < 365
    )
  )

cat(sprintf(
  "  S8032-only screening claims: %d patients screened by the full definition vs. %d by G0297/71271 alone (%+d)\n",
  sum(cohort$screened), sum(cohort$screened_primary_only),
  sum(cohort$screened) - sum(cohort$screened_primary_only)))


cat("Deriving peri-diagnosis brain MRI and PET...\n")

# CPT claims from IMG_PRE_DAYS before to IMG_POST_DAYS after dx_date_end
imaging_win <- cohort %>%
  select(PATIENT_ID, dx_date_end) %>%
  inner_join(
    proc_claims %>% filter(HCPCS_CD %in% c(BRAIN_MRI_CPTS, PET_CPTS)),
    by = "PATIENT_ID"
  ) %>%
  mutate(
    modality     = if_else(HCPCS_CD %in% BRAIN_MRI_CPTS, "brain_mri", "pet"),
    days_from_dx = as.numeric(proc_date - dx_date_end)
  ) %>%
  filter(days_from_dx >= -IMG_PRE_DAYS, days_from_dx <= IMG_POST_DAYS)

summarize_modality <- function(mod) {
  imaging_win %>%
    filter(modality == mod) %>%
    group_by(PATIENT_ID) %>%
    summarize(
      n_scans      = n(),
      n_pre        = sum(days_from_dx <  0),
      n_post       = sum(days_from_dx >= 0),
      first_date   = min(proc_date),
      last_date    = max(proc_date),
      days_from_dx = min(days_from_dx),
      .groups = "drop"
    ) %>%
    setNames(c("PATIENT_ID",
               paste0(mod, "_pre12post3_count"),
               paste0(mod, "_n_pre"),
               paste0(mod, "_n_post"),
               paste0(mod, "_first_date"),
               paste0(mod, "_last_date"),
               paste0(mod, "_days_from_dx")))
}

brain_mri_prior <- imaging_win %>%
  filter(HCPCS_CD %in% BRAIN_MRI_CPTS_PRIOR) %>%
  distinct(PATIENT_ID) %>%
  mutate(brain_mri_pre12post3_prior = 1L)

cohort <- cohort %>%
  left_join(summarize_modality("brain_mri"), by = "PATIENT_ID") %>%
  left_join(summarize_modality("pet"),       by = "PATIENT_ID") %>%
  left_join(brain_mri_prior,                 by = "PATIENT_ID") %>%
  mutate(
    across(c(brain_mri_pre12post3_count, brain_mri_n_pre, brain_mri_n_post,
             pet_pre12post3_count,       pet_n_pre,       pet_n_post),
           ~ as.integer(replace_na(.x, 0L))),
    brain_mri_pre12post3       = as.integer(brain_mri_pre12post3_count >= 1L),
    brain_mri_pre12post3_prior = as.integer(replace_na(brain_mri_pre12post3_prior, 0L)),
    pet_pre12post3             = as.integer(pet_pre12post3_count       >= 1L)
  )

# console diagnostics
cat("Peri-diagnosis imaging (12 months before to 3 months after diagnosis):\n")
print(
  tibble(
    modality = c("Brain MRI", "PET / PET-CT"),
    n        = c(sum(cohort$brain_mri_pre12post3), sum(cohort$pet_pre12post3)),
    n_pre_dx = c(sum(cohort$brain_mri_n_pre  >= 1L), sum(cohort$pet_n_pre  >= 1L)),
    n_post_dx = c(sum(cohort$brain_mri_n_post >= 1L), sum(cohort$pet_n_post >= 1L))
  ) %>%
    mutate(pct = round(100 * n / nrow(cohort), 1))
)
cat(sprintf("  Brain MRI incl. 70552: %d (%.1f%%) | prior study's 70551+70553 only: %d (%.1f%%) | difference: %d\n",
            sum(cohort$brain_mri_pre12post3),
            100 * mean(cohort$brain_mri_pre12post3),
            sum(cohort$brain_mri_pre12post3_prior),
            100 * mean(cohort$brain_mri_pre12post3_prior),
            sum(cohort$brain_mri_pre12post3) - sum(cohort$brain_mri_pre12post3_prior)))

img_ok <- cohort %>% filter(tx_fup_analytic == 1L)
cat(sprintf("\n  Restricted to complete 12-month post-dx coverage: n = %d (%.1f%% of cohort)\n",
            nrow(img_ok), 100 * nrow(img_ok) / nrow(cohort)))
print(
  tibble(
    modality  = c("Brain MRI", "PET / PET-CT"),
    n         = c(sum(img_ok$brain_mri_pre12post3), sum(img_ok$pet_pre12post3)),
    n_pre_dx  = c(sum(img_ok$brain_mri_n_pre  >= 1L), sum(img_ok$pet_n_pre  >= 1L)),
    n_post_dx = c(sum(img_ok$brain_mri_n_post >= 1L), sum(img_ok$pet_n_post >= 1L))
  ) %>%
    mutate(pct = round(100 * n / nrow(img_ok), 1))
)
cat("\n  Brain MRI and PET by screening status, complete-coverage group only:\n")
print(img_ok %>%
  group_by(scan_type) %>%
  summarize(n = n(),
            pct_brain_mri = round(100 * mean(brain_mri_pre12post3 == 1L), 1),
            pct_pet       = round(100 * mean(pet_pre12post3       == 1L), 1),
            .groups = "drop"))
cat("  (chemotherapy and immunotherapy in the same subgroup are reported below,\n")
cat("   after those variables are derived)\n")


cat("Deriving post-diagnosis systemic therapy...\n")

# treatment window: day 0 through day 365 from dx_date_first
systemic_post <- cohort %>%
  select(PATIENT_ID, dx_date_first) %>%
  inner_join(
    proc_claims %>%
      filter(HCPCS_CD %in% c(CHEMO_CPTS, CHEMO_ADMIN_RANGE, CHEMO_DRUG_RANGE,
                             IMMUNO_HCPCS_ALL, RADIATION_CPT, RESECTION_CPT)),
    by = "PATIENT_ID"
  ) %>%
  mutate(days_after = as.numeric(proc_date - dx_date_first)) %>%
  filter(days_after >= 0, days_after <= 365)

summarize_therapy <- function(codes, prefix) {
  systemic_post %>%
    filter(HCPCS_CD %in% codes) %>%
    group_by(PATIENT_ID) %>%
    summarize(
      n_claims   = n(),
      first_date = min(proc_date),
      days_after = min(days_after),
      .groups = "drop"
    ) %>%
    setNames(c("PATIENT_ID",
               paste0(prefix, "_count"),
               paste0(prefix, "_first_date"),
               paste0(prefix, "_days_after")))
}

cohort <- cohort %>%
  # HCPCS-only flags; radiation_post12 is widened with the MEDPAR sources below
  left_join(summarize_therapy(c(CHEMO_ADMIN_RANGE, CHEMO_DRUG_RANGE),
                              "chemo_post12_hcpcs"),                by = "PATIENT_ID") %>%
  left_join(summarize_therapy(CHEMO_CPTS,   "chemo_post12_enum"),   by = "PATIENT_ID") %>%
  left_join(summarize_therapy(CHEMO_ADMIN_RANGE, "chemo_post12_admin"), by = "PATIENT_ID") %>%
  left_join(summarize_therapy(CHEMO_DRUG_RANGE,  "chemo_post12_drug"),  by = "PATIENT_ID") %>%
  left_join(summarize_therapy(IMMUNO_HCPCS_SCLC, "immuno_post12"),  by = "PATIENT_ID") %>%
  left_join(summarize_therapy(RADIATION_DELIVERY_CPT, "radiation_post12"), by = "PATIENT_ID") %>%
  left_join(summarize_therapy(RADIATION_CPT, "radiation_post12_any"), by = "PATIENT_ID") %>%
  left_join(summarize_therapy(RESECTION_CPT, "resection_post12_hcpcs"), by = "PATIENT_ID") %>%
  mutate(
    across(c(chemo_post12_hcpcs_count, chemo_post12_enum_count,
             chemo_post12_admin_count, chemo_post12_drug_count,
             immuno_post12_count,
             radiation_post12_count, radiation_post12_any_count,
             resection_post12_hcpcs_count),
           ~ as.integer(replace_na(.x, 0L))),
    chemo_post12_admin      = as.integer(chemo_post12_admin_count      >= 1L),
    chemo_post12_drug       = as.integer(chemo_post12_drug_count       >= 1L),
    chemo_post12_enumerated = as.integer(chemo_post12_enum_count       >= 1L),
    immuno_post12 = as.integer(immuno_post12_count >= 1L),
    radiation_post12        = as.integer(radiation_post12_count        >= 1L),
    radiation_post12_any     = as.integer(radiation_post12_any_count    >= 1L),
    resection_post12_hcpcs  = as.integer(resection_post12_hcpcs_count  >= 1L)
  )


cat("  Reading MEDPAR inpatient procedure/diagnosis codes...\n")

# MEDPAR procedure and diagnosis columns melted once; the claim date is the discharge date
medpar_px <- tibble(PATIENT_ID = character(0), claim_date = as.Date(character(0)),
                    days_after = numeric(0), days_from_dx = numeric(0),
                    px_code = character(0))
medpar_dx <- tibble(PATIENT_ID = character(0), claim_date = as.Date(character(0)),
                    days_after = numeric(0), dx_code = character(0))

if (file.exists(medpar_path <- file.path(EXTRACTED, "medpar.parquet"))) {
  mp_cols    <- names(open_dataset(medpar_path))
  px_cols    <- mp_cols[str_detect(mp_cols, "^SRGCL_PRCDR_[0-9]+_CD$")]
  mp_dx_cols <- mp_cols[str_detect(mp_cols, "^DGNS_[0-9]+_CD$")]

  medpar_base <- open_dataset(medpar_path) %>%
    select(PATIENT_ID, ADMSN_DT, DSCHRG_DT, all_of(c(px_cols, mp_dx_cols))) %>%
    filter(PATIENT_ID %in% pids) %>%
    collect() %>%
    left_join(cohort %>% select(PATIENT_ID, dx_date_first, dx_date_end),
              by = "PATIENT_ID") %>%
    mutate(claim_date   = ymd(DSCHRG_DT),
           days_after   = as.numeric(claim_date - dx_date_first),
           days_from_dx = as.numeric(claim_date - dx_date_end)) %>%
    filter(!is.na(claim_date))

  if (length(px_cols) > 0) {
    medpar_px <- medpar_base %>%
      select(PATIENT_ID, claim_date, days_after, days_from_dx, all_of(px_cols)) %>%
      pivot_longer(all_of(px_cols), names_to = NULL, values_to = "px_code") %>%
      mutate(px_code = str_to_upper(str_remove_all(str_trim(px_code), "\\."))) %>%
      filter(!is.na(px_code), px_code != "")
  }
  if (length(mp_dx_cols) > 0) {
    medpar_dx <- medpar_base %>%
      select(PATIENT_ID, claim_date, days_after, all_of(mp_dx_cols)) %>%
      pivot_longer(all_of(mp_dx_cols), names_to = NULL, values_to = "dx_code") %>%
      mutate(dx_code = str_to_upper(str_remove_all(str_trim(dx_code), "\\."))) %>%
      filter(!is.na(dx_code), dx_code != "")
  }
  rm(medpar_base)
}

first_claim_date <- function(df, out_name) {
  df %>%
    group_by(PATIENT_ID) %>%
    summarize(!!out_name := min(claim_date), .groups = "drop")
}

is_radiation_pcs <- function(px) {
  str_detect(px, RADIATION_PX_ICD10_SECTION) &
    !(str_sub(px, 3, 3) == "Y" & str_sub(px, 5, 5) %in% RADIATION_PX_ICD10_EXCL_MQ)
}


cat("  Adding ICD-based chemotherapy (Z51.11 / V58.11, inpatient administration)...\n")

chemo_icd_dx <- tibble(PATIENT_ID = character(0), chemo_icd_dx_date = as.Date(character(0)))
if (file.exists(dx_path_chemo <- file.path(EXTRACTED, "lung_diagnosis.parquet"))) {
  dx_cols <- names(open_dataset(dx_path_chemo))
  dx_code_cols_all <- dx_cols[str_detect(dx_cols, "ICD_DGNS_CD|PRNCPAL_DGNS")]

  chemo_icd_dx <- open_dataset(dx_path_chemo) %>%
    select(PATIENT_ID, CLM_THRU_DT, all_of(dx_code_cols_all)) %>%
    filter(PATIENT_ID %in% pids) %>%
    collect() %>%
    left_join(cohort %>% select(PATIENT_ID, dx_date_first), by = "PATIENT_ID") %>%
    mutate(claim_date = claim_date_of(CLM_THRU_DT),
           days_after = as.numeric(claim_date - dx_date_first)) %>%
    filter(days_after >= 0, days_after <= 365) %>%
    pivot_longer(all_of(dx_code_cols_all), names_to = NULL, values_to = "dx_code") %>%
    mutate(dx_code = str_to_upper(str_remove_all(str_trim(dx_code), "\\."))) %>%
    filter(dx_code %in% c(CHEMO_DX_ICD10, CHEMO_DX_ICD9)) %>%
    group_by(PATIENT_ID) %>%
    summarize(chemo_icd_dx_date = min(claim_date), .groups = "drop")
}

medpar_px_tx <- medpar_px %>% filter(days_after >= 0, days_after <= 365)
medpar_dx_tx <- medpar_dx %>% filter(days_after >= 0, days_after <= 365)

chemo_icd_px <- medpar_px_tx %>%
  filter(px_code %in% c(CHEMO_PX_ICD9, CHEMO_PX_ICD10_NEWTECH) |
           str_detect(px_code, CHEMO_PX_ICD10_PATTERN)) %>%
  first_claim_date("chemo_icd_px_date")

cat("  Adding ICD-based lung resection (MEDPAR inpatient procedure codes)...\n")

resection_px_pattern <- paste0("^(", paste(RESECTION_PX_ICD9_PREFIX, collapse = "|"), ")")
resection_icd_px <- medpar_px_tx %>%
  filter(px_code %in% RESECTION_PX_ICD10 |
           str_detect(px_code, resection_px_pattern)) %>%
  first_claim_date("resection_icd_px_date")

# anatomic subset derived here while medpar_px_tx is in memory; descriptive only
resection_anatomic_px_pattern <- paste0(
  "^(", paste(RESECTION_ANATOMIC_PX_ICD9_PREFIX, collapse = "|"), ")")
resection_anatomic_icd_px <- medpar_px_tx %>%
  filter(px_code %in% RESECTION_ANATOMIC_PX_ICD10 |
           str_detect(px_code, resection_anatomic_px_pattern)) %>%
  first_claim_date("resection_anatomic_icd_px_date")

cat("  Adding ICD-based radiation therapy (MEDPAR inpatient procedure codes)...\n")

radiation_px_pattern <- paste0("^(", paste(RADIATION_PX_ICD9_PREFIX, collapse = "|"), ")")
radiation_icd_px <- medpar_px_tx %>%
  filter(is_radiation_pcs(px_code) | str_detect(px_code, radiation_px_pattern)) %>%
  first_claim_date("radiation_icd_px_date")

radiation_icd_dx <- tibble(PATIENT_ID = character(0),
                           radiation_icd_dx_date = as.Date(character(0)))
if (file.exists(dx_path_chemo)) {
  radiation_icd_dx <- open_dataset(dx_path_chemo) %>%
    select(PATIENT_ID, CLM_THRU_DT, all_of(dx_code_cols_all)) %>%
    filter(PATIENT_ID %in% pids) %>%
    collect() %>%
    left_join(cohort %>% select(PATIENT_ID, dx_date_first), by = "PATIENT_ID") %>%
    mutate(claim_date = claim_date_of(CLM_THRU_DT),
           days_after = as.numeric(claim_date - dx_date_first)) %>%
    filter(days_after >= 0, days_after <= 365) %>%
    pivot_longer(all_of(dx_code_cols_all), names_to = NULL, values_to = "dx_code") %>%
    mutate(dx_code = str_to_upper(str_remove_all(str_trim(dx_code), "\\."))) %>%
    filter(dx_code %in% c(RADIATION_DX_ICD10, RADIATION_DX_ICD9)) %>%
    group_by(PATIENT_ID) %>%
    summarize(radiation_icd_dx_date = min(claim_date), .groups = "drop")
}
radiation_icd_dx <- bind_rows(
    radiation_icd_dx,
    medpar_dx_tx %>%
      filter(dx_code %in% c(RADIATION_DX_ICD10, RADIATION_DX_ICD9)) %>%
      first_claim_date("radiation_icd_dx_date")
  ) %>%
  group_by(PATIENT_ID) %>%
  summarize(radiation_icd_dx_date = min(radiation_icd_dx_date), .groups = "drop")

cat("  Adding ICD-based immunotherapy (MEDPAR, class-level, sensitivity only)...\n")

immuno_icd_px <- medpar_px_tx %>%
  filter(str_detect(px_code, IMMUNO_PX_ICD10_PATTERN)) %>%
  first_claim_date("immuno_icd_px_date")

cat("  Adding ICD-based brain MRI / PET (MEDPAR, component columns only)...\n")

medpar_px_img <- medpar_px %>%
  filter(days_from_dx >= -IMG_PRE_DAYS, days_from_dx <= IMG_POST_DAYS)

brain_mri_icd_px <- medpar_px_img %>%
  filter(str_detect(px_code, paste0("^", BRAIN_MRI_PX_ICD10_PREFIX)) |
           px_code %in% BRAIN_MRI_PX_ICD9) %>%
  first_claim_date("brain_mri_icd_px_date")

pet_icd_px <- medpar_px_img %>%
  filter(str_detect(px_code, PET_PX_ICD10_PATTERN)) %>%
  first_claim_date("pet_icd_px_date")

cohort <- cohort %>%
  left_join(chemo_icd_dx,     by = "PATIENT_ID") %>%
  left_join(chemo_icd_px,     by = "PATIENT_ID") %>%
  left_join(resection_icd_px, by = "PATIENT_ID") %>%
  left_join(resection_anatomic_icd_px, by = "PATIENT_ID") %>%
  left_join(radiation_icd_px, by = "PATIENT_ID") %>%
  left_join(radiation_icd_dx, by = "PATIENT_ID") %>%
  left_join(immuno_icd_px,    by = "PATIENT_ID") %>%
  left_join(brain_mri_icd_px, by = "PATIENT_ID") %>%
  left_join(pet_icd_px,       by = "PATIENT_ID") %>%
  mutate(
    chemo_post12_icd_dx = as.integer(!is.na(chemo_icd_dx_date)),
    chemo_post12_icd_px = as.integer(!is.na(chemo_icd_px_date)),
    chemo_post12 = as.integer(
      chemo_post12_hcpcs_count >= 1L |
        chemo_post12_icd_dx == 1L |
        chemo_post12_icd_px == 1L
    ),
    chemo_post12_count = chemo_post12_hcpcs_count,
    chemo_post12_first_date = pmin(chemo_post12_hcpcs_first_date, chemo_icd_dx_date,
                                   chemo_icd_px_date, na.rm = TRUE),

    resection_post12_icd_px = as.integer(!is.na(resection_icd_px_date)),
    resection_post12 = as.integer(
      resection_post12_hcpcs_count >= 1L |
        resection_post12_icd_px == 1L
    ),
    resection_post12_count = resection_post12_hcpcs_count,
    resection_post12_first_date = pmin(resection_post12_hcpcs_first_date,
                                       resection_icd_px_date, na.rm = TRUE),

    # HCPCS-only value captured before radiation_post12 becomes the three-source union
    radiation_post12_hcpcs  = radiation_post12,
    radiation_post12_icd_px = as.integer(!is.na(radiation_icd_px_date)),
    radiation_post12_icd_dx = as.integer(!is.na(radiation_icd_dx_date)),
    radiation_post12 = as.integer(
      radiation_post12_hcpcs  == 1L |
        radiation_post12_icd_px == 1L |
        radiation_post12_icd_dx == 1L
    ),
    radiation_post12_any = as.integer(
      radiation_post12_any    == 1L |
        radiation_post12_icd_px == 1L |
        radiation_post12_icd_dx == 1L
    ),
    radiation_post12_first_date = pmin(radiation_post12_first_date,
                                       radiation_icd_px_date,
                                       radiation_icd_dx_date, na.rm = TRUE),
    radiation_post12_days_after = as.numeric(radiation_post12_first_date -
                                               dx_date_first),

    immuno_post12_icd_px = as.integer(!is.na(immuno_icd_px_date)),
    immuno_post12_any = as.integer(
      immuno_post12 == 1L | immuno_post12_icd_px == 1L
    ),

    # inpatient imaging: component columns only
    brain_mri_inpatient_icd_px = as.integer(!is.na(brain_mri_icd_px_date)),
    pet_inpatient_icd_px       = as.integer(!is.na(pet_icd_px_date))
  )

rm(medpar_px, medpar_dx, medpar_px_tx, medpar_dx_tx, medpar_px_img)
invisible(gc())

# console diagnostics
cat("Post-diagnosis systemic therapy (12 months after diagnosis):\n")
print(
  tibble(
    therapy = c("Chemotherapy (NCI range-based, primary)",
                "  - administration codes (96400-96549, Q0083-5)",
                "  - drug codes (J9000-J9999, excl. immunotherapy)",
                "  - chemotherapy encounter ICD (Z51.11 / V58.11)",
                "  - INPATIENT administration ICD procedure (MEDPAR)",
                "Chemotherapy (prior enumerated definition)",
                "Immunotherapy (atezo/durva/nivo/pembro J-codes, primary)",
                "  - INPATIENT antineoplastic mAb ICD px (class-level only)",
                "  - either of the above (immuno_post12_any, sensitivity)",
                "Radiation therapy (delivery: carrier + outpatient + INPATIENT)",
                "  - CPT/HCPCS delivery, carrier+outpatient only (was primary)",
                "  - INPATIENT radiation ICD procedure (MEDPAR)",
                "  - radiation encounter ICD (Z51.0 / V58.0)",
                "  - + planning/management/guidance (radiation_post12_any)",
                "Lung resection",
                "  - CPT/HCPCS (outpatient)",
                "  - INPATIENT ICD procedure (MEDPAR)"),
    n       = c(sum(cohort$chemo_post12),
                sum(cohort$chemo_post12_admin),
                sum(cohort$chemo_post12_drug),
                sum(cohort$chemo_post12_icd_dx),
                sum(cohort$chemo_post12_icd_px),
                sum(cohort$chemo_post12_enumerated),
                sum(cohort$immuno_post12),
                sum(cohort$immuno_post12_icd_px),
                sum(cohort$immuno_post12_any),
                sum(cohort$radiation_post12),
                sum(cohort$radiation_post12_hcpcs),
                sum(cohort$radiation_post12_icd_px),
                sum(cohort$radiation_post12_icd_dx),
                sum(cohort$radiation_post12_any),
                sum(cohort$resection_post12),
                sum(cohort$resection_post12_hcpcs_count >= 1L),
                sum(cohort$resection_post12_icd_px))
  ) %>%
    mutate(pct = round(100 * n / nrow(cohort), 1))
)
cat(sprintf(
  "  Radiation: adding inpatient ICD px + Z51.0 moved it %+d patients (%.1f%% -> %.1f%%)\n",
  sum(cohort$radiation_post12) - sum(cohort$radiation_post12_hcpcs),
  100 * mean(cohort$radiation_post12_hcpcs),
  100 * mean(cohort$radiation_post12)))
cat(sprintf("    of which inpatient ICD px alone: %+d; Z51.0/V58.0 alone: %+d\n",
            sum(cohort$radiation_post12_icd_px == 1L &
                  cohort$radiation_post12_hcpcs == 0L),
            sum(cohort$radiation_post12_icd_dx == 1L &
                  cohort$radiation_post12_hcpcs == 0L &
                  cohort$radiation_post12_icd_px == 0L)))
cat(sprintf(
  "  Inpatient imaging (component columns, NOT in the primary flags): brain MRI %d, PET %d\n",
  sum(cohort$brain_mri_inpatient_icd_px), sum(cohort$pet_inpatient_icd_px)))
cat(sprintf("  Range-based vs enumerated: %+d patients (%.1f%% -> %.1f%%)\n",
            sum(cohort$chemo_post12) - sum(cohort$chemo_post12_enumerated),
            100 * mean(cohort$chemo_post12_enumerated),
            100 * mean(cohort$chemo_post12)))


# the six reported measures in the analytic subgroup
tx_ok <- cohort %>% filter(tx_fup_analytic == 1L)

cat(sprintf(
  "\nANALYTIC SUBGROUP: complete 12-month post-dx coverage, n = %d of %d (%.1f%%)\n",
  nrow(tx_ok), nrow(cohort), 100 * nrow(tx_ok) / nrow(cohort)))
cat(sprintf("  Excluded: %d with a coverage gap, %d with no observable month at all\n",
            sum(cohort$tx_fup_complete == 0L & cohort$tx_fup_observable == 1L),
            sum(cohort$tx_fup_observable == 0L)))
cat(sprintf("  Of the %d who died inside the diagnosis month, retained: %d\n",
            sum(!is.na(cohort$death_date) &
                  floor_date(cohort$death_date, "month") == cohort$dx_date_first),
            sum(!is.na(cohort$death_date) &
                  floor_date(cohort$death_date, "month") == cohort$dx_date_first &
                  cohort$tx_fup_analytic == 1L)))

four_measures <- function(d) {
  tibble(
    measure = c("Chemotherapy (12 mo post-dx)", "Immunotherapy (12 mo post-dx)",
                "Radiation therapy (12 mo post-dx)", "Lung resection (12 mo post-dx)",
                "Brain MRI (12 mo pre - 3 mo post)", "PET / PET-CT (12 mo pre - 3 mo post)"),
    n       = c(sum(d$chemo_post12), sum(d$immuno_post12),
                sum(d$radiation_post12), sum(d$resection_post12),
                sum(d$brain_mri_pre12post3), sum(d$pet_pre12post3))
  ) %>%
    mutate(pct = round(100 * n / nrow(d), 1))
}

cat("\n  Full cohort:\n");        print(four_measures(cohort))
cat("\n  Analytic subgroup:\n");  print(four_measures(tx_ok))

cat("\n  By screening status, analytic subgroup:\n")
print(tx_ok %>%
  group_by(scan_type) %>%
  summarize(n = n(),
            pct_chemo     = round(100 * mean(chemo_post12         == 1L), 1),
            pct_immuno    = round(100 * mean(immuno_post12        == 1L), 1),
            pct_radiation = round(100 * mean(radiation_post12     == 1L), 1),
            pct_resection = round(100 * mean(resection_post12     == 1L), 1),
            pct_brain_mri = round(100 * mean(brain_mri_pre12post3 == 1L), 1),
            pct_pet       = round(100 * mean(pet_pre12post3       == 1L), 1),
            .groups = "drop"))


# =============================================================================
# VII. NCCN guideline-concordant treatment
# =============================================================================
cat("Deriving NCCN guideline-concordant treatment...\n")

# signed gap = radiation date minus chemotherapy date; concurrent -30 to +60 days, tight within 14
GCC_CONCURRENT_BEFORE_DAYS <- 30L
GCC_CONCURRENT_AFTER_DAYS  <- 60L

GCC_CONCURRENT_TIGHT_DAYS <- 14L

# any_of is OR over modality combinations and AND within one; tokens map to 0/1 columns
GCC_MODALITY_TOKENS <- c(
  chemo     = "chemo_post12",
  radiation = "radiation_post12",
  resection = "resection_post12",
  immuno    = "immuno_post12"
)

# the seven versions state the same modality-level rules
GCC_GUIDELINE_VERSIONS <- list(
  `NCCN 2.2015` = list(
    citation = paste(
      "NCCN Guidelines Version 2.2015, Small Cell Lung Cancer (released",
      "2015-03/04; in force at mid-2015). No JNCCN companion exists for any",
      "2014-2017 version; nearest full JNCCN statement of the era algorithm:",
      "Kalemkerian GP, Akerley W, Bogner P, et al. J Natl Compr Canc Netw.",
      "2013;11(1):78-98 (PMID 23307984)."
    ),
    stages = list(
      `I-IIA (cT1-2N0)` = list(any_of = list(c("resection", "chemo"),
                                             c("chemo", "radiation"))),
      `IIB-III`         = list(any_of = list(c("chemo", "radiation"))),
      `IV`              = list(any_of = list(c("chemo")))
    )
  ),

  `NCCN 1.2016` = list(
    stages = list(
      `I-IIA (cT1-2N0)` = list(any_of = list(c("resection", "chemo"),
                                             c("chemo", "radiation"))),
      `IIB-III`         = list(any_of = list(c("chemo", "radiation"))),
      `IV`              = list(any_of = list(c("chemo")))
    )
  ),

  `NCCN 3.2017` = list(
    stages = list(
      `I-IIA (cT1-2N0)` = list(any_of = list(c("resection", "chemo"),
                                             c("chemo", "radiation"))),
      `IIB-III`         = list(any_of = list(c("chemo", "radiation"))),
      `IV`              = list(any_of = list(c("chemo")))
    )
  ),

  `NCCN 2.2018` = list(
    stages = list(
      `I-IIA (cT1-2N0)` = list(any_of = list(c("resection", "chemo"),
                                             c("chemo", "radiation"))),
      `IIB-III`         = list(any_of = list(c("chemo", "radiation"))),
      `IV`              = list(any_of = list(c("chemo")))
    )
  ),

  `NCCN 2.2019` = list(
    stages = list(
      `I-IIA (cT1-2N0)` = list(any_of = list(c("resection", "chemo"),
                                             c("chemo", "radiation"))),
      `IIB-III`         = list(any_of = list(c("chemo", "radiation"))),
      `IV`              = list(any_of = list(c("chemo")))
    )
  ),

  `NCCN 3.2020` = list(
    stages = list(
      `I-IIA (cT1-2N0)` = list(any_of = list(c("resection", "chemo"),
                                             c("chemo", "radiation"))),
      `IIB-III`         = list(any_of = list(c("chemo", "radiation"))),
      `IV`              = list(any_of = list(c("chemo")))
    )
  ),

  `NCCN 3.2021` = list(
    stages = list(
      `I-IIA (cT1-2N0)` = list(any_of = list(c("resection", "chemo"),
                                             c("chemo", "radiation"))),
      `IIB-III`         = list(any_of = list(c("chemo", "radiation"))),
      `IV`              = list(any_of = list(c("chemo")))
    )
  )
)

GCC_GUIDELINE_BY_YEAR <- c(
  "2015" = "NCCN 2.2015",
  "2016" = "NCCN 1.2016",
  "2017" = "NCCN 3.2017",
  "2018" = "NCCN 2.2018",
  "2019" = "NCCN 2.2019",
  "2020" = "NCCN 3.2020",
  "2021" = "NCCN 3.2021"
)


# load-time validation of the guideline tables
gcc_validate_guidelines <- function() {
  stage_keys <- c("I-IIA (cT1-2N0)", "IIB-III", "IV")

  yrs <- as.character(START_YEAR:END_YEAR)
  missing_yr <- setdiff(yrs, names(GCC_GUIDELINE_BY_YEAR))
  if (length(missing_yr) > 0) {
    stop(sprintf("GCC_GUIDELINE_BY_YEAR has no entry for diagnosis year(s): %s",
                 paste(missing_yr, collapse = ", ")))
  }
  unknown_ver <- setdiff(unname(GCC_GUIDELINE_BY_YEAR), names(GCC_GUIDELINE_VERSIONS))
  if (length(unknown_ver) > 0) {
    stop(sprintf("GCC_GUIDELINE_BY_YEAR points at undefined version(s): %s",
                 paste(unknown_ver, collapse = ", ")))
  }

  for (v in names(GCC_GUIDELINE_VERSIONS)) {
    spec <- GCC_GUIDELINE_VERSIONS[[v]]
    if (!setequal(names(spec$stages), stage_keys)) {
      stop(sprintf("Guideline '%s' must define exactly these stage keys: %s",
                   v, paste(stage_keys, collapse = " | ")))
    }
    for (s in stage_keys) {
      combos <- spec$stages[[s]]$any_of
      if (!is.list(combos) || length(combos) == 0L) {
        stop(sprintf("Guideline '%s' stage '%s' has an empty any_of", v, s))
      }
      bad <- setdiff(unlist(combos), names(GCC_MODALITY_TOKENS))
      if (length(bad) > 0) {
        stop(sprintf(paste0("Guideline '%s' stage '%s' uses unknown modality ",
                            "token(s): %s. Valid tokens: %s"),
                     v, s, paste(bad, collapse = ", "),
                     paste(names(GCC_MODALITY_TOKENS), collapse = ", ")))
      }
    }
  }
  invisible(TRUE)
}
gcc_validate_guidelines()

gcc_eval_rule <- function(df, rule) {
  Reduce(`|`, lapply(rule$any_of, function(combo) {
    Reduce(`&`, lapply(combo, function(tok) {
      col <- GCC_MODALITY_TOKENS[[tok]]
      if (!col %in% names(df)) {
        stop(sprintf("Modality token '%s' maps to column '%s', which does not exist",
                     tok, col))
      }
      df[[col]] == 1L
    }))
  }))
}

gcc_score <- function(df) {
  out <- rep(NA_integer_, nrow(df))
  ver  <- unname(GCC_GUIDELINE_BY_YEAR[as.character(df$yr)])
  for (v in unique(stats::na.omit(ver))) {
    rules <- GCC_GUIDELINE_VERSIONS[[v]]$stages
    for (s in names(rules)) {
      sel <- which(ver == v & !is.na(df$gcc_rule) & df$gcc_rule == s)
      if (length(sel) == 0L) next
      out[sel] <- as.integer(gcc_eval_rule(df[sel, , drop = FALSE], rules[[s]]))
    }
  }
  out
}

gcc_rule_text <- function(rule) {
  paste(vapply(rule$any_of, function(c) paste(c, collapse = " + "), character(1)),
        collapse = " OR ")
}

regimen_first_date <- function(codes, out_name) {
  systemic_post %>%
    filter(HCPCS_CD %in% codes) %>%
    group_by(PATIENT_ID) %>%
    summarize(!!out_name := min(proc_date), .groups = "drop")
}

cohort <- cohort %>%
  # drug-level flags are HCPCS-only: MEDPAR records antineoplastic, never the agent
  left_join(regimen_first_date(PLATINUM_HCPCS,   "platinum_first_date"),   by = "PATIENT_ID") %>%
  left_join(regimen_first_date(ETOPOSIDE_HCPCS,  "etoposide_first_date"),  by = "PATIENT_ID") %>%
  left_join(regimen_first_date(IRINOTECAN_HCPCS, "irinotecan_first_date"), by = "PATIENT_ID") %>%
  left_join(regimen_first_date(IO_FIRSTLINE_ES_HCPCS,
                               "io_firstline_first_date"),               by = "PATIENT_ID") %>%
  left_join(regimen_first_date(RESECTION_ANATOMIC_CPT,
                               "resection_anatomic_hcpcs_first_date"),   by = "PATIENT_ID") %>%
  mutate(
    # regimen and extent detail, descriptive only
    platinum_post12   = as.integer(!is.na(platinum_first_date)),
    etoposide_post12  = as.integer(!is.na(etoposide_first_date)),
    irinotecan_post12 = as.integer(!is.na(irinotecan_first_date)),
    platinum_etoposide_post12 = as.integer(platinum_post12 == 1L &
                                             etoposide_post12 == 1L),
    platinum_doublet_es_post12 = as.integer(
      platinum_post12 == 1L & (etoposide_post12 == 1L | irinotecan_post12 == 1L)
    ),
    io_firstline_post12 = as.integer(!is.na(io_firstline_first_date)),

    resection_anatomic_post12 = as.integer(
      !is.na(resection_anatomic_hcpcs_first_date) |
        !is.na(resection_anatomic_icd_px_date)
    ),
    resection_anatomic_first_date = pmin(resection_anatomic_hcpcs_first_date,
                                         resection_anatomic_icd_px_date,
                                         na.rm = TRUE),

    chemort_gap_days = as.numeric(radiation_post12_first_date -
                                    chemo_post12_first_date),
    chemort_any_post12 = as.integer(chemo_post12 == 1L & radiation_post12 == 1L),
    no_treatment_post12 = as.integer(
      chemo_post12 == 0L & radiation_post12 == 0L &
        resection_post12 == 0L & immuno_post12 == 0L
    ),
    chemort_concurrent_post12 = as.integer(
      chemort_any_post12 == 1L &
        !is.na(chemort_gap_days) &
        chemort_gap_days >= -GCC_CONCURRENT_BEFORE_DAYS &
        chemort_gap_days <=  GCC_CONCURRENT_AFTER_DAYS
    ),
    chemort_concurrent_tight_post12 = as.integer(
      chemort_any_post12 == 1L &
        !is.na(chemort_gap_days) &
        abs(chemort_gap_days) <= GCC_CONCURRENT_TIGHT_DAYS
    ),
    chemort_sequential_post12 = as.integer(
      chemort_any_post12 == 1L & chemort_concurrent_post12 == 0L
    ),

    # six-month modality flags for gcc_180d
    chemo_post6 = as.integer(!is.na(chemo_post12_first_date) &
      as.numeric(chemo_post12_first_date - dx_date_first) <= 180),
    radiation_post6 = as.integer(!is.na(radiation_post12_first_date) &
      as.numeric(radiation_post12_first_date - dx_date_first) <= 180),
    resection_post6 = as.integer(!is.na(resection_post12_first_date) &
      as.numeric(resection_post12_first_date - dx_date_first) <= 180),
    immuno_post6 = as.integer(!is.na(immuno_post12_first_date) &
      as.numeric(immuno_post12_first_date - dx_date_first) <= 180),

    # stage tiers: the AJCC reconstruction first, summary stage as the fallback
    nccn_stage = case_when(
      stage_simple %in% c("Stage IA", "Stage IB")                   ~ "I",
      stage_simple %in% c("Stage IIA", "Stage IIB")                 ~ "II",
      stage_simple %in% c("Stage IIIA", "Stage IIIB", "Stage IIIC") ~ "III",
      stage_simple == "Stage IV"                                    ~ "IV",
      TRUE                                                          ~ NA_character_
    ),

    # cT1-2N0 is AJCC 8th stage I plus the T2bN0 cell of stage IIA
    nccn_surgical_eligible_stage = as.integer(
      stage_simple %in% c("Stage IA", "Stage IB", "Stage IIA")
    ),

    sclc_stage_group = case_when(
      nccn_stage %in% c("I", "II", "III") ~ "Limited",
      nccn_stage == "IV"                  ~ "Extensive",
      staging %in% c(1L, 2L)              ~ "Limited",
      staging == 3L                       ~ "Extensive",
      TRUE                                ~ NA_character_
    ),

    gcc_stage_basis = case_when(
      !is.na(nccn_stage)         ~ "AJCC TNM (SEER registry-reported)",
      staging %in% c(1L, 2L, 3L) ~ "SEER Combined Summary Stage",
      TRUE                       ~ "Unstageable"
    ),

    # evaluated top-down: Stage IIA meets the surgical-eligible rule before the nccn_stage II test
    gcc_rule = case_when(
      nccn_surgical_eligible_stage == 1L               ~ "I-IIA (cT1-2N0)",
      !is.na(nccn_stage) & nccn_stage %in% c("II","III") ~ "IIB-III",
      !is.na(nccn_stage) & nccn_stage == "IV"          ~ "IV",
      is.na(nccn_stage) & staging == 1L                ~ "I-IIA (cT1-2N0)",
      is.na(nccn_stage) & staging == 2L                ~ "IIB-III",
      is.na(nccn_stage) & staging == 3L                ~ "IV",
      TRUE                                             ~ NA_character_
    )
  )

cohort$gcc_guideline <- unname(GCC_GUIDELINE_BY_YEAR[as.character(cohort$yr)])
cohort$gcc           <- gcc_score(cohort)

# the same engine with each modality column swapped for its six-month counterpart
cohort$gcc_180d <- local({
  has_6mo <- c("chemo", "radiation", "resection", "immuno")
  used <- unique(unlist(lapply(GCC_GUIDELINE_VERSIONS, function(v)
    lapply(v$stages, function(s) unlist(s$any_of)))))
  no_6mo <- setdiff(used, has_6mo)
  if (length(no_6mo) > 0) {
    warning(sprintf(paste0("gcc_180d: rules use token(s) with no 180-day ",
                           "counterpart (%s); those keep their 12-month value, so ",
                           "gcc_180d mixes windows. Derive <token>_post6 above to fix."),
                    paste(no_6mo, collapse = ", ")))
  }
  df6 <- cohort
  df6$chemo_post12     <- df6$chemo_post6
  df6$radiation_post12 <- df6$radiation_post6
  df6$resection_post12 <- df6$resection_post6
  df6$immuno_post12    <- df6$immuno_post6
  gcc_score(df6)
})

cohort <- cohort %>%
  mutate(
    # gcc_eligible and the early-death flags are applied as denominators in 02
    gcc_category = case_when(
      is.na(gcc)                 ~ NA_character_,
      gcc == 1L                  ~ "Guideline-concordant",
      no_treatment_post12 == 1L  ~ "No treatment",
      TRUE                       ~ "Less than recommended"
    ),

    gcc_eligible = as.integer(!is.na(gcc) & tx_fup_analytic == 1L),

    early_death_30d = as.integer(!is.na(death_date) &
                                   as.numeric(death_date - dx_date_first) <= 30),
    early_death_60d = as.integer(!is.na(death_date) &
                                   as.numeric(death_date - dx_date_first) <= 60)
  )

# console diagnostics
cat("\n  GUIDELINE IN FORCE (from GCC_GUIDELINE_VERSIONS / GCC_GUIDELINE_BY_YEAR):\n")
for (v in unique(unname(GCC_GUIDELINE_BY_YEAR))) {
  yrs <- names(GCC_GUIDELINE_BY_YEAR)[GCC_GUIDELINE_BY_YEAR == v]
  spec <- GCC_GUIDELINE_VERSIONS[[v]]
  cat(sprintf("    %s  ->  diagnosis year%s %s (n = %d)\n",
              v, if (length(yrs) > 1) "s" else "",
              paste(unique(range(yrs)), collapse = "-"),
              sum(cohort$gcc_guideline == v)))
  for (s in names(spec$stages)) {
    cat(sprintf("      %-18s concordant if: %s\n", s,
                gcc_rule_text(spec$stages[[s]])))
  }
}

cat(sprintf("\n  Stage basis: AJCC TNM %d (%.1f%%) | SEER Summary Stage %d (%.1f%%) | unstageable %d (%.1f%%)\n",
            sum(cohort$gcc_stage_basis == "AJCC TNM (SEER registry-reported)"),
            100 * mean(cohort$gcc_stage_basis == "AJCC TNM (SEER registry-reported)"),
            sum(cohort$gcc_stage_basis == "SEER Combined Summary Stage"),
            100 * mean(cohort$gcc_stage_basis == "SEER Combined Summary Stage"),
            sum(cohort$gcc_stage_basis == "Unstageable"),
            100 * mean(cohort$gcc_stage_basis == "Unstageable")))

cat("\n  Guideline-concordant treatment by NCCN rule:\n")
print(cohort %>%
  filter(!is.na(gcc_rule)) %>%
  group_by(gcc_rule, gcc_stage_basis) %>%
  summarize(n            = n(),
            pct_gcc      = round(100 * mean(gcc == 1L), 1),
            pct_untreated = round(100 * mean(no_treatment_post12 == 1L), 1),
            .groups = "drop") %>%
  arrange(gcc_rule, gcc_stage_basis), n = Inf)

cat("\n  Treatment detail (descriptive; none of these enters gcc):\n")
print(tibble(
  component = c("Platinum (cisplatin or carboplatin)",
                "Etoposide",
                "Irinotecan (extensive-stage platinum partner only)",
                "Platinum + etoposide doublet",
                "Platinum + (etoposide or irinotecan), extensive-stage doublet",
                "First-line checkpoint inhibitor (atezolizumab/durvalumab)",
                "Any checkpoint inhibitor (immuno_post12, for contrast)",
                "Anatomic resection (lobectomy or greater, or segmentectomy)",
                "Any resection (resection_post12, for contrast)",
                "Chemotherapy + radiation, any interval",
                "  - concurrent, Corso window (RT -30 to +60 d from chemo)",
                "  - concurrent, tight (|gap| <= 14 d, Bobbili)",
                "  - sequential (outside the Corso window)"),
  n = c(sum(cohort$platinum_post12), sum(cohort$etoposide_post12),
        sum(cohort$irinotecan_post12),
        sum(cohort$platinum_etoposide_post12),
        sum(cohort$platinum_doublet_es_post12),
        sum(cohort$io_firstline_post12),
        sum(cohort$immuno_post12), sum(cohort$resection_anatomic_post12),
        sum(cohort$resection_post12), sum(cohort$chemort_any_post12),
        sum(cohort$chemort_concurrent_post12),
        sum(cohort$chemort_concurrent_tight_post12),
        sum(cohort$chemort_sequential_post12))
) %>% mutate(pct = round(100 * n / nrow(cohort), 1)))

cat("\n  Three-level version (Blom 2020 reporting convention):\n")
print(cohort %>% filter(!is.na(gcc_category)) %>% count(gcc_category) %>%
        mutate(pct = round(100 * n / sum(n), 1)))

cat(sprintf("\n  Overall: gcc = %d of %d scorable (%.1f%%) | NA (unstageable) %d\n",
            sum(cohort$gcc == 1L, na.rm = TRUE), sum(!is.na(cohort$gcc)),
            100 * mean(cohort$gcc == 1L, na.rm = TRUE),
            sum(is.na(cohort$gcc))))
cat(sprintf("  6-month window (gcc_180d, Caprario 2013): %.1f%% vs %.1f%% at 12 months (%+.1f pp)\n",
            100 * mean(cohort$gcc_180d == 1L, na.rm = TRUE),
            100 * mean(cohort$gcc == 1L, na.rm = TRUE),
            100 * (mean(cohort$gcc_180d == 1L, na.rm = TRUE) -
                     mean(cohort$gcc == 1L, na.rm = TRUE))))
cat(sprintf("  Restricted to complete 12-month coverage (gcc_eligible): %.1f%% concordant of %d\n",
            100 * mean(cohort$gcc[cohort$gcc_eligible == 1L] == 1L),
            sum(cohort$gcc_eligible == 1L)))
cat(sprintf("  Excluding deaths within 60 days of the window anchor: %.1f%% concordant of %d\n",
            100 * mean(cohort$gcc[cohort$gcc_eligible == 1L &
                                    cohort$early_death_60d == 0L] == 1L),
            sum(cohort$gcc_eligible == 1L & cohort$early_death_60d == 0L)))


# =============================================================================
# VIII. Comorbidity flags
# =============================================================================
cat("Deriving comorbidities...\n")

# pre-diagnosis code prefixes; the Charlson-Deyo scores used for adjustment are computed in 02
dx_prefixes <- list(
  asthma = list(
    icd9  = c("493"),
    icd10 = c("J45")
  ),
  ischemic_heart_disease = list(
    icd9  = c("410","411","412","413","414","4292"),
    icd10 = c("I20","I21","I22","I23","I24","I25")
  ),
  cvd = list(
    icd9  = c("394","395","396","397","424","425"),
    icd10 = c("I05","I06","I07","I08","I09","I34","I35","I36","I37","I38","I42")
  ),
  mi = list(
    icd9  = c("410","412"),
    icd10 = c("I21","I22")
  ),
  copd = list(
    icd9  = c("490","491","492","494","495","496"),
    icd10 = c("J40","J41","J42","J43","J44","J47")
  ),
  chf = list(
    icd9  = c("428","4254","4255","4257","4258","4259"),
    icd10 = c("I50","I110","I130","I132")
  ),
  pvd = list(
    icd9  = c("440","441","4431","4432","4433","4434","4435","4436","4437",
              "4438","4439","4471","5571","5579"),
    icd10 = c("I70","I71","I731","I738","I739","I771","I790","I792",
              "K551","K558","K559","Z958","Z959")
  ),
  cerebrovascular_disease = list(
    icd9  = c("430","431","432","433","434","435","436","437","438"),
    icd10 = c("G45","G46","H340","I60","I61","I62","I63","I64","I65",
              "I66","I67","I68","I69")
  ),
  diab_uc = list(
    icd9  = c("2500","2501","2502","2503","2508","2509"),
    icd10 = c("E100","E101","E109","E110","E111","E119","E120","E121",
              "E129","E130","E131","E139","E140","E141","E149")
  ),
  interstitial_lung = list(
    icd9  = c("515","5160","5161","5162","5163","5164","5168","5169",
              "5171","5178","5179"),
    icd10 = c("J70","J84")
  ),
  o2 = list(
    icd9  = c("V462"),
    icd10 = c("Z9981")
  ),
  kidney_diseases = list(
    icd9  = c("582","583","585","586","5880","V420","V451","V56"),
    icd10 = c("N00","N01","N02","N03","N04","N05","N06","N07","N08",
              "N10","N11","N12","N13","N14","N15","N16","N17","N18","N19","Q61")
  )
)

dx_path <- file.path(EXTRACTED, "lung_diagnosis.parquet")

if (file.exists(dx_path)) {
  dx_wide <- read_parquet(dx_path) %>%
    filter(PATIENT_ID %in% pids) %>%
    left_join(cohort %>% select(PATIENT_ID, dx_date_end), by = "PATIENT_ID") %>%
    filter(claim_date_of(CLM_THRU_DT) < dx_date_end)

  dx_code_cols <- names(dx_wide)[
    str_detect(names(dx_wide), "ICD_DGNS|PRNCPAL_DGNS|FST_DGNS")
  ]
  dx_long <- dx_wide %>%
    select(PATIENT_ID, CLM_THRU_DT, all_of(dx_code_cols)) %>%
    pivot_longer(all_of(dx_code_cols), names_to = NULL, values_to = "dx_code") %>%
    filter(!is.na(dx_code), dx_code != "") %>%
    mutate(
      claim_date = claim_date_of(CLM_THRU_DT),
      icd_ver    = if_else(claim_date >= ICD10_START, "icd10", "icd9"),
      dx_code    = str_to_upper(str_remove_all(str_trim(dx_code), "\\."))
    ) %>%
    distinct(PATIENT_ID, claim_date, dx_code, icd_ver)

  make_pattern <- function(prefixes) {
    paste0("^(", paste(str_replace_all(prefixes, "\\.", "\\\\."), collapse = "|"), ")")
  }

  comorbidity_counts <- cohort %>% select(PATIENT_ID)

  for (cx in names(dx_prefixes)) {
    pat9  <- make_pattern(dx_prefixes[[cx]]$icd9)
    pat10 <- make_pattern(dx_prefixes[[cx]]$icd10)

    counts <- dx_long %>%
      mutate(hit = (icd_ver == "icd9"  & str_detect(dx_code, pat9)) |
                   (icd_ver == "icd10" & str_detect(dx_code, pat10))) %>%
      group_by(PATIENT_ID) %>%
      summarize(!!sprintf("%s_pre_count", cx) := sum(hit, na.rm = TRUE),
                .groups = "drop")

    comorbidity_counts <- comorbidity_counts %>% left_join(counts, by = "PATIENT_ID")
  }

  comorbidity_counts <- comorbidity_counts %>%
    mutate(across(ends_with("_pre_count"), ~ replace_na(.x, 0L)))

  comorbidity_flags <- comorbidity_counts %>%
    mutate(across(
      ends_with("_pre_count"),
      ~ as.integer(.x >= 1L),
      .names = "{str_remove(.col, '_pre_count')}"
    ))

  cohort <- cohort %>% left_join(comorbidity_flags, by = "PATIENT_ID")

  cat("Comorbidity prevalence (% of cohort with pre-dx claim):\n")
  pct <- sapply(names(dx_prefixes), function(cx) {
    round(100 * mean(cohort[[cx]] == 1L, na.rm = TRUE), 1)
  })
  print(pct)
} else {
  cat("  Warning: lung_diagnosis.parquet not found — comorbidities set to NA\n")
}


# =============================================================================
# IX. Outputs
# =============================================================================

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
write_parquet(cohort, file.path(OUTPUT_DIR, "sclc_cohort.parquet"))


cat("Writing code definition table...\n")

# code definition table; the notes column comes from code_definition_notes.R
lcs_code_claims <- cohort %>%
  select(PATIENT_ID, dx_date_end) %>%
  inner_join(proc_claims %>% filter(HCPCS_CD %in% LCS_CPTS), by = "PATIENT_ID") %>%
  filter(as.numeric(dx_date_end - proc_date) >= 0)

LCS_WINDOW     <- "Any claim on or before the end of the diagnosis month"
IMG_WINDOW     <- "12 months before through 3 months after diagnosis"
THERAPY_WINDOW <- "First day of the diagnosis month through 365 days after"

code_definitions <- tribble(
  ~variable, ~codes, ~window,

  "Lung cancer screening (LDCT)",
    paste(LCS_CPTS, collapse = ", "),
    LCS_WINDOW,

  "Brain MRI (CPT, carrier + outpatient)",
    paste(BRAIN_MRI_CPTS, collapse = ", "),
    IMG_WINDOW,

  "Brain MRI (ICD, MEDPAR INPATIENT)",
    paste0("ICD-10-PCS ", BRAIN_MRI_PX_ICD10_PREFIX, "* (MRI of Brain); ",
           "ICD-9-CM ", paste(BRAIN_MRI_PX_ICD9, collapse = ", ")),
    IMG_WINDOW,

  "PET / PET-CT (CPT, carrier + outpatient)",
    paste(PET_CPTS, collapse = ", "),
    IMG_WINDOW,

  "PET / PET-CT (ICD, MEDPAR INPATIENT)",
    paste0("ICD-10-PCS ", PET_PX_ICD10_PATTERN,
           " (Nuclear Medicine, root type 3 = PET, any body system); no ICD-9-CM equivalent exists"),
    IMG_WINDOW,

  "Chemotherapy (CPT/HCPCS, carrier + outpatient)",
    paste(
      "CPT 96400-96549, Q0083-Q0085",
      "HCPCS J9000-J9999 (excluding the immunotherapy codes below)",
      sep = "; "),
    THERAPY_WINDOW,

  "Chemotherapy (ICD-CM encounter diagnosis)",
    paste(c(CHEMO_DX_ICD10, CHEMO_DX_ICD9), collapse = ", "),
    THERAPY_WINDOW,

  "Chemotherapy (ICD, MEDPAR INPATIENT)",
    paste(
      paste0("ICD-9-CM ", paste(CHEMO_PX_ICD9, collapse = ", ")),
      paste0("ICD-10-PCS ", CHEMO_PX_ICD10_PATTERN,
             " (substance 0 = Antineoplastic, qualifier 0/5, any body part/approach)"),
      paste0("ICD-10-PCS New Technology, enumerated: ",
             paste(CHEMO_PX_ICD10_NEWTECH, collapse = ", ")),
      sep = "; "),
    THERAPY_WINDOW,

  "Immunotherapy (HCPCS J-codes, primary)",
    paste0("Unrestricted, no date windowing: ", paste(IMMUNO_HCPCS_SCLC, collapse = ", ")),
    THERAPY_WINDOW,

  "Immunotherapy (ICD, MEDPAR INPATIENT)",
    paste0("ICD-10-PCS ", IMMUNO_PX_ICD10_PATTERN,
           " (substance 0 = Antineoplastic; qualifier M = Monoclonal Antibody, 2/3 = high/low-dose Interleukin-2)"),
    THERAPY_WINDOW,

  "Radiation therapy (CPT/HCPCS, carrier + outpatient)",
    paste(
      paste0("Delivery: ", paste(RADIATION_DELIVERY_CPT, collapse = ", ")),
      paste0("Supportive (planning/simulation/management/guidance, secondary only): ",
             paste(RADIATION_SUPPORTIVE_CPT, collapse = ", ")),
      sep = "; "),
    THERAPY_WINDOW,

  "Radiation therapy (ICD, MEDPAR INPATIENT)",
    paste(
      paste0("ICD-9-CM ", paste(RADIATION_PX_ICD9_PREFIX, collapse = "/"), "x"),
      paste0("ICD-10-PCS section ", RADIATION_PX_ICD10_SECTION,
             " excluding modality-Y qualifiers ",
             paste(RADIATION_PX_ICD10_EXCL_MQ, collapse = "/")),
      sep = "; "),
    THERAPY_WINDOW,

  "Radiation therapy (ICD-CM encounter diagnosis)",
    paste(c(RADIATION_DX_ICD10, RADIATION_DX_ICD9), collapse = ", "),
    THERAPY_WINDOW,

  "Lung resection (CPT, carrier + outpatient)",
    paste0("CPT ", paste(RESECTION_CPT, collapse = ", ")),
    THERAPY_WINDOW,

  "Lung resection (ICD, MEDPAR INPATIENT)",
    paste(
      paste0("ICD-9-CM ", paste(RESECTION_PX_ICD9_PREFIX, collapse = "/"), "*"),
      paste0("ICD-10-PCS ", paste(RESECTION_PX_ICD10, collapse = ", ")),
      sep = "; "),
    THERAPY_WINDOW,

  "NCCN concordance -- platinum (HCPCS)",
    paste(PLATINUM_HCPCS, collapse = ", "),
    THERAPY_WINDOW,

  "NCCN concordance -- etoposide (HCPCS)",
    paste(ETOPOSIDE_HCPCS, collapse = ", "),
    THERAPY_WINDOW,

  "NCCN concordance -- irinotecan (HCPCS)",
    paste(IRINOTECAN_HCPCS, collapse = ", "),
    THERAPY_WINDOW,

  "NCCN concordance -- first-line checkpoint inhibitor (HCPCS)",
    paste(IO_FIRSTLINE_ES_HCPCS, collapse = ", "),
    THERAPY_WINDOW,

  "NCCN concordance -- anatomic resection (CPT + ICD)",
    paste(
      paste0("CPT ", paste(RESECTION_ANATOMIC_CPT, collapse = ", ")),
      paste0("ICD-9-CM ", paste(RESECTION_ANATOMIC_PX_ICD9_PREFIX, collapse = "/"), "*"),
      paste0("ICD-10-PCS ", paste(RESECTION_ANATOMIC_PX_ICD10, collapse = ", "),
             " (root operation T = Resection, i.e. removal of an entire lobe or lung)"),
      sep = "; "),
    THERAPY_WINDOW,

  # generated from whichever guideline tables are in force at run time
  "NCCN guideline-concordant treatment (derived, not a code list)",
    paste(vapply(names(GCC_GUIDELINE_VERSIONS[[GCC_GUIDELINE_BY_YEAR[[1]]]]$stages),
                 function(s) sprintf("%s: %s", s,
                   gcc_rule_text(GCC_GUIDELINE_VERSIONS[[GCC_GUIDELINE_BY_YEAR[[1]]]]$stages[[s]])),
                 character(1)), collapse = "; "),
    THERAPY_WINDOW,

  "Chemoradiation timing (derived, descriptive)",
    "chemo_post12 = 1 AND radiation_post12 = 1, split by the signed chemo-to-RT gap",
    THERAPY_WINDOW,

  "NCCN concordance, 6-month sensitivity (derived)",
    "Same rules as gcc, every modality restricted to the first 180 days",
    "0-180 days from the first of the diagnosis month",

  "NCCN concordance, three-level (derived)",
    "Guideline-concordant / Less than recommended / No treatment",
    THERAPY_WINDOW
) %>%
  mutate(notes = unname(CODE_DEFINITION_NOTES[variable]))

if (anyNA(code_definitions$notes)) {
  stop(sprintf("code_definition_notes.R has no entry for: %s",
               paste(code_definitions$variable[is.na(code_definitions$notes)],
                     collapse = "; ")))
}
write_csv(code_definitions, file.path(OUTPUT_DIR, "code_definitions.csv"))


# guideline specification, one row per version and stage group
collapse_years <- function(y) {
  if (length(y) == 0L) return("(defined, not applied)")
  y <- sort(unique(as.integer(y)))
  brk <- cumsum(c(1L, diff(y) != 1L))
  paste(vapply(split(y, brk), function(g) {
    if (length(g) == 1L) as.character(g) else paste0(min(g), "-", max(g))
  }, character(1)), collapse = "; ")
}

guideline_spec <- bind_rows(lapply(names(GCC_GUIDELINE_VERSIONS), function(v) {
  spec <- GCC_GUIDELINE_VERSIONS[[v]]
  yrs  <- names(GCC_GUIDELINE_BY_YEAR)[GCC_GUIDELINE_BY_YEAR == v]
  bind_rows(lapply(names(spec$stages), function(s) {
    tibble(
      guideline_version   = v,
      applied_to_dx_years = collapse_years(yrs),
      stage_group         = s,
      concordant_if       = gcc_rule_text(spec$stages[[s]]),
      citation            = spec$citation
    )
  }))
}))
write_csv(guideline_spec, file.path(OUTPUT_DIR, "guideline_definitions.csv"))
cat(sprintf("  Guideline spec -> %s (%d stage rules)\n",
            file.path(OUTPUT_DIR, "guideline_definitions.csv"),
            nrow(guideline_spec)))


# patients per HCPCS code within each variable's window
observed <- bind_rows(
  lcs_code_claims %>%
    transmute(variable = "Lung cancer screening (LDCT)", HCPCS_CD, PATIENT_ID),
  imaging_win %>%
    transmute(variable = if_else(modality == "brain_mri", "Brain MRI", "PET / PET-CT"),
              HCPCS_CD, PATIENT_ID),
  systemic_post %>%
    transmute(variable = case_when(
                HCPCS_CD %in% IMMUNO_HCPCS_ALL   ~ "Immunotherapy",
                HCPCS_CD %in% RADIATION_CPT      ~ "Radiation therapy",
                HCPCS_CD %in% RESECTION_CPT      ~ "Lung resection",
                HCPCS_CD %in% CHEMO_ADMIN_RANGE ~ "Chemotherapy - administration",
                TRUE                            ~ "Chemotherapy - drug"),
              HCPCS_CD, PATIENT_ID)
) %>%
  distinct(variable, HCPCS_CD, PATIENT_ID) %>%
  count(variable, code = HCPCS_CD, name = "n_patients") %>%
  filter(n_patients > 0) %>%
  mutate(
    description = unname(CODE_LABELS[code]),
    description = if_else(is.na(description),
                          "(matched by range; not individually enumerated)",
                          description),
    pct_cohort  = round(100 * n_patients / nrow(cohort), 2)
  ) %>%
  arrange(variable, desc(n_patients))

write_csv(observed, file.path(OUTPUT_DIR, "observed_codes.csv"))

# every enumerated code must have a CODE_LABELS entry
enumerated <- unique(c(LCS_CPTS, BRAIN_MRI_CPTS, PET_CPTS,
                       CHEMO_ADMIN_CPTS, CHEMO_DRUG_HCPCS, IMMUNO_HCPCS_ALL,
                       RADIATION_CPT, RESECTION_CPT,
                       PLATINUM_HCPCS, ETOPOSIDE_HCPCS, IRINOTECAN_HCPCS,
                       IO_FIRSTLINE_ES_HCPCS, RESECTION_ANATOMIC_CPT))
undescribed <- setdiff(enumerated, names(CODE_LABELS))
if (length(undescribed) > 0) {
  stop(sprintf("CODE_LABELS has no description for: %s",
               paste(undescribed, collapse = ", ")))
}

cat(sprintf("  Specification -> %s (%d rows)\n",
            file.path(OUTPUT_DIR, "code_definitions.csv"), nrow(code_definitions)))
cat(sprintf("  Observed codes -> %s (%d codes with >=1 patient)\n",
            file.path(OUTPUT_DIR, "observed_codes.csv"), nrow(observed)))
print(observed %>% group_by(variable) %>%
        summarize(codes_with_claims = n(),
                  reached_only_by_range = sum(str_starts(description, "\\(matched")),
                  .groups = "drop"))


cat(sprintf(
  "\nExtraction complete. Final N: %d  |  Screened: %d  |  Missed Scan: %d  |  Never Screened: %d\n",
  nrow(cohort),
  sum(cohort$scan_type == "Screened",       na.rm = TRUE),
  sum(cohort$scan_type == "Missed Scan",    na.rm = TRUE),
  sum(cohort$scan_type == "Never Screened", na.rm = TRUE)
))

# closing console summary
cat("\n--- Chemotherapy (SEER recode) ---\n")
print(cohort %>% count(chemotherapy))
cat("\n--- Radiation (SEER recode) ---\n")
print(cohort %>% count(radiation) %>% arrange(desc(n)))
cat("\n--- Radiation collapsed ---\n")
print(cohort %>% count(radiation_yes_no_unk))
cat("\n--- Radiation therapy, claims-based (radiation_post12) ---\n")
print(cohort %>% count(radiation_post12))
cat("\n--- Lung resection, claims-based (resection_post12) ---\n")
print(cohort %>% count(resection_post12))
cat("\n--- Peri-diagnosis imaging by screening status ---\n")
print(cohort %>%
  group_by(scan_type) %>%
  summarize(
    n             = n(),
    pct_brain_mri = round(100 * mean(brain_mri_pre12post3 == 1L), 1),
    pct_pet       = round(100 * mean(pet_pre12post3       == 1L), 1),
    .groups = "drop"
  ))
cat("  (imaging window: 12 months before to 3 months after diagnosis)\n")
cat("\n--- SEER T (all eras, clinical or pathologic as SEER reports it) ---\n")
print(cohort %>% count(seer_t_group))
cat("\n--- SEER N ---\n")
print(cohort %>% count(seer_n_group))
cat("\n--- SEER stage group (raw, by source field) ---\n")
print(cohort %>% count(seer_stage_source, seer_stage_group), n = Inf)
cat("\n--- SEER T/N availability by dx year (TNM-reporting registries only) ---\n")
print(cohort %>%
  filter(tnm_reported_registry == 1L) %>%
  group_by(yr, seer_stage_source) %>%
  summarize(n = n(),
            pct_seer_t = round(100 * mean(seer_t_group != "Unknown/not available"), 1),
            pct_seer_n = round(100 * mean(seer_n_group != "Unknown/not available"), 1),
            .groups = "drop"))

# agreement between the registry-reported T/N groups and the derived tnm_t/tnm_n
tn_check <- cohort %>%
  transmute(
    t_ok = (seer_t_group == "Unknown/not available" & is.na(tnm_t)) |
           (seer_t_group == paste0("T", tnm_t)),
    n_ok = (seer_n_group == "Unknown/not available" & is.na(tnm_n)) |
           (seer_n_group == paste0("N", tnm_n))
  )
cat(sprintf("\nseer_t_group vs tnm_t agreement: %.2f%% | seer_n_group vs tnm_n: %.2f%%\n",
            100 * mean(tn_check$t_ok, na.rm = TRUE),
            100 * mean(tn_check$n_ok, na.rm = TRUE)))

cat(sprintf(
  "\nRegistries not reporting TNM: %d of %d cohort patients (%.1f%%) — seer_t/seer_n and stage_simple are structurally missing for these.\n",
  sum(cohort$tnm_reported_registry == 0L), nrow(cohort),
  100 * mean(cohort$tnm_reported_registry == 0L)))
cat(sprintf("Registries not reporting radiation: %d of %d (%.1f%%)\n",
  sum(cohort$radiation_reported == 0L), nrow(cohort),
  100 * mean(cohort$radiation_reported == 0L)))

yearly <- cohort %>%
  group_by(yr) %>%
  summarize(
    n_total      = n(),
    n_screened   = sum(screened == 1L, na.rm = TRUE),
    pct_screened = round(100 * n_screened / n_total, 1),
    .groups = "drop"
  ) %>%
  arrange(yr)

cat("\n--- Cohort and Screening by Diagnosis Year ---\n")
cat(sprintf("%-6s  %8s  %10s  %12s\n", "Year", "N Total", "N Screened", "% Screened"))
cat(strrep("-", 42), "\n")
for (i in seq_len(nrow(yearly))) {
  cat(sprintf("%-6d  %8d  %10d  %11.1f%%\n",
              yearly$yr[i], yearly$n_total[i],
              yearly$n_screened[i], yearly$pct_screened[i]))
}
cat(strrep("-", 42), "\n")
cat(sprintf("%-6s  %8d  %10d  %11.1f%%\n",
            "Total", nrow(cohort),
            sum(cohort$screened, na.rm = TRUE),
            round(100 * mean(cohort$screened, na.rm = TRUE), 1)))
