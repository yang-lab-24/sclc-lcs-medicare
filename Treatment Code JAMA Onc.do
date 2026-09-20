*Dx date 
gen dxdate = date(dx_date_end, "YMD")
gen surgdate_hcpcs = date(resection_post12_hcpcs_first_dat, "YMD")
gen surgdate_icd = date(resection_icd_px_date, "YMD")

*Chemotherapy
destring chemo_post12_hcpcs_days_after chemo_post12_enum_days_after chemo_post12_admin_days_after chemo_post12_drug_days_after, replace ignore("NA")

egen minval_chemo = rowmin(chemo_post12_hcpcs_days_after chemo_post12_enum_days_after chemo_post12_admin_days_after chemo_post12_drug_days_after)

*Surgery 
gen surg_days_after_hcpcs = surgdate_hcpcs - dxdate 
gen surg_days_after_icd = surgdate_icd - dxdate 

egen minval_surg = rowmin(surg_days_after_hcpcs surg_days_after_icd)
replace minval_surg = 0 if minval_surg<0

gen timing_of_chemo = ""
replace timing_of_chemo = "Adjuvant" if minval_surg<minval_chemo & minval_surg!=. & minval_chemo!=.
replace timing_of_chemo = "Neoadjuvant" if minval_surg>minval_chemo & minval_surg!=. & minval_chemo!=.
replace timing_of_chemo = "Other" if minval_chemo == minval_surg & minval_surg!=. & minval_chemo!=.

gen diff_chemo_surg = minval_chemo - minval_surg if minval_surg!=. & minval_chemo!=.

gen guideline_concordant_treatment = 0
replace guideline_concordant_treatment = 1 if staging == 1 & ((timing_of_chemo == "Adjuvant" & resection_post12 == 1) | (chemo_post12 == 1 & radiation_post12_any == 1))
replace guideline_concordant_treatment = 1 if staging == 2 & (chemo_post12 == 1 & radiation_post12_any == 1)
replace guideline_concordant_treatment = 1 if staging == 3 & chemo_post12 ==1 

gen no_treatment = 0 
replace no_treatment = 1 if chemo_post12 == 0 & resection_post12 == 0 & radiation_post12 == 0 & immuno_post12_any ==0 