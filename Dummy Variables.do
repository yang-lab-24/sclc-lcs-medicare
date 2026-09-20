*Clean Staging Variables 
gen tnm_n_clean = tnm_n if inlist(tnm_n, "0", "1", "2", "3")
destring tnm_n_clean, replace force

replace t_stage_8 = "" if t_stage_8 == "NA"
encode t_stage_8, gen (t_stage_8_num)


*Define Dummy Variables 
gen race1 = 0
replace race1 = 1 if race_new == 1

gen sex1 = 0
replace sex1 = 1 if sex_n == 1

gen dxyear2015 = 0
replace dxyear2015 = 1 if dxyear == 2015

gen dxyear_simple_base  = 0
replace dxyear_simple_base = 1 if dx_year_simple == "2015-2016"

gen income_group1 = 0
replace income_group1 = 1 if income_group == 1

gen yost_quintile1 = 0
replace yost_quintile1 = 1 if yost_acs_2013_2017_quintile == 1

gen rural = 0
replace rural = 1 if rural_urban == "Rural"

gen state1 = 0
replace state1 = 1 if state == 6


gen tumorlocationbilateral = 0
replace tumorlocationbilateral = 1 if tumorlocation_clean == 1

gen no_acute_mi = 0
replace no_acute_mi = 1 if acute_mi_binary == 0

gen no_history_mi = 0
replace no_history_mi = 1 if history_mi_binary == 0

gen no_chf = 0
replace no_chf = 1 if chf_binary == 0

gen no_pvd = 0
replace no_pvd = 1 if pvd_binary == 0

gen no_cvd = 0
replace no_cvd = 1 if cvd_binary == 0

gen no_copd = 0
replace no_copd = 1 if copd_binary == 0

gen no_dementia = 0
replace no_dementia = 1 if dementia_binary == 0

gen no_paralysis = 0
replace no_paralysis = 1 if paralysis_binary == 0

gen no_diabetes = 0
replace no_diabetes = 1 if diabetes_binary == 0

gen no_diabetes_comp = 0
replace no_diabetes_comp = 1 if diabetes_comp_binary == 0

gen no_renal_disease = 0
replace no_renal_disease = 1 if renal_disease_binary == 0

gen no_mild_liver_disease = 0
replace no_mild_liver_disease = 1 if mild_liver_disease_binary == 0

gen no_liver_disease = 0
replace no_liver_disease = 1 if liver_disease_binary == 0

gen no_ulcers = 0
replace no_ulcers = 1 if ulcers_binary == 0

gen no_rheum_disease = 0
replace no_rheum_disease = 1 if rheum_disease_binary == 0

gen no_aids = 0
replace no_aids = 1 if aids_binary == 0

gen t_stage_8_0 = 0
replace t_stage_8_0 = 1 if t_stage_8_num == 1 

gen tnm_n_clean_0 = 0 
replace tnm_n_clean_0 = 1 if tnm_n_clean_0 == 0 

gen brain_mri_pre12post3_none = 0 
replace brain_mri_pre12post3_none = 1 if brain_mri_pre12post3 == 0 

gen pet_pre12post3_none = 0 
replace pet_pre12post3_none = 1 if pet_pre12post3 == 0 