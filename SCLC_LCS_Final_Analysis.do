**********************************************************
*Title: Low-dose Computed Tomography Screening is Associated with Earlier Stage Diagnoses and Improved Survival for Small Cell Lung Cancer 
*Dataset: SEER-Medicare
*Date Last Revised: 3/8/2025	
**********************************************************

**Load in data (see R code for how dataset was cleaned)
clear
import delimited "C:\Users\cjy2\Dropbox\STUDIES\SCLC and Lung Cancer Screening\Data\SCLC.screening.data.final.csv"


*Notes:
/*
The code is structured into sections, including:

1. Section I: Defining Variables
2. Section II: Defining Inclusion Criteria 
3. Section III: Unadjusted Outcomes 
4. Section IV: Propensity Score Matched Analysis 
5. Section V: Assessing Standardized Differences of the Propensity Score Match 
*/


******************************
* Section I: Defining Variables
******************************

******Tumor Characteristics Variables******

*AJCC Staging Variable (Only available for patients in select registrations) 
gen ajcc_staging = stage_simple 
replace ajcc_staging = "Missing/Unknown" if strpos(stage_simple, "Missing Data") | strpos(stage_simple, "Data available")

*Missing AJCC Staging Data 
gen missing_ajcc = 0 
replace missing_ajcc = 1 if seer_registry>=61 & seer_registry<=66


*Tumor Size
gen tumor_size_ = size
replace tumor_size_ = "" if tumor_size_ == "NA"
destring tumor_size_, gen(tumor_size_clean)


*Tumor Location
gen tumorlocation_ = ""
replace tumorlocation_ = "Right upper lobe" if primary_site == "C341" & laterality == 1
replace tumorlocation_ = "Right middle lobe" if primary_site == "C342" & laterality == 1
replace tumorlocation_ = "Right lower lobe" if primary_site == "C343" & laterality == 1
replace tumorlocation_ = "Left upper lobe" if primary_site == "C341" & laterality == 2
replace tumorlocation_ = "Left lower lobe" if primary_site == "C343" & laterality == 2
replace tumorlocation_ = "Bilateral" if laterality == 4
encode tumorlocation_, gen(tumorlocation_clean)

*Other Demographics
encode sex, gen(sex_n)
encode race, gen(race_new)
encode rural_urban, gen(urbanrural)

*Comorbidities 
gen asthma_new = .
replace asthma_new = 0 if asthma_pre_count == 0
replace asthma_new = 1 if asthma_pre_count >= 1

gen ischemic_heart_disease_new = .
replace ischemic_heart_disease_new = 0 if ischemic_heart_disease_pre_count == 0
replace ischemic_heart_disease_new = 1 if ischemic_heart_disease_pre_count >= 1

gen cvd_new = .
replace cvd_new = 0 if cvd_pre_count == 0
replace cvd_new = 1 if cvd_pre_count >= 1

gen mi_new = .
replace mi_new = 0 if mi_pre_count == 0
replace mi_new = 1 if mi_pre_count >= 1

gen copd_new = .
replace copd_new = 0 if copd_and_allied_conditions_pre_c == 0
replace copd_new = 1 if copd_and_allied_conditions_pre_c >= 1

gen chf_new = .
replace chf_new = 0 if chf_pre_count == 0
replace chf_new = 1 if chf_pre_count >= 1

gen pvd_new = .
replace pvd_new = 0 if pvd_pre_count == 0
replace pvd_new = 1 if pvd_pre_count >= 1

gen cerebrovascular_disease_new = .
replace cerebrovascular_disease_new = 0 if cerebrovascular_disease_pre_coun == 0
replace cerebrovascular_disease_new = 1 if cerebrovascular_disease_pre_coun >= 1

gen diab_uc_new = .
replace diab_uc_new = 0 if diab_uc_pre_count == 0
replace diab_uc_new = 1 if diab_uc_pre_count >= 1

gen interstitial_lung_new = .
replace interstitial_lung_new = 0 if interstitial_lung_pre_count == 0
replace interstitial_lung_new = 1 if interstitial_lung_pre_count >= 1

gen o2_new = .
replace o2_new = 0 if o2_pre_count == 0
replace o2_new = 1 if o2_pre_count >= 1

gen kidney_diseases_new = .
replace kidney_diseases_new = 0 if kidney_diseases_pre_count == 0
replace kidney_diseases_new = 1 if kidney_diseases_pre_count >= 1

**Tabulate all comorbidities to assess frequency
tab1 asthma_new ischemic_heart_disease_new cvd_new mi_new copd_new chf_new pvd_new cerebrovascular_disease_new diab_uc_new interstitial_lung_new o2_new kidney_diseases_new

******************************
* Section II: Defining Survival Information 
******************************
replace survival_months = 0.01 if survival_months == 0

gen fupdate = deathdate
replace fupdate = "2020-12-31" if deathdate == "NA"


gen death_date_new = date(fupdate, "YMD")
format death_date_new %td

gen dxdate_new = date(dxdate, "YMD")
format dxdate_new %td

gen fup_days = death_date_new - dxdate_new
summ fup_days
replace fup_days = 0.01 if fup_days <= 0
gen fup_months = fup_days / 30
summ fup_months 

*Overall Survival
gen vitalstatus = 1 if deathdate != "NA"
replace vitalstatus = 0 if deathdate == "NA"

/*
******************************
* Section III: Unadjusted Analyses  
******************************

*Overall Numbers
tab screened 
tab dxyear screened, row 

*Stage Distribution (Figure 1)
tab staging screened, col 
tab ajcc_staging screened if missing_ajcc!=1, col 

*Overall Survival (Figure 2)
stset fup_months, fail(vitalstatus==1) 

sts graph, by(screened) xlabel(0(12)36) graphregion(color(white)) ylabel(,nogrid) risktable(,order(1 "Unscreened" 2 "Screened"))  tmax(36) ytitle("Overall Survival (Probability)") xtitle("Time (months)") legend(order(1 "Unscreened" 2 "Screened"))

sts test screened 

sts list, by(screened) at (12 24 36)

*Evaluating the time between the last LDCT scan and the SCLC diagnosis (Supplementary Figure 2)
replace timediff = "" if timediff=="NA"
destring timediff, replace force 
histogram timediff if screened==1, xtitle("Days Between the LDCT Screening Exam and Diagnosis Date") ytitle("Proportion")

*Baseline Characteristics 
table1, by(screened) vars(age conts\ sex_n cat\ race_recode_w_b_ai_api cat\ dxyear cat\ asthma_new bin\ ischemic_heart_disease_new bin\ cvd_new bin\ mi_new bin\ copd_new bin\ chf_new bin\ pvd_new bin\ cerebrovascular_disease_new bin\ diab_uc_new bin\ interstitial_lung_new bin\ o2_new bin\ kidney_diseases_new bin \ medianhouseholdincomeinflationa cat\ yost_acs_2013_2017_quintile cat\ rural_urban cat\ state cat), if sex_n!=5, format(%2.1f) onecol missing saving(tables1_baselinecharacteristics_overall_sclc_lcs_final.xls, replace)


******************************
* Section IV: Propensity Score Matched Analyses 
******************************

*Pre-match checks 
mdesc screened age sex_n race_recode_w_b_ai_api dxyear yost_acs_2013_2017_quintile medianhouseholdincomeinflationa urbanrural state asthma_new ischemic_heart_disease_new cvd_new mi_new copd_new chf_new pvd_new cerebrovascular_disease_new diab_uc_new interstitial_lung_new o2_new kidney_diseases_new

replace yost_acs_2013_2017_quintile=9 if yost_acs_2013_2017_quintile==. 
// run this line to ensure that patients with missing data for this variable are not excluded 


*Propensity Score-Match 
set seed 123456
gen u=uniform()
sort u

psmatch2 screened age sex_n i.race_recode_w_b_ai_api i.dxyear i.yost_acs_2013_2017_quintile i.medianhouseholdincomeinflationa i.state urbanrural asthma_new ischemic_heart_disease_new cvd_new mi_new copd_new chf_new pvd_new cerebrovascular_disease_new diab_uc_new interstitial_lung_new o2_new kidney_diseases_new, cal(0.01) noreplacement common logit

*This code creates a variable with an "ID" for each matched pair
gen pair = _id if screened ==0 
replace pair = _n1 if screened ==1

*This code creates a variable called paircount which counts the number of patients who have the same match "ID." Patients with a paircount !=2 were not included in the match. 
bysort pair: egen paircount = count(pair)

*Remove patients with a paircount not equal to 2
drop if paircount !=2


********Analyses********

*Stage Distribution
tab staging screened, col 
tab ajcc_staging screened if missing_ajcc!=1, col


*Overall Survival 
*Set survival time
stset fup_months, fail(vitalstatus==1) id(patient_id)

*Kaplan Meier Survival Code
sts graph, by(screened) xlabel(0(12)36) graphregion(color(white)) ylabel(,nogrid) risktable(,order(1 "Not Screened" 2 "Screened"))  tmax(36) ytitle("Overall Survival (Probability)") xtitle("Time (months)") legend(order(1 "Not Screened" 2 "Screened"))

*Stratified Log-rank Test
sts test screened, strata(pair) 

*List the survivor function at 5-years 
sts list, by (screened) at (12 24 36)

*Baseline Characteristics 
table1, by(screened) vars(age conts\ sex_n cat\ race_recode_w_b_ai_api cat\ dxyear cat\ asthma_new bin\ ischemic_heart_disease_new bin\ cvd_new bin\ mi_new bin\ copd_new bin\ chf_new bin\ pvd_new bin\ cerebrovascular_disease_new bin\ diab_uc_new bin\ interstitial_lung_new bin\ o2_new bin\ kidney_diseases_new bin \ medianhouseholdincomeinflationa cat\ yost_acs_2013_2017_quintile cat\ rural_urban cat\ state cat), if sex_n!=5, format(%2.1f) onecol missing saving(tables1_baselinecharacteristics_psmatch.xls, replace)


******************************
* Section V: Assessing Balance in the Propensity Score-Matched Cohort  
******************************

set seed 123456
gen u=uniform()
sort u

replace yost_acs_2013_2017_quintile=9 if yost_acs_2013_2017_quintile==.

psmatch2 screened age sex_n i.race_recode_w_b_ai_api i.dxyear i.yost_acs_2013_2017_quintile i.medianhouseholdincomeinflationa i.urbanrural i.state i.asthma_new i.ischemic_heart_disease_new i.cvd_new i.mi_new i.copd_new i.chf_new i.pvd_new i.cerebrovascular_disease_new i.diab_uc_new i.interstitial_lung_new i.o2_new i.kidney_diseases_new, cal(0.01) noreplacement common logit


*Visualize the distribution of propensity scores in the unmatched and matched cohorts 
twoway (kdensity _pscore if _treated==1) (kdensity _pscore if _treated==0, lpattern(dash)), legend( label( 1 "treated") label( 2 "control" )) xtitle("propensity scores BEFORE matching") 

gen match=_n1
replace match=_id if match==.
duplicates tag match, gen(dup)
twoway (kdensity _pscore if _treated==1) (kdensity _pscore if _treated==0 & dup>0, lpattern(dash)), legend( label( 1 "treated") label( 2 "control" )) 


*Define Dummy Variables 
gen race1 = 0
replace race1 = 1 if race_recode_w_b_ai_api== 1

gen sex1 = 0
replace sex1 = 1 if sex_n == 1

gen dxyear2015 = 0
replace dxyear2015 = 1 if dxyear == 2015

gen medianhouseholdincome1 = 0
replace medianhouseholdincome1 = 1 if medianhouseholdincomeinflationa == 1

gen yost_quintile1 = 0
replace yost_quintile1 = 1 if yost_acs_2013_2017_quintile == 1

gen rural = 0
replace rural = 1 if rural_urban == "Rural"

gen state1 = 0
replace state1 = 1 if state == 6

gen staging0local = 0
replace staging0local = 1 if staging0 == 1

gen tumorlocationbilateral = 0
replace tumorlocationbilateral = 1 if tumorlocation_clean == 1

gen no_asthma = 0
replace no_asthma = 1 if asthma_new == 0

gen no_ischemic_heart_disease = 0
replace no_ischemic_heart_disease = 1 if ischemic_heart_disease_new == 0

gen no_cvd = 0
replace no_cvd = 1 if cvd_new == 0

gen no_mi = 0
replace no_mi = 1 if mi_new == 0

gen no_copd = 0
replace no_copd = 1 if copd_new == 0

gen no_chf = 0
replace no_chf = 1 if chf_new == 0

gen no_pvd = 0
replace no_pvd = 1 if pvd_new == 0

gen no_cerebrovascular_disease = 0
replace no_cerebrovascular_disease = 1 if cerebrovascular_disease_new == 0

gen no_diab_uc = 0
replace no_diab_uc = 1 if diab_uc_new == 0

gen no_interstitial_lung = 0
replace no_interstitial_lung = 1 if interstitial_lung_new == 0

gen no_o2 = 0
replace no_o2 = 1 if o2_new == 0

gen no_kidney_diseases = 0
replace no_kidney_diseases = 1 if kidney_diseases_new == 0

tab race_recode_w_b_ai_api, gen(race_new)


pstest age sex_n race1 i.race_recode_w_b_ai_api dxyear2015 i.dxyear yost_quintile1 i.yost_acs_2013_2017_quintile medianhouseholdincome1 i.medianhouseholdincomeinflationa rural urbanrural state1 i.state i.asthma_new i.ischemic_heart_disease_new i.cvd_new i.mi_new i.copd_new i.chf_new i.pvd_new i.cerebrovascular_disease_new i.diab_uc_new i.interstitial_lung_new i.o2_new i.kidney_diseases_new, treated(screened) both
