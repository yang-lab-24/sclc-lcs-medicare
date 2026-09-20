**********************************************************
*Title: Low-dose Computed Tomography Screening is Associated with Earlier Stage Diagnoses and Improved Survival for Small Cell Lung Cancer 
*Dataset: SEER-Medicare
*Date Last Revised: 9/5/2026	
**********************************************************

**Load in data (see R code for how dataset was cleaned)
clear
import delimited "C:\Users\alp5\Dropbox\STUDIES\SCLC and Lung Cancer Screening\Data\final_analysis_cohort_info_9.5.2026.csv"

do "C:\Users\alp5\Dropbox\STUDIES\SCLC and Lung Cancer Screening\Code\Treatment Code JAMA Onc.do" //load in treatment definitions

*Drop pre-cleaned variables already in the new CSV so the do-file can recreate them
*capture drop ajcc_staging missing_ajcc sex_n no_cod_info fup_months tumor_size_clean screened 


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

*Renaming Relevant Variables 
rename year_of_diagnosis dxyear

gen dx_year_simple = "2015-2016" if inlist(dxyear, 2015, 2016)
replace dx_year_simple = "2017" if dxyear == 2017
replace dx_year_simple = "2018" if dxyear == 2018
replace dx_year_simple = "2019" if dxyear == 2019
replace dx_year_simple = "2020" if dxyear == 2020 
replace dx_year_simple = "2021" if dxyear == 2021 

encode dx_year_simple, gen (dx_year_simple_num) 

*AJCC Staging Variable (Only available for patients in select registrations) 
gen ajcc_staging = stage_simple 
replace ajcc_staging = "Missing/Unknown" if strpos(stage_simple, "Missing Data") | strpos(stage_simple, "Data available")

gen ajcc_staging_simple = ""
replace ajcc_staging_simple = "Stage I" if inlist(ajcc_staging, "Stage IA", "Stage IB")
replace ajcc_staging_simple = "Stage II" if inlist(ajcc_staging, "Stage IIA", "Stage IIB")
replace ajcc_staging_simple = ajcc_staging if strpos(ajcc_staging, "III")
replace ajcc_staging_simple = ajcc_staging if ajcc_staging == "Stage IV"

*Missing AJCC Staging Data 
gen missing_ajcc = 0 
replace missing_ajcc = 1 if seer_registry>=61 & seer_registry<=66


*Screened 
gen screened = 0 
replace screened = 1 if scan_type == "Screened"

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

**Tabulate all comorbidities to assess frequency
tab1 acute_mi_binary history_mi_binary chf_binary pvd_binary cvd_binary copd_binary dementia_binary paralysis_binary diabetes_binary diabetes_comp_binary renal_disease_binary mild_liver_disease_binary liver_disease_binary ulcers_binary rheum_disease_binary aids_binary

*Median Household Income - Consolidate 16 SEER categories into 4 quartiles
gen income_group = .
replace income_group = 1 if medianhouseholdincomeinflation >= 1 & medianhouseholdincomeinflation <= 4
replace income_group = 2 if medianhouseholdincomeinflation >= 5 & medianhouseholdincomeinflation <= 8
replace income_group = 3 if medianhouseholdincomeinflation >= 9 & medianhouseholdincomeinflation <= 12
replace income_group = 4 if medianhouseholdincomeinflation >= 13 & medianhouseholdincomeinflation <= 16

label define income_grp 1 "< $55,000" 2 "$55,000 - $74,999" 3 "$75,000 - $94,999" 4 ">= $95,000"
label values income_group income_grp

tab income_group screened, col

*Yost Index clean 
gen yost_index_clean = yost_acs_2013_2017_quintile
replace yost_index_clean = 6 if yost_acs_2013_2017_quintile == .

******************************
* Section II: Defining Survival Information 
******************************
destring survival_months, replace force
replace survival_months = 0.01 if survival_months == 0

gen fupdate = death_date
replace fupdate = "2022-12-31" if death_date == "NA"


gen death_date_new = date(fupdate, "YMD")
format death_date_new %td

gen dxdate_new = date(dx_date, "YMD")
format dxdate_new %td

gen fup_days = death_date_new - dxdate_new
summ fup_days
replace fup_days = 0.01 if fup_days <= 0
gen fup_months = fup_days / 30
summ fup_months 

*Overall Survival
gen vitalstatus = 1 if death_date != "NA"
replace vitalstatus = 0 if death_date == "NA"

*Cause of death information available 
gen no_cod_info = 1 if seer_registry==61 | seer_registry==63

*Lung Cancer Death 
gen lungcancerdeath = 0 
replace lungcancerdeath = 1 if seercausespecificdeathclassific == 1

*Other Death 
gen otherdeath = 0
replace otherdeath = 1 if vitalstatus==1  & seerothercauseofdeathclassifica  == 1


******************************
* Section III: Unadjusted Analyses  
******************************

*Overall Numbers
tab screened 
tab dxyear screened, row 

*Stage Distribution (Figure 1)
tab staging screened, col 
tab ajcc_staging screened if missing_ajcc!=1, col 
tab ajcc_staging_simple screened if missing_ajcc!=1, col missing

*Overall Survival (Figure 2)
stset fup_months, fail(vitalstatus==1) 
sts graph, by(screened) xlabel(0(12)36) graphregion(color(white)) ylabel(,nogrid) risktable(,order(1 "Unscreened" 2 "Screened"))  tmax(36) ytitle("Overall Survival (Probability)") xtitle("Time (months)") legend(order(1 "Unscreened" 2 "Screened"))
sts test screened 
sts list, by(screened) at (12 24 36)


*Evaluating the time between the last LDCT scan and the SCLC diagnosis (Supplementary Figure 2)
replace time_diff = "" if time_diff=="NA"
destring time_diff, replace force 
histogram time_diff if screened==1, xtitle("Days Between the LDCT Screening Exam and Diagnosis Date") ytitle("Proportion")

*Baseline Characteristics 
table1, by(screened) vars(age conts\ sex_n cat\ race_new cat\ dx_year_simple cat\ acute_mi_binary bin\ history_mi_binary bin\ chf_binary bin\ pvd_binary bin\ cvd_binary bin\ copd_binary bin\ dementia_binary bin\ paralysis_binary bin\ diabetes_binary bin\ diabetes_comp_binary bin\ renal_disease_binary bin\ mild_liver_disease_binary bin\ liver_disease_binary bin\ ulcers_binary bin\ rheum_disease_binary bin\ aids_binary bin\ income_group cat\ yost_index_clean cat\ rural_urban cat\ state cat\ brain_mri_pre12post3 cat\ pet_pre12post3 cat\ guideline_concordant_treatment cat), if sex_n!=5, format(%2.1f) onecol missing saving(tables1_baselinecharacteristics_overall_sclc_lcs_final.xls, replace)

