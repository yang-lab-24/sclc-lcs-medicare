*Title: Outcomes Analysis of Screen-Detected SCLC
*Dataset: SEER-Medicare
*Date Last Revised: 2/17/2025	

clear

import delimited "C:\Users\sonee\Downloads\finalSCLCstudyDATA.csv"


*Notes:
/*
The code is structured into sections, including:

1. Section I: Defining Variables
2. Section II: Defining Survival
3. Section III: Evaluating the Proportion and Absolute Number of Screen-Detected SCLC
4. Section IV: Age and Sex Adjusted Propensity Score Matched Cohorts
5. Section V: Survival Outcomes
6. Section VI: Baseline Characteristics

*/


******************************
* Section I: Defining Variables
******************************

******Tumor Characteristics Variables******
*Stage 
gen staging = . 
replace staging = combined_summary_stage_2004

label define stagelabel 0 "In Situ" 1 "Localized" 2 "Regional" 7 "Distant" 9 "Unknown" 
label values staging stagelabel 


*tumor size
gen tumor_size_ = size
replace tumor_size_ = "" if tumor_size_ == "NA"
destring tumor_size_, gen(tumor_size_clean)


gen tumorlocation_ = ""
replace tumorlocation_ = "Right upper lobe" if primary_site == "C341" & laterality == 1
replace tumorlocation_ = "Right middle lobe" if primary_site == "C342" & laterality == 1
replace tumorlocation_ = "Right lower lobe" if primary_site == "C343" & laterality == 1
replace tumorlocation_ = "Left upper lobe" if primary_site == "C341" & laterality == 2
replace tumorlocation_ = "Left lower lobe" if primary_site == "C343" & laterality == 2
replace tumorlocation_ = "Bilateral" if laterality == 4
encode tumorlocation_, gen(tumorlocation_clean)

encode sex, gen(sex_n)

encode race, gen(race_new)

******Screening Variables******
gen unscreened_new = .
replace unscreened_new = 1 if screened == 0 & age >= 66.5 & age <= 79

gen grouped_screening = .
replace grouped_screening = 0 if unscreened_new == 1
replace grouped_screening = 1 if screened == 1


******************************
* Section II: Defining Survival & Inclusion Criteria
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


******************************
* Section III: Evaluating the Proportion and Absolute Number of Screen-Detected SCLC Overall and By Stage
******************************

tab grouped_screening
tab year_of_diagnosis grouped_screening, row 
tab staging grouped_screening, col


******************************
* Section IV: Evaluating the Receipt of Immunotherapy by Immunotherapy
******************************
tab nivo12months dxyear, col
tab pembro12months dxyear, col
tab atezo12months dxyear, col
tab durva12months dxyear, col
tab lurb12months dxyear, col
tab immuno12months dxyear, col

******************************
* Section IV: Age and Sex Adjusted Propensity Score Matched Cohort
******************************


*Overall Survival
clear

import delimited "C:\Users\sonee\Downloads\finalSCLCstudyDATA.csv"


set seed 123456
gen u=uniform()
sort u

psmatch2 screened age sex_n dxyear, cal(0.01) noreplacement common logit

*This code creates a variable with an "ID" for each matched pair
gen pair = _id if grouped_screening ==0 
replace pair = _n1 if grouped_screening ==1

*This code creates a variable called paircount which counts the number of patients who have the same match "ID." Patients with a paircount !=2 were not included in the match. 
bysort pair: egen paircount = count(pair)

*Remove patients with a paircount not equal to 2
drop if paircount !=2

*Save your paired dataset (this part of the code is not absolutely necessary)
save paired, replace

**Survival Analyses Accounting for the Paired Nature of the Dataset 
*Set survival time
stset fup_months, fail(vitalstatus==1) id(patient_id)

*Kaplan Meier Survival Code
sts graph, by(grouped_screening) xlabel(0(12)36) graphregion(color(white)) ylabel(,nogrid) risktable(,order(1 "Not Screened" 2 "Screened"))  tmax(36) ytitle("Overall Survival (Probability)") xtitle("Time (months)") legend(order(1 "Not Screened" 2 "Screened"))

*Stratified Log-rank Test
sts test grouped_screening, strata(pair) 

*List the survivor function at 5-years 
sts list, by (grouped_screening) at (12 36 60)

*Propensity Score Test

clear

import delimited "C:\Users\sonee\Downloads\finalSCLCstudyDATA.csv"

set seed 123456
gen u=uniform()
sort u

psmatch2 grouped_screening age sex_n dxyear, cal(0.01) noreplacement common logit
pstest age sex_n, treated(grouped_screening) both, if lc_specific_analysis == 1

*Baseline Characteristics of the Cohort
table1, by(grouped_screening) vars(age conts\ sex_n cat\ race_new cat\ medianhouseholdincomeinflationa conts\ yost_acs_2013_2017_quintile cat\ seer_registry cat\ staging cat\ tumorlocation_clean cat\ tumor_size_clean conts), if sex_n!=5, format(%2.1f) onecol missing saving(tables_psm_overallcohort_sclc_lcs_final.xls, replace)


******************************
* Section V: Survival Outcomes
******************************


***Overall Survival 
stset fup_months, fail(vitalstatus==1)

*Three-Year Overall Survival Screened vs. Not Screened 
sts graph, by(grouped_screening) xlabel(0(12)36) graphregion(color(white)) ylabel(,nogrid) risktable(,order(1 "Not Screened" 2 "Screened"))  tmax(36) ytitle("Overall Survival (Probability)") xtitle("Time (months)") legend(order(1 "Not Screened" 2 "Screened"))

sts list if sclc==1, by(grouped_screening) at(12 24 36)

*Three-Year Overall Survival By Stage Among Patients with Screen-Detected SCLC 
sts graph, by(staging) xlabel(0(12)36) graphregion(color(white)) ylabel(,nogrid) risktable(,order(1 "Localized" 2 "Regional" 3 "Distant"))  tmax(36) ytitle("Overall Survival (Probability)") xtitle("Time (months)") legend(order(1 "Localized" 2 "Regional" 3 "Distant")), if grouped_screening==1 & staging<9

sts list if sclc==1 & grouped_screening==1 & staging<9, by(staging) at(12 24 36)


*survival table - overall survival_months
stset fup_months, fail(vitalstatus==1)
sts list, by(grouped_screening staging) at(12 24 36)



******************************
* Section VI: Baseline Characteristics
******************************

*Baseline Characteristics Table
table1, by(grouped_screening) vars(age conts\ sex_n cat\ race_new cat\ medianhouseholdincomeinflationa conts\ yost_acs_2013_2017_quintile cat\ seer_registry cat\ staging cat\ tumorlocation_clean cat\ tumor_size_clean conts), if sex_n!=5, format(%2.1f) onecol missing saving(tables1_baselinecharacteristics_overall_sclc_lcs_final.xls, replace)