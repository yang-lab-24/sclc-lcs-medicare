clear

*--------------------------------------------------------*
* Import and clean original data
*--------------------------------------------------------*
    do "C:\Users\alp5\Dropbox\STUDIES\SCLC and Lung Cancer Screening\Code\SCLC_LCS_Final_Analysis_FINAL.do"
    do "C:\Users\alp5\Dropbox\STUDIES\SCLC and Lung Cancer Screening\Code\Dummy Variables.do"
	
******************************
* Section IV: Propensity Score Matched Analyses 
******************************

*Pre-match checks 
mdesc screened age sex_n race_new yost_acs_2013_2017_quintile income_group urbanrural state acute_mi_binary history_mi_binary chf_binary pvd_binary cvd_binary copd_binary dementia_binary paralysis_binary diabetes_binary diabetes_comp_binary renal_disease_binary mild_liver_disease_binary liver_disease_binary ulcers_binary rheum_disease_binary aids_binary

*Propensity Score-Match 
set seed 123456
gen u=uniform()
sort u

psmatch2 screened age sex_n i.race_new i.dxyear i.yost_index_clean i.income_group i.state urbanrural acute_mi_binary history_mi_binary chf_binary pvd_binary cvd_binary copd_binary dementia_binary paralysis_binary diabetes_binary diabetes_comp_binary renal_disease_binary mild_liver_disease_binary liver_disease_binary ulcers_binary rheum_disease_binary aids_binary, cal(0.01) noreplacement common logit

*This code creates a variable with an "ID" for each matched pair
gen pair = _id if screened ==0 
replace pair = _n1 if screened ==1

*This code creates a variable called paircount which counts the number of patients who have the same match "ID." Patients with a paircount !=2 were not included in the match. 
bysort pair: egen paircount = count(pair)

*Visualize the distribution of propensity scores in the unmatched and matched cohorts 
twoway (kdensity _pscore if _treated==1) (kdensity _pscore if _treated==0, lpattern(dash)), legend( label( 1 "treated") label( 2 "control" )) xtitle("propensity scores BEFORE matching") 

gen match=_n1
replace match=_id if match==.
duplicates tag match, gen(dup)
twoway (kdensity _pscore if _treated==1) (kdensity _pscore if _treated==0 & dup>0, lpattern(dash)), legend( label( 1 "treated") label( 2 "control" )) 

pstest age sex_n race1 i.race_new dxyear_simple_base i.dx_year_simple_num yost_quintile1 i.yost_index_clean income_group1 i.income_group rural urbanrural state1 i.state i.acute_mi_binary i.history_mi_binary i.chf_binary i.pvd_binary i.cvd_binary i.copd_binary i.dementia_binary i.paralysis_binary i.diabetes_binary i.diabetes_comp_binary i.renal_disease_binary i.mild_liver_disease_binary i.liver_disease_binary i.ulcers_binary i.rheum_disease_binary i.aids_binary, treated(screened) both

*Remove patients with a paircount not equal to 2
drop if paircount !=2

********Analyses********

*Stage Distribution
tab staging screened, col 
tab ajcc_staging screened if missing_ajcc!=1, col
tab ajcc_staging_simple screened if missing_ajcc!=1, col missing


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
table1, by(screened) vars(age conts\ sex_n cat\ race_new cat\ dx_year_simple cat\ acute_mi_binary bin\ history_mi_binary bin\ chf_binary bin\ pvd_binary bin\ cvd_binary bin\ copd_binary bin\ dementia_binary bin\ paralysis_binary bin\ diabetes_binary bin\ diabetes_comp_binary bin\ renal_disease_binary bin\ mild_liver_disease_binary bin\ liver_disease_binary bin\ ulcers_binary bin\ rheum_disease_binary bin\ aids_binary bin\ income_group cat\ yost_index_clean cat\ rural_urban cat\ state cat\ brain_mri_pre12post3 cat\ pet_pre12post3 cat\ guideline_concordant_treatment cat), if sex_n!=5, format(%2.1f) onecol missing saving(tables1_baselinecharacteristics_psmatch.xls, replace)
