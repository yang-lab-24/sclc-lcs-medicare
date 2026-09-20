clear

do "C:\Users\alp5\Dropbox\STUDIES\SCLC and Lung Cancer Screening\Code\SCLC_LCS_Final_Analysis_FINAL.do" // Import data and clean
do "C:\Users\alp5\Dropbox\STUDIES\SCLC and Lung Cancer Screening\Code\Dummy Variables.do" // Define Dummy Variables 


*Describe Missing Data 
mdesc screened age sex_n race_new yost_acs_2013_2017_quintile income_group urbanrural state acute_mi_binary history_mi_binary chf_binary pvd_binary cvd_binary copd_binary dementia_binary paralysis_binary diabetes_binary diabetes_comp_binary renal_disease_binary mild_liver_disease_binary liver_disease_binary ulcers_binary rheum_disease_binary aids_binary  brain_mri_pre12post3 pet_pre12post3

replace yost_acs_2013_2017_quintile=9 if yost_acs_2013_2017_quintile==. 

set seed 123456
gen u=uniform()
sort u

psmatch2 screened age sex_n i.race_new i.dxyear i.yost_index_clean i.income_group i.state urbanrural acute_mi_binary history_mi_binary chf_binary pvd_binary cvd_binary copd_binary dementia_binary paralysis_binary diabetes_binary diabetes_comp_binary renal_disease_binary mild_liver_disease_binary liver_disease_binary ulcers_binary rheum_disease_binary aids_binary guideline_concordant_treatment, cal(0.01) noreplacement common logit

pstest age sex_n race1 i.race_new dxyear_simple_base i.dx_year_simple_num yost_quintile1 i.yost_index_clean income_group1 i.income_group rural urbanrural state1 i.state i.acute_mi_binary i.history_mi_binary i.chf_binary i.pvd_binary i.cvd_binary i.copd_binary i.dementia_binary i.paralysis_binary i.diabetes_binary i.diabetes_comp_binary i.renal_disease_binary i.mild_liver_disease_binary i.liver_disease_binary i.ulcers_binary i.rheum_disease_binary i.aids_binary guideline_concordant_treatment, treated(screened) both

gen pair = _id if screened ==0 
replace pair = _n1 if screened ==1
bysort pair: egen paircount = count(pair)
drop if paircount !=2

tab screened staging, row chi2 

stset fup_months, fail(vitalstatus == 1)

sts graph, by(screened) xlabel(0(12)36) graphregion(color(white)) ylabel(,nogrid) risktable(,order(1 "Not Screened" 2 "Screened"))  tmax(36) ytitle("Overall Survival (Probability)") xtitle("Time (months)") legend(order(1 "Not Screened" 2 "Screened"))

sts list, by(screened) at (12 24 36)
sts test screened, strata(pair)
