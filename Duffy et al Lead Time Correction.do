clear

do "C:\Users\alp5\Dropbox\STUDIES\SCLC and Lung Cancer Screening\Code\SCLC_LCS_Final_Analysis_FINAL.do" // Overall Study Code 

do "C:\Users\alp5\Downloads\Propensity Score Match Overall.do" // Propensity Score-Matched Analysis 

* Assumed mean sojourn time = 3, 6, 12, 18 months

local mst = 15

gen k = 1/`mst'

gen expected_lt = 0

* Screened patients who died
replace expected_lt = ///
    (1-exp(-k*fup_months)-k*fup_months*exp(-k*fup_months)) / ///
    (k*(1-exp(-k*fup_months))) ///
    if screened==1 & vital_status==1

* Screened patients who were censored/alive
replace expected_lt = ///
    (1-exp(-k*fup_months))/k ///
    if screened==1 & vital_status==0

* Adjust only screened patients
gen surv_adjusted = fup_months
replace surv_adjusted = fup_months-expected_lt if screened==1

* Check the correction
summ fup_months expected_lt surv_adjusted if screened==1
assert expected_lt <= fup_months if screened==1 & vital_status==1
assert surv_adjusted >= 0 if !missing(surv_adjusted)

kdensity surv_adjusted if screened==1, ///
    addplot(kdensity fup_months if screened==1)

stset surv_adjusted, failure(vital_status==1)

sts graph, by(screened) tmax(36)
sts list, by(screened) at(12 24 36)
sts test screened
sts test screened, strata(pair)