clear

*============================================================*
* Loop through stage groups
*============================================================*

foreach stagegroup in stage4  { //stage3 stage4

    clear

    *--------------------------------------------------------*
    * Import and clean original data
    *--------------------------------------------------------*
    do "C:\Users\alp5\Dropbox\STUDIES\SCLC and Lung Cancer Screening\Code\SCLC_LCS_Final_Analysis_FINAL.do"S
    do "C:\Users\alp5\Dropbox\STUDIES\SCLC and Lung Cancer Screening\Code\Dummy Variables.do"

    *--------------------------------------------------------*
    * Restrict to stage group
    *--------------------------------------------------------*
    if "`stagegroup'" == "stage12" {
        keep if inlist(ajcc_staging, ///
            "Stage IA", "Stage IB", "Stage IIA", "Stage IIB")
        local stagelabel "Stage I-II"
    }

    if "`stagegroup'" == "stage3" {
        keep if inlist(ajcc_staging, ///
            "Stage IIIA", "Stage IIIB", "Stage IIIC")
        local stagelabel "Stage III"
    }

    if "`stagegroup'" == "stage4" {
        keep if ajcc_staging == "Stage IV"
        local stagelabel "Stage IV"
    }

    display "===================================================="
    display " ANALYSIS FOR `stagelabel'"
    display "===================================================="

    *--------------------------------------------------------*
    * Unmatched survival
    *--------------------------------------------------------*
    stset fup_months, fail(vitalstatus == 1)

    sts graph, by(screened) ///
        xlabel(0(12)36) ///
        graphregion(color(white)) ///
        ylabel(,nogrid) ///
        risktable(,order(1 "Not Screened" 2 "Screened")) ///
        tmax(36) ///
        ytitle("Overall Survival (Probability)") ///
        xtitle("Time (months)") ///
        legend(order(1 "Not Screened" 2 "Screened")) ///
        title("`stagelabel'")

    tab screened guideline_concordant_treatment, row

    *--------------------------------------------------------*
    * Missingness
    *--------------------------------------------------------*
    mdesc screened age sex_n race_new ///
        yost_acs_2013_2017_quintile income_group ///
        urbanrural state acute_mi_binary history_mi_binary ///
        chf_binary pvd_binary cvd_binary copd_binary ///
        dementia_binary paralysis_binary diabetes_binary ///
        diabetes_comp_binary renal_disease_binary ///
        mild_liver_disease_binary liver_disease_binary ///
        ulcers_binary rheum_disease_binary aids_binary ///
        t_stage_8_num tnm_n_clean tumor_size_clean ///
        brain_mri_pre12post3 pet_pre12post3 resection_post12 ///
		chemo_post12 radiation_post12_any

    *--------------------------------------------------------*
    * Propensity score matching
    *--------------------------------------------------------*
    replace yost_acs_2013_2017_quintile = 6 ///
        if yost_acs_2013_2017_quintile == .

    set seed 123456
    gen u = uniform()
    sort u

    psmatch2 screened ///
        age sex_n ///
        i.race_new ///
        i.dxyear ///
        i.yost_acs_2013_2017_quintile ///
        i.income_group ///
        i.state ///
        urbanrural ///
        acute_mi_binary ///
        history_mi_binary ///
        chf_binary ///
        pvd_binary ///
        cvd_binary ///
        copd_binary ///
        dementia_binary ///
        paralysis_binary ///
        diabetes_binary ///
        diabetes_comp_binary ///
        renal_disease_binary ///
        mild_liver_disease_binary ///
        liver_disease_binary ///
        ulcers_binary ///
        rheum_disease_binary ///
        aids_binary ///
        i.t_stage_8_num ///
        i.tnm_n_clean ///
        tumor_size_clean ///
        brain_mri_pre12post3 ///
        pet_pre12post3, ///
        cal(0.01) noreplacement common logit

    *--------------------------------------------------------*
    * Assess covariate balance
    *--------------------------------------------------------*
    pstest age sex_n race1 i.race_new ///
        dxyear_simple_base i.dx_year_simple_num ///
        yost_quintile1 i.yost_acs_2013_2017_quintile ///
        income_group1 i.income_group ///
        rural urbanrural ///
        state1 i.state ///
        i.acute_mi_binary ///
        i.history_mi_binary ///
        i.chf_binary ///
        i.pvd_binary ///
        i.cvd_binary ///
        i.copd_binary ///
        i.dementia_binary ///
        i.paralysis_binary ///
        i.diabetes_binary ///
        i.diabetes_comp_binary ///
        i.renal_disease_binary ///
        i.mild_liver_disease_binary ///
        i.liver_disease_binary ///
        i.ulcers_binary ///
        i.rheum_disease_binary ///
        i.aids_binary ///
        tumor_size_clean ///
        t_stage_8_0 i.t_stage_8_num ///
        tnm_n_clean_0 i.tnm_n_clean ///
        brain_mri_pre12post3 ///
        pet_pre12post3 ///
        resection_post12 ///
        chemo_post12 ///
        radiation_post12_any ///
		immuno_post12, ///
        treated(screened) both

    *--------------------------------------------------------*
    * Identify matched pairs
    *--------------------------------------------------------*
    gen pair = _id if screened == 0
    replace pair = _n1 if screened == 1

    bysort pair: egen paircount = count(pair)
    drop if paircount != 2

    *--------------------------------------------------------*
    * Matched survival analysis
    *--------------------------------------------------------*
    stset fup_months, fail(vitalstatus == 1)

    sts graph, by(screened) ///
        xlabel(0(12)36) ///
        graphregion(color(white)) ///
        ylabel(,nogrid) ///
        risktable(,order(1 "Not Screened" 2 "Screened")) ///
        tmax(36) ///
        ytitle("Overall Survival (Probability)") ///
        xtitle("Time (months)") ///
        legend(order(1 "Not Screened" 2 "Screened")) ///
        title("`stagelabel' – Propensity Score Matched")

    sts list, by(screened) at(12 24 36)

    sts test screened, strata(pair)

    *--------------------------------------------------------*
    * Baseline characteristics table
    *--------------------------------------------------------*

dtable age tumor_size_clean i.sex_n i.race_new i.dx_year_simple_num i.acute_mi_binary i.history_mi_binary i.chf_binary i.pvd_binary i.cvd_binary i.copd_binary i.dementia_binary i.paralysis_binary i.diabetes_binary i.diabetes_comp_binary i.renal_disease_binary i.mild_liver_disease_binary i.liver_disease_binary i.ulcers_binary i.rheum_disease_binary i.aids_binary i.income_group i.yost_acs_2013_2017_quintile i.urbanrural i.state i.t_stage_8_num i.tnm_n_clean i.pet_pre12post3 i.brain_mri_pre12post3 i.chemo_post12 i.radiation_post12_any i.resection_post12 i.immuno_post12 if sex_n!=5, by(screened, nototal) continuous(age tumor_size_clean, statistics(median iqr)) factor(, statistics(fvfrequency fvpercent)) export(tables1_baselinecharacteristics_psmatch2_`stagegroup'.xlsx, replace)

}