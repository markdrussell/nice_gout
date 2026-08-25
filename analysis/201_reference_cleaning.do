version 16

/*==============================================================================
DO FILE NAME:			Clean primary cohort using dataset definition file
PROJECT:				OpenSAFELY NICE 
AUTHOR:					M Russell								
DATASETS USED:			Primary dataset defintion
USER-INSTALLED ADO: 	 
  (place .ado file(s) in analysis folder)						
==============================================================================*/

*Set filepaths
/*
global projectdir "C:\Users\k1754142\OneDrive\PhD Project\OpenSAFELY NICE\nice_gout"
global running_locally = 1 // Running on local machine
*/

global projectdir `c(pwd)'
global running_locally = 0 // Running on OpenSAFELY console

capture mkdir "$projectdir/output/data"
capture mkdir "$projectdir/output/figures"
capture mkdir "$projectdir/output/tables"

*Open log file
global logdir "$projectdir/logs"
cap log close
log using "$logdir/cohort_cleaning_ref.log", replace
log off

*Set Ado file path
adopath + "$projectdir/analysis/extra_ados"

*Run checks
global checks 0

*Set disease list, characteristics of interest, and study dates (passed from yaml)
global arglist studystart_date studyend_date studyfup_date comorbidities bloods
args $arglist

if $running_locally ==0 {
	foreach var of global arglist {
		local `var' : subinstr local `var' "|" " ", all
		global `var' "``var''"
		di "$`var'"
	}
}

if $running_locally ==1 {
	global studystart_date "2016-07-01"
	global studyend_date "2025-06-30"
	global studyfup_date "2026-06-30"
	global comorbidities "chd diabetes cva ckd hypertension depression heart_failure liver_disease transplant alcohol"
	global bloods "urate creatinine cholesterol hba1c"
}

set scheme plotplainblind

*Import primary dataset
import delimited "$projectdir/output/dataset_incidence_ref.csv", clear

*Conversion for dates ====================================================*

**Convert format for variables containing dates that are in string format
ds *date*, has(type string) //check list of variables is appropriate
local string_dates `r(varlist)'

foreach var of local string_dates {
	gen double `var'_num = daily(`var', "YMD")
	quietly count if !missing(`var'_num)
	if r(N) {
		format `var'_num %td
		order `var'_num, after(`var')
		drop `var'
		rename `var'_num `var'
	}
	else {
		drop `var'_num
	}
}

**Convert format for variables containing dates that are in numeric format
ds *date*, has(type numeric) //check list of variables is appropriate
capture ds *date*, has(type numeric)
if !_rc & "`r(varlist)'" != "" {
    foreach var of varlist `r(varlist)' {
        format `var' %td
    }
}

*Clean and label demographic and comorbidity variables ====================================*/

**Check cohort entry date
gen reference = 1 if cohort_entry_date!=.
codebook cohort_entry_date

**Age
lab var age "Age at diagnosis"
tabstat age, stat(n mean sd p50 p25 p75)

***Define 10-year age bands
recode age 18/29.9999 = 1 /// 
		   30/39.9999 = 2 ///
           40/49.9999 = 3 ///
		   50/59.9999 = 4 ///
	       60/69.9999 = 5 ///
		   70/79.9999 = 6 ///
		   80/max = 7, gen(agegroup_broad) 

label define agegroup_broad	1 "18 to 29" ///
							2 "30 to 39" ///
							3 "40 to 49" ///
							4 "50 to 59" ///
							5 "60 to 69" ///
							6 "70 to 79" ///
							7 "80 or above", modify
						
label values agegroup_broad agegroup_broad
lab var agegroup_broad "Age group"
order agegroup_broad, after(age)
tab agegroup_broad, missing

***Define broader age bands
recode age 18/39.9999 = 1 /// 
		   40/59.9999 = 2 ///
           60/79.9999 = 3 ///
		   80/max = 4, gen(agegroup) 

label define agegroup 	1 "18 to 39" ///
						2 "40 to 59" ///
						3 "60 to 79" ///
						4 "80 or above", modify
						
label values agegroup agegroup
lab var agegroup "Age group"
order agegroup, after(age)
tab agegroup, missing

**Sex
rename sex sex_s
encode sex_s, gen(sex)
lab def sex 1 "Female" 2 "Male", modify
lab val sex sex
lab var sex "Sex"
drop sex_s
tab sex, missing

**Ethnicity
gen ethnicity_n = 1 if ethnicity == "White"
replace ethnicity_n = 2 if ethnicity == "Asian or Asian British"
replace ethnicity_n = 3 if ethnicity == "Black or Black British"
replace ethnicity_n = 4 if ethnicity == "Mixed"
replace ethnicity_n = 5 if ethnicity == "Chinese or Other Ethnic Groups"
replace ethnicity_n = 9 if ethnicity == "Unknown"

label define ethnicity_n	1 "White"  	///
							2 "Asian"	///
							3 "Black"  	///
							4 "Mixed"	///
							5 "Chinese or Other" ///
							9 "Not known", modify
							
label val ethnicity_n ethnicity_n
lab var ethnicity_n "Ethnicity"
drop ethnicity
rename ethnicity_n ethnicity 
tab ethnicity, missing

**Practice region (at cohort entry)
replace region="Not known" if region==""
replace region="Yorkshire Humber" if region=="Yorkshire and The Humber"
encode region, gen(nuts_region)
replace nuts_region = 9 if region=="Not known"
drop region
rename nuts_region region
lab var region "Region"
tab region, missing

**Index of multiple deprivation (at cohort entry)
gen imd = 1 if imd_quintile == "1 (most deprived)"
replace imd = 2 if imd_quintile == "2"
replace imd = 3 if imd_quintile == "3"
replace imd = 4 if imd_quintile == "4"
replace imd = 5 if imd_quintile == "5 (least deprived)"
replace imd = 9 if imd_quintile == "Unknown"

label define imd 1 "1 most deprived" 2 "2" 3 "3" 4 "4" 5 "5 least deprived" 9 "Not known", modify
label val imd imd 
lab var imd "Index of multiple deprivation"
drop imd_quintile
tab imd, missing

**Body Mass Index
***Recode values that are more likely to be erroneous
replace bmi_value = . if !inrange(bmi_value, 10, 80)

***Restrict to BMI recorded within 10 years of primary diagnosis date and aged > 16 years old
gen bmi_time = (cohort_entry_date - bmi_date)/365.25
gen bmi_age = age - bmi_time
replace bmi_value = . if bmi_age < 16 
replace bmi_value = . if bmi_time > 10 & bmi_time !=. 
replace bmi_value = . if bmi_date == . 
replace bmi_date = . if bmi_value == . 
replace bmi_time = . if bmi_value == . 
replace bmi_age = . if bmi_value == . 

***Create BMI categories
gen bmicat = .
recode bmicat . = 1 if bmi_value < 18.5
recode bmicat . = 2 if bmi_value < 25
recode bmicat . = 3 if bmi_value < 30
recode bmicat . = 4 if bmi_value < 35
recode bmicat . = 5 if bmi_value < 40
recode bmicat . = 6 if bmi_value >= 40 & bmi_value!=.
replace bmicat = 9 if bmi_value == .

label define bmicat 1 "Underweight (<18.5)" 	///
					2 "Normal (18.5-24.9)"		///
					3 "Overweight (25-29.9)"	///
					4 "Obese I (30-34.9)"		///
					5 "Obese II (35-39.9)"		///
					6 "Obese III (40+)"			///
					9 "Not known"
					
label values bmicat bmicat
lab var bmicat "BMI"
order bmicat, after (bmi_value)
drop bmi_age bmi_time
tab bmicat, missing

**Smoking status
gen smoke = 1 if smoking_status == "N"
replace smoke = 2 if smoking_status == "E"
replace smoke = 3 if smoking_status == "S"
replace smoke = 9 if smoking_status == "M"
replace smoke = 9 if smoking_status == "" 
label define smoke 1 "Never" 2 "Former" 3 "Current" 9 "Not known"
label values smoke smoke
lab var smoke "Smoking status"
drop smoking_status ever_smoked most_recent_smoking_code
tab smoke, missing

***Create non-missing 3-category variable for current smoking (assumes missing smoking is never smoking)
recode smoke 9 = 1, gen(smoke_nomiss)
order smoke_nomiss, after(smoke)
label values smoke_nomiss smoke
lab var smoke_nomiss "Smoking status"
tab smoke_nomiss, missing

**Clinical comorbidities at baseline and after diagnosis (using recorded codes only for now)
foreach comorbidity in $comorbidities {
    local lbl : subinstr local comorbidity "_" " ", all
	local lbl = strproper("`lbl'")
	di "`lbl'"
	gen `comorbidity'_bl = 1 if (`comorbidity'_date <= cohort_entry_date) & `comorbidity'_date!=.
	recode `comorbidity'_bl .=0
	lab define `comorbidity'_bl 0 "No" 1 "Yes", modify
	lab val `comorbidity'_bl `comorbidity'_bl
	lab var `comorbidity'_bl "`lbl'"
	order `comorbidity'_bl, after(`comorbidity'_date)
	tab `comorbidity'_bl, missing
}

***Re-label variables with acronyms (amend manually)
lab var ckd_bl "CKD"
lab var diabetes_bl "T2DM"
lab var chd_bl "CHD"
lab var cva_bl "Stroke/TIA"
lab var liver_disease_bl "Chronic liver disease"
lab var transplant_bl "Solid organ transplant"
lab var alcohol_bl "Excess alcohol"

**Specify drug of interest: diuretics 
foreach drug in diuretic {

	local Drug = strproper("`drug'") //first letter capitalised for labelling
	di "`drug'"

	***Drug use at baseline (prescribed within 6 months prior to index diagnosis)
	gen `drug'_bl = 1 if (cohort_entry_date - `drug'_bl_date) <= 183 & `drug'_bl_date != .
	recode `drug'_bl .=0
	lab define `drug'_bl 0 "No" 1 "Yes", modify
	lab val `drug'_bl `drug'_bl
	lab var `drug'_bl "`Drug' at baseline"
}

**Relevant blood tests (passed from yaml, but amend thresholds as necessary)

foreach blood in $bloods {
	
	**Set thresholds for plausible low and high values (amend as necessary)
	if "`blood'" == "creatinine" {
		local low 20
		local high 3000
	}
	else if "`blood'" == "hba1c" {
		local low 10
		local high 200
	}
	else if "`blood'" == "urate" {
		local low 0.05
		local high 3000
	}
	else if "`blood'" == "cholesterol" {
		local low 0.5
		local high 20
	}
	else {
		local low 0
		local high 100000
	}

	***Set implausible values to missing (as defined above)
	tabstat `blood'_bl_value, stats(n mean sd p50 p25 p75)
	summ `blood'_bl_value if inrange(`blood'_bl_value, `low', `high')
	replace `blood'_bl_value = . if !inrange(`blood'_bl_value, `low', `high')
	replace `blood'_bl_value = . if `blood'_bl_date == . 
	replace `blood'_bl_date = . if `blood'_bl_value == . 
	tabstat `blood'_bl_date, stats(n mean sd p50 p25 p75)

	gen had_`blood'_bl = 1 if `blood'_bl_value!=.
	recode had_`blood'_bl .=0
	lab var had_`blood'_bl "Serum `blood' performed at baseline"
	lab define had_`blood'_bl 0 "No" 1 "Yes", modify
	lab val had_`blood'_bl had_`blood'_bl
	tab had_`blood'_bl
}

**Generate eGFR from serum creatinine (using CKD-EPI formula with no ethnicity) ============================*/ 
gen SCr_adj = creatinine_bl_value/88.4

gen min = .
replace min = SCr_adj/0.7 if sex==1
replace min = SCr_adj/0.9 if sex==2
replace min = min^-0.329  if sex==1
replace min = min^-0.411  if sex==2
replace min = 1 if min<1

gen max = .
replace max=SCr_adj/0.7 if sex==1
replace max=SCr_adj/0.9 if sex==2
replace max=max^-1.209
replace max=1 if max>1

gen egfr_bl_value=min*max*141
replace egfr_bl_value=egfr_bl_value*(0.993^age)
replace egfr_bl_value=egfr_bl_value*1.018 if sex==1
label var egfr_bl_value "eGFR"
drop min max SCr_adj

gen egfr_bl_date = creatinine_bl_date
format egfr_bl_date %td
lab var egfr_bl_value "Baseline eGFR at diagnosis"
lab var egfr_bl_date "Date of eGFR at diagnosis"

***Categorise baseline eGFR into CKD stages
gen egfr_bl_cat = .
recode egfr_bl_cat . = 3 if egfr_bl_value < 30
recode egfr_bl_cat . = 2 if egfr_bl_value < 60
recode egfr_bl_cat . = 1 if egfr_bl_value < .
replace egfr_bl_cat = 9 if egfr_bl_value >= .

label define egfr_bl_cat 	1 ">=60" 		///
							2 "30-59"		///
							3 "<30"			///
							9 "Not known"
					
label val egfr_bl_cat egfr_bl_cat
lab var egfr_bl_cat "eGFR at baseline"
tab egfr_bl_cat, missing

***Categorise baseline eGFR into more granular CKD stages
gen egfr_bl_finecat = .
recode egfr_bl_finecat . = 6 if egfr_bl_value < 15
recode egfr_bl_finecat . = 5 if egfr_bl_value < 30
recode egfr_bl_finecat . = 4 if egfr_bl_value < 45
recode egfr_bl_finecat . = 3 if egfr_bl_value < 60
recode egfr_bl_finecat . = 2 if egfr_bl_value < 90
recode egfr_bl_finecat . = 1 if egfr_bl_value < .
replace egfr_bl_finecat = 9 if egfr_bl_value >= .

label define egfr_bl_finecat 	1 ">=90" 		///
								2 "60-89"		///
								3 "45-59"		///
								4 "30-44"		///
								5 "15-29"		///
								6 "<15"			///
								9 "Not known"
					
label val egfr_bl_finecat egfr_bl_finecat
lab var egfr_bl_finecat "eGFR at baseline"
tab egfr_bl_finecat, missing

***Generate baseline CKD code that combines CKD coding (stages 3-5) + eGFR (stages 3-5)
gen ckd_comb_bl = 0
replace ckd_comb_bl = 1 if egfr_bl_value != . & egfr_bl_value < 60
replace ckd_comb_bl = 1 if ckd_bl == 1
label define ckd_comb_bl 0 "No" 1 "Yes"
label val ckd_comb_bl ckd_comb_bl
label var ckd_comb_bl "CKD"
tab ckd_comb_bl, missing

**Categorise HbA1c at baseline ============================*/
gen hba1c_bl_cat = 0 if hba1c_bl_value < 58
replace hba1c_bl_cat = 1 if hba1c_bl_value >= 58 & hba1c_bl_val !=.
replace hba1c_bl_cat = 9 if hba1c_bl_cat ==. 
label define hba1c_bl_cat 0 "HbA1c <58mmol/mol" 1 "HbA1c >=58mmol/mol" 9 "Not known"
label val hba1c_bl_cat hba1c_bl_cat
lab var hba1c_bl_cat "HbA1c at baseline"
tab hba1c_bl_cat, missing

***Create combined diabetes code that combines diabetes coding + HBA1c
gen diab_bl_cat = 1 if diabetes_bl==0
replace diab_bl_cat = 2 if diabetes_bl==1 & hba1c_bl_cat==0
replace diab_bl_cat = 3 if diabetes_bl==1 & hba1c_bl_cat==1
replace diab_bl_cat = 4 if diabetes_bl==1 & hba1c_bl_cat==9

label define diab_bl_cat 1 "No diabetes" 			///
						2 "Diabetes with HbA1c <58mmol/mol"		///
						3 "Diabetes with HbA1c >58mmol/mol" 	///
						4 "Diabetes with no recorded HbA1c"
label values diab_bl_cat diab_bl_cat
lab var diab_bl_cat "Type 2 diabetes mellitus with HbA1c categorisation"
tab diab_bl_cat, missing

save "$projectdir/output/data/cohort_processed_ref.dta", replace

log close
