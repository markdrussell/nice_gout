version 16

/*==============================================================================
DO FILE NAME:			Produce rounded and redacted summary tables using cleaned referencec cohort
PROJECT:				OpenSAFELY NICE 
AUTHOR:					M Russell								
DATASETS USED:			Cleaned primary cohort
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
log using "$logdir/summary_tables_ref.log", replace

*Set Ado file path
adopath + "$projectdir/analysis/extra_ados"

*Set disease list, characteristics of interest, and study dates (passed from yaml)
global arglist comorbidities bloods
args $arglist

if $running_locally ==0 {
	foreach var of global arglist {
		local `var' : subinstr local `var' "|" " ", all
		global `var' "``var''"
		di "$`var'"
	}
}

if $running_locally ==1 {
	global comorbidities "chd diabetes cva ckd hypertension depression heart_failure liver_disease transplant alcohol"
	global bloods "urate creatinine cholesterol hba1c"
}

set type double

set scheme plotplainblind

*Function to round and redact categorical variables ======================*/
program define rounded_categorical
    syntax varlist(min=1 max=1), outfile(string)
	local outcome `varlist'
	preserve 
		contract `outcome'
		local outcome_desc : variable label `outcome' 
		gen variable = `"`outcome_desc'"'
		decode `outcome', gen(categories)
		replace categories = "Missing" if categories == ""
		gen count = round(_freq, 5)
		egen total = total(count)
		gen percent = round((count/total)*100, 0.01)
		order total, before(percent)
		*egen mincount = min(count)
		replace percent =. if count<=7
		replace total =. if count<=7
		replace count =. if count<=7
		format percent %14.4f
		format count total %14.0f
		list variable categories count total percent
		keep variable categories count total percent
		capture append using `"`outfile'"'
		save `"`outfile'"', replace
    restore
end

*Function to round and redact continuous variables (amend to sum, rather than count, for ordinal variables) ======================*/
program define rounded_continuous
    syntax varlist(min=1 max=1), outfile(string)
	local outcome `varlist'
	preserve 
		local outcome_desc : variable label `outcome'
		collapse (count) count_un=`outcome' total_un=reference (mean) mean=`outcome' (sd) stdev=`outcome'
		gen count = round(count_un, 5)
		gen total = round(total_un, 5)
		replace stdev = . if count<=7
		replace mean = . if count<=7
		replace total = . if count<=7
		replace count = . if count<=7
		gen variable = `"`outcome_desc'"'
		order variable, first
		gen categories = "Not applicable"
		order categories, after(variable)
		order count, after(stdev)
		order total, after(count)
		format mean %14.2f
		format stdev %14.2f
		format count %14.0f
		list variable categories mean stdev count total
		keep variable categories mean stdev count total
		capture append using `"`outfile'"'
		save `"`outfile'"', replace
    restore
end

*Baseline summary table ========================*

**Store table name
local cohort "baseline"

**Erase any existing data file
capture erase "$projectdir/output/data/summary_table_ref_`cohort'.dta"

**Load processed dataset
use "$projectdir/output/data/cohort_processed_ref.dta", clear

**Store list of additional variables of interest
local comorbidity_vars_bl
foreach comorbidity in $comorbidities {
    local comorbidity_vars_bl `comorbidity_vars_bl' `comorbidity'_bl
}
di "`comorbidity_vars_bl'"

local blood_vars_value_bl
foreach blood in $bloods {
    local blood_vars_value_bl `blood_vars_value_bl' `blood'_bl_value
}
di "`blood_vars_value_bl'"

local blood_vars_test_bl
foreach blood in $bloods {
    local blood_vars_test_bl `blood_vars_test_bl' had_`blood'_bl
}
di "`blood_vars_test_bl'"

**Process catergorical outcomes of interest (other than outpatient-based variables)
foreach outcome of varlist diuretic_bl diab_bl_cat hba1c_bl_cat ckd_comb_bl egfr_bl_cat `blood_vars_test_bl' `comorbidity_vars_bl' bmicat smoke region imd ethnicity sex agegroup {
	rounded_categorical `outcome', outfile("$projectdir/output/data/summary_table_ref_`cohort'.dta")
}

**Process continuous outcomes of interest
foreach outcome of varlist `blood_vars_value_bl' age {
	rounded_continuous `outcome', outfile("$projectdir/output/data/summary_table_ref_`cohort'.dta")
}

**Export to CSV
use "$projectdir/output/data/summary_table_ref_`cohort'.dta", clear
export delimited using "$projectdir/output/tables/summary_table_ref_`cohort'.csv", datafmt replace

log close
