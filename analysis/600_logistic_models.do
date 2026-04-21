version 16

/*==============================================================================
DO FILE NAME:			Logistic regression analyses
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
log using "$logdir/logistic_models.log", replace

*Set Ado file path
adopath + "$projectdir/analysis/extra_ados"

*Set disease list (passed from yaml)
global arglist disease
args $arglist

if $running_locally ==0 {
	foreach var of global arglist {
		local `var' : subinstr local `var' "|" " ", all
		global `var' "``var''"
		di "$`var'"
	}
}

if $running_locally ==1 {
	global disease "gout"
}
di "$disease"

set type double

set scheme plotplainblind

**Define programme to run multi-level logistic models and output values of interest===================

program define melogit_model, rclass

	***Model arguments
	args outcome model_terms focal_predictor model_label outlabel inclusion measures do_icc

	***Run model
	capture noisily melogit `outcome' `model_terms' || practice_id:, or
	
	****Skip if estimation failed
	if _rc {
		di as txt "Skipping model (estimation failure): `model_terms'"
		return scalar model_ok = 0
		exit
	}

	****Skip if model doesn't converge
	capture confirm scalar e(converged)
	if !_rc {
		if e(converged) == 0 {
			di as txt "Skipping model (no convergence): `model_terms'"
			return scalar model_ok = 0
			exit
		}
	}
	
	***Check to ensure model ran ok
	return scalar model_ok = 1
		
	***Store column names for terms
	local cnames : colnames e(b)
	
	***Cycle through terms
	foreach term of local cnames {
		
		****Skip intercept and random-effect variance parameters
		if "`term'" == "_cons" continue
		if strpos("`term'", "var(") continue
		
		****Skip omitted terms
		if regexm("`term'", "^(o\.|[0-9]+o\.).+$") {
			
			if regexm("`term'", "^o\.(.+)$") {
				local varname "`=regexs(1)'"
				local category "Omitted"
			}
			else if regexm("`term'", "^([0-9]+)o\.(.+)$") {
				local levelnum "`=regexs(1)'"
				local varname  "`=regexs(2)'"
				
				local labname : value label `varname'
				if "`labname'" != "" {
					capture local category : label `labname' `levelnum'
					if _rc local category "`levelnum' (omitted)"
					else local category "`category' (omitted)"
				}
				else {
					local category "`levelnum' (omitted)"
				}
			}
			
			*****Store variable label
			local varlabel : variable label `varname'
			if "`varlabel'" == "" local varlabel "`varname'"
			
			*****Post output for omitted variables
			post `measures' ("`inclusion'") ("`outcome'") ("`outlabel'") ("`model_label'") ///
				("`varlabel'") ("`category'") (.) (.) (.) (.)
			continue
		}
				
		****Skip adjustment terms in age/sex-adjusted models unless they are the focal predictor
        if "`model_label'" == "Age/sex-adjusted" {
            if "`focal_predictor'" != "age_decile" & "`term'" == "age_decile" continue
            if "`focal_predictor'" != "i.sex" & regexm("`term'", "^[0-9]+b?\.sex$") continue
        }

		****Defaults for continuous variables
		local varname "`term'"
		local category "Continuous"

		****Handle base levels
		if regexm("`term'", "^([0-9]+)b\.(.+)$") {
			local levelnum "`=regexs(1)'"
			local varname  "`=regexs(2)'"
			
			*****Store factor level labels
			local labname : value label `varname'
			if "`labname'" != "" {
				capture local category : label `labname' `levelnum'
				if _rc local category "`levelnum'"
			}
			else {
				local category "`levelnum'"
			}
			
			*****Skip adjustment base sex terms
            if "`model_label'" == "Age/sex-adjusted" {
                if "`focal_predictor'" != "i.sex" & "`varname'" == "sex" continue
            }
			
			*****Store variable label
			local varlabel : variable label `varname'
			if "`varlabel'" == "" local varlabel "`varname'"

			*****Post default output for base factor levels
			local oddsratio = 1
			local lower95 = .
			local upper95 = .
			local pvalue = .

			*****Post results for base factor levels
			post `measures' ("`inclusion'") ("`outcome'") ("`outlabel'") ("`model_label'") ("`varlabel'") ("`category'") (`oddsratio') (`lower95') (`upper95') (`pvalue')
			
			continue
		}

		****Output estimates for terms of interest
		capture scalar b = _b[`outcome':`term']
		if _rc continue

		capture scalar se = _se[`outcome':`term']
		if _rc continue
		if missing(se) continue

		****Calculate OR, CI, p-values
		local oddsratio = round(exp(_b[`outcome':`term']), 0.001)
		local lower95 = round(exp(_b[`outcome':`term'] - invnormal(0.975)*_se[`outcome':`term']), 0.001)
		local upper95 = round(exp(_b[`outcome':`term'] + invnormal(0.975)*_se[`outcome':`term']), 0.001)
		local pvalue = round(2*normal(-abs(_b[`outcome':`term']/_se[`outcome':`term'])), 0.0001)
/*
		
		local j = colnumb(T, "`term'")
		if `j' >= . continue

		local oddsratio = round(T[1,`j'], 0.001)
		local lower95   = round(T[5,`j'], 0.001)
		local upper95   = round(T[6,`j'], 0.001)
		local pvalue    = round(T[4,`j'], 0.0001)
*/
		****If factor-variable term, split level and variable
		if regexm("`term'", "^([0-9]+)\.(.+)$") {
			local levelnum "`=regexs(1)'"
			local varname  "`=regexs(2)'"

			*****Store factor level labels
			local labname : value label `varname'
			if "`labname'" != "" {
				capture local category : label `labname' `levelnum'
				if _rc local category "`levelnum'"
			}
			else {
				local category "`levelnum'"
			}
		}

		****Store variable label
		local varlabel : variable label `varname'
		if "`varlabel'" == "" local varlabel "`varname'"
			
		****Post model results
		post `measures' ("`inclusion'") ("`outcome'") ("`outlabel'") ("`model_label'") ("`varlabel'") ("`category'") (`oddsratio') (`lower95') (`upper95') (`pvalue')
	}
	
	****Output ICC (Proportion of the total variation in the outcome attributable to differences between practices)
    
	****Passed argument to run ICC
	if "`do_icc'" == "1" {
	
		*****Estimate ICC
		capture noisily estat icc

		*****Check estimation/model hasn't failed
		if !_rc {
			capture confirm scalar r(icc2)
			if !_rc {
				local icc = round(r(icc2),0.001)
				
				capture matrix CI = r(ci2)
				if !_rc {
					local icc_lo = round(el(CI,1,1), 0.001)
					local icc_hi = round(el(CI,1,2), 0.001)
				}
				else {
					local icc_lo = .
					local icc_hi = .
				}
				
				******Post ICC results
				post `measures' ("`inclusion'") ("`outcome'") ("`outlabel'") ("`model_label'") ///
					("Intra-class correlation") ("Practice-level ICC") (`icc') (`icc_lo') (`icc_hi') (.)
			}
		}
	}
end

**Define programme to run standard logistic models and output values of interest ===================

program define logistic_model

	args outcome model_terms focal_predictor model_label outlabel inclusion measures

	***Run logistic model
	capture noisily logistic `outcome' `model_terms', vce(cluster practice_id)
	
	****Skip if estimation failed
	if _rc {
		di as txt "Skipping model (estimation failure): `model_terms'"
		exit
	}
		
	****Store outputs from model
	matrix T = r(table)
	matrix list T
	local cnames : colnames T
	
	****Strip factor prefix from focal predictor
    local focalvar "`focal_predictor'"
    local focalvar = subinstr("`focalvar'", "i.", "", .)
    local focalvar = subinstr("`focalvar'", "c.", "", .)
	
	***Cycle through column names
	foreach term of local cnames {
		
		****Skip intercept
		if "`term'" == "_cons" continue
		
		****Defaults for continuous variables
		local varname "`term'"
		local category "Continuous"
		
		*****If factor-variable term, split level and variable
		if regexm("`term'", "^([0-9]+)([a-z]*)\.(.+)$") {
			local levelnum "`=regexs(1)'"
			local varname  "`=regexs(3)'"

			local labname : value label `varname'
			if "`labname'" != "" {
				local category : label `labname' `levelnum'
			}
			else {
				local category "`levelnum'"
			}
		}
		
		*****Select which variables to output, depending on model
        if "`model_label'" == "Multivariable" {
        }
        else {
            if "`varname'" != "`focalvar'" continue

			****Skip adjustment terms in age/sex-adjusted models unless they are the focal predictor
            if "`model_label'" == "Age/sex-adjusted" {
                if "`focal_predictor'" != "age_decile"   & "`term'" == "age_decile" continue
                if "`focal_predictor'" != "i.sex" & "`varname'" == "sex" continue
            }
        }

		*****Keep outputs of relevance
		local j = colnumb(T, "`term'")
		if missing(`j') continue

		local oddsratio = round(T[1,`j'], 0.001)
		local lower95 = round(T[5,`j'], 0.001)
		local upper95 = round(T[6,`j'], 0.001)
		local pvalue = round(T[4,`j'], 0.0001)

		****Store variable label
		local varlabel : variable label `varname'
		if "`varlabel'" == "" local varlabel "`varname'"
			
		****Post model results
		post `measures' ("`inclusion'") ("`outcome'") ("`outlabel'") ("`model_label'") ("`varlabel'") ("`category'") (`oddsratio') (`lower95') (`upper95') (`pvalue')
	}	
end

*Run multi-level logistic models ================

**Generate temporary file to store outputs
tempname melogit_measures
postfile `melogit_measures' str30(inclusion) str30(outcome) str30(outcome_label) str30(model) str30(variable) str30(category) double oddsratio lower95 upper95 pvalue ///
    using "$projectdir/output/data/melogit_summary.dta", replace

*Select parameters for models ==========================================

**Define inclusion criteria
local inclusions has_12m_fup has_12m_fup_ult

**Loop through inclusion criteria
foreach inclusion of local inclusions {

    di "Inclusion criterion: `inclusion'"

	**Import cleaned/processed cohort
	use "$projectdir/output/data/cohort_processed.dta", clear
	
    keep if `inclusion' == 1

    **Define outcomes (specific to inclusion criteria)
	if "`inclusion'" == "has_12m_fup" {
		*local outcomes ult_12m
		local outcomes ckd_comb_bl
    }
    else if "`inclusion'" == "has_12m_fup_ult" {
		*local outcomes urate_12m_ult
        local outcomes chd_bl
    }
    else {
        di as error "No outcomes defined for `inclusion'"
        continue
    }

    **Define patient-level predictors (specific to inclusion criteria, if required); also need to think about interactions between i.post_nice and variables - Nb. programmes about won't output factors properly
    local patient_predictors ///
        age_decile i.sex i.imd i.ethnicity i.diabetes_bl i.heart_failure_bl i.cva_bl i.bmicat i.diuretic_bl i.alcohol_bl i.smoke
		*i.ckd_comb_bl i.chd_bl //removed as this is temporary outcome
		*could add address urban vs. rural as well


    if "`inclusion'" == "has_12m_fup" {
        local patient_predictors `patient_predictors' i.post_nice_diag
    }
    else if "`inclusion'" == "has_12m_fup_ult" {
        local patient_predictors `patient_predictors' i.post_nice_ult
    }

	**Define practice-level predictors
    local practice_predictors practice_list_n practice_${disease}_ratio
		
	***Loop through outcomes
	foreach outcome of local outcomes {
		
		****Store outcome variable name
		local outlabel : variable label `outcome'
		if "`outlabel'" == "" local outlabel "`outcome'"
				
		****Estimate practice-level variation/ICC without predictors ============
		melogit_model `outcome' `""' `""' `"Empty model"' `"`outlabel'"' `"`inclusion'"' `melogit_measures' 1

		****Univariable models with patient-level and practice-level predictors ==============
		
		*****Store predictors
		local predictors `patient_predictors' `practice_predictors'
		
		*****Loop through predictors
		foreach predictor of local predictors {
			melogit_model `outcome' `"`predictor'"' `"`predictor'"' `"Univariable"' `"`outlabel'"' `"`inclusion'"' `melogit_measures' 0
		}
		
		****Age and sex-adjusted models (Nb. don't usually need to present age/sex-adjusted practice-level variables) ============
		
		*****Store predictors
		local predictors `patient_predictors' 
		
		*****Loop through predictors
		foreach predictor of local predictors {
			
			******Define age and sex-adjustment
			local model_terms `predictor' age_decile i.sex
			if "`predictor'" == "age_decile" local model_terms age_decile i.sex
			if "`predictor'" == "i.sex" local model_terms i.sex age_decile
			
			******Run model
			melogit_model `outcome' `"`model_terms'"' `"`predictor'"' `"Age/sex-adjusted"' `"`outlabel'"' `"`inclusion'"' `melogit_measures' 0
		}
		
		****Multivariable model with patient-level predictors only ===========		
		melogit_model `outcome' `"`predictors'"' `""' `"Multivariable patient"' `"`outlabel'"' `"`inclusion'"'  `melogit_measures' 1
		 
		****Multivariable model with patient-level predictors and practice-level predictors ============
		
		*****Store predictors
		local predictors `patient_predictors' `practice_predictors'
		
		*****Run model
		melogit_model `outcome' `"`predictors'"' `""' `"Multivariable pt/practice"' `"`outlabel'"' `"`inclusion'"' `melogit_measures' 1
		 
	}
}

postclose `melogit_measures'

*Output postfiles to csv
use "$projectdir/output/data/melogit_summary.dta", clear
format oddsratio lower95 upper95 %9.3f
format pvalue %9.4f

export delimited using "$projectdir/output/tables/melogit_summary.csv", replace

*Run logistic regression models =====================================

**Generate temporary file to store outputs
tempname logistic_measures
postfile `logistic_measures' str30(inclusion) str30(outcome) str30(outcome_label) str20(model) str30(variable) str30(category) double oddsratio lower95 upper95 pvalue ///
    using "$projectdir/output/data/logistic_summary.dta", replace
	
**Define inclusion criteria
local inclusions has_12m_fup has_12m_fup_ult

**Loop through inclusion criteria
foreach inclusion of local inclusions {

    di "Inclusion criterion: `inclusion'"

	**Import cleaned/processed cohort
	use "$projectdir/output/data/cohort_processed.dta", clear
	
    keep if `inclusion' == 1

    **Define outcomes (specific to inclusion criteria)
	if "`inclusion'" == "has_12m_fup" {
		*local outcomes ult_12m
		local outcomes ckd_comb_bl
    }
    else if "`inclusion'" == "has_12m_fup_ult" {
		*local outcomes urate_12m_ult
        local outcomes chd_bl
    }
    else {
        di as error "No outcomes defined for `inclusion'"
        continue
    }

    **Define patient-level predictors (specific to inclusion criteria, if required); also need to think about interactions between i.post_nice and variables - Nb. programmes about won't output factors properly
    local patient_predictors ///
        age_decile i.sex i.imd i.ethnicity i.diabetes_bl i.heart_failure_bl i.cva_bl i.bmicat i.diuretic_bl i.alcohol_bl i.smoke
		*i.ckd_comb_bl i.chd_bl //removed as this is temporary outcome

    if "`inclusion'" == "has_12m_fup" {
        local patient_predictors `patient_predictors' i.post_nice_diag
    }
    else if "`inclusion'" == "has_12m_fup_ult" {
        local patient_predictors `patient_predictors' i.post_nice_ult
    }

	**Loop through outcomes
	foreach outcome of local outcomes {
		
		***Store outcome variable name
		local outlabel : variable label `outcome'
		if "`outlabel'" == "" local outlabel "`outcome'"
		
		***Store predictors
		local predictors `patient_predictors'
		
		***Run univariable models ===============

		***Loop through predictors
		foreach predictor of local predictors {
			
			logistic_model `outcome' `"`predictor'"' `"`predictor'"' `"Univariable"' `"`outlabel'"' `"`inclusion'"' `logistic_measures'
		}
		
		***Run age and sex-adjusted models ================
		foreach predictor of local predictors {
			
			local model_terms `predictor' age_decile i.sex
			if "`predictor'" == "age_decile" local model_terms age_decile i.sex
			if "`predictor'" == "i.sex" local model_terms i.sex age_decile
			
			logistic_model `outcome' `"`model_terms'"' `"`predictor'"' `"Age/sex-adjusted"' `"`outlabel'"' `"`inclusion'"' `logistic_measures'
		}
		
		***Run multivariable model with all selected predictors ============
		logistic_model `outcome' `"`predictors'"' `""' `"Multivariable"' `"`outlabel'"' `"`inclusion'"' `logistic_measures'
		
	}
}

*Close tempfile
postclose `logistic_measures'

*Output postfiles to csv
use "$projectdir/output/data/logistic_summary.dta", clear
format oddsratio lower95 upper95 %9.3f
format pvalue %9.4f

export delimited using "$projectdir/output/tables/logistic_summary.csv", replace

log close
