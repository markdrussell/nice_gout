version 16

/*==============================================================================
DO FILE NAME:			Produce temporal plots
PROJECT:				OpenSAFELY NICE 
AUTHOR:					M Russell								
DATASETS USED:			Rounded and redacted data tables
USER-INSTALLED ADO: 	 
  (place .ado file(s) in analysis folder)						
==============================================================================*/

*Set filepaths
/*
global projectdir "C:\Users\k1754142\OneDrive\PhD Project\OpenSAFELY NICE\nice_gout"
global projectdir "C:\Users\Mark\OneDrive\PhD Project\OpenSAFELY NICE\nice_gout"
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
log using "$logdir/temporal_plots.log", replace

*Set Ado file path
adopath + "$projectdir/analysis/extra_ados"

*Set disease list, characteristics of interest, and study dates (passed from yaml)
global arglist disease demographic studystart_date studyend_date studyfup_date intervention_date_2
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
	global demographic "agegroup sex ethnicity imd region"
	global studystart_date "2016-07-01"
	global studyend_date "2025-06-30"
	global studyfup_date "2025-12-31"
	global intervention_date_2 "2022-06-01"
}

di "$disease"
di "$studystart_date"
di "$studyend_date"
di "$studyfup_date"
di "$intervention_date_2"

global intervention_date "$intervention_date_2"
di "$intervention_date"

*Start year, end year, intervention date (derived from above)
foreach date in studystart studyend studyfup intervention {
	local year = real(substr("$`date'_date", 1, 4))
	local month = real(substr("$`date'_date", 6, 2))
	local year_month = ym(`year', `month')
	global `date' `year_month'
	global `date'_year `year'
	global `date'_tm `year'm`month'
	di %tm $`date'
	di %ty $`date'_year
}

set type double

set scheme plotplainblind

*Single line figures for full cohort ==================================*/

**Loop through data tables with different inclusion criteria

**Baseline data (no additional inclusion criteria)
**Individuals with X months minimum duration of follow-up post-diagnosis
**Individuals with X months minimum duration of follow-up post-ULT initiation, assuming ULT initiation within X months of diagnosis
**Individuals who should have been offered ULT at diagnosis or subsequently on the basis of risk factors

foreach table in ultrisk postult postdiagnosis baseline {
	
	**Import rounded and redacted data tables
	import delimited "$projectdir/output/tables/data_table_`table'.csv", clear
	di "`table'"

	**Extract outcomes of interest from data table
	levelsof outcome_name, local(outcomes)
	di `outcomes'

	**Loop through outcomes of interest for full cohort
	foreach outcome in `outcomes' {

		preserve
		
		**Keep outcome of interest
		keep if outcome_name == "`outcome'"
		di "`outcome'"
		
		**Reshape to long format
		reshape long count_ total_ prop_, i(month_year outcome_name outcome_desc) j(demographic) string
		gen demog_group = substr(demographic, 1, 3)
		gen demog_level = substr(demographic, 5, .)
		replace demog_level = subinstr(demog_level, "_", " ", .)
		rename count_ count
		rename total_ total
		rename prop_ prop
		drop demographic
		
		**Keep full cohort only
		keep if demog_group == "all"
		
		**Convert date format
		rename month_year month_year_s
		gen month_year = monthly(month_year_s, "MY") 
		format month_year %tmMon-CCYY
		drop month_year_s
		order month_year, after(outcome_desc)
		
		**Generate 3-monthly moving averages for proportions
		sort month_year
		gen prop_ma = (prop[_n-1]+prop[_n]+prop[_n+1])/3

		**Set y-axis format and title
		quietly summarize prop, meanonly
		local smin = r(min)
		local format "format(%03.2f)"
		local ytitle "Proportion with outcome"

		***Set x-axis format and title
		local xlabel ""
		forvalues y = $studystart_year(1)`= $studyend_year + 1' {
			local m = ym(`y', 1)
			local xlabel `xlabel' `m' "`y'" 
		}
		di as txt `"`xlabel'"'
		local xtitle "Date of diagnosis"
		
		***Set title
		local outcome_desc = outcome_desc

		***Temporal plot over study period (scatter with moving average)
		twoway scatter prop month_year, ytitle("`ytitle'", size(medsmall)) color(emerald%20) msymbol(circle) || line prop_ma month_year, lcolor(emerald) lstyle(solid) ylabel(, `format' nogrid labsize(small)) xtitle("`xtitle'", size(medsmall) margin(medsmall)) xlabel(`xlabel', nogrid labsize(small)) title("`outcome_desc'", size(medium) margin(b=2)) xline($intervention) legend(off) name("`outcome'", replace) saving("$projectdir/output/figures/temporal_plot_`outcome'.gph", replace)
		*graph export "$projectdir/output/figures/temporal_plot_`table'_`outcome'.png", replace
		graph export "$projectdir/output/figures/temporal_plot_`table'_`outcome'.svg", replace
/*
		***Set time-series for ITSA graphs
		tsset month_year

		***Extract minimum follow-up duration from variable name
		if regexm("`outcome'", "([0-9]+)m") {
			local fup_months = regexs(1)
		}
		else {
			local fup_months = 0
		}
		di "`fup_months'"

		***Skip if there are any gaps in monthly series or missing data
		local start = $studystart
		local end = ($studyfup - `fup_months')	
		quietly count if inrange(month_year, `start', `end') & missing(prop)
		local n_missing = r(N)
		di "`n_missing'"
		if (`n_missing' > 0) {
			di "Skipping ITSA for `outcome': missing data in study window"
		}
		else {
			***ITSA (Newey Standard Errors with 5 lags)
			
			di as txt `"`xlabel'"'
			*macro list xlabel

			itsa prop if inrange(month_year, `start', `end'), single trperiod($intervention_tm) lag(5) replace figure(title("`outcome_desc'", size(small)) subtitle("", size(medsmall)) ytitle("`ytitle'", size(medsmall) margin(small)) xlabel(`xlabel', nogrid labsize(small)) ylabel(, nogrid `format' labsize(small)) xtitle("`xtitle'", size(medsmall) margin(medsmall)) note("", size(v.small)) legend(off)) posttrend 
			*graph export "$projectdir/output/figures/ITSA_`table'_`outcome'.png", replace
			graph export "$projectdir/output/figures/ITSA_`table'_`outcome'.svg", replace
			
			actest, lag(18)		
		}
		*/
		restore
	}
}

*Multi-line figures by demographic characteristics ==================================*/

**Loop through data tables
foreach table in ultrisk postult postdiagnosis baseline {
	
	***Import rounded and redacted data tables
	import delimited "$projectdir/output/tables/data_table_`table'.csv", clear
	di "`table'"
	
	tempfile table_base
	save `table_base', replace
		
	**Extract outcomes of interest from data table
	levelsof outcome_name, local(outcomes)
	di `outcomes'

	**Loop through outcomes of interest
	foreach outcome in `outcomes' {
		
		use `table_base', clear
	
		keep if outcome_name == "`outcome'"
		di "`outcome'"
		
		**Reshape to long format
		reshape long count_ total_ prop_, i(month_year outcome_name outcome_desc) j(demographic) string
		gen demog_group = substr(demographic, 1, 3)
		gen demog_level = substr(demographic, 5, .)
		replace demog_level = subinstr(demog_level, "_", " ", .)
		replace demog_level = proper(demog_level)
		rename count_ count
		rename total_ total
		rename prop_ prop
		drop demographic
		
		**Remove not known categories
		drop if regexm(demog_level, "Not Known")
		
		**Save temporary file for that outcome
		tempfile outcome_base
		save `outcome_base', replace
		
		**Loop through demographic variables of interest
		foreach demog_var in $demographic {
			
			use `outcome_base', clear
			
			di "`demog_var'"
			keep if demog_group == substr("`demog_var'", 1, 3)
			
			***Skip if no observations
			count
			if r(N)==0 continue
	
			***Convert date format
			rename month_year month_year_s
			gen month_year = monthly(month_year_s, "MY") 
			format month_year %tmMon-CCYY
			drop month_year_s
			order month_year, after(outcome_desc)
				
			***Generate 3-monthly moving averages for proportions
			bys demog_level (month_year): gen prop_ma = (prop[_n-1]+prop[_n]+prop[_n+1])/3

			***Set y-axis format and title
			local format "format(%03.2f)"
			local ytitle "Proportion with outcome"

			***Set x-axis format and title
			local xlabel ""
			forvalues y = $studystart_year(1)`= $studyend_year + 1' {
				local m = ym(`y', 1)
				local xlabel `xlabel' `m' "`y'"
			}
			di as txt `"`xlabel'"'
			local xtitle "Date of diagnosis"
			
			***Set title
			local outcome_desc = outcome_desc
			
			***Choose colour palette based on demographic variable (change as needed)
			if "`demog_var'" == "sex" {
				local colours "red blue"
			}
			else if "`demog_var'" == "agegroup" {
				local colours "ltblue eltblue midblue ebblue blue navy black"
			}
			else if inlist("`demog_var'", "imd", "ethnicity") {
				local colours "ltblue eltblue ebblue blue navy"
			}
			else {
				local colours "emerald orange red blue dkgreen cranberry navy maroon teal sienna purple"
			}

			***Store plots and legend labels
			local plots ""
			local legorder ""
			local leglabels ""
			
			***Extract variables of interest from data table
			levelsof demog_level, local(demog_subset)

			local i = 0
			foreach subset of local demog_subset {
				di as txt `"`subset'"'
				local ++i
				local colour : word `i' of `colours'
				if "`colour'"=="" local colour "black" // fallback if more outcomes than colours

				****Two plots per subset: scatter and moving average line
				local plots `plots' ///
					(scatter prop month_year if demog_level=="`subset'", mcolor(`colour'%20) msymbol(circle)) ///
					(line prop_ma month_year if demog_level=="`subset'", lcolor(`colour') lpattern(solid))
				di as txt `"`plots'"'

				****Legend: keep only the moving average lines (2,4,6,...)
				local lineidx = 2*`i'
				local legorder `legorder' `lineidx'
				local outcome_disp : subinstr local subset "_" " " , all   //
				local leglabels `leglabels' label(`lineidx' "`outcome_disp'")
				di as txt `"`leglabels'"'
			}
			
			***Shorten name if long variable
			local gname = strtoname("`table'_`outcome'_`demog_var'")
			local gname = substr("`gname'", 1, 32)

			***Build plots
			twoway `plots', ytitle("`ytitle'", size(medsmall)) ylabel(, `format' nogrid labsize(small)) xtitle("`xtitle'", size(medsmall) margin(medsmall)) xlabel(`xlabel', nogrid labsize(small)) title("`outcome_desc'", size(medium) margin(b=2)) xline($intervention) legend(order(`legorder') `leglabels') name(`gname', replace) saving("$projectdir/output/figures/temporal_plot_`table'_`outcome'_`demog_var'.gph", replace)
			*graph export "$projectdir/output/figures/temporal_plot_`table'_`outcome'_`demog_var'.png", replace
			graph export "$projectdir/output/figures/temporal_plot_`table'_`outcome'_`demog_var'.svg", replace
		}	
	}	
}

*Multi-line figures by treatment ==================================*/
foreach table in flares {
	
	***Import rounded and redacted data tables
	import delimited "$projectdir/output/tables/data_table_`table'.csv", clear
	di "`table'"

	***Convert date format
	rename month_year month_year_s
	gen month_year = monthly(month_year_s, "MY") 
	format month_year %tmMon-CCYY
	drop month_year_s
	order month_year, after(outcome_desc)
		
	***Generate proportions from rounded data
	gen prop = count/total
	
	***Generate 3-monthly moving averages for proportions
	bys outcome_name (month_year): gen prop_ma = (prop[_n-1]+prop[_n]+prop[_n+1])/3

	***Set y-axis format and title
	local format "format(%03.2f)"
	local ytitle "Proportion with treatment" //change if needed

	***Set x-axis format and title
	local xlabel ""
	forvalues y = $studystart_year(1)`= $studyend_year + 1' {
		local m = ym(`y', 1)
		local xlabel `xlabel' `m' "`y'"
	}
	local xtitle "Date of diagnosis"
	
	***Extract variables of interest from data table
	levelsof outcome_name, local(outcomes)
	di as txt `"`outcomes'"'
	
	***Colour palette to cycle through
    local colours "emerald orange blue dkgreen cranberry navy maroon teal sienna purple"

    ***Store plots and legend labels
    local plots ""
    local legorder ""
    local leglabels ""

    local i = 0
    foreach outcome of local outcomes {
        local ++i
        local colour : word `i' of `colours'
        if "`colour'"=="" local colour "black" // fallback if more outcomes than colours

        ****Two plots per outcome: scatter and moving average line
        local plots `plots' ///
            (scatter prop month_year if outcome_name=="`outcome'", mcolor(`colour'%20) msymbol(circle)) ///
            (line prop_ma month_year if outcome_name=="`outcome'", lcolor(`colour') lpattern(solid))
		di as txt `"`plots'"'

        ****Legend: keep only the moving average lines (2,4,6,...)
        local lineidx = 2*`i'
        local legorder `legorder' `lineidx'
		local outcome_disp : subinstr local outcome "_" " " , all   //
        local leglabels `leglabels' label(`lineidx' "`outcome_disp'")
		di as txt `"`leglabels'"'
    }

	***Build plots
    twoway `plots', ytitle("`ytitle'", size(medsmall)) ylabel(, `format' nogrid labsize(small)) xtitle("`xtitle'", size(medsmall) margin(medsmall)) xlabel(`xlabel', nogrid labsize(small)) title("", size(medium) margin(b=2)) xline($intervention) legend(order(`legorder') `leglabels') name("`table'", replace) saving("$projectdir/output/figures/temporal_plot_`table'.gph", replace)
	*graph export "$projectdir/output/figures/temporal_plot_`table'.png", replace
	graph export "$projectdir/output/figures/temporal_plot_`table'.svg", replace
}

log close
