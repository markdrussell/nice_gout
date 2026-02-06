version 16

/*==============================================================================
DO FILE NAME:			Box plots
PROJECT:				Inflammatory Rheum OpenSAFELY project
DATE: 					20/05/2025
AUTHOR:					M Russell
DESCRIPTION OF FILE:	Box plots
DATASETS USED:			main data file
DATASETS CREATED: 		Box plots and outputs
OTHER OUTPUT: 			logfiles, printed to folder $Logdir
USER-INSTALLED ADO: 	 
  (place .ado file(s) in analysis folder)						
==============================================================================*/

**Set filepaths
*global projectdir "C:\Users\k1754142\OneDrive\PhD Project\OpenSAFELY NEIAA\inflammatory_rheum"
*global projectdir "C:\Users\Mark\OneDrive\PhD Project\OpenSAFELY NEIAA\inflammatory_rheum"
global projectdir `c(pwd)'

capture mkdir "$projectdir/output/data"
capture mkdir "$projectdir/output/figures"
capture mkdir "$projectdir/output/tables"

global logdir "$projectdir/logs"

**Open a log file
cap log close
log using "$logdir/box_plots.log", replace

**Set Ado file path
adopath + "$projectdir/analysis/extra_ados"

**Use cleaned data from previous step
use "$projectdir/output/data/file_eia_all.dta", clear

set scheme plotplainblind

**Set index dates ===========================================================*/

global index_date = "01/04/2016"
global start_date = "01/04/2019" //for outpatient analyses, data only be available from April 2019
global end_date = "31/03/2025"
global fup_date = "30/09/2025"
global base_year = year(date("$start_date", "DMY"))
global end_year = year(date("$end_date", "DMY"))
global max_year = $end_year - $base_year

*Rheum ref to appt performance by region, all years, using HES referral received date below==========================================================================*/

gen qs2_0 =1 if time_hes_rheum_appt<=21 & time_hes_rheum_appt!=.
recode qs2_0 .=0 if time_hes_rheum_appt!=.
gen qs2_1 =1 if time_hes_rheum_appt>21 & time_hes_rheum_appt<=42 & time_hes_rheum_appt!=.
recode qs2_1 .=0 if time_hes_rheum_appt!=.
gen qs2_2 = 1 if time_hes_rheum_appt>42 & time_hes_rheum_appt!=.
recode qs2_2 .=0 if time_hes_rheum_appt!=.

expand=2, gen(copy)
replace nuts_region = 0 if copy==1  
lab define nuts_region 0 "National" 9 "Yorkshire/Humber", modify
lab val nuts_region nuts_region

graph hbar (mean) qs2_0 (mean) qs2_1 (mean) qs2_2, over(nuts_region, relabel(1 "National")) stack ytitle(Percentage of patients) ytitle(, size(small)) ylabel(0.0 "0" 0.2 "20" 0.4 "40" 0.6 "60" 0.8 "80" 1.0 "100") legend(order(1 "Within 3 weeks" 2 "Within 6 weeks" 3 "More than 6 weeks")) title("Time from referral to rheumatology assessment") name(regional_qs2_bar, replace)
graph export "$projectdir/output/figures/regional_qs2_bar_overall.svg", replace

*By individuals years, merged into one graph

*keep if appt_year==1 | appt_year==2 | appt_year==3 | appt_year==4 | appt_year==5
*lab define appt_year 1 "Year 1" 2 "Year 2" 3 "Year 3" 4 "Year 4" 5 "Year 5", modify
*lab val appt_year appt_year

graph hbar (mean) qs2_0 (mean) qs2_1 (mean) qs2_2, over(appt_year, gap(20) label(labsize(*0.65))) over(nuts_region, gap(60) label(labsize(*0.8))) stack ytitle(Percentage of patients)  ytitle(, size(small)) ylabel(0.0 "0" 0.2 "20" 0.4 "40" 0.6 "60" 0.8 "80" 1.0 "100") legend(order(1 "Within 3 weeks" 2 "Within 6 weeks" 3 "More than 6 weeks")) title("Time from referral to rheumatology assessment") name(regional_qs2_bar_merged, replace) 
graph export "$projectdir/output/figures/regional_qs2_bar_GP_merged.svg", width(12in) replace

*By individual years
forvalues i = 1/$max_year {
    local start = $base_year + `i' - 1
	di "`start'"
    local end = `start' + 1

	preserve
	keep if appt_year==`i'
	
	graph hbar (mean) qs2_0 (mean) qs2_1 (mean) qs2_2, over(nuts_region, relabel(1 "National")) stack ytitle(Percentage of patients) ytitle(, size(small)) ylabel(0.0 "0" 0.2 "20" 0.4 "40" 0.6 "60" 0.8 "80" 1.0 "100") legend(order(1 "Within 3 weeks" 2 "Within 6 weeks" 3 "More than 6 weeks")) title("Time from referral to rheumatology assessment `i'") name(regional_qs2_bar_`i', replace)
	graph export "$projectdir/output/figures/regional_qs2_bar_`i'.svg", replace
	restore
}
/*
forvalues i = 1/$max_year {
    local start = $base_year + `i' - 1
	di "`start'"
    local end = `start' + 1

	preserve
	keep if appt_year==`i'
	gen qs2_0 =1 if time_hes_rheum_appt<=21 & time_hes_rheum_appt!=.
	recode qs2_0 .=0 if time_hes_rheum_appt!=.
	gen qs2_1 =1 if time_hes_rheum_appt>21 & time_hes_rheum_appt<=42 & time_hes_rheum_appt!=.
	recode qs2_1 .=0 if time_hes_rheum_appt!=.
	gen qs2_2 = 1 if time_hes_rheum_appt>42 & time_hes_rheum_appt!=.
	recode qs2_2 .=0 if time_hes_rheum_appt!=.

	expand=2, gen(copy)
	replace nuts_region = 0 if copy==1
	lab define nuts_region 0 "National" 9 "Yorkshire/Humber", modify
	lab val nuts_region nuts_region

	graph hbar (mean) qs2_0 (mean) qs2_1 (mean) qs2_2, over(nuts_region, relabel(1 "National")) stack ytitle(Percentage of patients) ytitle(, size(small)) ylabel(0.0 "0" 0.2 "20" 0.4 "40" 0.6 "60" 0.8 "80" 1.0 "100") legend(order(1 "Within 3 weeks" 2 "Within 6 weeks" 3 "More than 6 weeks")) title("Time from referral to rheumatology assessment") name(regional_qs2_bar_`i', replace)
	graph export "$projectdir/output/figures/regional_qs2_bar_`i'.svg", replace
	restore
}
*/
//for output checking tables for boxplot - see output/tables/referral_byregion_rounded.csv

/*
/*GP referral performance by ethnicity, merged===========================================================================*/

preserve
keep if appt_year==1 | appt_year==2 | appt_year==3 | appt_year==4
lab define appt_year 1 "Year 1" 2 "Year 2" 3 "Year 3" 4 "Year 4", modify
lab val appt_year appt_year

gen qs2_0 =1 if time_gp_rheum_appt<=21 & time_gp_rheum_appt!=.
recode qs2_0 .=0 if time_gp_rheum_appt!=.
gen qs2_1 =1 if time_gp_rheum_appt>21 & time_gp_rheum_appt<=42 & time_gp_rheum_appt!=.
recode qs2_1 .=0 if time_gp_rheum_appt!=.
gen qs2_2 = 1 if time_gp_rheum_appt>42 & time_gp_rheum_appt!=.
recode qs2_2 .=0 if time_gp_rheum_appt!=.

expand=2, gen(copy)
replace ethnicity = 0 if copy==1  

lab define ethnicity 0 "Overall" 2 "Asian", modify
lab val ethnicity ethnicity

graph hbar (mean) qs2_0 (mean) qs2_1 (mean) qs2_2, over(appt_year, gap(20) label(labsize(*0.9))) over(ethnicity, gap(60) label(labsize(*0.8))) stack ytitle(Percentage of patients)  ytitle(, size(small)) ylabel(0.0 "0" 0.2 "20" 0.4 "40" 0.6 "60" 0.8 "80" 1.0 "100") legend(order(1 "Within 3 weeks" 2 "Within 6 weeks" 3 "More than 6 weeks")) title("Time from referral to rheumatology assessment, by ethnicity") name(regional_qs2_bar_GP_ethnicity, replace)
graph export "$projectdir/output/figures/regional_qs2_bar_GP_ethnicity.svg", replace
restore

/*GP referral performance by IMD quintile, merged===========================================================================*/

preserve
keep if appt_year==1 | appt_year==2 | appt_year==3 | appt_year==4
lab define appt_year 1 "Year 1" 2 "Year 2" 3 "Year 3" 4 "Year 4", modify
lab val appt_year appt_year

gen qs2_0 =1 if time_gp_rheum_appt<=21 & time_gp_rheum_appt!=.
recode qs2_0 .=0 if time_gp_rheum_appt!=.
gen qs2_1 =1 if time_gp_rheum_appt>21 & time_gp_rheum_appt<=42 & time_gp_rheum_appt!=.
recode qs2_1 .=0 if time_gp_rheum_appt!=.
gen qs2_2 = 1 if time_gp_rheum_appt>42 & time_gp_rheum_appt!=.
recode qs2_2 .=0 if time_gp_rheum_appt!=.

expand=2, gen(copy)
replace imd = 0 if copy==1  

lab define imd 0 "Overall" 1 "1st Quintile" 2 "2nd Quintile" 3 "3rd Quintile" 4 "4th Quintile" 5 "5th Quintile", modify
lab val imd imd

graph hbar (mean) qs2_0 (mean) qs2_1 (mean) qs2_2, over(appt_year, gap(20) label(labsize(*0.9))) over(imd, gap(60) label(labsize(*0.8)) relabel(2 `" "1st Quintile" "(most deprived)" "' 6 `" "5th Quintile" "(least deprived)" "')) stack ytitle(Percentage of patients)  ytitle(, size(small)) ylabel(0.0 "0" 0.2 "20" 0.4 "40" 0.6 "60" 0.8 "80" 1.0 "100") legend(order(1 "Within 3 weeks" 2 "Within 6 weeks" 3 "More than 6 weeks")) title("Time from referral to rheumatology assessment, by IMD quintile") name(regional_qs2_bar_GP_imd, replace)
graph export "$projectdir/output/figures/regional_qs2_bar_GP_imd.svg", replace
restore

*Time from first rheumatology appointment to first csDMARD shared care by region ======================================================*/

use "$projectdir/output/data/file_eia_all.dta", clear

**For RA, PsA and Undiff IA patients combined (not including AxSpA - low counts)
keep if (eia_diagnosis==1 | eia_diagnosis==2 | eia_diagnosis==4)

gen csdmard_0 =1 if time_to_csdmard<=90 & time_to_csdmard!=.
recode csdmard_0 .=0
gen csdmard_1 =1 if time_to_csdmard>90 & time_to_csdmard<=180 & time_to_csdmard!=.
recode csdmard_1 .=0
gen csdmard_2 = 1 if time_to_csdmard>180 | time_to_csdmard==.
recode csdmard_2 .=0 

expand=2, gen(copy)
replace nuts_region = 0 if copy==1  
lab define nuts_region 0 "National" 9 "Yorkshire/Humber", modify
lab val nuts_region nuts_region

graph hbar (mean) csdmard_0 (mean) csdmard_1 (mean) csdmard_2, over(nuts_region, relabel(1 "National")) stack ytitle(Percentage of patients) ytitle(, size(small)) ylabel(0.0 "0" 0.2 "20" 0.4 "40" 0.6 "60" 0.8 "80" 1.0 "100") legend(order(1 "Within 3 months" 2 "Within 6 months" 3 "None within 6 months")) title("Time to first csDMARD in primary care, overall") name(regional_csdmard_bar, replace)
graph export "$projectdir/output/figures/regional_csdmard_bar_overall.svg", replace

*By individuals years, merged into one graph

*keep if appt_year==1 | appt_year==2 | appt_year==3 | appt_year==4 | appt_year==5
*lab define appt_year 1 "Year 1" 2 "Year 2" 3 "Year 3" 4 "Year 4" 5 "Year 5", modify
*lab val appt_year appt_year

graph hbar csdmard_0 (mean) csdmard_1 (mean) csdmard_2, over(appt_year, gap(20) label(labsize(*0.65))) over(nuts_region, gap(60) label(labsize(*0.8))) stack ytitle(Percentage of patients) ytitle(, size(small)) ylabel(0.0 "0" 0.2 "20" 0.4 "40" 0.6 "60" 0.8 "80" 1.0 "100") legend(order(1 "Within 3 months" 2 "Within 6 months" 3 "None within 6 months")) title("Time to first csDMARD in primary care") name(regional_csdmard_bar_merged, replace)
graph export "$projectdir/output/figures/regional_csdmard_bar_merged.svg", width(12in) replace

*By individual years
forvalues i = 1/$max_year {
    local start = $base_year + `i' - 1
	di "`start'"
    local end = `start' + 1

	preserve
	keep if appt_year==`i'
	
	graph hbar (mean) csdmard_0 (mean) csdmard_1 (mean) csdmard_2, over(nuts_region, relabel(1 "National")) stack ytitle(Percentage of patients) ytitle(, size(small)) ylabel(0.0 "0" 0.2 "20" 0.4 "40" 0.6 "60" 0.8 "80" 1.0 "100") legend(order(1 "Within 3 months" 2 "Within 6 months" 3 "None within 6 months")) title("Time to first csDMARD in primary care `i'") name(regional_csdmard_bar_`i', replace)
graph export "$projectdir/output/figures/regional_csdmard_bar_`i'.svg", replace
	restore
}

//for output checking table for boxplot - see output/tables/drug_byyearandregion_rounded.csv

/*
/*csDMARD shared care performance by ethnicity, merged===========================================================================*/

preserve
keep if appt_year==1 | appt_year==2 | appt_year==3 | appt_year==4
lab define appt_year 1 "Year 1" 2 "Year 2" 3 "Year 3" 4 "Year 4", modify
lab val appt_year appt_year

gen csdmard_0 =1 if time_to_csdmard<=90 & time_to_csdmard!=.
recode csdmard_0 .=0
gen csdmard_1 =1 if time_to_csdmard>90 & time_to_csdmard<=180 & time_to_csdmard!=.
recode csdmard_1 .=0
gen csdmard_2 = 1 if time_to_csdmard>180 | time_to_csdmard==.
recode csdmard_2 .=0 

expand=2, gen(copy)
replace ethnicity = 0 if copy==1  

lab define ethnicity 0 "Overall" 2 "Asian", modify
lab val ethnicity ethnicity

graph hbar (mean) csdmard_0 (mean) csdmard_1 (mean) csdmard_2, over(appt_year, gap(20) label(labsize(*0.9))) over(ethnicity, gap(60) label(labsize(*0.8))) stack ytitle(Percentage of patients)  ytitle(, size(small)) ylabel(0.0 "0" 0.2 "20" 0.4 "40" 0.6 "60" 0.8 "80" 1.0 "100") legend(order(1 "Within 3 months" 2 "Within 6 months" 3 "None within 6 months")) title("Time to first csDMARD in primary care, by ethnicity") name(regional_csdmard_bar_ethnicity, replace)
graph export "$projectdir/output/figures/regional_csdmard_bar_ethnicity.svg", replace
restore

/*csDMARD shared care performance by IMD quintile, merged===========================================================================*/

preserve
keep if appt_year==1 | appt_year==2 | appt_year==3 | appt_year==4
lab define appt_year 1 "Year 1" 2 "Year 2" 3 "Year 3" 4 "Year 4", modify
lab val appt_year appt_year

gen csdmard_0 =1 if time_to_csdmard<=90 & time_to_csdmard!=.
recode csdmard_0 .=0
gen csdmard_1 =1 if time_to_csdmard>90 & time_to_csdmard<=180 & time_to_csdmard!=.
recode csdmard_1 .=0
gen csdmard_2 = 1 if time_to_csdmard>180 | time_to_csdmard==.
recode csdmard_2 .=0 

expand=2, gen(copy)
replace imd = 0 if copy==1  

lab define imd 0 "Overall" 1 "1st Quintile" 2 "2nd Quintile" 3 "3rd Quintile" 4 "4th Quintile" 5 "5th Quintile", modify
lab val imd imd

graph hbar (mean) csdmard_0 (mean) csdmard_1 (mean) csdmard_2, over(appt_year, gap(20) label(labsize(*0.9))) over(imd, gap(60) label(labsize(*0.8)) relabel(2 `" "1st Quintile" "(most deprived)" "' 6 `" "5th Quintile" "(least deprived)" "')) stack ytitle(Percentage of patients)  ytitle(, size(small)) ylabel(0.0 "0" 0.2 "20" 0.4 "40" 0.6 "60" 0.8 "80" 1.0 "100") legend(order(1 "Within 3 months" 2 "Within 6 months" 3 "None within 6 months")) title("Time to first csDMARD in primary care, by IMD quintile") name(regional_csdmard_bar_imd, replace)
graph export "$projectdir/output/figures/regional_csdmard_bar_imd.svg", replace
restore

*/

/*GP referral performance by region, all years; low capture of rheum referrals presently, therefore using last GP appt as proxy measure currently - see below===========================================================================*

***Restrict all analyses below to patients with rheum appt, GP appt and 6m follow-up and registration (changed from 12m requirement, for purposes of OpenSAFELY report)
keep if has_6m_post_appt==1

preserve
gen qs1_0 =1 if time_gp_rheum_ref_appt<=3 & time_gp_rheum_ref_appt!=.
recode qs1_0 .=0 if time_gp_rheum_ref_appt!=.
gen qs1_1 =1 if time_gp_rheum_ref_appt>3 & time_gp_rheum_ref_appt<=7 & time_gp_rheum_ref_appt!=.
recode qs1_1 .=0 if time_gp_rheum_ref_appt!=.
gen qs1_2 = 1 if time_gp_rheum_ref_appt>7 & time_gp_rheum_ref_appt!=.
recode qs1_2 .=0 if time_gp_rheum_ref_appt!=.

expand=2, gen(copy)
replace nuts_region = 0 if copy==1  

graph hbar (mean) qs1_0 (mean) qs1_1 (mean) qs1_2, over(nuts_region, relabel(1 "National")) stack ytitle(Percentage of patients) ytitle(, size(small)) ylabel(0.0 "0" 0.2 "20" 0.4 "40" 0.6 "60" 0.8 "80" 1.0 "100") legend(order(1 "Within 3 days" 2 "Within 7 days" 3 "More than 7 days")) title("Time to rheumatology referral") name(regional_qs1_bar, replace)
graph export "$projectdir/output/figures/regional_qs1_bar_overall.svg", replace
restore

/*GP referral performance by region, Apr 2019 to Apr 2020; low capture of rheum referrals presently, therefore using last GP appt as proxy measure currently - see below===========================================================================*/

preserve
keep if appt_year==1
gen qs1_0 =1 if time_gp_rheum_ref_appt<=3 & time_gp_rheum_ref_appt!=.
recode qs1_0 .=0 if time_gp_rheum_ref_appt!=.
gen qs1_1 =1 if time_gp_rheum_ref_appt>3 & time_gp_rheum_ref_appt<=7 & time_gp_rheum_ref_appt!=.
recode qs1_1 .=0 if time_gp_rheum_ref_appt!=.
gen qs1_2 = 1 if time_gp_rheum_ref_appt>7 & time_gp_rheum_ref_appt!=.
recode qs1_2 .=0 if time_gp_rheum_ref_appt!=.

expand=2, gen(copy)
replace nuts_region = 0 if copy==1  

graph hbar (mean) qs1_0 (mean) qs1_1 (mean) qs1_2, over(nuts_region, relabel(1 "National")) stack ytitle(Percentage of patients) ytitle(, size(small)) ylabel(0.0 "0" 0.2 "20" 0.4 "40" 0.6 "60" 0.8 "80" 1.0 "100") legend(order(1 "Within 3 days" 2 "Within 7 days" 3 "More than 7 days")) title("Time to rheumatology referral") name(regional_qs1_bar, replace)
graph export "$projectdir/output/figures/regional_qs1_bar_2019.svg", replace
restore

/*GP referral performance by region, Apr 2020 to Apr 2021; low capture of rheum referrals presently, therefore using last GP appt as proxy measure currently - see below===========================================================================*/

preserve
keep if appt_year==2
gen qs1_0 =1 if time_gp_rheum_ref_appt<=3 & time_gp_rheum_ref_appt!=.
recode qs1_0 .=0 if time_gp_rheum_ref_appt!=.
gen qs1_1 =1 if time_gp_rheum_ref_appt>3 & time_gp_rheum_ref_appt<=7 & time_gp_rheum_ref_appt!=.
recode qs1_1 .=0 if time_gp_rheum_ref_appt!=.
gen qs1_2 = 1 if time_gp_rheum_ref_appt>7 & time_gp_rheum_ref_appt!=.
recode qs1_2 .=0 if time_gp_rheum_ref_appt!=.

expand=2, gen(copy)
replace nuts_region = 0 if copy==1  

graph hbar (mean) qs1_0 (mean) qs1_1 (mean) qs1_2, over(nuts_region, relabel(1 "National")) stack ytitle(Percentage of patients) ytitle(, size(small)) ylabel(0.0 "0" 0.2 "20" 0.4 "40" 0.6 "60" 0.8 "80" 1.0 "100") legend(order(1 "Within 3 days" 2 "Within 7 days" 3 "More than 7 days")) title("Time to rheumatology referral") name(regional_qs1_bar, replace)
graph export "$projectdir/output/figures/regional_qs1_bar_2020.svg", replace
restore

/*GP referral performance by region, Apr 2021 to Apr 2022; low capture of rheum referrals presently, therefore using last GP appt as proxy measure currently - see below===========================================================================*/

preserve
keep if appt_year==3
gen qs1_0 =1 if time_gp_rheum_ref_appt<=3 & time_gp_rheum_ref_appt!=.
recode qs1_0 .=0 if time_gp_rheum_ref_appt!=.
gen qs1_1 =1 if time_gp_rheum_ref_appt>3 & time_gp_rheum_ref_appt<=7 & time_gp_rheum_ref_appt!=.
recode qs1_1 .=0 if time_gp_rheum_ref_appt!=.
gen qs1_2 = 1 if time_gp_rheum_ref_appt>7 & time_gp_rheum_ref_appt!=.
recode qs1_2 .=0 if time_gp_rheum_ref_appt!=.

expand=2, gen(copy)
replace nuts_region = 0 if copy==1  

graph hbar (mean) qs1_0 (mean) qs1_1 (mean) qs1_2, over(nuts_region, relabel(1 "National")) stack ytitle(Percentage of patients) ytitle(, size(small)) ylabel(0.0 "0" 0.2 "20" 0.4 "40" 0.6 "60" 0.8 "80" 1.0 "100") legend(order(1 "Within 3 days" 2 "Within 7 days" 3 "More than 7 days")) title("Time to rheumatology referral") name(regional_qs1_bar, replace)
graph export "$projectdir/output/figures/regional_qs1_bar_2021.svg", replace
restore

*/

log off