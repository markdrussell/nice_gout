from ehrql import create_dataset, days, months, years, case, when, get_parameter, minimum_of, maximum_of
from ehrql.tables.tpp import patients, medications, practice_registrations, clinical_events, apcs, addresses, ethnicity_from_sus 
from ehrql.codes import ICD10Code
from datetime import date, datetime
from functools import reduce
import codelists_ehrQL as codelists
import json

# Read parameters from project.yaml
studystart_date = get_parameter("studystart_date")
studyend_date = get_parameter("studyend_date")
studyfup_date = get_parameter("studyfup_date")

# Read other parameters and ensure they are lists
def ensure_list(x):
    if isinstance(x, list):
        return x
    if x is None:
        return []
    if isinstance(x, str):
        s = x.strip()
        if not s:
            return []
        if s[0] == "[" and s[-1] == "]":
            try:
                v = json.loads(s)
                if isinstance(v, list):
                    return [str(t).strip() for t in v if str(t).strip()]
            except Exception:
                s = s[1:-1]
        if "," in s:
            return [p.strip().strip('\'"') for p in s.split(",") if p.strip().strip('\'"')]
        if " " in s:
            return [p for p in s.split() if p]
        return [s]
    s = str(x).strip()
    return [s] if s else []

comorbidities_list = ensure_list(get_parameter("comorbidities_list"))
bloods_list = ensure_list(get_parameter("bloods_list"))

# Medications list abbreviated for reference population (amend as necessary or pull in from project.yaml)
medications_list = ensure_list("diuretic")

dataset = create_dataset()
dataset.configure_dummy_data(population_size=1000)

# Function to identify code on or before cohort entry date (can define as first or last code in period)
def code_before_entry(dx_codelist):
    return clinical_events.where(
        clinical_events.snomedct_code.is_in(dx_codelist)
    ).where(
        clinical_events.date.is_on_or_before(dataset.cohort_entry_date)
    ).sort_by(
        clinical_events.date
    )

# Function to identify last code in primary care up to X [specify] months before cohort entry date
def code_before_entry_window(dx_codelist, pre_time_window):
    return clinical_events.where(
        clinical_events.snomedct_code.is_in(dx_codelist)
    ).where(
        clinical_events.date.is_on_or_between((dataset.cohort_entry_date - months(pre_time_window)), dataset.cohort_entry_date)
    ).sort_by(
        clinical_events.date
    ).last_for_patient()

# Function to identify recurrent blood tests within a specified time period (only those with associated numeric values)
def blood_before_entry(dx_codelist, time_window):
    return clinical_events.where(
        clinical_events.snomedct_code.is_in(dx_codelist)
    ).where(
        clinical_events.date.is_on_or_between((dataset.cohort_entry_date - months(time_window)), (dataset.cohort_entry_date + months(time_window)))
    ).where(
        clinical_events.numeric_value.is_not_null() & (clinical_events.numeric_value > 0) & (clinical_events.numeric_value < 3000)
    ).sort_by(
        clinical_events.date
    ).last_for_patient()

# Function to identify last prescription issued up to X (specify) months before cohort entry date
def medication_before_entry(dmd_codelist, pre_time_window):
    return medications.where(
        medications.dmd_code.is_in(dmd_codelist)
    ).where(
        medications.date.is_on_or_between((dataset.cohort_entry_date - months(pre_time_window)), dataset.cohort_entry_date)
    ).sort_by(
        medications.date
    ).last_for_patient()

# Any practice registration before study end date
registration_in_window = (
    practice_registrations.where(practice_registrations.start_date <= studyend_date)
    .except_where(practice_registrations.end_date < studystart_date)
)

# Keep only overlapping episodes that (observably) last at least 12 months
registration_in_window_12m = (registration_in_window
    .where(
        practice_registrations.start_date + months(12) <= studyend_date
    )
    .except_where(
        practice_registrations.end_date < (practice_registrations.start_date + months(12))
    )
)

dataset.any_registration_12m = registration_in_window_12m.exists_for_patient()

dataset.first_registration_12m = (
    registration_in_window_12m.sort_by(practice_registrations.start_date)
    .first_for_patient()
    .start_date
)

dataset.region = (
    registration_in_window_12m.sort_by(practice_registrations.start_date)
    .first_for_patient()
    .practice_nuts1_region_name
)

# Cohort entry date
dataset.cohort_entry_date = maximum_of(studystart_date, (dataset.first_registration_12m + months(12)))

# Define sex
dataset.sex = patients.sex

# Date of death
dataset.date_of_death = patients.date_of_death

dataset.alive_at_entry = (
    patients.date_of_death.is_null() | patients.date_of_death.is_after(dataset.cohort_entry_date)
)

# Age at index date
dataset.age = patients.age_on(dataset.cohort_entry_date)

# Define patient ethnicity
latest_ethnicity_code = (
    clinical_events.where(clinical_events.snomedct_code.is_in(codelists.ethnicity_codes))
    .where(clinical_events.date.is_on_or_before(studyend_date))
    .sort_by(clinical_events.date)
    .last_for_patient().snomedct_code.to_category(codelists.ethnicity_codes)
)

# Extract ethnicity from SUS records if it isn't present in primary care data 
ethnicity_sus = ethnicity_from_sus.code

dataset.ethnicity = case(
    when((latest_ethnicity_code == "1") | ((latest_ethnicity_code.is_null()) & (ethnicity_sus.is_in(["A", "B", "C"])))).then("White"),
    when((latest_ethnicity_code == "2") | ((latest_ethnicity_code.is_null()) & (ethnicity_sus.is_in(["D", "E", "F", "G"])))).then("Mixed"),
    when((latest_ethnicity_code == "3") | ((latest_ethnicity_code.is_null()) & (ethnicity_sus.is_in(["H", "J", "K", "L"])))).then("Asian or Asian British"),
    when((latest_ethnicity_code == "4") | ((latest_ethnicity_code.is_null()) & (ethnicity_sus.is_in(["M", "N", "P"])))).then("Black or Black British"),
    when((latest_ethnicity_code == "5") | ((latest_ethnicity_code.is_null()) & (ethnicity_sus.is_in(["R", "S"])))).then("Chinese or Other Ethnic Groups"),
    otherwise="Unknown", 
) 

# Define patient IMD at cohort entry
imd = addresses.for_patient_on(dataset.cohort_entry_date).imd_rounded

dataset.imd_quintile = case(
    when((imd >= 0) & (imd < int(32844 * 1 / 5))).then("1 (most deprived)"),
    when(imd < int(32844 * 2 / 5)).then("2"),
    when(imd < int(32844 * 3 / 5)).then("3"),
    when(imd < int(32844 * 4 / 5)).then("4"),
    when(imd < int(32844 * 5 / 5)).then("5 (least deprived)"),
    otherwise="Unknown",
)

# Define population
dataset.define_population(
    ((dataset.age >= 18) & (dataset.age <= 110))
    & dataset.sex.is_in(["male", "female"])
    & dataset.alive_at_entry
    & dataset.any_registration_12m
    & dataset.cohort_entry_date.is_not_null()
)

# Date of diagnosis for comorbidities (first recorded code before study end date)
for comorbidity in comorbidities_list:
    comorbidity_codelist = getattr(codelists, f"{comorbidity}_codes")
    dataset.add_column(f"{comorbidity}_date", code_before_entry(comorbidity_codelist).first_for_patient().date)

# Body mass index (last recorded code up to 120 [specify] months before entry date)
dataset.bmi_value = code_before_entry_window(codelists.bmi_codes, 120).numeric_value
dataset.bmi_date = code_before_entry_window(codelists.bmi_codes, 120).date

# Smoking status (last recorded code before primary diagnosis)
dataset.most_recent_smoking_code=clinical_events.where(
        clinical_events.ctv3_code.is_in(codelists.clear_smoking_codes)
    ).where(
        clinical_events.date.is_on_or_before(dataset.cohort_entry_date)
    ).sort_by(
        clinical_events.date
    ).last_for_patient().ctv3_code.to_category(codelists.clear_smoking_codes)

def filter_codes_by_category(codelist, include):
    return {k:v for k,v in codelist.items() if v in include}

dataset.ever_smoked=clinical_events.where(
        clinical_events.ctv3_code.is_in(filter_codes_by_category(codelists.clear_smoking_codes, include=["S", "E"]))
    ).where(
        clinical_events.date.is_on_or_before(dataset.cohort_entry_date)
    ).exists_for_patient()

dataset.smoking_status=case(
    when(dataset.most_recent_smoking_code == "S").then("S"),
    when((dataset.most_recent_smoking_code == "E") | ((dataset.most_recent_smoking_code == "N") & (dataset.ever_smoked == True))).then("E"),
    when((dataset.most_recent_smoking_code == "N") & (dataset.ever_smoked == False)).then("N"),
    otherwise="M"
)

# Blood test performed within X [specify] months before or after cohort entry date 
for blood in bloods_list:
    blood_codelist = getattr(codelists, f"{blood}_codes")
    ## Last blood test in the X (specify) months before cohort entry date
    dataset.add_column(f"{blood}_bl_date", blood_before_entry(blood_codelist, 2).date)
    dataset.add_column(f"{blood}_bl_value", blood_before_entry(blood_codelist, 2).numeric_value)


# Medication use within X [specify] months before cohort entry date
for medication in medications_list:
    medication_codelist = getattr(codelists, f"{medication}_codes")
    ## Last prescription in the X (specify) months before cohort entry date
    dataset.add_column(f"{medication}_bl_date", medication_before_entry(medication_codelist, 6).date)