#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import sys

if __name__ == "__main__":
    subject_path = sys.argv[1]
    visit_path = sys.argv[2]
    stool_path = sys.argv[3]
    preparation_path = sys.argv[4]
    mgen_path = sys.argv[5]
    mtab_path = sys.argv[6]

    meta = pd.read_csv(sys.stdin, low_memory=False)

    # Extract table: 'subject'.

    column_rename = {
        "Participant ID": "subject_id",
        "Age at diagnosis": "ibd_diagnosis_at_age",
        "consent_age": "age",
        "site_name": "site",
        "Has the subject had a prior abdominal surgery (other)?": "history_other_abdominal_surgery",
        "Education Level": "education",
        "Occupation": "occupation",
        "Has the subject had a tonsillectomy?": "history_tonsillectomy",
        "diagnosis": "ibd_diagnosis",
        "Did you grow up on a farm?": "history_childhood_farm",
        "Did you attend daycare as a child?": "history_childhood_daycare",
        "Were you exposed to cigarette smoke as a child?": "history_childhood_environmental_tobacco_smoke",
        "Were you born prematurely (more than 3 weeks early)?": "history_birth_premature",
        "Were you born in a hospital?": "history_birth_hospital",
        "Were you born via C-section?": "history_birth_csection",
        "Were you breastfed as an infant?": "history_childhood_breastfed",
        "Were you treated with antibiotics before the age of one?": "history_childhood_antibiotics_pre_age_one",
        "Were you hospitalized before the age of five?": "history_childhood_antibiotics_pre_age_five",
        "Did you have pets growing up?": "history_childhood_pets",
        "race": "race",
        "sex": "sex",
        "smoking status": "status_smoker",
        "Number years smoked": "status_smoker_years",
        "Age when started smoking": "history_smoker_age_at_start",
        "How many cigarettes/cigars/etc. do you smoke per day?": "status_smoker_number_per_day",
        "Height": "baseline_height",
        "Weight.1": "baseline_weight",
    }

    column_replace = {
        "history_other_abdominal_surgery": {
            "Yes": True,
            "No": False,
            "Not Sure": np.nan,
        },
        "education": {"Unknown/Not Reported": np.nan},
        "occupation": {"Unknown/Not Reported": np.nan},
        "history_tonsillectomy": {"Yes": True, "No": False, "Not Sure": np.nan},
        "history_childhood_farm": {"Yes": True, "No": False, "Not Sure": np.nan},
        "history_childhood_daycare": {"Yes": True, "No": False, "Not Sure": np.nan},
        "history_childhood_environmental_tobacco_smoke": {
            "Yes": True,
            "No": False,
            "Not Sure": np.nan,
        },
        "history_birth_premature": {"Yes": True, "No": False, "Not Sure": np.nan},
        "history_birth_hospital": {
            "Yes": True,
            "No": False,
            "Not sure": np.nan,
        },  # Notice the different capitalization of 'Not sure'.
        "history_birth_csection": {
            "Yes": True,
            "No": False,
            "Not sure": np.nan,
        },  # Notice the different capitalization of 'Not sure'.
        "history_childhood_breastfed": {"Yes": True, "No": False, "Not Sure": np.nan},
        "history_childhood_antibiotics_pre_age_one": {
            "Yes": True,
            "No": False,
            "Not sure": np.nan,
        },  # Notice the different capitalization of 'Not sure'.
        "history_childhood_antibiotics_pre_age_five": {
            "Yes": True,
            "No": False,
            "Not Sure": np.nan,
        },
        "history_childhood_pets": {"Yes": True, "No": False, "Not Sure": np.nan},
    }

    _subject = meta[column_rename.keys()].rename(  # Select columns of interest.
        columns=column_rename
    )  # Rename columns.

    # Replace malformed values
    for col in column_replace:
        _subject[col].replace(column_replace[col], inplace=True)

    # Fill NaNs when other valid values.
    _subject = _subject.groupby("subject_id").apply(
        lambda d: d.fillna(method="bfill").fillna(method="ffill")
    )

    _subject["history_smoker"] = _subject.status_smoker.replace(
        {"Current smoker": True, "Former smoker": True, "Never smoked": False}
    ).astype(float)
    _subject["status_smoker"] = _subject.status_smoker.replace(
        {"Current smoker": True, "Former smoker": False, "Never smoked": False}
    ).astype(float)

    n_unique_entries = _subject.drop_duplicates().shape[0]
    n_unique_subjects = _subject.subject_id.unique().shape[0]
    assert n_unique_subjects == n_unique_entries, (n_unique_subjects, n_unique_entries)

    print(f"Found {n_unique_subjects} unique subjects.")

    subject = _subject.drop_duplicates().set_index("subject_id")
    assert subject.index.is_unique

    # Extract table: 'visit'.

    column_rename = {
        "site_sub_coll": "visit_id",
        "Participant ID": "subject_id",
        "week_num": "week_number",
        "date_of_receipt": "visit_date",
        "visit_num": "visit_number",  # Visit number 1 CAN be after other visits; this describes the procedures done.
        "Soft drinks, tea or coffee with sugar (corn syrup, maple syrup, cane sugar, etc)": "diet_sugary_beverage",
        "Diet soft drinks, tea or coffee with sugar (Stevia, Equal, Splenda etc)": "diet_artificial_sweetener_beverage",
        "Fruit juice (orange, apple, cranberry, prune etc.)": "diet_fruit_juice",
        "Water": "diet_water",
        "Alcohol (beer, brandy, spirits, hard liquor, wine, aperitif, etc.)": "diet_alcohol",
        "Yogurt or other foods containing active bacterial cultures (kefir, sauerkraut)": "diet_live_bacteria",
        "Dairy (milk, cream, ice cream, cheese, cream cheese)": "diet_dairy",
        "Probiotic": "diet_probiotic",
        "Fruits (no juice) (Apples, raisins, bananas, oranges, strawberries, blueberries": "diet_whole_fruits",
        "Vegetables (salad, tomatoes, onions, greens, carrots, peppers, green beans, etc)": "diet_whole_vegetables",
        "Beans (tofu, soy, soy burgers, lentils, Mexican beans, lima beans etc)": "diet_beans",
        "Whole grains (wheat, oats, brown rice, rye, quinoa, wheat bread, wheat pasta)": "diet_whole_grains",
        "Starch (white rice, bread, pizza, potatoes, yams, cereals, pancakes, etc.)": "diet_starch",
        "Eggs": "diet_eggs",
        "Processed meat (other red or white meat such as lunch meat, ham, salami, bologna": "diet_meats_processed",
        "Red meat (beef, hamburger, pork, lamb)": "diet_meats_red",
        "White meat (chicken, turkey, etc.)": "diet_meats_white",
        "Shellfish (shrimp, lobster, scallops, etc.)": "diet_meats_shellfish",
        "Fish (fish nuggets, breaded fish, fish cakes, salmon, tuna, etc.)": "diet_meats_fish",
        "Sweets (pies, jam, chocolate, cake, cookies, etc.)": "diet_sweets",
        "Antibiotics": "status_antibiotics",
        "Chemotherapy": "status_chemotherapy",
        "Immunosuppressants (e.g. oral corticosteroids)": "status_immunosuppressants",
        "2) In the past 2 weeks, have you undergone a colonoscopy or other procedure": "status_recent_colonoscopy",
        "3) In the past 2 weeks, have you used an oral contrast": "status_oral_contrast",
        "4) In the past 2 weeks, have you had diarrhea?": "status_diarrhea",
        "5) In the past 2 weeks, have you been hospitalized?": "status_recent_hospitalization",
        "6) Have you ever had bowel surgery?": "status_ever_had_bowel_surgery",  # This is with visits because some individuals had surgery during the study.
        "Tea or coffee no sugar and no sugar replacement": "diet_tea_or_coffee",
        "General wellbeing": "status_general_wellbeing",
        "sccai": "status_sccai",
        "Urgency of defecation": "status_defecation_urgency",
        "Blood in the stool": "status_blood_in_stool",
        "Abdominal pain": "status_abdominal_pain",
        "Number of liquid or very soft stools in the past 24 hours:": "status_number_soft_stools",
        "Abdominal mass": "status_abdominal_mass",
        "Arthralgia": "status_arthralgia",
        "hbi": "status_hbi",
        "CRP (mg/L)": "status_crp",
        "ESR (mm/hr)": "status_esr",
        "Weight": "status_weight",
    }

    column_replace = {
        "status_antibiotics": {"Yes": True, "No": False},
        "status_chemotherapy": {"Yes": True, "No": False},
        "status_immunosuppressants": {"Yes": True, "No": False},
        "status_recent_colonoscopy": {"Yes": True, "No": False},
        "status_oral_contrast": {"Yes": True, "No": False},
        "status_diarrhea": {"Yes": True, "No": False},
        "status_recent_hospitalization": {"Yes": True, "No": False},
        "status_ever_had_bowel_surgery": {"Yes": True, "No": False},
        "status_arthralgia": {"Yes": True, "No": False},
        "status_weight": {999: np.nan},
        "status_crp": {999: np.nan},
        "status_esr": {999: np.nan},
        "status_number_soft_stools": {999: np.nan},
    }

    _visit = meta[column_rename.keys()].rename(  # Select columns of interest.
        columns=column_rename
    )  # Rename columns.

    # Replace malformed values
    for col in column_replace:
        _visit[col].replace(column_replace[col], inplace=True)

    # Fix week numbers
    # Week numbers <= zero are sometimes mixed in with other numbers for a single visit number.
    def _ifthen(bool, iftrue, iffalse):
        if bool:
            return iftrue
        else:
            return iffalse

    _visit["week_number"] = _visit["week_number"].map(
        lambda x: _ifthen(x <= 0, np.nan, x)
    )

    # Fill NaNs when other valid values.
    _visit = _visit.groupby("visit_id").apply(
        lambda d: d.fillna(method="bfill").fillna(method="ffill")
    )

    n_unique_entries = _visit.drop_duplicates().shape[0]
    n_unique_visits = _visit.visit_id.unique().shape[0]
    assert n_unique_visits == n_unique_entries, (n_unique_visits, n_unique_entries)

    print(f"Found {n_unique_visits} unique visits.")

    visit = _visit.drop_duplicates().set_index("visit_id")
    assert visit.index.is_unique
    assert visit.subject_id.isin(subject.index).all()

    # Extract table: 'stool'.

    column_rename = {
        "External ID": "stool_id",
        "site_sub_coll": "visit_id",
        "fecalcal": "fecal_calprotectin",
    }

    _stool = meta[
        meta.data_type.isin(  # Exclude host_transcriptomics, methylome, host_genome, serology
            [
                "viromics",
                "metabolomics",
                "metagenomics",
                "proteomics",
                "metagenomics",
                "metatranscriptomics",
                "stool_16S",
            ]
        )
    ][
        column_rename.keys()
    ].rename(  # Select columns of interest.
        columns=column_rename
    )  # Rename columns.

    # Remove suffixes _P, _TR.  What do these mean?
    _stool.stool_id = _stool.stool_id.str.replace("_TR", "").str.replace("_P", "")

    # Fill NaNs when other valid values.
    _stool = _stool.groupby("stool_id").apply(
        lambda d: d.fillna(method="bfill").fillna(method="ffill")
    )

    n_unique_entries = _stool.drop_duplicates().shape[0]
    n_unique_stools = _stool.stool_id.unique().shape[0]
    assert n_unique_stools == n_unique_entries, (n_unique_stools, n_unique_entries)

    print(f"Found {n_unique_stools} unique stools.")

    stool = _stool.drop_duplicates().set_index("stool_id")
    assert stool.index.is_unique
    assert stool.visit_id.isin(visit.index).all()

    # Extract table: 'preparation'.  # Because there are sometimes multiple libraries for the same target library.

    column_rename = {
        "Project": "project_name",
        "External ID": "external_id",
        "data_type": "library_type",
    }

    _preparation = meta[
        meta.data_type.isin(  # Exclude host_transcriptomics, methylome, host_genome
            [
                "viromics",
                "metabolomics",
                "metagenomics",
                "proteomics",
                "metagenomics",
                "metatranscriptomics",
                "stool_16S",
            ]
        )
    ][
        column_rename.keys()
    ].rename(  # Select columns of interest.
        columns=column_rename
    )  # Rename columns.

    # Remove suffixes _P, _TR.  What do these mean?
    _preparation["stool_id"] = _preparation.external_id.str.replace(
        "_TR", ""
    ).str.replace("_P", "")
    # Construct the primary key
    _preparation["preparation_id"] = (
        _preparation.stool_id + "_" + _preparation.project_name
    )

    # Drop the external_id, this'll be used in a later table
    _preparation.drop(columns=["external_id"], inplace=True)

    # Fill NaNs when other valid values.
    _preparation = _preparation.groupby("preparation_id").apply(
        lambda d: d.fillna(method="bfill").fillna(method="ffill")
    )

    n_unique_entries = _preparation.drop_duplicates().shape[0]
    n_unique_preparations = _preparation.preparation_id.unique().shape[0]
    assert n_unique_preparations == n_unique_entries, (
        n_unique_preparations,
        n_unique_entries,
    )

    print(f"Found {n_unique_preparations} unique preparations.")

    preparation = _preparation.drop_duplicates().set_index("preparation_id")
    assert preparation.index.is_unique
    assert preparation.stool_id.isin(stool.index).all()

    # Extract table: 'library'.

    column_rename = {
        "Project": "project_name",
        "External ID": "external_id",
        "data_type": "library_type",
        "# Lanes in Aggregation": "number_of_lanes_aggregated",
        "reads_raw": "sequenced_reads",
        #    'reads_filtered': 'sequenced_reads_filtered',
        #    'reads_qc_fail': 'sequenced_reads_qc_fail',
        #    'reads_human': 'sequenced_reads_human',
        #    'reads_ribosomal': 'sequenced_reads_ribosomal',
        #    'reads_viral': 'sequenced_reads_viral',
    }

    _library = meta[
        meta.data_type.isin(  # Exclude host_transcriptomics, methylome, host_genome
            [
                "viromics",
                "metabolomics",
                "metagenomics",
                "proteomics",
                "metagenomics",
                "metatranscriptomics",
                "stool_16S",
            ]
        )
    ][
        column_rename.keys()
    ].rename(  # Select columns of interest.
        columns=column_rename
    )  # Rename columns.

    # Remove suffixes _P, _TR.  What do these mean?
    _library["stool_id"] = _library.external_id.str.replace("_TR", "").str.replace(
        "_P", ""
    )
    # Construct the primary key
    _library["library_id"] = _library.external_id + "_" + _library.project_name
    # Construct the foreign key
    _library["preparation_id"] = _library.stool_id + "_" + _library.project_name

    _library.drop(columns=["stool_id", "project_name"], inplace=True)

    # Fill NaNs when other valid values.
    _library = _library.groupby("library_id").apply(
        lambda d: d.fillna(method="bfill").fillna(method="ffill")
    )

    n_unique_entries = _library.drop_duplicates().shape[0]
    n_unique_librarys = _library.library_id.unique().shape[0]
    assert n_unique_librarys == n_unique_entries, (n_unique_librarys, n_unique_entries)

    print(f"Found {n_unique_librarys} unique libraries.")

    library = _library.drop_duplicates().set_index("library_id")
    assert library.index.is_unique
    assert library.preparation_id.isin(preparation.index).all()

    mgen = library[library.library_type == "metagenomics"]
    mtab = library[library.library_type == "metabolomics"]

    for path, table in [
        (subject_path, subject),
        (visit_path, visit),
        (stool_path, stool),
        (preparation_path, preparation),
        (mgen_path, mgen),
        (mtab_path, mtab),
    ]:
        table.to_csv(path, sep="\t")
