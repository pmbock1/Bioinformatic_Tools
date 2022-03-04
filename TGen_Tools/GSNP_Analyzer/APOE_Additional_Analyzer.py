#!/usr/bin/env python
import pandas as pd
import sys

import APOE
import GSNP_Analyzer

APOE_wt = 0


def SNP_combiner(APOE_csv, other_csv, output):
    global APOE_wt
    # Asks for benign variants and initializes list of df
    benign_check = input("Do you want to include clinVar Benign Variants? (Y/N)")
    df_list = []

    # Sets benign check equal to the GSNP benign check so filtering can occur
    GSNP_Analyzer.benign_check = benign_check

    # Call GSNP analyzer for secondary analysis, APOE for main
    other_df = GSNP_Analyzer.gene_SNP_combiner(other_csv)
    APOE_list = APOE.APOE_Additional(APOE_csv)

    # Seperates list into parts, first is dataframe second is number of APOE3/APOE2
    APOE_df = APOE_list[0]
    APOE_wt = APOE_list[1]

    # Combines dataframes, all patients with either variant are now present
    final_df = pd.concat([other_df, APOE_df], axis=1)
    final_df.columns = ["Number of Variants", "Genotype", "SNP", "Biomarker", "Diplotype"]
    final_df.rename_axis("Patient", inplace=True)

    # Fills NAs for each column caused by concat, only Number of Variants is important
    final_df["Number of Variants"].fillna("0", inplace=True)
    final_df["Genotype"].fillna("wild-type", inplace=True)
    final_df["SNP"].fillna("None", inplace=True)
    final_df["Biomarker"].fillna("None", inplace=True)
    final_df["Diplotype"].fillna("APOE3/APOE3", inplace=True)

    # calls grouping method
    grouped = grouper(final_df)

    # Adds dfs to list and writes to output csv
    df_list.append(grouped)
    df_list.append(final_df)

    with open(output, mode='w') as output:
        for df in df_list:
            df.to_csv(output, na_rep="Unknown")
            output.write("\n")
        df_list.clear()


def grouper(dataframe):
    # Resets index, drops extra columns for grouping improvement
    new_df = dataframe.set_index("Diplotype")
    new_df.drop(["Genotype", "SNP", "Biomarker"], axis=1, inplace=True)

    # Ensures variant column is float data type (not string)
    new_df["Number of Variants"] = new_df["Number of Variants"].astype("float64")

    # Groups by diplotype, turns into dataframe. This df sees how many of each diplotype there is
    saved_df = new_df.groupby(by="Diplotype").size()
    saved_df = pd.DataFrame(saved_df)
    saved_df.columns = ["Size"]
    # Since all APOE3/APOE3 wont be here (only those with variants) , updates the number to all APOE3 wild types
    saved_df.at["APOE3/APOE3", "Size"] = APOE_wt

    # Groups original df by diplotype and SUMS column (so we have # of variants for each genotype group)
    grouped_df = new_df.groupby(by="Diplotype").sum()
    grouped_df_t = pd.DataFrame(grouped_df)

    # Divides # of variants by saved_df, variants/patients = variant average per patient. also rounds column
    grouped_df_t["Number of Variants"] = grouped_df_t["Number of Variants"] / saved_df["Size"]
    grouped_df_t['Number of Variants'] = grouped_df_t['Number of Variants'].astype(float).round(2)

    return grouped_df_t


if __name__ == "__main__":
    SNP_combiner(sys.argv[1], sys.argv[2], sys.argv[3])
