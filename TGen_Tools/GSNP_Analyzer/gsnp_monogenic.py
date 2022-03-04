#!/usr/bin/env python
import pandas
import pandas as pd
import sys
import csv

import GSNP_Analyzer

pd.set_option('display.max_rows', None)

df_to_output = []
txt = 0

def monogenic(indir, inexcel, outfile):
    global df_to_output
    global txt

    # Determines input type using GSNP analyzer
    dataframes = GSNP_Analyzer.input_checker(indir)
    txt = GSNP_Analyzer.txt

    # Reads in excel, keeps only BMI information and ID
    excel_input = pandas.read_excel(inexcel)
    excel_input.drop(excel_input.columns.difference(['BMI (earliest value)', "BMI category", "Avatar ID", "Ancestry"]),
                     1, inplace=True)

    # Drops duplicates (multiple patients) and sets id as index
    excel_input.drop_duplicates(subset="Avatar ID", inplace=True)
    excel_input.set_index("Avatar ID", inplace=True)

    # Combines all CSV dataframes into one large dataframe
    combined_df = pd.concat(dataframes)
    combined_df = pd.DataFrame(combined_df)

    # Filters dataframe
    combined_df = GSNP_Analyzer.data_filtering(combined_df).reset_index()

    # Grabs count of each gene, adds to output to be written
    gene_list = combined_df["Gene (HGNC)"].value_counts()
    gene_list = pd.DataFrame(gene_list)
    gene_list.columns = ["Number of Rare/Private Variants"]

    # Table builder
    table = table_builder(combined_df)

    # Combined tables so BMI information is preset. Adds Yes/No column
    final_table = pd.concat([table, excel_input], axis=1)
    final_table["Variant in Genes"] = final_table.apply(lambda row: categorise(row), axis=1)

    # Adds homozygous variant column
    final_table["genotype"] = final_table["genotype"].astype(str)
    final_table["Homozygous Variant?"] = final_table.apply(lambda row: categorise2(row), axis=1)

    # Weighted Frequency and Ancestry/BMI information
    bmi_table = bmi_calc(final_table)
    ancestry_final = ancestry_calc(final_table)

    # Appendings dataframes to list
    df_to_output.append(ancestry_final)
    df_to_output.append(bmi_table)
    df_to_output.append(gene_list)

    # Reorders columns
    final_table = final_table[
        ["BMI category", "BMI (earliest value)", "Variant in Genes", "Homozygous Variant?", "Ancestry", "Gene (HGNC)",
         "dbSNP Identifier (142)",
         "biomarker", "genotype", "Tumor Type"]]

    df_to_output.append(final_table)

    with open(outfile, mode='w') as output:
        writer = csv.writer(output, dialect="excel")
        for df in df_to_output:
            df.to_csv(output, na_rep="")
            output.write("\n")
        df_to_output.clear()


def table_builder(dataframe):
    dataframe["dbSNP Identifier (142)"].fillna("Unknown(Private)", inplace=True)
    dataframe["dbSNP Identifier (142)"] = dataframe["dbSNP Identifier (142)"].astype(str)

    # Combines Variants into one row per patient
    counter = dataframe.groupby('ProjectRun').agg({"Gene (HGNC)": ', '.join, 'genotype': ', '.join,
                                                   'dbSNP Identifier (142)': ', '.join,
                                                   'biomarker': ', '.join}).reset_index()
    avatar_col = counter["ProjectRun"].str.split("_", n=3, expand=True)
    counter["Avatar ID"] = avatar_col[2]
    counter["Tumor Type"] = avatar_col[0]

    counter.drop("ProjectRun", axis=1, inplace=True)
    counter.set_index("Avatar ID", inplace=True)

    return counter


def categorise(row):
    if pd.isnull(row['Gene (HGNC)']) == False:
        return 'Yes'
    return 'No'


def categorise2(row):
    if "1/1" in row["genotype"]:
        return "Yes"
    return "No"

def ancestry_calc(dataframe):
    global txt
    categories = ["AFR", "EAS", "EUR", "Mix(EUR+AFR)","Mix(EUR+EAS)","Mix(EUR+NAT)","NAT","SAS"]
    homo_dictionary = {}
    hetero_dictionary = {}
    for cats in categories:
        homo_dictionary.update({cats : []})
        hetero_dictionary.update({cats : []})

    category_df = dataframe["Ancestry"].value_counts()

    for values in categories:
        homozygous_count = 0
        heterozygous_count = 0
        dataframe_g = dataframe[dataframe["Ancestry"] == values].reset_index()
        dataframe_g.drop(dataframe_g.columns.difference(['Ancestry', 'genotype']), axis=1, inplace=True)
        for index, row in dataframe_g.iterrows():
            genotype = row["genotype"]
            homozygous_count += genotype.count("1/1")
            heterozygous_count += genotype.count("0/1")

        homo_dictionary[values].append(homozygous_count)
        hetero_dictionary[values].append(heterozygous_count)

    homo_df = pd.DataFrame.from_dict(homo_dictionary, orient='index', columns=["number of homozygous variants"])
    hetero_df = pd.DataFrame.from_dict(hetero_dictionary, orient='index', columns=["number of heterozygous variants"])

    merged_df = pd.concat([category_df, homo_df, hetero_df], axis=1)
    merged_df.columns = ["total patients", "homozygous variants", "heterozygous variants"]

    merged_df["weighted average"] = (merged_df["homozygous variants"] + merged_df["heterozygous variants"]) / merged_df[
        "total patients"]
    merged_df["weighted average"] = round(merged_df["weighted average"], 3)
    return merged_df

def bmi_calc(dataframe):
    global txt
    categories = ["Normal weight", "Obese", "Overweight","Underweight"]
    homo_dictionary = {"Normal weight": [], "Obese": [], "Overweight": [],"Underweight": []}
    hetero_dictionary = {"Normal weight": [], "Obese": [], "Overweight": [], "Underweight": []}

    category_df = dataframe["BMI category"].value_counts()

    for values in categories:
        homozygous_count = 0
        heterozygous_count = 0
        dataframe_g = dataframe[dataframe["BMI category"] == values].reset_index()
        dataframe_g.drop(dataframe_g.columns.difference(['BMI category', 'genotype']), axis=1, inplace=True)
        for index, row in dataframe_g.iterrows():
            genotype = row["genotype"]
            homozygous_count += genotype.count("1/1")
            heterozygous_count += genotype.count("0/1")
        homo_dictionary[values].append(homozygous_count)
        hetero_dictionary[values].append(heterozygous_count)

    homo_df = pd.DataFrame.from_dict(homo_dictionary, orient='index', columns=["number of homozygous variants"])
    hetero_df = pd.DataFrame.from_dict(hetero_dictionary, orient='index', columns=["number of heterozygous variants"])

    merged_df = pd.concat([category_df,homo_df,hetero_df], axis=1)
    merged_df.columns = ["total patients","homozygous variants", "heterozygous variants"]

    merged_df["weighted average"] = (merged_df["homozygous variants"] + merged_df["heterozygous variants"]) / merged_df["total patients"]
    merged_df["weighted average"] = round(merged_df["weighted average"],3)

    merged_df.drop(index="not available",inplace=True)
    return merged_df

if __name__ == "__main__":
    monogenic(sys.argv[1], sys.argv[2], sys.argv[3])
