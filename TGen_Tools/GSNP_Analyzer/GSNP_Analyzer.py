#!/usr/bin/env python
import pandas as pd
import sys
import csv
import os

# Lists of created dataframes

df_list = []

# Initializes variables
txt = 0
benign_check = ""
avg_variant = 0
avg_weighted = 0

avg_weighted_common = 0
avg_weighted_rare = 0
avg_weighted_private = 0


# main function
def merger(indir, outfile):  # outcsv):
    # Declares global variables
    global txt
    global benign_check
    global avg_variant
    global avg_weighted
    global avg_weighted_common
    global avg_weighted_rare
    global avg_weighted_private

    # Initializes avg variant variable
    avg_variant = 0

    # Runs input checker to determine what kind of input was provided
    CSV_dfs = input_checker(indir)

    # Unique to MERGER script. Combines all CSV dataframes into one large dataframe
    combined_df = pd.concat(CSV_dfs)

    # Drops duplicate Data (same patient that has multiple germline samples that found the same SNP)
    combined_df = data_filtering(combined_df)

    # Calls Dataframe builder
    dataframe_builder(combined_df)

    # Row to be written at start of CSV
    fixed_average = round(avg_variant, 3)
    fixed_weighted = round(avg_weighted, 3)
    avg_row = ["Average Number of Consequential Variants Per Person: ", str(fixed_average), ""
        , "Average Number of Variants (Weighted):", str(fixed_weighted)]

    avg_common_row = ["Average Number of Common Variants (Weighted):", str(avg_weighted_common)]

    avg_rare_row = ["Average Number of Rare Variants (Weighted):", str(avg_weighted_rare)]

    avg_private_row = ["Average Number of Private Variants (Weighted):", str(avg_weighted_private)]

    weighted_category_list = [avg_common_row, avg_rare_row, avg_private_row]

    # CSV Writing
    with open(outfile, mode='w') as output:
        writer = csv.writer(output, dialect="excel")
        writer.writerow(avg_row)

        for averages in weighted_category_list:
            writer.writerow(averages)

        writer.writerow(("\n"))

        for df in df_list:
            df.to_csv(output, float_format='{:,.2%}'.format, na_rep="Unknown")
            output.write("\n")
        df_list.clear()


def input_checker(indir):
    global txt
    global benign_check

    CSV_dfs = []
    # Checks if input is a tuple pair (This means GVDF is using the script)
    if isinstance(indir, tuple):
        for files in indir:
            # Reads in CSV as a dataframe and adds to list
            datafile = pd.read_csv(files)
            CSV_dfs.append(datafile)
        return CSV_dfs

    # Reads in CSV as DataFrame, if input is a directory
    elif os.path.isfile(indir) == False:
        txt = input("Input is detected as a directory: How many total patients were analyzed?: ")
        benign_check = input("Do you want to include clinVar Benign Variants? (Y/N)")
        for root, dirs, files in os.walk(indir):
            for csvs in files:
                # Joins path so file can be accessed
                files_to_read = os.path.join(root, csvs)
                # Reads in CSV as a dataframe and adds to list
                datafile = pd.read_csv(files_to_read)
                CSV_dfs.append(datafile)
        return CSV_dfs

    # Reads in CSV as DataFrame, if input is a file
    elif os.path.isfile(indir) == True:
        txt = input("Input is detected as a file: How many total patients were analyzed?: ")
        benign_check = input("Do you want to include clinVar Benign Variants? (Y/N)")
        datafile = pd.read_csv(indir)
        CSV_dfs.append(datafile)
        return CSV_dfs


def data_filtering(dataframe):
    global txt
    global benign_check

    # Drop repeat header lines caused by multiple pages
    dataframe.drop(index=dataframe[dataframe['ProjectRun'] == 'ProjectRun'].index, inplace=True)

    # Drops Tumor samples
    new_dataframe = dataframe[~dataframe.filename.str.contains('Whole_T')]

    # Drop samples that have the same ProjectRun (meaning multiple germline samples in same patient)
    new_dataframe.drop_duplicates(subset=['ProjectRun', 'biomarker'], inplace=True)

    # Drop samples that have the same sample name (dropping patients part of multiple tumor groups)
    new_dataframe.drop_duplicates(subset=['sample name', 'biomarker'], inplace=True)

    # Drops Benign and likely benign clinVar variants if wanted
    new_dataframe["clinvar Significance"].fillna("Unknown", inplace=True)
    new_dataframe.reset_index(inplace=True)

    if benign_check == "N":
        print("dropping benigns...")
        new_dataframe.drop(new_dataframe[new_dataframe['clinvar Significance'] == "Benign"].index, inplace=True)
        new_dataframe.drop(new_dataframe[new_dataframe['clinvar Significance'] == "Benign/Likely_benign"].index,
                           inplace=True)

    return new_dataframe


def dataframe_builder(combined_df):
    global avg_variant
    global txt
    global avg_weighted

    # Calls specific analysis functions
    variant_table = variant_frequency_table(combined_df)
    gene_table = gene_information(combined_df)
    snp_table = snp_information(combined_df)

    # Adds tables to returned df_list
    df_list.append(variant_table)
    df_list.append(snp_table)
    df_list.append(gene_table)

    # Columns to run normal analysis on
    default_analysis_columns = ["genotype", "clinvar Significance", "filetype", "frequencyCategory",
                                "Functional Impact"]

    for dcolumns in default_analysis_columns:
        # Basic Analysis for all other columns. Normalizes frequency and combines frequency + number of occurence into dataframe
        counter = pd.concat([combined_df[dcolumns].value_counts(normalize=True, dropna=False),
                             combined_df[dcolumns].value_counts(dropna=False)], axis=1)
        counter.rename_axis(dcolumns, inplace=True)
        counter.columns = ["Percent Frequency", "Number of Variants"]
        df_list.append(counter)


def variant_frequency_table(combined_df):
    global avg_variant
    global txt
    global avg_weighted
    snp_avg_by_category(combined_df)

    # Counts unique patients in file, if there are more than 1 it means there is multiple SNPs affecting patient
    counter_value = combined_df["ProjectRun"].value_counts()

    # Counts how many patients there are with variants
    patient_count = counter_value.count()

    # Turns counter value into a DataFrame, counts up occurences of 1 variant patients, 2 variant patients, etc
    df_val = pd.DataFrame(counter_value.value_counts())

    # Determines patients with 0 variants and adds to DataFrame
    zero_variants = int(txt) - patient_count
    df_val.loc[0] = zero_variants
    df_val.rename_axis("Number of Consequential Variants", inplace=True)
    df_val.columns = ["Number of People"]

    # Determines average of variants per patient
    avg = df_val.index * df_val["Number of People"]
    avg_variant = (avg.sum() / int(txt))

    # Calculates weighted average
    num_homo = weighted_avg_calculator(combined_df)
    avg_weighted = (avg.sum() + num_homo) / int(txt)

    return df_val


def snp_information(combined_df):
    # Creates dataframe by grouping the dbSNP and maximim population frequency, adds number of variants
    counter = pd.DataFrame({'Number of Variants': combined_df.groupby(
        ["dbSNP Identifier (142)", "Maximum Population Allele Freq", "biomarker"],
        dropna=False).size()}).reset_index()
    # Names Columns
    counter.columns = ["RSID", "Maximum Population Frequency", "biomarker", "Number of Variants"]
    # Converts Population Frequency into float data type
    counter["Maximum Population Frequency"] = counter["Maximum Population Frequency"].astype("float64")
    # Orders dataframe by number of variants
    counter.sort_values("Number of Variants", ascending=False, inplace=True)

    # Sets and reorders the columns
    order_of_columns = ["biomarker", "RSID", "Number of Variants",
                        "Maximum Population Frequency"]
    counter = counter.reindex(columns=order_of_columns)

    # Sets index as biomarker so concat can take placex
    counter.set_index("biomarker", inplace=True)

    # Retrieves frequency and concats dataframe
    frequency = allele_frequency_calculator(combined_df)
    final_counter = pd.concat([counter, frequency], axis=1).reset_index()
    final_counter.columns = ["Biomarker", "RSID", "Number of People with Alt SNP (at-least 1)",
                             "Maximum Population Frequency", "Frequency of Alt SNP in Cohort"]

    # Final Ordering of dataframe
    final_counter.reindex(columns=order_of_columns)
    final_counter.set_index("Biomarker", inplace=True)

    return final_counter


def gene_information(combined_df):
    # Determines the number of unique variants per gene found in cohort
    # Groups each gene,biomarker combination
    counter = pd.DataFrame({'Number of Variants': combined_df.groupby(
        ["Gene (HGNC)", "biomarker"], dropna=False).size()}).reset_index()

    # Counts how many times each gene is found
    counter_value = counter["Gene (HGNC)"].value_counts().reset_index()

    # Sets and renames index to be the list of genes
    counter_value.set_index("index", inplace=True)

    # Creates DataFrame of how many genes have variants, combines with other dataframe
    other_info = pd.concat([combined_df["Gene (HGNC)"].value_counts(dropna=False),
                            combined_df["Gene (HGNC)"].value_counts(normalize=True, dropna=False)], axis=1)
    final_df = pd.concat([counter_value, other_info], axis=1)

    # Renames column and adds DataFrame to list
    final_df.rename_axis("Gene", inplace=True)
    final_df.columns = ["Number of Unique Variants", "Total Number of Variants",
                        "Percent of Total Variants"]
    return final_df


def allele_frequency_calculator(dataframe):
    # Hidden dataframe for frequency calculation
    freq_column = []

    # Counts the number of occurences of each SNP + Genotype combination
    frequency_df = pd.DataFrame({'Number of Variants': dataframe.groupby(
        ["biomarker", "genotype"], dropna=False).size()}).reset_index()

    # Iterates through df and calculates # of alleles based off genotype (2 for homozygous)
    for index, row in frequency_df.iterrows():
        if row[1] == "0/1":
            number_of_alleles = row[2] * 1
            freq_column.append(number_of_alleles)
        if row[1] == "1/1":
            number_of_alleles = row[2] * 2
            freq_column.append(number_of_alleles)
        if row[1] == "1/2":
            number_of_alleles = row[2] * 2
            freq_column.append(number_of_alleles)

    # Adds number of alleles as column
    frequency_df["Number of Alt Alleles"] = freq_column

    # Combines identical SNPs and names column
    freq_df = frequency_df.groupby("biomarker").agg({"Number of Alt Alleles": "sum"})
    freq_df.columns = ["Frequency of Alt SNP in Cohort"]

    # Calculates frequency by dividing # of each alleles by 2 * patients
    frequency = freq_df["Frequency of Alt SNP in Cohort"] / (int(txt) * 2)
    return frequency


def snp_avg_by_category(dataframe):
    global txt
    global avg_weighted_common
    global avg_weighted_rare
    global avg_weighted_private
    categories = ["common", "rare", "private"]
    dictionary = {"common": [], "rare": [], "private": []}

    category_df = dataframe["frequencyCategory"].value_counts()

    for values in categories:
        dataframe_g = dataframe[dataframe["frequencyCategory"] == values].reset_index()
        dataframe_g.drop(dataframe_g.columns.difference(['frequencyCategory', 'genotype']), axis=1, inplace=True)
        dataframe_g.drop(dataframe_g[dataframe_g["genotype"] != "1/1"].index, inplace=True)
        dictionary[values].append(len(dataframe_g))

    homo_df = pd.DataFrame.from_dict(dictionary, orient='index', columns=["number of homozygous variants"])
    merged_df = pd.concat([homo_df, category_df], axis=1)
    merged_df.columns = ["homozygous variants", "total variants"]

    merged_df["weighted average"] = (merged_df["homozygous variants"] + merged_df["total variants"]) / int(txt)

    avg_weighted_common = round(merged_df.at["common", "weighted average"], 3)
    avg_weighted_rare = round(merged_df.at["rare", "weighted average"], 3)
    avg_weighted_private = round(merged_df.at["private", "weighted average"], 3)


def weighted_avg_calculator(dataframe):
    # Counts how many homozygous genotypes we have
    try:
        genotype_counter = dataframe["genotype"].value_counts()
        homozygous = genotype_counter.loc["1/1"]
        return homozygous
    except:
        print("No homozygous variants found")
        return 0


def gene_SNP_combiner(indir):
    global benign_check
    CSV_dfs = []
    # Checks if input is a tuple pair (This means GVDF is using the script)
    if isinstance(indir, tuple):
        for files in indir:
            # Reads in CSV as a dataframe and adds to list
            datafile = pd.read_csv(files)
            CSV_dfs.append(datafile)

    # Reads in CSV as DataFrame, if input is a directory
    elif os.path.isfile(indir) == False:
        for root, dirs, files in os.walk(indir):
            for csvs in files:
                # Joins path so file can be accessed
                files_to_read = os.path.join(root, csvs)
                # Reads in CSV as a dataframe and adds to list
                datafile = pd.read_csv(files_to_read)
                CSV_dfs.append(datafile)
    # Reads in CSV as DataFrame, if input is a file
    elif os.path.isfile(indir) == True:
        datafile = pd.read_csv(indir)
        CSV_dfs.append(datafile)

    # Combined filtered dataframes
    combined_df = pd.concat(CSV_dfs)

    # Drops duplicate Data (same patient that has multiple germline samples that found the same SNP)
    combined_df = data_filtering(combined_df)

    counter_value = combined_df["ProjectRun"].value_counts()

    return counter_value


if __name__ == "__main__":
    merger(sys.argv[1], sys.argv[2])
