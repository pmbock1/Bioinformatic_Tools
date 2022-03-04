#!/usr/bin/env python
import pandas as pd
import sys
import os
from itertools import combinations

import GSNP_Analyzer

# Lists of created dataframes
df_list = []


# main function
def merger(indir1, indir2, outdir):  # outcsv):
    # Determines if benign variants should be included
    benign_check = input("Do you want to include clinVar Benign Variants? (Y/N)")

    GSNP_Analyzer.benign_check = benign_check
    # Asks user for how many patients total were looked at to determine frequency
    txt = input("How many total patients were analyzed (first directory) ?: ")

    # Columns to include in analysis
    columns_to_include = ["ProjectRun", "Maximum Population Allele Freq", "genotype", "clinvar Significance",
                          "dbSNP Identifier (142)", "filetype", "frequencyCategory", "Functional Impact", "Gene (HGNC)"]

    # Empty lists where data will be stored
    CSV_dfs = []
    combs_dir1 = []
    combs_dir2 = []
    files_list = []
    avg_list_1 = []
    avg_list_2 = []

    # Reads in CSV from first directory as DataFrame
    for root, dirs, files in os.walk(indir1):
        for csvs in sorted(files):
            # Joins path so file can be accessed
            files_to_read = os.path.join(root, csvs)
            files_list.append(files_to_read)

        # Creates combinations of available CSVs of all lengths
        for i in range(1, len(files_list) + 1):
            els = combinations(files_list, i)
            combs_dir1.extend(els)

    # Checks all combination of files available for greatest difference in average variants
    for combination_of_files in combs_dir1:
        for filez in combination_of_files:
            datafile = pd.read_csv(filez)
            CSV_dfs.append(datafile)

        # Calls GSNP Analzyer to filter data
        GSNP_Analyzer.txt = txt
        GSNP_Analyzer.benign_check = benign_check

        combined_df = pd.concat(CSV_dfs)
        combined_df = GSNP_Analyzer.data_filtering(combined_df)

        print("Calculating Combinations for Directory 1")

        for (dcolumns, data) in combined_df.iteritems():
            # Only includes columns noted in "columns_to_include" list
            if dcolumns in columns_to_include:
                if dcolumns == "ProjectRun":
                    # Counts unique patients in file, if there are more than 1 it means there is multiple SNPs affecting patient
                    counter_value = combined_df[dcolumns].value_counts()

                    # Counts how many patients there are with variants
                    patient_count = counter_value.count()

                    # Turns counter value into a DataFrame, counts up occurences of 1 variant patients, 2 variant patients, etc
                    df_val = pd.DataFrame(counter_value.value_counts())

                    # Determines patients with 0 variants and adds to DataFrame
                    zero_variants = int(txt) - patient_count
                    df_val.loc[0] = zero_variants
                    df_val.rename_axis("Number of Consequential Variants", inplace=True)
                    df_val.columns = ["Number of People"]

                    # Calculates weighted average
                    num_homo = GSNP_Analyzer.weighted_avg_calculator(combined_df)

                    # Determines average of variants per patient
                    avg = df_val.index * df_val["Number of People"]
                    avg_variant = round((avg.sum() + num_homo) / int(txt), 3)
                    print("Dir1 avg: " + str(avg_variant))
                    avg_list_1.append(avg_variant)
                    # Resets dataframe list after every combination's average is calculated
                    CSV_dfs.clear()

    # Clears lists so second directory can be analyzed
    files_list.clear()

    txt2 = input("How many total patients were analyzed (second directory) ?: ")
    # Reads in CSV from second directory as DataFrame
    for root, dirs, files in os.walk(indir2):
        for csvs in sorted(files):
            # Joins path so file can be accessed
            files_to_read = os.path.join(root, csvs)
            files_list.append(files_to_read)
        # Creates combinations of available CSVs of all lengths
        for i in range(1, len(files_list) + 1):
            els = combinations(files_list, i)
            combs_dir2.extend(els)

    # Checks all combination of files available for greatest difference in average variants
    for combination_of_files in combs_dir2:
        for filez in combination_of_files:
            datafile = pd.read_csv(filez)
            CSV_dfs.append(datafile)

        # Calls GSNP Analzyer to filter data
        GSNP_Analyzer.txt = txt2
        GSNP_Analyzer.benign_check = benign_check

        combined_df = pd.concat(CSV_dfs)
        combined_df = GSNP_Analyzer.data_filtering(combined_df)

        print("Calculating Combinations for Directory 2")

        for (dcolumns, data) in combined_df.iteritems():
            # Only includes columns noted in "columns_to_include" list
            if dcolumns in columns_to_include:
                if dcolumns == "ProjectRun":
                    # Counts unique patients in file, if there are more than 1 it means there is multiple SNPs affecting patient
                    counter_value = combined_df[dcolumns].value_counts()

                    # Counts how many patients there are with variants
                    patient_count = counter_value.count()

                    # Turns counter value into a DataFrame, counts up occurences of 1 variant patients, 2 variant patients, etc
                    df_val = pd.DataFrame(counter_value.value_counts())

                    # Determines patients with 0 variants and adds to DataFrame
                    zero_variants = int(txt) - patient_count
                    df_val.loc[0] = zero_variants
                    df_val.rename_axis("Number of Consequential Variants", inplace=True)
                    df_val.columns = ["Number of People"]

                    # Calculates weighted average
                    num_homo = GSNP_Analyzer.weighted_avg_calculator(combined_df)

                    # Determines average of variants per patient
                    avg = df_val.index * df_val["Number of People"]
                    avg_variant = round((avg.sum() + num_homo) / int(txt2), 3)
                    print("Dir2 avg: " + str(avg_variant))
                    avg_list_2.append(avg_variant)
                    # Resets dataframe list after every combination's average is calculated
                    CSV_dfs.clear()

    # Zips averages into tuples so comparison can be made
    zipped_avg = zip(avg_list_1, avg_list_2)

    # Calculate difference between each tuple pair
    temp = [abs(b - a) for a, b in zipped_avg]
    print(temp)

    # Determines max difference in variant average and retrieves files
    max_value = max(temp)
    max_index = temp.index(max(temp))

    # Retrieves maxes for both directories
    dir1_files = combs_dir1[max_index]
    dir2_files = combs_dir2[max_index]

    # Prints out largest variant difference and the files responsible, asks user if they want script run on these files
    print('\n')
    print(
        "The largest average variant difference is " + str(round(max_value,3)) + " for files " + str(dir1_files) + " and " + str(dir2_files))
    checker = input("Would you like to run the gSNP analzyer script on these files? (Y/N) ")

    # If Y is selected
    if checker == "Y":
        # Runs GSNP_Analyzer on max valued files in dir1
        outfile = os.path.join(outdir + "dir1_analysis.csv")
        GSNP_Analyzer.txt = txt
        GSNP_Analyzer.merger(dir1_files, outfile)

        # Runs GSNP_Analyzer on max valued files in dir2
        outfile = os.path.join(outdir + "dir2_analysis.csv")
        GSNP_Analyzer.txt = txt2
        GSNP_Analyzer.merger(dir2_files, outfile)

        print("CSVs Created")

    # Other checker options
    elif checker == "N":
        quit()
    else:
        print("Invalid Option")


if __name__ == "__main__":
    merger(sys.argv[1], sys.argv[2], sys.argv[3])
