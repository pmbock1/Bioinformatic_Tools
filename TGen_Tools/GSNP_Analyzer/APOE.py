#!/usr/bin/env python
import pandas as pd
import sys
import csv

# Lists of created dataframes
df_list = []
APOE4_freq = 0
APOE3_freq = 0
APOE2_freq = 0
APOE_wt = 0
txt = ""

df_genotyped = []

# List of diplotypes
list_of_diplotypes = []
# Columns to include in analysis
columns_to_include = ["ProjectRun", "dbSNP Identifier (142)", "genotype", "biomarker"]


# main function
def merger(indir, outfile):  # outcsv):
    # sets variables as global
    global APOE2_freq
    global APOE3_freq
    global APOE4_freq
    global APOE_wt
    global txt

    # Asks user for how many patients total were looked at to determine frequency
    txt = input("How many total patients were analyzed?: ")

    # Empty lists where data will be stored
    CSV_dfs = []

    # Reads in CSV file
    datafile = pd.read_csv(indir)
    CSV_dfs.append(datafile)

    # Drops duplicate Data (same patient that has multiple germline samples that found the same SNP)
    data_filtering(CSV_dfs)

    # Combines all CSV dataframes into one large dataframe (should only be 1)
    combined_df = pd.concat(CSV_dfs)
    # Renames columns
    combined_df.columns = columns_to_include
    # Fixes SNP column to be string instead of float
    combined_df["dbSNP Identifier (142)"] = combined_df["dbSNP Identifier (142)"].astype(str)

    # Combines Variants into one row per patient
    counter = combined_df.groupby('ProjectRun').agg({'genotype': ', '.join,
                                                     'dbSNP Identifier (142)': ', '.join,
                                                     'biomarker': ', '.join}).reset_index()

    # Genotyper Function calling
    genotyper(counter)

    # Creates a new column "diplotypee" and assigns the diplotype list to it
    counter["diplotype"] = list_of_diplotypes

    # Calls reportbuilder to finish the analysis
    reportbuilder(counter, outfile)


def genotyper(dataframe):
    global APOE2_freq
    global APOE3_freq
    global APOE4_freq
    global APOE_wt
    global list_of_diplotypes
    # Initialized Variables
    APOE = []
    # Intiates Count of Haplotypes
    APOE4 = 0
    APOE3 = 0
    APOE2 = 0

    # APOE Diplotyping. Goes through every combination of variants defining APOE alleles, adds to list
    for index, row in dataframe.iterrows():
        if row[2] == "rs429358" and row[1] == "0/1":
            list_of_diplotypes.append("APOE3/APOE4")
            APOE4 += 1
            APOE3 += 1
        elif row[2] == "rs429358" and row[1] == "1/1":
            list_of_diplotypes.append("APOE4/APOE4")
            APOE4 += 2
        elif row[2] == "rs7412" and row[1] == "0/1":
            list_of_diplotypes.append("APOE3/APOE2")
            APOE3 += 1
            APOE2 += 1
        elif row[2] == "rs7412" and row[1] == "1/1":
            list_of_diplotypes.append("APOE2/APOE2")
            APOE2 += 2
        elif row[2] == "rs7412, rs429358" or row[2] == "rs429358, rs7412":
            list_of_diplotypes.append("APOE2/APOE4")
            APOE4 += 1
            APOE2 += 1
        # If row[2] is not in APOE allele format
        elif row[2] not in APOE:
            if row[1] == "0/1":
                list_of_diplotypes.append("APOE3/" + row[2])
                APOE3 += 1
            elif row[1] == "1/1":
                list_of_diplotypes.append(row[2] + "/" + row[2])
            else:
                list_of_diplotypes.append("Unknown")

    # Calculates how many people had wild-tyoe diplotype (APOE3/APOE3)
    APOE_wt = int(txt) - dataframe["ProjectRun"].count()
    APOE3 += (APOE_wt * 2)

    # Calculates and Stores APOE Allele Frequencies. Divides allele count by total alleles (2 * people)
    APOE4_freq = APOE4 / (int(txt) * 2)
    APOE3_freq = APOE3 / (int(txt) * 2)
    APOE2_freq = APOE2 / (int(txt) * 2)

    return APOE_wt


def reportbuilder(dataframe, outfile):
    # Makes new DataFrame to count each occurence of each diplotype
    diplotype_counter = pd.DataFrame(dataframe["diplotype"].value_counts()).reset_index()
    diplotype_counter.columns = ["Diplotype", "Number of People"]

    # Adds in the "normal" patients. People having no APOE variants are wildtype
    diplotype_counter.loc[len(diplotype_counter.index)] = ["APOE3/APOE3", APOE_wt]

    # Sorts by number of people and resets index as the diplotypes
    diplotype_counter.sort_values("Number of People", ascending=False, inplace=True)
    diplotype_counter.set_index("Diplotype", inplace=True)
    diplotype_counter["Frequency in Cohort"] = diplotype_counter["Number of People"] / int(txt)

    # Sets index of original dataframe as patient ID
    dataframe.set_index("ProjectRun", inplace=True)

    # Adds Dataframes to be written out to CSV file
    df_list.append(diplotype_counter)
    df_list.append(dataframe)

    # CSV Writing
    with open(outfile, mode='w') as output:
        # Creates writer object
        writer = csv.writer(output, dialect="excel")

        # Round Frequencies
        APOE4_freq_fixed = round(APOE4_freq, 3)
        APOE3_freq_fixed = round(APOE3_freq, 3)
        APOE2_freq_fixed = round(APOE2_freq, 3)

        # Writes frequency row
        writer.writerow(["APOE4 Frequency is: ", str(APOE4_freq_fixed), " ", "Worldwide Frequency:", "0.137"])
        writer.writerow(["APOE3 Frequency is: ", str(APOE3_freq_fixed), " ", "Worldwide Frequency:", "0.779"])
        writer.writerow(["APOE2 Frequency is: ", str(APOE2_freq_fixed), " ", "Worldwide Frequency:", "0.084"])

        writer.writerow(("\n"))

        # Writes out all dataframes in df_list
        for df in df_list:
            df.to_csv(output, float_format='{:,.2%}'.format, na_rep="Unknown")
            output.write("\n")
        df_list.clear()


def data_filtering(list_of_dataframes):
    global txt
    # Drops duplicate Data (same patient that has multiple germline samples that found the same SNP)
    for csvs in list_of_dataframes:
        # Drop samples that have the same ProjectRun (meaning multiple germline samples in same patient)
        csvs.drop_duplicates(subset=['ProjectRun', 'biomarker'], inplace=True)

        # Drop samples that have the same sample name (dropping patients part of multiple tumor groups)
        # Also updates the number of patients to accomadate (same sample counted as multiple patients)
        txt = int(txt) - (len(csvs['ProjectRun']) - len(csvs.drop_duplicates(subset=['sample name', 'biomarker'])))
        csvs.drop_duplicates(subset=['sample name', 'biomarker'], inplace=True)

        # Drop repeat header lines caused by multiple pages
        csvs.drop(index=csvs[csvs['ProjectRun'] == 'ProjectRun'].index, inplace=True)

        # Drops all other columns not wanted for analysis
        csvs.drop(csvs.columns.difference(columns_to_include), axis=1, inplace=True)
    print("The number of patients has been reduced to " + str(txt) + " after filtering")


def APOE_Additional(indir):
    # sets variables as global
    global APOE2_freq
    global APOE3_freq
    global APOE4_freq
    global txt

    # Asks user for how many patients total were looked at to determine frequency
    txt = input("How many total patients were analyzed?: ")

    # Empty lists where data will be stored
    CSV_dfs = []

    # Reads in CSV file
    datafile = pd.read_csv(indir)
    CSV_dfs.append(datafile)

    # Drops duplicate Data (same patient that has multiple germline samples that found the same SNP)
    data_filtering(CSV_dfs)

    # Combines all CSV dataframes into one large dataframe (should only be 1)
    combined_df = pd.concat(CSV_dfs)
    # Renames columns
    combined_df.columns = columns_to_include
    # Fixes SNP column to be string instead of float
    combined_df["dbSNP Identifier (142)"] = combined_df["dbSNP Identifier (142)"].astype(str)

    # Combines Variants into one row per patient
    counter = combined_df.groupby('ProjectRun').agg({'genotype': ', '.join,
                                                     'dbSNP Identifier (142)': ', '.join,
                                                     'biomarker': ', '.join}).reset_index()

    # Genotyper Function calling
    APOE_wt = genotyper(counter)

    # Creates a new column "diplotypee" and assigns the diplotype list to it
    counter["diplotype"] = list_of_diplotypes
    counter.set_index("ProjectRun", inplace=True)

    # Returns list containing dataframe and number APOE wild-type
    return_list = [counter, APOE_wt]
    return return_list


if __name__ == "__main__":
    merger(sys.argv[1], sys.argv[2])
