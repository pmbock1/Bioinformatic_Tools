#!/usr/bin/env python3
import json
import sys
import os
import csv

import analysis

def masterbuilder(direc, outdirec):

    # sets prefix to input direc
    prefix_path = direc
    ext = '.report.json'  # IMPORTANT : change if different JSON is desired

    # lists to hold sample name and subsequent retrieved information
    phase_list = []

    genotype_list = []
    genotype_list_final = []

    phenotype_list = []
    phenotype_list_final = []

    #Count used for CYP2D6 fix
    count = 6

    #walks down directory, looking for files with extension
    for root, dirs, files in os.walk(direc, topdown=True):
        for file in files:
            dirs.sort()
            filename = root + os.sep + file #ensures filename is retrieved properly

            #checks to see if extension is present
            if filename.endswith(ext):
                #opens file to read
                with open(os.path.join(prefix_path, filename),'r') as data:

                    #loads json files
                    obj = json.loads(data.read())
                    fname = os.path.basename(os.path.normpath(file))
                    name = os.path.splitext(fname)[0]

                    #prints name of file as it goes
                    print(name)

                    #lists to hold sample name and subsequent retrieved information
                    phase_list.append(name)
                    genotype_list = [name]
                    phenotype_list = [name]

                    for keys, values in obj.items():  # pulls out each dictionary in json

                        # Pulls Phasing Information
                        if keys == "geneCalls":  # only considers geneCalls field
                            for genes in values:
                                for subfields, phase in genes.items():
                                    if phase == "G6PD" or phase == "MT-RNR1":
                                        break
                                    else:
                                        if subfields == "phaseStatus":  # Final dictionary that contains phaseStatus of each gene
                                            output = str(phase)  # removes brackets and quotations
                                            phase_list.append(output)

                        # Pulls genotype/phenotype information
                        if keys == "genotypes":  # only considers genotype field
                            #Resets gene_holder list
                            gene_holder = []
                            for genes in values:
                                for subfields, genotypes in genes.items():
                                    if subfields == "phenotype":  # Final dictionary that contains phenotype of each gene
                                        if "," in str(genotypes):
                                            stringtosplit = str(genotypes)
                                            splitstring = stringtosplit.split(",", 1)[0]
                                            fixedoutput = splitstring[2:-1]
                                            phenotype_list.append(fixedoutput)
                                        else:
                                            fixedoutput = str(genotypes)[2:-2]  # removes brackets and quotations
                                            phenotype_list.append(fixedoutput)

                                    if subfields == "calls":  # Final dictionary that contains genotype of each gene
                                            fixedoutput = str(genotypes)[2:-2]  # removes brackets and quotations
                                            genotype_list.append(fixedoutput)
                                    if subfields == "gene":
                                            fixedoutput = str(genotypes)  # removes brackets and quotations
                                            gene_holder.append(fixedoutput)

                            #Checks if CYP2D6 call is absent and fixes output for consistency
                            if "CYP2D6" not in gene_holder:
                                genotype_list.insert(count,"Not Available")
                                phenotype_list.insert(count, "N/A")

                            #Stores "fixed" data into final list for table creation
                            genotype_list_final.extend(genotype_list)
                            phenotype_list_final.extend(phenotype_list)




    #Data Combining Section
    pheno_phase_list = []

    #Adds phasing status to genotype for phenotype-phase table
    for (phenotypes,phase) in zip(phenotype_list_final,phase_list):
        if phenotypes.endswith("report"):
            pheno_phase_list.append(phenotypes)
        else:
            pheno_phase_list.append(phenotypes + "(" + phase + ")")

    #sets up path for output
    out_path = outdirec

    #Set number of genes to +1 of genes (number of columns = 1 (sample) +  # of genes)
    num_genes = 17

    #Header row to be written for each table
    gene_list = ['sample', 'CACNA1S', 'CFTR', 'CYP2B6', 'CYP2C19', 'CYP2C9', 'CYP2D6','CYP3A5', 'CYPF2', 'DPYD','IFNL3', 'NUDT15', 'RYR1', 'SLCO1B1', 'TPMT', 'UGT1A1', 'VKORC1']

    # Makes directory if it does not exist
    if not os.path.exists(out_path):
        os.makedirs(out_path)

    # Makes directories for organization
    os.makedirs(out_path + "/analysis")
    os.makedirs(out_path + "/tables")

    #sets folders
    analysis_folder = out_path + "/analysis"
    tables_folder = out_path + "/tables"

    # Sets up file output names and paths
    filename_geno = os.path.join(tables_folder, 'genotype.csv')
    filename_pheno = os.path.join(tables_folder, 'phenotype.csv')
    filename_phase = os.path.join(tables_folder, 'phase_status.csv')
    filename_pheno_phase = os.path.join(tables_folder, 'phenotype_with_phase.csv')

    filename_analysisp = os.path.join(analysis_folder, 'analysis_phenotype.csv')
    filename_analysisg = os.path.join(analysis_folder, 'analysis_genotype.csv')
    filename_analysis_phase = os.path.join(analysis_folder, 'analysis_phased.csv')

    #Writes genotype table
    with open(filename_geno, mode='w') as datafile:
        datawrite = csv.writer(datafile, dialect='excel', delimiter=',',
                               quoting=csv.QUOTE_ALL)  # opens up CSV file to write
        datawrite.writerow(gene_list)  # COLUMN HEADERS

        counttemp = 0
        while counttemp < len(
                genotype_list_final):  # MAIN CSV writing. Change 'counttemp + COLUMN#' below to match number of data column (Genes + 1)
            datawrite.writerow(genotype_list_final[counttemp:counttemp + num_genes])
            counttemp = counttemp + num_genes


    #Writes phenotype table
    with open(filename_pheno, mode='w') as datafile:
        datawrite = csv.writer(datafile, dialect='excel', delimiter=',',
                               quoting=csv.QUOTE_ALL)  # opens up CSV file to write
        datawrite.writerow(gene_list)  # COLUMN HEADERS

        counttemp = 0
        while counttemp < len(
                phenotype_list_final):  # MAIN CSV writing. Change 'counttemp + COLUMN#' below to match number of data column (Genes + 1)
            datawrite.writerow(phenotype_list_final[counttemp:counttemp + num_genes])
            counttemp = counttemp + num_genes

    # Writes phase status table
    with open(filename_phase, mode='w') as datafile:
        datawrite = csv.writer(datafile, dialect='excel', delimiter=',',
                               quoting=csv.QUOTE_ALL)  # opens up CSV file to write
        datawrite.writerow(gene_list)  # COLUMN HEADERS

        counttemp = 0
        while counttemp < len(
                phase_list):  # MAIN CSV writing. Change 'counttemp + COLUMN#' below to match number of data column (Genes + 1)
            datawrite.writerow(phase_list[counttemp:counttemp + num_genes])
            counttemp = counttemp + num_genes


    #writes phenotype table that includes phase status
    with open(filename_pheno_phase, mode='w') as datafile:
        datawrite = csv.writer(datafile, dialect='excel', delimiter=',',
                               quoting=csv.QUOTE_ALL)  # opens up CSV file to write
        datawrite.writerow(gene_list)  # COLUMN HEADERS

        counttemp = 0
        while counttemp < len(
                pheno_phase_list):  # MAIN CSV writing. Change 'counttemp + COLUMN#' below to match number of data column (Genes + 1)
            datawrite.writerow(pheno_phase_list[counttemp:counttemp + num_genes])
            counttemp = counttemp + num_genes


    #Calls analysis script. Outputs in designated directory
    analyzer = analysis.analysis
    analyzer(filename_geno, filename_analysisg)
    analyzer(filename_pheno, filename_analysisp)
    analyzer(filename_phase, filename_analysis_phase)


if __name__ == "__main__":
    # calls function(first input, second input).
    masterbuilder(sys.argv[1], sys.argv[2])