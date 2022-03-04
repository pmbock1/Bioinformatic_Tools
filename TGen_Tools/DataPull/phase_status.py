#!/usr/bin/env python3
import json
import sys
import os
import csv


# USAGE: Used for large scale retrieval of phasing stasus from PharmCat JSON output
# Useful for cohort analysis of metabolism levels
# How to use: In command line : python3 ./phase_status.py /path/to/Directory_to_Search /path/to/filetowrite.csv

def phaser(direc, cfile):
    # sets prefix to input direc
    prefix_path = direc
    ext = '.report.json'  # IMPORTANT : change if different JSON is desired
    with open(cfile, mode='w') as datafile:
        datawrite = csv.writer(datafile, dialect='excel', delimiter=',',
                               quoting=csv.QUOTE_ALL)  # opens up CSV file to write
        datawrite.writerow(['sample', 'CACNA1S', 'CFTR', 'CYP2B6', 'CYP2C19'
                               , 'CYP2C9','CYP2D6', 'CYP3A5', 'CYPF2', 'DPYD',
                            'IFNL3', 'NUDT15', 'RYR1', 'SLCO1B1', 'TPMT', 'UGT1A1', 'VKORC1'])  # COLUMN HEADERS

        # WALKS DOWN (top to bottom) directory and subdirectories looking for files ending with EXTENSION designated
        for root, dirs, files in os.walk(direc, topdown=True):
            for file in files:
                dirs.sort()
                filename = root + os.sep + file
                if filename.endswith(ext):
                    with open(os.path.join(prefix_path, filename),
                              'r') as data:  # ensures script can reach file in directory
                        obj = json.loads(data.read())
                        fname = os.path.basename(os.path.normpath(file))
                        name = os.path.splitext(fname)[0]
                        print(name)
                        dictlist = [name]
                        for keys, values in obj.items():  # pulls out each dictionary in json
                            if keys == "geneCalls":  # only considers geneCalls field
                                for genes in values:
                                    for subfields, phase in genes.items():
                                        if phase == "G6PD" or phase =="MT-RNR1":
                                            break
                                        else:
                                            if subfields == "phaseStatus":  # Final dictionary that contains phaseStatus of each gene
                                                fixedoutput = str(phase)  # removes brackets and quotations
                                                dictlist.append(fixedoutput)
                        counttemp = 0
                    while counttemp < len(
                            dictlist):  # MAIN CSV writing. Change 'counttemp + COLUMN#' below to match number of data column (Genes + 1)
                        datawrite.writerow(dictlist[counttemp:counttemp + 17])
                        counttemp = counttemp + 17


if __name__ == "__main__":
    # calls function(first input, second input).
    phaser(sys.argv[1], sys.argv[2])
