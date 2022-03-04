import csv
import sys
import json
import os

import catergorizer


def Reporter(input, output):
    # Instance Variables (Data Values that will be retrieved and stored)
    listy = []
    phenotype = ""
    missingFlag = ""
    gene = ""
    reportable = ""
    recommendationClass = ""
    guideline = ""
    population = ""
    reporty = []
    dicy = []
    dicy_info = []
    drugName = ""
    geneholder = []

    #Sets up CSV Output Directory
    path = output
    if not os.path.exists(path):
        os.makedirs(path)

    guideline_csv = os.path.join(path, 'guideline.csv')
    metabolism_csv = os.path.join(path, 'metabolism.csv')

    # Generic Dats Outputs List (All of these are considered generic and not actionable, should not be included)
    Generic = [
        'Initiate therapy with recommended starting dose. In accordance with the prescribing information, use the lowest effective dosage for shortest duration consistent with individual patient treatment goals.',
        "Based on genotype, there is no indication to change dose or therapy. Use label-recommended dosage and administration.",
        "Clinical findings, family history, further genetic testing and other laboratory data should guide use of halogenated volatile anesthetics or depolarizing muscle relaxants.",
        "Clopidogrel: label-recommended dosage and administration",
        "There is no need to avoid prescribing of atazanavir based on UGT1A1 genetic test result",
        "Initiate efavirenz with standard dosing (600 mg/day)"]

    # DATA RETRIEVAL
    # JSON File Reading. Extracts Information for Report
    with open(input, "r") as infile:
        obj = json.loads(infile.read())
        for layer, extra in obj.items():
            # Collects information regarding gene function and phenotype
            if layer == "genotypes":
                for info in extra:
                    for key, value in info.items():
                        if key == "phenotype":
                            phenotype = value
                        if key == "gene":
                            gene = value
                        if key == "missingVariants":
                            missingFlag = value
                    dictInfo = dict({"Gene": gene, "Phenotype": phenotype, "Missing Variants": missingFlag})
                    dicy.append(dictInfo)

            # Retrieves information regarding guidelines and medication recommendations
            if layer == "guidelines":
                for info in extra:
                    for key, value in info.items():
                        if key == "geneCalls":
                            for calls in value:
                                for geneinfo, values in calls.items():
                                    if geneinfo == "gene":
                                        geneholder.append(values)
                            gene = geneholder
                            # print(gene)
                        if key == "drugs":
                            drugName = value
                            # print(drugName)
                        if key == "reportable":
                            reportable = value
                        if key == "groups":
                            for comments in value:
                                for name, guidelines in comments.items():
                                    for drugs in guidelines:
                                        for storage, recommendations in drugs.items():
                                            if recommendations == "Recommendation":
                                                guideline = list(drugs.items())[0][1]
                                            if recommendations == "Population":
                                                population = list(drugs.items())[0][1]
                                            if recommendations == "Classification of Recommendation":
                                                recommendationClass = list(drugs.items())[0][1]

                                    # BUILDS DICTIONARY from retrieved values
                                    dictDrugs = dict(
                                        {"Drug": drugName, "Guidelines": guideline, "Population": population,
                                         "Class of Recommendation": recommendationClass, "Gene(s) Responsible": gene})
                                    # print(dictDrugs)
                                    dictDrugWrap = dict({drugName: dictDrugs, "reportable": reportable})
                                    listy.append(dictDrugWrap)
                        # Resets gene holder list
                        geneholder = []

                # DATA FILTERING
                # List builder for drug guideline information (dictDrugWrap list)
                for items in listy:
                    # checks if the drug guideline is reportable
                    if True in items.values():
                        for drugs, info in items.items():
                            # ignores reportable value
                            if info != True:
                                for guidelines, values in info.items():
                                    reporty.append(values)
                # List builder for phenotype information (dictInfo list)
                for items in dicy:
                    for drugs, info in items.items():
                        dicy_info.append(info)

                # Removes "Normal" Dosage recommendations so only actionable drug changes are reported
                valcounter = len(reporty)
                for values in reversed(reporty):
                    valcounter -= 1
                    for messages in Generic:
                        if str(values).startswith(messages):
                            del reporty[valcounter-1:valcounter+4]

                # Removes brackets and quotes from data outputs (if present)
                valcounter2 = 0
                for values in reporty:
                    valcounter2 += 1
                    if isinstance(values, list):
                        fixedvalue = str(values)
                        translated = fixedvalue.replace("'", "")
                        reporty[valcounter2 - 1] = translated[1:-1]

                valcounter3 = 0
                for values in dicy_info:
                    valcounter3 += 1
                    if isinstance(values, list):
                        fixedvalue = str(values)
                        translated = fixedvalue.replace("'", "")
                        dicy_info[valcounter3 - 1] = translated[1:-1]

    #Calls Drug Categorizer script
    table_list = catergorizer.catergorizer(reporty,path)


    # CSV Writing
    #Guideline CSV
    with open(guideline_csv, "w") as outfile:
        datawrite = csv.writer(outfile, dialect='excel', delimiter=',',
                               quoting=csv.QUOTE_ALL)  # opens up CSV file to write
        datawrite.writerow(['Actionable Drug', 'Guideline', "Population", 'Recommendation Strength',
                            'Gene(s) Responsible'])  # COLUMN HEADERS

        counttemp = 0
        # Writes information regarding drug guidelines
        while counttemp < len(
                reporty):
            datawrite.writerow(reporty[counttemp:counttemp + 5])
            counttemp = counttemp + 5
    #Metabolism CSV
    with open(metabolism_csv, "w") as outfile:
        datawrite = csv.writer(outfile, dialect='excel', delimiter=',',
                                   quoting=csv.QUOTE_ALL)  # opens up CSV file to write
        datawrite.writerow(['Gene', "Metabolism", "Missing Data Input"])  # COLUMN HEADERS
        counttemp2 = 0
        # Writes information regarding gene metabolism
        while counttemp2 < len(
                dicy_info):
            datawrite.writerow(dicy_info[counttemp2:counttemp2 + 3])
            counttemp2 = counttemp2 + 3

    #Returns list of tables created by "categorizer" script
    return table_list

#Calls script
if __name__ == "__main__":
    Reporter(sys.argv[1], sys.argv[2])
