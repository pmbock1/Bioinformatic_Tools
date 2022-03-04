import csv
import sys
import os

def catergorizer(in_list,outdirectory):
    # Drug Categories
    nsaid = ["nsaid", "celexocib", "flurbiprofen", "ibuprofen", "lornoxicam", "meloxicam", "piroxicam", "tenoxicam"]
    anti = ["antidepressants", "amitriptyline", "clomipramine", "doxepin", "imipramine", "trimipramine", "desipramine",
            "nortriptyline", "citalopram", "escitalopram", "sertraline", "fluvoxamine", "paroxetine"]
    anes = ["anesthetics", "desflurane", "enflurane", "halothane", "isoflurane", "methoxyflurane", "sevoflurane",
            "succinylcholine"]
    opioids = ["opioid","codeine","hydrocodone","tramadol"]
    chemotherapy = ["chemotherapy","capecitabine","fluouracil","thioguanine","mercaptopurine"]

    master = [nsaid, anti, anes, opioids, chemotherapy]

    #Makes Category folder
    path = outdirectory
    pathfolder = path + "/drug_categories"

    if not os.path.exists(pathfolder):
        os.makedirs(pathfolder)

    #Empty Lists to fill drugs if present
    nsaid_present = ["NSAIDs"]
    anti_present = ["Antidepressants (tricylics and SSRIs)"]
    anes_present = ["General Anesthetics"]
    op_present = ["Opiate Analgesics"]
    chemo_present = ["Chemotherapy (antimetabolites)"]

    # Data Classifier: Filtered Data is separated into drug categories
    for indx, values in enumerate(in_list):
        for categories in master:
            for drugs in categories:
                if values == drugs and categories[0] == "nsaid":
                    nsaid_present.extend(in_list[indx:indx+5])
                if values == drugs and categories[0] == "antidepressants":
                    anti_present.extend(in_list[indx:indx+5])
                if values == drugs and categories[0] == "anesthetics":
                    anes_present.extend(in_list[indx:indx+5])
                if values == drugs and categories[0] == "opioids":
                    op_present.extend(in_list[indx:indx+5])
                if values == drugs and categories[0] == "chemotherapy":
                    chemo_present.extend(in_list[indx:indx + 5])

    #Category Writing
    antidepressant_table = os.path.join(pathfolder,"antidepressants.csv")
    nsaid_table = os.path.join(pathfolder, "nsaid.csv")
    anesthetics_table = os.path.join(pathfolder,"anesthetics.csv")
    opioids_table = os.path.join(pathfolder, "opioids.csv")
    chemotherapy_table = os.path.join(pathfolder, "chemotherapy.csv")

    # Catches drugs not included in categories and writes them to the "other" table
    other_table = os.path.join(pathfolder, "other.csv")
    other_present = ["Other"]

    #checks every 5 elements of list (a drug and its additional information) to see if its part of a category
    counter = 0
    for drugs in in_list[::5]:
        counter += 5
        #if not part of a category, add it to other list
        if drugs not in (item for catergories in master for item in catergories):
            other_present.extend(in_list[counter - 5:counter])

    #Dictionary defining each table to each drug category created based off patient data
    table_dictionary = {antidepressant_table: anti_present, nsaid_table: nsaid_present, anesthetics_table: anes_present
                        ,opioids_table: op_present, chemotherapy_table: chemo_present, other_table: other_present}

    #Creates list of tables present for patient
    tables_present = {}
    for keys, values in table_dictionary.items():
        if len(values) != 1:
            tables_present.update([(keys,values)])
            with open(keys,'w') as dfile:
                datawrite = csv.writer(dfile, dialect='excel', delimiter=',',
                                   quoting=csv.QUOTE_ALL)  # opens up CSV file to write
                datawrite.writerow(['Actionable Drug', 'Guideline', "Population", 'Recommendation Strength',
                                'Gene(s) Responsible'])  # COLUMN HEADERS

                counttemp = 1
                # Writes information regarding drug guidelines
                while counttemp < len(values):
                    datawrite.writerow(values[counttemp:counttemp + 5])
                    counttemp = counttemp + 5

    #Returns List of what tables are present
    #print(tables_present)
    return tables_present

# Calls script
if __name__ == "__main__":
    catergorizer(sys.argv[1],sys.argv[2])