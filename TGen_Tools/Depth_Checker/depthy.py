import pandas
import sys
import os
import numpy as np

def depth(input,outfile):
    dataframes = []
    df_list = []
    for root, dirs, files in os.walk(input):
        for csvs in files:
            # Joins path so file can be accessed
            files_to_read = os.path.join(root, csvs)
            # Reads in CSV as a dataframe and adds to list
            datafile = pandas.read_table(files_to_read)
            dataframes.append(datafile)

    combined = pandas.concat(dataframes,axis=1)
    statdf =stats(combined)

    #num_positions_atleast_one = statdf["Individuals with at-least 10 Depth"][statdf["Individuals with at-least 10 Depth"] > 0].count()
    #print(num_positions_atleast_one)

    df_list.append(statdf)
    with open(outfile, mode='w') as output:
        for df in df_list:
            df.to_csv(output, na_rep="Unknown")
            output.write("\n")
        df_list.clear()


def stats(df):

    stats = round(df.mean(axis=1),2)
    statsdf = pandas.DataFrame(stats)
    statsdf.set_axis(['Average Depth of Position'], axis=1, inplace=True)

    statsdf["max"] = df.max(axis=1)
    statsdf["min"] = df.min(axis=1)
    statsdf["Depth Range"] = statsdf["min"].astype(str) + " - " + statsdf["max"].astype(str)
    statsdf.drop(["max","min"],axis=1,inplace=True)

    statsdf["Individuals with at-least 10 Depth"] = df.select_dtypes(np.number).gt(9).sum(axis=1)
    statsdf["Percent of Individuals with at-least 10 Depth"] = round(statsdf["Individuals with at-least 10 Depth"] / len(df.columns),2) * 100
    statsdf["Percent of Individuals with at-least 10 Depth"] = statsdf["Percent of Individuals with at-least 10 Depth"].astype(str) + "%"
    statsdf = statsdf.reindex(columns=["Individuals with at-least 10 Depth","Percent of Individuals with at-least 10 Depth",
                             'Average Depth of Position',"Depth Range"])

    return statsdf



if __name__ == "__main__":
    depth(sys.argv[1],sys.argv[2])