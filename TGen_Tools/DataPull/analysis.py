#!/usr/bin/env python
import pandas as pd
import sys

# Lists of created dataframes
df_list = []


# main function
def analysis(incsv, outfile):  # outcsv):
    # Read column names from file
    cols = list(pd.read_csv(incsv, nrows=1))
    # Use list comprehension to remove the unwanted column in **usecol**
    datafile = pd.read_csv(incsv, usecols=[i for i in cols if i != 'sample'])
    # iterates through csv file
    for (dcolumns, data) in datafile.iteritems():
        counter = pd.concat([datafile[dcolumns].value_counts(normalize=True, dropna=False),
                             datafile[dcolumns].value_counts(dropna=False)], axis=1)
        counter.rename_axis(dcolumns, inplace=True)
        counter.columns = ["Percent Frequency", "Number of Patients"]
        df_list.append(counter)

    with open(outfile, mode='w') as output:
        for df in df_list:
            df.to_csv(output, float_format='{:.2%}'.format, na_rep="Unknown")
            output.write("\n")
        df_list.clear()


if __name__ == "__main__":
    analysis(sys.argv[1], sys.argv[2])
