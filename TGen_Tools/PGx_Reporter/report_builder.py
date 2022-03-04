import pandas
import sys
import os
from datetime import date

import reporter


def builder(infile, outdir):
    #prepares name of output file
    remove_path = infile.split("/")[-1]
    name = remove_path.split(".")[0]
    out_html = "report_" + name + ".html"

    output = os.path.join(outdir,out_html)

    #Retrieves Date Information
    today = date.today()
    format_today = today.strftime("%B %d, %Y")

    #Calls reporter to make CSV files
    table_return = reporter.Reporter(infile,outdir)

    #Grabs all CSV tables
    table_list = []
    table_name_list = ["Metabolism Ability for Important Pharmacogenomics Genes"]

    table_list.append(outdir + "/metabolism.csv")

    for keys, values in table_return.items():
        table_list.append(keys)
        table_name_list.append(values[0])

    #turns CSV files into HTML tables
    html_tables = []
    for tables in table_list:
        temp = pandas.read_csv(tables)
        html_tables.append(temp.to_html(index=False, justify="left", na_rep="Not Available", table_id="tableStyle"))

    #HTML coding
    #header section. name, date, title
    html_str = "<header> <h1> PGx Report </h1> <p> " \
               "Name: " + name + "</p>" + \
               "<p> Date Created: " + format_today + "</p> </header>"

    #Information String Displayed at top of file
    information_str = """
    <span>
        Disclosures: This report is for RESEARCH PURPOSES only. All content is sourced from the CPIC database. This
        report was generated utilizing a variety of softwares to determine drug metabolism ability based off genomic data.
    </span>
    """
    # CSS Style - Internal
    html_str += """
    <style>
    h1 {
     border-radius: 25px;
     border: 2px solid #A7B0B9;
     padding: 20px;
     }
    span{
     display:inline-block;
     border: 2px solid #A7B0B9;
     width: 100%;
     margin-bottom: 30px;
     margin-top: 15px;
     }
    #tableStyle{
      margin-bottom: 30px;
      margin-right: 30px;
      }
    h2 {
     border-radius: 15px;
     border: 2px solid #A7B0B9;
     padding: 10px;
     font-size: 20px;
     width: 50%;
     }
    </style>
    """

    with open(output,"w") as out:
        out.write(html_str)
        out.write(information_str)

        for names, tables in zip(table_name_list,html_tables):
            out.write("<h2>" + names + "</h2>")
            out.write(tables)

if __name__ == "__main__":
    builder(sys.argv[1], sys.argv[2])