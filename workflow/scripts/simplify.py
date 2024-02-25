"""Removes default values from data file 

Code adapted from 
https://github.com/ClimateCompatibleGrowth/clicSANDMac/blob/main/cloud_converter/sand_filter_v2.py
"""

import os, sys
from collections import defaultdict
import re


def main(datafile_in, datafile_out):

    lines = []
    param_line = []
    parsing = False
    param_dr = False
    index_line = []
    year_line = []
    index_tag = []
    year_tag = []

    with open(datafile_in, 'r') as file_in:
        for line in file_in:
            line = line.strip('\t ').replace('\t', ' ')
            line = re.sub(' +', ' ', line)
            line = line.replace(' \n', '\n')

            if line.startswith('set YEAR'):
                start_year = line.split(':=')[1].split(' ')[1]

            if line.startswith('param'):
                parsing = True

                if line.startswith('param DiscountRate'):
                    param_dr = True

                if line.startswith('param ResultsPath'):
                    continue

            if parsing:
                line_values = []
                if line.startswith('param'):
                    line_elements = list(line.split(' '))
                    line_elements = [i.strip("\n:=") for i in line_elements]
                    lines.append(line)

                    if 'default' in line_elements:
                        default_index = line_elements.index('default')  # Find position of 'default'
                        default_value = line_elements[default_index+1]  # Extract default value

                elif line.startswith('['):
                    index_line = line
                    param_reset = True
                    index_tag = True

                elif line.startswith(str(start_year)):
                    year_line = []
                    year_line = line
                    param_reset = True
                    year_tag = True

                else:
                    line_values = list(set(line.rstrip('\n').split(' ')[1:-1]))

                    if param_dr:
                        if not line.startswith(';'):
                            lines.append(line)
                        param_dr = False

                    if line_values != list(default_value):
                        if line_values != []:
                            param_line = line

                    if len(param_line) > 0:
                        if param_reset:
                            if index_tag:
                                if index_line != []:
                                    lines.append(index_line)
                            if year_tag:
                                if year_line != []:
                                    lines.append(year_line)
                        param_reset = False
                        if not param_reset:
                            lines.append(param_line)
                            param_line = []
            if line.startswith(';'):
                if lines[-1].startswith("["):
                    lines.pop()
                # this is not a very robust solution :( 
                if ("default" in lines[-2].split(" ")) and (len(lines[-1].split(" ")) > 2):
                    lines.pop()
                if lines[-1].startswith("param"):
                    lines[-1] = lines[-1].replace(" :\n", " :=")
                lines.append(line)
                parsing = False
                index_tag = False
                year_tag = False
            elif not parsing:
                lines.append(line)

    with open(datafile_out, 'w') as file_out:
        file_out.writelines(lines)


if __name__ == "__main__":

    if "snakemake" in globals():
        main(snakemake.input.txt, snakemake.output.txt)
    else:
        if len(sys.argv) != 3:
            msg = "Usage: python {} <data_in> <data_out>"
            print(msg.format(sys.argv[0]))
            sys.exit(1)
        else:
            data_in=sys.argv[1]
            data_out = sys.argv[2]
            main(data_in, data_out)

