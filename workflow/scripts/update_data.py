"""Updates SDK data"""

from typing import Dict
import pandas as pd
import sys
from otoole import read, write

def correct_data(data: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
    """Corrects Discount Rate"""
    discount_rate = 0.05
    df = data["DiscountRate"]
    df["VALUE"] = discount_rate
    data["DiscountRate"] = df
    return data

if __name__ == "__main__":
    if "snakemake" in globals():
        config = snakemake.inputs.config
        data_in = snakemake.inputs.txt
        data_out = snakemake.outputs.txt
        data, defaults = read(config, 'datafile', data_in) 
        data = correct_data(data)
        write(config, 'datafile', data_out, data, defaults)
    else:
        if len(sys.argv) != 0: # replace 0 with the number of arguments you expect 
            msg = "Usage: python {}" # correct usage statement 
            print(msg.format(sys.argv[0]))
            sys.exit(1)
        else:
            # add parsing of arguments 
            data, defaults = read() # add arguments to read function 
            data = correct_data(data)
            write() # add arguments to write function