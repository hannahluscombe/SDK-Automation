"""Updates SDK data"""

from typing import Dict
import pandas as pd
import sys
from otoole import read, write

def correct_data(data: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Corrects Discount Rate"""
    # add logic here
    pass

if __name__ == "__main__":
    if "snakemake" in globals():
        data, defults = read() # add arguments to read function 
        data = correct_data(data)
        write() # add arguments to write function
    else:
        if len(sys.argv) != 0: # replace 0 with the number of arguments you expect 
            msg = "Usage: python {}" # correct usage statement 
            print(msg.format(sys.argv[0]))
            sys.exit(1)
        else:
            # add parsing of arguments 
            data, defults = read() # add arguments to read function 
            data = correct_data(data)
            write() # add arguments to write function