"""Module to download Zenodo SDK xlsx files"""

import requests
import sys

def download_file(url, file_path):
    """Downloads a file from a URL"""

    response = requests.get(url)
    
    # Check if the request was successful (status code 200)
    if response.status_code == 200:
        with open(file_path, 'wb') as file:
            file.write(response.content)
        print(f"File downloaded successfully to {file_path}")
    else:
        print(f"Failed to download file. Status code: {response.status_code}")

if __name__ == "__main__":

    if "snakemake" in globals():
        download_file(snakemake.params.url, snakemake.output.xlsm)
    else:
        if len(sys.argv) != 3:
            msg = "Usage: python {} <url> <xlsm_file>"
            print(msg.format(sys.argv[0]))
            sys.exit(1)
        else:
            url = sys.argv[1]
            xlsm = sys.argv[2]
            download_file(url, xlsm)