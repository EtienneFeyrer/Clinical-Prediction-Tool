# download_hpo_files.py
import requests
import os

def download_hpo_files():
    """Download HPO files to local directory."""
    files = {
        'hp.obo': 'https://github.com/obophenotype/human-phenotype-ontology/releases/download/v2025-09-01/hp.obo',
        'phenotype_to_genes.txt': 'https://github.com/obophenotype/human-phenotype-ontology/releases/download/v2025-09-01/phenotype_to_genes.txt'
    }

    for filename, url in files.items():
        print(f"Downloading {filename}...")
        response = requests.get(url)
        with open(filename, 'wb') as f:
            f.write(response.content)
        print(f"Downloaded {filename}")

if __name__ == "__main__":
    download_hpo_files()