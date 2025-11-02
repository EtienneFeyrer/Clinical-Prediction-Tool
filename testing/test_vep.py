"""
    Annotates a VCF file using Ensembl VEP API and saves results in a JSON file.

    Args:
        input_vcf_path (str): Path to input VCF or VCF.GZ file.
        output_json_path (str): Path to output JSON file.
        fields (list): List of VEP fields to request.
    """


import gzip
import requests
import sys


def annotate_vcf_with_vep(input_vcf_path, output_json_path):

    if not output_json_path.lower().endswith(".json"):
        raise ValueError("Output file must have a .json extension")

    opener = gzip.open if input_vcf_path.endswith(('.gz', '.bgz')) else open
    variants = []
    with opener(input_vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            chrom, pos, ref, alt = cols[0], cols[1], cols[3], cols[4]
            for a in alt.split(","):
                variants.append(f"{chrom} {pos} {pos} {ref}/{a} +")

    payload = {
        "variants": variants,
        "REVEL" : True,
        "CADD" : True,
        "SpliceAI": True,
        "protein": True,
        "gencode_basic": True,
        "LoF": True,
        "mane": True,
        "hgvs": True,
        "dbNSFP": "clinvar_OMIM_id,GERP++_RS",
    }

    headers = {
        "Content-Type": "application/json",
        "Accept": "application/json"
    }
    vep_url = "https://rest.ensembl.org/vep/human/region"

    print("Sending request to VEP API...")
    response = requests.post(vep_url, headers=headers, json=payload)

    if response.ok:
        with open(output_json_path, "w") as out_file:
            out_file.write(response.text)
        print(f" Annotation complete. Results saved to {output_json_path}")
    else:
        print(f"API Error {response.status_code}: {response.text}")
        sys.exit(1)




def compress_vcf(input_vcf, output_vcf_gz):
    with open(input_vcf, 'rb') as f_in:
        with gzip.open(output_vcf_gz, 'wb') as f_out:
            f_out.writelines(f_in)


# Example usage:
if __name__ == "__main__":
    #compress_vcf("src/short.vcf", "short.vcf.gz")
    annotate_vcf_with_vep(
        input_vcf_path="rare_coding_variants.vcf",
        output_json_path="annotations.json",
    )