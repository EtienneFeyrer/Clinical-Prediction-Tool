import gzip
import json
import requests
import sys

#default annotation is always with GERP++

class VEPAnnotator:
    def __init__(self):
        self.vcf_path = None
        self.output_json = None
        self.fields = []
        self.vep_url = "https://rest.ensembl.org/vep/human/region"

    def set_vcf(self, vcf_path):
        self.vcf_path = vcf_path
#the output path is the path to the output JSON file
    def set_output(self, output_json):
        self.output_json = output_json
        if not self.output_json:
            raise ValueError("Output path not set.")
        if self.output_json.endswith(".json") is False:
            raise ValueError("Output file must have a .json extension.")
# the field parameter is a list of fields to request from VEP
        # e.g. GnomADg_AF, GENCODE_basic and so on
        
    def set_fields(self, fields):
        self.fields = fields
        #activate the fields that are requested annotations
        if not self.fields:
            print("No fields set for VEP annotation. Defaulting to GnomADg_AF and GENCODE_basic.")
            self.fields = ["GnomADg_AF", "GENCODE_basic"]
        else:
            # Ensure the URL includes activation
           fields_params = []
        for i, field in enumerate(self.fields):
            fields_params.append(f"{field}=1")
        self.vep_url += "?" + "&".join(fields_params)
        self.vep_url += "&dbNSFP=GERP%2B%2B_RS"

    def extract_variants(self):
        if not self.vcf_path:
            raise ValueError("VCF path not set.")
        opener = gzip.open if self.vcf_path.endswith(".gz") else open
        variants = []
        with opener(self.vcf_path, "rt", encoding="latin-1") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                cols = line.strip().split('\t')
                chrom, pos, ref, alt = cols[0], cols[1], cols[3], cols[4]
                for a in alt.split(","):  # handle multi-allelic sites
                    start = int(pos)
                    end = start + len(ref) - 1
                    variants.append(f"{chrom} {start} {end} {ref}/{a} +")
        return variants

    def run(self):
        if not (self.vcf_path and self.output_json and self.fields):
            raise ValueError("VCF path, output path, or fields not set.")
        variants = self.extract_variants()
        payload = {
            "variants": variants,
            "fields": self.fields
        }
        print("Payload to be sent to VEP API:")
        print(json.dumps(payload, indent=2))
        headers = {"Content-Type": "application/json"}
        response = requests.post(self.vep_url, headers=headers, data=json.dumps(payload))
        if response.status_code == 200:
            with open(self.output_json, "w") as out_f:
                out_f.write(response.text)
            print(f"VEP API annotation saved to {self.output_json}")
        else:
            print(f"Error: {response.status_code} {response.text}")
            sys.exit(1)

# Optional: Command-line interface for backward compatibility
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Annotate VCF using Ensembl VEP REST API")
    parser.add_argument("--input_vcf", required=True, help="Path to input VCF.gz file")
    parser.add_argument("--output_json", required=True, help="Path to output JSON file")
    parser.add_argument("--fields", required=True, nargs="+", help="Fields to request from VEP (space-separated)")
    args = parser.parse_args()

    annotator = VEPAnnotator()
    annotator.set_vcf(args.input_vcf)
    annotator.set_output(args.output_json)
    annotator.set_fields(args.fields)
    annotator.run()