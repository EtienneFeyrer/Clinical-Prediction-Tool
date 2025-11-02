import gzip
import os


class ClinvarSplit:
    """
    Class to handle splitting of ClinVar VCF data into benign and pathogenic files.
    """

    def __init__(self):
        """
        Initialize with paths for input VCF and output files.
        """
        self.vcf_path = None
        self.output_benign = None
        self.output_pathogenic = None

    def set_vcf_path(self, vcf_path):
        self.vcf_path = vcf_path

    def set_output_benign(self, output_benign):
        self.output_benign = output_benign

    def set_output_pathogenic(self, output_pathogenic):
        self.output_pathogenic = output_pathogenic

    def split_vcf(self):
        """
        Splits the ClinVar VCF file into two separate files based on significance.
        Writes benign/likely_benign and pathogenic/likely_pathogenic variants.
        Handles both .vcf and .vcf.gz files using gzip if needed.
        Header lines are collected and written to both output files, just like in Batches.py.
        Uses robust INFO field parsing (dictionary lookup) for CLNSIG.
        """
        if not self.vcf_path or not self.output_benign or not self.output_pathogenic:
            raise ValueError("VCF path and output paths must be set before splitting.")

        # Use gzip.open for .vcf.gz, else open for .vcf
        open_func = gzip.open if self.vcf_path.endswith('.gz') else open
        headers = []
        benign_records = []
        pathogenic_records = []
        with open_func(self.vcf_path, 'rt') as vcf_file:
            for line in vcf_file:
                if line.startswith('#'):
                    headers.append(line)
                else:
                    fields = line.strip().split('\t')  # Split by semicolon for VCF INFO field
                    info = fields[7]
                    info_items = info.split(';')
                    info_dict = dict(item.split('=', 1) for item in info_items if '=' in item)
                    clnsig = info_dict.get('CLNSIG', None)
                    if clnsig:
                        if any(sig in clnsig for sig in ["benign", "likely_benign"]):
                            benign_records.append(line)
                        elif any(sig in clnsig for sig in ["pathogenic", "likely_pathogenic"]):
                            pathogenic_records.append(line)

        with open(self.output_benign, 'w') as benign_file:
            benign_file.writelines(headers)
            benign_file.writelines(benign_records)
        with open(self.output_pathogenic, 'w') as pathogenic_file:
            pathogenic_file.writelines(headers)
            pathogenic_file.writelines(pathogenic_records)
        print(f"Wrote {len(benign_records)} benign and {len(pathogenic_records)} pathogenic records.")

# This approach now robustly parses the INFO field for CLNSIG, regardless of its position.
