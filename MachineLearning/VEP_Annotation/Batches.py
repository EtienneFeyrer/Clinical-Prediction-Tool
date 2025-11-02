import gzip
import os

class Batch:
    def __init__(self, header, records):
        self.header = header  # List of header lines
        self.records = records  # List of record lines

    def write_to_file(self, filepath):
        with open(filepath, "w") as f:
            f.writelines(self.header)
            f.writelines(self.records)

class VCFBatchSplitter:
    def __init__(self):
        self.vcf_path = None
        self.batch_size = 1000  # Default batch size

    def set_vcf(self, vcf_path):
        self.vcf_path = vcf_path

    def set_batch_size(self, batch_size):
        self.batch_size = batch_size

    def split(self):
        if not self.vcf_path:
            raise ValueError("VCF path must be set before splitting.")

        batch_dir = self.vcf_path.replace(".vcf.gz", "_batches").replace(".vcf", "_batches")
        os.makedirs(batch_dir, exist_ok=True)

        # Read header lines and records
        headers = []
        records = []
        open_func = gzip.open if self.vcf_path.endswith('.gz') else open
        with open_func(self.vcf_path, "rt", encoding="latin-1") as f:
            for line in f:
                if line.startswith("#"):
                    headers.append(line)
                else:
                    records.append(line)

        # Split into batches using the Batch class
        total_batches = (len(records) + self.batch_size - 1) // self.batch_size
        for i in range(total_batches):
            batch_records = records[i*self.batch_size : (i+1)*self.batch_size]
            batch = Batch(headers, batch_records)
            batch_file = os.path.join(
                batch_dir, f"{os.path.basename(self.vcf_path).replace('.vcf.gz','').replace('.vcf','')}_batch_{i+1}.vcf"
            )
            batch.write_to_file(batch_file)
            if (i + 1) % 10000 == 0:
                print(f"Saved {i + 1} batches so far...")

# Example usage:
# splitter = VCFBatchSplitter()
# splitter.set_vcf("../Data/Main/clinvar.vcf.gz")
# splitter.set_batch_size(1000)
# splitter.split()
# splitter.set_vcf("../Data/Main/clinvar.vcf")
# splitter.split()