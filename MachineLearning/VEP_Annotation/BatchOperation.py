import os
import sys
import importlib
# Ensure the parent directory is in the system path to import VEPAnnotator
sys.path.append(os.path.abspath('..'))
import VEP_Annotation.APIrequest
importlib.reload(VEP_Annotation.APIrequest)
from VEP_Annotation import APIrequest
from concurrent.futures import ThreadPoolExecutor, as_completed
from VEP_Annotation import VEP_extract

def annotate_vcf(file, path_batches, field_list):
    annotator = APIrequest.VEPAnnotator()
    annotator.set_vcf(os.path.join(path_batches, file))
    annotator.set_output(os.path.join(path_batches, file.replace(".vcf", "_annotated.json")))
    annotator.set_fields(field_list)
    annotator.run()
    return file



def annotate_vcf_batches(path_batches, field_list, print_progress=True, batch_label="batches", max_workers=15):
    """
    Annotate all VCF files in the given directory with gnomAD fields using VEPAnnotator, in parallel.
    """
    vcf_files = [f for f in os.listdir(path_batches) if f.endswith(".vcf")]

    counter = 0
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(annotate_vcf, file, path_batches, field_list): file for file in vcf_files}
        for future in as_completed(futures):
            counter += 1
            if print_progress and counter % 10 == 0:
                print(f"Annotated {counter} {batch_label} in the folder '{path_batches}'.")

def extract_single(file, path_batches):
    """
    Convert the annotated JSON file to a CSV file.
    """
    json_path = os.path.join(path_batches, file)
    output_csv = json_path.replace("_annotated.json", "_reduced.csv")
    VEP_Annotation.VEP_extract.main(input_json=json_path, output_csv=output_csv)


def extract_annotated_csv(path_batches, max_workers=15):
    """
    Convert the annotated JSON files from a given directory to CSV files.
    """
    json_files = [f for f in os.listdir(path_batches) if f.endswith("_annotated.json")]
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(extract_single, file, path_batches): file for file in json_files}
        for future in as_completed(futures):
            file = futures[future]
            try:
                future.result()
                print(f"Converted {file} to CSV.")
            except Exception as e:
                print(f"Failed to convert {file}: {e}")


if __name__ == "__main__":
    # Example usage:
    # annotate_vcf_batches("/path/to/batches", field_list=["GnomADg_AF"])
    pass
