import json
import csv
import sys
import os
import pandas as pd
#this script processes a JSON file returned my the VEP rest API and extracts the following Data:
# - variant Input string (chromosome, position, reference and alternate alleles)
# - GnomAd-allele frequency
# - GnomadMaxallelfrequency
# - Transcript either the mane if it is present or the first GENCODE basic transcript
# - for this transcript return the SpliceAI score there you sum up all the DS scores
# - for this transcript return the GERP++ score
# - for this transcript return the PolyPhen score if applicable else return 0
# - for this transcript return the REVEL score
# - for this transcript return the CADD score
# - for this transcript return the LoF score
# - for this transcript return the clinvar label


def gnomadFinder(variant):
        colocated = variant.get("colocated_variants", [])
        gnomADAF = 0.0  # Default value if no AF is found
        gnomADmaxAF = "0.0"  # Default value if no max AF is found
        for col in colocated:
            freqs = col.get("frequencies", {})
            for allele, freq_dict in freqs.items():
                # check for the maximum allele frequency
                for item in freq_dict.items():
                    float_value = float(item[1])
                    if float_value > float(gnomADmaxAF):
                        # if the current allele frequency is greater than the current maximum allele frequency  
                        gnomADmaxAF = item[1]
                #print(f"Processing allele: {allele}")
                if "af" in freq_dict:
                    gnomADAF = freq_dict["af"]
                    break
            if 'af' in locals(): 
                #if sf not in the current allel it will continue
                # with the next allele if it exists
                break
            
        #print(f"GnomAD AF found: {gnomADAF}")
        return gnomADAF, gnomADmaxAF

def TranscriptFinder(variant):
    #finds the transcript that is either a MANE transcript or the first GENCODE basic transcript 
    # (similar to the first transcript sice the gencode_basic filter is already applied)
    variant_transcripts = variant.get("transcript_consequences", [])
    maintranscript = ""
    if not variant_transcripts:
        print("No transcripts found for this variant.")
        return None
    # Check for MANE transcripts first
    for transcript in variant_transcripts:
        if "MANE_Select" in transcript.get("mane", []):
        # This is a MANE Select transcript
            #print(f"Found MANE transcript: {transcript.get('transcript_id')}")
            maintranscript = transcript 
            break
    if not maintranscript:  
        maintranscript = variant_transcripts[0] # Return the first transcript if no MANE transcript is found
    #print(f"Selected main transcript: {maintranscript.get('transcript_id')}")
    return maintranscript

def SpliceAIFinder(trancript):
    splice_ai_ds = 0
    #finds the SpliceAI score for a give trancript
    # returns the the maximum of all DS scores
    splice_ai_scores = trancript.get("spliceai", {})
    if splice_ai_scores:
        for score in splice_ai_scores.get("ds", []):
            if score > splice_ai_ds:
                splice_ai_ds = score
    return splice_ai_ds
def GERPFinder(transcript):
    """
    Finds the GERP++ score for a given transcript.
    Returns 0 if not available or not convertible to float.
    Prints a message if the GERP++ score does not exist.
    """
    if "gerp++_rs" not in transcript:
        print(f"for this transcript there is no GERP++ score: {transcript}")
        return 0
    try:
        gerp_score = float(transcript["gerp++_rs"])
    except (TypeError, ValueError):
        print(f"for this transcript the GERP++ score is not convertible: {transcript}")
        gerp_score = 0
    return gerp_score

def PolyPhenFinder(transcript):
    """
    Finds the PolyPhen score for a given transcript.
    Returns 0.5 if not available or not convertible to float.
    """
    if "polyphen_score" not in transcript:
        print("for this transcript there is no polyphen score: " + str(transcript))
        return 0.5
    # Default value for PolyPhen score

    polyphen_score = 0.5 # Default value for PolyPhen score
    try:
        # PolyPhen scores are often under "polyphen_score" or "polyphen"
        polyphen_value = transcript.get("polyphen_score", 0)
        polyphen_score = float(polyphen_value)
    except (TypeError, ValueError):
        print("for this transcript there is no polyphen score: " + transcript)
        polyphen_score = 0.5 
    return polyphen_score

def REVELFinder(transcript):
    """
    Finds the REVEL score for a given transcript.
    Returns 0.5 if not available or not convertible to float.
    Because in REVEL scores 0.5 is the default for no constraint
    """
    if "revel" not in transcript:
        print("for this transcript there is no REVEL score: " + str(transcript))
        return 0.5
    # Default value for REVEL score
    try:
        revel_value = transcript.get("revel", 0)
        revel_score = float(revel_value)
    except (TypeError, ValueError):
        print("for this transcript there is no REVEL score: " + transcript)
        revel_score = 0.5
    return revel_score

def CADDFinder(transcript):
    """
    Finds the CADD phred score for a given transcript.
    Returns 0 if not available or not convertible to float.
    """
    if "cadd_phred" not in transcript:
        print("for this transcript there is no CADD score: " + str(transcript))
        return 0
    try:
        cadd_value = transcript.get("cadd_phred", 0)
        cadd_score = float(cadd_value)
    except (TypeError, ValueError):
        print("for this transcript there is no CADD score: " + transcript)
        cadd_score = 0
    return cadd_score

def ClinvarFinder(variant):
    """
    Finds the ClinVar label (clin_sig) for a given variant.
    Returns empty string if not available.
    """
    colocated = variant.get("colocated_variants", [])
    for col in colocated:
        clinic_score = col.get("clin_sig")
        if clinic_score and isinstance(clinic_score, list) and len(clinic_score) > 0:
            return clinic_score[0]
    return ""

    # Check for ClinVar significance in colocated variants
    
    clinic_score = variant.get("clin_sig", [])
    if not clinic_score:
        return ""       
    # Return the first ClinVar significance label if available
    return clinic_score[0] if isinstance(clinic_score, list) else clinic_score

def LoFFinder(transcript):
    """
    Finds the LoF (Loss-of-Function) label for a given transcript.
    Returns 0 if not available.
    """
    return transcript.get("lof", 0)

def ConsequenceFinder(transcript):
    """
    Finds the consequence type for a given transcript.
    Returns the first consequence type if available, otherwise returns an empty string.
    """
    consequences = transcript.get("consequence_terms", [])
    if consequences:
        return consequences[0]  # Return the first consequence type
    return 0  # Return empty string if no consequences found

def ImpactFinder(transcript):
    """
    Finds the impact of a given transcript.
    Returns the impact value if available, otherwise returns an empty string.
    """
    return transcript.get("impact", 0)  # Return impact or empty string if not found


def process_json_to_dataframe(json_data):
    """
    Process loaded JSON data and return a pandas DataFrame.
    
    Args:
        json_data: Already loaded JSON data (list of variant dictionaries)
        
    Returns:
        pandas.DataFrame: DataFrame with extracted variant features
    """
    # List of fields in the order you want them in the DataFrame
    fieldnames = [
        "variant",
        "gnomad_af",
        "gnomadg_max_af",
        "transcript_id",
        "splice_ai",
        "gerp_score",
        "polyphen_score",
        "revel_score",
        "cadd_score",
        "lof_score",
        "consequence",
        "impact",
        "clinvar_label"
    ]
    
    # List to store all rows
    rows = []
    
    for v in json_data:
        row = {}
        row["variant"] = v.get("input", v.get("id", ""))

        # Try to find the global gnomAD genome AF
        gnomad_af, gnomadg_max_af = gnomadFinder(v)
        row["gnomad_af"] = gnomad_af
        row["gnomadg_max_af"] = gnomadg_max_af

        # Find the main transcript 
        transcript = TranscriptFinder(v)
        if not transcript:
            print(f"No transcript found for variant {row['variant']}. Not Skipping.")
            #continue
            row["transcript_id"] = ""
            row["splice_ai"] = 0
            row["gerp_score"] = 0
            row["polyphen_score"] = 0.5
            row["revel_score"] = 0.5
            row["cadd_score"] = 0
            row["lof_score"] = 0
            row["consequence"] = None
            row["impact"] = None
        else:
            row["transcript_id"] = transcript.get('transcript_id', "")

            # Extract the SpliceAI score
            row["splice_ai"] = SpliceAIFinder(transcript)

            # Extract the GERP++ score  
            row["gerp_score"] = GERPFinder(transcript)

            # Extract the PolyPhen score
            row["polyphen_score"] = PolyPhenFinder(transcript)

            # Extract the REVEL score
            row["revel_score"] = REVELFinder(transcript)

            # Extract the CADD score
            row["cadd_score"] = CADDFinder(transcript)

            # Extract the LoF field
            row["lof_score"] = LoFFinder(transcript)

            # Extract the consequence term
            row["consequence"] = ConsequenceFinder(transcript)

            # Extract the impact
            row["impact"] = ImpactFinder(transcript)

        # Extract the ClinVar label
        row["clinvar_label"] = ClinvarFinder(v)

        rows.append(row)
    
    # Create DataFrame from the list of rows
    df = pd.DataFrame(rows, columns=fieldnames)
    print(f"Processed {len(df)} variants into DataFrame")
    
    return df


def main(input_json, output_csv):
    with open(input_json) as f:
        data = json.load(f)

    # List of fields in the order you want them in the CSV
    fieldnames = [
        "variant",
        "gnomad_af",
        "gnomadg_max_af",
        "transcript_id",
        "splice_ai",
        "gerp_score",
        "polyphen_score",
        "revel_score",
        "cadd_score",
        "lof_score",
        "consequence",
        "impact",
        "clinvar_label"
    ]

    with open(output_csv, "w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter=";")
        writer.writeheader()

        for v in data:
            row = {}
            row["variant"] = v.get("input", v.get("id", ""))

            # Try to find the global gnomAD genome AF
            gnomad_af, gnomadg_max_af = gnomadFinder(v)
            row["gnomad_af"] = gnomad_af
            row["gnomadg_max_af"] = gnomadg_max_af

            # Find the main transcript 
            transcript = TranscriptFinder(v)
            if not transcript:
                print(f"No transcript found for variant {row['variant']}. Skipping.")
                continue
            row["transcript_id"] = transcript.get('transcript_id', "")

            # Extract the SpliceAI score
            row["splice_ai"] = SpliceAIFinder(transcript)

            # Extract the GERP++ score  
            row["gerp_score"] = GERPFinder(transcript)

            # Extract the PolyPhen score
            row["polyphen_score"] = PolyPhenFinder(transcript)

            # Extract the REVEL score
            row["revel_score"] = REVELFinder(transcript)

            # Extract the CADD score
            row["cadd_score"] = CADDFinder(transcript)

            # Extract the LoF field
            row["lof_score"] = LoFFinder(transcript)

            # Extract the consequence term
            row["consequence"] = ConsequenceFinder(transcript)

            # Extract the impact
            row["impact"] = ImpactFinder(transcript)

            # Extract the ClinVar label
            row["clinvar_label"] = ClinvarFinder(v)

            writer.writerow(row)

    print(f"Saved: {output_csv}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: python {os.path.basename(__file__)} <input_json> <output_csv>")
        sys.exit(1)
    input_json = sys.argv[1]
    output_csv = sys.argv[2]
    main(input_json, output_csv)