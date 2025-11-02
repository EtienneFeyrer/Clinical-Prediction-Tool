#This Parser is used to parse the VEP string values into numerical values
#required for the following annotaions:
#All unique impact classes: ['HIGH', 'LOW', 'MODERATE', 'MODIFIER']
#All unique LOF classes: ['0', 'HC', 'LC']
#All unique consequence classes: ['3_prime_UTR_variant', '5_prime_UTR_variant', 'coding_sequence_variant', 'downstream_gene_variant', 'frameshift_variant', 'inframe_deletion', 'inframe_insertion', 'intron_variant', 'mature_miRNA_variant', 'missense_variant', 'non_coding_transcript_exon_variant', 'protein_altering_variant', 'splice_acceptor_variant', 'splice_donor_5th_base_variant', 'splice_donor_region_variant', 'splice_donor_variant', 'splice_polypyrimidine_tract_variant', 'splice_region_variant', 'start_lost', 'stop_gained', 'stop_lost', 'stop_retained_variant', 'synonymous_variant', 'transcript_ablation', 'upstream_gene_variant']
#All unique ClinVar classes: ['Benign/Likely_benign', 'Benign/Likely_benign|drug_response', 'Benign/Likely_benign|drug_response|other', 'Benign/Likely_benign|other', 'Benign/Likely_benign|other|risk_factor', 'Benign/Likely_benign|risk_factor', 'Conflicting_classifications_of_pathogenicity', 'Conflicting_classifications_of_pathogenicity|Affects', 'Conflicting_classifications_of_pathogenicity|association', 'Conflicting_classifications_of_pathogenicity|association|risk_factor', 'Conflicting_classifications_of_pathogenicity|drug_response', 'Conflicting_classifications_of_pathogenicity|other', 'Conflicting_classifications_of_pathogenicity|other|risk_factor', 'Conflicting_classifications_of_pathogenicity|risk_factor', 'Likely_benign', 'Likely_benign|association', 'Likely_benign|drug_response', 'Likely_benign|other', 'Likely_benign|risk_factor', 'Likely_pathogenic', 'Likely_pathogenic/Likely_risk_allele', 'Likely_pathogenic|Affects', 'Likely_pathogenic|association', 'Likely_pathogenic|drug_response', 'Likely_pathogenic|other', 'Likely_pathogenic|risk_factor', 'Pathogenic/Likely_pathogenic', 'Pathogenic/Likely_pathogenic/Likely_risk_allele', 'Pathogenic/Likely_pathogenic/Pathogenic', 'Pathogenic/Likely_pathogenic|association', 'Pathogenic/Likely_pathogenic|other', 'Pathogenic/Likely_pathogenic|risk_factor']
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import OneHotEncoder

def parse_loftee(loftee_str):
    """
    Parses the LOFTEE string and returns a numerical value.
    Returns 0 if not available.
    """
    if not loftee_str:
        return 0
    loftee_str = loftee_str.strip()
    if loftee_str == "HC":
        return 1  # High confidence
    elif loftee_str == "LC":
        return 0.5  # Low confidence
    else:
        return 0  # Not a loss-of-function variant

def parse_impact(impact_str):
    """
    Parses the impact string and returns a numerical value.
    Returns 0 if not available.
    """
    if not impact_str:
        return 0
    impact_str = impact_str.strip().upper()
    if impact_str == "HIGH":
        return 1
    elif impact_str == "MODERATE":
        return 0.5
    elif impact_str == "LOW":
        return 0.25
    elif impact_str == "MODIFIER":
        return 0.1
    else:
        return 0  # Not a recognized impact type
    
def parse_hot_encoding(column):
    print("Parsing hot encoding for column:")
    split_consequences = column.str.split(',', expand=False)
    #print first 10 rows of split_consequences
    #print(split_consequences.head(10))
    return column.str.get_dummies(sep=',')

#recieves a column and a list of classes(possible values for the specified 
#column) and returns a dataframe with one-hot encoded columns for each class
#returns a dataframe with one-hot encoded columns for each class
def single_hot_encoding(column, classes):
    """
    Creates one-hot encoded columns for a single-value column.
    
    Args:
        column: pandas Series with single values (not comma-separated)
        classes: list of all possible class values
        
    Returns:
        pandas DataFrame with one-hot encoded columns for each class
    """
    # Create a DataFrame with columns for each class, initialized with 0
    df = pd.DataFrame(0, index=column.index, columns=classes)
    
    # For each row, set the corresponding class column to 1
    for idx, value in column.items():
        if value in classes:
            df.loc[idx, value] = 1
        # If value is not in classes, all columns remain 0 (unknown/other category)
    
    return df

def parse_consequence(consequence_str):
    """
    Parses the consequence string and returns a numerical value.
    Returns 0 if not available.
    """
    if not consequence_str:
        return 0
    consequence_str = consequence_str.strip().lower()
    consequence_map = {
        "3_prime_utr_variant": 1,
        "5_prime_utr_variant": 1,
        "coding_sequence_variant": 2,
        "downstream_gene_variant": 0.5,
        "frameshift_variant": 3,
        "inframe_deletion": 2.5,
        "inframe_insertion": 2.5,
        "intron_variant": 0.1,
        "mature_mirna_variant": 1,
        "missense_variant": 2,
        "non_coding_transcript_exon_variant": 1,
        "protein_altering_variant": 2,
        "splice_acceptor_variant": 3,
        "splice_donor_5th_base_variant": 3,
        "splice_donor_region_variant": 3,
        "splice_donor_variant": 3,
        "splice_polypyrimidine_tract_variant": 3,
        "splice_region_variant": 3,
        "start_lost": 3,
        "stop_gained": 3,
        "stop_lost": 3,
        "stop_retained_variant": 1,
        "synonymous_variant": 0.1,
        "transcript_ablation": 4, 
        "upstream_gene_variant": 0.5
    }
    return consequence_map.get(consequence_str, 0)  # Default to 0 if not found

def parse_clinvar(clinvar_str):
    """
    Parses the ClinVar string and returns a numerical value.
    Returns 0 if not available.
    """
    if not clinvar_str:
        return 0
    clinvar_str = clinvar_str.strip().lower()
    clinvar_map = {
        "benign/likely_benign": 0,
        "benign/likely_benign|drug_response": 0,
        "benign/likely_benign|drug_response|other": 0,
        "benign/likely_benign|other": 0,
        "benign/likely_benign|other|risk_factor": 0,
        "benign/likely_benign|risk_factor": 0,
        "likely_benign": 0.25,
        "likely_benign|association": 0.25,
        "likely_benign|drug_response": 0.25,
        "likely_benign|other": 0.25,
        "likely_benign|risk_factor": 0.25,
        "likely_pathogenic": 0.75,
        "likely_pathogenic/likely_risk_allele": 0.75,
        "likely_pathogenic|affects": 0.75,
        "likely_pathogenic|association": 0.75,
        "likely_pathogenic|drug_response": 0.75,
        "likely_pathogenic|other": 0.75,
        "likely_pathogenic|risk_factor": 0.75,
        "pathogenic/likely_pathogenic": 1,
        "pathogenic/likely_pathogenic/likely_risk_allele": 1,
        "pathogenic/likely_pathogenic/pathogenic": 1,
        "pathogenic/likely_pathogenic|association": 1,
        "pathogenic/likely_pathogenic|other": 1,
        "pathogenic/likely_pathogenic|risk_factor": 1
    }
    return clinvar_map.get(clinvar_str, 0)  # Default to 0 if not found

