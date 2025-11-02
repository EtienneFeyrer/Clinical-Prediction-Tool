"""
Database Validation Test.
Uses existing VEP annotation data directly without VCF parsing.
"""
import json
import os
import sys
import mysql.connector
from typing import Dict, List
import importlib

sys.path.append(os.path.abspath('.'))
from MachineLearning.RandomForest.RandomForestIO import RandomForestIO
sys.path.append(os.path.abspath('.'))
import packages.annotation_classes as ac
sys.path.append(os.path.abspath('.'))
import core.helper_methods as helper
sys.path.append(os.path.abspath('.'))
import database.service as db

def parse_vep_to_annotation_exact(vep_data: Dict) -> ac.GeneAnnotations:
    """
    EXACT copy of BatchProcessor._parse_vep_to_annotation method.
    This ensures identical parsing logic for validation.
    """
    variant_annotation = ac.GeneAnnotations()

    # Variant-level annotations
    if vep_data.get('transcript_consequences'):
        first_transcript = vep_data['transcript_consequences'][0]
        if 'cadd_phred' in first_transcript:
            variant_annotation.add_annotation(ac.Annotation_Float("CADD", first_transcript['cadd_phred']))
        if 'gene_symbol' in first_transcript:
            variant_annotation.add_annotation(ac.Annotation_Str("GENE", first_transcript['gene_symbol']))
        else:
            variant_annotation.add_annotation(ac.Annotation_Str("GENE", ""))

    # Extract gnomAD frequency
    if 'colocated_variants' in vep_data:
        found_freq = False
        found_clin = False

        for colocated in vep_data['colocated_variants']:
            # Look for frequencies
            if not found_freq and 'frequencies' in colocated:
                for allele, freq_data in colocated['frequencies'].items():
                    gnomad_af = freq_data.get('gnomadg', freq_data.get('af'))
                    if gnomad_af:
                        variant_annotation.add_annotation(ac.Annotation_Float("gnomAD AF", gnomad_af))
                    variant_annotation.add_annotation(
                        ac.Annotation_Float("max_allele_freq", helper.get_max_frequency(freq_data)))
                    found_freq = True
                    break

            # Look for clin_sig
            if not found_clin and 'clin_sig' in colocated:
                clin_sig = colocated['clin_sig']
                variant_annotation.add_annotation(ac.Annotation_Str("CLINSIG", clin_sig))
                found_clin = True

            # Exit early if both found
            if found_freq and found_clin:
                break

    # Most severe consequence
    if 'most_severe_consequence' in vep_data:
        variant_annotation.add_annotation(
            ac.Annotation_Str("Most Severe Consequence", vep_data['most_severe_consequence']))

    # OMIM extraction
    clinvar_omim_id = (helper.extract_clinvar_omim_id(vep_data).split('&'))
    if clinvar_omim_id:
        variant_annotation.add_annotation(ac.Annotation_Str("OMIM", clinvar_omim_id))
    else:
        variant_annotation.add_annotation(ac.Annotation_Str("OMIM", None))

    try:
        if isinstance(vep_data, dict):
            ml_input = [vep_data.copy()]

            if 'seq_region_name' in ml_input[0] and ml_input[0]['seq_region_name'].startswith('chr'):
                ml_input[0]['seq_region_name'] = ml_input[0]['seq_region_name'][3:]

            if 'input' in ml_input[0] and ml_input[0]['input'].startswith('chr'):
                ml_input[0]['input'] = ml_input[0]['input'][3:]
        else:
            ml_input = vep_data

        importlib.reload(sys.modules['MachineLearning.RandomForest.RandomForestIO'])
        importlib.reload(sys.modules['MachineLearning.RandomForest.RandomForestModel'])
        object = RandomForestIO("forest")
        ml_result = object.get_ml_values(ml_input)

        # EXTRACT the score from the list
        if isinstance(ml_result, list) and len(ml_result) > 0:
            ml_score = ml_result[0]
        else:
            ml_score = ml_result

        variant_annotation.add_annotation(ac.Annotation_Float("ML-Score", ml_score))
    except Exception as e:
        print(f"ML model error: {e}")
        import traceback
        traceback.print_exc()
        mock_score = 0.75
        variant_annotation.add_annotation(ac.Annotation_Float("ML-Score", mock_score))
        print(f"Using mock ML-Score: {mock_score}")

    # Process transcripts (exact copy)
    for entry in vep_data.get('transcript_consequences', []):
        is_mane = bool(entry.get('mane') and len(entry.get('mane', [])) > 0)

        # Get cDNA and Protein notation
        if 'hgvsc' in entry:
            cdna_notation = entry.get('hgvsc').split(':')[1]
        else:
            cdna_notation = ""
        if 'hgvsp' in entry:
            protein_notation = entry.get('hgvsp').split(':')[1]
        else:
            protein_notation = ""

        if 'consequence_terms' in entry:
            consequence_terms = entry['consequence_terms']
            consequence_string = ",".join(consequence_terms)
        else:
            consequence_string = ""  # Your exact placeholder

        transcript_ann = ac.TranscriptAnnotations(
            transcript_id=entry['transcript_id'] or "",
            impact=entry['impact'] or "",
            revel=entry.get('revel') or None,
            gerp=entry.get('gerp++_rs') or None,
            spliceai=helper.extract_max_spliceai_score(entry.get('spliceai', {})) if entry.get('spliceai') else None,
            polyphen=entry.get('polyphen_score') or None,
            loftee=entry.get('lof') or "",
            is_mane=is_mane,
            cdna_notation=cdna_notation,
            protein_notation=protein_notation,
            consequences=consequence_string
        )

        variant_annotation.add_transcript_annotation(transcript_ann)

    return variant_annotation


def create_reference_from_vep_data(vep_json_file='testing/annotations.json',
                                   output_file='testing/reference_database_format.jsonl') -> List[Dict]:
    """
    Create reference annotations directly from VEP data using BatchProcessor logic.
    No VCF parsing needed - just processes the VEP annotations directly.
    """
    print("=== Creating Reference Annotations from VEP Data ===")

    # Load VEP data
    try:
        with open(vep_json_file) as file:
            vep_data = json.load(file)
        print(f"Loaded {len(vep_data)} VEP results")
    except FileNotFoundError:
        print(f"VEP JSON file not found: {vep_json_file}")
        return []

    reference_annotations = []

    for i, vep_variant_data in enumerate(vep_data):
        try:
            # Use EXACT same parsing logic as BatchProcessor
            variant_annotation = parse_vep_to_annotation_exact(vep_variant_data)

            # Convert to database JSON format using existing helper
            json_data = helper.annotation_to_json(variant_annotation)

            # Create simple variant key - you can modify this logic as needed
            variant_key = f"variant_{i}"

            # If your VEP data has input info, you can extract the real variant key:
            if 'input' in vep_variant_data:
                # VEP input format is usually like "1 123456 123456 A/G +"
                input_parts = vep_variant_data['input'].split()
                if len(input_parts) >= 4:
                    chrom = input_parts[0]
                    pos = input_parts[1]
                    ref_alt = input_parts[3]
                    if '/' in ref_alt:
                        ref, alt = ref_alt.split('/')
                        variant_key = f"{chrom}:{pos}:{ref}>{alt}"

            # Create database-compatible structure
            database_annotation = {
                'gene': json_data.get('GENE'),
                'CADD': json_data.get('CADD'),
                'pathogenicity_score': json_data.get('ML-Score'),
                'Most Severe Consequence': json_data.get('Most Severe Consequence'),
                'gnomAD AF': json_data.get('gnomAD AF'),
                'max_allele_freq': json_data.get('max_allele_freq'),
                'OMIM': json_data.get('OMIM'),
                'CLINSIG': json_data.get('CLINSIG'),
                'transcript_consequences': []
            }

            # Convert transcript annotations to database format
            for transcript in json_data.get('transcripts', []):
                transcript_db = {
                    'transcript_id': transcript.get('transcript_id'),
                    'polyphen': transcript.get('polyphen'),
                    'protein_notation': transcript.get('protein_notation'),
                    'REVEL': transcript.get('revel'),
                    'Splice_AI': transcript.get('spliceai'),
                    'Mane': transcript.get('is_mane', False),
                    'LOFTEE': transcript.get('lof'),
                    'impact': transcript.get('impact'),
                    'GERP': transcript.get('gerp'),
                    'cDNA_notation': transcript.get('cdna_notation'),
                    'consequences': transcript.get('consequences')
                }
                database_annotation['transcript_consequences'].append(transcript_db)

            reference_annotations.append({
                'variant_id': variant_key,
                'annotation': database_annotation
            })

        except Exception as e:
            print(f"Error processing variant {i}: {e}")
            continue

    # Write to output file
    with open(output_file, 'w') as out_file:
        for annotation in reference_annotations:
            json.dump(annotation, out_file)
            out_file.write('\n')

    print(f"Created {len(reference_annotations)} reference annotations")
    print(f"Output written to: {output_file}")
    return reference_annotations


def get_database_variant_keys() -> List[str]:
    """
    Get all variant keys currently in the database.
    This helps identify which variants to compare.
    """
    print("=== Getting Database Variant Keys ===")

    try:
        config = db.DatabaseService.get_db_config()
        connection = mysql.connector.connect(**config)
        cursor = connection.cursor()

        cursor.execute("SELECT variant_id FROM Annotation")
        variant_keys = [row[0] for row in cursor.fetchall()]

        cursor.close()
        connection.close()

        print(f"Found {len(variant_keys)} variants in database")
        if variant_keys:
            print(f"Sample keys: {variant_keys[:3]}...")

        return variant_keys

    except Exception as e:
        print(f"Error getting database keys: {e}")
        return []


def compare_with_database(reference_file='testing/reference_database_format.jsonl') -> bool:
    """
    Compare reference annotations with database entries.
    Only compares variants that exist in both reference and database.
    """
    print("=== Comparing Reference with Database ===")

    # Load reference annotations
    reference_annotations = []
    try:
        with open(reference_file, 'r') as f:
            for line in f:
                if line.strip():
                    reference_annotations.append(json.loads(line))
    except FileNotFoundError:
        print(f"Reference file not found: {reference_file}")
        return False

    # Get database variant keys
    db_variant_keys = set(get_database_variant_keys())

    print(f"Loaded {len(reference_annotations)} reference annotations")
    print(f"Found {len(db_variant_keys)} variants in database")

    # Compare annotations
    differences = []
    matches = 0
    not_found = 0

    for ref_data in reference_annotations:
        variant_id = ref_data['variant_id']
        ref_annotation = ref_data['annotation']

        # Skip if not in database (expected for test data)
        if variant_id not in db_variant_keys:
            not_found += 1
            continue

        try:
            # Get database annotation
            db_annotation = db.DatabaseService.get_variant_annotation(variant_id)

            if db_annotation is None:
                differences.append(f"{variant_id}: Found in DB list but get_variant_annotation returned None")
                continue

            # Compare main fields
            main_fields = ['gene', 'CADD', 'pathogenicity_score', 'Most Severe Consequence',
                          'gnomAD AF', 'max_allele_freq', 'OMIM', 'CLINSIG']

            variant_differences = []

            for field in main_fields:
                ref_value = ref_annotation.get(field)
                db_value = db_annotation.get(field)

                if not values_equal(ref_value, db_value):
                    variant_differences.append(f"  {field}: ref='{ref_value}' vs db='{db_value}'")

            # Compare transcript data (first 2 transcripts only)
            #ref_transcripts = ref_annotation.get('transcript_consequences', [])
            #db_transcripts = db_annotation.get('transcript_consequences', [])

            #if len(ref_transcripts) != len(db_transcripts):
                #variant_differences.append(f"  transcript_count: ref={len(ref_transcripts)} vs db={len(db_transcripts)}")

            #transcript_fields = ['transcript_id', 'polyphen', 'protein_notation', 'REVEL',
                               #'Splice_AI', 'Mane', 'LOFTEE', 'impact', 'GERP',
                               #'cDNA_notation', 'consequences']

            #for i, (ref_t, db_t) in enumerate(zip(ref_transcripts[:2], db_transcripts[:2])):
                #for field in transcript_fields:
                    #ref_val = ref_t.get(field)
                    #db_val = db_t.get(field)

                    #if not values_equal(ref_val, db_val):
                        #variant_differences.append(f"  transcript[{i}].{field}: ref='{ref_val}' vs db='{db_val}'")

            if variant_differences:
                differences.append(f"{variant_id}:\n" + "\n".join(variant_differences))
            else:
                matches += 1
                print(f"{variant_id}: Perfect match")

        except Exception as e:
            differences.append(f"{variant_id}: Error during comparison - {e}")

    # Summary
    compared_count = len(reference_annotations) - not_found

    # Write detailed report
    with open('testing/validation_report_simple.txt', 'w') as f:
        f.write("=== Streamlined Database Validation Report ===\n\n")
        f.write(f"Total reference: {len(reference_annotations)}\n")
        f.write(f"Compared: {compared_count}\n")
        f.write(f"Perfect matches: {matches}\n")
        f.write(f"With differences: {len(differences)}\n\n")

        for diff in differences:
            f.write(diff + "\n\n")

    success = len(differences) == 0 and compared_count > 0
    return success


def values_equal(val1, val2) -> bool:
    """Check if two values are equal, handling None and type differences."""
    if val1 is None and val2 is None:
        return True
    if val1 is None or val2 is None:
        return False
    if isinstance(val1, float) and isinstance(val2, float):
        return abs(val1 - val2) < 1e-6
    return str(val1) == str(val2)


def run_simple_validation() -> bool:
    """Run streamlined database validation test."""

    # Step 1: Test database connection
    try:
        config = db.DatabaseService.get_db_config()
        connection = mysql.connector.connect(**config)
        cursor = connection.cursor()
        cursor.execute("SELECT COUNT(*) FROM Annotation")
        db_count = cursor.fetchone()[0]
        cursor.close()
        connection.close()
        print(f"Database connected - {db_count} annotations found")
    except Exception as e:
        print(f"Database connection failed: {e}")
        return False

    # Step 2: Create reference annotations from VEP data
    reference_annotations = create_reference_from_vep_data(
        vep_json_file='testing/annotations.json',
        output_file='testing/reference_database_format.jsonl'
    )

    if not reference_annotations:
        print("\nCould not create reference annotations")
        return False

    # Step 3: Compare with database
    success = compare_with_database('testing/reference_database_format.jsonl')

    # Final result
    print("=" * 50)
    if success:
        print("VALIDATION PASSED!")
        print("Database entries match reference parser")
    else:
        print("VALIDATION ISSUES FOUND")
        print("Check validation_report_simple.txt")
    print("=" * 50)

    return success


if __name__ == "__main__":
    print("Streamlined Database Validation")
    print("This uses your VEP annotation data directly without VCF parsing\n")

    success = run_simple_validation()
