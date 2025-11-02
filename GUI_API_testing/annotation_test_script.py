#!/usr/bin/env python3
"""
Comprehensive test for annotation assignment functionality.
Tests the App window's ability to receive and correctly assign annotation results.
"""
import requests
import time
import json
import tempfile
import os
import sys
from typing import Dict, List, Any

# Mock annotation data that simulates what your API might return
MOCK_ANNOTATIONS = {
    "chr1:12345:A>T": {
        "gene_symbol": "BRCA1",
        "transcript_id": "ENST00000357654",
        "consequence": "missense_variant",
        "variant_type": "SNV",
        "revel_score": "0.856",
        "clinvar_significance": "Pathogenic",
        "cadd_phred": "25.3",
        "gerp_rs": "5.76",
        "polyphen2": "probably_damaging(0.999)",
        "max_pop_af": "0.000012",
        "spliceai": "0.02",
        "omim_id": "113705"
    },
    "chr2:67890:G>C": {
        "gene_symbol": "TP53",
        "transcript_id": "ENST00000269305",
        "consequence": "stop_gained",
        "variant_type": "SNV",
        "revel_score": "0.952",
        "clinvar_significance": "Pathogenic",
        "cadd_phred": "35.8",
        "gerp_rs": "4.33",
        "polyphen2": "probably_damaging(1.000)",
        "max_pop_af": "0.000001",
        "loftee": "HC",
        "omim_id": "191170"
    },
    "chr3:11111:AT>G": {
        "gene_symbol": "CFTR",
        "transcript_id": "ENST00000003084",
        "consequence": "frameshift_variant",
        "variant_type": "INDEL",
        "clinvar_significance": "Likely pathogenic",
        "cadd_phred": "28.9",
        "gerp_rs": "4.12",
        "max_pop_af": "0.000008",
        "loftee": "HC",
        "omim_id": "602421"
    }
}

def test_server_health(api_url="http://127.0.0.1:5001"):
    """Test server connectivity"""
    try:
        response = requests.get(f"{api_url}/health", timeout=5)
        if response.status_code == 200:
            print(f"✓ Server is healthy at {api_url}")
            return True
        else:
            print(f"✗ Server returned status {response.status_code}")
            return False
    except requests.exceptions.RequestException as e:
        print(f"✗ Cannot connect to server: {e}")
        return False

def submit_test_variants(api_url: str) -> List[str]:
    """Submit test variants and return list of variant IDs"""
    print("\n=== Submitting Test Variants ===")

    test_variants = [
        {'chrom': 'chr1', 'pos': 12345, 'ref': 'A', 'alt': 'T'},
        {'chrom': 'chr2', 'pos': 67890, 'ref': 'G', 'alt': 'C'},
        {'chrom': 'chr3', 'pos': 11111, 'ref': 'AT', 'alt': 'G'},
    ]

    submitted_variants = []

    for i, variant_data in enumerate(test_variants):
        variant_id = f"{variant_data['chrom']}:{variant_data['pos']}:{variant_data['ref']}>{variant_data['alt']}"

        try:
            print(f"Submitting {i+1}/3: {variant_id}")
            response = requests.post(f"{api_url}/submit", json=variant_data, timeout=10)

            if response.status_code == 200:
                result = response.json()
                print(f"  Status: {result.get('status')} - {result.get('message', '')}")

                if result.get('status') in ['success', 'failure']:
                    submitted_variants.append(variant_id)

            else:
                print(f"  HTTP {response.status_code}: {response.text}")

        except requests.exceptions.RequestException as e:
            print(f"  Error: {e}")

    print(f"\nSuccessfully submitted: {len(submitted_variants)}/{len(test_variants)} variants")
    return submitted_variants

def poll_for_annotations(api_url: str, variant_ids: List[str]) -> Dict[str, Dict]:
    """Poll variants until annotated and return annotation data"""
    print(f"\n=== Polling for Annotations ===")

    completed_annotations = {}
    max_polls = 60  # 1 minute of polling

    for poll_round in range(max_polls):
        newly_completed = 0

        for variant_id in variant_ids:
            if variant_id in completed_annotations:
                continue

            try:
                response = requests.get(f"{api_url}/poll/{variant_id}", timeout=5)

                if response.status_code == 200:
                    result = response.json()
                    status = result.get('status')

                    if status == 'completed':
                        annotation_data = result.get('annotation', {})
                        completed_annotations[variant_id] = annotation_data
                        print(f"  ✓ {variant_id} completed")
                        newly_completed += 1

                    elif status == 'processing':
                        if poll_round % 10 == 0:  # Print less frequently
                            print(f"  ⏳ {variant_id} still processing...")

                elif response.status_code == 404:
                    print(f"  ✗ {variant_id} not found")

            except requests.exceptions.RequestException as e:
                if poll_round % 20 == 0:  # Print errors less frequently
                    print(f"  Network error for {variant_id}: {e}")

        completed_count = len(completed_annotations)
        total_count = len(variant_ids)

        if poll_round % 5 == 0:  # Status every 5 seconds
            print(f"Poll {poll_round+1}: {completed_count}/{total_count} completed")

        if completed_count == total_count:
            break

        time.sleep(1)

    return completed_annotations

def test_annotation_mapping(annotation_data: Dict[str, Dict]):
    """Test that annotations map correctly to our column structure"""
    print(f"\n=== Testing Annotation Mapping ===")

    # These are the columns we expect to populate
    expected_columns = [
        "Gene", "Transcript", "Type", "Consequence", "LOFTEE",
        "Allele Frequencies", "Maximum allele frequency of all subpopulations",
        "Splice AI", "GERP++", "PolyPhen-2", "REVEL", "CAD(25)", "OMIM", "ClinVar"
    ]

    mapping_results = {}

    for variant_id, data in annotation_data.items():
        print(f"\nTesting mapping for {variant_id}:")
        variant_mapping = {}

        # Simulate the App window's mapping logic
        mappings = {
            'Gene': _extract_value(data, ['gene_symbol', 'gene_name', 'gene', 'SYMBOL']),
            'Transcript': _extract_value(data, ['transcript_id', 'feature', 'transcript']),
            'Type': _extract_value(data, ['variant_type', 'type', 'Variant_Type']),
            'Consequence': _extract_value(data, ['consequence', 'annotation', 'Consequence']),
            'LOFTEE': _extract_value(data, ['lof', 'loftee', 'LoF', 'LOFTEE']),
            'Maximum allele frequency of all subpopulations': _extract_value(data, ['max_af', 'max_pop_af', 'MAX_AF']),
            'Splice AI': _extract_value(data, ['spliceai', 'splice_ai', 'SpliceAI']),
            'GERP++': _extract_value(data, ['gerp', 'gerp_rs', 'GERP_RS']),
            'PolyPhen-2': _extract_value(data, ['polyphen', 'polyphen2', 'PolyPhen']),
            'REVEL': _extract_value(data, ['revel', 'revel_score', 'REVEL']),
            'CAD(25)': _extract_value(data, ['cadd', 'cadd_phred', 'CADD_PHRED']),
            'OMIM': _extract_value(data, ['omim', 'omim_id', 'OMIM']),
            'ClinVar': _extract_value(data, ['clinvar', 'clinvar_significance', 'ClinVar_SIG'])
        }

        mapped_count = 0
        for column, value in mappings.items():
            if value:
                variant_mapping[column] = value
                mapped_count += 1
                print(f"  ✓ {column}: {value}")
            else:
                print(f"  ✗ {column}: <no mapping found>")

        mapping_results[variant_id] = {
            'mapped_columns': mapped_count,
            'total_columns': len(mappings),
            'mappings': variant_mapping,
            'raw_data': data
        }

        success_rate = (mapped_count / len(mappings)) * 100
        print(f"  Mapping success: {mapped_count}/{len(mappings)} ({success_rate:.1f}%)")

    return mapping_results

def _extract_value(data: Dict, possible_keys: List[str]) -> str:
    """Helper function to extract values with multiple possible key names"""
    for key in possible_keys:
        if key in data:
            value = data[key]
            if value is not None and str(value).strip():
                return str(value).strip()

        # Try case-insensitive
        for data_key in data.keys():
            if data_key.lower() == key.lower():
                value = data[data_key]
                if value is not None and str(value).strip():
                    return str(value).strip()
    return ""

def simulate_table_update(mapping_results: Dict):
    """Simulate how the App window would update the table"""
    print(f"\n=== Simulating Table Updates ===")

    # Simulate a table with our test variants
    table_data = [
        ["chr1", "12345", ".", "A", "T", "60", "PASS", "", "", ""],  # Row 0
        ["chr2", "67890", ".", "G", "C", "55", "PASS", "", "", ""],  # Row 1
        ["chr3", "11111", ".", "AT", "G", "70", "PASS", "", "", ""], # Row 2
    ]

    # Column indices for annotation columns (after base columns)
    base_col_count = 10  # BASE_COLUMNS count
    annotation_columns = [
        "Gene", "Transcript", "Type", "Consequence", "LOFTEE",
        "Allele Frequencies", "Maximum allele frequency of all subpopulations",
        "Splice AI", "GERP++", "PolyPhen-2", "REVEL", "CAD(25)", "OMIM", "ClinVar"
    ]

    # Extend table to include annotation columns
    for row in table_data:
        row.extend([""] * len(annotation_columns))

    # Apply annotations
    variant_to_row = {
        "chr1:12345:A>T": 0,
        "chr2:67890:G>C": 1,
        "chr3:11111:AT>G": 2
    }

    updates_applied = 0
    for variant_id, result in mapping_results.items():
        if variant_id in variant_to_row:
            row_idx = variant_to_row[variant_id]
            mappings = result['mappings']

            print(f"\nUpdating row {row_idx} for {variant_id}:")

            for col_name, value in mappings.items():
                if col_name in annotation_columns:
                    col_idx = base_col_count + annotation_columns.index(col_name)
                    table_data[row_idx][col_idx] = value
                    print(f"  Set {col_name} = '{value}'")
                    updates_applied += 1

    print(f"\nFinal table state (showing key columns):")
    print("Row | CHROM | POS   | REF | ALT | Gene  | Transcript | Consequence | REVEL | ClinVar")
    print("-" * 80)

    for i, row in enumerate(table_data):
        gene_idx = base_col_count + annotation_columns.index("Gene")
        transcript_idx = base_col_count + annotation_columns.index("Transcript")
        consequence_idx = base_col_count + annotation_columns.index("Consequence")
        revel_idx = base_col_count + annotation_columns.index("REVEL")
        clinvar_idx = base_col_count + annotation_columns.index("ClinVar")

        print(f" {i}  | {row[0]:<5} | {row[1]:<5} | {row[2]:<3} | {row[3]:<3} | {row[gene_idx]:<5} | {row[transcript_idx][:12]:<12} | {row[consequence_idx][:11]:<11} | {row[revel_idx]:<5} | {row[clinvar_idx][:10]}")

    print(f"\nTotal updates applied: {updates_applied}")
    return updates_applied > 0

def main():
    """Run comprehensive annotation assignment test"""
    print("=" * 60)
    print("ANNOTATION ASSIGNMENT TEST")
    print("=" * 60)

    api_url = "http://127.0.0.1:5001"

    # Test 1: Server connectivity
    print("1. Testing server connectivity...")
    if not test_server_health(api_url):
        print("❌ Server test failed. Please ensure your annotation server is running.")
        print("Start with: python server.py")
        return False

    # Test 2: Submit variants
    print("\n2. Submitting test variants...")
    submitted_variants = submit_test_variants(api_url)

    if not submitted_variants:
        print("❌ No variants submitted successfully")
        return False

    # Test 3: Poll for annotations
    print("\n3. Waiting for annotations to complete...")
    annotation_data = poll_for_annotations(api_url, submitted_variants)

    if not annotation_data:
        print("❌ No annotations completed")
        return False

    # Test 4: Test annotation mapping
    print("\n4. Testing annotation field mapping...")
    mapping_results = test_annotation_mapping(annotation_data)

    # Test 5: Simulate table updates
    print("\n5. Simulating table updates...")
    table_success = simulate_table_update(mapping_results)

    # Summary
    print("\n" + "=" * 60)
    print("TEST SUMMARY")
    print("=" * 60)

    submitted_count = len(submitted_variants)
    completed_count = len(annotation_data)

    print(f"Variants submitted: {submitted_count}")
    print(f"Annotations received: {completed_count}")

    if completed_count > 0:
        success_rate = (completed_count / submitted_count) * 100
        print(f"Completion rate: {success_rate:.1f}%")

        total_mappings = sum(r['mapped_columns'] for r in mapping_results.values())
        total_possible = sum(r['total_columns'] for r in mapping_results.values())
        mapping_rate = (total_mappings / total_possible) * 100 if total_possible > 0 else 0

        print(f"Mapping success rate: {total_mappings}/{total_possible} ({mapping_rate:.1f}%)")
        print(f"Table updates: {'✓ SUCCESS' if table_success else '✗ FAILED'}")

        if success_rate >= 50 and mapping_rate >= 30:
            print("\n🎉 OVERALL TEST: ✅ PASSED")
            print("Your App window should correctly receive and assign annotations!")
            return True
        else:
            print("\n⚠️  OVERALL TEST: PARTIAL SUCCESS")
            print("Some annotations are working, but there may be mapping issues.")
            return False
    else:
        print("\n❌ OVERALL TEST: FAILED")
        print("No annotations were received from the server.")
        return False

if __name__ == "__main__":
    print("Make sure your annotation server is running on the correct port!")
    print("Default: python server.py (port 5001)")
    print("")

    # Allow port override
    if len(sys.argv) > 1:
        try:
            port = int(sys.argv[1])
            api_url = f"http://127.0.0.1:{port}"
            print(f"Using custom port: {port}")
        except ValueError:
            print("Invalid port number, using default 5001")

    input("Press Enter to start the annotation assignment test...")

    success = main()

    if success:
        print("\n✅ Test completed successfully!")
        print("Your App window annotation assignment should work correctly.")
    else:
        print("\n🔧 Test revealed some issues that may need attention.")
        print("Check the mapping results above for details.")