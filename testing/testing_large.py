"""
Simple test setup for annotation service with VCF file support

This tests the annotation classes, batch processor, and API
WITH DEBUGGING TO FIND MISSING VARIANTS
"""
import requests
import time

import gzip


def parse_vcf(vcf_path):
    """Parse VCF and return list of variants."""
    variants = []
    opener = gzip.open if vcf_path.endswith('.gz') else open

    with opener(vcf_path, "rt") as f:
        for line_num, line in enumerate(f, 1):
            if line.startswith("#"):
                continue

            cols = line.strip().split("\t")

            # Check if line has enough columns
            if len(cols) < 5:
                print(f"   Skipping malformed line {line_num}: {line.strip()}")
                continue

            try:
                chrom, pos, ref, alt = cols[0], int(cols[1]), cols[3], cols[4]
                for a in alt.split(","):
                    variants.append({
                        'chrom': chrom,
                        'pos': pos,
                        'ref': ref,
                        'alt': a
                    })
            except (ValueError, IndexError) as e:
                print(f"   Skipping invalid line {line_num}: {e}")
                continue

    return variants


def test_api_basic():
    """Test basic API functionality"""
    base_url = "http://localhost:5001"

    print("=== Testing Annotation Service WITH DEBUGGING ===\n")

    # Test 1: Health Check
    print("1. Testing health endpoint...")
    try:
        response = requests.get(f"{base_url}/health")
        print(f"   Health: {response.json()}")
    except Exception as e:
        print(f"   Error: {e}")
        return

    # Test 2: Load variants from VCF and submit them
    print("\n2. Loading variants from VCF...")
    vcf_file = "rare_coding_variants.vcf"
    variants = parse_vcf(vcf_file)
    print(f"   Found {len(variants)} variants")

    # Limit variants for testing (you can change this number)
    test_variants = variants
    print(f"   Testing first {len(test_variants)} variants")

    print("\n3. Submitting all variants...")
    submitted_variants = []
    submit_failed = 0
    variant_id_mapping = {}  # Track what we sent vs what server returned

    # Submit all variants first
    for i, variant_data in enumerate(test_variants):
        # Create variant_id the same way we always do
        client_variant_id = f"{variant_data['chrom']}:{variant_data['pos']}:{variant_data['ref']}>{variant_data['alt']}"

        try:
            response = requests.post(f"{base_url}/submit", json=variant_data)
            submit_result = response.json()

            print(f"   SERVER RESPONSE: {submit_result}")

            # Check what variant_id the server returned
            server_variant_id = submit_result.get('variant_id')
            print(f"   SERVER RETURNED ID: {server_variant_id}")

            # Check if IDs match
            if server_variant_id != client_variant_id:
                print(f"   ⚠️  ID MISMATCH! Client: {client_variant_id}, Server: {server_variant_id}")
                variant_id_mapping[client_variant_id] = server_variant_id

            if submit_result.get('status') in ['success', 'failure']:  # Both success and "already annotated" are OK
                submitted_variants.append(client_variant_id)  # Use our client ID for now
                status_msg = "success" if submit_result.get('status') == 'success' else "already_exists"
                print(f"   [{i+1:2d}/{len(test_variants)}] {client_variant_id} {status_msg}")
            else:
                print(f"   [{i+1:2d}/{len(test_variants)}] {client_variant_id} ERROR: {submit_result}")
                submit_failed += 1

        except Exception as e:
            print(f"   [{i+1:2d}/{len(test_variants)}] {client_variant_id} EXCEPTION: {e}")
            submit_failed += 1

    print(f"\n   === SUBMISSION SUMMARY ===")
    print(f"   Total variants: {len(test_variants)}")
    print(f"   Successfully submitted: {len(submitted_variants)}")
    print(f"   Failed submissions: {submit_failed}")
    print(f"   ID mismatches found: {len(variant_id_mapping)}")

    if variant_id_mapping:
        print(f"   ID Mapping issues:")
        for client_id, server_id in variant_id_mapping.items():
            print(f"     {client_id} → {server_id}")

    if not submitted_variants:
        print("No variants were submitted successfully")
        return

    print(f"\n4. Polling database completion status...")

    # Poll to see how many are completed in database
    max_polls = 60  # Poll for up to 1 minute
    missing_variants = []
    # Initialize counters so they're available outside the loop
    completed_count = 0
    processing_count = 0
    not_found_count = 0
    retry_count = 0
    failed_count = 0

    for poll_round in range(max_polls):
        # Reset counters for this round
        completed_count = 0
        processing_count = 0
        not_found_count = 0
        retry_count = 0
        failed_count = 0
        checked_variants = []
        poll_results = {}

        print(f"\n Poll round {poll_round+1}, checking {len(submitted_variants)} variants")

        for variant_id in submitted_variants:
            checked_variants.append(variant_id)

            # Use server ID if we have a mapping, otherwise use client ID
            poll_id = variant_id_mapping.get(variant_id, variant_id)

            try:
                poll_response = requests.get(f"{base_url}/poll/{poll_id}")
                poll_result = poll_response.json()

                poll_results[variant_id] = poll_result

                if poll_result['status'] == 'completed':
                    completed_count += 1
                elif poll_result['status'] == 'processing':
                    processing_count += 1
                elif poll_result['status'] == 'not_found':
                    not_found_count += 1
                elif poll_result['status'] == 'retry_available':
                    retry_count += 1
                    print(f"   RETRY: {variant_id} is available for retry")
                elif poll_result['status'] == 'failed':
                    failed_count += 1
                    print(f"   FAILED: {variant_id} failed processing")
                else:
                    print(f"   UNKNOWN STATUS for {variant_id}: {poll_result['status']}")

            except Exception as e:
                not_found_count += 1
                poll_results[variant_id] = f"EXCEPTION: {e}"
                print(f"   POLLING EXCEPTION for {variant_id}: {e}")

        # Show progress
        total_variants = len(submitted_variants)
        accounted_for = completed_count + processing_count + not_found_count + retry_count + failed_count
        missing_count = total_variants - accounted_for

        progress_bar = "█" * (completed_count * 20 // total_variants) + "▒" * (20 - (completed_count * 20 // total_variants))

        print(f"\r   Poll {poll_round+1:2d}: [{progress_bar}] {completed_count:2d}/{total_variants} completed, {processing_count:2d} processing, {retry_count:2d} retry, {failed_count:2d} failed, {not_found_count:2d} not found", end="", flush=True)

        # DEBUG: Check for missing variants
        if missing_count > 0:
            print(f"\n MISSING VARIANTS: {missing_count}")
            print(f"   Total submitted: {total_variants}")
            print(f"   Accounted for: completed({completed_count}) + processing({processing_count}) + retry({retry_count}) + failed({failed_count}) + not_found({not_found_count}) = {accounted_for}")

            # Find which variants are missing
            missing_this_round = []
            for variant_id in submitted_variants:
                if variant_id not in poll_results:
                    missing_this_round.append(variant_id)
                elif isinstance(poll_results[variant_id], str) and "EXCEPTION:" in poll_results[variant_id]:
                    missing_this_round.append(variant_id)

            print(f"   Missing variants this round: {missing_this_round}")

        # If all are completed or failed (final states), break
        if (completed_count + failed_count) == total_variants:
            print(f"\n All {total_variants} variants processed! ({completed_count} completed, {failed_count} failed)")
            break

        # If nothing is processing anymore, break
        if processing_count == 0 and retry_count == 0 and poll_round > 10:  # Give at least 10 seconds
            print(f"\n  No more variants processing. Final: {completed_count}/{total_variants} completed, {failed_count} failed")

            # Find the missing variants
            missing_variants = []
            for variant_id in submitted_variants:
                result = poll_results.get(variant_id)
                if not result or (isinstance(result, str) and 'EXCEPTION' in result):
                    missing_variants.append(variant_id)
                elif result.get('status') not in ['completed', 'processing', 'not_found', 'retry_available', 'failed']:
                    missing_variants.append(variant_id)

            if missing_variants:
                print(f"\n   COMPLETELY MISSING VARIANTS:")
                for missing_id in missing_variants:
                    poll_id = variant_id_mapping.get(missing_id, missing_id)
                    print(f"     - {missing_id} (polling as: {poll_id})")

                    # Try one more direct poll to see what happens
                    try:
                        direct_response = requests.get(f"{base_url}/poll/{poll_id}")
                        direct_result = direct_response.json()
                        print(f"       Direct poll result: {direct_result}")
                    except Exception as e:
                        print(f"       Direct poll failed: {e}")

            break

        time.sleep(1)

    print(f"\n\n=== FINAL RESULTS ===")
    print(f"   Total submitted: {len(submitted_variants)}")
    print(f"   Completed: {completed_count}")
    print(f"   Still processing: {processing_count}")
    print(f"   Retry available: {retry_count}")
    print(f"   Failed: {failed_count}")
    print(f"   Not found: {not_found_count}")
    print(f"   Completely missing: {len(missing_variants)}")

    if len(submitted_variants) > 0:
        # Success rate only counts completed variants as successful
        success_rate = completed_count/len(submitted_variants)*100
        print(f"   Success rate: {success_rate:.1f}% (only completed variants)")

        # Show processing rate (completed + processing + retry)
        processing_rate = (completed_count + processing_count + retry_count)/len(submitted_variants)*100
        print(f"   Processing rate: {processing_rate:.1f}% (completed + processing + retry)")


    # Test 6: Statistics
    print("\n6. Checking final statistics...")
    try:
        response = requests.get(f"{base_url}/statistics")
        stats = response.json()
        print(f"   Statistics: {stats}")
    except Exception as e:
        print(f"   Error: {e}")


if __name__ == "__main__":
    print("Make sure the annotation service is running first!")
    print("Start it with: python server.py or docker compose up")
    print("Then run this test with: python testing_large.py\n")

    # Wait a moment for user to read
    input("Press Enter to start testing...")

    # Run tests
    test_api_basic()