"""
Minimal test for retry functionality.

"""
import requests
import time


def test_retry_logic():
    """Test retry functionality with a variant that will fail."""
    base_url = "http://localhost:5001"

    # Use a test variant
    test_variant = {
        'chrom': 'chr1',
        'pos': 99999999,
        'ref': 'GCG',
        'alt': 'C'
    }
    variant_id = f"{test_variant['chrom']}:{test_variant['pos']}:{test_variant['ref']}>{test_variant['alt']}"

    print(f"Testing retry logic with: {variant_id}")

    # Submit variant multiple times to trigger retry limit
    for attempt in range(7):
        print(f"\nAttempt {attempt + 1}:")

        # Submit variant
        submit_response = requests.post(f"{base_url}/submit", json=test_variant)
        submit_result = submit_response.json()
        print(f"  Submit: {submit_result['status']} - {submit_result.get('message', '')}")

        # Wait a bit for processing
        time.sleep(3)

        # Poll for result
        poll_response = requests.get(f"{base_url}/poll/{variant_id}")
        poll_result = poll_response.json()
        print(f"  Poll: {poll_result['status']} - {poll_result.get('message', '')}")

        # Show retry info if available
        if 'retry_info' in poll_result:
            retry = poll_result['retry_info']
            print(f"  Retries: {retry['current_retries']}/{retry['max_retries']}")

        # Stop if exceeded retry limit
        if submit_result['status'] == 'failure' and 'retry limit' in submit_result.get('message', ''):
            print("  Retry limit reached!")
            break


if __name__ == "__main__":
    test_retry_logic()