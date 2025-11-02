"""
Simple test setup for annotation service

This tests the annotation classes, batch processor, and API without needing a database.
"""
import requests
import time
import json
import threading


def test_api_basic():
    """Test basic API functionality"""
    base_url = "http://localhost:5001"

    print("=== Testing Annotation Service ===\n")

    # Test 1: Health Check
    print("1. Testing health endpoint...")
    try:
        response = requests.get(f"{base_url}/health")
        print(f"   Health: {response.json()}")
    except Exception as e:
        print(f"   Error: {e}")
        return

    # Test 2: Submit a variant
    print("\n2. Submitting test variant...")
    variant_data = {
        'chrom': 'chr2',
        'pos':  162279995,
        'ref': 'C',
        'alt': 'G'
    }

    try:
        response = requests.post(f"{base_url}/submit", json=variant_data)
        submit_result = response.json()
        print(f"   Submit result: {submit_result}")

        if submit_result.get('status') == 'success':
            variant_id = submit_result['variant_id']

            # Test 3: Poll for result
            print(f"\n3. Polling for variant {variant_id}...")
            for i in range(10):
                time.sleep(2)
                poll_response = requests.get(f"{base_url}/poll/{variant_id}")
                poll_result = poll_response.json()
                print(f"   Poll {i + 1}: {poll_result['status']}")

                if poll_result['status'] == 'completed':
                    print(f"Annotation completed!")
                    print(f"Result: {json.dumps(poll_result['annotation'], indent=2)}")
                    break
                elif poll_result['status'] == 'processing':
                    continue
                else:
                    print(f"   Status: {poll_result}")
                    break

    except Exception as e:
        print(f"   Error: {e}")

    # Test 4: Statistics
    print("\n4. Checking statistics...")
    try:
        response = requests.get(f"{base_url}/statistics")
        stats = response.json()
        print(f"   Statistics: {stats}")
    except Exception as e:
        print(f"   Error: {e}")


if __name__ == "__main__":
    import sys
    
    print("Make sure the annotation service is running first!")
    print("Start it with: python app.py or python -m api.server")
    print("Then run this test with: python test_service.py\n")

    # Only wait for input if running interactively (not in CI)
    if sys.stdin.isatty():
        input("Press Enter to start testing...")
    else:
        print("Running in non-interactive mode (CI), starting tests automatically...")

    # Run tests
    test_api_basic()
