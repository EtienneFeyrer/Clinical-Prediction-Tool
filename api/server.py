"""
Flask API endpoints for polling-based annotation system.
In-progress tracking in local memory only.
"""
from flask import Flask, request, jsonify
import os
import sys
import signal
import re


sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from core.batch_processor import batch_processor
import database.service as db

app = Flask(__name__)


def shutdown_handler(signum, frame):
    """Handle Ctrl+C shutdown"""
    print("\nShutting down...")
    batch_processor.shutdown(wait=True)
    sys.exit(0)


@app.route('/submit', methods=['POST'])
def submit_variant():
    """
    Submit variant for annotation.
    Returns success + variant_id if new, failure if already exists.

    In-progress tracking in local memory only.
    """
    try:
        data = request.get_json()

        # Validate input
        required_fields = ['chrom', 'pos', 'ref', 'alt']
        if not all(field in data for field in required_fields):
            return jsonify({
                'status': 'error',
                'message': 'Missing required fields: chrom, pos, ref, alt'
            }), 400

        chrom = data['chrom']
        pos = int(data['pos'])
        ref = data['ref']
        alt = data['alt']

        chrom_re = re.compile(r'^chr(?:[1-9]|1[0-9]|2[0-3]|X|Y|M|MT)$', re.IGNORECASE)
        base_re = re.compile(r'^[ACGTacgt]+$')

        if not chrom_re.match(chrom):
            return jsonify({
                'status': 'error',
                'message': f"Invalid chromosome format: '{chrom}'for a human."
            }), 400

        if not base_re.match(ref):
            return jsonify({
                'status': 'error',
                'message': f"Invalid REF allele: '{ref}'"
            }), 400

        if not base_re.match(alt):
            return jsonify({
                'status': 'error',
                'message': f"Invalid ALT allele: '{alt}'"
            }), 400

        # Create variant key
        variant_key = f"{chrom}:{pos}:{ref}>{alt}"

        # Check if variant already exists in database
        if db.DatabaseService.variant_exists(variant_key):
            return jsonify({
                'status': 'failure',
                'message': 'Variant already annotated',
                'variant_id': variant_key
            }), 200

            # Check retry status
        current_retries, max_retries, exceeded_limit = batch_processor.get_retry_info(variant_key)

        if exceeded_limit:
            return jsonify({
                'status': 'failure',
                'message': f'Variant exceeded retry limit ({current_retries}/{max_retries} attempts)',
                'variant_id': variant_key,
                'retry_info': {
                    'current_retries': current_retries,
                    'max_retries': max_retries,
                    'exceeded_limit': True
                }
            }), 200

        # Check if currently processing
        if batch_processor.is_variant_in_progress(variant_key):
            return jsonify({
                'status': 'success',
                'message': 'Variant already in progress',
                'variant_id': variant_key,
                'retry_info': {
                    'current_retries': current_retries,
                    'max_retries': max_retries,
                    'exceeded_limit': False
                }
            }), 200

        # Add to batch processor (non-blocking)
        variant_id = batch_processor.add_variant(variant_key, chrom, pos, ref, alt)

        return jsonify({
            'status': 'success',
            'message': f'Variant added to queue (attempt {current_retries + 1}/{max_retries})',
            'variant_id': variant_id,
            'retry_info': {
                'current_retries': current_retries,
                'max_retries': max_retries,
                'exceeded_limit': False
            }
        }), 200

    except ValueError as e:
        return jsonify({
            'status': 'error',
            'message': f'Invalid input: {str(e)}'
        }), 400
    except Exception as e:
        return jsonify({
            'status': 'error',
            'message': f'Server error: {str(e)}'
        }), 500


@app.route('/poll/<variant_id>', methods=['GET'])
def poll_variant(variant_id: str):
    """
    Poll for variant annotation status.
    Returns annotation data if ready, or 'processing' if still in progress.
    """
    try:
        # Check if annotation is ready in database
        annotation_data = db.DatabaseService.get_variant_annotation(variant_id)

        if annotation_data:
            return jsonify({
                'status': 'completed',
                'source': 'cache',
                'variant_id': variant_id,
                'annotation': annotation_data
            }), 200

        # Check if variant is currently being processed (memory check)
        if batch_processor.is_variant_in_progress(variant_id):
            return jsonify({
                'status': 'processing',
                'message': 'Annotation in progress',
                'variant_id': variant_id
            }), 202

        # Check retry status
        current_retries, max_retries, exceeded_limit = batch_processor.get_retry_info(variant_id)

        if exceeded_limit:
            return jsonify({
                'status': 'failed',
                'message': f'Annotation failed after {max_retries} attempts',
                'variant_id': variant_id,
                'retry_info': {
                    'current_retries': current_retries,
                    'max_retries': max_retries,
                    'exceeded_limit': True
                }
            }), 200

        elif current_retries > 0:
            return jsonify({
                'status': 'retry_available',
                'message': f'Previous attempts failed ({current_retries}/{max_retries}). Can retry.',
                'variant_id': variant_id,
                'retry_info': {
                    'current_retries': current_retries,
                    'max_retries': max_retries,
                    'exceeded_limit': False
                }
            }), 404

        else:
            return jsonify({
                'status': 'not_found',
                'message': 'Variant not found. Please submit first.',
                'variant_id': variant_id
            }), 404

    except Exception as e:
        return jsonify({
            'status': 'error',
            'message': f'Server error: {str(e)}'
        }), 500


@app.route('/health', methods=['GET'])
def health_check():
    """Health check endpoint."""
    return jsonify({
        'status': 'healthy',
        'service': 'annotation-api',
        'in_progress_count': len(batch_processor.in_progress)
    }), 200


@app.route('/statistics', methods=['GET'])
def get_statistics():
    """Get service statistics."""
    return jsonify({
        'in_progress_count': len(batch_processor.in_progress),
        'batch_size_limit': batch_processor.max_batch_size,
        'batch_time_limit': batch_processor.max_wait_time,
        'in_progress_variants': list(batch_processor.in_progress_keys)
    }), 200


if __name__ == '__main__':
    signal.signal(signal.SIGINT, shutdown_handler)

    host = os.getenv('API_HOST', '0.0.0.0')
    port = int(os.getenv('API_PORT', '5001'))
    debug = os.getenv('FLASK_DEBUG', 'false').lower() == 'true'
    app.run(debug=debug, host=host, port=port)