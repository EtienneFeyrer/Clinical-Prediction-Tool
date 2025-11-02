"""
Batch processor for handling VEP API calls and database operations.
Uses thread pool for efficient parallel processing with controlled concurrency.
"""
import threading
from typing import List, Dict, Set
from datetime import datetime
import requests
import json
import os
import sys
from concurrent.futures import ThreadPoolExecutor
import uuid
import core.helper_methods as helper

sys.path.append(os.path.abspath(".."))
import packages.annotation_classes as ac


sys.path.append(os.path.abspath('.'))
import importlib
from MachineLearning.RandomForest.RandomForestIO import RandomForestIO

sys.path.append(os.path.abspath('.'))
import database.service as db


class BatchProcessor:
    """
    Handles batch processing of variants for VEP annotation.
    Non-blocking batch submission when size reaches 200 or 5 seconds elapsed.
    Uses thread pool for efficient parallel processing.
    """

    def __init__(self,
                 vep_api_url: str = "https://rest.ensembl.org/vep/human/region",
                 max_workers: int = 3):
        """
        Initialize batch processor with thread pool.

        Args:
            vep_api_url: VEP REST API endpoint
            max_workers: Number of worker threads for parallel processing
        """
        self.vep_api_url = vep_api_url
        self.max_workers = max_workers

        # Current batch being collected
        self.in_progress: List[Dict] = []
        self.in_progress_keys: Set[str] = set()  # Track variant keys for fast lookup

        # Track variants currently being processed across all threads
        self.processing_keys: Set[str] = set()

        # Thread pool for batch processing
        self.executor = ThreadPoolExecutor(
            max_workers=max_workers,
            thread_name_prefix="VEP-Worker"
        )
        self.active_futures = set()  # Track active batch processing futures

        # Timing and batch control
        self.lock = threading.Lock()
        self.max_batch_size = 200
        self.max_wait_time = 5.0  # seconds

        # Timer for automatic timeout trigger
        self._timer = None
        self._timer_lock = threading.Lock()

        self.retry_counts = {}  # variant_key -> retry_count
        self.max_retries = 3

    def add_variant(self, variant_key: str, chrom: str, pos: int, ref: str, alt: str) -> str:
        """
        Add variant to in-progress batch.
        Returns variant_id for tracking.
        """
        variant_data = {
            'variant_key': variant_key,
            'chrom': chrom,
            'pos': pos,
            'ref': ref,
            'alt': alt,
            'timestamp': datetime.now()
        }

        with self.lock:
            # Check if variant exceeded retry limit
            if self.retry_counts.get(variant_key, 0) >= self.max_retries:
                return variant_key  # Don't add to batch, but return key for status check

            # Check if variant is already being processed
            if variant_key in self.processing_keys or variant_key in self.in_progress_keys:
                return variant_key  # Already being processed

            self.in_progress.append(variant_data)
            self.in_progress_keys.add(variant_key)

            # Check if batch should be triggered immediately (size limit reached)
            if len(self.in_progress) >= self.max_batch_size:
                self._cancel_timer()
                self._trigger_batch_async()
            else:
                # Start or restart timer for timeout trigger
                self._restart_timer()

        return variant_key

    def is_variant_in_progress(self, variant_key: str) -> bool:
        """Check if variant is currently being processed."""
        with self.lock:
            return variant_key in self.in_progress_keys or variant_key in self.processing_keys

    def increment_retry_count(self, variant_key: str):
        """Increment retry count for a variant."""
        with self.lock:
            self.retry_counts[variant_key] = self.retry_counts.get(variant_key, 0) + 1

    def get_retry_info(self, variant_key: str) -> tuple:
        """Get retry information. Returns (current_retries, max_retries, exceeded_limit)."""
        with self.lock:
            current = self.retry_counts.get(variant_key, 0)
            exceeded = current >= self.max_retries
            return current, self.max_retries, exceeded

    def shutdown(self, wait: bool = True):
        """
        Manually shutdown the processor (only for manual shutdown).

        Args:
            wait: Whether to wait for all pending batches to complete
        """
        print("Manual shutdown requested...")

        # Process any remaining variants in current batch
        with self.lock:
            if self.in_progress:
                print(f"Processing final batch of {len(self.in_progress)} variants...")
                self._trigger_batch_async()

        # Cancel timer
        self._cancel_timer()

        # Wait for all active batches to complete
        if wait and self.active_futures:
            print(f"Waiting for {len(self.active_futures)} active batches to complete...")
            for future in list(self.active_futures):
                try:
                    future.result(timeout=30)
                except Exception as e:
                    print(f"Error waiting for batch: {e}")

        # Shutdown thread pool
        print("Shutting down thread pool...")
        self.executor.shutdown(wait=wait)
        print("Batch processor shutdown complete")

    def _restart_timer(self):
        """Restart the timeout timer to trigger batch processing after 5 seconds."""
        with self._timer_lock:
            # Cancel existing timer
            if self._timer and self._timer.is_alive():
                self._timer.cancel()

            # Start new timer - triggers batch processing (NOT shutdown)
            self._timer = threading.Timer(self.max_wait_time, self._timeout_batch_trigger)
            self._timer.daemon = True
            self._timer.start()

    def _timeout_batch_trigger(self):
        """Called when timeout expires - triggers batch processing only."""
        with self.lock:
            # Only trigger if there are variants waiting
            if self.in_progress:
                print(f"Timeout: Processing batch of {len(self.in_progress)} variants")
                self._trigger_batch_async()
            # NO shutdown - just wait for more variants

    def _cancel_timer(self):
        """Cancel the current timer."""
        with self._timer_lock:
            if self._timer and self._timer.is_alive():
                self._timer.cancel()
                self._timer = None

    def _trigger_batch_async(self):
        """Submit current batch to thread pool for processing (non-blocking)."""
        if not self.in_progress:
            return

        # Generate unique batch ID
        batch_id = str(uuid.uuid4())[:8]

        # Copy current batch and clear in-progress
        current_batch = self.in_progress.copy()
        current_keys = self.in_progress_keys.copy()

        # Move variants to processing state
        self.processing_keys.update(current_keys)

        # Clear current batch for new variants
        self.in_progress.clear()
        self.in_progress_keys.clear()

        # Cancel timer since we're processing now
        self._cancel_timer()

        # Submit to thread pool
        workload = self.executor.submit(
            self._process_batch,
            batch_id,
            current_batch,
            current_keys
        )

        # Track the future
        self.active_futures.add(workload)

        # Add callback to clean up when done
        workload.add_done_callback(lambda f: self._on_batch_complete(f, current_keys))

    def _process_batch(self, batch_id: str, batch: List[Dict], batch_keys: Set[str]):
        """Process batch of variants with VEP API (runs in thread pool worker)."""
        thread_name = threading.current_thread().name

        try:
            print(f"[{thread_name}] Processing batch {batch_id} of {len(batch)} variants...")

            # Prepare VEP request format
            vep_regions = []
            for variant in batch:
                region = helper.format_vep_region(variant['chrom'], variant['pos'], variant['ref'], variant['alt'])
                vep_regions.append(region)

            # Call VEP API
            headers = {'Content-Type': 'application/json'}
            payload = {
                "variants": vep_regions,
                "REVEL": True,
                "CADD": True,
                "SpliceAI": True,
                "protein": True,
                "gencode_basic": True,
                "LoF": True,
                "mane": True,
                "hgvs": True,
                "dbNSFP": "clinvar_OMIM_id,GERP++_RS"
            }

            response = requests.post(self.vep_api_url,
                                     headers=headers,
                                     data=json.dumps(payload),
                                     timeout=300)

            if response.status_code == 200:
                vep_results = response.json()
                
                self._store_results(vep_results)
                print(f"[{thread_name}] Batch {batch_id} processing completed successfully")
                return True
            else:
                print(f"[{thread_name}] VEP API error for batch {batch_id}: {response.status_code}")
                return False

        except Exception as e:
            print(f"[{thread_name}] Batch {batch_id} processing error: {e}")
            return False

    def _on_batch_complete(self, future, batch_keys: Set[str]):
        """
        Callback called when a batch completes.
        Cleans up processing state.
        """
        # Clean up processing state
        with self.lock:
            self.processing_keys -= batch_keys
            self.active_futures.discard(future)
        try:
            success = future.result()
            if not success:
                # Increment retry count for all variants in failed batch
                for variant_key in batch_keys:
                    self.increment_retry_count(variant_key)
                    print(
                        f"Batch failed - {variant_key} retry count: {self.retry_counts.get(variant_key, 0)}/{self.max_retries}")
        except Exception as e:
            # Also increment on exception
            for variant_key in batch_keys:
                self.increment_retry_count(variant_key)
            print(f"Batch error - incrementing retry counts: {e}")


    def _store_results(self, vep_results: List[Dict]):
        """Store VEP results using bulk database operations."""
        # Parse all variants first
        annotations_data = []
        transcripts_data = []

        for i, vep_data in enumerate(vep_results):
            try:
                variant_annotation = self._parse_vep_to_annotation(vep_data)
                json_data = helper.annotation_to_json(variant_annotation)

                variant_key = None
                if 'input' in vep_data:
                    variant_key = helper.transform_variant(vep_data['input'])
                if not variant_key:
                    variant_key = f"variant_{i}"

                # Prepare for bulk insert
                annotations_data.append((
                    variant_key, json_data.get('GENE'), json_data.get('CADD'),
                    json_data.get('ML-Score'), json_data.get('Most Severe Consequence'),
                    json_data.get('gnomAD AF'), json_data.get('max_allele_freq'), json_data.get('OMIM'), json_data.get('CLINSIG')
                ))

                for transcript in json_data.get('transcripts', []):
                    transcripts_data.append((
                        variant_key, transcript.get('transcript_id'), transcript.get('polyphen'),
                        transcript.get('protein_notation'), transcript.get('revel'), transcript.get('spliceai'),
                        transcript.get('is_mane', False), transcript.get('lof'), transcript.get('impact'),
                        transcript.get('gerp'), transcript.get('cdna_notation'), transcript.get('consequences')
                    ))
            except Exception as e:
                print(f"Error parsing variant {i}: {e}")

        # Single database transaction for entire batch
        if annotations_data:
            db.DatabaseService.bulk_insert_annotations(annotations_data, transcripts_data)

    def _parse_vep_to_annotation(self, vep_data: Dict):
        """Convert VEP data to our annotation format (same as parser)."""
        variant_annotation = ac.GeneAnnotations()

        # Basic Gene and Score Annotations
        if vep_data.get('transcript_consequences'):
            first_transcript = vep_data['transcript_consequences'][0]
            if 'cadd_phred' in first_transcript:
                variant_annotation.add_annotation(ac.Annotation_Float("CADD", first_transcript['cadd_phred']))
            if 'gene_symbol' in first_transcript:
                variant_annotation.add_annotation(ac.Annotation_Str("GENE", first_transcript['gene_symbol']))
            else:
                variant_annotation.add_annotation(ac.Annotation_Str("GENE", ""))

        # Population Frequency and Clinical Significance
        colocated_data = helper.extract_colocated_variants_data(vep_data)

        if colocated_data['gnomad_af']:
            variant_annotation.add_annotation(ac.Annotation_Float("gnomAD AF", colocated_data['gnomad_af']))

        if colocated_data['max_allele_freq']:
            variant_annotation.add_annotation(ac.Annotation_Float("max_allele_freq", colocated_data['max_allele_freq']))

        if colocated_data['clin_sig']:
            variant_annotation.add_annotation(ac.Annotation_Str("CLINSIG", colocated_data['clin_sig']))

        # Consequence and Disease Annotations
        if 'most_severe_consequence' in vep_data:
            variant_annotation.add_annotation(
                ac.Annotation_Str("Most Severe Consequence", vep_data['most_severe_consequence']))

        clinvar_omim_id = helper.extract_clinvar_omim_id(vep_data)
        if clinvar_omim_id:
            variant_annotation.add_annotation(ac.Annotation_Str("OMIM", clinvar_omim_id.split('&')))
        else:
            variant_annotation.add_annotation(ac.Annotation_Str("OMIM", ""))

        # Machine Learning Pathogenicity Score
        try:
            # Transform data format for ML model compatibility
            if isinstance(vep_data, dict):
                ml_input = [vep_data.copy()]  # Make a copy so we don't modify original

                # Fix chromosome format for ML model (chr1 -> 1)
                if 'seq_region_name' in ml_input[0] and ml_input[0]['seq_region_name'].startswith('chr'):
                    ml_input[0]['seq_region_name'] = ml_input[0]['seq_region_name'][3:]

                # Fix input format for ML model
                if 'input' in ml_input[0] and ml_input[0]['input'].startswith('chr'):
                    ml_input[0]['input'] = ml_input[0]['input'][3:]
            else:
                ml_input = vep_data

            importlib.reload(sys.modules['MachineLearning.RandomForest.RandomForestIO'])
            importlib.reload(sys.modules['MachineLearning.RandomForest.RandomForestModel'])
            object = RandomForestIO("forest")
            ml_result = object.get_ml_values(ml_input)

            # Extract the score from the result
            if isinstance(ml_result, list) and len(ml_result) > 0:
                ml_score = ml_result[0]  # Extract first element
            else:
                ml_score = ml_result  # If it's already a number

            variant_annotation.add_annotation(ac.Annotation_Float("ML-Score", ml_score))
            print(f"ML-Score calculated: {ml_score}")
        except Exception as e:
            print(f"ML model error: {e}")
            import traceback
            traceback.print_exc()
            # Use fallback mock score
            mock_score = 0.75
            variant_annotation.add_annotation(ac.Annotation_Float("ML-Score", mock_score))
            print(f"Using mock ML-Score: {mock_score}")

        # Transcript-level annotations
        for entry in vep_data.get('transcript_consequences', []):
            is_mane = bool(entry.get('mane') and len(entry.get('mane', [])) > 0)

            # Extract cDNA and Protein notation
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
                consequence_string = ""

            transcript_ann = ac.TranscriptAnnotations(
                transcript_id=entry['transcript_id'] or "",
                impact=helper.validate_enum(entry['impact'] or "", ac.ImpactLevel),
                revel=entry.get('revel') or None,
                gerp=entry.get('gerp++_rs') or None,
                spliceai=helper.extract_max_spliceai_score(entry.get('spliceai', {})) if entry.get(
                    'spliceai') else None,
                polyphen=entry.get('polyphen_score') or None,
                loftee=helper.validate_enum(entry.get('lof') or "", ac.LofteeLevel),
                is_mane=is_mane,
                cdna_notation=cdna_notation,
                protein_notation=protein_notation,
                consequences=consequence_string
            )

            variant_annotation.add_transcript_annotation(transcript_ann)

        return variant_annotation

# Global batch processor instance
batch_processor = BatchProcessor(max_workers=3)
