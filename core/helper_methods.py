"""
Helper methods for NGS variant annotation processing.
Contains utility functions for VEP data extraction, frequency calculations,
and annotation format conversions.
"""

import enum
import sys
import os

sys.path.append(os.path.abspath(".."))
import packages.annotation_classes as ac


def format_vep_region(chrom: str, pos: int, ref: str, alt: str) -> str:
    """Format variant for VEP API with correct position calculation."""
    if len(ref) == len(alt) == 1:
        # SNV: same start and end position
        return f"{chrom} {pos} {pos} {ref}/{alt} +"

    elif len(ref) > len(alt):
        # Deletion: end position is start + length of reference - 1
        end_pos = pos + len(ref) - 1
        return f"{chrom} {pos} {end_pos} {ref}/{alt} +"

    elif len(alt) > len(ref):
        # Insertion: VEP expects same start/end position for simple insertions, under the assumption that there is only one anchor base
        return f"{chrom} {pos} {pos} {ref}/{alt} +"

    else:
        # Same length substitution
        end_pos = pos + len(ref) - 1
        return f"{chrom} {pos} {end_pos} {ref}/{alt} +"


def extract_max_spliceai_score(spliceai_data: dict) -> float:
    """Extract maximum delta score from SpliceAI data."""
    if not spliceai_data:
        return None
    try:
        delta_scores = []
        for key in ['DS_AG', 'DS_AL', 'DS_DG', 'DS_DL']:
            if key in spliceai_data and spliceai_data[key] is not None:
                delta_scores.append(abs(float(spliceai_data[key])))
        return max(delta_scores) if delta_scores else None
    except (ValueError, TypeError):
        return None


def extract_clinvar_omim_id(vep_data: dict) -> str:
    """Extract ClinVar OMIM ID from VEP data."""
    for transcript in vep_data.get('transcript_consequences', []):
        clinvar_omim = transcript.get('clinvar_omim_id')
        if clinvar_omim:
            return str(clinvar_omim)
    return ''


def extract_colocated_variants_data(vep_data: dict) -> dict:
    """Extract gnomAD frequency and Clinical significance from VEP colocated variants data."""
    result = {
        'gnomad_af': None,
        'max_allele_freq': None,
        'clin_sig': None
    }

    if 'colocated_variants' not in vep_data:
        return result

    found_freq = False
    found_clin = False

    for colocated in vep_data['colocated_variants']:
        # Look for frequencies
        if not found_freq and 'frequencies' in colocated:
            for allele, freq_data in colocated['frequencies'].items():
                gnomad_af = freq_data.get('gnomadg', freq_data.get('af'))
                if gnomad_af:
                    result['gnomad_af'] = gnomad_af
                result['max_allele_freq'] = get_max_frequency(freq_data)
                found_freq = True
                break

        # Look for clin_sig
        if not found_clin and 'clin_sig' in colocated:
            result['clin_sig'] = colocated['clin_sig']
            found_clin = True

        # Exit early if both found
        if found_freq and found_clin:
            break

    return result


def get_max_frequency(freq_data: dict) -> float:
    """Get the maximum frequency from all populations."""
    all_frequencies = []

    # Get all frequency values from the dictionary
    for key, value in freq_data.items():
        if value is not None and isinstance(value, (int, float)):
            all_frequencies.append(float(value))

    # Return the highest frequency, or None if no frequencies found
    return max(all_frequencies) if all_frequencies else None


def annotation_to_json(annotation_data):
    """Convert annotation to JSON format for database storage."""
    d = {}

    # Add variant-level annotations
    for entry in annotation_data.get_all_annotations():
        # Use the specialized method for string annotations
        if isinstance(entry, ac.Annotation_Str):
            d[entry.origin] = entry.get_as_string()
        else:
            d[entry.origin] = entry.data

    # Add transcript-level annotations
    transcripts = []
    for transcript in annotation_data.get_transcript_annotations():
        transcript_dict = {
            'transcript_id': transcript.transcript_id,
            'impact': transcript.impact,
            'revel': transcript._revel,
            'gerp': transcript._gerp,
            'spliceai': transcript._spliceai,
            'polyphen': transcript._polyphen,
            'lof': transcript._loftee,
            'is_mane': transcript._is_mane,
            'cdna_notation': transcript._cdna_notation,
            'protein_notation': transcript._protein_notation,
            'consequences': transcript._consequences
        }
        transcripts.append(transcript_dict)

    d['transcripts'] = transcripts
    return d


def transform_variant(variant: str) -> str:
    """
    Transform variant string from 'chr2 148483494 148483494 C/A +'
    into 'chr2:148483494:C>A'
    """
    chrom, start, end, alleles, strand = variant.strip().split()
    ref, alt = alleles.split("/")
    return f"{chrom}:{start}:{ref}>{alt}"


def validate_enum(value: str, enum_class: type[enum.Enum]) -> str:
    """
    General enum validator - works for any enum class.
    Returns validated string or empty string if invalid.
    """
    if not value:
        return ""
    try:
        return enum_class(value.upper()).value
    except (ValueError, AttributeError):
        return ""



def create_cdna_notation(transcript_data):
    """Create cDNA HGVS notation"""
    if not transcript_data.get('cdna_start') or not transcript_data.get('codons'):
        return ""

    position = transcript_data['cdna_start']
    codons = transcript_data['codons']

    # Parse codons like "gCt/gTt" to get the changed base
    if '/' in codons:
        ref_codon, alt_codon = codons.split('/')
        # Find the different base
        for i, (ref_base, alt_base) in enumerate(zip(ref_codon.upper(), alt_codon.upper())):
            if ref_base != alt_base:
                return f"c.{position + i}{ref_base}>{alt_base}"

    return ""


def create_protein_notation(transcript_data):
    """Create protein HGVS notation"""
    if not transcript_data.get('protein_start') or not transcript_data.get('amino_acids'):
        return ""

    position = transcript_data['protein_start']
    amino_acids = transcript_data['amino_acids']

    # Parse amino acids like "A/V"
    if '/' in amino_acids:
        ref_aa, alt_aa = amino_acids.split('/')
        # Convert single letter to 3-letter code
        aa_dict = {
            'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys',
            'Q': 'Gln', 'E': 'Glu', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
            'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro',
            'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val'
        }

        ref_long = aa_dict.get(ref_aa, ref_aa)
        alt_long = aa_dict.get(alt_aa, alt_aa)

        return f"p.{ref_long}{position}{alt_long}"

    return ""
