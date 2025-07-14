import os
import uuid
import time
import atexit
import logging
from io import StringIO
from urllib.parse import quote # Still useful for other special characters
from typing import Dict, List, Any

from flask import Flask, request, jsonify, send_from_directory
from flask_cors import CORS
from apscheduler.schedulers.background import BackgroundScheduler

from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, CompoundLocation

# Configure logging for better feedback
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s: %(message)s')

# --- Configuration ---
TEMP_DIR = 'temp_files'
MAX_CONTENT_LENGTH = 20 * 1024 * 1024  # 20 MB, adjust as needed for large GenBank files
FILE_LIFETIME_SECONDS = 300  # 5 minutes

# Create the temporary directory if it doesn't exist
if not os.path.exists(TEMP_DIR):
    os.makedirs(TEMP_DIR)

# --- Flask App Initialization ---
app = Flask(__name__)
app.config['MAX_CONTENT_LENGTH'] = MAX_CONTENT_LENGTH

# --- CORS Configuration ---
# IMPORTANT: Adjust allowed_origins to your actual frontend domain(s)
allowed_origins = [
    "https://sciencecodons.com", # Replace with your production domain
    "http://sciencecodons.com",  # Replace with your production domain
    "http://localhost:3000",            # For local frontend development
    "http://127.0.0.1:5637"             # Another common local dev port
]
CORS(app, resources={r"/api/*": {"origins": allowed_origins}})

# --- Temporary File Cleanup Scheduler ---
scheduler = BackgroundScheduler(daemon=True)

def delete_old_files():
    """Deletes files in the temp directory older than the configured lifetime."""
    try:
        current_time = time.time()
        for filename in os.listdir(TEMP_DIR):
            file_path = os.path.join(TEMP_DIR, filename)
            if os.path.isfile(file_path):
                file_age = current_time - os.path.getmtime(file_path)
                if file_age > FILE_LIFETIME_SECONDS:
                    os.remove(file_path)
                    logging.info(f"Deleted old file: {filename}")
    except Exception as e:
        logging.error(f"Error during file cleanup: {e}")

scheduler.add_job(func=delete_old_files, trigger="interval", minutes=5)
scheduler.start()
atexit.register(lambda: scheduler.shutdown())

# --- GenBank Conversion Logic (Adapted from your original script) ---

def _format_gff3_attribute_value(value: str) -> str:
    """
    Formats a string value for GFF3 attributes.
    Replaces internal double quotes with single quotes.
    URL-encodes problematic GFF3 characters (;,=,%,,,\\t,\\n,\\r) but *preserves spaces*.
    """
    # Replace internal double quotes with single quotes to avoid parsing issues
    value = value.replace('"', "'")

    # URL-encode other problematic characters that are GFF3 delimiters or special,
    # but explicitly add space to 'safe' characters to prevent %20 encoding.
    # The 'safe' set includes alphanumeric characters, and common symbols like /:.-_
    # Adding ' ' to safe will prevent spaces from being encoded.
    return quote(value, safe="/:.,-_ ")

def _get_common_qualifiers(feature, excluded_qualifiers: List[str] = None) -> Dict[str, Any]:
    """
    Extracts common and frequently used qualifiers from a Biopython SeqFeature.
    Returns a dictionary with default None or empty list values if a qualifier is not found.
    """
    if excluded_qualifiers is None:
        excluded_qualifiers = []

    qualifiers = feature.qualifiers
    common = {
        "gene": qualifiers.get("gene", [None])[0],
        "locus_tag": qualifiers.get("locus_tag", [None])[0],
        "product": qualifiers.get("product", [None])[0],
        "protein_id": qualifiers.get("protein_id", [None])[0],
        "transcript_id": qualifiers.get("transcript_id", [None])[0],
        "db_xref": qualifiers.get("db_xref", []),
        "note": qualifiers.get("note", [None])[0],
        "codon_start": qualifiers.get("codon_start", [None])[0],
        "transl_table": qualifiers.get("transl_table", [None])[0],
        "experiment": qualifiers.get("experiment", [None])[0],
        "function": qualifiers.get("function", [None])[0],
        "inference": qualifiers.get("inference", [None])[0],
        "old_locus_tag": qualifiers.get("old_locus_tag", [None])[0],
        "organism": qualifiers.get("organism", [None])[0],
        "mol_type": qualifiers.get("mol_type", [None])[0],
        "collection_date": qualifiers.get("collection_date", [None])[0],
        "country": qualifiers.get("country", [None])[0],
        "isolation_source": qualifiers.get("isolation_source", [None])[0],
        "host": qualifiers.get("host", [None])[0],
        "strain": qualifiers.get("strain", [None])[0],
        "segment": qualifiers.get("segment", [None])[0],
        "map": qualifiers.get("map", [None])[0],
        "allele": qualifiers.get("allele", [None])[0],
        "citation": qualifiers.get("citation", [None])[0],
        "compare": qualifiers.get("compare", [None])[0],
        "direction": qualifiers.get("direction", [None])[0],
        "exception": qualifiers.get("exception", [None])[0],
        "gap": qualifiers.get("gap", [None])[0],
        "label": qualifiers.get("label", [None])[0],
        "mobile_element_type": qualifiers.get("mobile_element_type", [None])[0],
        "rpt_family": qualifiers.get("rpt_family", [None])[0],
        "rpt_type": qualifiers.get("rpt_type", [None])[0],
        "standard_name": qualifiers.get("standard_name", [None])[0],
        "transl_except": qualifiers.get("transl_except", []),
        "translation": qualifiers.get("translation", [None])[0],
        "pseudo": "pseudo" in qualifiers, # Check for presence of 'pseudo' flag
        "partial": "partial" in qualifiers, # Check for presence of 'partial' flag
        "exception": qualifiers.get("exception", [None])[0]
    }

    # Add any other qualifiers not explicitly handled, unless they are in excluded_qualifiers
    for key, value in qualifiers.items():
        if key not in common and key not in excluded_qualifiers and value:
            common[key] = value[0] if isinstance(value, list) and len(value) == 1 else value

    return common

def _generate_fallback_id(seqname: str, feature_type: str, index: int, is_protein_id: bool = False) -> str:
    """Generates a consistent fallback ID."""
    prefix = "protein" if is_protein_id else feature_type
    return f"{prefix}_{seqname}_{index}"

def _calculate_gtf_frame(loc: FeatureLocation) -> str:
    """Calculates the GTF frame (0, 1, 2) for a CDS feature."""
    # GTF frame is 0-based, relative to the start of the CDS, indicating
    # the number of bases to skip at the start of the feature to reach the first complete codon.
    # Biopython's codon_start is 1-based.
    # frame = (3 - (codon_start - 1) % 3) % 3
    # If codon_start is 1, frame is 0. If codon_start is 2, frame is 2. If codon_start is 3, frame is 1.
    codon_start = loc.qualifiers.get("codon_start", ["1"])[0] # Default to 1 if not present
    try:
        cs = int(codon_start)
        return str((3 - (cs - 1) % 3) % 3)
    except ValueError:
        return "." # Invalid codon_start

def _calculate_gff3_phase(loc: FeatureLocation) -> str:
    """Calculates the GFF3 phase (0, 1, 2) for a CDS feature."""
    # GFF3 phase is 0-based, relative to the start of the CDS, indicating
    # the number of bases to skip at the start of the feature to reach the first complete codon.
    # It's the same as GTF frame.
    codon_start = loc.qualifiers.get("codon_start", ["1"])[0] # Default to 1 if not present
    try:
        cs = int(codon_start)
        return str((cs - 1) % 3)
    except ValueError:
        return "." # Invalid codon_start

def _get_partiality_attributes_gff3(loc: FeatureLocation) -> Dict[str, str]:
    """
    Determines GFF3 partiality attributes based on Biopython FeatureLocation.
    """
    attrs = {}
    if loc.start_original_position != 0: # Indicates 5' partiality
        attrs["five_prime_partial"] = "true"
    if loc.end_original_position != None and loc.end_original_position != loc.length: # Indicates 3' partiality
        attrs["three_prime_partial"] = "true"
    return attrs

def _process_feature_gtf(feature, seqname: str, source: str, feature_index: int,
                         gtf_gene_ids: Dict[str, str], gtf_transcript_ids: Dict[str, str],
                         record_length: int) -> str | None:
    """Processes a single Biopython SeqFeature for GTF output."""
    feature_type = feature.type
    # GTF typically only includes gene, transcript, exon, CDS, start_codon, stop_codon
    if feature_type not in ["gene", "mRNA", "CDS", "tRNA", "rRNA", "exon"]:
        logging.debug(f"Skipping GTF conversion for feature type: {feature_type}")
        return None

    # GTF coordinates are 1-based, inclusive
    start = int(feature.location.start) + 1
    end = int(feature.location.end)

    # Handle compound locations (e.g., for CDS across introns)
    if isinstance(feature.location, CompoundLocation):
        # For GTF, compound locations (like multi-exon CDS) are typically broken into individual exons/CDS segments.
        # This simplified approach will just take the overall start/end.
        # A more robust GTF converter would iterate through sub_features.
        # For now, we'll output a single line for the overall feature.
        pass # Start and end are already calculated from overall location

    strand = "+" if feature.location.strand == 1 else "-" if feature.location.strand == -1 else "."
    score = "."

    common_qualifiers = _get_common_qualifiers(feature, excluded_qualifiers=[
        "codon_start", "transl_table", "note", "db_xref", "experiment", "function",
        "inference", "old_locus_tag", "organism", "mol_type", "collection_date",
        "country", "isolation_source", "host", "strain", "segment", "map",
        "allele", "citation", "compare", "direction", "gap", "label",
        "mobile_element_type", "rpt_family", "rpt_type", "standard_name",
        "transl_except", "translation", "pseudo", "partial", "exception"
    ])

    attributes = []
    gene_name = common_qualifiers.get("gene") or common_qualifiers.get("locus_tag")
    product = common_qualifiers.get("product")
    protein_id = common_qualifiers.get("protein_id")
    transcript_id_qual = common_qualifiers.get("transcript_id")

    # Gene ID handling for GTF
    if gene_name:
        gene_id = gtf_gene_ids.setdefault(gene_name, f"gene_{gene_name}")
    else:
        gene_id = gtf_gene_ids.setdefault(f"{seqname}_{feature_index}_gene", f"gene_{seqname}_{feature_index}")
    attributes.append(f'gene_id "{gene_id}"')

    # Transcript ID handling for GTF
    if feature_type == "mRNA":
        if transcript_id_qual:
            transcript_id = gtf_transcript_ids.setdefault(transcript_id_qual, f"transcript_{transcript_id_qual}")
        elif gene_name:
            transcript_id = gtf_transcript_ids.setdefault(f"{gene_name}_transcript", f"transcript_{gene_name}_{feature_index}")
        else:
            transcript_id = gtf_transcript_ids.setdefault(f"{seqname}_{feature_index}_transcript", f"transcript_{seqname}_{feature_index}")
        attributes.append(f'transcript_id "{transcript_id}"')
    elif feature_type == "CDS" or feature_type == "exon":
        # For CDS/exon, try to link to an mRNA/transcript if available
        # This is a simplification; a full GTF conversion would need a pre-pass or more complex logic
        # to ensure CDS/exons are correctly nested under transcripts.
        # For now, we'll try to derive a transcript_id.
        parent_gene_name = common_qualifiers.get("gene") or common_qualifiers.get("locus_tag")
        if parent_gene_name and parent_gene_name in gtf_gene_ids:
            # Try to find an associated mRNA/transcript ID if one was generated for this gene
            # This is heuristic and might not always be correct without explicit parentage in GenBank
            derived_transcript_id = gtf_transcript_ids.get(f"{parent_gene_name}_transcript") or \
                                    gtf_transcript_ids.get(f"{seqname}_{feature_index}_transcript") # Fallback
            if derived_transcript_id:
                attributes.append(f'transcript_id "{derived_transcript_id}"')
            else: # Fallback if no derived transcript ID found
                attributes.append(f'transcript_id "{_generate_fallback_id(seqname, "transcript", feature_index)}"')
        else:
            attributes.append(f'transcript_id "{_generate_fallback_id(seqname, "transcript", feature_index)}"')


    if gene_name:
        attributes.append(f'gene_name "{_format_gff3_attribute_value(gene_name)}"')
    if product:
        attributes.append(f'product "{_format_gff3_attribute_value(product)}"')
    if protein_id:
        attributes.append(f'protein_id "{_format_gff3_attribute_value(protein_id)}"')

    frame = "."
    if feature_type == "CDS":
        frame = _calculate_gtf_frame(feature.location)

    # GTF: seqname source feature start end score strand frame attributes
    return f"{seqname}\t{source}\t{feature_type}\t{start}\t{end}\t{score}\t{strand}\t{frame}\t{'; '.join(attributes)};"

def _process_feature_gff3(feature, seqname: str, source: str, feature_index: int, record_length: int,
                          gff3_id_map: Dict[str, str], gene_id_bases: Dict[Any, str], mRNA_id_bases: Dict[Any, str],
                          current_source_id: str) -> str | None:
    """Processes a single Biopython SeqFeature for GFF3 output."""
    feature_type = feature.type

    # GFF3 coordinates are 1-based, inclusive
    start = int(feature.location.start) + 1
    end = int(feature.location.end)

    strand = "+" if feature.location.strand == 1 else "-" if feature.location.strand == -1 else "."
    score = "."
    phase = "." # For CDS features

    # Exclude qualifiers that are handled as specific GFF3 attributes or are not typically included
    excluded_qualifiers = [
        "codon_start", "transl_table", "translation", "pseudo", "partial", "exception",
        "experiment", "function", "inference", "old_locus_tag", "organism", "mol_type",
        "collection_date", "country", "isolation_source", "host", "strain", "segment",
        "map", "allele", "citation", "compare", "direction", "gap", "label",
        "mobile_element_type", "rpt_family", "rpt_type", "standard_name", "transl_except"
    ]
    if feature_type == "CDS": # NCBI GFF3 often excludes 'note' for CDS
        excluded_qualifiers.append("note")

    common_qualifiers = _get_common_qualifiers(feature, excluded_qualifiers)

    attributes = []
    feature_id = ""
    parent_ids = []

    # Handle feature IDs and Parent relationships
    if feature_type == "gene":
        gene_id_base = common_qualifiers.get("gene") or common_qualifiers.get("locus_tag")
        if gene_id_base:
            feature_id = gene_id_bases.get((seqname, gene_id_base), _generate_fallback_id(seqname, "gene", feature_index))
        else:
            feature_id = gene_id_bases.get((seqname, feature_index), _generate_fallback_id(seqname, "gene", feature_index))
        gff3_id_map[f"{seqname}_gene_{feature_index}"] = feature_id # Map for later reference
        parent_ids.append(current_source_id) # Gene is child of the sequence region

    elif feature_type == "mRNA":
        mrna_id_base = common_qualifiers.get("product")
        if mrna_id_base:
            feature_id = mRNA_id_bases.get((seqname, mrna_id_base), _generate_fallback_id(seqname, "mRNA", feature_index))
        else:
            feature_id = mRNA_id_bases.get((seqname, feature_index), _generate_fallback_id(seqname, "mRNA", feature_index))
        gff3_id_map[f"{seqname}_mRNA_{feature_index}"] = feature_id # Map for later reference

        # mRNA's parent is the gene
        gene_qual = common_qualifiers.get("gene") or common_qualifiers.get("locus_tag")
        if gene_qual and (seqname, gene_qual) in gene_id_bases:
            parent_ids.append(gene_id_bases[(seqname, gene_qual)])
        else:
            # Fallback: try to find a gene parent by proximity or a general gene ID for this sequence
            # This is a simplification; robust linking needs more context
            parent_ids.append(current_source_id) # Default to source if gene parent not found easily

    elif feature_type == "CDS":
        # CDS ID
        protein_id = common_qualifiers.get("protein_id")
        if protein_id:
            feature_id = f"cds_{_format_gff3_attribute_value(protein_id)}"
        else:
            feature_id = _generate_fallback_id(seqname, "CDS", feature_index)
        gff3_id_map[f"{seqname}_CDS_{feature_index}"] = feature_id

        # CDS parent is mRNA (if available) or gene
        mrna_product = common_qualifiers.get("product")
        gene_qual = common_qualifiers.get("gene") or common_qualifiers.get("locus_tag")

        if mrna_product and (seqname, mrna_product) in mRNA_id_bases:
            parent_ids.append(mRNA_id_bases[(seqname, mrna_product)])
        elif gene_qual and (seqname, gene_qual) in gene_id_bases:
            parent_ids.append(gene_id_bases[(seqname, gene_qual)])
        else:
            # Fallback: link to the sequence region
            parent_ids.append(current_source_id)

        phase = _calculate_gff3_phase(feature.location)

    elif feature_type == "exon":
        feature_id = _generate_fallback_id(seqname, "exon", feature_index)
        gff3_id_map[f"{seqname}_exon_{feature_index}"] = feature_id

        # Exon parent is mRNA
        mrna_product = common_qualifiers.get("product")
        if mrna_product and (seqname, mrna_product) in mRNA_id_bases:
            parent_ids.append(mRNA_id_bases[(seqname, mrna_product)])
        else:
            # Fallback: link to the sequence region
            parent_ids.append(current_source_id)

    elif feature_type in ["tRNA", "rRNA"]:
        # For tRNA/rRNA, they might be top-level or children of a gene.
        # For simplicity, we'll make them children of the gene if a gene qualifier exists,
        # otherwise children of the source.
        gene_qual = common_qualifiers.get("gene") or common_qualifiers.get("locus_tag")
        if gene_qual and (seqname, gene_qual) in gene_id_bases:
            parent_ids.append(gene_id_bases[(seqname, gene_qual)])
        else:
            parent_ids.append(current_source_id)
        feature_id = _generate_fallback_id(seqname, feature_type, feature_index)
        gff3_id_map[f"{seqname}_{feature_type}_{feature_index}"] = feature_id

    else:
        # For other feature types, make them children of the source region
        feature_id = _generate_fallback_id(seqname, feature_type, feature_index)
        parent_ids.append(current_source_id)


    # Add ID and Parent attributes
    if feature_id:
        attributes.append(f"ID={feature_id}")
    if parent_ids:
        attributes.append(f"Parent={','.join(parent_ids)}")

    # Add Name attribute if available and not already used as ID
    name = common_qualifiers.get("gene") or common_qualifiers.get("locus_tag") or common_qualifiers.get("product")
    if name and name != feature_id: # Avoid redundant ID/Name if they are the same
        attributes.append(f"Name={_format_gff3_attribute_value(name)}")

    # Add Dbxref
    db_xrefs = common_qualifiers.get("db_xref", [])
    if feature_type == "CDS" and common_qualifiers.get("protein_id"):
        db_xrefs.append(f"NCBI_GP:{common_qualifiers['protein_id']}")
    if db_xrefs:
        attributes.append(f"Dbxref={','.join([_format_gff3_attribute_value(x) for x in db_xrefs])}")

    # Add other common qualifiers as GFF3 attributes
    for key, value in common_qualifiers.items():
        if value is not None and key not in ["gene", "locus_tag", "product", "protein_id", "transcript_id", "db_xref", "note", "codon_start", "transl_table"] and key not in excluded_qualifiers:
            if isinstance(value, list): # Handle list values (e.g., transl_except)
                attributes.append(f"{key}={','.join([_format_gff3_attribute_value(str(v)) for v in value])}")
            else:
                attributes.append(f"{key}={_format_gff3_attribute_value(str(value))}")

    # Add partiality attributes
    partiality_attrs = _get_partiality_attributes_gff3(feature.location)
    for k, v in partiality_attrs.items():
        attributes.append(f"{k}={v}")

    # Add note if present and not excluded
    note = common_qualifiers.get("note")
    if note and "note" not in excluded_qualifiers:
        attributes.append(f"Note={_format_gff3_attribute_value(note)}")

    # GFF3: seqname source feature start end score strand phase attributes
    return f"{seqname}\t{source}\t{feature_type}\t{start}\t{end}\t{score}\t{strand}\t{phase}\t{';'.join(attributes)}"


def convert_genbank_data_to_annotation_format(genbank_data: str, output_format: str) -> tuple[str | None, str | None]:
    """
    Converts GenBank data (as a string) to GTF or GFF3 format.
    Returns the converted string and an error message (if any).
    """
    if output_format not in ["gtf", "gff3"]:
        return None, "Invalid output format specified. Must be 'gtf' or 'gff3'."

    converted_lines = []
    gff3_id_map = {} # To store IDs for GFF3 parent-child relationships

    # Pre-processing for GFF3 parent-child relationships (gene and mRNA IDs)
    # This helps ensure Parent attributes can reference existing IDs.
    gene_id_bases = {}
    mRNA_id_bases = {}

    try:
        # Use StringIO to treat the string data as a file
        genbank_handle = StringIO(genbank_data)
        records = list(SeqIO.parse(genbank_handle, "genbank"))
        genbank_handle.close() # Close StringIO handle

        if not records:
            return None, "No valid GenBank records found in the provided data."

        for record_index, record in enumerate(records):
            seqname = record.id if record.id else record.name if record.name else f"sequence_{record_index}"
            record_length = len(record.seq)
            current_source_id = f"region_{seqname}" # Unique ID for the source region of this record

            # GFF3 Header and Source Region feature
            if output_format == "gff3":
                if record_index == 0: # Only add GFF version once at the very beginning
                    converted_lines.append("##gff-version 3")

                source_qualifiers = {}
                for f in record.features:
                    if f.type == "source":
                        source_qualifiers = {k: v[0] for k, v in f.qualifiers.items() if v}
                        break

                record_description = source_qualifiers.get("organism", seqname)
                # Adding a source region feature
                converted_lines.append(
                    f"{seqname}\tGenBank\tregion\t1\t{record_length}\t.\t.\t.\tID={current_source_id};Name={seqname};accession={record.id};Dbxref=taxon:{source_qualifiers.get('db_xref', [''])[0].replace('taxon:', '')};description={_format_gff3_attribute_value(record_description)}"
                )

            # Pre-pass for GFF3 to identify gene and mRNA IDs for parent-child linking
            if output_format == "gff3":
                for feature_index, feature in enumerate(record.features):
                    if feature.type == "gene":
                        gene_qual_val = feature.qualifiers.get("gene", [None])[0] or feature.qualifiers.get("locus_tag", [None])[0]
                        if gene_qual_val:
                            gene_id_bases[(seqname, gene_qual_val)] = f"gene_{seqname}_{_format_gff3_attribute_value(gene_qual_val)}"
                        else:
                            # Fallback if no gene/locus_tag qualifier
                            gene_id_bases[(seqname, feature_index)] = _generate_fallback_id(seqname, "gene", feature_index)
                    elif feature.type == "mRNA":
                        mrna_qual_val = feature.qualifiers.get("product", [None])[0]
                        if mrna_qual_val:
                            mRNA_id_bases[(seqname, mrna_qual_val)] = f"mrna_{seqname}_{_format_gff3_attribute_value(mrna_qual_val)}"
                        else:
                            # Fallback if no product qualifier
                            mRNA_id_bases[(seqname, feature_index)] = _generate_fallback_id(seqname, "mRNA", feature_index)


            gtf_gene_ids = {} # For GTF gene_id tracking (per record)
            gtf_transcript_ids = {} # For GTF transcript_id tracking (per record)

            # Sort features to ensure parents appear before children in GFF3
            # (e.g., gene before mRNA before CDS/exon)
            # This is a heuristic sort, not guaranteed to be perfect for all GenBank files.
            feature_sort_order = {"source": 0, "gene": 1, "mRNA": 2, "tRNA": 3, "rRNA": 4, "CDS": 5, "exon": 6}
            sorted_features = sorted(record.features, key=lambda f: feature_sort_order.get(f.type, 99))

            for feature_index, feature in enumerate(sorted_features):
                line = None
                if output_format == "gtf":
                    line = _process_feature_gtf(feature, seqname, "GenBank", feature_index, gtf_gene_ids, gtf_transcript_ids, record_length)
                elif output_format == "gff3":
                    line = _process_feature_gff3(
                        feature, seqname, "GenBank", feature_index, record_length,
                        gff3_id_map, gene_id_bases, mRNA_id_bases, current_source_id
                    )

                if line:
                    converted_lines.append(line)

            if output_format == "gff3" and record_index < len(records) - 1:
                converted_lines.append("###") # Delimiter between sequences in GFF3

    except Exception as e:
        logging.error(f"Conversion error: {e}", exc_info=True)
        return None, f"An error occurred during conversion: {str(e)}"

    return "\n".join(converted_lines), None


# --- API Endpoints ---
@app.route('/api/convert', methods=['POST'])
def handle_genbank_conversion():
    """Handles the file upload/paste, conversion, and response."""

    # Check for file or raw data
    if 'file' not in request.files and 'genbank_data' not in request.form:
        return jsonify({"success": False, "error": "No GenBank data provided. Please upload a file or paste data."}), 400

    # Determine output format
    output_format = request.form.get('output_format', 'gff3').lower()
    if output_format not in ['gtf', 'gff3']:
        return jsonify({"success": False, "error": "Invalid output format. Choose 'gtf' or 'gff3'."}), 400

    genbank_data = ""
    try:
        # Read data from file upload
        if 'file' in request.files and request.files['file'].filename != '':
            genbank_file = request.files['file']
            # Basic file type validation
            if not genbank_file.filename.lower().endswith(('.gb', '.gbk', '.genbank')):
                return jsonify({"success": False, "error": "Invalid file type. Please upload a .gb, .gbk, or .genbank file."}), 400
            genbank_data = genbank_file.read().decode('utf-8')
        # Read data from raw text input
        elif 'genbank_data' in request.form:
            genbank_data = request.form['genbank_data']
        else:
            return jsonify({"success": False, "error": "No GenBank data found in the request."}), 400

        if not genbank_data.strip():
            return jsonify({"success": False, "error": "Submitted GenBank data is empty."}), 400

    except Exception as e:
        logging.error(f"Error processing input: {e}", exc_info=True)
        return jsonify({"success": False, "error": f"Error processing input: {e}"}), 400

    # Perform the conversion
    converted_result, error_message = convert_genbank_data_to_annotation_format(genbank_data, output_format)

    if error_message:
        return jsonify({"success": False, "error": error_message}), 400

    # Generate a unique filename for the output
    unique_filename = f"{uuid.uuid4()}.{output_format}"
    output_path = os.path.join(TEMP_DIR, unique_filename)

    try:
        with open(output_path, "w") as f:
            f.write(converted_result)
    except IOError as e:
        logging.error(f"Error writing file to temporary storage: {e}", exc_info=True)
        return jsonify({"success": False, "error": "A server error occurred while preparing the file for download."}), 500

    # Construct download URL
    download_url = request.host_url.rstrip('/') + f"/api/download/{unique_filename}"

    return jsonify({
        "success": True,
        "converted_data": converted_result, # Optionally return data directly for smaller files
        "download_url": download_url
    })

@app.route('/api/download/<filename>')
def download_converted_file(filename):
    """Serves the converted file for download."""
    # Security: Use secure_filename to prevent directory traversal
    from werkzeug.utils import secure_filename
    safe_filename = secure_filename(filename)

    # Ensure the filename matches the secured version to prevent malicious paths
    if safe_filename != filename:
        return jsonify({"success": False, "error": "Invalid filename."}), 400

    # Ensure the file exists and is within the TEMP_DIR
    full_file_path = os.path.join(TEMP_DIR, safe_filename)
    if not os.path.exists(full_file_path):
        logging.warning(f"Attempted to download non-existent file: {full_file_path}")
        return jsonify({"success": False, "error": "File not found or has expired."}), 404

    return send_from_directory(TEMP_DIR, safe_filename, as_attachment=True)

@app.errorhandler(413)
def request_entity_too_large(error):
    """Handles requests where the payload size exceeds MAX_CONTENT_LENGTH."""
    return jsonify({"success": False, "error": f"The file is too large. The maximum allowed size is {MAX_CONTENT_LENGTH / 1024 / 1024:.0f}MB."}), 413

@app.route('/')
def serve_index():
    """Simple health check endpoint."""
    return "GenBank to GTF/GFF3 Converter API is running. Use /api/convert to interact."

if __name__ == '__main__':
    # For local development, run on port 5000
    app.run(host='0.0.0.0', port=5000, debug=True)
