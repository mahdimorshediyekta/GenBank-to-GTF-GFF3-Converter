import os
import uuid
import time
import atexit
import logging
from io import StringIO
from urllib.parse import quote
from typing import Dict, List, Any

from flask import Flask, request, jsonify, send_from_directory
from flask_cors import CORS
from apscheduler.schedulers.background import BackgroundScheduler

from Bio import SeqIO
from Bio.SeqFeature import (
    FeatureLocation, CompoundLocation,
    BeforePosition, AfterPosition, UnknownPosition
)

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
# For production, replace 'http://localhost:3000' and 'http://127.0.0.1:5637'
# with your actual hosted frontend domains (e.g., 'https://yourfrontend.com').
# Avoid using "*" in production for security reasons.
allowed_origins = [
    "http://localhost:3000",            # For local frontend development
    "http://127.0.0.1:5637",            # Another common local dev port
    # "https://your-genbank-website.com", # Replace with your production domain
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

# --- GenBank Conversion Logic ---

def _format_gff3_attribute_value(value: str) -> str:
    """
    Formats a string value for GFF3 attributes.
    Replaces internal double quotes with single quotes.
    URL-encodes problematic GFF3 characters (;,=,%,,\t,\n,\r) but *preserves spaces*.
    """
    value = value.replace('"', "'")
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
        "pseudo": "pseudo" in qualifiers,
        "partial": "partial" in qualifiers,
    }

    for key, value in qualifiers.items():
        if value is not None and key not in common and key not in excluded_qualifiers:
            common[key] = value[0] if isinstance(value, list) and len(value) == 1 else value

    return common

def _generate_fallback_id(seqname: str, feature_type: str, index: int) -> str:
    """Generates a consistent fallback ID."""
    # Ensure ID parts are safe for GFF3/GTF by replacing spaces
    safe_seqname = seqname.replace(' ', '_')
    safe_feature_type = feature_type.replace(' ', '_')
    return f"{safe_feature_type}_{safe_seqname}_{index}"

def _calculate_gtf_frame(feature) -> str:
    """Calculates the GTF frame (0, 1, 2) for a CDS feature."""
    codon_start = feature.qualifiers.get("codon_start", ["1"])[0]
    try:
        cs = int(codon_start)
        return str((3 - (cs - 1) % 3) % 3)
    except ValueError:
        return "."

def _calculate_gff3_phase(feature) -> str:
    """Calculates the GFF3 phase (0, 1, 2) for a CDS feature."""
    codon_start = feature.qualifiers.get("codon_start", ["1"])[0]
    try:
        cs = int(codon_start)
        return str((cs - 1) % 3)
    except ValueError:
        return "."

def _get_partiality_attributes_gff3(loc: FeatureLocation) -> Dict[str, str]:
    """
    Determines GFF3 partiality attributes based on Biopython FeatureLocation.
    Checks for specific partial position types.
    """
    attrs = {}
    if isinstance(loc.start, (BeforePosition, UnknownPosition)):
        attrs["five_prime_partial"] = "true"
    if isinstance(loc.end, (AfterPosition, UnknownPosition)):
        attrs["three_prime_partial"] = "true"
    return attrs

def _process_feature_gtf(feature, seqname: str, source: str, feature_index: int,
                         gtf_gene_ids: Dict[str, str], gtf_transcript_ids: Dict[str, str]) -> str | None:
    """Processes a single Biopython SeqFeature for GTF output."""
    feature_type = feature.type
    if feature_type not in ["gene", "mRNA", "CDS", "tRNA", "rRNA", "exon"]:
        logging.debug(f"Skipping GTF conversion for feature type: {feature_type}")
        return None

    start = int(feature.location.start) + 1
    end = int(feature.location.end)

    strand = "+" if feature.location.strand == 1 else "-" if feature.location.strand == -1 else "."
    score = "."

    common_qualifiers = _get_common_qualifiers(feature, excluded_qualifiers=[
        "codon_start", "transl_table", "note", "db_xref", "translation", "pseudo", "partial", "exception"
    ])

    attributes = []
    gene_name = common_qualifiers.get("gene") or common_qualifiers.get("locus_tag")
    product = common_qualifiers.get("product")
    protein_id = common_qualifiers.get("protein_id")
    transcript_id_qual = common_qualifiers.get("transcript_id")

    # Gene ID handling for GTF
    if gene_name:
        gene_id = gtf_gene_ids.setdefault(gene_name, f"gene_{_format_gff3_attribute_value(gene_name.replace(' ', '_'))}")
    else:
        gene_id = gtf_gene_ids.setdefault(f"{seqname}_{feature_index}_gene", _generate_fallback_id(seqname, "gene", feature_index))
    attributes.append(f'gene_id "{gene_id}"')

    # Transcript ID handling for GTF
    current_transcript_id = None
    if feature_type == "mRNA":
        if transcript_id_qual:
            current_transcript_id = gtf_transcript_ids.setdefault(transcript_id_qual, f"transcript_{_format_gff3_attribute_value(transcript_id_qual.replace(' ', '_'))}")
        elif gene_name:
            current_transcript_id = gtf_transcript_ids.setdefault(f"{gene_name}_transcript", f"transcript_{_format_gff3_attribute_value(gene_name.replace(' ', '_'))}_{feature_index}")
        else:
            current_transcript_id = gtf_transcript_ids.setdefault(f"{seqname}_{feature_index}_mRNA_transcript", _generate_fallback_id(seqname, "mRNA_transcript", feature_index))
    elif feature_type in ["CDS", "exon"]:
        # Try to link to an mRNA/transcript based on gene_name first
        if gene_name and gene_name in gtf_gene_ids:
            # Check if an mRNA transcript was already generated for this gene
            derived_transcript_id = gtf_transcript_ids.get(f"{gene_name}_transcript")
            if derived_transcript_id:
                current_transcript_id = derived_transcript_id
            elif protein_id:
                current_transcript_id = f"transcript_{_format_gff3_attribute_value(protein_id.replace(' ', '_'))}"
        if not current_transcript_id: # Fallback if no specific transcript ID derived
             current_transcript_id = _generate_fallback_id(seqname, f"{feature_type}_transcript", feature_index)
    elif feature_type in ["tRNA", "rRNA"]:
        if gene_name:
            current_transcript_id = f"transcript_{_format_gff3_attribute_value(gene_name.replace(' ', '_'))}_{feature_type}"
        else:
            current_transcript_id = _generate_fallback_id(seqname, f"{feature_type}_transcript", feature_index)


    if current_transcript_id:
        attributes.append(f'transcript_id "{current_transcript_id}"')

    if gene_name:
        attributes.append(f'gene_name "{_format_gff3_attribute_value(gene_name)}"')
    if product:
        attributes.append(f'product "{_format_gff3_attribute_value(product)}"')
    if protein_id:
        attributes.append(f'protein_id "{_format_gff3_attribute_value(protein_id)}"')
    if common_qualifiers.get("pseudo"):
        attributes.append('pseudo "true"')
    if common_qualifiers.get("exception"):
        attributes.append(f'exception "{_format_gff3_attribute_value(common_qualifiers["exception"])}"')
    if common_qualifiers.get("transl_table") and feature_type == "CDS":
        attributes.append(f'transl_table "{common_qualifiers["transl_table"]}"')
    
    # Add other remaining qualifiers
    for qual_key, qual_value in common_qualifiers.items():
        if qual_key not in ["gene", "locus_tag", "product", "protein_id", "transcript_id", "pseudo", "exception", "transl_table", "db_xref", "note", "codon_start", "translation", "partial"]:
            if qual_value is not None:
                if isinstance(qual_value, list): # For qualifiers like transl_except
                    attributes.append(f'{qual_key} "{",".join([_format_gff3_attribute_value(str(v)) for v in qual_value])}"')
                else:
                    attributes.append(f'{qual_key} "{_format_gff3_attribute_value(str(qual_value))}"')

    frame = "."
    if feature_type == "CDS":
        frame = _calculate_gtf_frame(feature)

    return f"{seqname}\t{source}\t{feature_type}\t{start}\t{end}\t{score}\t{strand}\t{frame}\t{'; '.join(attributes)};"


def _process_feature_gff3(feature, seqname: str, source: str, feature_index: int,
                          gff3_id_map: Dict[str, str], gene_id_bases: Dict[Any, str], mRNA_id_bases: Dict[Any, str],
                          current_source_id: str) -> str | None:
    """Processes a single Biopython SeqFeature for GFF3 output."""
    feature_type = feature.type

    start = int(feature.location.start) + 1
    end = int(feature.location.end)

    strand = "+" if feature.location.strand == 1 else "-" if feature.location.strand == -1 else "."
    score = "."
    phase = "."

    excluded_qualifiers_gff3 = [
        "codon_start", "transl_table", "translation", "pseudo", "partial", "exception",
        "gene", "locus_tag", "product", "protein_id", "transcript_id", "db_xref", "note"
    ]
    if feature_type == "CDS":
        excluded_qualifiers_gff3.append("note") # NCBI GFF3 often excludes 'note' for CDS

    common_qualifiers = _get_common_qualifiers(feature, excluded_qualifiers_gff3)

    attributes = []
    feature_id = ""
    parent_ids = []

    # Handle feature IDs and Parent relationships
    if feature_type == "gene":
        gene_id_base_key = common_qualifiers.get("gene") or common_qualifiers.get("locus_tag")
        if gene_id_base_key:
            feature_id = gene_id_bases.get((seqname, gene_id_base_key))
        if not feature_id: # Fallback if not found in pre-pass map
            feature_id = _generate_fallback_id(seqname, "gene", feature_index)
        gff3_id_map[f"{seqname}_gene_{feature_index}"] = feature_id # Map for later reference
        
        # Genes are typically top-level features, but in some cases, they might be children of a larger region.
        # For GenBank conversion, we'll make them children of the sequence region.
        parent_ids.append(current_source_id)

    elif feature_type == "mRNA":
        mrna_id_base_key = common_qualifiers.get("product") or common_qualifiers.get("transcript_id") or common_qualifiers.get("gene") or common_qualifiers.get("locus_tag")
        if mrna_id_base_key:
            feature_id = mRNA_id_bases.get((seqname, mrna_id_base_key))
        if not feature_id: # Fallback if not found in pre-pass map
            feature_id = _generate_fallback_id(seqname, "mRNA", feature_index)
        gff3_id_map[f"{seqname}_mRNA_{feature_index}"] = feature_id

        # mRNA's parent is the gene
        gene_qual_for_parent = common_qualifiers.get("gene") or common_qualifiers.get("locus_tag")
        if gene_qual_for_parent and (seqname, gene_qual_for_parent) in gene_id_bases:
            parent_ids.append(gene_id_bases[(seqname, gene_qual_for_parent)])
        elif current_source_id: # Fallback to source region
            parent_ids.append(current_source_id)

    elif feature_type == "CDS":
        protein_id = common_qualifiers.get("protein_id")
        if protein_id:
            feature_id = f"cds-{_format_gff3_attribute_value(protein_id)}"
        else:
            feature_id = _generate_fallback_id(seqname, "CDS", feature_index)
        gff3_id_map[f"{seqname}_CDS_{feature_index}"] = feature_id

        # CDS parent is mRNA (preferable) or gene
        mrna_parent_key = common_qualifiers.get("product") or common_qualifiers.get("transcript_id") or common_qualifiers.get("gene") or common_qualifiers.get("locus_tag")
        if mrna_parent_key and (seqname, mrna_parent_key) in mRNA_id_bases:
            parent_ids.append(mRNA_id_bases[(seqname, mrna_parent_key)])
        else:
            gene_parent_key = common_qualifiers.get("gene") or common_qualifiers.get("locus_tag")
            if gene_parent_key and (seqname, gene_parent_key) in gene_id_bases:
                parent_ids.append(gene_id_bases[(seqname, gene_parent_key)])
            elif current_source_id: # Fallback to source region
                parent_ids.append(current_source_id)
        
        phase = _calculate_gff3_phase(feature)

    elif feature_type == "exon":
        # Exon ID - often based on mRNA ID
        mrna_parent_key = common_qualifiers.get("product") or common_qualifiers.get("transcript_id") or common_qualifiers.get("gene") or common_qualifiers.get("locus_tag")
        if mrna_parent_key and (seqname, mrna_parent_key) in mRNA_id_bases:
            parent_mrna_id = mRNA_id_bases[(seqname, mrna_parent_key)]
            feature_id = f"{parent_mrna_id.replace('rna-', 'exon-')}_part{feature_index}" # Derive from mRNA ID
        else:
            feature_id = _generate_fallback_id(seqname, "exon", feature_index)
        gff3_id_map[f"{seqname}_exon_{feature_index}"] = feature_id

        # Exon parent is mRNA
        if mrna_parent_key and (seqname, mrna_parent_key) in mRNA_id_bases:
            parent_ids.append(mRNA_id_bases[(seqname, mrna_parent_key)])
        elif current_source_id: # Fallback to source region
            parent_ids.append(current_source_id)

    elif feature_type in ["tRNA", "rRNA"]:
        gene_qual_for_parent = common_qualifiers.get("gene") or common_qualifiers.get("locus_tag")
        if gene_qual_for_parent and (seqname, gene_qual_for_parent) in gene_id_bases:
            parent_ids.append(gene_id_bases[(seqname, gene_qual_for_parent)])
        elif current_source_id:
            parent_ids.append(current_source_id)
        
        id_base = common_qualifiers.get("product") or common_qualifiers.get("gene") or common_qualifiers.get("locus_tag")
        if id_base:
            feature_id = f"{feature_type}-{_format_gff3_attribute_value(id_base.replace(' ', '_'))}"
        else:
            feature_id = _generate_fallback_id(seqname, feature_type, feature_index)
        gff3_id_map[f"{seqname}_{feature_type}_{feature_index}"] = feature_id
    else:
        # For other feature types, link to the source region
        feature_id = _generate_fallback_id(seqname, feature_type, feature_index)
        parent_ids.append(current_source_id)


    if feature_id:
        attributes.append(f"ID={_format_gff3_attribute_value(feature_id)}")
    if parent_ids:
        attributes.append(f"Parent={','.join([_format_gff3_attribute_value(pid) for pid in parent_ids])}")

    name = common_qualifiers.get("gene") or common_qualifiers.get("locus_tag") or common_qualifiers.get("product")
    if name: # GFF3 Name is often a display name, not necessarily unique
        attributes.append(f"Name={_format_gff3_attribute_value(name)}")

    # Add Dbxref
    db_xrefs = common_qualifiers.get("db_xref", [])
    if feature_type == "CDS" and common_qualifiers.get("protein_id"):
        db_xrefs.insert(0, f"NCBI_GP:{common_qualifiers['protein_id']}") # Add protein_id as Dbxref
    if db_xrefs:
        attributes.append(f"Dbxref={','.join([_format_gff3_attribute_value(x) for x in db_xrefs])}")

    # Add note if present and not excluded for this feature type
    note = common_qualifiers.get("note")
    if note and "note" not in excluded_qualifiers_gff3:
        attributes.append(f"Note={_format_gff3_attribute_value(note)}")

    # Add gbkey
    gbkey_map = {
        "gene": "Gene", "mRNA": "mRNA", "CDS": "CDS", "exon": "mRNA",
        "tRNA": "tRNA", "rRNA": "rRNA", "source": "Src" # Added source
    }
    if feature_type in gbkey_map:
        attributes.append(f"gbkey={gbkey_map[feature_type]}")

    # Add gene_biotype for gene features
    if feature_type == "gene":
        attributes.append(f'gene_biotype=protein_coding') # Assuming protein_coding for genes from GenBank

    # Add pseudo attribute
    if common_qualifiers.get("pseudo"):
        attributes.append(f'pseudo=true')

    # Add exception attribute
    if common_qualifiers.get("exception"):
        attributes.append(f'exception={_format_gff3_attribute_value(common_qualifiers["exception"])}')
    
    # Add transl_table for CDS
    if common_qualifiers.get("transl_table") and feature_type == "CDS":
        attributes.append(f'transl_table={common_qualifiers["transl_table"]}')

    # Add gene and product attributes for CDS (as seen in NCBI GFF3)
    if feature_type == "CDS":
        if common_qualifiers.get("gene"):
            attributes.append(f'gene={_format_gff3_attribute_value(common_qualifiers["gene"])}')
        if common_qualifiers.get("product"):
            attributes.append(f'product={_format_gff3_attribute_value(common_qualifiers["product"])}')
        if common_qualifiers.get("protein_id"):
            attributes.append(f'protein_id={_format_gff3_attribute_value(common_qualifiers["protein_id"])}')


    # Add partiality attributes (five_prime_partial, three_prime_partial)
    partiality_attrs = _get_partiality_attributes_gff3(feature.location)
    for k, v in partiality_attrs.items():
        attributes.append(f"{k}={v}")

    # Add any other remaining qualifiers
    for qual_key, qual_value in common_qualifiers.items():
        if qual_key not in excluded_qualifiers_gff3:
            if qual_value is not None:
                if isinstance(qual_value, list):
                    attributes.append(f"{qual_key}={','.join([_format_gff3_attribute_value(str(v)) for v in qual_value])}")
                else:
                    attributes.append(f"{qual_key}={_format_gff3_attribute_value(str(qual_value))}")

    return f"{seqname}\t{source}\t{feature_type}\t{start}\t{end}\t{score}\t{strand}\t{phase}\t{';'.join(attributes)}"


def convert_genbank_data_to_annotation_format(genbank_data: str, output_format: str) -> tuple[str | None, str | None]:
    """
    Converts GenBank data (as a string) to GTF or GFF3 format.
    Returns the converted string and an error message (if any).
    """
    if output_format not in ["gtf", "gff3"]:
        return None, "Invalid output format specified. Must be 'gtf' or 'gff3'."

    converted_lines = []

    # Dictionaries for ID tracking across features and records
    # GFF3 needs to know parent IDs from previous features.
    gff3_id_map = {} # Generic map for any GFF3 ID generated
    gene_id_bases = {} # Stores base IDs for gene features (seqname, qual_value) -> actual_id
    mRNA_id_bases = {} # Stores base IDs for mRNA features (seqname, qual_value) -> actual_id

    try:
        genbank_handle = StringIO(genbank_data)
        records = list(SeqIO.parse(genbank_handle, "genbank"))
        genbank_handle.close()

        if not records:
            return None, "No valid GenBank records found in the provided data."

        if output_format == "gtf":
            converted_lines.append("##gff-version 2.2") # GTF often uses this header, though not strictly part of its spec.
        elif output_format == "gff3":
            converted_lines.append("##gff-version 3")
            converted_lines.append("##gff-spec-version 1.21")


        for record_index, record in enumerate(records):
            seqname = record.id if record.id else record.name if record.name else f"sequence_{record_index}"
            record_length = len(record.seq)
            current_source_id = f"region-{_format_gff3_attribute_value(seqname)}"

            # --- GFF3 Specific Header and Region Feature ---
            if output_format == "gff3":
                converted_lines.append(f"##sequence-region {seqname} 1 {record_length}")

                source_qualifiers = {}
                for f in record.features:
                    if f.type == "source":
                        source_qualifiers = {k: v[0] for k, v in f.qualifiers.items() if v}
                        break
                
                # Add ##species line if taxon is found
                taxon_ids = [x.split(':')[-1] for x in source_qualifiers.get("db_xref", []) if "taxon:" in x]
                if taxon_ids:
                    converted_lines.append(f"##species https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id={taxon_ids[0]}") # Only first taxon ID

                # Region feature
                region_attributes = [f'ID={current_source_id}', 'gbkey=Src', 'mol_type=genomic DNA']
                region_attributes.append(f'Name={_format_gff3_attribute_value(seqname)}')
                if record.id:
                    region_attributes.append(f'accession={_format_gff3_attribute_value(record.id)}')
                if source_qualifiers.get("organism"):
                    region_attributes.append(f'organism={_format_gff3_attribute_value(source_qualifiers["organism"])}')
                if source_qualifiers.get("strain"):
                    region_attributes.append(f'strain={_format_gff3_attribute_value(source_qualifiers["strain"])}')
                
                # Add all db_xrefs from source feature to region feature
                source_db_xrefs = source_qualifiers.get("db_xref", [])
                if source_db_xrefs:
                    region_attributes.append(f'Dbxref={",".join([_format_gff3_attribute_value(x) for x in source_db_xrefs])}')
                
                converted_lines.append(
                    f"{seqname}\tGenBank\tregion\t1\t{record_length}\t.\t.\t.\t{';'.join(region_attributes)}"
                )

            # Pre-pass for GFF3 to identify gene and mRNA IDs for parent-child linking
            if output_format == "gff3":
                for feature_index, feature in enumerate(record.features):
                    if feature.type == "gene":
                        gene_qual_val = feature.qualifiers.get("gene", [None])[0] or feature.qualifiers.get("locus_tag", [None])[0]
                        if gene_qual_val:
                            gene_id_bases[(seqname, gene_qual_val)] = f"gene-{_format_gff3_attribute_value(gene_qual_val.replace(' ', '_'))}"
                        else:
                            gene_id_bases[(seqname, feature_index)] = _generate_fallback_id(seqname, "gene", feature_index)
                    elif feature.type == "mRNA":
                        # mRNA ID can be based on product, transcript_id, gene, or locus_tag
                        mrna_qual_val = feature.qualifiers.get("product", [None])[0] or \
                                        feature.qualifiers.get("transcript_id", [None])[0] or \
                                        feature.qualifiers.get("gene", [None])[0] or \
                                        feature.qualifiers.get("locus_tag", [None])[0]
                        if mrna_qual_val:
                            mRNA_id_bases[(seqname, mrna_qual_val)] = f"rna-{_format_gff3_attribute_value(mrna_qual_val.replace(' ', '_'))}"
                        else:
                            mRNA_id_bases[(seqname, feature_index)] = _generate_fallback_id(seqname, "mRNA", feature_index)

            # Sort features to ensure parents appear before children in GFF3
            # and to maintain a somewhat logical order (source, gene, mRNA, tRNA, rRNA, exon, CDS)
            feature_sort_order = {"source": 0, "gene": 1, "mRNA": 2, "tRNA": 3, "rRNA": 4, "exon": 5, "CDS": 6}
            sorted_features = sorted(record.features, key=lambda f: feature_sort_order.get(f.type, 99))

            for feature_index, feature in enumerate(sorted_features):
                line = None
                if output_format == "gtf":
                    line = _process_feature_gtf(feature, seqname, "GenBank", feature_index, gtf_gene_ids, gtf_transcript_ids)
                elif output_format == "gff3":
                    # For GFF3, if the feature is 'source', it's handled as the sequence-region line
                    # and often duplicated as a "region" feature. We've handled the "region" explicitly above.
                    # Other features are processed here.
                    if feature.type == "source":
                        continue # Skip processing 'source' features here as they are handled by the 'region' line
                    line = _process_feature_gff3(
                        feature, seqname, "GenBank", feature_index,
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

    if 'file' not in request.files and 'genbank_data' not in request.form:
        return jsonify({"success": False, "error": "No GenBank data provided. Please upload a file or paste data."}), 400

    output_format = request.form.get('output_format', 'gff3').lower()
    if output_format not in ['gtf', 'gff3']:
        return jsonify({"success": False, "error": "Invalid output format. Choose 'gtf' or 'gff3'."}), 400

    genbank_data = ""
    try:
        if 'file' in request.files and request.files['file'].filename != '':
            genbank_file = request.files['file']
            if not genbank_file.filename.lower().endswith(('.gb', '.gbk', '.genbank')):
                return jsonify({"success": False, "error": "Invalid file type. Please upload a .gb, .gbk, or .genbank file."}), 400
            genbank_data = genbank_file.read().decode('utf-8')
        elif 'genbank_data' in request.form:
            genbank_data = request.form['genbank_data']
        else:
            return jsonify({"success": False, "error": "No GenBank data found in the request."}), 400

        if not genbank_data.strip():
            return jsonify({"success": False, "error": "Submitted GenBank data is empty."}), 400

    except Exception as e:
        logging.error(f"Error processing input: {e}", exc_info=True)
        return jsonify({"success": False, "error": f"Error processing input: {e}"}), 400

    converted_result, error_message = convert_genbank_data_to_annotation_format(genbank_data, output_format)

    if error_message:
        return jsonify({"success": False, "error": error_message}), 400

    unique_filename = f"{uuid.uuid4()}.{output_format}"
    output_path = os.path.join(TEMP_DIR, unique_filename)

    try:
        with open(output_path, "w") as f:
            f.write(converted_result)
    except IOError as e:
        logging.error(f"Error writing file to temporary storage: {e}", exc_info=True)
        return jsonify({"success": False, "error": "A server error occurred while preparing the file for download."}), 500

    download_url = request.host_url.rstrip('/') + f"/api/download/{unique_filename}"

    return jsonify({
        "success": True,
        "converted_data": converted_result,
        "download_url": download_url
    })

@app.route('/api/download/<filename>')
def download_converted_file(filename):
    """Serves the converted file for download."""
    from werkzeug.utils import secure_filename
    safe_filename = secure_filename(filename)

    if safe_filename != filename:
        return jsonify({"success": False, "error": "Invalid filename."}), 400

    full_file_path = os.path.join(TEMP_DIR, safe_filename)
    if not os.path.exists(full_file_path):
        logging.warning(f"Attempted to download non-existent file: {full_file_path}")
        return jsonify({"success": False, "error": "File not found or has expired."}), 404

    return send_from_directory(TEMP_DIR, safe_filename, as_attachment=True)

@app.errorhandler(413)
def request_entity_too_large(error):
    """Handles requests where the payload size exceeds MAX_CONTENT_LENGTH."""
    return jsonify({"success": False, "error": f"The file is too large. The maximum allowed size is {MAX_CONTENT_LENGTH / (1024 * 1024):.0f} MB."}), 413

@app.route('/')
def serve_index():
    """Simple health check endpoint."""
    return "GenBank to GTF/GFF3 Converter API is running. Use /api/convert to interact."

if __name__ == '__main__':
    # For local development, run on port 5000.
    # IMPORTANT: Set debug=False in production for security and performance.
    app.run(host='0.0.0.0', port=5000, debug=True)
