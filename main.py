import sys
import os
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, CompoundLocation
import logging
from urllib.parse import quote
from flask import Flask, request, jsonify # Import Flask

# Configure logging for better feedback
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

app = Flask(__name__) # Initialize Flask app

def _format_gff3_attribute_value(value: str) -> str:
    """
    Formats a string value for GFF3 attributes.
    Replaces internal double quotes with single quotes.
    URL-encodes problematic GFF3 characters (;,=,%,,,\t,\n,\r) but *preserves spaces*.
    """
    # Replace internal double quotes with single quotes to avoid parsing issues
    value = value.replace('"', "'")

    # URL-encode other problematic characters that are GFF3 delimiters or special,
    # but explicitly add space to 'safe' characters to prevent %20 encoding.
    # The 'safe' set includes alphanumeric characters, and common symbols like /:.-_
    # Adding ' ' to safe will prevent spaces from being encoded.
    return quote(value, safe="/:.,-_ ")

def _get_common_qualifiers(feature) -> dict:
    """
    Extracts common and frequently used qualifiers from a Biopython SeqFeature.
    Returns a dictionary with default None or empty list values if a qualifier is not found.
    """
    qualifiers = feature.qualifiers
    return {
        "gene": qualifiers.get("gene", [None])[0],
        "locus_tag": qualifiers.get("locus_tag", [None])[0],
        "product": qualifiers.get("product", [None])[0],
        "protein_id": qualifiers.get("protein_id", [None])[0],
        "transcript_id": qualifiers.get("transcript_id", [None])[0],
        "db_xrefs": qualifiers.get("db_xref", []),
        "note": qualifiers.get("note", [None])[0],
        "codon_start": int(qualifiers.get("codon_start", ["1"])[0]) if qualifiers.get("codon_start") else 1,
        "pseudo": "true" if "pseudo" in qualifiers else None,
        "exception": qualifiers.get("exception", [None])[0],
        "transl_table": int(qualifiers.get("transl_table", ["1"])[0]) if qualifiers.get("transl_table") else 1
    }

def _generate_fallback_id(prefix: str, seqname: str, feature_index: int, part_index: int = None) -> str:
    """
    Generates a unique fallback ID when no suitable identifier is available from GenBank qualifiers.
    """
    if part_index is not None:
        return f"{prefix}-{seqname}_{feature_index}_{part_index}"
    return f"{prefix}-{seqname}_{feature_index}"

def _get_partiality_attributes_gff3(feature_location, start, end) -> list:
    """
    Determines if a feature is partial at its start or end and returns
    the appropriate attributes for GFF3.
    """
    attributes = []
    is_partial_start = False
    is_partial_end = False

    if hasattr(feature_location, 'start_original') and feature_location.start_original.startswith("<"):
        is_partial_start = True
    if hasattr(feature_location, 'end_original') and feature_location.end_original.startswith(">"):
        is_partial_end = True

    if is_partial_start:
        attributes.append(f'partial=true')
        attributes.append(f'start_range=.,{start}')
    if is_partial_end:
        if 'partial=true' not in attributes:
            attributes.append(f'partial=true')
        attributes.append(f'end_range={end},.')
    return attributes

def _calculate_gff3_phase(cumulative_cds_length: int, initial_codon_start: int) -> str:
    """
    Calculates the GFF3 'phase' (frame) for a CDS segment.
    Phase is 0, 1, or 2, indicating the number of bases that should be
    removed from the beginning of this feature to reach the first base of the next codon.
    """
    phase = (initial_codon_start - 1 + cumulative_cds_length) % 3
    return str(phase)


def _process_feature_gtf(
    seqname: str,
    feature_index: int,
    feature,
    part_index: int,
    loc,
    common_ids: dict,
    cumulative_cds_length: int,
    gene_id_map: dict,
    gtf_mRNA_transcript_id_map: dict,
) -> tuple:
    """
    Processes a single Biopython SeqFeature (or part of a compound feature)
    and returns the GTF formatted line and updated cumulative CDS length.
    """
    start = int(loc.start) + 1
    end = int(loc.end)
    strand_char = "+" if loc.strand == 1 else ("-" if loc.strand == -1 else ".")
    score = "."
    frame = "."

    current_genbank_feature_type = feature.type
    genbank_gene_qual = common_ids["gene"]
    genbank_locus_tag_qual = common_ids["locus_tag"]
    genbank_product_qual = common_ids["product"]
    genbank_protein_id_qual = common_ids["protein_id"]
    genbank_transcript_id_qual = common_ids["transcript_id"]
    initial_codon_start = common_ids["codon_start"]

    if current_genbank_feature_type == "CDS":
        current_frame = (cumulative_cds_length + initial_codon_start - 1) % 3
        frame = str(current_frame)
        cumulative_cds_length += (end - start + 1)

    attributes = []
    gtf_feature_type_output = current_genbank_feature_type
    if gtf_feature_type_output == "mRNA":
        gtf_feature_type_output = "transcript" # GTF standard uses 'transcript'

    # Gene ID
    current_gene_identifier = genbank_gene_qual or genbank_locus_tag_qual
    if not current_gene_identifier:
        current_gene_identifier = _generate_fallback_id("gene", seqname, feature_index)

    if (seqname, current_gene_identifier) not in gene_id_map:
        gene_id_map[(seqname, current_gene_identifier)] = current_gene_identifier

    attributes.append(f'gene_id "{gene_id_map[(seqname, current_gene_identifier)]}"')

    # Transcript ID
    current_transcript_id = None

    if current_genbank_feature_type == "mRNA":
        transcript_id_base = genbank_transcript_id_qual or genbank_product_qual or genbank_gene_qual or genbank_locus_tag_qual
        if not transcript_id_base:
            transcript_id_base = _generate_fallback_id("transcript", seqname, feature_index)
        current_transcript_id = f"{transcript_id_base.replace(' ', '_')}"
        gtf_mRNA_transcript_id_map[(seqname, current_gene_identifier)] = current_transcript_id
    elif current_genbank_feature_type in ["CDS", "exon"]:
        if (seqname, current_gene_identifier) in gtf_mRNA_transcript_id_map:
            current_transcript_id = gtf_mRNA_transcript_id_map[(seqname, current_gene_identifier)]
        else:
            transcript_id_base = genbank_protein_id_qual or genbank_gene_qual or genbank_locus_tag_qual
            if not transcript_id_base:
                transcript_id_base = _generate_fallback_id("transcript", seqname, feature_index, part_index)
            current_transcript_id = f"{transcript_id_base.split('.')[0] if genbank_protein_id_qual else transcript_id_base}"
    elif current_genbank_feature_type in ["tRNA", "rRNA"]:
        transcript_id_base = genbank_gene_qual or genbank_locus_tag_qual or genbank_product_qual
        if not transcript_id_base:
            transcript_id_base = _generate_fallback_id("transcript", seqname, feature_index)
        current_transcript_id = f"{transcript_id_base}_{current_genbank_feature_type}"

    if not current_transcript_id:
        current_transcript_id = _generate_fallback_id("transcript", seqname, feature_index, part_index)

    attributes.append(f'transcript_id "{current_transcript_id}"')

    # Other GTF attributes
    if genbank_gene_qual:
        attributes.append(f'gene_name "{genbank_gene_qual}"')
    if genbank_product_qual:
        attributes.append(f'product "{genbank_product_qual.replace('"', '""')}"')
    if current_genbank_feature_type == "CDS" and genbank_protein_id_qual:
        attributes.append(f'protein_id "{genbank_protein_id_qual}"')
    if common_ids["pseudo"]:
        attributes.append(f'pseudo "true"')
    if common_ids["exception"]:
        attributes.append(f'exception "{common_ids["exception"].replace('"', '""')}"')
    if common_ids["transl_table"] and current_genbank_feature_type == "CDS":
        attributes.append(f'transl_table "{common_ids["transl_table"]}"')

    excluded_qualifiers = [
        "gene", "locus_tag", "transcript_id", "product", "protein_id",
        "codon_start", "translation", "db_xref", "note", "pseudo", "exception",
        "organism", "strain", "mol_type", "transl_table"
    ]
    for qual_key, qual_values in feature.qualifiers.items():
        if qual_key not in excluded_qualifiers:
            if qual_values and isinstance(qual_values, list) and qual_values[0]:
                escaped_value = qual_values[0].replace('"', '""')
                attributes.append(f'{qual_key} "{escaped_value}"')

    gtf_attributes_str = "; ".join(attributes) + ";"
    output_line = f"{seqname}\tGenbank\t{gtf_feature_type_output}\t{start}\t{end}\t{score}\t{strand_char}\t{frame}\t{gtf_attributes_str}\n"
    return output_line, cumulative_cds_length


def _process_feature_gff3(
    seqname: str,
    feature_index: int,
    feature,
    part_index: int,
    loc,
    common_ids: dict,
    cumulative_cds_length: int,
    gff3_id_map: dict,
) -> tuple:
    """
    Processes a single Biopython SeqFeature (or part of a compound feature)
    and returns the GFF3 formatted line and updated cumulative CDS length.
    """
    start = int(loc.start) + 1
    end = int(loc.end)
    strand_char = "+" if loc.strand == 1 else ("-" if loc.strand == -1 else ".")
    score = "."
    frame = "."

    current_genbank_feature_type = feature.type
    genbank_gene_qual = common_ids["gene"]
    genbank_locus_tag_qual = common_ids["locus_tag"]
    genbank_product_qual = common_ids["product"]
    genbank_protein_id_qual = common_ids["protein_id"]
    db_xrefs = common_ids["db_xrefs"]
    note = common_ids["note"]
    initial_codon_start = common_ids["codon_start"]

    if current_genbank_feature_type == "CDS":
        frame = _calculate_gff3_phase(cumulative_cds_length, initial_codon_start)
        cumulative_cds_length += (end - start + 1)

    attributes = []
    gff3_feature_type_output = current_genbank_feature_type
    gff3_source = "Genbank"

    # 1. ID (unique identifier for the feature)
    current_gff3_id = None
    id_base = None

    if gff3_feature_type_output == "gene":
        id_base = genbank_gene_qual or genbank_locus_tag_qual
        if id_base:
            current_gff3_id = f"gene-{_format_gff3_attribute_value(id_base.replace(' ', '_'))}"
        else:
            current_gff3_id = _generate_fallback_id("gene", seqname, feature_index)
    elif gff3_feature_type_output == "mRNA":
        id_base = genbank_gene_qual or genbank_locus_tag_qual or genbank_product_qual
        if id_base:
            current_gff3_id = f"rna-{_format_gff3_attribute_value(id_base.replace(' ', '_'))}"
        else:
            current_gff3_id = _generate_fallback_id("rna", seqname, feature_index)
    elif gff3_feature_type_output == "CDS":
        if genbank_protein_id_qual:
            current_gff3_id = f"cds-{_format_gff3_attribute_value(genbank_protein_id_qual)}"
        else:
            current_gff3_id = _generate_fallback_id("cds", seqname, feature_index, part_index)
    elif gff3_feature_type_output == "exon":
        mRNA_parent_id_key = (seqname, "mRNA", genbank_gene_qual or genbank_locus_tag_qual or genbank_product_qual)
        parent_mRNA_id = gff3_id_map.get(mRNA_parent_id_key)
        if parent_mRNA_id:
            current_gff3_id = f"exon-{_format_gff3_attribute_value(parent_mRNA_id.replace('rna-', ''))}-{part_index+1}"
        else:
            current_gff3_id = _generate_fallback_id("exon", seqname, feature_index, part_index)
    elif gff3_feature_type_output in ["tRNA", "rRNA"]:
        id_base = genbank_gene_qual or genbank_locus_tag_qual or genbank_product_qual
        if id_base:
            current_gff3_id = f"{gff3_feature_type_output}-{_format_gff3_attribute_value(id_base.replace(' ', '_'))}"
        else:
            current_gff3_id = _generate_fallback_id(gff3_feature_type_output, seqname, feature_index)

    if not current_gff3_id:
        current_gff3_id = _generate_fallback_id("feature", seqname, feature_index, part_index)

    attributes.append(f'ID={current_gff3_id}')
    gff3_id_map[(seqname, gff3_feature_type_output, current_gff3_id)] = current_gff3_id


    # 2. Parent attribute (for hierarchical features)
    parent_gff3_id = None
    if gff3_feature_type_output == "mRNA":
        parent_id_key = (seqname, "gene", genbank_gene_qual or genbank_locus_tag_qual)
        parent_gff3_id = gff3_id_map.get(parent_id_key)
    elif gff3_feature_type_output in ["exon", "CDS", "tRNA", "rRNA"]:
        mRNA_parent_id_key = (seqname, "mRNA", genbank_gene_qual or genbank_locus_tag_qual or genbank_product_qual)
        parent_gff3_id = gff3_id_map.get(mRNA_parent_id_key)

        if not parent_gff3_id and gff3_feature_type_output not in ["tRNA", "rRNA"]:
            gene_parent_id_key = (seqname, "gene", genbank_gene_qual or genbank_locus_tag_qual)
            parent_gff3_id = gff3_id_map.get(gene_parent_id_key)

    if parent_gff3_id:
        attributes.append(f'Parent={_format_gff3_attribute_value(parent_gff3_id)}')


    # 3. Name (often gene name, product, or protein_id)
    if genbank_gene_qual and gff3_feature_type_output == "gene":
        attributes.append(f'Name={_format_gff3_attribute_value(genbank_gene_qual)}')
    # NCBI example mRNA does NOT have a Name attribute.
    elif genbank_protein_id_qual and gff3_feature_type_output == "CDS":
        attributes.append(f'Name={_format_gff3_attribute_value(genbank_protein_id_qual)}')
    elif (genbank_product_qual or genbank_gene_qual) and gff3_feature_type_output in ["tRNA", "rRNA"]:
        # For tRNA/rRNA, use product or gene name as Name, replacing spaces with underscores for ID-like names
        attributes.append(f'Name={_format_gff3_attribute_value(genbank_product_qual or genbank_gene_qual)}')


    # 4. Dbxref (cross-references)
    for xref in db_xrefs:
        if current_genbank_feature_type == "CDS" and "protein_id" in xref and genbank_protein_id_qual:
            attributes.append(f'Dbxref=NCBI_GP:{_format_gff3_attribute_value(genbank_protein_id_qual)}')
        else:
            attributes.append(f'Dbxref={_format_gff3_attribute_value(xref)}')

    # 5. Note (general notes) - Exclude for CDS as per NCBI example
    if note and current_genbank_feature_type != "CDS":
        attributes.append(f'Note={_format_gff3_attribute_value(note)}')

    # 6. gbkey (NCBI specific feature key mapping)
    gbkey_map = {
        "gene": "Gene", "mRNA": "mRNA", "CDS": "CDS", "exon": "mRNA",
        "tRNA": "tRNA", "rRNA": "rRNA"
    }
    if gff3_feature_type_output in gbkey_map:
        attributes.append(f'gbkey={gbkey_map[gff3_feature_type_output]}')

    # 7. gene_biotype (NCBI specific for gene)
    if gff3_feature_type_output == "gene":
        attributes.append(f'gene_biotype=protein_coding')

    # 8. Pseudo attribute
    if common_ids["pseudo"]:
        attributes.append(f'pseudo=true')

    # 9. Exception attribute
    if common_ids["exception"]:
        attributes.append(f'exception={_format_gff3_attribute_value(common_ids["exception"])}')

    # 10. transl_table for CDS
    if common_ids["transl_table"] and current_genbank_feature_type == "CDS":
        attributes.append(f'transl_table={common_ids["transl_table"]}')

    # 11. Partial attributes (start_range, end_range, partial=true)
    attributes.extend(_get_partiality_attributes_gff3(feature.location, start, end))

    # 12. Add gene and product attributes to CDS features (as seen in NCBI example)
    if current_genbank_feature_type == "CDS":
        if genbank_gene_qual:
            attributes.append(f'gene={_format_gff3_attribute_value(genbank_gene_qual)}')
        if genbank_product_qual:
            attributes.append(f'product={_format_gff3_attribute_value(genbank_product_qual)}')
        if genbank_protein_id_qual:
            attributes.append(f'protein_id={_format_gff3_attribute_value(genbank_protein_id_qual)}')

    # 13. Other general qualifiers
    excluded_qualifiers = [
        "gene", "locus_tag", "transcript_id", "product", "protein_id",
        "codon_start", "translation", "db_xref", "note", "pseudo", "exception",
        "organism", "strain", "mol_type", "old-name", "transl_table"
    ]
    for qual_key, qual_values in feature.qualifiers.items():
        if qual_key not in excluded_qualifiers:
            if qual_values and isinstance(qual_values, list) and qual_values[0]:
                attributes.append(f'{qual_key}={_format_gff3_attribute_value(qual_values[0])}')

    gff3_attributes_str = ";".join(attributes)
    output_line = f"{seqname}\t{gff3_source}\t{gff3_feature_type_output}\t{start}\t{end}\t{score}\t{strand_char}\t{frame}\t{gff3_attributes_str}\n"
    return output_line, cumulative_cds_length


def convert_genbank_to_annotation_format(
    genbank_file_path: str,
    output_file_path: str,
    output_format: str = "gtf"
) -> None:
    """
    Converts a GenBank file to either GTF or GFF3 format.

    This function parses a GenBank file, extracts relevant features, and formats
    them into either GTF or GFF3 lines based on the specified output_format.
    It handles common GenBank qualifiers to populate attributes like gene_id,
    transcript_id (for GTF), or ID and Parent (for GFF3).

    Args:
        genbank_file_path (str): The path to the input GenBank file (.gb, .gbk).
        output_file_path (str): The path where the output annotation file will be saved.
        output_format (str): The desired output format ("gtf" (default) or "gff3").

    Raises:
        FileNotFoundError: If the input GenBank file does not exist.
        ValueError: If an unsupported output_format is specified.
        IOError: If there's an issue reading the GenBank file or writing
                 the output file.
        Exception: For any other unexpected errors during conversion.
    """
    if not os.path.exists(genbank_file_path):
        raise FileNotFoundError(f"Input GenBank file not found: {genbank_file_path}")

    output_format_lower = output_format.lower()
    if output_format_lower not in ["gtf", "gff3"]:
        raise ValueError(f"Unsupported output format: '{output_format}'. Choose 'gtf' or 'gff3'.")

    logging.info(f"Starting conversion of '{genbank_file_path}' to '{output_format_lower.upper()}' format...")

    try:
        with open(genbank_file_path, "r") as gb_handle, \
             open(output_file_path, "w") as out_handle:

            if output_format_lower == "gtf":
                out_handle.write("##gff-version 2.2\n")
            elif output_format_lower == "gff3":
                out_handle.write("##gff-version 3\n")
                out_handle.write("##gff-spec-version 1.21\n")

            gene_id_map = {}
            gtf_mRNA_transcript_id_map = {}
            gff3_id_map = {}

            for seq_record in SeqIO.parse(gb_handle, "genbank"):
                seqname = seq_record.id if seq_record.id else seq_record.name
                if not seqname:
                    logging.warning(f"Skipping a record with no ID or name in {genbank_file_path}")
                    continue

                logging.info(f"Processing record: {seqname} (Length: {len(seq_record.seq)} bp)")

                # --- GFF3 Specific Headers and Region Feature ---
                if output_format_lower == "gff3":
                    out_handle.write(f"##sequence-region {seqname} 1 {len(seq_record.seq)}\n")

                    source_feature_qualifiers = None
                    organism_qual = None
                    strain_qual = None
                    db_xref_qual = []

                    for feature in seq_record.features:
                        if feature.type == "source":
                            source_feature_qualifiers = feature.qualifiers
                            organism_qual = source_feature_qualifiers.get("organism", [None])[0]
                            strain_qual = source_feature_qualifiers.get("strain", [None])[0]
                            db_xref_qual = source_feature_qualifiers.get("db_xref", [])
                            break

                    for xref in db_xref_qual:
                        if "taxon:" in xref:
                            out_handle.write(f"##species https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id={xref.split(':')[-1]}\n")
                            break

                    region_id = f"{seqname}:1..{len(seq_record.seq)}"
                    gff3_id_map[(seqname, "region", region_id)] = region_id

                    region_attributes = [
                        f'ID={_format_gff3_attribute_value(region_id)}',
                        'gbkey=Src',
                        'mol_type=genomic DNA'
                    ]
                    if organism_qual:
                        region_attributes.append(f'old-name={_format_gff3_attribute_value(organism_qual)}')
                    if strain_qual:
                        region_attributes.append(f'strain={_format_gff3_attribute_value(strain_qual)}')
                    for xref in db_xref_qual:
                        region_attributes.append(f'Dbxref={_format_gff3_attribute_value(xref)}')

                    region_attributes_str = ";".join(region_attributes)
                    out_handle.write(f"{seqname}\tGenbank\tregion\t1\t{len(seq_record.seq)}\t.\t+\t.\t{region_attributes_str}\n")


                # --- Pre-processing for GFF3: Identify and store IDs for parent features (gene, mRNA) ---
                for feature_index, feature in enumerate(seq_record.features):
                    if feature.type not in ["gene", "mRNA", "CDS", "tRNA", "rRNA", "exon"]:
                        continue

                    common_ids = _get_common_qualifiers(feature)
                    genbank_gene_qual = common_ids["gene"]
                    genbank_locus_tag_qual = common_ids["locus_tag"]
                    genbank_product_qual = common_ids["product"]

                    if output_format_lower == "gff3":
                        if feature.type == "gene":
                            gene_id_base = genbank_gene_qual or genbank_locus_tag_qual
                            gene_gff3_id = f"gene-{_format_gff3_attribute_value(gene_id_base.replace(' ', '_'))}" if gene_id_base else _generate_fallback_id("gene", seqname, feature_index)
                            gff3_id_map[(seqname, "gene", gene_id_base or gene_gff3_id)] = gene_gff3_id

                        elif feature.type == "mRNA":
                            mRNA_id_base = genbank_gene_qual or genbank_locus_tag_qual or genbank_product_qual
                            mRNA_gff3_id = f"rna-{_format_gff3_attribute_value(mRNA_id_base.replace(' ', '_'))}" if mRNA_id_base else _generate_fallback_id("rna", seqname, feature_index)
                            gff3_id_map[(seqname, "mRNA", mRNA_id_base or mRNA_gff3_id)] = mRNA_gff3_id

                # Second pass: Process all features, including individual parts for CDS/Exon
                sorted_features_for_output = []
                for feature_index, feature in enumerate(seq_record.features):
                    if feature.type not in ["gene", "mRNA", "CDS", "tRNA", "rRNA", "exon"]:
                        continue
                    sorted_features_for_output.append((feature_index, feature))

                sorted_features_for_output.sort(key=lambda x: (
                    0 if x[1].type == "gene" else
                    1 if x[1].type == "mRNA" else
                    2 if x[1].type == "exon" else
                    3 if x[1].type == "CDS" else
                    4,
                    x[1].location.start
                ))

                for feature_index, feature in sorted_features_for_output:
                    common_ids = _get_common_qualifiers(feature)
                    strand_char = "+" if feature.location.strand == 1 else ("-" if feature.location.strand == -1 else ".")

                    locations = []
                    if isinstance(feature.location, CompoundLocation):
                        locations.extend(feature.location.parts)
                    else:
                        locations.append(feature.location)

                    cumulative_cds_length = 0

                    if output_format_lower == "gff3":
                        if feature.type == "gene":
                            gene_id_base = common_ids["gene"] or common_ids["locus_tag"]
                            gene_gff3_id = gff3_id_map.get((seqname, "gene", gene_id_base))

                            gene_start = int(feature.location.start) + 1
                            gene_end = int(feature.location.end)
                            if isinstance(feature.location, CompoundLocation):
                                gene_start = min([int(p.start)+1 for p in feature.location.parts])
                                gene_end = max([int(p.end) for p in feature.location.parts])

                            attributes = [f'ID={gene_gff3_id}']
                            if common_ids["gene"]:
                                attributes.append(f'Name={_format_gff3_attribute_value(common_ids["gene"])}')
                            attributes.append(f'gbkey=Gene')
                            attributes.append(f'gene_biotype=protein_coding')
                            attributes.extend(_get_partiality_attributes_gff3(feature.location, gene_start, gene_end))
                            if common_ids["pseudo"]: attributes.append(f'pseudo=true')
                            if common_ids["exception"]: attributes.append(f'exception={_format_gff3_attribute_value(common_ids["exception"])}')

                            for xref in common_ids["db_xrefs"]:
                                attributes.append(f'Dbxref={_format_gff3_attribute_value(xref)}')

                            out_handle.write(f"{seqname}\tGenbank\tgene\t{gene_start}\t{gene_end}\t.\t{strand_char}\t.\t{';'.join(attributes)}\n")
                            continue

                        elif feature.type == "mRNA":
                            mRNA_id_base = common_ids["gene"] or common_ids["locus_tag"] or common_ids["product"]
                            mRNA_gff3_id = gff3_id_map.get((seqname, "mRNA", mRNA_id_base))

                            mRNA_start = int(feature.location.start) + 1
                            mRNA_end = int(feature.location.end)
                            if isinstance(feature.location, CompoundLocation):
                                mRNA_start = min([int(p.start)+1 for p in feature.location.parts])
                                mRNA_end = max([int(p.end) for p in feature.location.parts])

                            attributes = [f'ID={mRNA_gff3_id}']

                            parent_id_key = (seqname, "gene", common_ids["gene"] or common_ids["locus_tag"])
                            parent_gff3_id = gff3_id_map.get(parent_id_key)
                            if parent_gff3_id:
                                attributes.append(f'Parent={_format_gff3_attribute_value(parent_gff3_id)}')

                            # NCBI example mRNA does NOT have a Name attribute.
                            attributes.append(f'gbkey=mRNA')
                            attributes.extend(_get_partiality_attributes_gff3(feature.location, mRNA_start, mRNA_end))
                            if common_ids["gene"]: attributes.append(f'gene={_format_gff3_attribute_value(common_ids["gene"])}')
                            if common_ids["product"]: attributes.append(f'product={_format_gff3_attribute_value(common_ids["product"])}')
                            if common_ids["pseudo"]: attributes.append(f'pseudo=true')
                            if common_ids["exception"]: attributes.append(f'exception={_format_gff3_attribute_value(common_ids["exception"])}')

                            for xref in common_ids["db_xrefs"]:
                                attributes.append(f'Dbxref={_format_gff3_attribute_value(xref)}')

                            out_handle.write(f"{seqname}\tGenbank\tmRNA\t{mRNA_start}\t{mRNA_end}\t.\t{strand_char}\t.\t{';'.join(attributes)}\n")
                            continue

                        elif feature.type in ["tRNA", "rRNA"]:
                            current_gff3_id = f"{feature.type}-{_format_gff3_attribute_value((common_ids['gene'] or common_ids['locus_tag'] or seqname).replace(' ', '_'))}_{feature_index}"
                            attributes = [f'ID={current_gff3_id}']

                            parent_id_key = (seqname, "gene", common_ids["gene"] or common_ids["locus_tag"])
                            parent_gff3_id = gff3_id_map.get(parent_id_key)
                            if parent_gff3_id:
                                attributes.append(f'Parent={_format_gff3_attribute_value(parent_gff3_id)}')

                            if common_ids["product"]:
                                attributes.append(f'Name={_format_gff3_attribute_value(common_ids["product"])}')
                            attributes.append(f'gbkey={feature.type}')
                            attributes.extend(_get_partiality_attributes_gff3(feature.location, int(feature.location.start)+1, int(feature.location.end)))
                            if common_ids["pseudo"]: attributes.append(f'pseudo=true')
                            if common_ids["exception"]: attributes.append(f'exception={_format_gff3_attribute_value(common_ids["exception"])}')
                            for xref in common_ids["db_xrefs"]:
                                attributes.append(f'Dbxref={_format_gff3_attribute_value(xref)}')

                            out_handle.write(f"{seqname}\tGenbank\t{feature.type}\t{int(feature.location.start)+1}\t{int(feature.location.end)}\t.\t{strand_char}\t.\t{';'.join(attributes)}\n")
                            continue


                    for part_index, loc in enumerate(locations):
                        output_line = ""
                        if output_format_lower == "gtf":
                            output_line, cumulative_cds_length = _process_feature_gtf(
                                seqname, feature_index, feature, part_index, loc,
                                common_ids, cumulative_cds_length, gene_id_map,
                                gtf_mRNA_transcript_id_map
                            )
                            out_handle.write(output_line)
                        elif output_format_lower == "gff3":
                            if feature.type == "exon":
                                current_gff3_id = f"exon-{_format_gff3_attribute_value((common_ids['gene'] or common_ids['locus_tag'] or seqname).replace(' ', '_'))}-{part_index+1}"
                                attributes = [f'ID={current_gff3_id}']

                                mRNA_parent_id_key = (seqname, "mRNA", common_ids["gene"] or common_ids["locus_tag"] or common_ids["product"])
                                parent_mRNA_id = gff3_id_map.get(mRNA_parent_id_key)
                                if parent_mRNA_id:
                                    attributes.append(f'Parent={_format_gff3_attribute_value(parent_mRNA_id)}')

                                attributes.append(f'gbkey=mRNA')
                                attributes.extend(_get_partiality_attributes_gff3(loc, int(loc.start)+1, int(loc.end)))
                                if common_ids["gene"]: attributes.append(f'gene={_format_gff3_attribute_value(common_ids["gene"])}')
                                if common_ids["product"]: attributes.append(f'product={_format_gff3_attribute_value(common_ids["product"])}')
                                if common_ids["pseudo"]: attributes.append(f'pseudo=true')
                                if common_ids["exception"]: attributes.append(f'exception={_format_gff3_attribute_value(common_ids["exception"])}')

                                for xref in common_ids["db_xrefs"]:
                                    attributes.append(f'Dbxref={_format_gff3_attribute_value(xref)}')

                                output_line = f"{seqname}\tGenbank\texon\t{int(loc.start)+1}\t{int(loc.end)}\t.\t{strand_char}\t.\t{';'.join(attributes)}\n"
                                out_handle.write(output_line)

                            elif feature.type == "CDS":
                                output_line, cumulative_cds_length = _process_feature_gff3(
                                    seqname, feature_index, feature, part_index, loc,
                                    common_ids, cumulative_cds_length, gff3_id_map
                                )
                                out_handle.write(output_line)

        logging.info("Conversion complete!")

    except FileNotFoundError as e:
        logging.error(f"Error: {e}")
        raise
    except ValueError as e:
        logging.error(f"Configuration error: {e}")
        raise
    except IOError as e:
        logging.error(f"File I/O error during conversion: {e}")
        raise
    except Exception as e:
        logging.error(f"An unexpected error occurred: {e}")
        raise

# API Endpoint
@app.route('/convert', methods=['POST'])
def convert_file():
    if 'genbank_file' not in request.files:
        return jsonify({"error": "No genbank_file part in the request"}), 400

    genbank_file = request.files['genbank_file']
    output_format = request.form.get('format', 'gtf').lower()

    if output_format not in ['gtf', 'gff3']:
        return jsonify({"error": "Invalid output format. Choose 'gtf' or 'gff3'."}), 400

    if genbank_file.filename == '':
        return jsonify({"error": "No selected file"}), 400

    if genbank_file:
        # Save the uploaded file temporarily
        input_filepath = "/tmp/input.gb"
        output_filepath = "/tmp/output.annotation"
        genbank_file.save(input_filepath)

        try:
            convert_genbank_to_annotation_format(input_filepath, output_filepath, output_format)
            with open(output_filepath, "r") as f:
                output_content = f.read()
            return output_content, 200, {'Content-Type': 'text/plain'} # Return as plain text

        except Exception as e:
            logging.error(f"Conversion API error: {e}")
            return jsonify({"error": str(e)}), 500
        finally:
            # Clean up temporary files
            if os.path.exists(input_filepath):
                os.remove(input_filepath)
            if os.path.exists(output_filepath):
                os.remove(output_filepath)

# To run the Flask app
if __name__ == "__main__":
    # Cloud Run expects the application to listen on PORT environment variable
    port = int(os.environ.get("PORT", 8080))
    app.run(host="0.0.0.0", port=port)
