import argparse
import logging
import pandas as pd
from collections import defaultdict
from Bio import SeqIO
from typing import Dict, Set, List, Tuple, Optional

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def process_species_name_from_fasta_header(header: str) -> Optional[str]:
    """Extract a normalized species name from a FASTA header string.

    Attempts to parse a "Genus species" pattern from various FASTA header formats.

    Args:
        header (str): FASTA header line (with or without leading '>').

    Returns:
        Optional[str]: The extracted species name (e.g., "Genus species"), or None if parsing fails.
    """
    try:
        if header.startswith('>'):
            header = header[1:]

        parts = header.split('.')[0].split('_')
        if len(parts) >= 2:
            species_candidate = " ".join(parts[:2])
            genus_species = []
            temp_name = ""
            for char in species_candidate:
                if char.isalpha() or char == ' ':
                    temp_name += char
                else:
                    break
            if temp_name.strip():
                name_parts = temp_name.strip().split()
                if len(name_parts) == 2 and name_parts[0][0].isupper() and name_parts[1].islower():
                    return " ".join(name_parts)
                elif len(name_parts) == 1 and name_parts[0][0].isupper():
                    simple_name = header.split(' ')[0].replace('_', ' ').strip()
                    simple_name_parts = simple_name.split()
                    if len(simple_name_parts) >= 2:
                        return " ".join(simple_name_parts[:2])
                    return simple_name_parts[0] if simple_name_parts else None

        simple_name = header.split(' ')[0].replace('_', ' ').strip()
        simple_name_parts = simple_name.split()
        if len(simple_name_parts) >= 2 :
            return " ".join(simple_name_parts[:2])
        elif simple_name_parts:
            return simple_name_parts[0]

    except Exception as e:
        logging.warning(f"Could not parse species name from header '{header}': {e}")
    return None

def load_taxonomy_data(
    names_file: str, nodes_file: str
) -> Tuple[Dict[str, Dict[str, str]], Dict[str, str], Dict[str, List[str]]]:
    """Load taxonomy data from NCBI names.dmp and nodes.dmp files.

    Builds a unified map of tax IDs to their name, rank, and parent ID.
    Also creates a mapping of species names to tax IDs and a parent-to-children hierarchy.

    Args:
        names_file (str): Path to NCBI names.dmp file.
        nodes_file (str): Path to NCBI nodes.dmp file.

    Returns:
        Tuple[Dict[str, Dict[str, str]], Dict[str, str], Dict[str, List[str]]]: A tuple containing:
            - tax_id_to_info: Mapping of tax_id to {"parent_id", "rank", "name"}.
            - species_name_to_tax_id: Mapping of normalized species names to tax IDs.
            - parent_to_children: Mapping of parent tax IDs to lists of child tax IDs.
    """
    logging.info("Loading taxonomy data...")
    tax_id_to_info: Dict[str, Dict[str, str]] = {}
    species_name_to_tax_id: Dict[str, str] = {}
    parent_to_children: Dict[str, List[str]] = defaultdict(list)

    logging.info(f"Parsing nodes file: {nodes_file}...")
    try:
        with open(nodes_file, 'r') as f_nodes:
            for line_num, line in enumerate(f_nodes):
                parts = line.strip().split('\t|\t')
                if len(parts) < 3:
                    logging.warning(f"Malformed line {line_num + 1} in {nodes_file} (expected at least 3 parts separated by '|'): {line.strip()}")
                    continue

                tax_id = parts[0].strip()
                parent_id = parts[1].strip()
                rank_from_nodes = parts[2].strip().lower()

                tax_id_to_info[tax_id] = {
                    "parent_id": parent_id,
                    "rank": rank_from_nodes,
                    "name": ""
                }

                if parent_id != tax_id and parent_id != "0":
                    parent_to_children[parent_id].append(tax_id)

    except FileNotFoundError:
        logging.error(f"Nodes file not found: {nodes_file}")
        raise
    logging.info(f"Loaded basic hierarchy and ranks for {len(tax_id_to_info)} tax_ids from {nodes_file}.")

    logging.info(f"Parsing names file: {names_file}...")
    names_assigned_from_file = 0
    species_names_mapped = 0
    try:
        with open(names_file, 'r') as f_names:
            for line_num, line in enumerate(f_names):
                parts = line.strip().split('\t')
                if len(parts) < 3:
                    logging.warning(f"Malformed line {line_num + 1} in {names_file} (expected 3 tab-separated parts): {line.strip()}")
                    continue

                current_tax_id = parts[0].strip()
                rank_from_name_file = parts[1].strip().lower()
                name_txt = parts[2].strip()

                if current_tax_id in tax_id_to_info:
                    tax_id_to_info[current_tax_id]["name"] = name_txt
                    names_assigned_from_file += 1

                    if rank_from_name_file == "species":
                        norm_name_parts = name_txt.split()[:2]
                        if len(norm_name_parts) == 2:
                            normalized_name = " ".join(norm_name_parts)
                            species_name_to_tax_id[normalized_name] = current_tax_id
                            species_names_mapped += 1
                        else:
                            logging.debug(f"Species name '{name_txt}' for TaxID {current_tax_id} from {names_file} line {line_num+1} could not be normalized to 'Genus species' format. Not added to species_name_to_tax_id map.")
                else:
                    logging.debug(f"TaxID {current_tax_id} ('{name_txt}') from {names_file} line {line_num+1} was not found in pre-loaded data from {nodes_file}. It will be ignored for lineage building.")

    except FileNotFoundError:
        logging.error(f"Names file not found: {names_file}")
        raise
    logging.info(f"Assigned names to {names_assigned_from_file} tax_ids using {names_file}.")
    logging.info(f"Mapped {species_names_mapped} species names to TaxIDs for FASTA matching.")

    unnamed_count = 0
    for tax_id, info in tax_id_to_info.items():
        if not info.get("name") or info["name"] == "":
            info["name"] = f"Unnamed_taxid_{tax_id}"
            unnamed_count += 1
    if unnamed_count > 0:
        logging.warning(f"{unnamed_count} tax_ids (from {nodes_file}) did not receive a name from {names_file} and were assigned placeholder names.")

    logging.info(f"Taxonomy data loading complete. Total TaxIDs in info map: {len(tax_id_to_info)}.")
    return tax_id_to_info, species_name_to_tax_id, parent_to_children

def get_lineage_for_taxid(tax_id: str, tax_id_to_info: Dict[str, Dict[str, str]]) -> List[Tuple[str, str, str]]:
    """Build the full taxonomic lineage for a given tax ID.

    Traverses parent pointers up to the root, returning the lineage
    ordered from highest (root) to lowest (species).

    Args:
        tax_id (str): The taxonomic ID to build lineage for.
        tax_id_to_info (Dict[str, Dict[str, str]]): Mapping of tax_id to metadata.

    Returns:
        List[Tuple[str, str, str]]: List of (name, tax_id, rank) tuples ordered from root to terminal.
    """
    lineage = []
    current_tax_id = tax_id
    visited_ids_for_this_lineage = set()

    while current_tax_id and current_tax_id in tax_id_to_info and current_tax_id not in visited_ids_for_this_lineage:
        visited_ids_for_this_lineage.add(current_tax_id)
        info = tax_id_to_info[current_tax_id]
        name = info.get("name", f"Unnamed_{current_tax_id}")
        rank = info.get("rank", "no rank")
        lineage.append((name, current_tax_id, rank))

        parent_id = info.get("parent_id")
        if not parent_id or parent_id == current_tax_id or parent_id == "0" or parent_id == "1":
            break
        current_tax_id = parent_id

    return lineage[::-1]


def filter_hits(input_file, hit_threshold):
    """Filter rows with HitRatio greater than hit_threshold from the input TSV file.

    Args:
        input_file (str): Path to the input TSV file.
        hit_threshold (float): Threshold value for HitRatio filtering.

    Returns:
        pd.DataFrame: Filtered DataFrame, or an empty DataFrame if column is missing.
    """
    logging.info(f"Filtering hits from {input_file} with threshold {hit_threshold}...")
    df = pd.read_csv(input_file, sep='\t')
    if 'HitRatio' not in df.columns:
        logging.error(f"'HitRatio' column not found in {input_file}. Columns are: {df.columns.tolist()}")
        return pd.DataFrame(columns=df.columns)

    filtered_df = df[df['HitRatio'] > hit_threshold]
    logging.info(f"Filtered to {len(filtered_df)} rows.")
    return filtered_df


def classify_probes(filtered_df, classification_level_rank: str):
    """Classify probes into groups based on the specified taxonomic rank.

    Groups probes by their PhyloLevel name where the corresponding rank matches
    the specified classification level.

    Args:
        filtered_df (pd.DataFrame): DataFrame with 'Rank' and 'PhyloLevel' columns.
        classification_level_rank (str): The taxonomic rank to classify by (e.g., 'class', 'order').

    Returns:
        defaultdict: Mapping of PhyloLevel names to lists of ProbeIDs.
    """
    logging.info(f"Classifying probes by PhyloLevel names whose rank is '{classification_level_rank}'...")
    probe_sets = defaultdict(list)

    for _, row in filtered_df.iterrows():
        if row['Rank'].lower() == classification_level_rank.lower():
            probe_sets[row['PhyloLevel']].append(row['ProbeID'])

    logging.info(f"Classification results: {len(probe_sets)} groups based on rank '{classification_level_rank}'.")
    for level_name, probes in probe_sets.items():
        logging.debug(f"  Group '{level_name}' (rank: {classification_level_rank}) has {len(probes)} probes.")
    return probe_sets


def determine_species_and_lineages(
    fasta_file: str,
    probe_sets: Dict[str, List[str]],
    species_name_to_tax_id: Dict[str, str],
    tax_id_to_info: Dict[str, Dict[str, str]]
) -> Tuple[Dict[str, Set[str]], Dict[str, int], Dict[str, str]]:
    """Determine species and lineages from FASTA, mapping sequences to probe set PhyloLevels.

    Parses FASTA headers for species names, resolves their tax IDs and lineages,
    and associates each sequence with the target PhyloLevel groups from probe_sets.

    Args:
        fasta_file (str): Path to the background FASTA file.
        probe_sets (Dict[str, List[str]]): Mapping of PhyloLevel names to lists of ProbeIDs.
        species_name_to_tax_id (Dict[str, str]): Mapping of species names to tax IDs.
        tax_id_to_info (Dict[str, Dict[str, str]]): Mapping of tax_id to metadata.

    Returns:
        Tuple[Dict[str, Set[str]], Dict[str, int], Dict[str, str]]: A tuple containing:
            - species_dict: Mapping of PhyloLevel names to sets of associated sequence IDs.
            - sequence_lengths: Mapping of sequence IDs to their lengths.
            - full_sequences: Mapping of sequence IDs to full sequence strings.
    """
    logging.info("Determining species, TaxIDs, lineages from FASTA and mapping to probe set PhyloLevels...")

    species_dict = defaultdict(set)
    sequence_lengths = {}
    full_sequences = defaultdict(str)

    target_phylo_level_names = set(probe_sets.keys())
    logging.info(f"Target PhyloLevel names from probe_sets: {target_phylo_level_names}")

    processed_fasta_entries = 0
    matched_to_target_phylolevel = 0

    for record in SeqIO.parse(fasta_file, "fasta"):
        processed_fasta_entries += 1
        current_seq_id = record.id

        fasta_header_species_name = process_species_name_from_fasta_header(record.description)

        if not fasta_header_species_name:
            logging.debug(f"Could not parse species name for sequence ID '{current_seq_id}'. Header: '{record.description}'")
            continue

        species_tax_id = species_name_to_tax_id.get(fasta_header_species_name)
        if not species_tax_id:
            simple_parsed_name = " ".join(fasta_header_species_name.split()[:2])
            species_tax_id = species_name_to_tax_id.get(simple_parsed_name)

        if not species_tax_id:
            logging.debug(f"No TaxID found for parsed species name '{fasta_header_species_name}' (from seq '{current_seq_id}').")
            continue

        lineage = get_lineage_for_taxid(species_tax_id, tax_id_to_info)
        if not lineage:
            logging.debug(f"No lineage found for TaxID '{species_tax_id}' (species: '{fasta_header_species_name}', seq: '{current_seq_id}').")
            continue

        found_match = False
        for phylo_name_in_lineage, _, _ in lineage:
            if phylo_name_in_lineage in target_phylo_level_names:
                species_dict[phylo_name_in_lineage].add(current_seq_id)
                if current_seq_id not in sequence_lengths:
                    sequence_lengths[current_seq_id] = 0
                    full_sequences[current_seq_id] = ""
                found_match = True

        if found_match:
            matched_to_target_phylolevel +=1
            sequence_lengths[current_seq_id] = len(record.seq)
            full_sequences[current_seq_id] = str(record.seq)

    logging.info(f"Processed {processed_fasta_entries} FASTA entries.")
    logging.info(f"Associated {matched_to_target_phylolevel} FASTA entries with target PhyloLevels based on lineage.")
    logging.info(f"Final species_dict contains {len(species_dict)} PhyloLevel groups.")
    for level, seq_ids in species_dict.items():
        logging.debug(f"  PhyloLevel '{level}' has {len(seq_ids)} associated FASTA sequences.")
    if not any(species_dict.values()):
         logging.warning("species_dict is empty or all its lists are empty. This will lead to no coverage calculation or empty output.tsv.")
         logging.warning("Check: 1. FASTA headers and parsing (process_species_name_from_fasta_header).")
         logging.warning("       2. species_name_to_tax_id mapping (names.dmp, species names).")
         logging.warning("       3. Lineage construction and matching against target_phylo_level_names.")

    return species_dict, sequence_lengths, full_sequences


def calculate_coverage(probe_sets, blast_file, sequence_lengths):
    """Calculate the coverage of probe sets within the classification level.

    Reads BLAST results, maps probes to their target sequences, and computes
    per-sequence coverage percentages for each phylogenetic level.

    Args:
        probe_sets (dict): Mapping of PhyloLevel names to lists of ProbeIDs.
        blast_file (str): Path to the BLAST results file.
        sequence_lengths (dict): Mapping of target sequence IDs to their lengths.

    Returns:
        tuple: A tuple containing:
            - coverage_results: List of (PhyloLevel, target_id, coverage_percent) tuples.
            - empty_dict: An empty dictionary (placeholder, retained for compatibility).
    """
    logging.info("Calculating coverage...")
    coverage_results = []
    covered_regions = defaultdict(lambda: defaultdict(set))

    with open(blast_file, 'r') as f:
        for line in f:
            columns = line.strip().split()
            if len(columns) < 10:
                logging.warning(f"Skipping malformed BLAST line (less than 10 columns): {line.strip()}")
                continue

            probe_id = columns[0]
            target_id = columns[1]

            try:
                sstart = int(columns[8])
                send = int(columns[9])
            except ValueError:
                logging.warning(f"Could not parse sstart/send from BLAST line: {line.strip()}")
                continue

            probe_phylo_level = None
            for phylo_level, probes_in_level in probe_sets.items():
                if probe_id in probes_in_level:
                    probe_phylo_level = phylo_level
                    break

            if probe_phylo_level and target_id in sequence_lengths:
                start_pos = min(sstart, send)
                end_pos = max(sstart, send)
                for pos in range(start_pos, end_pos + 1):
                    covered_regions[probe_phylo_level][target_id].add(pos)

    for phylo_level, targets_in_level in covered_regions.items():
        for target_id, covered_positions in targets_in_level.items():
            if target_id in sequence_lengths and sequence_lengths[target_id] > 0:
                total_covered_bases = len(covered_positions)
                coverage_percent = (total_covered_bases / sequence_lengths[target_id]) * 100
                coverage_results.append((phylo_level, target_id, coverage_percent))
            elif target_id not in sequence_lengths:
                 logging.warning(f"Target ID '{target_id}' from BLAST results not found in sequence_lengths (derived from FASTA). Skipping coverage calculation for it.")
            elif sequence_lengths[target_id] == 0:
                 logging.warning(f"Target ID '{target_id}' has zero length in sequence_lengths. Skipping coverage calculation.")

    logging.info(f"Coverage calculation completed. Found coverage data for {len(coverage_results)} (phylo_level, target) pairs.")
    return coverage_results, {}


def save_probes_by_phylo_level(probe_sets, blast_file, output_probe_prefix, probe_fasta_file):
    """Save probe sequences (FASTA) belonging to specified PhyloLevels into separate files.

    Args:
        probe_sets (dict): Mapping of PhyloLevel names to lists of ProbeIDs.
        blast_file (str): Path to the BLAST results file.
        output_probe_prefix (str): Prefix for output FASTA filenames.
        probe_fasta_file (str): Path to the FASTA file containing all probe sequences.
    """
    logging.info("Writing probe sequences to separate FASTA files by PhyloLevel...")

    probe_sequences_db = {}
    try:
        for record in SeqIO.parse(probe_fasta_file, "fasta"):
            probe_sequences_db[record.id] = record
    except FileNotFoundError:
        logging.error(f"Probe FASTA file not found: {probe_fasta_file}")
        return

    if not probe_sequences_db:
        logging.warning(f"No sequences loaded from probe FASTA file: {probe_fasta_file}")
        return

    for level_name, probe_id_list in probe_sets.items():
        filename = f"{output_probe_prefix}_{level_name.replace(' ', '_').replace('/', '_')}_probes.fasta"
        probes_to_write_for_this_level = []
        for probe_id in sorted(list(set(probe_id_list))):
            if probe_id in probe_sequences_db:
                probes_to_write_for_this_level.append(probe_sequences_db[probe_id])
            else:
                logging.warning(f"ProbeID '{probe_id}' (for PhyloLevel '{level_name}') not found in probe FASTA file '{probe_fasta_file}'.")

        if probes_to_write_for_this_level:
            with open(filename, 'w') as probe_out:
                SeqIO.write(probes_to_write_for_this_level, probe_out, "fasta")
            logging.info(f"Probe FASTA file for {level_name} written to {filename} with {len(probes_to_write_for_this_level)} sequences.")
        else:
            logging.info(f"No probe sequences to write for PhyloLevel '{level_name}' (either no probes in set or probes not in FASTA db).")


def main():
    """Parse command-line arguments and run the enhanced grouped coverage analysis workflow."""
    parser = argparse.ArgumentParser(description='Analyze probe data and calculate coverage using NCBI taxonomy.')
    parser.add_argument('filter_output_file', type=str, help='Input TSV file from previous filtering (e.g., probe_filtered.tsv). Columns: ProbeID, PhyloLevel, TaxID, Rank, Total, Hit, HitRatio, SelfSpeciesName')
    parser.add_argument('hit_threshold', type=float, help='HitRatio threshold (e.g., 0.0 to filter only by rank from filter_output_file, or higher to pre-filter probes)')
    parser.add_argument('classification_level_rank', type=str, help='The taxonomic *Rank* to classify probes by (e.g., "class", "order"). Probes will be grouped by their PhyloLevel NAME if that name corresponds to this rank.')
    parser.add_argument('fasta_file', type=str, help='Background FASTA file path (sequences of target species). Headers should be parsable for species names.')
    parser.add_argument('blast_file', type=str, help='BLAST results file path (probe vs background FASTA).')
    parser.add_argument('names_dmp_file', type=str, help='Path to NCBI names.dmp file.')
    parser.add_argument('nodes_dmp_file', type=str, help='Path to NCBI nodes.dmp file.')
    parser.add_argument('probe_fasta_file', type=str, help='FASTA file containing all probe sequences.')
    parser.add_argument('output_tsv_file', type=str, help='Output TSV file path for coverage summary.')
    parser.add_argument('output_probe_prefix', type=str, help='Output prefix for probe FASTA files per PhyloLevel.')

    args = parser.parse_args()
    logging.info("Started processing with arguments: %s", args)

    tax_id_to_info, species_name_to_tax_id, _ = load_taxonomy_data(args.names_dmp_file, args.nodes_dmp_file)
    if not tax_id_to_info or not species_name_to_tax_id:
        logging.error("Taxonomy data could not be loaded. Exiting.")
        return

    initial_filtered_df = filter_hits(args.filter_output_file, args.hit_threshold)

    probe_sets = {}
    if initial_filtered_df.empty:
        logging.warning(f"No probes passed initial HitRatio > {args.hit_threshold} from {args.filter_output_file}.")
    else:
        probe_sets = classify_probes(initial_filtered_df, args.classification_level_rank)

    if not probe_sets:
        logging.warning(f"No probe groups were formed (either no initial probes passed filter, or no filtered probes matched rank '{args.classification_level_rank}').")
        logging.warning(f"The output TSV file ({args.output_tsv_file}) will not be created.")
    else:
        species_dict, sequence_lengths, _ = determine_species_and_lineages(
            args.fasta_file,
            probe_sets,
            species_name_to_tax_id,
            tax_id_to_info
        )

        if not species_dict or not any(species_dict.values()):
            logging.warning("No sequences from FASTA file could be associated with the target PhyloLevel groups from probe_sets.")
            logging.warning("Coverage calculation will likely yield no results. The output TSV will be generated but may contain NA values for coverage.")

        coverage_results, _ = calculate_coverage(probe_sets, args.blast_file, sequence_lengths)

        logging.info(f"Writing final summary results to {args.output_tsv_file}...")
        with open(args.output_tsv_file, 'w') as output:
            output.write("PhyloLevelName\tTotalProbeCountInGroup\tCoverageQ1\tCoverageQ2\tCoverageQ3\n")

            for phylo_level_name, probes_in_group in probe_sets.items():
                total_probes_for_this_level_group = len(probes_in_group)

                relevant_fasta_seq_ids = species_dict.get(phylo_level_name, set())
                if not relevant_fasta_seq_ids:
                    logging.debug(f"No FASTA sequences were associated with PhyloLevel '{phylo_level_name}'. Coverage will be NA.")

                coverage_values_for_level = [
                    cov_percent for res_level, res_target_id, cov_percent in coverage_results
                    if res_level == phylo_level_name and res_target_id in relevant_fasta_seq_ids
                ]

                if coverage_values_for_level:
                    series = pd.Series(coverage_values_for_level)
                    q1 = series.quantile(0.25)
                    q2 = series.quantile(0.50)
                    q3 = series.quantile(0.75)
                else:
                    q1 = q2 = q3 = "NA"
                    if relevant_fasta_seq_ids:
                         logging.debug(f"PhyloLevel '{phylo_level_name}' has associated FASTA sequences ({len(relevant_fasta_seq_ids)}), but no coverage results from BLAST. Check BLAST file or alignment parameters.")
                output.write(f"{phylo_level_name}\t{total_probes_for_this_level_group}\t{q1}\t{q2}\t{q3}\n")

    save_probes_by_phylo_level(probe_sets, args.blast_file, args.output_probe_prefix, args.probe_fasta_file)
    logging.info("Processing completed.")

if __name__ == "__main__":
    main()

