import argparse
import logging
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import os
import pandas as pd
from typing import List, Dict, Set, Tuple, Any

# Logging setup
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - [%(filename)s:%(lineno)d] - %(message)s'
)

# Define a type for a lineage item: (level_name, level_tax_id, level_rank)
LineageItem = Tuple[str, str, str]
# Define a type for a full lineage: List of LineageItem
Lineage = List[LineageItem]

RANK_PRIORITY_LIST = [
    "kingdom", "subkingdom", "phylum", "subphylum",
    "class", "subclass", "order", "suborder",
    "family", "subfamily", "genus", "species" # Added species for completeness if needed
]
RANK_PRIORITY_SET = set(RANK_PRIORITY_LIST)
RANK_ORDER_MAP = {rank: i for i, rank in enumerate(RANK_PRIORITY_LIST)}


def process_species_name(species_name: str) -> str:
    """Normalize species name to retain only genus and species level.

    Args:
        species_name (str): Raw species name potentially containing underscores, hyphens, or extra parts.

    Returns:
        str: Normalized species name with only the first two words (genus and species).
    """
    return " ".join(species_name.replace('_', ' ').replace('-', ' ').split()[:2])

def parse_taxonomy_data(taxonomy_file: str, taxonomy_nodes_file: str) -> Dict[str, Dict[str, str]]:
    """Parse taxonomy and node files to build a comprehensive map from tax_id to its details.

    Reads the taxonomy names file and the NCBI-style nodes file to construct
    a unified dictionary of tax IDs with their name, rank, and parent ID.

    Args:
        taxonomy_file (str): Path to taxonomy names file (TSV: tax_id, rank, name).
        taxonomy_nodes_file (str): Path to taxonomy nodes file (TSV: tax_id | parent_id | rank | ...).

    Returns:
        Dict[str, Dict[str, str]]: Mapping of tax_id to {"name": str, "rank": str, "parent_id": str}.
    """
    logging.info("Parsing taxonomy information and parent-child relationships...")
    tax_id_to_details: Dict[str, Dict[str, Any]] = defaultdict(dict)

    # Read taxonomy_file: tax_id \t rank \t name
    try:
        with open(taxonomy_file, 'r') as f:
            for i, line in enumerate(f):
                parts = line.strip().split('\t')
                if len(parts) < 3:
                    logging.warning(f"Malformed line {i+1} in {taxonomy_file}: {line.strip()}")
                    continue
                tax_id, rank, name = parts[0], parts[1], parts[2]
                tax_id_to_details[tax_id].update({"name": name, "rank": rank.lower()})
    except FileNotFoundError:
        logging.error(f"File not found: {taxonomy_file}")
        raise
    logging.info(f"Loaded {len(tax_id_to_details)} entries from {taxonomy_file}")

    # Read taxonomy_nodes_file: tax_id \t|\t parent_id \t|\t rank \t|\t ...
    # This file establishes hierarchy and can supplement rank information.
    found_in_nodes = 0
    updated_from_nodes = 0
    try:
        with open(taxonomy_nodes_file, 'r') as f:
            for i, line in enumerate(f):
                parts = line.strip().split("\t|\t")
                if len(parts) < 3: # tax_id, parent_id, rank are essential
                    logging.warning(f"Malformed line {i+1} in {taxonomy_nodes_file}: {line.strip()}")
                    continue
                tax_id, parent_id, rank_from_nodes = parts[0].strip(), parts[1].strip(), parts[2].strip().lower()

                found_in_nodes += 1
                # Ensure entry exists
                entry = tax_id_to_details[tax_id] # defaultdict creates if not exists
                entry["parent_id"] = parent_id
                # Prefer rank from nodes.dmp if taxonomy_file didn't provide one or it's generic
                if not entry.get("rank") or entry.get("rank") == "no rank":
                    entry["rank"] = rank_from_nodes
                    updated_from_nodes +=1
                # If name is missing, it's an issue, but store tax_id as name placeholder
                if "name" not in entry:
                    entry["name"] = f"Unnamed_taxid_{tax_id}"

    except FileNotFoundError:
        logging.error(f"File not found: {taxonomy_nodes_file}")
        raise
    logging.info(f"Processed {found_in_nodes} entries from {taxonomy_nodes_file}, updated/added rank for {updated_from_nodes} entries.")
    logging.info(f"Total unique TaxIDs in combined taxonomy: {len(tax_id_to_details)}")
    return tax_id_to_details

def parse_fasta_species(fasta_file: str) -> Set[str]:
    """Extract and process species names from a FASTA file.

    Args:
        fasta_file (str): Path to the FASTA file.

    Returns:
        Set[str]: Set of unique processed species names extracted from FASTA headers.
    """
    logging.info(f"Parsing species names from FASTA file {fasta_file}...")
    species_names: Set[str] = set()
    try:
        with open(fasta_file, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    # Example header: >species.name.identifier other_stuff
                    header_part = line.split('.', 1)[0].replace('>', '').strip()
                    if not header_part: # Handle cases like ">" or "> "
                         header_part = line.split(' ',1)[0].replace('>', '').strip()

                    processed_name = process_species_name(header_part)
                    if processed_name:
                        species_names.add(processed_name)
    except FileNotFoundError:
        logging.error(f"File not found: {fasta_file}")
        raise
    logging.info(f"Extracted {len(species_names)} unique processed species names from FASTA.")
    return species_names

def parse_alignment_file(alignment_file: str) -> Dict[str, Set[str]]:
    """Parse alignment results to map probe_id to a set of hit species names.

    Args:
        alignment_file (str): Path to the alignment results file (TSV).

    Returns:
        Dict[str, Set[str]]: Mapping of probe IDs to sets of hit species names.
    """
    logging.info(f"Parsing alignment file {alignment_file}...")
    probe_to_hit_species: Dict[str, Set[str]] = defaultdict(set)
    try:
        with open(alignment_file, 'r') as f:
            for i, line in enumerate(f):
                parts = line.strip().split('\t')
                if len(parts) < 2:
                    logging.warning(f"Malformed line {i+1} in alignment file: {line.strip()}")
                    continue
                probe_id, species_name_raw = parts[0], parts[1]
                processed_name = process_species_name(species_name_raw)
                if processed_name:
                    probe_to_hit_species[probe_id].add(processed_name)
    except FileNotFoundError:
        logging.error(f"File not found: {alignment_file}")
        raise
    logging.info(f"Parsed {len(probe_to_hit_species)} probes from alignment file.")
    return probe_to_hit_species

def build_lineages_and_totals(
    all_species_names_from_fasta: Set[str],
    tax_id_to_details: Dict[str, Dict[str, str]]
) -> Tuple[Dict[str, str], Dict[str, Lineage], Dict[LineageItem, Set[str]]]:
    """Build species name to tax_id map, cache lineages, and precompute total species per phylogenetic level.

    Args:
        all_species_names_from_fasta (Set[str]): Set of all processed species names from the FASTA file.
        tax_id_to_details (Dict[str, Dict[str, str]]): Mapping of tax_id to name, rank, and parent_id.

    Returns:
        Tuple[Dict[str, str], Dict[str, Lineage], Dict[LineageItem, Set[str]]]: A tuple containing:
            - species_name_to_tax_id: Mapping of species names to their tax IDs.
            - lineage_cache: Cached lineages for each tax ID.
            - phylo_level_totals: Mapping of phylogenetic level items to sets of species.
    """
    logging.info("Building species-to-tax_id map, lineages, and phylogenetic totals...")

    species_name_to_tax_id: Dict[str, str] = {}
    for tax_id, details in tax_id_to_details.items():
        if details.get("rank", "").lower() == "species" and "name" in details:
            species_name_to_tax_id[details["name"]] = tax_id
            processed_canonical_name = process_species_name(details["name"])
            if processed_canonical_name not in species_name_to_tax_id:
                 species_name_to_tax_id[processed_canonical_name] = tax_id


    logging.info(f"Built map for {len(species_name_to_tax_id)} species names to TaxIDs.")

    lineage_cache: Dict[str, Lineage] = {}
    def get_lineage(current_tax_id: str) -> Lineage:
        if not current_tax_id or current_tax_id == "0": # "0" sometimes used for no parent
            return []
        if current_tax_id in lineage_cache:
            return lineage_cache[current_tax_id]

        details = tax_id_to_details.get(current_tax_id)
        if not details or "name" not in details: # Should not happen if data is consistent
            logging.warning(f"TaxID {current_tax_id} missing in tax_id_to_details or lacks a name. Stopping lineage.")
            lineage_cache[current_tax_id] = []
            return []

        name = details["name"]
        rank = details.get("rank", "unknown rank")
        parent_id = details.get("parent_id")

        current_lineage: Lineage = []
        # Recurse for parent, avoid self-recursion if parent_id is same as current_tax_id (e.g. root "1")
        if parent_id and parent_id != current_tax_id and parent_id != "0" and parent_id != "1": # "1" is NCBI root
            current_lineage.extend(get_lineage(parent_id))

        # Add current node to lineage
        # Only add if rank is in our priority list or it's species (always include species itself)
        if rank in RANK_PRIORITY_SET or rank == "species":
             current_lineage.append((name, current_tax_id, rank))

        lineage_cache[current_tax_id] = current_lineage

        return current_lineage

    phylo_level_totals: Dict[LineageItem, Set[str]] = defaultdict(set)
    processed_species_count = 0
    for species_name in all_species_names_from_fasta:
        species_tax_id = species_name_to_tax_id.get(species_name)
        if not species_tax_id:
            # Attempt with a more broadly processed name from FASTA if direct match failed
            # This part can be tricky if FASTA names differ significantly from taxonomy names
            logging.debug(f"Species '{species_name}' from FASTA not directly mapped to TaxID. Check name consistency.")
            continue

        processed_species_count +=1
        lineage = get_lineage(species_tax_id)
        for phylo_level_item in lineage: # (name, tax_id, rank)
            # We only care about totals for ranks in our priority list (excluding species itself for 'total')
            if phylo_level_item[2] in RANK_PRIORITY_SET: # rank is phylo_level_item[2]
                 phylo_level_totals[phylo_level_item].add(species_name)

    logging.info(f"Built lineages for {len(lineage_cache)} TaxIDs based on {processed_species_count} FASTA species.")
    logging.info(f"Precomputed totals for {len(phylo_level_totals)} phylogenetic levels.")
    return species_name_to_tax_id, lineage_cache, phylo_level_totals


def analyze_probe_task(
    probe_id: str,
    hit_species_names_for_probe: Set[str],
    lineage_cache: Dict[str, Lineage],
    species_name_to_tax_id: Dict[str, str],
    phylo_level_totals: Dict[LineageItem, Set[str]],
    self_species_name: str,
    ratio_threshold: float
) -> List[Tuple[str, str, str, str, int, int, float, str]]:
    """Analyze a single probe, calculating hit counts and ratios for relevant phylogenetic levels.

    For each phylogenetic level, computes the hit count and hit ratio based on the
    probe's hit species and the total species at that level.

    Args:
        probe_id (str): The probe identifier.
        hit_species_names_for_probe (Set[str]): Set of hit species names for this probe.
        lineage_cache (Dict[str, Lineage]): Cached lineages for tax IDs.
        species_name_to_tax_id (Dict[str, str]): Mapping of species names to tax IDs.
        phylo_level_totals (Dict[LineageItem, Set[str]]): Precomputed totals per phylogenetic level.
        self_species_name (str): The name of the species this probe originates from.
        ratio_threshold (float): Minimum hit ratio threshold for inclusion.

    Returns:
        List[Tuple[str, str, str, str, int, int, float, str]]: List of tuples containing
            (ProbeID, PhyloLevel, TaxID, Rank, Total, Hit, HitRatio, SelfSpeciesName).
    """
    probe_results: List[Tuple[str, str, str, str, int, int, float, str]] = []

    # Determine which phylogenetic levels this probe's hits belong to
    probe_hits_at_phylo_level: Dict[LineageItem, Set[str]] = defaultdict(set)
    for hit_species_name in hit_species_names_for_probe:
        hit_species_tax_id = species_name_to_tax_id.get(hit_species_name)
        if not hit_species_tax_id:
            logging.debug(f"Probe {probe_id}: Hit species '{hit_species_name}' not found in taxonomy map.")
            continue

        # It's possible a hit species tax_id might not have a pre-cached lineage if it wasn't in FASTA
        # This implies get_lineage might need to be callable here, or ensure all hit species are processed by build_lineages_and_totals
        # For simplicity, assuming lineage_cache is comprehensive for all hit species.
        # If not, a call to get_lineage(hit_species_tax_id) using tax_id_to_details (if passed) would be needed.
        lineage = lineage_cache.get(hit_species_tax_id, [])

        for phylo_level_item in lineage: # (name, tax_id, rank)
            # Only consider levels that are in our precomputed totals (i.e., valid ranks)
            if phylo_level_item[2] in RANK_PRIORITY_SET:
                probe_hits_at_phylo_level[phylo_level_item].add(hit_species_name)

    # Calculate stats for each relevant phylogenetic level
    for phylo_level_item, total_species_set in phylo_level_totals.items():
        level_name, level_tax_id, level_rank = phylo_level_item
        total_count = len(total_species_set)
        if total_count == 0: # Should not happen if phylo_level_totals is built correctly
            continue

        hit_species_set = probe_hits_at_phylo_level.get(phylo_level_item, set())
        hit_count = len(hit_species_set)
        hit_ratio = hit_count / total_count if total_count > 0 else 0.0

        if hit_ratio >= ratio_threshold:
            probe_results.append((
                probe_id, level_name, level_tax_id, level_rank,
                total_count, hit_count, hit_ratio, self_species_name
            ))
    return probe_results

def process_all_probes(
    alignment_file: str, fasta_file: str, taxonomy_file: str,
    taxonomy_nodes_file: str, ratio_threshold: float, num_workers: int
) -> List[Tuple[str, str, str, str, int, int, float, str]]:
    """Main processing orchestrator.

    Loads taxonomy data, parses FASTA species names and alignment results,
    builds lineages, and analyzes all probes in parallel.

    Args:
        alignment_file (str): Path to alignment results file (TSV).
        fasta_file (str): Path to FASTA file with sequences of all species.
        taxonomy_file (str): Path to taxonomy names file (TSV).
        taxonomy_nodes_file (str): Path to taxonomy nodes file.
        ratio_threshold (float): Minimum hit ratio threshold.
        num_workers (int): Number of worker processes.

    Returns:
        List[Tuple[str, str, str, str, int, int, float, str]]: List of analysis result tuples.
    """

    tax_id_to_details = parse_taxonomy_data(taxonomy_file, taxonomy_nodes_file)
    all_species_names_from_fasta = parse_fasta_species(fasta_file)
    probe_to_hit_species = parse_alignment_file(alignment_file)

    species_name_to_tax_id, lineage_cache, phylo_level_totals = build_lineages_and_totals(
        all_species_names_from_fasta, tax_id_to_details
    )

    all_analysis_results: List[Tuple[str, str, str, str, int, int, float, str]] = []
    futures = []

    logging.info(f"Starting probe analysis with {num_workers} workers...")
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        for probe_id, hit_species_names in probe_to_hit_species.items():
            self_species_name = process_species_name(probe_id.split('.')[0]) # Assuming probe ID format like 'Genus_species.probe123'

            future = executor.submit(
                analyze_probe_task,
                probe_id,
                hit_species_names,
                lineage_cache,
                species_name_to_tax_id,
                phylo_level_totals,
                self_species_name,
                ratio_threshold
            )
            futures.append(future)

        for i, future in enumerate(as_completed(futures)):
            try:
                probe_results = future.result()
                all_analysis_results.extend(probe_results)
                if (i + 1) % 1000 == 0: # Log progress
                    logging.info(f"Completed analysis for {i+1}/{len(futures)} probes.")
            except Exception as e:
                logging.error(f"Error processing a probe task: {e}", exc_info=True)

    logging.info(f"Finished analysis. Total results generated: {len(all_analysis_results)}")
    return all_analysis_results

def filter_and_write_highest_taxonomic_level(
    results: List[Tuple[str, str, str, str, int, int, float, str]],
    output_file: str,
    total_threshold: int
):
    """Filter results to get the highest taxonomic level per probe and write to file.

    Groups results by probe, selects the highest-priority taxonomic rank per probe,
    and writes the result to a tab-separated file.

    Args:
        results (List[Tuple[str, str, str, str, int, int, float, str]]): List of analysis result tuples.
        output_file (str): Path to the output file.
        total_threshold (int): Minimum total species at a phylo_level to be considered.
    """
    if not results:
        logging.warning("No results to process for highest taxonomic level filtering.")
        if output_file: # Create empty file with header if output_file is specified
             with open(output_file, 'w') as out_f:
                out_f.write("ProbeID\tPhyloLevel\tTaxID\tRank\tTotal\tHit\tHitRatio\tSelfSpeciesName\n")
        return

    df = pd.DataFrame(results, columns=[
        "ProbeID", "PhyloLevel", "TaxID", "Rank",
        "Total", "Hit", "HitRatio", "SelfSpeciesName"
    ])
    logging.info(f"Initial records for filtering: {len(df)}")

    df = df[df['Total'] >= total_threshold]
    logging.info(f"Records after Total >= {total_threshold} filter: {len(df)}")

    if df.empty:
        logging.warning(f"No records remain after 'Total >= {total_threshold}' filter.")
        with open(output_file, 'w') as out_f: # Create empty file with header
            out_f.write("ProbeID\tPhyloLevel\tTaxID\tRank\tTotal\tHit\tHitRatio\tSelfSpeciesName\n")
        return

    # Use the predefined rank order for sorting (lower index means higher rank)
    df['RankPriority'] = df['Rank'].str.lower().map(RANK_ORDER_MAP).fillna(len(RANK_ORDER_MAP))

    # Get the row with the minimum RankPriority for each ProbeID
    # idxmin returns the index of the first occurrence of the minimum value
    highest_rank_df = df.loc[df.groupby('ProbeID')['RankPriority'].idxmin()]
    highest_rank_df = highest_rank_df.drop(columns=['RankPriority'])

    logging.info(f"Writing {len(highest_rank_df)} highest taxonomic level results to {output_file}")
    highest_rank_df.to_csv(output_file, sep='\t', index=False, header=True)


def main():
    """Parse command-line arguments and run the probe analysis workflow."""
    parser = argparse.ArgumentParser(description="Analyze probe alignment results and generate reports.")
    parser.add_argument("alignment_file", type=str, help="Path to alignment results file (TSV: probe_id, species_name_raw).")
    parser.add_argument("fasta_file", type=str, help="Path to FASTA file containing sequences of all species.")
    parser.add_argument("taxonomy_file", type=str, help="Path to taxonomy file (TSV: tax_id, rank, name).")
    parser.add_argument("taxonomy_nodes_file", type=str, help="Path to NCBI-like taxonomy nodes file (TSV: tax_id | parent_id | rank | ...).")
    parser.add_argument("output_file", type=str, help="Path to output file for all results passing ratio_threshold.")
    parser.add_argument("filter_output_file", type=str, help="Path to output file for highest taxonomic level results after filtering.")
    parser.add_argument("--total_threshold", type=int, default=0, help="Minimum total species at a phylo_level to be considered for the filtered output (default: 0).")
    parser.add_argument("--ratio_threshold", type=float, default=0.7, help="Minimum HitRatio for a result to be included (default: 0.7).")
    parser.add_argument("--num_workers", type=int, default=os.cpu_count(), help="Number of worker processes (default: all CPU cores).")

    args = parser.parse_args()
    logging.info("Program started with arguments: %s", args)

    all_results = process_all_probes(
        args.alignment_file, args.fasta_file, args.taxonomy_file,
        args.taxonomy_nodes_file, args.ratio_threshold, args.num_workers
    )

    # Write all results that passed the ratio threshold
    logging.info(f"Writing {len(all_results)} raw results (passing ratio threshold) to {args.output_file}...")
    header = "ProbeID\tPhyloLevel\tTaxID\tRank\tTotal\tHit\tHitRatio\tSelfSpeciesName\n"
    with open(args.output_file, 'w') as output_f:
        output_f.write(header)
        for result_tuple in all_results:
            output_f.write("\t".join(map(str, result_tuple)) + "\n")
    logging.info(f"Raw results written to {args.output_file}.")

    # Filter for highest taxonomic level and write
    filter_and_write_highest_taxonomic_level(all_results, args.filter_output_file, args.total_threshold)
    logging.info(f"Filtered highest taxonomic level results written to {args.filter_output_file}.")
    logging.info("Program finished successfully.")

if __name__ == "__main__":
    main()


