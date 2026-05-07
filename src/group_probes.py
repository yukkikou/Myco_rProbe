import argparse
import logging
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
import os
import pandas as pd

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


def clean_probe_id(probe_id):
    """Clean the probe id by removing the sliding part if exists.

    Args:
        probe_id (str): Raw probe identifier potentially containing a '_sliding:' segment.

    Returns:
        str: Cleaned probe ID without the sliding coordinate suffix.
    """
    if "_sliding:" in probe_id:
        return probe_id.split("_sliding:")[0]
    return probe_id


def parse_phylogenetic_info_with_filter(alignment_file, fasta_file):
    """Parse phylogenetic info from the FASTA file while filtering by alignment file query IDs.

    Args:
        alignment_file (str): File containing probe alignment results.
        fasta_file (str): Original FASTA file containing phylogenetic information.

    Returns:
        tuple: A tuple containing:
            - species_to_phylo: Mapping of species ID to phylogenetic levels.
            - species_name_map: Mapping of species ID to species name (last phylogenetic level).
    """
    logging.info(f"Parsing phylogenetic info from {fasta_file} with filtering based on {alignment_file}...")

    # Step 1: Extract valid IDs from alignment file
    valid_ids = set()
    with open(alignment_file, 'r') as f:
        for line in f:
            query_id = line.strip().split()[0]
            clean_id = clean_probe_id(query_id)
            valid_ids.add(clean_id)

    logging.info(f"Collected {len(valid_ids)} valid IDs from alignment file.")

    # Step 2: Parse the FASTA file with filtering
    species_to_phylo = {}
    species_name_map = {}
    with open(fasta_file, 'r') as f:
        current_species_id = None
        current_phylo_levels = []
        for line in f:
            if line.startswith('>'):  # Header line
                parts = line.strip().split(' ', 1)
                species_id = parts[0][1:]  # Skip the '>'
                if species_id in valid_ids:  # Only process if ID is valid
                    current_species_id = species_id
                    current_phylo_levels = parts[1].split(';')
                    species_to_phylo[species_id] = current_phylo_levels
                    species_name_map[species_id] = current_phylo_levels[-1] if current_phylo_levels else "Unknown"
            else:
                continue  # Skip sequence lines (not required for filtering)

    logging.info(f"Loaded phylogenetic information for {len(species_to_phylo)} valid species.")
    return species_to_phylo, species_name_map


def parse_alignment_results(alignment_file):
    """Parse alignment results from the alignment file.

    Args:
        alignment_file (str): Path to the alignment results file.

    Returns:
        tuple: A tuple containing:
            - probe_to_species_RAW: Mapping of raw probe IDs to sets of hit species IDs.
            - probe_to_self_species: Mapping of raw probe IDs to their self-species ID.
    """
    logging.info(f"Parsing alignment results from {alignment_file}...")
    probe_to_species_RAW = defaultdict(set)
    probe_to_self_species = {}

    with open(alignment_file, 'r') as f:
        for line in f:
            columns = line.strip().split()
            raw_probe_id = columns[0]
            species_id = columns[1]
            probe_to_species_RAW[raw_probe_id].add(species_id)
            # remove probe_to_self_species
            if raw_probe_id not in probe_to_self_species:
                probe_to_self_species[raw_probe_id] = species_id

    logging.info(f"Parsed alignment results for {len(probe_to_species_RAW)} probes.")
    return probe_to_species_RAW, probe_to_self_species


def parse_taxonomy_info(taxonomy_file):
    """Parse taxonomy information from the taxonomy file.

    Args:
        taxonomy_file (str): Path to the taxonomy information file.

    Returns:
        dict: Dictionary mapping taxonomic level names to their tax ID and rank.
    """
    logging.info(f"Parsing taxonomy information from {taxonomy_file}...")
    tax_info = {}
    with open(taxonomy_file, 'r') as f:
        next(f)
        for line in f:
            tax_id, rank, name = line.strip().split('\t')
            tax_info[name] = {'taxid': tax_id, 'rank': rank}
    logging.info(f"Loaded taxonomy information for {len(tax_info)} phylogenetic levels.")
    return tax_info


def precompute_phylo_total(species_to_phylo):
    """Precompute phylo_total information mapping phylogenetic levels to their species sets.

    Args:
        species_to_phylo (dict): Mapping of species IDs to lists of phylogenetic levels.

    Returns:
        defaultdict: Mapping of phylogenetic level names to sets of species IDs.
    """
    phylo_total = defaultdict(set)
    for species, levels in species_to_phylo.items():
        for level in levels[:-1]:
            phylo_total[level].add(species)
    return phylo_total


def analyze_probe(raw_probe_id, species_set, probe_to_self_species, phylo_total, species_to_phylo, species_name_map, taxonomy_info, ratio_threshold):
    """Analyze a single probe for its taxonomic implications.

    Args:
        raw_probe_id (str): Raw probe identifier.
        species_set (set): Set of species IDs hit by this probe.
        probe_to_self_species (dict): Mapping of raw probe IDs to self-species IDs.
        phylo_total (defaultdict): Precomputed phylogenetic level to species mapping.
        species_to_phylo (dict): Mapping of species IDs to phylogenetic level lists.
        species_name_map (dict): Mapping of species IDs to species names.
        taxonomy_info (dict): Dictionary of taxonomic level metadata.
        ratio_threshold (float): Hit ratio threshold for inclusion.

    Returns:
        list: List of result tuples (ProbeID, PhyloLevel, TaxID, Rank, Total, Hit, HitRatio, SelfSpeciesName).
    """
    clean_id = clean_probe_id(raw_probe_id)
    logging.debug(f"Analyzing probe {clean_id}...")

    phylo_hit = defaultdict(set)
    species_names = set()
    self_species_name = species_name_map.get(clean_id, ".")

    for species in species_set:
        species_names.add(species_name_map.get(species, "Unknown"))
        levels = species_to_phylo.get(species, [])
        for level in levels[:-1]:
            phylo_hit[level].add(species)

    results = []
    for level in phylo_total:
        taxid = taxonomy_info.get(level, {}).get("taxid", "Unknown")
        rank = taxonomy_info.get(level, {}).get("rank", "Unknown")

        total = len(phylo_total[level])
        hit = len(phylo_hit.get(level, set()))
        hit_ratio = hit / total if total > 0 else 0

        if hit_ratio >= ratio_threshold:
            results.append((
                raw_probe_id,
                level,
                taxid,
                rank,
                total,
                hit,
                hit_ratio,
                self_species_name
            ))

    return results


def process_results(alignment_file, fasta_file, taxonomy_file, ratio_threshold, num_workers):
    """Process alignment results and apply taxonomic analysis.

    Args:
        alignment_file (str): Path to the alignment results file.
        fasta_file (str): Path to the FASTA file with phylogenetic information.
        taxonomy_file (str): Path to the taxonomy information file.
        ratio_threshold (float): Hit ratio threshold for inclusion.
        num_workers (int): Number of worker processes for parallel processing.

    Returns:
        list: List of result tuples from probe analysis.
    """
    logging.info("Starting analysis process...")
    species_to_phylo, species_name_map = parse_phylogenetic_info_with_filter(alignment_file, fasta_file)
    probe_to_species_RAW, probe_to_self_species = parse_alignment_results(alignment_file)
    taxonomy_info = parse_taxonomy_info(taxonomy_file)

    phylo_total = precompute_phylo_total(species_to_phylo)

    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = {
            executor.submit(analyze_probe, raw_probe_id, species_set, probe_to_self_species, phylo_total, species_to_phylo, species_name_map, taxonomy_info, ratio_threshold)
            for raw_probe_id, species_set in probe_to_species_RAW.items()
        }
        all_results = []
        for future in futures:
            all_results.extend(future.result())

    logging.info("Analysis process completed.")
    return all_results


def get_highest_taxonomic_level(results, output_file, total_threshold):
    """Identify and write highest taxonomic levels of the results.

    Filters results by total threshold, ranks by taxonomic priority, and writes
    the highest rank per probe to a tab-separated file.

    Args:
        results (list): List of result tuples from analyze_probe.
        output_file (str): Path to the output file.
        total_threshold (int): Minimum total species required per tax level.
    """
    df = pd.DataFrame(results, columns=["ProbeID", "PhyloLevel", "TaxID", "Rank", "Total", "Hit", "HitRatio", "SelfSpeciesName"])

    df = df[df['Total'] >= total_threshold]
    if df.empty:
        logging.warning("No records remaining after total threshold filter.")
        return

    # Order by priority of taxonomic rank
    rank_priority = ["kingdom","subkingdom", "phylum","subphylum", "class","subclass", "order", "suborder","family","subfamily", "genus","subgenus", "species"]
    df['RankPriority'] = df['Rank'].map(lambda x: rank_priority.index(x.lower()) if x.lower() in rank_priority else len(rank_priority))

    highest_rank_df = df.loc[df.groupby('ProbeID')['RankPriority'].idxmin()].drop(columns=['RankPriority'])
    highest_rank_df.to_csv(output_file, sep='\t', index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Analyze probe alignment results and generate reports.')
    parser.add_argument('alignment_file', type=str, help='Path to the file containing probe alignment results.')
    parser.add_argument('fasta_file', type=str, help='Path to the FASTA file containing phylogenetic information.')
    parser.add_argument('taxonomy_file', type=str, help='Path to the taxonomy information file.')
    parser.add_argument('output_file', type=str, help='Path to the output file to write the results.')
    parser.add_argument('filter_output_file', type=str, help='Path to the filtered output file to write the highest taxonomic level results.')
    parser.add_argument('--total_threshold', type=int, default=0, help='Minimum total species required per tax level (default: 0).')
    parser.add_argument('--ratio_threshold', type=float, default=0.7, help='Hit ratio threshold.')
    parser.add_argument('--num_workers', type=int, default=os.cpu_count(), help='Number of worker processes for parallel processing.')

    args = parser.parse_args()

    logging.info("Program started with arguments: %s", args)

    # Step 1: Process results directly without intermediate filtered FASTA file
    results = process_results(args.alignment_file, args.fasta_file, args.taxonomy_file, args.ratio_threshold, args.num_workers)

    # Step 2: Write full results to file
    logging.info(f"Writing results to {args.output_file}...")
    with open(args.output_file, 'w') as output:
        output.write("ProbeID\tPhyloLevel\tTaxID\tRank\tTotal\tHit\tHitRatio\tSelfSpeciesName\n")
        for result in results:
            output.write("\t".join(map(str, result)) + "\n")
    logging.info("Results writing completed.")

    # Step 3: Write filtered highest taxonomic level results
    get_highest_taxonomic_level(results, args.filter_output_file, args.total_threshold)
    logging.info(f"Filtered highest taxonomic level results written to {args.filter_output_file}.")

