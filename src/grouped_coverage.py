import argparse
import logging
import pandas as pd
from collections import defaultdict
from Bio import SeqIO

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


def filter_hits(input_file, hit_threshold):
    """Filter rows with HitRatio greater than hit_threshold.

    Args:
        input_file (str): Path to the input TSV file containing HitRatio column.
        hit_threshold (float): Threshold value for HitRatio filtering.

    Returns:
        pd.DataFrame: Filtered DataFrame containing rows with HitRatio above the threshold.
    """
    logging.info("Filtering hits...")
    df = pd.read_csv(input_file, sep='\t')
    filtered_df = df[df['HitRatio'] > hit_threshold]
    logging.info(f"Filtered to {len(filtered_df)} rows.")
    return filtered_df


def classify_probes(filtered_df, classification_level):
    """Classify probes into groups based on the specified classification level.

    Args:
        filtered_df (pd.DataFrame): DataFrame containing probe data with 'Rank' and 'PhyloLevel' columns.
        classification_level (str): The taxonomic rank to classify by (e.g., 'subkingdom').

    Returns:
        defaultdict: Dictionary mapping PhyloLevel names to lists of ProbeIDs.
    """
    logging.info(f"Classifying probes by {classification_level}...")
    probe_sets = defaultdict(list)

    for _, row in filtered_df.iterrows():
        if row['Rank'] == classification_level:
            probe_sets[row['PhyloLevel']].append(row['ProbeID'])

    logging.info(f"Classification results: {len(probe_sets)} groups.")
    return probe_sets


def determine_species(fasta_file, probe_sets):
    """Determine species based on probe set phylogenetic levels and extract sequence lengths.

    Args:
        fasta_file (str): Path to the background FASTA file.
        probe_sets (defaultdict): Dictionary mapping PhyloLevel names to lists of ProbeIDs.

    Returns:
        tuple: A tuple containing:
            - species_dict: Mapping of PhyloLevel names to sets of sequence IDs.
            - sequence_lengths: Mapping of sequence IDs to their lengths.
            - full_sequences: Mapping of sequence IDs to their full sequence strings.
    """
    logging.info("Determining species based on probe sets and extracting sequence lengths...")

    species_dict = defaultdict(set)
    sequence_lengths = {}
    full_sequences = defaultdict(str)  # New dictionary to store complete sequences
    phylogy_levels = set(probe_sets.keys())
    logging.info(f"Identified phylogy levels: {phylogy_levels}")

    with open(fasta_file, 'r') as f:
        current_seq_id = None
        for line in f:
            if line.startswith('>'):
                current_seq_id = line[1:].split(' ', 1)[0].strip()
                header_parts = line[1:].strip().split(' ', 1)
                if len(header_parts) > 1:
                    phylo_levels = header_parts[1].split(';')
                    for level in phylogy_levels:
                        if level in phylo_levels:
                            species_dict[level].add(current_seq_id)
                            sequence_lengths[current_seq_id] = 0
                            full_sequences[current_seq_id] = ""  # Initialize sequence storage
                            break
            elif current_seq_id and current_seq_id in sequence_lengths:
                sequence_lengths[current_seq_id] += len(line.strip())
                full_sequences[current_seq_id] += line.strip()  # Save full sequence
    logging.info(f"Extracted sequence lengths for {len(sequence_lengths)} sequences.")
    return species_dict, sequence_lengths, full_sequences


def calculate_coverage(probe_sets, blast_file, sequence_lengths):
    """Calculate the coverage of probe sets within the classification level.

    Args:
        probe_sets (defaultdict): Mapping of PhyloLevel names to lists of ProbeIDs.
        blast_file (str): Path to the BLAST results file.
        sequence_lengths (dict): Mapping of target sequence IDs to their lengths.

    Returns:
        tuple: A tuple containing:
            - coverage_results: List of (PhyloLevel, target_id, coverage_percent) tuples.
            - probe_coverage_counts: Mapping of target IDs to counts of covering probes.
    """
    logging.info("Calculating coverage...")
    coverage_results = []
    covered_regions = defaultdict(lambda: defaultdict(set))
    target_probes = defaultdict(set)

    with open(blast_file, 'r') as f:
        for line in f:
            columns = line.strip().split()
            probe_id = columns[0]
            target_id = columns[1]
            for phylo_level, probes in probe_sets.items():
                if probe_id in probes:
                    target_probes[target_id].add(probe_id)

                    start, end = int(columns[8]), int(columns[9])
                    covered_range = range(min(start, end), max(start, end) + 1)
                    for pos in covered_range:
                        covered_regions[phylo_level][target_id].add(pos)
                    break

    probe_coverage_counts = {k: len(v) for k, v in target_probes.items()}

    for phylo_level, targets in covered_regions.items():
        for target_id in targets:
            if target_id in sequence_lengths and sequence_lengths[target_id] > 0:
                total_coverage = len(targets[target_id])
                coverage_percent = total_coverage / sequence_lengths[target_id] * 100
                coverage_results.append((phylo_level, target_id, coverage_percent))
    logging.info("Coverage calculation completed.")
    return coverage_results, probe_coverage_counts


def save_probes_by_phylo_level(probe_sets, blast_file, output_probe_prefix, probe_fasta):
    """Save probe sequences (FASTA) belonging to specified PhyloLevels into separate files.

    Args:
        probe_sets (defaultdict): Mapping of PhyloLevel names to lists of ProbeIDs.
        blast_file (str): Path to the BLAST results file.
        output_probe_prefix (str): Prefix for output FASTA filenames.
        probe_fasta (str): Path to the probe FASTA file containing all probe sequences.
    """
    logging.info("Writing probe sequences to separate FASTA files by PhyloLevel...")
    probe_hits_by_phylo_level = defaultdict(set)

    # Load probes from FASTA file
    probe_sequences = {
        record.id: record for record in SeqIO.parse(probe_fasta, "fasta")
    }

    # Read BLAST file to map probes to targets
    with open(blast_file, 'r') as f:
        for line in f:
            columns = line.strip().split()
            probe_id = columns[0]
            target_id = columns[1]
            for phylo_level, probes in probe_sets.items():
                if probe_id in probes:  # Check if the probe belongs to this phylo level
                    probe_hits_by_phylo_level[phylo_level].add(probe_id)
                    break

    # Write probe sequences to separate FASTA files by PhyloLevel
    for level, probes in probe_hits_by_phylo_level.items():
        filename = f"{output_probe_prefix}_{level}_probes.fasta"
        with open(filename, 'w') as probe_out:
            for probe_id in sorted(probes):
                if probe_id in probe_sequences:
                    SeqIO.write(probe_sequences[probe_id], probe_out, "fasta")
        logging.info(f"Probe FASTA file for {level} written to {filename}.")


def main():
    """Parse command-line arguments and run the grouped coverage analysis workflow."""
    parser = argparse.ArgumentParser(description='Analyze probe data and calculate coverage.')
    parser.add_argument('filter_output_file', type=str, help='Filter output file path')
    parser.add_argument('hit_threshold', type=float, help='Hit threshold (can be integer or float)')
    parser.add_argument('classification_level', type=str, help='Classification level (e.g., subkingdom)')
    parser.add_argument('fasta_file', type=str, help='Background FASTA file path')
    parser.add_argument('blast_file', type=str, help='BLAST results file path')
    parser.add_argument('output_file', type=str, help='Output TSV file path')
    parser.add_argument('probe_fasta_file', type=str, help='Probe FASTA file path')  # New parameter for probe FASTA file
    parser.add_argument('output_probe_prefix', type=str, help='Output probe file prefix')  # Modified parameter for probe output prefix

    args = parser.parse_args()

    logging.info("Started processing with arguments: %s", args)

    filtered_df = filter_hits(args.filter_output_file, args.hit_threshold)

    if filtered_df.empty:
        logging.warning("Filtered DataFrame is empty. No data to process. Exiting.")
        return

    probe_sets = classify_probes(filtered_df, args.classification_level)
    if not probe_sets:
        logging.warning("No probes were classified into any groups. Exiting.")
        return

    species_dict, sequence_lengths, full_sequences = determine_species(args.fasta_file, probe_sets)
    if not sequence_lengths:
        logging.warning("No sequences were found for classified probes. Exiting.")
        return

    coverage_results, probe_coverage_counts = calculate_coverage(probe_sets, args.blast_file, sequence_lengths)

    if filtered_df.empty or not coverage_results:
        logging.warning("No coverage results available to write. Skipping output.")
    else:
        logging.info(f"Writing final results to {args.output_file}...")
        with open(args.output_file, 'w') as output:
            output.write("PhyloLevelName\tTotalProbeCountInGroup\tCoverageQ1\tCoverageQ2\tCoverageQ3\n")
            for level, species_names in species_dict.items():
                logging.debug(f"Level is {level}, species_names is {species_names}")
                level_probe_count = len(probe_sets.get(level, set()))
                coverage_values = [
                    x[2] for x in coverage_results
                    if x[0] == level and x[1] in species_names
                ]
                if coverage_values:
                    q1 = pd.Series(coverage_values).quantile(0.25)
                    q2 = pd.Series(coverage_values).quantile(0.5)
                    q3 = pd.Series(coverage_values).quantile(0.75)
                else:
                    q1 = q2 = q3 = "NA"
                output.write(f"{level}\t{level_probe_count}\t{q1}\t{q2}\t{q3}\n")
    # Step 6: Save probe outputs for each PhyloLevel (if applicable)
    if not probe_sets:
        logging.warning("No probes classified into PhyloLevels. Skipping output of probe sequences.")
    else:
        save_probes_by_phylo_level(probe_sets, args.blast_file, args.output_probe_prefix, args.probe_fasta_file)
    logging.info("Processing completed.")

if __name__ == "__main__":
    main()

