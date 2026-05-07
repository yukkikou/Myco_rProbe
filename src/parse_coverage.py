import re
from collections import defaultdict
import argparse
import os
from Bio import SeqIO
import matplotlib.pyplot as plt

def load_species_info(species_file):
    """Load species information from a FASTA-style file with taxonomic lineages.

    Args:
        species_file (str): Path to the species information file.

    Returns:
        dict: A dictionary mapping sequence IDs to species names.
    """
    species_dict = {}
    try:
        with open(species_file, "r") as f:
            for line in f:
                if line.startswith(">"):
                    parts = line.strip().split()
                    if ";" in line:
                        seq_id = parts[0][1:]
                        species = " ".join(parts[1:])
                    else:
                        seq_id = parts[0][1:]
                        species = line.split("_5")[0]
                    species_dict[seq_id] = species
    except FileNotFoundError:
        print(f"Error: The file {species_file} was not found.")
    return species_dict

def extract_species(seq_id, species_dict_18S=None, species_dict_28S=None):
    """Extract species name from a sequence ID using available species dictionaries.

    Args:
        seq_id (str): The sequence identifier.
        species_dict_18S (dict, optional): Dictionary mapping 18S sequence IDs to species names.
        species_dict_28S (dict, optional): Dictionary mapping 28S sequence IDs to species names.

    Returns:
        str: The extracted species name, or "Unknown" if it cannot be determined.
    """
    if "::" in seq_id:
        return "_".join(re.split(r"[_.]", seq_id)[0:2])
    elif seq_id.startswith("sliva_"):
        if "18S" in seq_id and species_dict_18S:
            return species_dict_18S.get(seq_id, "Unknown")
        elif "28S" in seq_id and species_dict_28S:
            return species_dict_28S.get(seq_id, "Unknown")
    elif seq_id.startswith("unite_its_"):
        return seq_id.split("_")[2]
    return "Unknown"

def get_sequence_lengths(fasta_file):
    """Get sequence lengths from a FASTA file.

    Args:
        fasta_file (str): Path to the FASTA file.

    Returns:
        dict: A dictionary mapping sequence IDs to their lengths.
    """
    sequence_lengths = {}
    try:
        with open(fasta_file, "r") as file:
            for record in SeqIO.parse(file, "fasta"):
                sequence_id = record.id
                sequence_length = len(record.seq)
                sequence_lengths[sequence_id] = sequence_length
    except FileNotFoundError:
        print(f"Error: The file {fasta_file} was not found.")
    return sequence_lengths

def parse_blast_results(blast_file, fasta_file, species_dict_18S=None, species_dict_28S=None):
    """Parse BLAST results to compute coverage and region information per species.

    Args:
        blast_file (str): Path to the BLAST results file.
        fasta_file (str): Path to the FASTA file for sequence lengths.
        species_dict_18S (dict, optional): Dictionary for 18S species info.
        species_dict_28S (dict, optional): Dictionary for 28S species info.

    Returns:
        tuple: A tuple containing:
            - coverage_dict: Per-species per-sequence coverage data.
            - base_coverage_dict: Per-base coverage data.
            - covered_regions: Covered region intervals per species/sequence.
            - uncovered_regions: Uncovered region intervals per species/sequence.
    """
    sequence_lengths = get_sequence_lengths(fasta_file)
    coverage_dict = defaultdict(lambda: defaultdict(lambda: {"covered_positions": set(), "length": 0}))
    base_coverage_dict = defaultdict(lambda: defaultdict(list))

    try:
        with open(blast_file, "r") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 12:
                    print(f"Warning: Skipping invalid line: {line.strip()}")
                    continue
                target_id = parts[1]
                target_start, target_end = map(int, parts[8:10])
                species = extract_species(target_id, species_dict_18S, species_dict_28S)
                target_length = sequence_lengths.get(target_id, max(target_start, target_end))

                if not coverage_dict[species][target_id]["length"]:
                    coverage_dict[species][target_id]["length"] = target_length
                base_coverage = base_coverage_dict[species][target_id]
                if len(base_coverage) < target_length:
                    base_coverage.extend([0] * (target_length - len(base_coverage)))
                start, end = sorted((target_start, target_end))
                for pos in range(start, end + 1):
                    if 1 <= pos <= target_length:
                        coverage_dict[species][target_id]["covered_positions"].add(pos)
                        base_coverage[pos - 1] = 1
                    else:
                        print(f"Warning: pos {pos} over rnage，target_id {target_id} (length={target_length})")
    except FileNotFoundError:
        print(f"Error: The file {blast_file} was not found.")

    covered_regions = defaultdict(lambda: defaultdict(list))
    uncovered_regions = defaultdict(lambda: defaultdict(list))

    for species, seq_data in coverage_dict.items():
        for seq_id, info in seq_data.items():
            length = info["length"]
            covered_positions = info["covered_positions"]
            # Calculate covered regions
            start = None
            for pos in range(1, length + 1):
                if pos in covered_positions:
                    if start is None:
                        start = pos
                elif start is not None:
                    covered_regions[species][seq_id].append((start, pos - 1))
                    start = None
            if start is not None:
                covered_regions[species][seq_id].append((start, length))
            # Calculate uncovered regions
            start = None
            for pos in range(1, length + 1):
                if pos not in covered_positions:
                    if start is None:
                        start = pos
                elif start is not None:
                    uncovered_regions[species][seq_id].append((start, pos - 1))
                    start = None
            if start is not None:
                uncovered_regions[species][seq_id].append((start, length))

    return coverage_dict, base_coverage_dict, covered_regions, uncovered_regions

def save_regions_to_file(covered_regions, uncovered_regions, coverage_dict, file_path):
    """Save covered and uncovered region summaries to a tab-separated file.

    Args:
        covered_regions (dict): Covered region intervals per species/sequence.
        uncovered_regions (dict): Uncovered region intervals per species/sequence.
        coverage_dict (dict): Coverage data per species per sequence.
        file_path (str): Path to the output file.
    """
    with open(file_path, 'w') as f:
        f.write("rRNA_ID\tLength\tCovered_Regions_Length\tUncovered_Regions_Length\n")
        for species, seq_data in covered_regions.items():
            for seq_id, cov_regions in seq_data.items():
                total_covered_length = sum(end - start + 1 for start, end in cov_regions)
                total_uncovered_length = coverage_dict[species][seq_id]["length"] - total_covered_length
                f.write(f"{seq_id}\t{coverage_dict[species][seq_id]['length']}\t{total_covered_length}\t{total_uncovered_length}\n")

def calculate_coverage(coverage_dict):
    """Calculate coverage percentages for each species.

    Args:
        coverage_dict (dict): Coverage data per species per sequence.

    Returns:
        defaultdict: A dictionary mapping species names to lists of coverage percentages.
    """
    coverage_results = defaultdict(list)
    for species, seq_dict in coverage_dict.items():
        for seq_id, coverage_data in seq_dict.items():
            covered_positions = coverage_data["covered_positions"]
            total_length = coverage_data["length"]
            coverage_percentage = (len(covered_positions) / total_length) * 100 if total_length else 0
            coverage_results[species].append(coverage_percentage)
    return coverage_results

def write_total_coverage(coverage_dict, output_file):
    """Write total coverage percentages for each sequence to a file.

    Args:
        coverage_dict (dict): Coverage data per species per sequence.
        output_file (str): Path to the output file.
    """
    with open(output_file, "w") as outfile:
        outfile.write("Species\tSequence_ID\tTotal_Coverage(%)\n")
        for species, seq_dict in coverage_dict.items():
            for seq_id, coverage_data in seq_dict.items():
                total_length = coverage_data["length"]
                covered_positions = coverage_data["covered_positions"]
                coverage_percentage = (len(covered_positions) / total_length) * 100 if total_length else 0
                outfile.write(f"{species}\t{seq_id}\t{coverage_percentage:.2f}\n")

def plot_total_coverage(coverage_results, output_dir):
    """Plot a histogram of total coverage distribution across all sequences.

    Args:
        coverage_results (defaultdict): Dictionary mapping species names to lists of coverage percentages.
        output_dir (str): Directory to save the plot image.
    """
    all_coverage = [coverage for coverage_list in coverage_results.values() for coverage in coverage_list]
    plt.figure(figsize=(10, 6))
    plt.hist(all_coverage, bins=50, color="blue", edgecolor="black")
    plt.title("Total Coverage Distribution")
    plt.xlabel("Coverage (%)")
    plt.ylabel("Frequency")
    plt.grid(True)
    plt.savefig(os.path.join(output_dir, "total_coverage_distribution.png"))
    plt.close()

def write_per_base_coverage(base_coverage_dict, output_file):
    """Write per-base coverage data to a tab-separated file.

    Args:
        base_coverage_dict (dict): Per-base coverage data per species per sequence.
        output_file (str): Path to the output file.
    """
    with open(output_file, "w") as f:
        f.write("Species\tSequence_ID\tPosition\tCoverage\n")
        for species, seq_dict in base_coverage_dict.items():
            for seq_id, base_coverage in seq_dict.items():
                for pos, coverage in enumerate(base_coverage):
                    f.write(f"{species}\t{seq_id}\t{pos+1}\t{coverage}\n")

def plot_summary_coverage(file_path, output_dir):
    """Plot a histogram of coverage percentage from a summary file.

    Args:
        file_path (str): Path to the rRNA coverage summary file.
        output_dir (str): Directory to save the plot image.
    """
    try:
        coverage_percentages = []
        with open(file_path, "r") as f:
            # Skip the header
            next(f)
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 4:
                    try:
                        length = int(parts[1])
                        covered_length = int(parts[2])
                        coverage_percentage = (covered_length / length) * 100 if length else 0
                        coverage_percentages.append(coverage_percentage)
                    except ValueError:
                        print(f"Warning: Could not parse line: {line.strip()}")

        plt.figure(figsize=(10, 6))
        plt.hist(coverage_percentages, bins=50, color="green", edgecolor="black")
        plt.title("Covered Regions Length / Total Length (%) Distribution")
        plt.xlabel("Coverage (%)")
        plt.ylabel("Frequency")
        plt.grid(True)
        plt.savefig(os.path.join(output_dir, "regions_length_coverage_distribution.png"))
        plt.close()
    except FileNotFoundError:
        print(f"Error: The file {file_path} was not found.")

def main():
    """Parse command-line arguments and run the rRNA coverage analysis workflow."""
    parser = argparse.ArgumentParser(description="Analyze rRNA coverage from BLAST results.")
    parser.add_argument("--blast", required=True, help="Path to the BLAST results file.")
    parser.add_argument("--fasta", required=True, help="Path to the FASTA file for sequence lengths.")
    parser.add_argument("--species_18S", required=False, help="Path to the 18S species information file.")
    parser.add_argument("--species_28S", required=False, help="Path to the 28S species information file.")
    parser.add_argument("--outputdir", required=True, help="Directory to save output files and plots.")
    args = parser.parse_args()

    os.makedirs(args.outputdir, exist_ok=True)

    species_dict_18S = load_species_info(args.species_18S) if args.species_18S else None
    species_dict_28S = load_species_info(args.species_28S) if args.species_28S else None

    coverage_dict, base_coverage_dict, covered_regions, uncovered_regions = parse_blast_results(
        args.blast, args.fasta, species_dict_18S, species_dict_28S)

    summary_file_path = os.path.join(args.outputdir, "rRNA_coverage_summary.txt")
    save_regions_to_file(covered_regions, uncovered_regions, coverage_dict, os.path.join(args.outputdir, "rRNA_coverage_summary.txt"))

    #coverage_results = calculate_coverage(coverage_dict)

    #write_total_coverage(coverage_dict, os.path.join(args.outputdir, "total_coverage.txt"))

    #plot_total_coverage(coverage_results, args.outputdir)
    plot_summary_coverage(summary_file_path, args.outputdir)

    # Write per-base coverage results
    #per_base_output_file = os.path.join(args.outputdir, "per_base_coverage.txt")
    #write_per_base_coverage(base_coverage_dict, per_base_output_file)

if __name__ == "__main__":
    main()

