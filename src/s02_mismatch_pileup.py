import sys
from collections import defaultdict
import csv
import matplotlib.pyplot as plt
import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

IUPAC_AMBIGUITY_CODES = {
    ('A',): 'A',
    ('C',): 'C',
    ('G',): 'G',
    ('T',): 'T',
    ('A', 'G'): 'R',
    ('C', 'T'): 'Y',
    ('G', 'C'): 'S',
    ('A', 'T'): 'W',
    ('G', 'T'): 'K',
    ('A', 'C'): 'M',
    ('C', 'G', 'T'): 'B',
    ('A', 'G', 'T'): 'D',
    ('A', 'C', 'G'): 'H',
    ('A', 'C', 'T'): 'V',
    ('A', 'C', 'G', 'T'): 'N'
}


def analyze_pileup(pileup_file, mismatch_threshold=0.3):  # Added threshold
    """
    Analyzes a pileup file to differentiate between matching and mismatching
    bases, calculate statistics, and generate a consensus sequence with
    IUPAC ambiguity codes when mismatch rate exceeds a threshold.
    Args:
        pileup_file (str): Path to the pileup file.
        mismatch_threshold (float): Threshold for mismatch rate (0-1).

    Returns:
    A tuple containing:
        - dict: A dictionary of dictionaries (same as before).
        - str: The consensus sequence as a string.
    """


    pileup_results = defaultdict(dict)
    consensus_sequence = ""

    with open(pileup_file, "r") as pileup_file_handle:  # Using pileup_file_handle
      for line in pileup_file_handle:
        chrom, pos, ref, coverage, bases, qualities = line.strip().split('\t')
        pos = int(pos)

        pileup_results[pos]['reference_base'] = ref.upper()
        pileup_results[pos]['total_coverage'] = int(coverage)
        pileup_results[pos]['base_counts'] = defaultdict(int)  # Initialize counts
        total_quality = 0

        #Parse base and qualities (robust indel handling)
        base_index = 0
        bases_len = len(bases)
        qualities_len = len(qualities)

        while base_index < bases_len:
            base = bases[base_index]

            #Match to reference
            if base == ".":
                base = ref.upper()
            elif base == ",":
                base = ref.lower()
            # Handle insertions/deletions
            elif base in "+-":
                indel_length = 0
                indel_length_str = ""
                index = base_index + 1
                while index < bases_len and bases[index].isdigit():
                    indel_length_str += bases[index]
                    index += 1
                indel_length = int(indel_length_str)
                inserted_sequence = bases[index : index + indel_length]

                base_index = index + indel_length  # Index advances past indel
                continue #Go to next base

            elif base == "^": #Start of qualified segment (ignore)
                base_index+=2 #Skip the next char - mapping quality
                continue
            elif base == "$": #End of read segment
                base_index+=1 #Skip
                continue

            #Process base and quality
            pileup_results[pos]['base_counts'][base.upper()] += 1 #Aggregate to uppercase
            if qualities != "*" and base_index < qualities_len: # '*' means mapping quality is unavailable - but also no called base!
              total_quality += ord(qualities[base_index]) - 33 #Convert ASCII to quality score
            base_index += 1 #Advance base


        if int(coverage) != 0:
            pileup_results[pos]['average_quality'] = total_quality / int(coverage) #Avoid div/0
        else:
            pileup_results[pos]['average_quality'] = 0 #No coverage

        # Consensus base determination
        if int(coverage) > 0:
            base_counts = pileup_results[pos]['base_counts']
            total_reads = int(coverage)  # Use total coverage as total reads

            sorted_bases = sorted(base_counts.items(), key=lambda item: item[1], reverse=True)
            top_base = sorted_bases[0][0] if sorted_bases else ref.upper() # Default to reference base

            #Calculate top two bases
            if len(sorted_bases) > 1:
                second_base = sorted_bases[1][0]
                top_base_count = sorted_bases[0][1]
                second_base_count = sorted_bases[1][1]

                combined_rate = (top_base_count + second_base_count) / total_reads
                top_rate = top_base_count / total_reads

                #Ambiguous rules
                if top_rate < 1 - mismatch_threshold and combined_rate >= 1 - mismatch_threshold :
                    bases_to_combine = tuple(sorted([top_base,second_base]))
                    consensus_base = IUPAC_AMBIGUITY_CODES.get(bases_to_combine, 'N') #Fallback to N
                else:
                  consensus_base = top_base
            else:
                #More stringent rule for only one base present.
                rate = sorted_bases[0][1] / total_reads if sorted_bases else 0
                if rate < 1 - mismatch_threshold:
                    consensus_base = 'N' #If its very low and not clear, use N.
                else:
                    consensus_base = top_base #Otherwise just this base.

            consensus_sequence += consensus_base
        else:
            consensus_sequence += "N"  # No coverage - ambiguous

            #Remove ambiguity for a better N masking.

    return pileup_results, consensus_sequence


def load_depth_data(depth_file):
    """Loads depth information from a depth file."""
    depth_data = {}
    with open(depth_file, "r") as f:
        for line in f:
            chrom, pos, depth = line.strip().split('\t')
            pos = int(pos)
            depth = int(depth)
            depth_data[pos] = depth
    return depth_data


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python analyze_pileup.py <pileup_file> <depth_file>")
        sys.exit(1)

    pileup_file = sys.argv[1]
    depth_file = sys.argv[2]

    results, consensus_sequence = analyze_pileup(pileup_file, mismatch_threshold=0.3) #Call with threshold.
    depth_data = load_depth_data(depth_file)

    # Prepare data for TSV output
    output_data = []
    positions = []
    matching_rates = []
    mismatching_rates = []

    for position, data in results.items():
        reference_base = data['reference_base']
        total_coverage = data['total_coverage']
        base_counts = data['base_counts']
        average_quality = data['average_quality']

        matching_reads = data['base_counts'].get(data['reference_base'].upper(), 0) + data['base_counts'].get(data['reference_base'].lower(), 0)
        mismatching_reads = data['total_coverage'] - matching_reads

        if data['total_coverage'] > 0:
            matching_rate = (matching_reads / data['total_coverage']) * 100
            mismatching_rate = (mismatching_reads / data['total_coverage']) * 100
        else:
            matching_rate = 0
            mismatching_rate = 0

        output_data.append([
            position,
            reference_base,
            total_coverage,
            str(base_counts),  # Convert dict to string for TSV
            f"{average_quality:.2f}",  # Format to two decimal places
            f"{matching_rate:.2f}",
            f"{mismatching_rate:.2f}"
        ])

        positions.append(position)
        matching_rates.append(matching_rate)
        mismatching_rates.append(mismatching_rate)

    # Write to TSV file
    output_file = "pileup_analysis.tsv"
    with open(output_file, "w", newline="") as tsvfile:
        tsv_writer = csv.writer(tsvfile, delimiter="\t")
        tsv_writer.writerow([
            "Position",
            "Reference Base",
            "Total Coverage",
            "Base Counts",
            "Average Quality",
            "Matching Rate",
            "Mismatching Rate"
        ])  # Write header
        tsv_writer.writerows(output_data)

    print(f"Analysis complete. Results saved to {output_file}")

    # Create the matching rate plot
    plt.figure(figsize=(12, 6))
    plt.plot(positions, matching_rates)
    plt.xlabel("Genomic Position")
    plt.ylabel("Matching Rate (%)")
    plt.title("Matching Rate Across Genomic Positions")
    plt.grid(True)
    plt.ylim(0, 100)
    step = max(1, len(positions) // 20)
    plt.xticks(positions[::step], rotation=45, ha="right")
    plt.tight_layout()
    plot_file_matching = "matching_rate.png"
    plt.savefig(plot_file_matching)
    print(f"Matching rate plot saved to {plot_file_matching}")

    # Pileup Chart: Create the coverage and mismatch rate plot with shared X-axis
    fig, ax1 = plt.subplots(figsize=(12, 6))

    color = 'tab:blue'
    ax1.set_xlabel("Genomic Position")
    ax1.set_ylabel("Total Coverage", color=color)
    total_coverages_pileup = [results[p]['total_coverage'] for p in positions] #Extract coverage values from pileup results
    ax1.plot(positions, total_coverages_pileup, color=color)
    ax1.tick_params(axis='y', labelcolor=color)

    ax2 = ax1.twinx()  # Instantiate a second axes that shares the same x-axis

    color = 'tab:red'
    ax2.set_ylabel("Mismatching Rate (%)", color=color)  # We already handled the x-label with ax1
    ax2.plot(positions, mismatching_rates, color=color)
    ax2.tick_params(axis='y', labelcolor=color)

    fig.tight_layout()  # Otherwise the right y-label is slightly clipped
    plt.title("Pileup: Coverage and Mismatching Rate Across Genomic Positions")

    # Make sure x-axis labels are not overcrowded
    plt.xticks(positions[::step], rotation=45, ha="right") # use the positions and old step

    plot_file_pileup_coverage_mismatch = "pileup_coverage_vs_mismatch.png"
    plt.savefig(plot_file_pileup_coverage_mismatch)
    print(f"Pileup Coverage vs mismatch plot saved to {plot_file_pileup_coverage_mismatch}")
    plt.close('all')#Clean

    # DEPTH Chart: Create the coverage and mismatch rate plot with shared X-axis
    fig, ax1 = plt.subplots(figsize=(12, 6))

    color = 'tab:blue'
    ax1.set_xlabel("Genomic Position")
    ax1.set_ylabel("Total Coverage", color=color)

    # Use positions from pileup and depth data from depth file
    depth_positions = sorted(depth_data.keys()) # ensure plots in ascending order.

    total_coverages = [depth_data.get(pos, 0) for pos in depth_positions]
    ax1.plot(depth_positions, total_coverages, color=color)
    ax1.tick_params(axis='y', labelcolor=color)

    ax2 = ax1.twinx()  # Instantiate a second axes that shares the same x-axis

    color = 'tab:red'
    ax2.set_ylabel("Mismatching Rate (%)", color=color)  # We already handled the x-label with ax1
    #The depth position needs to be a shared value:
    mismatching_rates_aligned = [mismatching_rates[positions.index(pos)] if pos in positions else 0 for pos in depth_positions]
    ax2.plot(depth_positions, mismatching_rates_aligned, color=color)
    ax2.tick_params(axis='y', labelcolor=color)

    fig.tight_layout()  # Otherwise the right y-label is slightly clipped
    plt.title("Depth Coverage and Mismatching Rate Across Genomic Positions")

    # Make sure x-axis labels are not overcrowded
    plt.xticks(depth_positions[::step], rotation=45, ha="right") # use the new positions and old step

    plot_file_coverage_mismatch = "depth_coverage_vs_mismatch.png"
    plt.savefig(plot_file_coverage_mismatch)
    print(f"Depth Coverage vs mismatch plot saved to {plot_file_coverage_mismatch}")

    plt.close('all')#Clean

    # Write consensus sequence to FASTA file
    fasta_file = "consensus.fasta"
    record = SeqRecord(Seq(consensus_sequence), id="consensus", description=f"Generated from {pileup_file} with mismatch threshold {0.3}") #Pass it in.
    SeqIO.write(record, fasta_file, "fasta")
    print(f"Consensus sequence saved to {fasta_file}")

