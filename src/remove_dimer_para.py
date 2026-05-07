import sys
from Bio import SeqIO
from Bio.Seq import Seq
import primer3
from concurrent.futures import ThreadPoolExecutor, as_completed

def filter_sequence_by_n_content(sequence, max_n_proportion):
    """Filter sequence based on the maximum proportion of N content allowed.

    Args:
        sequence (SeqRecord): A BioPython SeqRecord object.
        max_n_proportion (float): Maximum allowed proportion of 'N' bases.

    Returns:
        bool: True if the N proportion is within the allowed threshold, False otherwise.
    """
    n_count = sequence.seq.count('N')
    n_proportion = n_count / len(sequence.seq)
    return n_proportion <= max_n_proportion

def generate_reverse_complements(probe_sequences, max_n_proportion):
    """Generate reverse complements for the given list of probe sequences.

    Args:
        probe_sequences (list): List of BioPython SeqRecord objects.
        max_n_proportion (float): Maximum allowed proportion of 'N' bases for filtering.

    Returns:
        list: List of reverse complement SeqRecord objects, filtered by N content.
    """
    reverse_complements = []
    for probe in probe_sequences:
        # Filter out sequences with too many 'N's
        if not filter_sequence_by_n_content(probe, max_n_proportion):
            continue
        rev_comp_seq = str(Seq(probe.seq).reverse_complement())
        rev_comp_record = SeqIO.SeqRecord(Seq(rev_comp_seq), id=probe.id, description="Reverse complement")
        reverse_complements.append(rev_comp_record)
    return reverse_complements

def check_dimerization(rev_comp_i, all_rev_complements, min_dimer_dG):
    """Check if a given reverse complement sequence is free of significant dimers.

    Args:
        rev_comp_i (SeqRecord): The reverse complement sequence to check.
        all_rev_complements (list): List of all reverse complement SeqRecord objects.
        min_dimer_dG (float): Minimum dimerization free energy threshold (kcal/mol).

    Returns:
        bool: True if no significant dimerization is found, False otherwise.
    """
    for rev_comp_j in all_rev_complements:
        if rev_comp_i.id != rev_comp_j.id:  # Ignore self-dimerization
            dimer_dG = primer3.calc_heterodimer(str(rev_comp_i.seq), str(rev_comp_j.seq), mv_conc=300).dg / 1000
            if dimer_dG < min_dimer_dG:
                return False
    return True

def calc_filter_dimers_parallel(rev_complements, min_dimer_dG, max_workers):
    """Filter reverse complement sequences based on dimer DG threshold using multithreading.

    Args:
        rev_complements (list): List of reverse complement SeqRecord objects.
        min_dimer_dG (float): Minimum dimerization free energy threshold (kcal/mol).
        max_workers (int): Maximum number of worker threads.

    Returns:
        list: Filtered list of reverse complement SeqRecord objects passing the dimer check.
    """
    filtered_complements = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(check_dimerization, rev_comp, rev_complements, min_dimer_dG): rev_comp for rev_comp in rev_complements}

        for future in as_completed(futures):
            rev_comp = futures[future]
            if future.result():
                filtered_complements.append(rev_comp)
    return filtered_complements

def main(input_fasta, output_fasta, min_dimer_dG, max_workers, max_n_proportion):
    """Main function to process the FASTA file and filter reverse complement sequences.

    Args:
        input_fasta (str): Path to the input FASTA file.
        output_fasta (str): Path to the output FASTA file.
        min_dimer_dG (float): Minimum dimerization free energy threshold (kcal/mol).
        max_workers (int): Maximum number of worker threads.
        max_n_proportion (float): Maximum allowed proportion of 'N' bases.
    """
    probe_sequences = list(SeqIO.parse(input_fasta, "fasta"))
    rev_complements = generate_reverse_complements(probe_sequences, max_n_proportion)
    filtered_complements = calc_filter_dimers_parallel(rev_complements, min_dimer_dG, max_workers)
    SeqIO.write(filtered_complements, output_fasta, "fasta")

    print(f"Filtered {len(rev_complements) - len(filtered_complements)} sequences due to dimerization issues.")
    print(f"Remaining {len(filtered_complements)} sequences written to {output_fasta}.")

if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python filter_revcomp_multithread.py <input_fasta> <output_fasta> <min_dimer_dG> <max_workers> <max_n_proportion>")
        sys.exit(1)

    input_fasta = sys.argv[1]
    output_fasta = sys.argv[2]
    min_dimer_dG = float(sys.argv[3])
    max_workers = int(sys.argv[4])
    max_n_proportion = float(sys.argv[5])

    main(input_fasta, output_fasta, min_dimer_dG, max_workers, max_n_proportion)

