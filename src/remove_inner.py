from Bio import SeqIO
import argparse

def parse_interval(header):
    """Extract interval from FASTA header.

    Args:
        header (str): FASTA header string containing a '_sliding:' segment.

    Returns:
        tuple: A tuple (start, end) of integers representing the interval.
    """
    sliding_part = header.split('_sliding:')[-1]
    start, end = map(int, sliding_part.split('-'))
    return start, end

def is_overlapping(current_start, current_end, prev_start, prev_end, next_start, next_end):
    """Check if the current interval overlaps with the previous or next interval.

    Args:
        current_start (int): Start of the current interval.
        current_end (int): End of the current interval.
        prev_start (int): Start of the previous interval.
        prev_end (int): End of the previous interval.
        next_start (int): Start of the next interval.
        next_end (int): End of the next interval.

    Returns:
        bool: True if the current interval overlaps with either neighbor, False otherwise.
    """
    return (current_start < prev_end and current_end > prev_start) or \
           (current_start < next_end and current_end > next_start)

def filter_probes(fasta_file, output_file):
    """Filter probes based on the described logic.

    Groups probes by source, then iteratively retains probes that do not overlap
    with previously retained probes or that skip over overlapping segments.

    Args:
        fasta_file (str): Path to the input FASTA file.
        output_file (str): Path to the output FASTA file.
    """
    records = list(SeqIO.parse(fasta_file, "fasta"))
    filtered_records = []

    # Group probes by source
    source_to_probes = {}
    for record in records:
        source = record.id.split('_sliding:')[0]
        if source not in source_to_probes:
            source_to_probes[source] = []
        source_to_probes[source].append(record)

    # Process probes for each source
    for source, probes in source_to_probes.items():
        if not probes:
            continue

        # Parse all intervals beforehand
        intervals = [parse_interval(probe.id) for probe in probes]

        # Always retain the first probe
        filtered_records.append(probes[0])
        prev_start, prev_end = intervals[0]

        # Process from the second probe
        i = 1
        while i < len(probes):
            current_start, current_end = intervals[i]

            # Check if the previous retained probe's end and the next probe's start overlap
            if i + 1 < len(probes):
                next_start, next_end = intervals[i + 1]
                if prev_end > next_start:
                    # If they overlap, proceed with the overlap logic
                    # Check if the current probe overlaps with the previous retained probe
                    if current_start < prev_end and current_end > prev_start:
                        # If it overlaps, check if it overlaps with the next probe
                        if i + 1 < len(probes):
                            next_start, next_end = intervals[i + 1]
                            if current_end <= next_start:
                                # Retain the current probe if it doesn't overlap with the next probe
                                filtered_records.append(probes[i])
                                # Update the previous interval
                                prev_start, prev_end = current_start, current_end
                        i += 1  # Skip the current probe if it overlaps with the next probe
                    else:
                        # Retain the current probe if it doesn't overlap with the previous retained probe
                        filtered_records.append(probes[i])
                        # Update the previous interval
                        prev_start, prev_end = current_start, current_end
                        i += 1
                else:
                    # If they do not overlap, retain the current probe
                    filtered_records.append(probes[i])
                    # Update the previous interval
                    prev_start, prev_end = current_start, current_end
                    i += 1
            else:
                # If there is no next probe, retain the current probe
                filtered_records.append(probes[i])
                i += 1

        # Always retain the last probe if it wasn't already retained
        if len(filtered_records) == 0 or filtered_records[-1].id != probes[-1].id:
            filtered_records.append(probes[-1])

    # Write filtered probes to the output file
    SeqIO.write(filtered_records, output_file, "fasta")

def main():
    """Parse command-line arguments and run the probe filtering process."""
    # Set up command-line argument parser
    parser = argparse.ArgumentParser(description="Filter overlapping probes from a FASTA file.")
    parser.add_argument("input_file", help="Path to the input FASTA file.")
    parser.add_argument("output_file", help="Path to the output FASTA file.")
    args = parser.parse_args()

    # Call the filter_probes function with the provided arguments
    filter_probes(args.input_file, args.output_file)

if __name__ == "__main__":
    main()

