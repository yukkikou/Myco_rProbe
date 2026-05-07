import os
import re
from pathlib import Path
from collections import defaultdict
import csv
import argparse
import sys # For redirecting stdout

# --- 1. Parse Taxonomy (Same as your previous working version) ---
def parse_taxonomy(nodes_file, info_file):
    """Parse taxonomy nodes and info files into a unified tax_data dictionary.

    Reads NCBI-style nodes.dmp and names.dmp files to build a comprehensive
    taxonomy tree with parent-child relationships.

    Args:
        nodes_file (str): Path to the taxonomy nodes file (nodes.dmp format).
        info_file (str): Path to the taxonomy info/names file (names.dmp format).

    Returns:
        tuple: A tuple containing:
            - tax_data: Dictionary mapping tax IDs to their name, rank, parent_id, and children.
            - name_to_tax_id: DefaultDict mapping taxonomic names to lists of matching tax IDs.
    """
    tax_data = {}
    name_to_tax_id = defaultdict(list)
    # print(f"Parsing taxonomy nodes from: {nodes_file}")
    with open(nodes_file, 'r') as f:
        for line in f:
            parts = [p.strip() for p in line.split('|')]
            if len(parts) < 3: continue
            tax_id, parent_id, rank = parts[0], parts[1], parts[2]
            tax_data[tax_id] = {'name': '', 'rank': rank, 'parent_id': parent_id, 'children': []}
    # print(f"Parsing taxonomy info from: {info_file}")
    with open(info_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if len(row) < 3: continue
            tax_id, _, name = row[0].strip(), row[1].strip(), row[2].strip()
            if tax_id in tax_data:
                tax_data[tax_id]['name'] = name
                name_to_tax_id[name].append(tax_id)
    for tax_id, data in tax_data.items():
        parent_id = data['parent_id']
        if parent_id in tax_data and tax_id != parent_id:
            tax_data[parent_id]['children'].append(tax_id)
    # print(f"Parsed {len(tax_data)} taxonomy entries.")
    return tax_data, name_to_tax_id

# --- 2. Parse Probe Conservation Data (Using CoverageQ2) ---
def parse_conservation_data(grouped_tsv_dir):
    """Parse probe conservation data from grouped TSV files in a directory.

    Reads all *_grouped.tsv files and extracts probe counts and coverage
    statistics (CoverageQ2) for each taxonomic rank/name combination.

    Args:
        grouped_tsv_dir (str or Path): Directory containing *_grouped.tsv files.

    Returns:
        dict: Dictionary mapping (rank, taxon_name) tuples to dicts with 'probes' and 'coverage' keys.
    """
    conservation_info = {}
    # print(f"Parsing conservation data from: {grouped_tsv_dir}")
    for tsv_file in Path(grouped_tsv_dir).glob("*_grouped.tsv"):
        match = re.search(r"_([a-zA-Z\s]+)_grouped\.tsv$", tsv_file.name)
        if not match:
            match_simple = re.search(r"_(class|order|subclass|subkingdom|subphylum|family|subfamily|suborder|genus|species)_grouped\.tsv$", tsv_file.name)
            if not match_simple:
                # print(f"Could not determine rank from filename: {tsv_file.name}")
                continue
            rank = match_simple.group(1)
        else:
            rank = match.group(1).replace("_", " ")

        with open(tsv_file, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            if 'CoverageQ2' not in reader.fieldnames:
                # print(f"Warning: 'CoverageQ2' column not found in {tsv_file.name}. Skipping coverage for this file.")
                has_q2 = False
            else:
                has_q2 = True

            for row in reader:
                taxon_name = row['PhyloLevelName'].strip()
                try:
                    probes = int(row['TotalProbeCountInGroup'])
                    coverage = float(row['CoverageQ2']) if has_q2 and row.get('CoverageQ2') else 0.0
                    conservation_info[(rank, taxon_name)] = {'probes': probes, 'coverage': coverage}
                except (ValueError, KeyError) as e:
                    print(f"Skipping row due to parsing error in {tsv_file.name} for {taxon_name}: {e}, Row: {row}", file=sys.stderr)
    # print(f"Parsed conservation data for {len(conservation_info)} taxon groups.")
    return conservation_info

# --- 3. Helper: Parse Species Map FASTA (New Function) ---
def parse_species_map_fasta(map_fasta_file):
    """Parse a FASTA file mapping sequence IDs to species names via their lineage headers.

    FASTA headers are expected to contain an ID followed by a semicolon-separated
    taxonomic lineage, with the species name as the last element.

    Args:
        map_fasta_file (str): Path to the species map FASTA file.

    Returns:
        dict: Dictionary mapping sequence IDs to refined species names (e.g., "Genus species").
    """
    id_to_species_map = {}
    if not map_fasta_file or not Path(map_fasta_file).exists():
        # print("No species map FASTA provided or file not found. Skipping.")
        return id_to_species_map

    # print(f"Parsing species map FASTA from: {map_fasta_file}")
    with open(map_fasta_file, 'r') as f:
        for line in f:
            if line.startswith(">"):
                header_content = line[1:].strip()
                parts = header_content.split(maxsplit=1) # Split only on the first space
                if len(parts) < 2:
                    # print(f"Warning: Malformed header in map FASTA: {header_content}")
                    continue

                seq_id = parts[0] # e.g., sliva_28S_CBZX020000007.2.3306
                lineage_str = parts[1] # e.g., Eukaryota;...;Candida parapsilosis

                # Extract species name from lineage (typically the last part)
                species_name_candidate = lineage_str.split(';')[-1].strip()

                # Further refine to get "Genus species" if possible
                # This helps match taxonomy which usually stores "Genus species"
                name_parts = species_name_candidate.split()
                if len(name_parts) >= 2 and name_parts[0][0].isupper() and name_parts[1].islower():
                    refined_species_name = f"{name_parts[0]} {name_parts[1]}"
                    id_to_species_map[seq_id] = refined_species_name
                # else:
                    # print(f"    Warning: Could not refine species name from '{species_name_candidate}' for ID {seq_id}. Using as is or might fail lookup.")
                    # id_to_species_map[seq_id] = species_name_candidate # Or skip if not confident
    # print(f"Parsed {len(id_to_species_map)} entries from species map FASTA.")
    return id_to_species_map

# --- 4. Parse Probe FASTA File (REVISED for dual input types) ---
def parse_probes_fasta(fasta_file, tax_data, name_to_tax_id, id_species_map=None):
    """Parse probe FASTA file and map probes to species tax IDs.

    Supports two types of probe headers: Type 1 (embedded Genus_species pattern)
    and Type 2 (mapped via an external species map).

    Args:
        fasta_file (str): Path to the probe FASTA file.
        tax_data (dict): Taxonomy data dictionary from parse_taxonomy().
        name_to_tax_id (defaultdict): Mapping of names to tax IDs.
        id_species_map (dict, optional): Optional mapping of sequence IDs to species names.

    Returns:
        tuple: A tuple containing:
            - species_probes: DefaultDict mapping tax IDs to lists of probe headers.
            - probes_sequences: Dictionary mapping probe headers to their sequences.
    """
    species_probes = defaultdict(list) # tax_id -> list of probe headers
    probes_sequences = {} # header -> sequence

    # print(f"Parsing probes FASTA from: {fasta_file}")
    current_header = None
    current_seq_parts = []

    species_name_to_tax_id_map = {} # Local cache for faster lookups
    for tax_id_val, data in tax_data.items():
        if data['rank'] == 'species' and data['name']:
            species_name_to_tax_id_map[data['name'].lower()] = tax_id_val

    def get_species_name_from_type1_header(header_line):
        # Type 1: e.g., Aspergillus_fumigatusa1163.ASM15014v1_5S_rRNA::...
        main_id_part = header_line.split("::")[0]
        # Attempt to extract "Genus_species" or "Genus species"
        # Looking for Genus_species pattern like Aspergillus_fumigatus
        match_gs_underscore = re.match(r"([A-Z][a-z]+)_([a-z]+)", main_id_part)
        if match_gs_underscore:
            return f"{match_gs_underscore.group(1)} {match_gs_underscore.group(2)}" # Convert to "Genus species"

        # If no underscore, try to parse from dot-separated parts, e.g. A.fumigatus might be Genus.species
        species_candidate_part = main_id_part.split('.')[0] if '.' in main_id_part else main_id_part
        species_candidate_part = species_candidate_part.replace('_', ' ')

        match_gs_space = re.match(r"([A-Z][a-z]+)\s+([a-z]+)", species_candidate_part)
        if match_gs_space:
            return f"{match_gs_space.group(1)} {match_gs_space.group(2)}"

        parts = species_candidate_part.split(' ')
        if len(parts) >= 2 and parts[0][0].isupper() and parts[1].islower():
            return f"{parts[0]} {parts[1]}"
        return None

    if not Path(fasta_file).exists():
        print(f"Warning: Probe FASTA file not found: {fasta_file}", file=sys.stderr)
        return species_probes, probes_sequences

    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line: continue
            if line.startswith(">"):
                if current_header and current_seq_parts:
                    probes_sequences[current_header] = "".join(current_seq_parts)

                    species_name_str = None
                    # Attempt 1: Parse Type 1 header
                    species_name_str = get_species_name_from_type1_header(current_header)

                    # Attempt 2: Use id_species_map if Type 1 failed and map exists
                    if not species_name_str and id_species_map:
                        # Extract the ID part that would match the map, e.g., "sliva_28S_CBZX020000007.2.3306"
                        # from "sliva_28S_CBZX020000007.2.3306_sliding:1-60"
                        probe_id_for_map = current_header.split('_sliding:')[0] # Common pattern
                        if not probe_id_for_map: # fallback if _sliding not present
                             probe_id_for_map = current_header.split('::')[0] # Another common pattern

                        species_name_str = id_species_map.get(probe_id_for_map)
                        # if species_name_str:
                        #     print(f"    Debug: Found '{species_name_str}' for ID '{probe_id_for_map}' using map.")
                        # elif probe_id_for_map in id_species_map: # Should not happen if .get used
                        #     print(f"    Debug: ID '{probe_id_for_map}' in map, but value is None or empty.")
                        # else:
                        #     print(f"    Debug: Probe ID '{probe_id_for_map}' (from header '{current_header}') not in id_species_map.")


                    if species_name_str:
                        species_tax_id = species_name_to_tax_id_map.get(species_name_str.lower())
                        if species_tax_id:
                            species_probes[species_tax_id].append(current_header)
                        # else:
                            # print(f"    Warning: Species '{species_name_str}' (from header/map) not found in taxonomy for: {current_header}", file=sys.stderr)
                    # else:
                        # print(f"    Warning: Could not determine species for header: {current_header}", file=sys.stderr)

                current_header = line[1:]
                current_seq_parts = []
            else:
                if current_header: current_seq_parts.append(line)

        if current_header and current_seq_parts: # Last sequence
            probes_sequences[current_header] = "".join(current_seq_parts)
            species_name_str = None
            species_name_str = get_species_name_from_type1_header(current_header)
            if not species_name_str and id_species_map:
                probe_id_for_map = current_header.split('_sliding:')[0]
                if not probe_id_for_map:
                     probe_id_for_map = current_header.split('::')[0]
                species_name_str = id_species_map.get(probe_id_for_map)
            if species_name_str:
                species_tax_id = species_name_to_tax_id_map.get(species_name_str.lower())
                if species_tax_id: species_probes[species_tax_id].append(current_header)
                # else:
                    # print(f"    Warning: Species '{species_name_str}' (from header/map) not found in taxonomy for: {current_header}", file=sys.stderr)
            # else:
                # print(f"    Warning: Could not determine species for header: {current_header}", file=sys.stderr)

    # print(f"Parsed {len(probes_sequences)} probes from FASTA. Mapped probes for {len(species_probes)} species.")
    return species_probes, probes_sequences

# --- 5. Build and Display Decision Tree (Same as your previous working version) ---
def display_decision_tree(tax_id, tax_data, conservation_info, species_probes,
                          indent_level=0, accumulated_probe_designs_from_groups=0,
                          probes_sequences=None): # probes_sequences is available but not used for seq printing
    """Recursively build and display a hierarchical taxonomic decision tree.

    Prints the taxonomic tree starting from the given tax ID, annotating nodes
    with probe counts and coverage information from group files and FASTA probes.

    Args:
        tax_id (str): Root tax ID to start the tree from.
        tax_data (dict): Taxonomy data dictionary.
        conservation_info (dict): Conservation/probe data per (rank, name).
        species_probes (defaultdict): Mapping of tax IDs to lists of probe headers.
        indent_level (int): Current indentation level for tree display.
        accumulated_probe_designs_from_groups (int): Running total of probe designs encountered.
        probes_sequences (dict, optional): Dictionary of probe sequences (unused for display, retained for compatibility).

    Returns:
        bool: True if this node or any of its descendants were printed, False otherwise.
    """
    node = tax_data.get(tax_id)
    if not node:
        return False

    node_name = node['name']
    node_rank = node['rank']

    if not node_name:
        any_child_printed_for_unnamed_node = False
        children_ids_unnamed = node.get('children', [])
        sorted_children_data_unnamed = []
        for child_id_unnamed in children_ids_unnamed:
            child_node_unnamed = tax_data.get(child_id_unnamed)
            if child_node_unnamed:
                 sorted_children_data_unnamed.append({'id': child_id_unnamed,
                                                      'name': child_node_unnamed.get('name', f"Unnamed_TaxID_{child_id_unnamed}")})
        sorted_children_data_unnamed.sort(key=lambda x: x['name'])
        for child_data_unnamed in sorted_children_data_unnamed:
            if display_decision_tree(child_data_unnamed['id'], tax_data, conservation_info, species_probes,
                                     indent_level,
                                     accumulated_probe_designs_from_groups,
                                     probes_sequences):
                any_child_printed_for_unnamed_node = True
        return any_child_printed_for_unnamed_node

    level_specific_group_info = conservation_info.get((node_rank, node_name))
    probes_from_this_group_file = 0
    if level_specific_group_info:
        probes_from_this_group_file = level_specific_group_info.get('probes', 0)

    is_species_with_fasta_probes = (node_rank == 'species' and tax_id in species_probes and species_probes[tax_id])
    current_total_accumulated_designs_from_groups = accumulated_probe_designs_from_groups + probes_from_this_group_file
    this_node_should_be_printed = bool(level_specific_group_info) or is_species_with_fasta_probes

    if this_node_should_be_printed:
        current_path_element = f"{node_name} ({node_rank} - TaxID: {tax_id})"
        print(f"{'  ' * indent_level}- {current_path_element}")

        if level_specific_group_info:
            coverage_str = f"{level_specific_group_info.get('coverage', 0.0):.2f}%"
            print(f"{'  ' * (indent_level + 1)}  Probes from its group file: {probes_from_this_group_file} designs, Group Coverage: {coverage_str}")

        print(f"{'  ' * (indent_level + 1)}  Accumulated designs from group files on this path (incl. this level): {current_total_accumulated_designs_from_groups}")

        if is_species_with_fasta_probes:
            probes_for_this_species_headers = species_probes[tax_id]
            print(f"{'  ' * (indent_level + 1)}  Actual FASTA probes for this species ({len(probes_for_this_species_headers)}):")
            for i, probe_header in enumerate(probes_for_this_species_headers):
                if i < 3:
                    print(f"{'  ' * (indent_level + 2)}  {probe_header}")
                elif i == 3:
                    print(f"{'  ' * (indent_level + 2)}  ...and {len(probes_for_this_species_headers) - 3} more probes.")
                    break

    any_child_printed = False
    children_ids = node.get('children', [])
    sorted_children_data = []
    for child_id in children_ids:
        child_node = tax_data.get(child_id)
        if child_node:
             sorted_children_data.append({'id': child_id, 'name': child_node.get('name', f"Unnamed_TaxID_{child_id}")})
    sorted_children_data.sort(key=lambda x: x['name'])
    child_indent_level = indent_level + 1 if this_node_should_be_printed else indent_level

    for child_data in sorted_children_data:
        if display_decision_tree(child_data['id'], tax_data, conservation_info, species_probes,
                                 child_indent_level,
                                 current_total_accumulated_designs_from_groups,
                                 probes_sequences):
            any_child_printed = True

    return this_node_should_be_printed or any_child_printed


# --- Main Execution (with argparse) ---
def main():
    """Parse command-line arguments and generate the hierarchical taxonomic report."""
    parser = argparse.ArgumentParser(description="Generate a hierarchical taxonomic report with probe information.")
    parser.add_argument("--nodes", required=True, help="Path to the taxonomy nodes.tsv file (e.g., nodes.dmp format).")
    parser.add_argument("--info", required=True, help="Path to the taxonomy info/names.tsv file (e.g., names.dmp format).")
    parser.add_argument("--grouped-dir", required=True, help="Path to the directory containing *_grouped.tsv files.")
    parser.add_argument("--probes-fasta", required=True, help="Path to the main FASTA file with probe sequences.")
    parser.add_argument("--species-map-fasta", help="Optional: Path to a FASTA file mapping generic probe IDs to species lineages (for Type 2 probe headers).")
    parser.add_argument("--fungi-root-name", default="Fungi", help="Name of the root taxon to start processing from (default: Fungi).")
    parser.add_argument("--output-file", help="Optional: Path to an output file. If not provided, output goes to stdout.")

    args = parser.parse_args()

    # Prepare output stream
    original_stdout = sys.stdout
    output_fh = None
    if args.output_file:
        try:
            output_fh = open(args.output_file, 'w')
            sys.stdout = output_fh
            print(f"Outputting to file: {args.output_file}")
        except IOError as e:
            sys.stdout = original_stdout # Revert to stdout on error
            print(f"Error: Could not open output file {args.output_file}: {e}", file=sys.stderr)
            # Decide if you want to exit or continue to stdout
            # exit(1)
            print("Continuing with output to stdout.", file=sys.stderr)


    print("Starting script with provided arguments...", file=sys.stderr if not args.output_file else original_stdout) # Progress to stderr or original stdout
    # Use Path objects for file paths
    path_to_taxonomy_nodes = Path(args.nodes)
    path_to_taxonomy_info = Path(args.info)
    path_to_grouped_tsv_dir = Path(args.grouped_dir)
    path_to_probes_fasta = Path(args.probes_fasta)
    path_to_species_map_fasta = Path(args.species_map_fasta) if args.species_map_fasta else None


    if not path_to_taxonomy_nodes.exists() or not path_to_taxonomy_info.exists():
        print(f"Error: Taxonomy files not found. Please check --nodes and --info paths.", file=sys.stderr)
        if output_fh: output_fh.close(); sys.stdout = original_stdout
        exit(1)

    tax_data, name_to_tax_id = parse_taxonomy(path_to_taxonomy_nodes, path_to_taxonomy_info)
    # print(f"Taxonomy parsed: {len(tax_data)} entries.", file=sys.stderr if not args.output_file else original_stdout)

    conservation_info = parse_conservation_data(path_to_grouped_tsv_dir)
    # if not conservation_info: print(f"Warning: No conservation data loaded from {path_to_grouped_tsv_dir.resolve()}", file=sys.stderr)
    # else: print(f"Conservation data parsed: {len(conservation_info)} groups.", file=sys.stderr if not args.output_file else original_stdout)

    id_species_map = {}
    if path_to_species_map_fasta and path_to_species_map_fasta.exists():
        id_species_map = parse_species_map_fasta(path_to_species_map_fasta)
    elif args.species_map_fasta : # if path provided but file doesn't exist
         print(f"Warning: Species map FASTA file not found: {path_to_species_map_fasta}", file=sys.stderr if not args.output_file else original_stdout)

    species_probes, probes_sequences = parse_probes_fasta(path_to_probes_fasta, tax_data, name_to_tax_id, id_species_map)
    # if not probes_sequences: print(f"Warning: No probes loaded from FASTA file: {path_to_probes_fasta.resolve()}", file=sys.stderr)
    # else: print(f"FASTA probes parsed: {len(probes_sequences)} total, mapped to {len(species_probes)} species.", file=sys.stderr if not args.output_file else original_stdout)

    fungal_root_tax_ids = name_to_tax_id.get(args.fungi_root_name)
    processed_root_ids = set()

    if fungal_root_tax_ids:
        fungal_root_id_found = None
        for tid in fungal_root_tax_ids:
            if tax_data.get(tid, {}).get('rank') in ['kingdom', 'superkingdom']:
                fungal_root_id_found = tid
                break
        if not fungal_root_id_found and fungal_root_tax_ids:
            fungal_root_id_found = fungal_root_tax_ids[0]

        if fungal_root_id_found:
            print(f"\n--- Starting Decision Tree from '{args.fungi_root_name}' (TaxID: {fungal_root_id_found}) ---")
            display_decision_tree(fungal_root_id_found, tax_data, conservation_info, species_probes,
                                  probes_sequences=probes_sequences)
            processed_root_ids.add(fungal_root_id_found)
            print("\n" + "="*50 + "\n")
        # else:
            # print(f"Could not find a suitable TaxID for Fungal root: '{args.fungi_root_name}'. IDs found: {fungal_root_tax_ids}", file=sys.stderr)
    # else:
        # print(f"Fungal root name '{args.fungi_root_name}' not found in taxonomy.", file=sys.stderr)


    # Fallback or additional processing for top-level taxa from conservation data
    # print("\n--- Checking for other top-level taxa from conservation data (if not covered by Fungi root) ---", file=sys.stderr)
    potential_starts = []
    for (rank, tax_name), data in conservation_info.items():
        ids = name_to_tax_id.get(tax_name)
        if ids:
            for tax_id_val in ids:
                if tax_data.get(tax_id_val, {}).get('rank') == rank:
                    potential_starts.append({'id': tax_id_val, 'name': tax_name, 'rank': rank})
                    break
    unique_potential_starts_dict = {(item['id']): item for item in potential_starts}
    unique_potential_starts = list(unique_potential_starts_dict.values())
    rank_order = ['kingdom', 'subkingdom', 'phylum', 'subphylum', 'class', 'subclass', 'order', 'suborder', 'family', 'subfamily', 'genus', 'species', 'species group'] # Added 'species group'
    unique_potential_starts.sort(key=lambda x: rank_order.index(x['rank']) if x['rank'] in rank_order else len(rank_order))

    actual_start_nodes = []
    for item in unique_potential_starts:
        is_child_of_already_selected = False
        # Check against already processed_root_ids (like Fungi)
        for proc_id in processed_root_ids:
            if proc_id == item['id']: is_child_of_already_selected = True; break
            path_trace = [item['id']]
            tc = item['id']
            while tax_data.get(tc,{}).get('parent_id') != tc and tax_data.get(tc,{}).get('parent_id') in tax_data:
                tc = tax_data[tc]['parent_id']
                path_trace.append(tc)
                if tc == proc_id: is_child_of_already_selected = True; break
            if is_child_of_already_selected: break
        if is_child_of_already_selected: continue

        for existing_start_node in actual_start_nodes:
            path_trace = [item['id']]
            tc = item['id']
            is_descendant = False
            while tax_data.get(tc,{}).get('parent_id') != tc and tax_data.get(tc,{}).get('parent_id') in tax_data :
                tc = tax_data[tc]['parent_id']
                path_trace.append(tc)
                if tc == existing_start_node['id']: is_descendant = True; break
            if is_descendant: is_child_of_already_selected = True; break
        if is_child_of_already_selected: continue

        if item['id'] not in processed_root_ids:
            actual_start_nodes.append(item)

    actual_start_nodes.sort(key=lambda x: x['name'])

    # if not actual_start_nodes and not processed_root_ids:
        # print("No suitable top-level taxa found in conservation data to start the tree, and Fungi root also not processed.", file=sys.stderr)
    # elif not actual_start_nodes and processed_root_ids:
        # print("All relevant taxa from conservation files seem to be under the Fungi root (if processed).", file=sys.stderr)

    for start_node_data in actual_start_nodes:
        if start_node_data['id'] not in processed_root_ids:
            print(f"\n--- Displaying Tree from: {start_node_data['name']} ({start_node_data['rank']} - TaxID: {start_node_data['id']}) ---")
            display_decision_tree(start_node_data['id'], tax_data, conservation_info, species_probes,
                                  probes_sequences=probes_sequences)
            processed_root_ids.add(start_node_data['id'])
            print("\n" + "="*50 + "\n")

    print("Decision tree generation complete.", file=sys.stderr if not args.output_file else original_stdout)

    # Close output file and restore stdout
    if output_fh:
        output_fh.close()
        sys.stdout = original_stdout

if __name__ == "__main__":
    main()

