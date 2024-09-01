#!/usr/bin/env python3

import re
import argparse

def extract_references(sam_file):
    references = {}

    with open(sam_file, 'r') as f:
        for line in f:
            if line.startswith('@SQ'):
                parts = line.strip().split('\t')
                ref_name = None
                ref_length = None
                for part in parts:
                    if part.startswith('SN:'):
                        ref_name = part.split(':')[1]
                    elif part.startswith('LN:'):
                        ref_length = int(part.split(':')[1])
                if ref_name and ref_length:
                    references[ref_name] = ref_length
            elif not line.startswith('@'):
                break  # Stop reading header once we reach the data section

    return references

def parse_cigar(cigar_string):
    operations = []
    current_op = ''
    for char in cigar_string:
        if char.isdigit():
            current_op += char
        else:
            operations.append((int(current_op), char))
            current_op = ''
    return operations

def get_last_n_characters(string, n):
    return string[-n:]

def parse_line(line, unique):
    columns = line.strip().split('\t')
    seq = columns[9]
    name = columns[0]
    ref = columns[2]
    start = int(columns[3])
    cigar = columns[5]
    AS = columns[11]
    XS = columns[12]

    if unique:
        top_score = None
        second_score = None

        match_as = re.search(r'AS:i:([\-\d]+)$', AS)
        if match_as:
            top_score = int(match_as.group(1))
        else:
            raise ValueError(f"ERROR: could not parse AS score value: {AS}")

        match_xs = re.search(r'XS:i:([\-\d]+)$', XS)
        if match_xs:
            second_score = int(match_xs.group(1))
            if top_score <= second_score:
                return None
        else:
            second_score = None  # If XS is not found, it might be None

    return {
        "seq": seq,
        "name": name,
        "ref": ref,
        "start": start,
        "cigar": cigar
    }

def calculate_sum_matches_gaps_deletions_insertions(sam_file, unique):
    references = extract_references(sam_file)
    results = []

    with open(sam_file, 'r') as f:
        for line in f:
            if line.startswith('@'):
                continue  # Skip header lines

            parsed_data = parse_line(line, unique)
            if not parsed_data:
                continue

            read_name = parsed_data["name"]
            ref_name = parsed_data["ref"]
            leftmost_position = parsed_data["start"]
            cigar_string = parsed_data["cigar"]
            operations = parse_cigar(cigar_string)
            cigar_sum = 0
            for length, op in operations:
                if op in ['M', 'N', 'D']:
                    cigar_sum += length

            # Add the sum of matches, gaps, and deletions to the leftmost position
            cigar_sum += leftmost_position - 1

            # Get tail
            tail = "other"
            seq_end = get_last_n_characters(parsed_data["seq"], 3)
            if cigar_sum == references[ref_name]:
                if seq_end == 'CCA':
                    tail = "CCA"
            elif cigar_sum == references[ref_name] - 1:
                if seq_end[1:] == 'CC':
                    tail = "CC"
            elif cigar_sum == references[ref_name] - 2:
                if seq_end[2:] == 'C':
                    tail = "C"
            elif cigar_sum < references[ref_name] - 7:
                tail = "tDR"

            results.append((read_name, ref_name, leftmost_position, cigar_sum, tail))

    return results

if __name__ == "__main__":
    import sys

    parser = argparse.ArgumentParser(description="Process SAM files and calculate read positions with optional unique alignment filtering.")
    parser.add_argument("input_sam_file", help="Input SAM file")
    parser.add_argument("output_tsv_file", help="Output TSV file")
    parser.add_argument("--unique", action="store_true", help="Filter for unique alignments where top_score > second_score")

    args = parser.parse_args()

    results = calculate_sum_matches_gaps_deletions_insertions(args.input_sam_file, args.unique)

    with open(args.output_tsv_file, 'w') as f:
        f.write("Read Name\tReference Name\tLeftmost Position\tRight Position\tTail\n")
        for result in results:
            f.write('\t'.join(map(str, result)) + '\n')
