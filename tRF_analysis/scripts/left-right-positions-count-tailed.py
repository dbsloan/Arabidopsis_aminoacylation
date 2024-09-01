#!/usr/bin/env python

import sys
import csv
from collections import defaultdict

def read_gene_lengths(file_path):
    """Reads the gene lengths file and returns a dictionary {gene_name: length}."""
    gene_lengths = {}
    with open(file_path, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            gene_name, length = row
            gene_lengths[gene_name] = int(length)
    return gene_lengths

def process_reads(file_path):
    """Processes the reads file and returns dictionaries for leftmost and rightmost positions grouped by tail types."""
    left_positions = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    right_positions = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    
    with open(file_path, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        next(reader)  # Skip header if present
        for row in reader:
            read_name, gene_name, left_pos, right_pos, tail = row
            left_pos = int(left_pos)
            right_pos = int(right_pos)
            tail_group = categorize_tail(tail)
            left_positions[gene_name][left_pos][tail_group] += 1
            right_positions[gene_name][right_pos][tail_group] += 1
            
    return left_positions, right_positions

def categorize_tail(tail):
    """Categorizes the tail into predefined groups."""
    if tail == 'C':
        return 'C'
    elif tail == 'CC':
        return 'CC'
    elif tail == 'CCA':
        return 'CCA'
    elif tail == 'tDR':
        return 'tDR'
    else:
        return 'other'

def write_output(file_path, gene_lengths, left_positions, right_positions):
    """Writes the output file with the specified format."""
    with open(file_path, 'w') as file:
        writer = csv.writer(file, delimiter='\t')
        header = ['GeneName', 'Position', 'Left_Other', 'Left_tDR', 'Left_C', 'Left_CC', 'Left_CCA', 'Right_Other','Right_tDR', 'Right_C', 'Right_CC', 'Right_CCA', 'Left_Total', 'Right_Total']
        writer.writerow(header)
        
        for gene_name, length in gene_lengths.items():
            for position in range(1, length + 1):
                left_counts = {
                    'other': left_positions[gene_name][position]['other'],
                    'tDR': left_positions[gene_name][position]['tDR'],
                    'C': left_positions[gene_name][position]['C'],
                    'CC': left_positions[gene_name][position]['CC'],
                    'CCA': left_positions[gene_name][position]['CCA'],
                }
                right_counts = {
                    'other': right_positions[gene_name][position]['other'],
                    'tDR': right_positions[gene_name][position]['tDR'],
                    'C': right_positions[gene_name][position]['C'],
                    'CC': right_positions[gene_name][position]['CC'],
                    'CCA': right_positions[gene_name][position]['CCA'],
                }
                left_total = sum(left_counts.values())
                right_total = sum(right_counts.values())
                row = [gene_name, position, 
                       left_counts['other'], left_counts['tDR'], left_counts['C'], left_counts['CC'], left_counts['CCA'], 
                       right_counts['other'], right_counts['tDR'], right_counts['C'], right_counts['CC'], right_counts['CCA'],
                       left_total, right_total]
                writer.writerow(row)

def main():
    if len(sys.argv) != 4:
        print("Usage: python script_name.py <gene_lengths_file> <reads_file> <output_file>")
        sys.exit(1)
    
    gene_lengths_file = sys.argv[1]
    reads_file = sys.argv[2]
    output_file = sys.argv[3]
    
    gene_lengths = read_gene_lengths(gene_lengths_file)
    left_positions, right_positions = process_reads(reads_file)
    write_output(output_file, gene_lengths, left_positions, right_positions)
    
    print(f"Output written to {output_file}")

if __name__ == "__main__":
    main()