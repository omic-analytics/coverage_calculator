#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import argparse
import os
import gzip

def count_bases(directory):
    bases_by_sample = {}
    for file_name in os.listdir(directory):
        if file_name.endswith('.fastq.gz'):
            file_path = os.path.join(directory, file_name)
            if '_' in file_name:
                sample_name, read_num = file_name.rsplit('_', 1)
                read_num = read_num.split('.')[0]  # Remove the file extension
            else:
                sample_name, read_num = file_name.split('.')[0], '1'  # Default to R1 if not specified
            if sample_name not in bases_by_sample:
                bases_by_sample[sample_name] = {'R1': 0, 'R2': 0}
            with gzip.open(file_path, 'rt') as file:
                for line_num, line in enumerate(file):
                    if line_num % 4 == 1:  # Sequence lines
                        if read_num == '1' or read_num == 'R1':
                            bases_by_sample[sample_name]['R1'] += len(line.strip())
                        elif read_num == '2' or read_num == 'R2':
                            bases_by_sample[sample_name]['R2'] += len(line.strip())
    return bases_by_sample

def count_bases_ont(directory):
    bases_by_sample = {}
    for file_name in os.listdir(directory):
        if file_name.endswith('.fastq.gz'):
            file_path = os.path.join(directory, file_name)
            sample_name = file_name.split('.')[0]  # Extract sample name without file extension
            if sample_name not in bases_by_sample:
                bases_by_sample[sample_name] = 0
            with gzip.open(file_path, 'rt') as file:
                for line_num, line in enumerate(file):
                    if line_num % 4 == 1:  # Sequence lines
                        bases_by_sample[sample_name] += len(line.strip())
    return bases_by_sample

def calculate_coverage(ref_genome_size, illumina_bases, nanopore_bases):
    illumina_total_bases = sum(sum(sample.values()) for sample in illumina_bases.values())
    nanopore_total_bases = sum(nanopore_bases.values()) if nanopore_bases else 0  # Handling empty nanopore_bases dictionary
    total_bases = illumina_total_bases + nanopore_total_bases
    coverage = total_bases / ref_genome_size
    return coverage

def process_illumina_samples(args, illumina_bases):
    for sample, bases in illumina_bases.items():
        illumina_coverage = calculate_coverage(args.ref_size, {sample: bases}, {})
        print(f'Illumina sample: {sample}')
        print(f'Illumina R1 bases sequenced: {bases["R1"]}')
        print(f'Illumina R2 bases sequenced: {bases["R2"]}')
        print(f'Illumina coverage: {illumina_coverage:.2f}x')
        print()

        # Writing Illumina statistics to a log file
        log_filename = f"{sample}_summary.txt"
        with open(log_filename, 'w') as log_file:
            log_file.write(f'Sequence platform: illumina\n')
            log_file.write(f'Reference genome size: {args.ref_size}\n')
            log_file.write(f'Illumina sample: {sample}\n')
            log_file.write(f'Illumina R1 bases sequenced: {bases["R1"]}\n')
            log_file.write(f'Illumina R2 bases sequenced: {bases["R2"]}\n')
            log_file.write(f'Illumina coverage: {illumina_coverage:.2f}x\n\n')

def process_nanopore_samples(args, nanopore_bases):
    for sample, bases in nanopore_bases.items():
        nanopore_coverage = calculate_coverage(args.ref_size, {}, {sample: bases})
        print(f'Nanopore sample: {sample}')
        print(f'Nanopore bases sequenced: {bases}')
        print(f'Nanopore coverage: {nanopore_coverage:.2f}x')
        print()

        # Writing Nanopore statistics to a log file
        log_filename = f"{sample}_summary.txt"
        with open(log_filename, 'w') as log_file:
            log_file.write(f'Sequence platform: nanopore\n')
            log_file.write(f'Reference genome size: {args.ref_size}\n')
            log_file.write(f'Nanopore sample: {sample}\n')
            log_file.write(f'Nanopore bases sequenced: {bases}\n')
            log_file.write(f'Nanopore coverage: {nanopore_coverage:.2f}x\n\n')

def process_hybrid_samples(args, illumina_bases, nanopore_bases):
    for sample in illumina_bases.keys() & nanopore_bases.keys():  # Get common sample names
        illumina_coverage = calculate_coverage(args.ref_size, {sample: illumina_bases[sample]}, {})
        nanopore_coverage = calculate_coverage(args.ref_size, {}, {sample: nanopore_bases[sample]})
        hybrid_coverage = illumina_coverage + nanopore_coverage
        print(f'Hybrid sample: {sample}')
        print(f'Illumina coverage: {illumina_coverage:.2f}x')
        print(f'Nanopore coverage: {nanopore_coverage:.2f}x')
        print(f'Hybrid coverage: {hybrid_coverage:.2f}x')
        print()

        # Writing Hybrid statistics to a log file
        log_filename = f"{sample}_summary.txt"
        with open(log_filename, 'w') as log_file:
            log_file.write(f'Sequence platform: hybrid\n')
            log_file.write(f'Reference genome size: {args.ref_size}\n')
            log_file.write(f'Hybrid sample: {sample}\n')
            log_file.write(f'Illumina coverage: {illumina_coverage:.2f}x\n')
            log_file.write(f'Nanopore coverage: {nanopore_coverage:.2f}x\n')
            log_file.write(f'Hybrid coverage: {hybrid_coverage:.2f}x\n\n')
            
def concatenate_summary_files():
    summary_files = [file for file in os.listdir() if file.endswith("_summary.txt")]
    if len(summary_files) >= 2:
        with open("all_samples_summary.txt", "w") as output_file:
            for file in summary_files:
                sample_name = file.split("_summary.txt")[0]
                output_file.write(f"Sample: {sample_name}\n\n")
                with open(file, "r") as input_file:
                    output_file.write(input_file.read())
                output_file.write("\n\n\n")
        print("All sample summaries have been concatenated into all_samples_summary.txt")
    else:
        print("There are less than 2 summary files. No concatenation needed.")

def main():
    parser = argparse.ArgumentParser(description='Calculate coverage of a reference genome.')
    parser.add_argument('-r', '--ref_size', type=int, help='Size of the reference genome')
    parser.add_argument('--illumina_dir', help='Directory containing Illumina reads in FASTQ.gz format')
    parser.add_argument('--nanopore_dir', help='Directory containing ONT reads in FASTQ.gz format')
    parser.add_argument('--platform', choices=['illumina', 'nanopore', 'hybrid'], help='Sequencing platform')
    args = parser.parse_args()

    if not args.ref_size:
        parser.error('Please provide the size of the reference genome.')

    if args.platform == 'illumina':
        if not args.illumina_dir:
            parser.error('Please provide the directory containing Illumina reads.')
        illumina_bases = count_bases(args.illumina_dir)
        process_illumina_samples(args, illumina_bases)

    if args.platform == 'nanopore':
        if not args.nanopore_dir:
            parser.error('Please provide the directory containing ONT reads.')
        nanopore_bases = count_bases_ont(args.nanopore_dir)
        process_nanopore_samples(args, nanopore_bases)

    if args.platform == 'hybrid':
        if not args.illumina_dir or not args.nanopore_dir:
            parser.error('Please provide both Illumina and Nanopore directories for hybrid mode.')
        illumina_bases = count_bases(args.illumina_dir)
        nanopore_bases = count_bases_ont(args.nanopore_dir)
        process_hybrid_samples(args, illumina_bases, nanopore_bases)

    # Concatenate summary files
    concatenate_summary_files()

if __name__ == '__main__':
    main()

