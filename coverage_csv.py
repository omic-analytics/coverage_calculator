import argparse
import os
import gzip
import csv

def count_bases(directory):
    bases_by_sample = {}
    reads_by_sample = {}
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
                reads_by_sample[sample_name] = {'R1': 0, 'R2': 0}
            with gzip.open(file_path, 'rt') as file:
                for line_num, line in enumerate(file):
                    if line_num % 4 == 1:  # Sequence lines
                        if read_num == '1' or read_num == 'R1':
                            bases_by_sample[sample_name]['R1'] += len(line.strip())
                            reads_by_sample[sample_name]['R1'] += 1
                        elif read_num == '2' or read_num == 'R2':
                            bases_by_sample[sample_name]['R2'] += len(line.strip())
                            reads_by_sample[sample_name]['R2'] += 1
    return bases_by_sample, reads_by_sample

def count_bases_ont(directory):
    bases_by_sample = {}
    reads_by_sample = {}
    total_length_by_sample = {}
    for file_name in os.listdir(directory):
        if file_name.endswith('.fastq.gz'):
            file_path = os.path.join(directory, file_name)
            sample_name = file_name.split('.')[0]  # Extract sample name without file extension
            if sample_name not in bases_by_sample:
                bases_by_sample[sample_name] = 0
                reads_by_sample[sample_name] = 0
                total_length_by_sample[sample_name] = 0
            with gzip.open(file_path, 'rt') as file:
                for line_num, line in enumerate(file):
                    if line_num % 4 == 1:  # Sequence lines
                        read_length = len(line.strip())
                        bases_by_sample[sample_name] += read_length
                        reads_by_sample[sample_name] += 1
                        total_length_by_sample[sample_name] += read_length
    average_length_by_sample = {sample: total_length_by_sample[sample] / reads_by_sample[sample] 
                                for sample in reads_by_sample if reads_by_sample[sample] > 0}
    return bases_by_sample, reads_by_sample, average_length_by_sample

def calculate_coverage(ref_genome_size, illumina_bases, nanopore_bases):
    illumina_total_bases = sum(sum(sample.values()) for sample in illumina_bases.values())
    nanopore_total_bases = sum(nanopore_bases.values()) if nanopore_bases else 0  # Handling empty nanopore_bases dictionary
    total_bases = illumina_total_bases + nanopore_total_bases
    coverage = total_bases / ref_genome_size
    return coverage

def process_illumina_samples(args, illumina_bases, illumina_reads):
    fieldnames = ['Sample', 'Platform', 'R1 Bases', 'R2 Bases', 'R1 Reads', 'R2 Reads', 'Coverage']
    for sample, bases in illumina_bases.items():
        illumina_coverage = calculate_coverage(args.ref_size, {sample: bases}, {})
        total_r1_reads = illumina_reads[sample]['R1']
        total_r2_reads = illumina_reads[sample]['R2']
        print(f'Illumina sample: {sample}')
        print(f'Illumina R1 bases sequenced: {bases["R1"]}')
        print(f'Illumina R2 bases sequenced: {bases["R2"]}')
        print(f'Illumina R1 reads: {total_r1_reads}')
        print(f'Illumina R2 reads: {total_r2_reads}')
        print(f'Illumina coverage: {illumina_coverage:.2f}x')
        print()

        # Writing Illumina statistics to a CSV file
        csv_filename = f"{sample}_summary.csv"
        with open(csv_filename, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerow({
                'Sample': sample,
                'Platform': 'illumina',
                'R1 Bases': bases['R1'],
                'R2 Bases': bases['R2'],
                'R1 Reads': total_r1_reads,
                'R2 Reads': total_r2_reads,
                'Coverage': f'{illumina_coverage:.2f}'
            })

def process_nanopore_samples(args, nanopore_bases, nanopore_reads, nanopore_avg_length):
    fieldnames = ['Sample', 'Platform', 'Bases', 'Reads', 'Average Length', 'Coverage']
    for sample, bases in nanopore_bases.items():
        nanopore_coverage = calculate_coverage(args.ref_size, {}, {sample: bases})
        total_reads = nanopore_reads[sample]
        avg_length = nanopore_avg_length[sample]
        print(f'Nanopore sample: {sample}')
        print(f'Nanopore bases sequenced: {bases}')
        print(f'Nanopore reads: {total_reads}')
        print(f'Nanopore average read length: {avg_length:.2f}')
        print(f'Nanopore coverage: {nanopore_coverage:.2f}x')
        print()

        # Writing Nanopore statistics to a CSV file
        csv_filename = f"{sample}_summary.csv"
        with open(csv_filename, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerow({
                'Sample': sample,
                'Platform': 'nanopore',
                'Bases': bases,
                'Reads': total_reads,
                'Average Length': f'{avg_length:.2f}',
                'Coverage': f'{nanopore_coverage:.2f}'
            })

def process_hybrid_samples(args, illumina_bases, illumina_reads, nanopore_bases, nanopore_reads, nanopore_avg_length):
    fieldnames = ['Sample', 'Platform', 'Illumina R1 Bases', 'Illumina R2 Bases', 'Illumina R1 Reads', 'Illumina R2 Reads', 'Nanopore Bases', 'Nanopore Reads', 'Nanopore Average Length', 'Illumina Coverage', 'Nanopore Coverage', 'Hybrid Coverage']
    for sample in illumina_bases.keys() & nanopore_bases.keys():  # Get common sample names
        illumina_coverage = calculate_coverage(args.ref_size, {sample: illumina_bases[sample]}, {})
        nanopore_coverage = calculate_coverage(args.ref_size, {}, {sample: nanopore_bases[sample]})
        hybrid_coverage = illumina_coverage + nanopore_coverage
        total_r1_reads = illumina_reads[sample]['R1']
        total_r2_reads = illumina_reads[sample]['R2']
        total_nano_reads = nanopore_reads[sample]
        avg_length = nanopore_avg_length[sample]
        print(f'Hybrid sample: {sample}')
        print(f'Illumina R1 bases sequenced: {illumina_bases[sample]["R1"]}')
        print(f'Illumina R2 bases sequenced: {illumina_bases[sample]["R2"]}')
        print(f'Illumina R1 reads: {total_r1_reads}')
        print(f'Illumina R2 reads: {total_r2_reads}')
        print(f'Nanopore bases sequenced: {nanopore_bases[sample]}')
        print(f'Nanopore reads: {total_nano_reads}')
        print(f'Nanopore average read length: {avg_length:.2f}')
        print(f'Illumina coverage: {illumina_coverage:.2f}x')
        print(f'Nanopore coverage: {nanopore_coverage:.2f}x')
        print(f'Hybrid coverage: {hybrid_coverage:.2f}x')
        print()

        # Writing Hybrid statistics to a CSV file
        csv_filename = f"{sample}_summary.csv"
        with open(csv_filename, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerow({
                'Sample': sample,
                'Platform': 'hybrid',
                'Illumina R1 Bases': illumina_bases[sample]['R1'],
                'Illumina R2 Bases': illumina_bases[sample]['R2'],
                'Illumina R1 Reads': total_r1_reads,
                'Illumina R2 Reads': total_r2_reads,
                'Nanopore Bases': nanopore_bases[sample],
                'Nanopore Reads': total_nano_reads,
                'Nanopore Average Length': f'{avg_length:.2f}',
                'Illumina Coverage': f'{illumina_coverage:.2f}',
                'Nanopore Coverage': f'{nanopore_coverage:.2f}',
                'Hybrid Coverage': f'{hybrid_coverage:.2f}'
            })

def concatenate_summary_files():
    summary_files = [file for file in os.listdir() if file.endswith("_summary.csv")]
    if len(summary_files) >= 2:
        fieldnames = ['Sample', 'Platform', 'Illumina R1 Bases', 'Illumina R2 Bases', 'Illumina R1 Reads', 'Illumina R2 Reads', 'Nanopore Bases', 'Nanopore Reads', 'Nanopore Average Length', 'Illumina Coverage', 'Nanopore Coverage', 'Hybrid Coverage']
        with open("all_samples_summary.csv", "w", newline='') as output_file:
            writer = csv.DictWriter(output_file, fieldnames=fieldnames)
            writer.writeheader()
            for file in summary_files:
                with open(file, "r") as input_file:
                    reader = csv.DictReader(input_file)
                    for row in reader:
                        writer.writerow(row)
        print("All sample summaries have been concatenated into all_samples_summary.csv")
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
        illumina_bases, illumina_reads = count_bases(args.illumina_dir)
        process_illumina_samples(args, illumina_bases, illumina_reads)

    if args.platform == 'nanopore':
        if not args.nanopore_dir:
            parser.error('Please provide the directory containing ONT reads.')
        nanopore_bases, nanopore_reads, nanopore_avg_length = count_bases_ont(args.nanopore_dir)
        process_nanopore_samples(args, nanopore_bases, nanopore_reads, nanopore_avg_length)

    if args.platform == 'hybrid':
        if not args.illumina_dir or not args.nanopore_dir:
            parser.error('Please provide both Illumina and Nanopore directories for hybrid mode.')
        illumina_bases, illumina_reads = count_bases(args.illumina_dir)
        nanopore_bases, nanopore_reads, nanopore_avg_length = count_bases_ont(args.nanopore_dir)
        process_hybrid_samples(args, illumina_bases, illumina_reads, nanopore_bases, nanopore_reads, nanopore_avg_length)

    # Concatenate summary files
    concatenate_summary_files()

if __name__ == '__main__':
    main()
