Calculate coverage of assembled genome.

options:
  -h, --help            show this help message and exit
  -r REF_SIZE, --ref_size REF_SIZE
                        Size of the assembly
  --illumina_dir ILLUMINA_DIR
                        Directory containing Illumina reads in FASTQ.gz format
  --nanopore_dir NANOPORE_DIR
                        Directory containing ONT reads in FASTQ.gz format
  --platform {illumina,nanopore,hybrid}
                        Sequencing platform
                        
example usage:                      
HYBRID
python coverage_calculator.py -r 4411532 --platform hybrid --nanopore_dir ~/nanopore_reads --illumina_dir ~/illumina_reads

ILLUMINA
python coverage_calculator.py -r 4411532 --platform illumina --illumina_dir ~/illumina_reads

NANOPORE
python coverage_calculator.py -r 4411532 --platform nanopore --nanopore_dir ~/nanopore_reads

Output will be saved in your current directory
