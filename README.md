Calculate coverage of assembled genome. \n

options: \n
  -h, --help            show this help message and exit \n
  -r REF_SIZE, --ref_size REF_SIZE \n 
                        Size of the assembly \n
  --illumina_dir ILLUMINA_DIR \n
                        Directory containing Illumina reads in FASTQ.gz format \n
  --nanopore_dir NANOPORE_DIR \n
                        Directory containing ONT reads in FASTQ.gz format \n
  --platform {illumina,nanopore,hybrid} \n 
                        Sequencing platform \n
                        
example usage:                      
HYBRID
python coverage_calculator.py -r 4411532 --platform hybrid --nanopore_dir ~/nanopore_reads --illumina_dir ~/illumina_reads

ILLUMINA
python coverage_calculator.py -r 4411532 --platform illumina --illumina_dir ~/illumina_reads

NANOPORE
python coverage_calculator.py -r 4411532 --platform nanopore --nanopore_dir ~/nanopore_reads

Output will be saved in your current directory
