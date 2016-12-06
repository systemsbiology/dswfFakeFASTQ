|Build Status|

=====
**dswfFakeFASTQ**
=====

Generate a FASTQ file for testing the Duplex-Sequencing workflow

Creates two paired-end FASTQ files or one FASTQ file for testing the Duplex-Sequencing workflow.

Only requires a FASTA file of the read sequence to use to generate the FASTQ file, but has
many options that allow customization of the FASTQ file(s) created.

usage::
usage: makeFakeFASTQ.py [-h] [--barcode_length BARCODE_LENGTH]
                        [--barcode_a BARCODE_A] [--barcode_b BARCODE_B]
                        [--spacer_length SPACER_LENGTH]
                        [--max_num_families MAX_NUM_FAMILIES]
                        [--min_num_families MIN_NUM_FAMILIES]
                        [--num_families NUM_FAMILIES]
                        [--max_num_reads MAX_NUM_READS]
                        [--min_num_reads MIN_NUM_READS]
                        [--num_reads NUM_READS] [--wobble_reads WOBBLE_READS]
                        [--read_length READ_LENGTH]
                        [--max_num_flipped MAX_NUM_FLIPPED]
                        [--min_num_flipped MIN_NUM_FLIPPED]
                        [--num_flipped NUM_FLIPPED] [--prefix PREFIX]
                        [--instrument INSTRUMENT] [--flow_cell FLOW_CELL]
                        [--x_min X_MIN] [--x_max X_MAX] [--y_min Y_MIN]
                        [--y_max Y_MAX] [--lane_min LANE_MIN]
                        [--lane_max LANE_MAX]
                        [--quality_type {high,medium,low}]
                        [--swathes SWATHES [SWATHES ...]]
                        [--tile_min TILE_MIN] [--tile_max TILE_MAX]
                        [--paired_end PAIRED_END] [--is_filtered IS_FILTERED]
                        --fasta [FASTA]
                        [--include_fasta_header_in_fastq_header INCLUDE_FASTA_HEADER_IN_FASTQ_HEADER]
                        [--include_barcode_in_fastq_header INCLUDE_BARCODE_IN_FASTQ_HEADER]
                        [--map_file MAP_FILE]
                        [--frequencies_file [FREQUENCIES_FILE]]
                        [--buffer_end BUFFER_END]
                        [--buffer_both_sides BUFFER_BOTH_SIDES]
                        [--truncate_end TRUNCATE_END]
                        [--truncate_start TRUNCATE_START]
                        [--truncate_both_sides TRUNCATE_BOTH_SIDES]
                        [--truncate_by_read TRUNCATE_BY_READ]
                        [--buffer_seq BUFFER_SEQ] [--quality QUALITY]

optional arguments:
  -h, --help            show this help message and exit
  --barcode_length BARCODE_LENGTH, -bcl BARCODE_LENGTH
                        Length of Barcode at beginning and end of sequence.
                        Default: 10
  --barcode_a BARCODE_A, -fwbc BARCODE_A
                        Forward barcode string to use. Default: None
  --barcode_b BARCODE_B, -rvbc BARCODE_B
                        Reverse barcode string to use. Default: None
  --spacer_length SPACER_LENGTH, -sl SPACER_LENGTH
                        Length of Spacer between Barcode and sequence.
                        Default: 1
  --max_num_families MAX_NUM_FAMILIES, -maxnf MAX_NUM_FAMILIES
                        Maximum number of families per molecule.
  --min_num_families MIN_NUM_FAMILIES, -minnf MIN_NUM_FAMILIES
                        Minimum number of families per molecule. Default: 1
  --num_families NUM_FAMILIES, -nf NUM_FAMILIES
                        Number of families per molecule.
  --max_num_reads MAX_NUM_READS, -maxnr MAX_NUM_READS
                        Maximum number of reads per family. Default: 10
  --min_num_reads MIN_NUM_READS, -minnr MIN_NUM_READS
                        Minimum number of reads per family. Default: 1
  --num_reads NUM_READS, -nr NUM_READS
                        Number of reads per family.
  --wobble_reads WOBBLE_READS, -wr WOBBLE_READS
                        Wobble reads to create a distribution
  --read_length READ_LENGTH, -rl READ_LENGTH
                        Length of sequence in FASTQ output. Default: 156
  --max_num_flipped MAX_NUM_FLIPPED, -maxnff MAX_NUM_FLIPPED
                        Maximum number of flipped families per molecule.
  --min_num_flipped MIN_NUM_FLIPPED, -minnff MIN_NUM_FLIPPED
                        Minimum number of flipped familes per molecule.
                        Default: 1
  --num_flipped NUM_FLIPPED, -nff NUM_FLIPPED
                        Number of flipped families per molecule.
  --prefix PREFIX, -p PREFIX
                        Prefix for FASTQ output files. Default: DSWF
  --instrument INSTRUMENT, -inst INSTRUMENT
                        Instrument name in FASTQ reads. Default: NS500770
  --flow_cell FLOW_CELL, -fc FLOW_CELL
                        Flow cell name in FASTQ reads. Default: H5VNJAFXX
  --x_min X_MIN         Minimum X coordinate for FASTQ read. Default: 1015
  --x_max X_MAX         Maximum X coordinate for FASTQ read. Default: 26894
  --y_min Y_MIN         Minimum Y coordinate for FASTQ read. Default: 1017
  --y_max Y_MAX         Maximum Y coordinate for FASTQ read. Default: 20413
  --lane_min LANE_MIN   Minimum lane number for FASTQ read
  --lane_max LANE_MAX   Maximum lane number for FASTQ read
  --quality_type {high,medium,low}, -q {high,medium,low}
                        The quality of sequence: high, medium, low
  --swathes SWATHES [SWATHES ...]
                        The swathes on Illumina chip for FASTQ record Default:
                        [111, 112, 113, 114, 115, 116, 211, 212, 213, 214,
                        215, 216]
  --tile_min TILE_MIN   Minimum tile number for FASTQ read
  --tile_max TILE_MAX   Maximum tile number for FASTQ read
  --paired_end PAIRED_END
                        Produce paired end output. Default: 1
  --is_filtered IS_FILTERED
                        Produce filtered output. List. Default: [N]
  --fasta [FASTA], -f [FASTA]
                        A FASTA file to use as sequence for the reads
  --include_fasta_header_in_fastq_header INCLUDE_FASTA_HEADER_IN_FASTQ_HEADER
                        Include the FASTA header in the FASTQ file after the
                        control
  --include_barcode_in_fastq_header INCLUDE_BARCODE_IN_FASTQ_HEADER
                        Include the family random barcode in the FASTQ file
                        after the control and FASTA header if also selected.
  --map_file MAP_FILE   Create a map file of molecules to number of families
                        to number of reads.
  --frequencies_file [FREQUENCIES_FILE], -ff [FREQUENCIES_FILE]
                        File of frequencies for families
  --buffer_end BUFFER_END, -be BUFFER_END
                        Add buffer sequence to end of FASTA line. Default: 1
  --buffer_both_sides BUFFER_BOTH_SIDES, -bb BUFFER_BOTH_SIDES
                        Add buffer sequence to both sides of FASTA line.
                        Default: 0
  --truncate_end TRUNCATE_END, -te TRUNCATE_END
                        Truncate sequence at the end of the FASTA line.
                        Default: 1
  --truncate_start TRUNCATE_START, -ts TRUNCATE_START
                        Truncate sequence at the start of the FASTA line
                        Default: 1
  --truncate_both_sides TRUNCATE_BOTH_SIDES, -tbs TRUNCATE_BOTH_SIDES
                        Truncate both sides of FASTA sequence line. Default: 0
  --truncate_by_read TRUNCATE_BY_READ, -tbr TRUNCATE_BY_READ
                        Truncate paired end 1 reads at end, truncate paired
                        end 2 reads at start. Default: 0
  --buffer_seq BUFFER_SEQ, -buffSeq BUFFER_SEQ
                        Buffer string to use. Default: None
  --quality QUALITY, -qual QUALITY
                        Quality string to use. Default: None

DSWFFakeFASTQ produces a set of paired end FASTQ files with barcodes and spacers as if the 
FASTQ file had been produced by the DSWF procedure - amplification of sequence, attaching 
barcodes and spacers.  

The DSWF procedure samples a certain number of amplified molecules from the source DNA sample.
These are 'family' members.  The procedure then separates each double stranded molecule and
sequences each molecule multiple times.  These are 'reads'. Each 'family' has a different barcode.

DSWFFakeFASTQ takes an input FASTA file.  The entries in the FASTA file should be greater than
the read length that you want DSWFFakeFASTQ to produce.  Entries of 100 bp or less are sometimes
difficult to match to the genome uniquely using bwa.  It is recommended that you provide FASTA
entries of 300 bp or more. If the sequence in the FASTA file is shorter than the read length
desired, DSWFFakeFASTQ will pad the sequence with randomly generated sequence.

For each sequence in the FASTA file, DSWFFakeFASTQ randomly creates a number of molecules 
that will be sequenced as 'Num Families'.  If you have two sequences in the FASTA file that
contain the same sequence except for a SNP near the beginning of the sequence, seq1:C and seq1:T,
DSWFFakeFASTQ will create a random number of 'families' for each sequence.  It creates two reads
for each sequence, one in each alignment pattern.
As a hypothetical example, DSWFFakeFASTQ creates 2 families for seq1:C and 4 for seq1:T.  Each
family gets assigned a unique barcode.  Then DSWFFakeFASTQ will create a random number of reads
for each family.  If DSWFFakeFASTQ creates 5 reads for family 1 of seq1:C with a barcode of
AACAAGCAGT, then there will be 10 FASTQ entries for seq1:C with barcode AACAAGCAGT.  If it creates
3 reads for family 2 of seq1:C with a barcode of GCGGCACATG, then there will be 6 FASTQ entries
for seq1:C with a barcode of GCGGCACATG.  The numbers of families and reads with associated
barcodes are stored in a map_file.txt produced when DSWFFakeFASTQ is run. Depending on the
options selected, the FASTQ header will include the FASTA file header and/or the barcode
information for troubleshooting.

FASTA file:

>seq1:C
GTGATAGAGTGGCATTAGAAATTCCAGATAGAGCTAAAACTGAAGCTTTCCTTATAGAGATTTATCCTAGTTAGTTTGCGGGGATACTGGTTGGGCCGAAATCCTTTTGAAACTGGTTAAAACTCTCAGGGGCCCTTCCATTTGGTTTTCTGCAGCTGTGGATTCCCAACCAACAGTCATTGTGATCTTCCAAGCCAGAATGTGCTCTGGGCTGGAGTGGCAGCCCCTTATTCTGGCATTCAAGAGCGTGGGCACCCTTTGGCTATTTCTAGCATTTGTCTGGTTAGCCTTTGGGAAACG
>seq1:T
GTGATAGAGTGGCATTAGAAATTCCAGATAGAGCTAAAACTGAAGCTTTCCTTATAGAGATTTATCCTAGTTAGTTTGCGGGGATACTGGTTGGGCCGAAATCCTTTTGAAACTGGTTAAAACTCTCAGGGGCCCTTCCATTTGGTTTTCTGCAGCTGTGGATTCCCAACCAACAGTCATTGTGATCTTCCAAGCCAGAATGTGCTCTGGGCTGGAGTGGCAGCCCCTTATTCTGGCATTCAAGAGCGTGGGCACCCTTTGGCTATTTCTAGCATTTGTCTGGTTAGCCTTTGGGAAACG

Map file:

VERSION 0.06
FASTA Header    Num Familes     Num Reads       Num Flipped     Barcode A       Barcode B       Full Barcode
>seq1:C 2       4       2       AGAGGTCCCC      AATTTGCTAA      AGAGGTCCCCAATTTGCTAA
>seq1:C 2       4       2       GCCGCGCAGT      GAAATCCAAT      GCCGCGCAGTGAAATCCAAT
>seq1:T 4       8       4       TCTCGTTCCT      GGTAAATCAC      TCTCGTTCCTGGTAAATCAC
>seq1:T 4       7       3       CTGCAACTTA      AACTGTCGAA      CTGCAACTTAAACTGTCGAA
>seq1:T 4       4       2       TGAATAGATC      TACTGTAGTA      TGAATAGATCTACTGTAGTA
>seq1:T 4       4       2       ATTTACAGGG      ACCCATTTTG      ATTTACAGGGACCCATTTTG

The Num Families information is duplicated on every line and indicates the total number of lines 
of families for the FASTA sequence.  The number of families generated depends on the max_num_families
value.  Each line should have a unique Barcode and a number of reads generated depending on the 
max_num_reads value. The Num Reads value will be duplicated because duplex sequencing requires
ab and ba reads for each sequence to be valid. So the first barcode has 8 reads in seq1 and 8 in seq2

.. |Build Status| image:: https://travis-ci.org/systemsbiology/dswfFakeFASTQ.svg?branch=master
   :target: https://travis-ci.org/systemsbiology/dswfFakeFASTQ
