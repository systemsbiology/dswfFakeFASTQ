|Build Status|

# dswfFakeFASTQ
Generate a FASTQ file for testing the Duplex-Sequencing workflow

Creates two paired-end FASTQ files or one FASTQ file for testing the Duplex-Sequencing workflow.

Only requires a FASTA file of the read sequence to use to generate the FASTQ file, but has
many options that allow customization of the FASTQ file(s) created.

```bash
usage: makeFakeFASTQ.py [-h] [--barcode_length BARCODE_LENGTH]
                        [--barcode BARCODE] [--spacer_length SPACER_LENGTH]
                        [--max_num_families MAX_NUM_FAMILIES]
                        [--max_num_reads MAX_NUM_READS]
                        [--read_length READ_LENGTH] [--prefix PREFIX]
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
                        [--map_file MAP_FILE] [--buffer_end BUFFER_END]
                        [--buffer_both_sides BUFFER_BOTH_SIDES]
                        [--truncate_end TRUNCATE_END]
                        [--truncate_both_sides TRUNCATE_BOTH_SIDES]
                        [--buffer_seq BUFFER_SEQ] [--quality QUALITY]
                        
optional arguments:
  -h, --help            show this help message and exit
  --barcode_length BARCODE_LENGTH, -bcl BARCODE_LENGTH
                        Length of Barcode at beginning and end of sequence.
                        Default: 10
  --barcode BARCODE, -bc BARCODE
                        Barcode string to use. Default: None
  --spacer_length SPACER_LENGTH, -sl SPACER_LENGTH
                        Length of Spacer between Barcode and sequence.
                        Default: 1
  --max_num_families MAX_NUM_FAMILIES, -nf MAX_NUM_FAMILIES
                        Maximum number of families per molecule. Default: 15
  --max_num_reads MAX_NUM_READS, -nr MAX_NUM_READS
                        Maximum number of reads per family. Default: 10
  --read_length READ_LENGTH, -rl READ_LENGTH
                        Length of sequence in FASTQ output. Default: 156
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
                        The swathes on this Illumina chip for FASTQ record.
                        Default: [111, 112, 113, 114, 115, 116, 211, 212, 213,
                        214, 215, 216]
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
                        after the control (and FASTA header if also selected).
  --map_file MAP_FILE   Create a map file of molecules to number of families
                        to number of reads.
  --buffer_end BUFFER_END, -be BUFFER_END
                        Add buffer sequence to end of FASTA line. Default: 1
  --buffer_both_sides BUFFER_BOTH_SIDES, -bbs BUFFER_BOTH_SIDES
                        Add buffer sequence to both sides of FASTA line.
                        Default: 0
  --truncate_end TRUNCATE_END, -te TRUNCATE_END
                        Truncate sequence at the end of the FASTA line.
                        Default: 1
  --truncate_both_sides TRUNCATE_BOTH_SIDES, -tbs TRUNCATE_BOTH_SIDES
                        Truncate both sides of FASTA sequence line. Default: 0
  --buffer_seq BUFFER_SEQ, -buffSeq BUFFER_SEQ
                        Buffer string to use. Default: None
  --quality QUALITY, -qual QUALITY
                        Quality string to use. Default: None
```

.. |Build Status| image:: https://travis-ci.org/systemsbiology/dswfFakeFASTQ.svg?branch=master
   :target: https://travis-ci.org/systemsbiology/dswfFakeFASTQ
