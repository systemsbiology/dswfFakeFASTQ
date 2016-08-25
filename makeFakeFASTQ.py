#!/usr/bin/env python
import os
import sys
import time
import argparse
import gzip
from random import randint, sample
from Bio.Seq import Seq # used for reverse_complement()

VERSION = 0.05

# To generate duplex sequencing simulated data, we start with N number of
# molecules of double stranded DNA.  These are represented by the number of
# entries in the FASTA file that feeds this program (N). We then generate
# 1-X families (-nf), composed of 1-Y reads (-nr).

# The barcode on each side of a molecule can be read in either direction
# by the sequencer, so some amount of the families must have reversed
# barcodes (Z, num_flipped).

# This program assumes that the FASTA file contains FORWARD strand
# sequence and thus generates the reverse complement.

# ab:1 is a - FORWARD - b:1
# ab:2 is a - REVERSE - b:2

# ba:1 is b - REVERSE - a:1
# ba:2 is b - FORWARD - a:2

# We combine ab:1 and ba:2 to get the FORWARD sequence
# We combine ab:2 and ba:1 to get the REVERSE sequence
# We then get one duplex consensus for FORWARD and one
# for REVERSE.
# Note: This can be flipped where ab:1/ba:2 contain REVERSE
# and then ab:2/ba:2 contain FORWARD in real data

#  If num_flipped is set to 0 then only produces ab:1/ab:2

# Example:  Starting with 360 entries in the FASTA file, with X = 15 and
#           Y = 10, and Z = 5. The maximum number of entries in the FASTQ
#           file is 360*15*10 = 54,000 with half of reads flipped.

# Valid values for quality_type are 'high', 'medium', 'low'

# There are two ways to specify frequencies in the FASTQ.
# 1) Input FASTA contains a number of fasta entries for frequency
#    and num_families and num_reads are specified in options.
# 2) Input FASTA contains one entry and a frequency file is supplied

# DEFAULT:
# If provided a FASTQ without a frequencies file, makeFakeFASTQ will produce
# 1-max_num_families families with 1-max_num_reads reads for each FASTA entry.
# Ex: A fasta contains three entries for rs10092491 with allele 1 sequence
# and provides 10 entries for rs10092491 with allele 2 sequence. makeFakeFASTQ
# randomly choose to make 5 families for the first allele 1 seq, 6 for the
# second and 2 for the third for a total of 13 families for allele1.  It
# randomly chooses to make 13, 2, 5, 8, 4, 1, 9, 5, 9, 2 familes for each of
# the allele 2 sequences for a total of 58 families for allele2.  Thus
# frequencies are not preserved from the fasta file and may significantly
# vary depending on randomness.

# Specified with options:
# makeFakeFASTQ using the num_families option will create a specified number of
# families for each entry in the FASTA file.
# Ex: A fasta contains 3 entries for rs10092491 with allele 1 sequence and 10
# entries for rs100092491 with allele 2 sequence.  num_families is specified as
# 5, so makeFakeFASTQ will create 3*5 = 15 families for rs10092491 allele 1 and
# 10*5 = 50 families for rs10092491 allele 2. If num_reads isn't specified then
# it will create 1-max_num_reads to support each of these families, but if it
# is specified as 5 then each family would have 5 reads supporting it for a
# total of 15*5+50*5 = 325 reads in the FASTQ file.

# Frequency file:
# A frequency file is used to specify the number of reads in SSCS, DCS, and
# families. The frequency file is in the format
#  <fasta header>\t<num_families>\t<num flipped>\t<num_reads>\t<quality type>

# Ex:
# >rs100092491:allele1	3	5 	high
# >rs100092491:allele2	10	5	high

# This can also be used to combine both methods of frequency generation and
# quality as follows:
# >rs100092491:allele1	3	2   5	high
# >rs100092491:allele1	2	1   5	medium
# >rs100092491:allele2	10	3   5	high
# >rs100092491:allele2	5	3   5	low
# This file would produce 3 families of 5 reads (with 3 ab and 2 ba) for
# allele1 line 1 and 2 families of 5 reads (with 1 ab and 1 ba) for allele1
# line 2.  As these FASTA sequences are, this would produce the same counts
# as having 5 families with 5 reads, but the quality distribution differs.

OPT_DEFAULTS = {
    'barcode_length': 10, 'spacer_length': 1,
    'min_num_families': 2, 'max_num_families': 15,
    'min_num_reads': 3, 'max_num_reads': 20,
    'min_num_flipped': 2, 'max_num_flipped': 10, # defaults to half of reads
    'read_length': 156, 'buffer_both_sides': 0,
    'buffer_end': 1, 'truncate_by_read': 1, 'rand_window': None,
    'truncate_both_sides': 0, 'truncate_end': 0, 'truncate_start': 0,
    'prefix': 'DSWF', 'instrument': 'NS500770', 'flow_cell': 'H5VNJAFXX',
    'x_min': 1015, 'x_max': 26894, 'y_min': 1017, 'y_max': 20413,
    'lane_min': 1, 'lane_max': 4, 'quality_type': 'high',
    'swathes': [111, 112, 113, 114, 115, 116, 211, 212, 213, 214, 215, 216],
    'tile_min': 1, 'tile_max': 12, 'paired_end': 1, 'is_filtered': ['N'],
    'include_fasta_header_in_fastq_header': 1,
    'include_barcode_in_fastq_header': 1,
    'map_file': 1, 'tag_file': 1  # DEBUG
}
DESCRIPTION = """Create FASTQ data."""

QUAL_SET = {'high': [35, 41], 'medium': [23, 30], 'low': [3, 15]}


# Design:
#  Multiple Duplex Consensus Sequences (DCS) form a clan
#  Multiple 5->3 and 3->5 Single Strand Consensus sequences form a family
#  An individual sequence is a read
#  There are always paired end reads, so one 5->3 and one 3->5
#  make_clan makes multiple families via calls to make_family
#  make_family makes multiple paired end reads
#  the 5->3 read for the top strand is called top_five_for_ds_read
#  the 3->5 read for the top strand is called top_three_for_ds_read
#  the 5->3 read for the bottom strand is called bottom_five_for_ds_read
#  the 3->5 read for the bottom strand is called bottom_three_for_ds_read

def main(argv):
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.set_defaults(**OPT_DEFAULTS)

    parser.add_argument('--barcode_length', '-bcl', type=int, required=False,
                        help='Length of Barcode at beginning and end of\
                        sequence. Default: 10')
    parser.add_argument('--fwbarcode', '-fwbc', type=int, required=False,
                        help='Forward barcode string to use. Default: None')
    parser.add_argument('--rvbarcode', '-rvbc', type=int, required=False,
                        help='Reverse barcode string to use. Default: None')
    parser.add_argument('--spacer_length', '-sl', type=int, required=False,
                        help='Length of Spacer between Barcode and sequence.\
                        Default: 1')
    parser.add_argument('--max_num_families', '-maxnf', type=int,
                        required=False, help='Maximum number of families per \
                        molecule.')
    parser.add_argument('--min_num_families', '-minnf', type=int,
                        required=False, help='Minimum number of families per \
                        molecule. Default: 1')
    parser.add_argument('--num_families', '-nf', type=int,
                        required=False, help='Number of families per \
                        molecule.')
    parser.add_argument('--max_num_reads', '-maxnr', type=int, required=False,
                        help='Maximum number of reads per family.\
                        Default: 10')
    parser.add_argument('--min_num_reads', '-minnr', type=int, required=False,
                        help='Minimum number of reads per family.\
                        Default: 1')
    parser.add_argument('--num_reads', '-nr', type=int, required=False,
                        help='Number of reads per family.')
    parser.add_argument('--read_length', '-rl', type=int, required=False,
                        help='Length of sequence in FASTQ output.\
                        Default: 156')
    parser.add_argument('--max_num_flipped', '-maxnff', type=int,
                        required=False, help='Maximum number of flipped \
                        families per molecule.')
    parser.add_argument('--min_num_flipped', '-minnff', type=int,
                        required=False, help='Minimum number of flipped \
                        familes per molecule. Default: 1')
    parser.add_argument('--num_flipped', '-nff', type=int,
                        required=False, help='Number of flipped families \
                        per molecule.')
    parser.add_argument('--prefix', '-p', type=str, required=False,
                        help='Prefix for FASTQ output files. Default: DSWF')
    parser.add_argument('--instrument', '-inst', type=str, required=False,
                        help='Instrument name in FASTQ reads.\
                        Default: NS500770')
    parser.add_argument('--flow_cell', '-fc', type=str, required=False,
                        help='Flow cell name in FASTQ reads.\
                        Default: H5VNJAFXX')
    parser.add_argument('--x_min', type=int, required=False,
                        help='Minimum X coordinate for FASTQ read.\
                        Default: 1015')
    parser.add_argument('--x_max', type=int, required=False,
                        help='Maximum X coordinate for FASTQ read.\
                        Default: 26894')
    parser.add_argument('--y_min', type=int, required=False,
                        help='Minimum Y coordinate for FASTQ read.\
                        Default: 1017')
    parser.add_argument('--y_max', type=int, required=False,
                        help='Maximum Y coordinate for FASTQ read.\
                        Default: 20413')
    parser.add_argument('--lane_min', type=int, required=False,
                        help='Minimum lane number for FASTQ read')
    parser.add_argument('--lane_max', type=int, required=False,
                        help='Maximum lane number for FASTQ read')
    parser.add_argument('--quality_type', '-q', type=str, required=False,
                        help="The quality of sequence: high, medium, low",
                        choices=["high", "medium", "low"])
    parser.add_argument('--swathes', nargs='+', type=int,
                        help='The swathes on Illumina chip for FASTQ record\
                        Default: [111, 112, 113, 114, 115, 116, 211, 212,\
                        213, 214, 215, 216]')
    parser.add_argument('--tile_min', type=int, required=False,
                        help='Minimum tile number for FASTQ read')
    parser.add_argument('--tile_max', type=int, required=False,
                        help='Maximum tile number for FASTQ read')
    parser.add_argument('--paired_end', type=int, required=False,
                        help='Produce paired end output.  Default: 1')
    parser.add_argument('--is_filtered', type=int, required=False,
                        help='Produce filtered output. List.  Default: [N]')
    parser.add_argument('--fasta', '-f', nargs='?', required=True,
                        help='A FASTA file to use as sequence for the reads')
    parser.add_argument('--include_fasta_header_in_fastq_header', type=int,
                        required=False, help='Include the FASTA header in the\
                        FASTQ file after the control ')
    parser.add_argument('--include_barcode_in_fastq_header', type=int,
                        required=False, help='Include the family random barcode\
                        in the FASTQ file after the control and FASTA header\
                        if also selected.')
    parser.add_argument('--map_file', type=int, required=False,
                        help='Create a map file of molecules to number of\
                        families to number of reads.')
    parser.add_argument('--frequencies_file', '-ff', nargs='?', required=False,
                        help='File of frequencies for families')

    # Sequence Buffering and Truncation Options
    parser.add_argument('--buffer_end', '-be', type=int, required=False,
                        help='Add buffer sequence to end of FASTA line.\
                        Default: 1')
    parser.add_argument('--buffer_both_sides', '-bb', type=int, required=False,
                        help='Add buffer sequence to both sides of FASTA line.\
                        Default: 0')
    parser.add_argument('--truncate_end', '-te', type=int, required=False,
                        help='Truncate sequence at the end of the FASTA line.\
                        Default: 1')
    parser.add_argument('--truncate_start', '-ts', type=int, required=False,
                        help='Truncate sequence at the start of the FASTA line\
                        Default: 1')
    parser.add_argument('--truncate_both_sides', '-tbs', type=int,
                        required=False, help='Truncate both sides of FASTA\
                        sequence line. Default: 0')
    parser.add_argument('--truncate_by_read', '-tbr', type=int, required=False,
                        help="Truncate paired end 1 reads at end, truncate\
                        paired end 2 reads at start. Default: 0")

    # Testing options
    parser.add_argument('--buffer_seq', '-buffSeq', type=int, required=False,
                        help='Buffer string to use. Default: None')
    parser.add_argument('--quality', '-qual', type=int, required=False,
                        help='Quality string to use. Default: None')

    args = parser.parse_args(argv[1:])
    # check that max_num_reads >= min_num_reads
    if args.min_num_reads > args.max_num_reads:
        args.max_num_reads = args.min_num_reads
    if args.min_num_families > args.max_num_families:
        args.max_num_families = args.min_num_families
    if args.min_num_flipped > args.max_num_flipped:
        args.max_num_flipped = args.min_num_flipped
    print("args type is {0}".format(type(args)))

    fasta = open(args.fasta)
    seq1_file = gzip.open(args.prefix + '_seq1.fastq.gz', 'wb')
    seq2_file = gzip.open(args.prefix + '_seq2.fastq.gz', 'wb')
    if args.map_file is 1:
        args.map_file = gzip.open(args.prefix + '_map.txt.gz', 'wb')
        args.map_file.write("VERSION\t"+str(VERSION)+"\n")
        args.map_file.write("\t".join(["FASTA Header", "Num Familes",
                                       "Num Reads", "Num Flipped", "Barcode"])
                                       + "\n")
    if args.tag_file is 1:
        args.tag_file = gzip.open(args.prefix + '_tags.txt.gz', 'wb')
        args.tag_file.write('VERSION\t'+str(VERSION)+"\n")
        args.tag_file.write("\t".join(["FASTA Header", "Barcode","Reads"]) 
                                      + "\n")

    print('opened file ' + args.fasta)
    while True:
        header = fasta.readline().rstrip('\r\n')
        line = fasta.readline()
        if not line:
            break
        seq = line.rstrip('\r\n').upper()
        (clan_seq1, clan_seq2) = make_clan(header, seq, args)
        seq1_file.write("\n".join(map("\n".join,clan_seq1)) + "\n")
        seq2_file.write("\n".join(map("\n".join,clan_seq2)) + "\n")
    if args.map_file:
        args.map_file.close()
    seq1_file.close()
    seq2_file.close()
    fasta.close()
    print("Finished generating FASTQ files")

    exit(0)

# make_clan creates a number of families per clan


def make_clan(header, seq, args):
    if args.min_num_families > args.max_num_families:
        raise TypeError("Incorrect value of min_num_families or max_num_families")
    if args.num_families == None:
        num_families = randint(args.min_num_families, args.max_num_families)
    print("making {0} families for {1}".format(num_families, header))
    clan_seq = []
    clan_seq2 = []
    for f in range(1, num_families + 1):
        family_seq1, family_seq2 = make_family(header, seq, args, num_families)
        clan_seq.append(family_seq1)
        clan_seq2.append(family_seq2)
    return (clan_seq, clan_seq2)

# add random sequence to the FASTA line to reach read length
# count is used to buffer a sequence that needs a barcode+spacer added on the
# front but since we may be buffer_both_sides, we can't include the
# barcode+spacer on front until the read is generated


def buffer_sequence(args, seq, count=None):
    seq_diff = count if count else args.read_length - len(seq)
    if seq_diff <= 0:
        return seq
    buffer_seq = args.buffer_seq if args.buffer_seq else random_sequence(
        seq_diff)
    if args.buffer_end:
        new_seq = ''.join(str(v) for v in [seq, buffer_seq[:seq_diff]])
    else:
        if args.buffer_both_sides:
            cnt_both_sides = int(seq_diff / 2)
            cnt_front = int(seq_diff % 2)
            end_buffer_seq = args.buffer_seq[
                :cnt_both_sides] or random_sequence(cnt_both_sides)
            front_buffer_seq = args.buffer_seq[
                :cnt_both_sides] or random_sequence(cnt_front)
            if cnt_front:
                joined_seq = ''.join(
                    [args.buffer_seq[:cnt_front], front_buffer_seq])
                front_buffer_seq = joined_seq
            new_seq = ''.join([front_buffer_seq, seq, end_buffer_seq])
    return new_seq

# remove sequence from the FASTA line to reach read length depending on config
# count is used to truncate a sequence that needs a barcode+spacer added on the
# front but since we may be truncate_both_sides, we can't include the
# barcode+spacer on front until the read is generated


def truncate_sequence(args, seq, count=None, read_type=None):
    seq_diff = count if count else len(seq) - args.read_length
    if seq_diff <= 0:
        return seq
    ts = args.truncate_start
    te = args.truncate_end
    ws = seq_diff  # window start
    we = len(seq)  # window end
    if args.truncate_by_read:
        if read_type == '5':
            ts = 0
            te = 1
        elif read_type == '3':
            ts = 1
            te = 0
            # taking the wobble out so that UnifiedConsensusMaker works 20160825
            # if we're making a 3' read we want to move the window randomly a
            # little bit to the 5' to make bwa happy with the read distribution
            #rand_num = randint(0, 15)
            #if args.rand_window is not None:
            #    rand_num = args.rand_window
            #print("rand num is {}".format(rand_num))
            #ws = ws - rand_num
            #we = we - rand_num
        else:
            raise TypeError("read_type not set in truncate_sequence")

    print("ts {} te {} ws {} we {} args.rand_window {} args.truncate_both_sides {}".format(ts, te, ws, we, args.rand_window, args.truncate_both_sides))
    if te:
        new_seq = seq[0:-seq_diff]
    elif ts:
        new_seq = seq[ws:we]
    elif args.truncate_both_sides:
            cnt_both_sides = int(seq_diff / 2)
            cnt_front = int(seq_diff % 2)
            new_seq = seq[cnt_both_sides:-cnt_both_sides]
            if cnt_front:
                new_seq = seq[cnt_both_sides + cnt_front:-cnt_both_sides]
    else:
        raise TypeError("no truncate type specified")
    return new_seq


def make_ds_read(args, seq, barcode, read_type=None):
    ds_spacer = 'T' * args.spacer_length
    ds_length = len(ds_spacer) + len(barcode)
    total_length = len(seq) + ds_length
    if total_length == args.read_length:
        return "{0}{1}{2}".format(barcode, ds_spacer, seq)
    if (total_length > args.read_length):
        ds_seq = truncate_sequence(args, seq, total_length - args.read_length,
                                   read_type)
    if (total_length < args.read_length):
        ds_seq = buffer_sequence(args, seq, args.read_length - total_length)
    read = "{0}{1}{2}".format(barcode, ds_spacer, ds_seq)
    return read

# FASTQ is:
# Line1 @sequence id
# Line2 raw sequence letters
# Line3 +sequence id or just +
# Line4 encoded quality values for sequence
# @NS500773:24:H5VNJAFXX:1:11101:13101:1024 1:N:0:CGATGT
# GGAAANGTGCTGCCCATTACAAAATTAAATCAAACTCAACCTACCACTCACCTGAAATGCCTATGGTTCAAAGTAATAAGATTCATGAGACTTCTCTAAAAGTGGATTATTATGATCAGAAAGAATATGATCCACATTGTATGGTTTTTAGGCACC
# +
# AAAAAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE<EEEEEEEEEEEAEEEEEEEEEAEEEEEEEEEEEEEEEEEEEE<EEEEEEEEEEEEEEEEEAEEEEE<EEEEEEEEEEEEEEEEEAE<<AA


def make_family(header, seq, args, num_families):
    print("START make_family header {} seq {} args {}".format(header, seq, args))
    if args.fwbarcode:
        fwbarcode = args.fwbarcode
    else:
        fwbarcode = random_sequence(args.barcode_length)
    if args.rvbarcode:
        rvbarcode = args.rvbarcode
    else:
        rvbarcode = random_sequence(args.barcode_length)

    while rvbarcode == fwbarcode:
        rvbarcode = random_sequence(args.barcode_length)

    if rvbarcode > fwbarcode:
        hold = rvbarcode
        rvbarcode = fwbarcode
        fwbarcode = hold

    fullbarcode = fwbarcode + rvbarcode
    # set up number of reads
    if args.min_num_reads > args.max_num_reads:
        raise TypeError("Incorrect value of min_num_reads or max_num_reads")
    num_reads = randint(args.min_num_reads, args.max_num_reads)
    if args.num_reads:
        num_reads = args.num_reads

    # set up number of flipped
    if args.min_num_flipped > args.max_num_flipped:
        raise TypeError("Incorrect value of min_num_flipped or max_num_flipped")
    num_flipped = randint(args.min_num_flipped, args.max_num_flipped)
    # must also be half or less of the num_reads
    if num_flipped > num_reads/2:
        print("changing num_flipped from {} to {}".format(num_flipped, int(num_reads/2)))
        num_flipped = int(num_reads/2)
    if args.num_flipped:
        num_flipped = args.num_flipped


    if args.map_file:
        args.map_file.write("{0}	{1}	{2}	{3}	{4}	{5}\n".format(
            header, num_families, num_reads, num_flipped, fwbarcode,
            rvbarcode))
    print("making {0} reads for fwbarcode {1} rvbarcode {2} header {3}"
     .format(num_reads, fwbarcode, rvbarcode, header))
    print("num_flipped {}".format(num_flipped))
    read_set = []
    read_set2 = []
    read_names_top = []
    read_names_bottom = []
    #  the 5->3 read for the top strand is called top_five_for_ds_read
    #  the 3->5 read for the top strand is called top_three_for_ds_read
    top_five_for_ds_read = make_ds_read(args, seq, fwbarcode, '5')
    top_three_for_ds_read = make_ds_read(args, seq, fwbarcode, '3')
    top_five_rev_ds_read = make_ds_read(args, seq, rvbarcode, '5')
    top_three_rev_ds_read = make_ds_read(args, seq, rvbarcode, '3')

    top_five_quality = args.quality if args.quality else fastq_quality(
        args, len(top_five_for_ds_read))
    top_three_quality = args.quality if args.quality else fastq_quality(
        args, len(top_three_for_ds_read))
    #  the 5->3 read for the bottom strand is called bottom_three_for_ds_read
    #  the 3->5 read for the bottom strand is called bottom_five_for_ds_read
    bSeq = Seq(seq)
    rev_seq = bSeq.reverse_complement()
    print("rev_seq {}".format(rev_seq))
    print("making bottom_for_ds_read")
    bottom_three_for_ds_read = make_ds_read(args, seq, fwbarcode, '3')
    bottom_five_for_ds_read = make_ds_read(args, seq, fwbarcode, '5')
    print("bottom_three_for_ds_read {} bottom_five_for_ds_read {}"
            .format(bottom_three_for_ds_read,bottom_five_for_ds_read))

    print("making bottom_rev_for_ds_read")
    bottom_five_rev_ds_read = make_ds_read(args, seq, rvbarcode, '5')
    bottom_three_rev_ds_read = make_ds_read(args, seq, rvbarcode, '3')
    print("bottom_three_rev_ds_read {} bottom_five_rev_ds_read {}"
            .format(bottom_three_rev_ds_read,bottom_five_rev_ds_read))
    bottom_three_quality = args.quality if args.quality else fastq_quality(
        args, len(bottom_three_for_ds_read))
    bottom_five_quality = args.quality if args.quality else fastq_quality(
        args, len(bottom_five_for_ds_read))

    # we have num_reads to produce the following:
    count_flipped = 0
    for i in range(1, num_reads + 1):
        # add a reversed read if we have flipped to add
        if (num_flipped > count_flipped):
            count_flipped =+ 1
            (bottom_fastq_header, bottom_paired_header) = fastq_entry_header(args,
                                                            header, fullbarcode)
            (bottom_fastq_header2, bottom_paired_header2) = fastq_entry_header(args,
                                                            header, fullbarcode)
            print("bot fastq_head {}, bot paired_head {}, head {} fullbarcode {}"
                    .format(bottom_fastq_header, bottom_paired_header, header, fullbarcode))
            print("bot fastq_head2 {}, bot paired_head2 {}, head {} fullbarcode {}"
                    .format(bottom_fastq_header2, bottom_paired_header2, header, fullbarcode))
            # in order to have the proper barcodes, we need family 1 to have for and family2 to have rev
            # read set 1 needs fastq_header and bottom_three_for_ds_read
            # read set 2 needs fastq_header and bottom_five_rev_ds_read
            # read set 1 needs fastq_header2 and bottom_five_for_ds_read
            # read set 2 needs fastq_header2 and bottom_three_rev_ds_read
            read_names_bottom.append(bottom_fastq_header)
            read_names_bottom.append(bottom_fastq_header2)
            read_names_bottom.append(bottom_paired_header)
            read_names_bottom.append(bottom_paired_header2)

            # each read set needs to have a bottom_three and a bottom_five
            # set up the ab reads
            read_set.append(bottom_fastq_header)
            read_set.append(bottom_three_for_ds_read)
            read_set.append("+")
            read_set.append(bottom_three_quality)

            read_set2.append(bottom_paired_header)
            read_set2.append(bottom_five_rev_ds_read)
            read_set2.append("+")
            read_set2.append(bottom_five_quality)
            print("forward:\n")
            print("{}\n{}".format(bottom_fastq_header, bottom_three_for_ds_read))
            print("{}\n{}".format(bottom_paired_header, bottom_five_for_ds_read))
            print("reverse:\n")
            print("{}\n{}".format(bottom_fastq_header2, bottom_five_rev_ds_read))
            print("{}\n{}".format(bottom_paired_header2, bottom_three_rev_ds_read))

            # set up ba reversed
            read_set.append(bottom_fastq_header2)
            read_set.append(bottom_three_rev_ds_read)
            read_set.append("+")
            read_set.append(bottom_three_quality)

            read_set2.append(bottom_paired_header2)
            read_set2.append(bottom_five_for_ds_read)
            read_set2.append("+")
            read_set2.append(bottom_five_quality)
        else:
            (fastq_header, paired_header) = fastq_entry_header(args, header,
                                                               fullbarcode)
            (fastq_header2, paired_header2) = fastq_entry_header(args, header,
                                                               fullbarcode)
            print("fastq_header {}, paired_header {}, header {} fullbarcode {}"
                    .format(fastq_header, paired_header, header, fullbarcode))
            read_names_top.append(fastq_header)
            read_names_top.append(paired_header)
            read_names_top.append(fastq_header2)
            read_names_top.append(paired_header2)

            # read set 1 needs fastq_header2 and top_five_for_ds_read
            # read set 2 needs fastq_header2 and top_three_rev_ds_read
            # read set 1 needs fastq_header and top_three_for_ds_read
            # read set 2 needs fastq_header and top_five_rev_ds_read
            # each read set needs to have a top_three and a top_five
            # set up the ab reads
            read_set.append(fastq_header)
            read_set.append(top_five_for_ds_read)
            read_set.append("+")
            read_set.append(top_five_quality)

            read_set2.append(paired_header)
            read_set2.append(top_three_rev_ds_read)
            read_set2.append("+")
            read_set2.append(top_three_quality)

            # set up the ba reads
            read_set.append(fastq_header2)
            read_set.append(top_five_rev_ds_read)
            read_set.append("+")
            read_set.append(top_five_quality)

            read_set2.append(paired_header2)
            read_set2.append(top_three_for_ds_read)
            read_set2.append("+")
            read_set2.append(top_three_quality)

    print("fullbarcode {} read_names_top {} read_names_bottom {}".format(fullbarcode, read_names_top, read_names_bottom))
    args.tag_file.write("\t".join([fullbarcode, ",".join(read_names_top), ",".join(read_names_bottom)]) + "\n")
    return read_set, read_set2


# quality scores are PHRED scores from 0-41 converted into a character +33
def fastq_quality(args, seq_len):
    qualities = []
    for i in range(1, seq_len + 1):
        qualities.append(randint(QUAL_SET[args.quality_type][
                         0], QUAL_SET[args.quality_type][1]))
    # avg_qual = sum(qualities) / len(qualities)
    # print("avg_qual ",avg_qual)
    qual_string = map(lambda x: chr(x + 33), qualities)
    return(''.join(qual_string))

# base composition is an array of each nucleotide
# length is the length of the random sequence to generate


def random_sequence(length):
    bases = ['A', 'G', 'T', 'C']
    random_sequence = ['0'] * length
    for i in range(length):
        random_sequence[i] = bases[randint(0, 3)]
    return ''.join(random_sequence)

# make a sequence id header
# @<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos> <read>:<is filtered>:<control number>:<sample number> # nopep8
# <instrument>     Characters allowed: a-z, A-Z, 0-9, underscore   Instrument ID # nopep8
# <run number>     Numerical                                       Run number on instrument # nopep8
# <flowcell ID>    Characters allowed: a-z, A-Z, 0-9               Flow cell ID. # nopep8
# <lane>           Numerical                                       Lane number # nopep8
# <tile>           Numerical                                       Tile number # nopep8
# <x_pos>          Numerical                                       X coordinate of cluster # nopep8
# <y_pos>          Numerical                                       Y coordinate of cluster # nopep8
# <read>           Numerical                                       Read number. 1 can be single read or Read 2 of paired-end # nopep8
# <is filtered>    Y or N                                          Y if the read is filtered (did not pass), N otherwise # nopep8
# <control number> Numerical                                       0 when none of the control bits are on, otherwise it is an even number. # nopep8
#                                                                  On HiSeq X and NextSeq systems, control specification is not performed # nopep8
#                                                                  and this number is always 0. # nopep8
# <sample number>  Numerical                                       Sample number from sample sheet # nopep8
# @NS500773:24:H5VNJAFXX:1:11101:13101:1024 1:N:0:CGATGT           Paired end #1 # nopep8
# @NS500773:24:H5VNJAFXX:1:11101:13101:1024 2:N:0:CGATGT           Paired end #2 # nopep8
#  there is one instrument per FASTQ file.
#  there is one run number per FASTQ file.
#  there is one flowcell ID per FASTQ file.
#  there are 4 lanes per FASTQ file (1-4)
#  there are 144 tiles per FASTQ file (1-21612 in increments of 12)
#     first 3 positions are swaths: 111, 112, 113, 114, 115, 116. 211, 212,
#                                   213, 214, 215, 216
#     last two are first three combined with 01-12 : 21601, 21602, 21603,
#                                   21604, 21605, 21606, 21607, 21608, 21609,
#                                   21610, 21611, 21612).  Called tiles
#  there are 25880 x_pos.  From 1015 to 26894.
#  there are 19397 y_pos.  From 1017 to 20413.
#  there is one read number per FASTQ file: 2 for paired end FASTQ files
#                                           or 1 otherwise
#  is filtered is either Y or N - there are no filtered reads in CODIS
#  control number is 0


def fastq_entry_header(args, fasta_header, barcode):
    x_pos = randint(args.x_min, args.x_max)
    y_pos = randint(args.y_min, args.y_max)
    lane = randint(args.lane_min, args.lane_max)
    swath = sample(args.swathes, 1)[0]
    tile = randint(args.tile_min, args.tile_max)
    full_tile = "{0}{1:02}".format(swath, tile)
    filter = sample(args.is_filtered, 1)[0]
    if args.include_fasta_header_in_fastq_header:
        if fasta_header.startswith(">"):
            c = '0' + fasta_header
            control = c.replace(">", ":")
        else:
            control = '0:' + fasta_header
    else:
        control = '0'
    if args.include_barcode_in_fastq_header:
        c = control
        control = c + ':' + barcode
    # always return both headers and allow downstream to decide what to output
    return ("@{0}:{1}:{2}:{3}:{4}:{5}:{6} {7}:{8}:{9}".format(args.instrument, '1', args.flow_cell, lane, full_tile, x_pos, y_pos, '1', filter, control),  # nopep8
            "@{0}:{1}:{2}:{3}:{4}:{5}:{6} {7}:{8}:{9}".format(args.instrument, '1', args.flow_cell, lane, full_tile, x_pos, y_pos, '2', filter, control))  # nopep8

if __name__ == '__main__':
    sys.exit(main(sys.argv))
