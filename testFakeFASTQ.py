import unittest
import makeFakeFASTQ
import re
from argparse import Namespace
from Bio.Seq import Seq  # used for reverse_complement()


# DSWF expects paired FASTQ files
# There needs to have 4 reads minimum per sequence location.  Two reads in the seq1 file
# and two reads in the seq2 file.
# Each read has a barcode constructed from the first 10 bp of the seq1 FASTQ file
# and the first 10bp of the seq2 FASTQ file.
# The barcodes can be in either direction seq1 -> seq2 (ab) or seq2 -> seq1 (ba)
# The sequence can be in either orientation seq1 -> seq2 (1) or seq2 -> seq1 (2)

# This leads to four reads:
# ab:1 seq1 contains the a barcode and seq2 contains the b barcode and the - sequence
# ab:2 seq1 contains the a barcode and seq2 contains the b barcode and the + sequence

# ba:1 seq1 contains the b barcode and seq2 contains the a barcode and the + sequence
# ba:2 seq1 contains the b barcode and seq2 contains the a barcode and the - sequence

# The pairs that should have matching sequence are:
# - sequence : ab:1 and ba:2
# + sequence : ab:2 and ba:1

# ex:
# barcode a AGAAAACGTC
# barcode b AAGATTTGTG

# ab1 AGAAAACGTC AACTAGGGCTTGTCCTAGCAGGGAAAGCTTTGACTTAAAAGCCTCTGCT
# ab2 AAGATTTGTG TGGGGGCATCTCTTATACTCATGAAATCAACAGAGGCTTGC

# ba1 AAGATTTGTG TGGGGGCATCTCTTATACTCATGAAATCAACAGAGGCTTGC
# ba2 AGAAAACGTC AACTAGGGCTTGTCCTAGCAGGGAAAGCTTTGACTTAAAAGCCTCTGCT

def avg_qual(seq):
    qual_array = list(map(lambda x: ord(x) - 33, list(seq)))
    avg = (sum(map(int, qual_array)) / len(qual_array))
    return avg


class FakeFASTQTest(unittest.TestCase):
    def setUp(self):
        self.fasta_header = "test1"
        # the default barcodes must be switched to make a full barcode
        self.barcode_a = "AGAAAACGTC"
        self.barcode_b = "AAGATTTGTG"
        self.fullbarcode = self.barcode_a + self.barcode_b
        self.seq1 = "AACTAGGGCTTGTCCTAGCAGGGAAAGCTTTGACTTAAAAGCCTCTGCT"
        self.seq2 = "TGGGGGCATCTCTTATACTCATGAAATCAACAGAGGCTTGC"

        self.rev_seq1 = "AGCAGAGGCTTTTAAGTCAAAGCTTTCCCTGCTAGGACAAGCCCTAGTT"

        self.read_ba1 = self.barcode_b + self.seq1
        self.read_ab2 = self.barcode_a + self.seq1

        self.read_ab1 = self.barcode_a + self.seq2
        self.read_ba2 = self.barcode_b + self.seq2

        # length of sequence only
        self.qual = "HFDEGEFHFIFDGIEHIFGEJJJHG"
        # length of sequence plus barcode (10) and spacer (1)
        self.full_qual = "HFDEGEFHFIFDGIEHIFGEJJJHGEFHFIFDGIE"
        self.buffer_seq = "GGGGG"
        self.args = Namespace(x_min=5, x_max=5, y_min=8, y_max=8,
                              lane_min=1, lane_max=1, swathes=[111],
                              buffer_end=1, buffer_both_sides=0, truncate_by_read=0,
                              wobble_read=None, rand_window=3,
                              truncate_end=1, truncate_both_sides=0, truncate_start=0,
                              tile_min=3, tile_max=3, is_filtered=['N'],
                              include_fasta_header_in_fastq_header=1, include_barcode_in_fastq_header=1,
                              paired_end=0, instrument="N5V", flow_cell="H5N", quality_type="high",
                              read_length=35, spacer_length=1,
                              max_num_reads=1, min_num_reads=1, num_reads=None,
                              max_num_flipped=1, min_num_flipped=1, num_flipped=None,
                              # min_num_families is expected to be 1 in all make_family() tests
                              max_num_families=1, min_num_families=1, num_families=None,
                              map_file=0, tag_file=0)

    # def tearDown(self):
    #    print("teardown")

    def test_reads(self):
        self.assertEqual(self.read_ab1[10:], self.read_ba2[10:], "read ab1 and \
        read ba2 should match after first 10 bp\nab1 {0}\nba2 {1}".format(
            self.read_ab1[10:], self.read_ba2[10:]))
        self.assertEqual(self.read_ab2[10:], self.read_ba1[10:], "read ab2 and \
        read ba1 should match after first 10 bp\nab2 {0}\nba1 {1}".format(
            self.read_ab2[10:], self.read_ba1[10:]))

    # test reverse complement
    def test_rev_seq(self):
        calc_rev_seq = Seq(self.seq1).reverse_complement()
        self.assertEqual(
            self.rev_seq1, calc_rev_seq,
            "self.rev_seq1 and Bio.Seq reverse_complement are the same\
            \nself.rev_seq1 {0}\nBio.seq revc {1}\n"
            .format(self.rev_seq1, calc_rev_seq))

    # TEST def random_sequence
    def test_random_sequence_length(self):
        # print("Testing that random sequence is specified length...")
        self.assertEqual(len(makeFakeFASTQ.random_sequence(10)), 10, "10 length random \
        sequence")

    def test_random_sequence_composition(self):
        # print("Testing that random sequence contains only AGTC...")
        reg = re.compile('^[AGTC]+$')
        seq = makeFakeFASTQ.random_sequence(10)
        self.assertTrue(reg.match(seq), "random_sequence contains only AGTC")

    #############################################################################################

    # TEST def fastq_quality(args, seq_len)
    #  requires args.quality_type where quality_type is high, medium, low
    #  returns a string of seq_len
    def test_fastq_quality_length(self):
        # print("Testing that fastq_quality is specified length...")
        self.assertEqual(len(makeFakeFASTQ.fastq_quality(self.args, 10)), 10, "10 length quality")

    def test_fastq_quality_high(self):
        # print("Testing that fastq_quality high works...")
        qual = makeFakeFASTQ.fastq_quality(self.args, 10)
        avg = avg_qual(qual)
        self.assertTrue(avg >= 35, "qual quality {0} is not greater than or equal to 35".format(avg))

    def test_fastq_quality_medium(self):
        # print("Testing that fastq_quality medium works...")
        args = Namespace(quality_type='medium')
        qual = makeFakeFASTQ.fastq_quality(args, 10)
        avg = avg_qual(qual)
        self.assertTrue(avg <= 30, "qual quality average {0} is less than\
                30".format(avg))
        self.assertTrue(avg >= 23, "qual quality average {0} is not \
                greater than or equal to 23".format(avg))

    def test_fastq_quality_low(self):
        # print("Testing that fastq_quality low works...")
        args = Namespace(quality_type='low')
        qual = makeFakeFASTQ.fastq_quality(args, 10)
        avg = avg_qual(qual)
        self.assertTrue(avg <= 15, "qual quality average {0} is less than\
                15".format(avg))
        self.assertTrue(avg >= 3, "qual quality average {0} is not \
                greater than or equal to 3".format(avg))

    def test_fastq_high_quality_letters(self):
        # print("Testing that fastq_quality result contains proper ascii letters...")
        qual = makeFakeFASTQ.fastq_quality(self.args, 10)
        reg = re.compile('^[A-J]+$')
        self.assertTrue(reg.match(qual), "fastq_quality high contains only ascii letters: {0}".format(qual))

    def test_fastq_low_quality_valid(self):
        # print("Testing that fastq_quality result contains valid values...")
        # valid  !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHI
        args = Namespace(quality_type="low")
        qual = makeFakeFASTQ.fastq_quality(args, 10)
        reg = re.compile('^[A-I0-9!\"#\$\%\&\'\(\)\*\+,=\./:;<=>?@-]+$')
        self.assertTrue(reg.match(qual), "fastq_quality low contains only valid values:{0}".format(qual))

    #############################################################################################

    # TEST def fastq_entry_header(args, fasta_header, barcode)
    def test_fastq_entry_header(self):
        # test header that doesn't include barcode or fasta header
        # print("Testing fastq entry header...")
        self.args.include_barcode_in_fastq_header = 0
        self.args.include_fasta_header_in_fastq_header = 0
        testHeader, testPaired = makeFakeFASTQ.fastq_entry_header(self.args, self.fasta_header, self.barcode_a)
        canonical = "@N5V:1:H5N:1:11103:5:8 1:N:0"
        self.assertEqual(canonical, testHeader, "fastq_entry_header returns proper header")

    def test_fastq_entry_header_fasta_header(self):
        self.args.include_barcode_in_fastq_header = 0
        testHeader, testPaired = makeFakeFASTQ.fastq_entry_header(self.args, self.fasta_header, self.barcode_a)
        canonical = "@N5V:1:H5N:1:11103:5:8 1:N:0:test1"
        self.assertEqual(canonical, testHeader, "fastq_entry_header returns proper header with fasta header no > {0}".format(testHeader))

    def test_fastq_entry_header_fasta_header_symbol(self):
        self.fasta_header = ">test1"
        self.args.include_barcode_in_fastq_header = 0
        testHeader, testPaired = makeFakeFASTQ.fastq_entry_header(self.args, self.fasta_header, self.barcode_a)
        canonical = "@N5V:1:H5N:1:11103:5:8 1:N:0:test1"
        self.assertEqual(canonical, testHeader, "fastq_entry_header returns proper header with fasta header >\ncanonical: {0}\ntest:      {1}".format(canonical, testHeader))

    def test_fastq_entry_header_barcode(self):
        self.args.include_fasta_header_in_fastq_header = 0
        testHeader, testPaired = makeFakeFASTQ.fastq_entry_header(self.args, self.fasta_header, self.barcode_a)
        canonical = "@N5V:1:H5N:1:11103:5:8 1:N:0:AGAAAACGTC"
        self.assertEqual(canonical, testHeader, "fastq_entry_header returns proper header with barcode canonical: {0} test:".format(canonical, testHeader))

    def test_fastq_entry_header_fasta_header_and_barcode(self):
        testHeader, testPaired = makeFakeFASTQ.fastq_entry_header(self.args, self.fasta_header, self.barcode_a)
        canonical = "@N5V:1:H5N:1:11103:5:8 1:N:0:test1:AGAAAACGTC"
        self.assertEqual(canonical, testHeader, "fastq_entry_header returns proper header with fasta header and barcode {0}".format(testHeader))

    def test_fastq_entry_header_fasta_header_symbol_and_barcode(self):
        self.fasta_header = ">test1"
        testHeader, testPaired = makeFakeFASTQ.fastq_entry_header(self.args, self.fasta_header, self.barcode_a)
        canonical = "@N5V:1:H5N:1:11103:5:8 1:N:0:test1:AGAAAACGTC"
        self.assertEqual(canonical, testHeader, "fastq_entry_header returns proper header with fasta header with > and barcode {0}".format(testHeader))

    #############################################################################################

    # TEST def buffer_sequence(args,  seq)
    def test_buffer_sequence_end_nothing(self):
        # defaults to buffer_end
        self.args.buffer_seq = self.buffer_seq
        testSeq = makeFakeFASTQ.buffer_sequence(self.args, self.seq1)
        self.assertEqual(len(testSeq), len(self.seq1), "buffer_sequence does nothing when it should args: {0} test:{1}".format(len(self.seq1), len(testSeq)))

    def test_buffer_sequence_end_nothing_rl(self):
        # defaults to buffer_end
        self.args.buffer_seq = self.buffer_seq
        self.args.read_length = len(self.seq1) + 5
        testSeq = makeFakeFASTQ.buffer_sequence(self.args, self.seq1)
        self.assertEqual(len(testSeq), self.args.read_length, "buffer_sequence does nothing when it should args: {0} test:{1}".format(self.args.read_length, len(testSeq)))

    def test_buffer_sequence_end(self):
        # defaults to buffer_end, should add 5 G's to end of self.seq1
        self.args.read_length = len(self.seq1) + 5
        self.args.buffer_seq = self.buffer_seq
        testSeq = makeFakeFASTQ.buffer_sequence(self.args, self.seq1)
        self.assertEqual(len(testSeq), self.args.read_length, "buffer sequence outputs correct read length args: {0} test:{1}".format(self.args.read_length, len(testSeq)))
        canonical = "AACTAGGGCTTGTCCTAGCAGGGAAAGCTTTGACTTAAAAGCCTCTGCTGGGGG"
        self.assertEqual(canonical, testSeq, "buffer sequence outputs correct read sequence".format(testSeq))

    def test_buffer_sequence_both_sides_nothing(self):
        # with buffer_both_sides set, it shouldn't do anything when the read length is correct
        self.args.buffer_end = 0
        self.args.buffer_both_sides = 1
        self.args.buffer_seq = self.buffer_seq
        self.args.read_length = len(self.seq1)
        testSeq = makeFakeFASTQ.buffer_sequence(self.args, self.seq1)
        self.assertEqual(len(testSeq), self.args.read_length, "buffer_sequence does nothing when it should args: {0} test:{1}".format(self.args.read_length, len(testSeq)))

    def test_buffer_sequence_both_sides_odd(self):
        # with buffer_both_sides set, it should add to both sides
        self.args.read_length = 54  # add 5 bp, 3 to left, 2 to right
        self.args.buffer_end = 0
        self.args.buffer_both_sides = 1
        self.args.buffer_seq = self.buffer_seq
        testSeq = makeFakeFASTQ.buffer_sequence(self.args, self.seq1)
        self.assertEqual(len(testSeq), self.args.read_length, "buffer_sequence proper length with odd addition\nLengths: args: {0} test:{1}\nSeqs:\ntest: {2}\n seq: {3}".format(self.args.read_length, len(testSeq), testSeq, self.seq1))
        canonical = "GGGAACTAGGGCTTGTCCTAGCAGGGAAAGCTTTGACTTAAAAGCCTCTGCTGG"
        self.assertEqual(testSeq, canonical, "buffer_sequence sequence is canonical\nLengths: args: {0} test:{1}\nSeqs:\ntest: {2}\n seq: {3}".format(self.args.read_length, len(testSeq), testSeq, canonical))

    def test_buffer_sequence_both_sides_even(self):
        self.args.read_length = len(self.seq1) + 4
        self.args.buffer_end = 0
        self.args.buffer_both_sides = 1
        self.args.buffer_seq = self.buffer_seq
        testSeq = makeFakeFASTQ.buffer_sequence(self.args, self.seq1)
        self.assertEqual(len(testSeq), self.args.read_length, "buffer_sequence proper length with even addition\nLengths: args: {0} test:{1}\nSeqs:\ntest: {2}\n seq: {3}".format(self.args.read_length, len(testSeq), testSeq, self.seq1))
        canonical = "GGAACTAGGGCTTGTCCTAGCAGGGAAAGCTTTGACTTAAAAGCCTCTGCTGG"
        self.assertEqual(testSeq, canonical, "buffer_sequence sequence is canonical\nLengths: args: {0} test:{1}\nSeqs:\ntest: {2}\n seq: {3}".format(self.args.read_length, len(testSeq), testSeq, canonical))

#############################################################################################

    # TEST def truncate_sequence(args,seq)
    def test_truncate_sequence_end_nothing(self):
        # defaults to truncate_end
        testSeq = makeFakeFASTQ.truncate_sequence(self.args, self.seq1)
        self.assertEqual(len(testSeq), self.args.read_length, "truncate_sequence does nothing when it should args: {0} test:{1}".format(self.args.read_length, len(testSeq)))

    def test_truncate_sequence_end(self):
        # defaults to truncate_end, should remove 5 BP from the end of self.seq1
        self.args.read_length = 20
        canonical = "AACTAGGGCTTGTCCTAGCA"
        testSeq = makeFakeFASTQ.truncate_sequence(self.args, self.seq1)
        self.assertEqual(len(testSeq), self.args.read_length, "truncate sequence outputs correct read length.\nLengths: args: {0} test:{1}\nSeqs:\ncanonical: {2}\n     test: {3}".format(self.args.read_length, len(testSeq), canonical, testSeq))
        self.assertEqual(canonical, testSeq, "truncate sequence outputs correct read sequence\ncanonical: {0}\n     test: {1}".format(canonical, testSeq))

    def test_truncate_sequence_start(self):
        self.args.read_length = 20
        self.args.truncate_end = 0
        self.args.truncate_start = 1
        canonical = "TTGACTTAAAAGCCTCTGCT"
        testSeq = makeFakeFASTQ.truncate_sequence(self.args, self.seq1)
        self.assertEqual(len(testSeq), self.args.read_length, "truncate sequence start outputs correct read length.\nLengths: args: {0} test:{1}\nSeqs:\ncanonical: {2}\n     test: {3}".format(self.args.read_length, len(testSeq), canonical, testSeq))
        self.assertEqual(canonical, testSeq, "truncate sequence start outputs correct read sequence\ncanonical: {0}\n     test: {1}".format(canonical, testSeq))

    def test_truncate_sequence_both_sides_nothing(self):
        # with truncate_both_sides set, it shouldn't do anything when the read length is correct
        self.args.truncate_end = 0
        self.args.truncate_both_sides = 1
        testSeq = makeFakeFASTQ.truncate_sequence(self.args, self.seq1)
        self.assertEqual(len(testSeq), self.args.read_length, "truncate_sequence does nothing when it should args: {0} test:{1}".format(self.args.read_length, len(testSeq)))

    def test_truncate_sequence_both_sides_odd(self):
        self.args.read_length = 20  # seq is 25 bp, if updated then canonical needs updating too
        self.args.truncate_end = 0
        self.args.truncate_both_sides = 1
        testSeq = makeFakeFASTQ.truncate_sequence(self.args, self.seq1)
        self.assertEqual(len(testSeq), self.args.read_length, "truncate_sequence proper length with odd addition\nLengths: args: {0} test:{1}\nSeqs:\ntest: {2}\n seq: {3}".format(self.args.read_length, len(testSeq), testSeq, self.seq1))
        canonical = "TAGCAGGGAAAGCTTTGACT"
        self.assertEqual(testSeq, canonical, "truncate_sequence sequence is canonical\nLengths: args: {0} test:{1}\nSeqs:\ntest: {2}\n seq: {3}".format(self.args.read_length, len(testSeq), testSeq, canonical))

    def test_truncate_sequence_both_sides_even(self):
        self.args.read_length = 21  # seq is 25 bp
        self.args.truncate_end = 0
        self.args.truncate_both_sides = 1
        testSeq = makeFakeFASTQ.truncate_sequence(self.args, self.seq1)
        self.assertEqual(len(testSeq), self.args.read_length, "truncate_sequence proper length with even addition\nLengths: args: {0} test:{1}\nSeqs:\ntest: {2}\n seq: {3}".format(self.args.read_length, len(testSeq), testSeq, self.seq1))
        canonical = "CTAGCAGGGAAAGCTTTGACT"
        self.assertEqual(testSeq, canonical, "truncate_sequence sequence is canonical\nLengths: args: {0} test:{1}\nSeqs:\ntest: {2}\n seq: {3}".format(self.args.read_length, len(testSeq), testSeq, canonical))

    def test_truncate_sequence_end_by_read(self):
        self.args.read_length = 20
        self.args.truncate_by_read = 1
        canonical = "AACTAGGGCTTGTCCTAGCA"
        testSeq = makeFakeFASTQ.truncate_sequence(self.args, self.seq1, len(self.seq1) - self.args.read_length, '5')
        self.assertEqual(len(testSeq), self.args.read_length, "truncate sequence end outputs correct read length.\nLengths: args: {0} test:{1}\nSeqs:\ncanonical: {2}\n     test: {3}".format(self.args.read_length, len(testSeq), canonical, testSeq))
        self.assertEqual(canonical, testSeq, "truncate sequence end outputs correct read sequence\ncanonical: {0}\n     test: {1}".format(canonical, testSeq))

    def test_truncate_sequence_start_by_read(self):
        self.args.read_length = 20
        self.args.truncate_by_read = 1
        self.args.rand_window = 0
        canonical = "TTGACTTAAAAGCCTCTGCT"
        testSeq = makeFakeFASTQ.truncate_sequence(self.args, self.seq1, len(self.seq1) - self.args.read_length, '3')
        self.assertEqual(len(testSeq), self.args.read_length, "truncate sequence start outputs correct read length.\nLengths: args: {0} test:{1}\nSeqs:\ncanonical: {2}\n     test: {3}".format(self.args.read_length, len(testSeq), canonical, testSeq))
        self.assertEqual(canonical, testSeq, "truncate sequence start outputs correct read sequence\ncanonical: {0}\n     test: {1}".format(canonical, testSeq))

    def test_truncate_sequence_start_by_read_window(self):
        self.args.read_length = 20
        self.args.truncate_by_read = 1
        canonical = "TTGACTTAAAAGCCTCTGCT"
        testSeq = makeFakeFASTQ.truncate_sequence(self.args, self.seq1, len(self.seq1) - self.args.read_length, '3')
        self.assertEqual(len(testSeq), self.args.read_length, "truncate sequence start outputs correct read length.\nLengths: args: {0} test:{1}\nSeqs:\ncanonical: {2}\n     test: {3}".format(self.args.read_length, len(testSeq), canonical, testSeq))
        self.assertEqual(canonical, testSeq, "truncate sequence start outputs correct read sequence\ncanonical: {0}\n     test: {1}".format(canonical, testSeq))

    def test_truncate_sequence_start_by_read_failure_read_type_not_set(self):
        self.args.read_length = 20
        self.args.truncate_end = 0
        self.args.truncate_start = 0
        self.args.truncate_by_read = 1
        try:
            testseq = makeFakeFASTQ.truncate_sequence(self.args, self.seq1, len(self.seq1) - self.args.read_length)
        except ValueError:
            pass
        except Exception as e:
            self.fail("Unexpected exception raised:", e)
        else:
            self.fail("Expected ValueError not raised")

    #############################################################################################

    # TEST def make_ds_read(args, seq, barcode, read_type)
    def test_make_ds_read_truncate(self):
        testRead = makeFakeFASTQ.make_ds_read(self.args, self.seq1, self.barcode_a)
        canonical = "AGAAAACGTCTAACTAGGGCTTGTCCTAGCAGGGA"
        self.assertEqual(len(testRead), self.args.read_length, "make_ds_read read length is proper. canonical: {0} test: {1}".format(self.args.read_length, len(testRead)))
        self.assertEqual(testRead, canonical, "make_ds_read sequence is canonical\ncanonical: {0}\n     test: {1}".format(canonical, testRead))

    def test_make_ds_read_equal(self):
        # test exactly equal read_length and barcode+seq+spacer
        self.args.read_length = 36
        testRead = makeFakeFASTQ.make_ds_read(self.args, self.seq1, self.barcode_a)
        # "AGAAAACGTC"+"T"+"AACTAGGGCTTGTCCTAGCAGGGAAAGCTTTGACTTAAAAGCCTCTGCT" => "AGAAAACGTC"+ "T" +"AACTAGGGCTTGTCCTAGCAGGGA"
        canonical = "AGAAAACGTCTAACTAGGGCTTGTCCTAGCAGGGAA"
        self.assertEqual(len(testRead), self.args.read_length, "make_ds_read read length is proper. canonical: {0} test: {1}".format(self.args.read_length, len(testRead)))
        self.assertEqual(testRead, canonical, "make_ds_read sequence is canonical\ncanonical: {0}\n     test: {1}".format(canonical, testRead))

    def test_make_ds_read_buffer(self):
        # test read_length > barcode+seq+spacer
        self.args.read_length = 40
        self.args.buffer_seq = self.buffer_seq
        testRead = makeFakeFASTQ.make_ds_read(self.args, self.seq1, self.barcode_a)
        # 10 + 1 + 29 = 40
        # "AGAAAACGTC"+"T"+"AACTAGGGCTTGTCCTAGCAGGGAAAGCTTTGACTTAAAAGCCTCTGCT" => "AGAAAACGTC"+ "T" +"AACTAGGGCTTGTCCTAGCAGGGAAAGCT"
        canonical = "AGAAAACGTCTAACTAGGGCTTGTCCTAGCAGGGAAAGCT"
        self.assertEqual(len(testRead), self.args.read_length, "make_ds_read read length is proper. args: {0} canonical: {1} test: {2}\ncanonical: {3}\n     test: {4}".format(self.args.read_length, len(canonical), len(testRead), canonical, testRead))
        self.assertEqual(testRead, canonical, "make_ds_read sequence is canonical\ncanonical: {0}\n     test: {1}".format(canonical, testRead))

    def test_make_ds_read_forward(self):
        # the 5 shouldn't do anything unless args.truncate_by_read is set to 1
        # forward reads should truncate the end
        self.args.truncate_end = 1
        self.args.truncate_start = 0
        testRead = makeFakeFASTQ.make_ds_read(self.args, self.seq1, self.barcode_a, '5')
        # "AGAAAACGTC"+"T"+"AACTAGGGCTTGTCCTAGCAGGGAAAGCTTTGACTTAAAAGCCTCTGCT" => "AGAAAACGTC"+ "T" +"AACTAGGGCTTGTCCTAGCAGGGA"
        canonical = "AGAAAACGTCTAACTAGGGCTTGTCCTAGCAGGGA"
        self.assertEqual(len(testRead), self.args.read_length, "make_ds_read read length is proper. canonical: {0} test: {1}".format(self.args.read_length, len(testRead)))
        self.assertEqual(testRead, canonical, "make_ds_read sequence is canonical\ncanonical: {0}\n     test: {1}".format(canonical, testRead))

    def test_make_ds_read_forward_by_read(self):
        # 5 should truncate the end
        self.args.truncate_end = 0
        self.args.truncate_by_read = 1
        testRead = makeFakeFASTQ.make_ds_read(self.args, self.seq1, self.barcode_a, '5')
        canonical = "AGAAAACGTCTAACTAGGGCTTGTCCTAGCAGGGA"
        self.assertEqual(len(testRead), self.args.read_length, "make_ds_read read length is proper. canonical: {0} test: {1}".format(self.args.read_length, len(testRead)))
        self.assertEqual(testRead, canonical, "make_ds_read sequence is canonical\ncanonical: {0}\n     test: {1}".format(canonical, testRead))

    def test_make_ds_read_forward_by_read_ignore_start(self):
        # 5 should truncate the end
        self.args.truncate_end = 0
        self.args.truncate_start = 1
        self.args.truncate_by_read = 1
        testRead = makeFakeFASTQ.make_ds_read(self.args, self.seq1, self.barcode_a, '5')
        canonical = "AGAAAACGTCTAACTAGGGCTTGTCCTAGCAGGGA"
        self.assertEqual(len(testRead), self.args.read_length, "read length is proper. canonical: {0} test: {1}".format(self.args.read_length, len(testRead)))
        self.assertEqual(testRead, canonical, "sequence is canonical\ncanonical: {0}\n     test: {1}".format(canonical, testRead))

    def test_make_ds_read_reverse_by_read(self):
        # 3 should truncate the start
        self.args.truncate_end = 0
        self.args.truncate_by_read = 1
        self.args.rand_window = 0
        testRead = makeFakeFASTQ.make_ds_read(self.args, self.seq1, self.barcode_a, '3')
        # "AGAAAACGTC"+"T"+"AACTAGGGCTTGTCCTAGCAGGGAAAGCTTTGACTTAAAAGCCTCTGCT" => "AGAAAACGTC"+"T"+"AGCTTTGACTTAAAAGCCTCTGCT"
        canonical = "AGAAAACGTCTAGCTTTGACTTAAAAGCCTCTGCT"
        self.assertEqual(len(testRead), self.args.read_length, "make_ds_read read length is proper. canonical: {0} test: {1}".format(self.args.read_length, len(testRead)))
        self.assertEqual(testRead, canonical, "make_ds_read sequence is canonical\ncanonical: {0}\n     test: {1}".format(canonical, testRead))

    def test_make_ds_read_reverse_by_read_window(self):
        # 3 should truncate the start
        self.args.truncate_end = 0
        self.args.truncate_by_read = 1
        testRead = makeFakeFASTQ.make_ds_read(self.args, self.seq1, self.barcode_a, '3')
        # "AGAAAACGTC"+"T"+"AACTAGGGCTTGTCCTAGCAGGGAAAGCTTTGACTTAAAAGCCTCTGCT" => "AGAAAACGTC"+"T"+"GCTTTGACTTAAAAGCCTCTGCT"
        canonical = "AGAAAACGTCTAGCTTTGACTTAAAAGCCTCTGCT"
        self.assertEqual(len(testRead), self.args.read_length, "make_ds_read read length is proper. canonical: {0} test: {1}".format(self.args.read_length, len(testRead)))
        self.assertEqual(testRead, canonical, "make_ds_read sequence is canonical\ncanonical: {0}\n     test: {1}".format(canonical, testRead))

    def test_make_ds_read_reverse_by_read_all_params(self):
        # 3 should truncate the start by the read and ignore the truncate_end and truncate_start params
        self.args.truncate_end = 1
        self.args.truncate_start = 1
        self.args.truncate_by_read = 1
        self.args.rand_window = 0
        testRead = makeFakeFASTQ.make_ds_read(self.args, self.seq1, self.barcode_a, '3')
        # "AGAAAACGTC"+"T"+"AACTAGGGCTTGTCCTAGCAGGGAAAGCTTTGACTTAAAAGCCTCTGCT" => "AGAAAACGTC"+"T"+"AGCTTTGACTTAAAAGCCTCTGCT"
        canonical = "AGAAAACGTCTAGCTTTGACTTAAAAGCCTCTGCT"
        self.assertEqual(len(testRead), self.args.read_length, "make_ds_read read length is proper. canonical: {0} test: {1}".format(self.args.read_length, len(testRead)))
        self.assertEqual(testRead, canonical, "make_ds_read sequence is canonical\ncanonical: {0}\n     test: {1}".format(canonical, testRead))

    def test_make_ds_read_reverse_by_read_all_params_window(self):
        # 3 should truncate the start by the read and ignore the truncate_end and truncate_start params
        self.args.truncate_end = 1
        self.args.truncate_start = 1
        self.args.truncate_by_read = 1
        testRead = makeFakeFASTQ.make_ds_read(self.args, self.seq1, self.barcode_a, '3')
        # "AGAAAACGTC"+"T"+"AACTAGGGCTTGTCCTAGCAGGGAAAGCTTTGACTTAAAAGCCTCTGCT" => "AGAAAACGTC"+"T"+"AGCTTTGACTTAAAAGCCTCTGCT"
        canonical = "AGAAAACGTCTAGCTTTGACTTAAAAGCCTCTGCT"
        self.assertEqual(len(testRead), self.args.read_length, "make_ds_read read length is proper. canonical: {0} test: {1}".format(self.args.read_length, len(testRead)))
        self.assertEqual(testRead, canonical, "make_ds_read sequence is canonical\ncanonical: {0}\n     test: {1}".format(canonical, testRead))

    def test_make_ds_read_no_truncates(self):
        # 3 should truncate the start
        self.args.truncate_end = 0
        self.args.truncate_start = 0
        self.args.truncate_by_read = 0
        try:
            testseq = makeFakeFASTQ.make_ds_read(self.args, self.seq1, self.barcode_a, '3')
        except ValueError:
            pass
        except Exception as e:
            self.fail("Unexpected exception raised:", e)
        else:
            self.fail("Expected ValueError not raised")

    def test_make_ds_read_end_three(self):
        # 3 should truncate the start, should ignore truncate_end
        self.args.truncate_end = 1
        self.args.truncate_start = 0
        self.args.truncate_by_read = 1
        self.args.rand_window = 0
        testRead = makeFakeFASTQ.make_ds_read(self.args, self.seq1, self.barcode_a, '3')
        # "AGAAAACGTC"+"T"+"AACTAGGGCTTGTCCTAGCAGGGAAAGCTTTGACTTAAAAGCCTCTGCT" => "AGAAAACGTC"+"T"+"AGCTTTGACTTAAAAGCCTCTGCT"
        canonical = "AGAAAACGTCTAGCTTTGACTTAAAAGCCTCTGCT"
        self.assertEqual(len(testRead), self.args.read_length, "make_ds_read read length is proper. canonical: {0} test: {1}".format(self.args.read_length, len(testRead)))
        self.assertEqual(testRead, canonical, "make_ds_read sequence is canonical\ncanonical: {0}\n     test: {1}".format(canonical, testRead))

    def test_make_ds_read_end_three_window(self):
        # 3 should truncate the start, should ignore truncate_end
        # rand_window is 3, so move it 3 to the left
        self.args.truncate_end = 1
        self.args.truncate_start = 0
        self.args.truncate_by_read = 1
        testRead = makeFakeFASTQ.make_ds_read(self.args, self.seq1, self.barcode_a, '3')
        # "AGAAAACGTC"+"T"+"AACTAGGGCTTGTCCTAGCAGGGAAAGCTTTGACTTAAAAGCCTCTGCT" => "AGAAAACGTC"+"T"+"AGCTTTGACTTAAAAGCCTCTGCT"
        canonical = "AGAAAACGTCTAGCTTTGACTTAAAAGCCTCTGCT"
        self.assertEqual(len(testRead), self.args.read_length, "make_ds_read read length is proper. canonical: {0} test: {1}".format(self.args.read_length, len(testRead)))
        self.assertEqual(testRead, canonical, "make_ds_read sequence is canonical\ncanonical: {0}\n     test: {1}".format(canonical, testRead))

    #############################################################################################

    # TEST def make_family(header, seq, args)
    def test_make_family_truncate(self):
        self.args.quality = self.qual
        self.args.barcode_a = self.barcode_a
        self.args.barcode_b = self.barcode_b
        (testFam1, testFam2) = makeFakeFASTQ.make_family(self.fasta_header, self.seq1, self.args, self.args.min_num_families)
        canonical1 = ['@N5V:1:H5N:1:11103:5:8 1:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AAGATTTGTGTAACTAGGGCTTGTCCTAGCAGGGA', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHG', '@N5V:1:H5N:1:11103:5:8 1:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AGAAAACGTCTAACTAGGGCTTGTCCTAGCAGGGA', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHG']
        self.assertEqual(' '.join(canonical1), ' '.join(testFam1), "make_family returns proper values for family 1.\ncanonical: {0}\n     test: {1} ".format(canonical1, testFam1))
        canonical2 = ['@N5V:1:H5N:1:11103:5:8 2:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AGAAAACGTCTAACTAGGGCTTGTCCTAGCAGGGA', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHG', '@N5V:1:H5N:1:11103:5:8 2:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AAGATTTGTGTAACTAGGGCTTGTCCTAGCAGGGA', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHG']
        self.assertEqual(canonical2, testFam2, "make_family returns proper values for family 2.\ncanonical: {0}\n     test: {1} ".format(canonical2, testFam2))

    def test_make_family_equal(self):
        self.args.quality = self.full_qual
        self.args.read_length = 36
        self.args.barcode_a = self.barcode_a
        self.args.barcode_b = self.barcode_b
        (testFam1, testFam2) = makeFakeFASTQ.make_family(self.fasta_header, self.seq1, self.args, self.args.min_num_families)
        canonical1 = ['@N5V:1:H5N:1:11103:5:8 1:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AAGATTTGTGTAACTAGGGCTTGTCCTAGCAGGGAA', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHGEFHFIFDGIE', '@N5V:1:H5N:1:11103:5:8 1:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AGAAAACGTCTAACTAGGGCTTGTCCTAGCAGGGAA', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHGEFHFIFDGIE']
        self.assertEqual(' '.join(canonical1), ' '.join(testFam1), "make_family returns proper values for family 1.\ncanonical: {0}\n     test: {1} ".format(canonical1, testFam1))
        canonical2 = ['@N5V:1:H5N:1:11103:5:8 2:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AGAAAACGTCTAACTAGGGCTTGTCCTAGCAGGGAA', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHGEFHFIFDGIE', '@N5V:1:H5N:1:11103:5:8 2:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AAGATTTGTGTAACTAGGGCTTGTCCTAGCAGGGAA', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHGEFHFIFDGIE']
        self.assertEqual(canonical2, testFam2, "make_family returns proper values for family 2.\ncanonical: {0}\n     test: {1} ".format(canonical2, testFam2))

    def test_make_family_buffer(self):
        self.args.quality = self.full_qual + 'GGGG'
        self.args.barcode_a = self.barcode_a
        self.args.barcode_b = self.barcode_b
        self.args.buffer_seq = self.buffer_seq
        self.args.read_length = 40
        (testFam1, testFam2) = makeFakeFASTQ.make_family(self.fasta_header, self.seq1, self.args, self.args.min_num_families)
        canonical1 = ['@N5V:1:H5N:1:11103:5:8 1:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AAGATTTGTGTAACTAGGGCTTGTCCTAGCAGGGAAAGCT', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHGEFHFIFDGIEGGGG', '@N5V:1:H5N:1:11103:5:8 1:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AGAAAACGTCTAACTAGGGCTTGTCCTAGCAGGGAAAGCT', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHGEFHFIFDGIEGGGG']
        self.assertEqual(' '.join(canonical1), ' '.join(testFam1), "make_family returns proper values for family 1.\ncanonical: {0}\n     test: {1} ".format(canonical1, testFam1))
        canonical2 = ['@N5V:1:H5N:1:11103:5:8 2:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AGAAAACGTCTAACTAGGGCTTGTCCTAGCAGGGAAAGCT', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHGEFHFIFDGIEGGGG', '@N5V:1:H5N:1:11103:5:8 2:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AAGATTTGTGTAACTAGGGCTTGTCCTAGCAGGGAAAGCT', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHGEFHFIFDGIEGGGG']
        self.assertEqual(canonical2, testFam2, "make_family returns proper values for family 2.\ncanonical: {0}\n     test: {1} ".format(canonical2, testFam2))

    def test_make_family_truncate_by_read(self):
        self.args.quality = self.qual
        self.args.barcode_a = self.barcode_a
        self.args.barcode_b = self.barcode_b
        self.args.truncate_end = 0
        self.args.read_length = 25
        self.args.truncate_by_read = 1
        (testFam1, testFam2) = makeFakeFASTQ.make_family(self.fasta_header, self.seq1, self.args, self.args.min_num_families)
        canonical1 = ['@N5V:1:H5N:1:11103:5:8 1:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AAGATTTGTGTAACTAGGGCTTGTC', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHG', '@N5V:1:H5N:1:11103:5:8 1:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AGAAAACGTCTTAAAAGCCTCTGCT', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHG']
        self.assertEqual(' '.join(canonical1), ' '.join(testFam1), "make_family returns proper values for family 1.\ncanonical: {0}\n     test: {1} ".format(canonical1, testFam1))
        # read pairs are a+end/b+start for canonical 2
        canonical2 = ['@N5V:1:H5N:1:11103:5:8 2:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AGAAAACGTCTTAAAAGCCTCTGCT', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHG', '@N5V:1:H5N:1:11103:5:8 2:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AAGATTTGTGTAACTAGGGCTTGTC', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHG']
        self.assertEqual(canonical2, testFam2, "make_family returns proper values for family 2.\ncanonical: {0}\n     test: {1} ".format(canonical2, testFam2))

    def test_make_family_truncate_by_read_wobble(self):
        # truncate by read means that one sequence gets the 5' of the sequence
        # and the other sequence gets the 3' of the sequence
        self.args.quality = self.qual
        self.args.barcode_a = self.barcode_a
        self.args.barcode_b = self.barcode_b
        self.args.read_length = 25
        self.args.truncate_end = 0
        self.args.truncate_by_read = 1
        self.args.rand_window = 3
        self.args.wobble_read = 1
        (testFam1, testFam2) = makeFakeFASTQ.make_family(self.fasta_header, self.seq1, self.args, self.args.min_num_families)
        # barcode a + start of seq
        # "AGAAAACGTC"+"T"+"AACTAGGGCTTGTCCTAGCAGGGAAAGCTTTGACTTAAAAGCCTCTGCT" => "AGAAAACGTC"+"T"+"AACTAGGGCTTGTC"
        # barcode a + end of seq
        # "AGAAAACGTC"+"T"+"AACTAGGGCTTGTCCTAGCAGGGAAAGCTTTGACTTAAAAGCCTCTGCT" => "AGAAAACGTC"+"T"+ "ACTTAAAAGCCTCT"
        # barcode b + start of seq
        # "AAGATTTGTG"+"T"+"AACTAGGGCTTGTCCTAGCAGGGAAAGCTTTGACTTAAAAGCCTCTGCT" => "AAGATTTGTG"+"T"+"AACTAGGGCTTGTC"
        # barcode b + end of seq
        # "AAGATTTGTG"+"T"+"AACTAGGGCTTGTCCTAGCAGGGAAAGCTTTGACTTAAAAGCCTCTGCT" => "AAGATTTGTG"+"T"+"ACTTAAAAGCCTCT"
        # the headers are expected to be the same in each set because the default parameters for testing are set to a single number. In practice the first entry in the array and the second entry in the array would have different read names with testFam1 having the 1:N:0 values and testFame2 having the 2:N:0 values
        # read pairs are b+start/a+end for canonical 1
        canonical1 = ['@N5V:1:H5N:1:11103:5:8 1:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AAGATTTGTGTAACTAGGGCTTGTC', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHG', '@N5V:1:H5N:1:11103:5:8 1:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AGAAAACGTCTACTTAAAAGCCTCT', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHG']
        self.assertEqual(' '.join(canonical1), ' '.join(testFam1), "make_family returns proper values for family 1.\ncanonical: {0}\n     test: {1} ".format(canonical1, testFam1))
        # read pairs are a+end/b+start for canonical 2
        canonical2 = ['@N5V:1:H5N:1:11103:5:8 2:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AGAAAACGTCTACTTAAAAGCCTCT', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHG', '@N5V:1:H5N:1:11103:5:8 2:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AAGATTTGTGTAACTAGGGCTTGTC', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHG']
        self.assertEqual(canonical2, testFam2, "make_family returns proper values for family 2.\ncanonical: {0}\n     test: {1} ".format(canonical2, testFam2))

    ###########################################################################
    # TEST make_clan
    def test_make_clan(self):
        self.args.quality = self.full_qual
        self.args.read_length = 36
        self.args.barcode_a = self.barcode_a
        self.args.barcode_b = self.barcode_b
        (testClan1, testClan2) = makeFakeFASTQ.make_clan(self.fasta_header, self.seq1, self.args)
        canonical1 = [['@N5V:1:H5N:1:11103:5:8 1:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AAGATTTGTGTAACTAGGGCTTGTCCTAGCAGGGAA', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHGEFHFIFDGIE', '@N5V:1:H5N:1:11103:5:8 1:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AGAAAACGTCTAACTAGGGCTTGTCCTAGCAGGGAA', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHGEFHFIFDGIE']]
        canonical2 = [['@N5V:1:H5N:1:11103:5:8 2:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AGAAAACGTCTAACTAGGGCTTGTCCTAGCAGGGAA', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHGEFHFIFDGIE', '@N5V:1:H5N:1:11103:5:8 2:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AAGATTTGTGTAACTAGGGCTTGTCCTAGCAGGGAA', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHGEFHFIFDGIE']]
        self.assertEqual(' '.join(map(" ".join, canonical1)), ' '.join(map(' '.join, testClan1)), "make_clan returns proper values for clan 1.\ncanonical: {0}\n     test: {1} ".format(canonical1, testClan1))
        self.assertEqual(canonical2, testClan2, "make_clan returns proper values for clan 2.\ncanonical: {0}\n     test: {1} ".format(canonical2, testClan2))

    def test_min_num_families_fails_without_max_num_families(self):
        self.args.quality = self.full_qual
        self.args.read_length = 36
        self.args.barcode_a = self.barcode_a
        self.args.barcode_b = self.barcode_b
        self.args.min_num_families = 3
        try:
            (testClan1, testClan2) = makeFakeFASTQ.make_clan(self.fasta_header, self.seq1, self.args)
        except ValueError:
            pass
        except Exception as e:
            self.fail("Unexpected exception raised:", e)
        else:
            self.fail("Expected ValueError not raised")

    # make clan should return 3 families with 2 reads each for each testClan
    def test_min_num_families(self):
        self.args.quality = self.full_qual
        self.args.read_length = 36
        self.args.barcode_a = self.barcode_a
        self.args.barcode_b = self.barcode_b
        self.args.min_num_families = 3
        self.args.max_num_families = 3
        (testClan1, testClan2) = makeFakeFASTQ.make_clan(self.fasta_header, self.seq1, self.args)
        canonical1 = [['@N5V:1:H5N:1:11103:5:8 1:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AAGATTTGTGTAACTAGGGCTTGTCCTAGCAGGGAA', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHGEFHFIFDGIE', '@N5V:1:H5N:1:11103:5:8 1:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AGAAAACGTCTAACTAGGGCTTGTCCTAGCAGGGAA', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHGEFHFIFDGIE'], ['@N5V:1:H5N:1:11103:5:8 1:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AAGATTTGTGTAACTAGGGCTTGTCCTAGCAGGGAA', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHGEFHFIFDGIE', '@N5V:1:H5N:1:11103:5:8 1:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AGAAAACGTCTAACTAGGGCTTGTCCTAGCAGGGAA', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHGEFHFIFDGIE'], ['@N5V:1:H5N:1:11103:5:8 1:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AAGATTTGTGTAACTAGGGCTTGTCCTAGCAGGGAA', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHGEFHFIFDGIE', '@N5V:1:H5N:1:11103:5:8 1:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AGAAAACGTCTAACTAGGGCTTGTCCTAGCAGGGAA', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHGEFHFIFDGIE']]
        canonical2 = [['@N5V:1:H5N:1:11103:5:8 2:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AGAAAACGTCTAACTAGGGCTTGTCCTAGCAGGGAA', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHGEFHFIFDGIE', '@N5V:1:H5N:1:11103:5:8 2:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AAGATTTGTGTAACTAGGGCTTGTCCTAGCAGGGAA', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHGEFHFIFDGIE'], ['@N5V:1:H5N:1:11103:5:8 2:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AGAAAACGTCTAACTAGGGCTTGTCCTAGCAGGGAA', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHGEFHFIFDGIE', '@N5V:1:H5N:1:11103:5:8 2:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AAGATTTGTGTAACTAGGGCTTGTCCTAGCAGGGAA', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHGEFHFIFDGIE'], ['@N5V:1:H5N:1:11103:5:8 2:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AGAAAACGTCTAACTAGGGCTTGTCCTAGCAGGGAA', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHGEFHFIFDGIE', '@N5V:1:H5N:1:11103:5:8 2:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AAGATTTGTGTAACTAGGGCTTGTCCTAGCAGGGAA', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHGEFHFIFDGIE']]
        self.assertEqual(' '.join(map(' '.join, canonical1)), ' '.join(map(' '.join, testClan1)), "make_clan returns proper values for clan 1.\ncanonical: {0}\n     test: {1} ".format(canonical1, testClan1))
        self.assertEqual(canonical2, testClan2, "make_clan returns proper values for clan 2.\ncanonical: {0}\n     test: {1} ".format(canonical2, testClan2))

    def test_make_clan_num_families(self):
        self.args.quality = self.full_qual
        self.args.read_length = 36
        self.args.barcode_a = self.barcode_a
        self.args.barcode_b = self.barcode_b
        self.args.num_families = 3
        (testClan1, testClan2) = makeFakeFASTQ.make_clan(self.fasta_header, self.seq1, self.args)
        self.assertEqual(len(testClan1), self.args.num_families, "make_clan returns proper number of lines")
        self.assertEqual(len(testClan2), self.args.num_families, "make_clan returns proper number of lines")
        canonical1 = [['@N5V:1:H5N:1:11103:5:8 1:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AAGATTTGTGTAACTAGGGCTTGTCCTAGCAGGGAA', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHGEFHFIFDGIE', '@N5V:1:H5N:1:11103:5:8 1:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AGAAAACGTCTAACTAGGGCTTGTCCTAGCAGGGAA', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHGEFHFIFDGIE'], ['@N5V:1:H5N:1:11103:5:8 1:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AAGATTTGTGTAACTAGGGCTTGTCCTAGCAGGGAA', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHGEFHFIFDGIE', '@N5V:1:H5N:1:11103:5:8 1:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AGAAAACGTCTAACTAGGGCTTGTCCTAGCAGGGAA', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHGEFHFIFDGIE'], ['@N5V:1:H5N:1:11103:5:8 1:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AAGATTTGTGTAACTAGGGCTTGTCCTAGCAGGGAA', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHGEFHFIFDGIE', '@N5V:1:H5N:1:11103:5:8 1:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AGAAAACGTCTAACTAGGGCTTGTCCTAGCAGGGAA', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHGEFHFIFDGIE']]
        canonical2 = [['@N5V:1:H5N:1:11103:5:8 2:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AGAAAACGTCTAACTAGGGCTTGTCCTAGCAGGGAA', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHGEFHFIFDGIE', '@N5V:1:H5N:1:11103:5:8 2:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AAGATTTGTGTAACTAGGGCTTGTCCTAGCAGGGAA', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHGEFHFIFDGIE'], ['@N5V:1:H5N:1:11103:5:8 2:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AGAAAACGTCTAACTAGGGCTTGTCCTAGCAGGGAA', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHGEFHFIFDGIE', '@N5V:1:H5N:1:11103:5:8 2:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AAGATTTGTGTAACTAGGGCTTGTCCTAGCAGGGAA', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHGEFHFIFDGIE'], ['@N5V:1:H5N:1:11103:5:8 2:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AGAAAACGTCTAACTAGGGCTTGTCCTAGCAGGGAA', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHGEFHFIFDGIE', '@N5V:1:H5N:1:11103:5:8 2:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AAGATTTGTGTAACTAGGGCTTGTCCTAGCAGGGAA', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHGEFHFIFDGIE']]
        self.assertEqual(' '.join(map(" ".join, canonical1)), ' '.join(map(' '.join, testClan1)), "make_clan returns proper values for clan 1.\ncanonical: {0}\n     test: {1} ".format(canonical1, testClan1))
        self.assertEqual(canonical2, testClan2, "make_clan returns proper values for clan 2.\ncanonical: {0}\n     test: {1} ".format(canonical2, testClan2))

    ###########################################################################
    # TEST reads options
    # this should fail because max_num_reads is too small
    def test_min_num_reads_fails_without_max_num_reads(self):
        self.args.quality = self.full_qual
        self.args.read_length = 36
        self.args.barcode_a = self.barcode_a
        self.args.barcode_b = self.barcode_b
        self.args.min_num_reads = 3
        try:
            (testFam1, testFam2) = makeFakeFASTQ.make_family(self.fasta_header, self.seq1, self.args, self.args.min_num_families)
        except ValueError:
            pass
        except Exception as e:
            self.fail("Unexpected exception raised:", e)
        else:
            self.fail("Expected ValueError not raised")

    def test_min_num_reads(self):
        self.args.quality = self.full_qual
        self.args.read_length = 36
        self.args.barcode_a = self.barcode_a
        self.args.barcode_b = self.barcode_b
        self.args.min_num_reads = 3
        self.args.max_num_reads = 3
        (testFam1, testFam2) = makeFakeFASTQ.make_family(self.fasta_header, self.seq1, self.args, self.args.min_num_families)
        canonical1 = ['@N5V:1:H5N:1:11103:5:8 1:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AGAAAACGTCTAACTAGGGCTTGTCCTAGCAGGGAA', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHGEFHFIFDGIE', '@N5V:1:H5N:1:11103:5:8 1:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AAGATTTGTGTAACTAGGGCTTGTCCTAGCAGGGAA', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHGEFHFIFDGIE', '@N5V:1:H5N:1:11103:5:8 1:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AAGATTTGTGTAACTAGGGCTTGTCCTAGCAGGGAA', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHGEFHFIFDGIE', '@N5V:1:H5N:1:11103:5:8 1:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AGAAAACGTCTAACTAGGGCTTGTCCTAGCAGGGAA', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHGEFHFIFDGIE', '@N5V:1:H5N:1:11103:5:8 1:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AAGATTTGTGTAACTAGGGCTTGTCCTAGCAGGGAA', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHGEFHFIFDGIE', '@N5V:1:H5N:1:11103:5:8 1:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AGAAAACGTCTAACTAGGGCTTGTCCTAGCAGGGAA', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHGEFHFIFDGIE']
        canonical2 = ['@N5V:1:H5N:1:11103:5:8 2:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AAGATTTGTGTAACTAGGGCTTGTCCTAGCAGGGAA', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHGEFHFIFDGIE', '@N5V:1:H5N:1:11103:5:8 2:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AGAAAACGTCTAACTAGGGCTTGTCCTAGCAGGGAA', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHGEFHFIFDGIE', '@N5V:1:H5N:1:11103:5:8 2:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AGAAAACGTCTAACTAGGGCTTGTCCTAGCAGGGAA', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHGEFHFIFDGIE', '@N5V:1:H5N:1:11103:5:8 2:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AAGATTTGTGTAACTAGGGCTTGTCCTAGCAGGGAA', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHGEFHFIFDGIE', '@N5V:1:H5N:1:11103:5:8 2:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AGAAAACGTCTAACTAGGGCTTGTCCTAGCAGGGAA', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHGEFHFIFDGIE', '@N5V:1:H5N:1:11103:5:8 2:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AAGATTTGTGTAACTAGGGCTTGTCCTAGCAGGGAA', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHGEFHFIFDGIE']
        self.assertEqual(len(testFam1), 8 * self.args.min_num_reads, "make_family returns proper number of lines testFam1 {0} canonical {1}".format(len(testFam1), 4 * self.args.min_num_reads))
        self.assertEqual(' '.join(canonical1), ' '.join(testFam1), "make_family returns proper values for forward.\ncanonical: {0}\n     test: {1} ".format(canonical1, testFam1))
        self.assertEqual(len(testFam2), 8 * self.args.min_num_reads, "make_family returns proper number of lines")
        self.assertEqual(canonical2, testFam2, "make_family returns proper values for reverse.\ncanonical: {0}\n     test: {1} ".format(canonical2, testFam2))

    def test_min_num_reads_overridden_by_num_reads(self):
        self.args.quality = self.full_qual
        self.args.read_length = 60
        self.args.barcode_a = self.barcode_a
        self.args.barcode_b = self.barcode_b
        self.args.min_num_reads = 3
        self.args.max_num_reads = 3
        self.args.num_reads = 1
        (testFam1, testFam2) = makeFakeFASTQ.make_family(self.fasta_header, self.seq1, self.args, self.args.min_num_families)
        canonical1 = ['@N5V:1:H5N:1:11103:5:8 1:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AAGATTTGTGTAACTAGGGCTTGTCCTAGCAGGGAAAGCTTTGACTTAAAAGCCTCTGCT', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHGEFHFIFDGIE', '@N5V:1:H5N:1:11103:5:8 1:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AGAAAACGTCTAACTAGGGCTTGTCCTAGCAGGGAAAGCTTTGACTTAAAAGCCTCTGCT', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHGEFHFIFDGIE']
        canonical2 = ['@N5V:1:H5N:1:11103:5:8 2:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AGAAAACGTCTAACTAGGGCTTGTCCTAGCAGGGAAAGCTTTGACTTAAAAGCCTCTGCT', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHGEFHFIFDGIE', '@N5V:1:H5N:1:11103:5:8 2:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AAGATTTGTGTAACTAGGGCTTGTCCTAGCAGGGAAAGCTTTGACTTAAAAGCCTCTGCT', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHGEFHFIFDGIE']
        # right now it returns 2 fastq entries per num_read, so 8 total lines per num_read
        self.assertEqual(len(testFam1), 8 * self.args.num_reads, "make_family returns proper number of lines testFam {0} canonical {1}".format(len(testFam1), 8 * self.args.num_reads))
        self.assertEqual(' '.join(canonical1), ' '.join(testFam1), "make_family returns proper values for forward.\ncanonical: {0}\n     test: {1} ".format(canonical1, testFam1))
        self.assertEqual(len(testFam2), 8 * self.args.num_reads, "make_family returns proper number of lines")
        self.assertEqual(canonical2, testFam2, "make_family returns proper values for reverse.\ncanonical: {0}\n     test: {1} ".format(canonical2, testFam2))

    ###########################################################################
    # TEST FREQUENCY FILE OPTIONS

    def test_frequency_file(self):
        self.args.quality = self.full_qual
        self.args.read_length = 36
        self.args.barcode_a = self.barcode_a
        self.args.barcode_b = self.barcode_b
        self.args.frequency = ">test1	3	5	high\n>test2	10	5	high"
        # Expecting 1 here for the self.args.min_num_families
        (testFam1, testFam2) = makeFakeFASTQ.make_family(self.fasta_header, self.seq1, self.args, self.args.min_num_families)
        canonical1 = ['@N5V:1:H5N:1:11103:5:8 1:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AAGATTTGTGTAACTAGGGCTTGTCCTAGCAGGGAA', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHGEFHFIFDGIE', '@N5V:1:H5N:1:11103:5:8 1:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AGAAAACGTCTAACTAGGGCTTGTCCTAGCAGGGAA', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHGEFHFIFDGIE']
        canonical2 = ['@N5V:1:H5N:1:11103:5:8 2:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AGAAAACGTCTAACTAGGGCTTGTCCTAGCAGGGAA', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHGEFHFIFDGIE', '@N5V:1:H5N:1:11103:5:8 2:N:0:test1:AGAAAACGTCAAGATTTGTG', 'AAGATTTGTGTAACTAGGGCTTGTCCTAGCAGGGAA', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHGEFHFIFDGIE']
        self.assertEqual(' '.join(canonical1), ' '.join(testFam1), "make_family returns proper values for family 1.\ncanonical: {0}\n     test: {1} ".format(canonical1, testFam1))
        self.assertEqual(canonical2, testFam2, "make_family returns proper values for family 2.\ncanonical: {0}\n     test: {1} ".format(canonical2, testFam2))

    ###########################################################################

if __name__ == '__main__':
    unittest.main()
