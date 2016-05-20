import unittest
import makeFakeFASTQ
import re
from argparse import Namespace

def avg_qual(seq):
    qual_array = map(lambda x: ord(x) - 33, list(seq))
    avg = (sum(map(int, qual_array)))/len(qual_array)
    return avg

class FakeFASTQTest(unittest.TestCase):
    def setUp(self):
        self.fasta_header="test1"
        self.barcode = "ATATAGGACA"
        self.seq =  "GCTAATACGAATTGACCCGGATAGA"
        self.qual = "HFDEGEFHFIFDGIEHIFGEJJJHG"
        self.args = Namespace(x_min=5,x_max=5,y_min=8,y_max=8,
                              lane_min=1,lane_max=1,swathes=[111],
                              tile_min=3,tile_max=3,is_filtered=['N'],
                              include_fasta_header_in_fastq_header=1,include_barcode_in_fastq_header=1,
                              paired_end=0,instrument="N5V",flow_cell="H5N",quality_type="high",
                              read_length=25, spacer_length=1, max_num_reads=1, max_num_families=1,
                              map_file=0)

    #def tearDown(self):
    #    print("teardown")

    #TEST def random_sequence
    def test_random_sequence_length(self):
        #print("Testing that random sequence is specified length...")
        self.assertEqual(len(makeFakeFASTQ.random_sequence(10)), 10, "10 length random \
        sequence")

    def test_random_sequence_composition(self):
        #print("Testing that random sequence contains only AGTC...")
        reg = re.compile('^[AGTC]+$')
        seq = makeFakeFASTQ.random_sequence(10)
        self.assertTrue(reg.match(seq),"random_sequence contains only AGTC")

    #TEST def fastq_quality(args, seq_len)
    # requires args.quality_type where quality_type is high, medium, low
    # returns a string of seq_len
    def test_fastq_quality_length(self):
        #print("Testing that fastq_quality is specified length...")
        self.assertEqual(len(makeFakeFASTQ.fastq_quality(self.args,10)), 10, "10 length quality")

    def test_fastq_quality_high(self):
        #print("Testing that fastq_quality high works...")
        qual = makeFakeFASTQ.fastq_quality(self.args,10)
        avg = avg_qual(qual)
        self.assertGreaterEqual(avg, 35, "qual quality {} is not greater than or equal to 35".format(avg))

    def test_fastq_quality_medium(self):
        #print("Testing that fastq_quality medium works...")
        args = Namespace(quality_type='medium')
        qual = makeFakeFASTQ.fastq_quality(args,10)
        avg = avg_qual(qual)
        self.assertLessEqual(avg, 30, "qual quality average {} is less than\
                30".format(avg))
        self.assertGreaterEqual(avg, 23, "qual quality average {} is not \
                greater than or equal to 23".format(avg))

    def test_fastq_quality_low(self):
        #print("Testing that fastq_quality low works...")
        args = Namespace(quality_type='low')
        qual = makeFakeFASTQ.fastq_quality(args,10)
        avg = avg_qual(qual)
        self.assertLessEqual(avg, 15, "qual quality average {} is less than\
                15".format(avg))
        self.assertGreaterEqual(avg, 3, "qual quality average {} is not \
                greater than or equal to 3".format(avg))

    def test_fastq_high_quality_letters(self):
        #print("Testing that fastq_quality result contains proper ascii letters...")
        qual = makeFakeFASTQ.fastq_quality(self.args,10)
        reg = re.compile('^[A-J]+$')
        self.assertTrue(reg.match(qual),"fastq_quality high contains only ascii letters: {}".format(qual))

    def test_fastq_low_quality_valid(self):
        #print("Testing that fastq_quality result contains valid values...")
        # valid  !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHI
        args = Namespace(quality_type="low")
        qual = makeFakeFASTQ.fastq_quality(args,10)
        reg = re.compile('^[A-I0-9!\"#\$\%\&\'\(\)\*\+,=\./:;<=>?@-]+$')
        self.assertTrue(reg.match(qual),"fastq_quality low contains only valid values:{}".format(qual))


    ### TEST def fastq_entry_header(args, fasta_header, barcode)
    def test_fastq_entry_header(self):
        # test header that doesn't include barcode or fasta header
        #print("Testing fastq entry header...")
        self.args.include_barcode_in_fastq_header = 0
        self.args.include_fasta_header_in_fastq_header = 0
        testHeader, testPaired = makeFakeFASTQ.fastq_entry_header(self.args,self.fasta_header,self.barcode)
        canonical = "@N5V:1:H5N:1:11103:5:8 1:N:0"
        self.assertEqual(canonical, testHeader, "fastq_entry_header returns proper header")

    def test_fastq_entry_header_fasta_header(self):
        self.args.include_barcode_in_fastq_header = 0
        testHeader, testPaired = makeFakeFASTQ.fastq_entry_header(self.args,self.fasta_header,self.barcode)
        canonical = "@N5V:1:H5N:1:11103:5:8 1:N:0:test1"
        self.assertEqual(canonical, testHeader, "fastq_entry_header returns proper header with fasta header no > {}".format(testHeader))

    def test_fastq_entry_header_fasta_header_symbol(self):
        self.fasta_header = ">test1"
        self.args.include_barcode_in_fastq_header = 0
        testHeader, testPaired = makeFakeFASTQ.fastq_entry_header(self.args,self.fasta_header,self.barcode)
        canonical = "@N5V:1:H5N:1:11103:5:8 1:N:0:test1"
        self.assertEqual(canonical, testHeader, "fastq_entry_header returns proper header with fasta header >\ncanonical: {}\ntest:      {}".format(canonical, testHeader))

    def test_fastq_entry_header_barcode(self):
        self.args.include_fasta_header_in_fastq_header = 0
        testHeader, testPaired = makeFakeFASTQ.fastq_entry_header(self.args,self.fasta_header,self.barcode)
        canonical = "@N5V:1:H5N:1:11103:5:8 1:N:0:ATATAGGACA"
        self.assertEqual(canonical, testHeader, "fastq_entry_header returns proper header with barcode canonical: {} test:".format(canonical,testHeader))

    def test_fastq_entry_header_fasta_header_and_barcode(self):
        testHeader, testPaired = makeFakeFASTQ.fastq_entry_header(self.args,self.fasta_header,self.barcode)
        canonical = "@N5V:1:H5N:1:11103:5:8 1:N:0:test1:ATATAGGACA"
        self.assertEqual(canonical, testHeader, "fastq_entry_header returns proper header with fasta header and barcode {}".format(testHeader))

    def test_fastq_entry_header_fasta_header_symbol_and_barcode(self):
        self.fasta_header = ">test1"
        testHeader,testPaired = makeFakeFASTQ.fastq_entry_header(self.args,self.fasta_header,self.barcode)
        canonical = "@N5V:1:H5N:1:11103:5:8 1:N:0:test1:ATATAGGACA"
        self.assertEqual(canonical, testHeader, "fastq_entry_header returns proper header with fasta header with > and barcode {}".format(testHeader))

    # TEST def buffer_sequence(args,seq)
    def test_buffer_sequence_front(self):
        testSeq = makeFakeFASTQ.buffer_sequence(self.args,self.seq)
        self.assertEqual(len(testSeq),self.args.read_length,"buffer sequence outputs correct read length args: {} test:{}".format(self.args.read_length, len(testSeq)))

    #TEST def make_family(header, seq, args)
    def test_make_family(self):
        self.args.quality = self.qual
        self.args.barcode = self.barcode
        (testFam1, testFam2) = makeFakeFASTQ.make_family(self.fasta_header, self.seq, self.args)
        print("testFam1 {} testFam2 {}".format(testFam1, testFam2))
        canonical1 = ['@N5V:1:H5N:1:11103:5:8 1:N:0:test1:ATATAGGACA', 'ATATAGGACATGCTAATACGAATTGACCCGGATAGATATATAGGACA', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHG']
        canonical2 = ['@N5V:1:H5N:1:11103:5:8 2:N:0:test1:CAATAACTAA', 'CAATAACTAATGCTAATACGAATTGACCCGGATAGATCAATAACTAA', '+', 'HFDEGEFHFIFDGIEHIFGEJJJHG']
        self.assertEqual(' '.join(canonical1), ' '.join(testFam1), "make_family returns proper values for family 1.\ncanonical: {}\n     test: {} ".format(canonical1, testFam1))
        self.assertEqual(canonical2, testFam2, "make_family returns proper values for family 2.\ncanonical: {}\n     test: {} ".format(canonical2, testFam2))




if __name__ == '__main__':
    unittest.main()
