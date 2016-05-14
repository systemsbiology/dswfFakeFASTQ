import unittest
import makeFakeFASTQ
import re
from argparse import Namespace

class FakeFASTQTest(unittest.TestCase):
    def setUp(self):
        print("setup")

    def tearDown(self):
        print("teardown")

    #TEST def random_sequence
    def test_random_sequence_length(self):
        print "Testing that random sequence is specified length..."
        self.assertEqual(len(makeFakeFASTQ.random_sequence(10)), 10, "10 length random \
        sequence")

    def test_random_sequence_composition(self):
        print "Testing that random sequence contains only AGTC..."
        reg = re.compile('^[AGTC]+$')
        seq = makeFakeFASTQ.random_sequence(10)
        self.assertTrue(reg.match(seq),"random_sequence contains only AGTC")

    #TEST def fastq_quality(args, seq_len)
    # requires args.quality_type where quality_type is high, medium, low
    # returns a string of seq_len
    def test_fastq_quality_length(self):
        print "Testing that fastq_quality is specified length..."
        args = Namespace(quality_type='high')
        test = makeFakeFASTQ.fastq_quality(args,10)
        self.assertEqual(len(makeFakeFASTQ.fastq_quality(args,10)), 10, "10 length quality")

    def test_fastq_quality_high(self):
        print "Testing that fastq_quality high works..."

    def test_fastq_quality_letters(self):
        print "Testing that fastq_quality result contains proper ascii \
        letters..."

if __name__ == '__main__':
    unittest.main()
