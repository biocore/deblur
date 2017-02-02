from unittest import TestCase, main
from os.path import join, abspath, dirname

from biom import load_table
from deblur.workflow import _system_call
from deblur.workflow import sequence_generator
from tempfile import mkdtemp


class TestScript(TestCase):

    def setUp(self):
        self.data_dir = join(dirname(abspath(__file__)), 'data')
        self.seqs_fp = join(self.data_dir, 'seqs_s5.fasta')
        self.orig_fp = join(self.data_dir, 'simset.s3.fasta')
        self.orig_one_seq_fp = join(self.data_dir, 'simset.s5.fasta')
        self.output_dir = mkdtemp()
        self.output_biom = join(self.output_dir, 'all.biom')
        self.output_seqs = join(self.output_dir, 'all.seqs.fa')
        self.trim_length = 150

    def validate_results(self, table_name, orig_fasta_name):
        res_table = load_table(table_name)

        res_seqs = list(res_table.ids(axis='observation'))
        exp_seqs = [item[1] for item in sequence_generator(orig_fasta_name)]
        exp_seqs = list(map(lambda x: x.upper()[:self.trim_length], exp_seqs))
        self.assertListEqual(res_seqs, exp_seqs)

    def test_workflow(self):
        # test default parameters, negative mode, single thread
        cmd = ["deblur", "workflow", "--seqs-fp", self.seqs_fp,
               "--output-dir", self.output_dir,
               "--trim-length", "150", '-w']
        sout, serr, res = _system_call(cmd)
        self.validate_results(self.output_biom, self.orig_one_seq_fp)

        # test default parameters, negative mode, multi thread
        cmd = ["deblur", "workflow", "--seqs-fp", self.seqs_fp,
               "--output-dir", self.output_dir,
               "--trim-length", "150", '-w',
               "-O", "2"]
        sout, serr, res = _system_call(cmd)
        self.validate_results(self.output_biom, self.orig_one_seq_fp)

        # test default parameters, positive mode, single thread
        cmd = ["deblur", "workflow", "--seqs-fp", self.seqs_fp,
               "--output-dir", self.output_dir,
               "--trim-length", "150", '-w']
        sout, serr, res = _system_call(cmd)
        self.validate_results(self.output_biom, self.orig_one_seq_fp)

        # test default parameters, positive mode, multi thread
        cmd = ["deblur", "workflow", "--seqs-fp", self.seqs_fp,
               "--output-dir", self.output_dir,
               "--trim-length", "150", '-w',
               "-O", "2"]
        sout, serr, res = _system_call(cmd)
        self.validate_results(self.output_biom, self.orig_one_seq_fp)

        # test default parameters except min-reads set to 0, negative mode, single thread
        cmd = ["deblur", "workflow", "--seqs-fp", self.seqs_fp,
               "--output-dir", self.output_dir,
               "--trim-length", "150", '-w', '--min-reads', '0']
        sout, serr, res = _system_call(cmd)
        self.validate_results(self.output_biom, self.orig_fp)

if __name__ == "__main__":
    main()
