from unittest import TestCase, main
from os.path import join, abspath, dirname

from biom import load_table
from deblur.workflow import _system_call
from deblur.workflow import sequence_generator
from tempfile import mkdtemp


class TestScript(TestCase):

    def setUp(self):
        self.data_dir = join(dirname(abspath(__file__)), 'data')
        self.seqs_fp = join(self.data_dir, 'seqs_s3.fasta')
        self.orig_fp = join(self.data_dir, 'simset.s3.fasta')
        self.output_dir = mkdtemp()
        self.output_biom = join(self.output_dir, 'final.biom')
        self.output_seqs = join(self.output_dir, 'final.seqs.fa')
        self.trim_length = 150

    def test_workflow(self):
        cmd = ["deblur", "workflow", "--seqs-fp", self.seqs_fp,
               "--output-dir", self.output_dir,
               "--trim-length", "150", '-w', '-n']
        sout, serr, res = _system_call(cmd)
        res_table = load_table(self.output_biom)

        res_seqs = list(res_table.ids(axis='observation'))
        exp_seqs = [item[1] for item in sequence_generator(self.orig_fp)]
        exp_seqs = list(map(lambda x: x.upper()[:self.trim_length], exp_seqs))
        self.assertListEqual(res_seqs, exp_seqs)


if __name__ == "__main__":
    main()
