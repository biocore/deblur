from unittest import TestCase, main
from os.path import join, abspath, dirname

from deblur.workflow import _system_call, create_otu_table
from tempfile import mkdtemp
from shutil import rmtree
from biom import load_table


class TestScript(TestCase):
    def setUp(self):
        self.data_dir = join(dirname(abspath(__file__)), 'data')
        self.output_dir = mkdtemp()

    def tearDown(self):
        rmtree(self.output_dir)

    def test_workflow(self):
        cmd = ["deblur", "workflow",
               "--seqs-fp", join(self.data_dir, 'usecase_mixedchars',
                                 'split'),
               "--output-dir", self.output_dir,
               "--error-dist", ("1, 0.06, 0.02, 0.02, 0.01, 0.005, 0.005, "
                                "0.005, 0.001, 0.001, 0.001, 0.0005"),
               '--indel-max', "3",
               '--indel-prob', "0.01",
               '--jobs-to-start', "1",
               '--mean-error', "0.005",
               '--min-reads', "0",
               '--min-size', "2",
               '--threads-per-sample', "1",
               '--trim-length', "100",
               '-w',
               ]
        # run deblur
        sout, serr, res = _system_call(cmd)

        # test that all feature sequences are all upper case
        with open(join(self.output_dir, 'all.seqs.fa'), 'r') as f:
            for line in f.readlines():
                if line.startswith('>'):
                    self.assertEqual(line[1:].strip(),
                                     line[1:].strip().upper())

    def test_create_otu_table(self):
        fp_table = join(self.output_dir, 'all.biom')
        create_otu_table(
            fp_table,
            [(join(self.data_dir, 'usecase_mixedchars',
                   ('1.SKB7.640196.fastq.trim.derep.'
                    'no_artifacts.msa.deblur.no_chimeras')), '1.SKB7.640196'),
             (join(self.data_dir, 'usecase_mixedchars',
                   ('1.SKB8.640193.fastq.trim.derep.'
                    'no_artifacts.deblur.no_chimeras')), '1.SKB8.640193')],
            outputfasta_fp=join(self.output_dir, 'all.seqs'), minreads=0)
        table = load_table(fp_table)

        # should produce a table with two samples and two features
        self.assertEqual(table.shape, (2, 2))

        # assert that counts from different case entries are collapsed
        self.assertTrue(list(table.to_dataframe().to_dense().loc[
            ('TACGGGGGGGGTTAGCGTTATTCAATGATATTTGGCGTAAAGTGCATGTAGATGGTGTTAC'
             'AAGTTAAAAAAATAAAAACTAAGGACAAATCTTTTCGTT'), :].values) == [60, 0])


if __name__ == "__main__":
    main()
