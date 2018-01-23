from unittest import TestCase, main
from os.path import join, abspath, dirname

from deblur.workflow import _system_call
from tempfile import mkdtemp
from shutil import rmtree


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


if __name__ == "__main__":
    main()
