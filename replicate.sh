set -x 
set -e

#conda create --yes -n deblur-hammer python=3.5 pip nose flake8 h5py libgcc # "numpy<1.14.0,>=1.9.2"
# source activate deblur-hammer
#conda install --yes -c bioconda "VSEARCH=2.8.3"
#conda install --yes -c bioconda MAFFT=7.310
#conda install --yes -c bioconda SortMeRNA=2.0
#pip install -U pip coveralls
#pip install --process-dependency-links .
#vsearch --version
#ldd $(which vsearch)

for i in $(seq 1 10)
do
  python deblur/test/test_workflow.py workflowTests.test_dereplicate_seqs
done
