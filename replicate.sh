set -x 
set -e

for i in $(seq 1 10)
do
  python deblur/test/test_workflow.py workflowTests.test_dereplicate_seqs
done
