70_otus.fasta - the 70% representative set from greengenes 13.8

simset.s1.fasta - 10 sequences from greengenes 13.8, selected randomly with >0.5 similarity to seed sequence
simset.s2.fasta - 10 sequences from greengenes 13.8, selected randomly with >0.9 similarity to seed sequence

seqs_s1.fa - a simulated dataset with 100 reads derived from simset.s1.fasta using art:

# simulate reads:
art_illumina -i simset.s1.fasta -l 150 -f 10 -amp -o simset.s1.artreads
# convert to fasta:
fastq_to_fasta.sh simset.s1.artreads.fq > simset.s1.artreads.fa
# add s1_ to every read fasta header (since we want a post split-libraries fasta)
~/scripts/AddStringToFastaHeader.py -i simset.s1.artreads.fa > seqs_s1.fa


seqs_s2.fa - a simulated dataset with 200 reads derived from simset.s2.fasta using grinder, with 10% chimeras, and additionally 35 phix reads:

# simulate reads and chimeras
grinder -rf simset.s2.fasta -rd 151 -un 1 -cp 10 -random_seed 100 -mutation_dist "poly4 3e-3 3.3e-8" -tr 200
mv grinder-reads.fa seqs_s2.grinder.fa
# simulate phix reads
~/bin/art_bin_VanillaIceCream/art_illumina -i ~/data/phix/PhiX.fasta -l 150 -o phix_sim2 -f 1
cp seqs_s2.grinder.fa seqs_s2.withphix.fa
cat phix_sim2.fa >> seqs_s2.withphix.fa
# add s1_ to every read fasta header (since we want a post split-libraries fasta)
~/scripts/AddStringToFastaHeader.py -i seqs_s2.withphix.fa > seqs_s2.fa



