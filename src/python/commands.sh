# convert the aligned sequences to logDet distances
#   You will have to change the paths for your setup
#   The script writes each PHYLIP distance matrix into the same folder as the fasta file
python msa.py --in-fasta "/media/sf_constraint_inc/data/100_seq/R0/rose.aln.true.fasta" --out "/media/sf_constraint_inc/data/100_seq/R0/logDet.dist" --os "ubuntu" --model "logDet" # --os "osx"

# You can also generate the newick trees inferred by INC or INC_NJ themselves.
python inc.py --in-phylip "/media/sf_constraint_inc/data/100_seq/R0/logDet.dist" --out "/media/sf_constraint_inc/data/100_seq/R0/inc_trees.txt" --inc --inc-nj