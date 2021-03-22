## Instructions
Python INC was developed and tested using `conda` environments.
There are three dependencies
```
conda install -c bioconda dendropy
conda install -c anaconda scikit-bio
conda install -c anaconda numpy
```

The command to convert an alignment (FASTA) to PHYLIP distance matrix is
```
python msa.py --in-fasta <fasta_path> --out <phylip_path> --os <ubuntu_or_osx> [--model <model>]
```
`logDet` is the default model, but `K2P`, `JC`, or `P` are acceptable. This tool depends on the (provided) paup* binary, so be sure to run the command from the `src/python` directory.

The command to create newick trees from a PHYLIP file, using INC and/or INC_NJ is
```
python inc.py --in-phylip <phylip_path> --out <newick_path> [--inc] [--inc-nj]
```
Consult the `commands.sh` file for examples.
For best results, specify the full paths.

Various other changes (for example to the voting scheme), are available in the code, but not accessible from the command line. Feel free to make changes to the code for your needs or exploration.

## Contact
For questions re: the python implementation, email raaronsy@berkeley.edu
