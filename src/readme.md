C++ utilities for parsing MSA and RNA secondary structure
```bash
a3m2msa     # convert a3m format MSA to fasta MSA without insertion states
fasta2pfam  # convert fasta to tab-eliminated table
fastaNA     # clean non-standard nucleotide in fasta
fastNf      # calculate length normalized number of effective sequence (Nf)
fixAlnX     # remove unknown residue type from MSA
pfam2fasta  # convert the output of fasta2pfam back to fasta
rFUpred     # FUpred domain partition algorithm for RNA secondary structure
RemoveNonQueryPosition # delete any position corresponding to gap in query
trimblastN  # trim sequence hits
```

Install the programs by
```bash
make
make install
```
