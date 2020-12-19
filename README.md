# rMSA #
RNA Multiple Sequence Alignment generation

## Install ##
```bash
git clone https://github.com/kad-ecoli/rMSA
cd rMSA
./database/script/update.sh    # Download RNAcentral and nt
```
rMSA is for 64bit Linux (x86-64) only.

## Run the program ##
```bash
./rMSA.pl seq.fasta
```
Run ``./rMSA.pl`` without command line argument to get full option list,
including alternative databases, temporary folders, number of CPU threads
(default is 1) and the secondary structure. Nucleotide U is converted to T.

## Third party programs ##
* nhmmer from HMMER 3.3
* cd-hit-est and cd-hit-est-2d from CD-HIT 4.8.1
* clustalo 1.2.4
* cmbuild, cmcalibrate, cmscan and cmsearch from INFERNAL 1.1.3
* hhfilter and reformat.pl from HH-suite 2.0.15
* RNAfold from ViennaRNA 2.4.14

The output format of nhmmer and cmsearch are modifed from
``eslMSAFILE_STOCKHOLM`` or ``eslMSAFILE_PFAM`` to ``eslMSAFILE_A2M``.

## License ##
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 2 of the License, or
(at your option) any later version.
