#!/bin/bash
FILE=`readlink -e $0`
bindir=`dirname $FILE`
rootdir=`dirname $bindir`
cd $rootdir

echo "download RNAcentral"
wget ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/sequences/rnacentral_species_specific_ids.fasta.gz -O rnacentral_species_specific_ids.fasta.gz 
wget ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/release_notes.txt -O release_notes.txt
wget ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/rfam/rfam_annotations.tsv.gz -O rfam_annotations.tsv.gz

echo "download rfam"
wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.full_region.gz -O Rfam.full_region.gz 
wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz -O Rfam.cm.gz 
wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/database_files/version.txt.gz -O version.txt.gz
gzip -d -f version.txt.gz
