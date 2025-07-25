#!/bin/bash

#Run from inside the folder where your nt database is held.
#Run to update nt, taxonomy dmp, and prep files for REVAMP.
#Requires blast+, wget, and taxonomy-tools (https://github.com/pmenzel/taxonomy-tools), all in PATH
echo "Date of last update to the nt database:" >> nt_lastupdate.txt
date >> nt_lastupdate.txt
echo >> nt_lastupdate.txt
echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" >> nt_lastupdate.txt
echo >> nt_lastupdate.txt

echo
echo "Updating nt BLAST database..."
update_blastdb.pl nt
for f in *.tar.gz; do tar -xzf $f; rm $f; done

echo
echo "Updating taxonomy database..."
rm -r taxdump
mkdir taxdump
cd taxdump
wget https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xzf taxdump.tar.gz
rm taxdump.tar.gz

echo
echo "Preparing files for REVAMP..."
grep "genbank common name" names.dmp > common_names.dmp

grep "scientific name" names.dmp | grep -i "uncultured" | cut -f1 > taxids_to_ignore.txt
grep "scientific name" names.dmp | grep -i "environmental samples" | cut -f1 >> taxids_to_ignore.txt
grep "scientific name" names.dmp | grep -i "metagenome" | cut -f1 >> taxids_to_ignore.txt
grep "scientific name" names.dmp | grep -i "unidentified" | cut -f1 >> taxids_to_ignore.txt
subtree -t nodes.dmp -i taxids_to_ignore.txt > taxid_exclusion_list_leavesinUnclassified.txt
grep "scientific name" names.dmp | grep -i "unclassified" | cut -f1 >> taxids_to_ignore.txt
subtree -t nodes.dmp -i taxids_to_ignore.txt > taxid_exclusion_list_removesUnclassified.txt
rm taxids_to_ignore.txt

#Removed due to new REVAMP use of the $TAXONKIT_DB env variable
#mkdir -p ~/.taxonkit
#cp names.dmp ~/.taxonkit/
#cp nodes.dmp ~/.taxonkit/
#cp delnodes.dmp ~/.taxonkit/
#cp merged.dmp ~/.taxonkit/
