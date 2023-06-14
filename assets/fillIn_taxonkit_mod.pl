#!/usr/bin/perl -w
use strict;
use Getopt::Std;

# - - - - - H E A D E R - - - - - - - - - - - - - - - - -
#Goals of script:
#Fill in blanks in taxonomic hierarchy up from Genus if present or between hierarchy otherwise.
#Add "Genus unknown__s" if genus present and species absent (rare).
#Add "unknown__s" if species and genus absent (rarer).


# - - - - - C O M M A N D    L I N E    O P T I O N S - - - - - - - -
my %options=();
getopts("i:h", \%options);

if ($options{h})
    {   print "\n\nHelp called:\nOptions:\n";
        print "-i = Input K/P/C/O/F/G/S file from taxonkit (reformatted_taxonkit_out.txt)\n";
        print "-h = This help message\n\n";
        die;
    }

# - - - - - G L O B A L  V A R I A B L E S  - - - - - -


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - M A I N - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
open(IN, "<$options{i}") or die "\n\nNADA $options{i} file!!\n\n";
my @data = <IN>; close(IN);

foreach my $line (@data)
    {   chomp($line);
        my $kingdom = "";
        my $phylum = "";
        my $class = "";
        my $order = "";
        my $family = "";
        my $genus = "";
        my $species = "";
        my @split_line = split('\t', $line);
        my $taxonID = $split_line[0];
        my $originaltaxonkitlineage = $split_line[1];
        my $taxonomyString = $split_line[2];
        my @splittaxonomy = split(';', $taxonomyString);
        if (exists $splittaxonomy[0]) {$kingdom = $splittaxonomy[0];}
        if (exists $splittaxonomy[1]) {$phylum = $splittaxonomy[1];}
        if (exists $splittaxonomy[2]) {$class = $splittaxonomy[2];}
        if (exists $splittaxonomy[3]) {$order = $splittaxonomy[3];}
        if (exists $splittaxonomy[4]) {$family = $splittaxonomy[4];}
        if (exists $splittaxonomy[5]) {$genus = $splittaxonomy[5];}
        if (exists $splittaxonomy[6]) {$species = $splittaxonomy[6];}
        if ($#splittaxonomy > 6)
            {   die "\n\nTaxonkit format is not correct\n\n";}
        if ($kingdom eq "")
            {   print "$line\n";}
        else
            {   if ($genus ne "")
                    {   if ($family eq "")
                            {   $family = $genus."__f";}
                        if ($order eq "")
                            {   $order = $family."__o";}
                        if ($class eq "")
                            {   $class = $order."__c";} 
                        if ($phylum eq "")
                            {   $phylum = $class."__p";}   
                    }
                if ($genus eq "")
                    {   if ($family ne "")
                            {   if ($order eq "")
                                    {   $order = $family."__o";}
                                if ($class eq "")
                                    {   $class = $order."__c";}
                                if ($phylum eq "")
                                    {   $phylum = $class."__p";}
                            }
                        if ($order ne "")
                            {   if ($class eq "")
                                    {   $class = $order."__c";}
                                if ($phylum eq "")
                                    {   $phylum = $class."__p";}
                            }
                        if ($class ne "")
                            {   if ($phylum eq "")
                                    {   $phylum = $class."__p";}
                            }
                    }
                my $mergetaxastring = $kingdom.";".$phylum.";".$class.";".$order.";".$family.";".$genus.";".$species;
                print "$taxonID\t$originaltaxonkitlineage\t$mergetaxastring\n";
            }
    }
    
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - S U B R O U T I N E S - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -
