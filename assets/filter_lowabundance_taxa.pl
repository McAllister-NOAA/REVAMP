#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use List::MoreUtils qw(uniq);
use List::Util qw( min max );


# - - - - - H E A D E R - - - - - - - - - - - - - - - - -
#Goals of script:
#Import ASVs_counts_NOUNKNOWNS_collapsedOnTaxonomy_percentabund.tsv & x_asvTaxonomyTable_NOUNKNOWNS.txt
#Output: Modified ASV taxonomy table that replaces taxonomy positions <x% with zzOther. Reassess at each level.

# - - - - - C O M M A N D    L I N E    O P T I O N S - - - - - - - -
my %options=();
getopts("a:t:p:h", \%options);

if ($options{h})
    {   print "\n\nHelp called:\nOptions:\n";
        print "-a = ASVs_counts_NOUNKNOWNS_collapsedOnTaxonomy_percentabund.tsv file\n";
        print "-t = ASV taxonomy table from processing script (x_asvTaxonomyTable_NOUNKNOWNS.txt)\n";
        print "-p = Filter percent (i.e. 5)\n";
        print "-h = This help message\n\n";
        die;
    }

# - - - - - G L O B A L  V A R I A B L E S  - - - - - -
my %ASV;
my %TAX;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - M A I N - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#IMPORT ASV Taxonomy table
open(IN2, "<$options{t}") or die "\n\nThere is no $options{t} file!!\n\n";
my @tax_dat = <IN2>; close(IN2);
my @uselater = @tax_dat;
shift(@tax_dat);
foreach my $line (@tax_dat)
	{	chomp($line);
		my @data = split('\t', $line);
		my $asv = $data[0]; chomp($asv);
        $TAX{"k_".$data[1]}{'asvs'} .= $asv.";";
        $TAX{"p_".$data[2]}{'asvs'} .= $asv.";";
        $TAX{"c_".$data[3]}{'asvs'} .= $asv.";";
        $TAX{"o_".$data[4]}{'asvs'} .= $asv.";";
        $TAX{"f_".$data[5]}{'asvs'} .= $asv.";";
        $TAX{"g_".$data[6]}{'asvs'} .= $asv.";";
        $TAX{"s_".$data[7]}{'asvs'} .= $asv.";";
	}
    
#IMPORT ASV count table info
open(IN1, "<$options{a}") or die "\n\nThere is no $options{a} file!!\n\n";
my @asv_dat = <IN1>; close(IN1);
my $sample_header = shift(@asv_dat);
my @file_sample_headers = split('\t', $sample_header);
foreach my $line (@asv_dat)
	{	chomp($line);
		my @data = split('\t', $line);
		my $asv = $data[0]; chomp($asv);
		foreach my $i (1..$#data)
			{	my $number = $data[$i]; chomp($number);
                my $samplehead = $file_sample_headers[$i]; chomp($samplehead);
				$ASV{$asv}{$samplehead} = $number; #Stores ASV count by sample
			}
	}

#Not all ASVs from taxonomy table are in ASV count table.

#Process Taxonomy dictionary for percent filter
#In ANY sample, is the sum of that taxa > filter percent
foreach my $i (sort keys %TAX)
    {   if ($i =~ m/^[kpcofgs]_NA$/)
            {   $TAX{$i}{'display'} = "NA";
            }
        else {
            my $goodtax = "FALSE";
            my $asvs = $TAX{$i}{'asvs'};
            my @asv_list = split(';', $asvs);
            foreach my $sample (1..$#file_sample_headers)
                {   my $samplehead = $file_sample_headers[$sample]; chomp($samplehead);
                    my $relabundsum = 0;
                    foreach my $j (@asv_list)
                        {   if (exists $ASV{$j})
                                {   $relabundsum += $ASV{$j}{$samplehead};
                                }
                        }
                    if ($relabundsum >= $options{p})
                        {   $goodtax = "TRUE";
                        }
                }
            if ($goodtax eq "TRUE")
                {   $TAX{$i}{'display'} = "TRUE";
                }
            else
                {   $TAX{$i}{'display'} = "FALSE";
                }
        }
    }


#Print out new taxonomy table.
my $header = shift(@uselater); chomp($header);
print "$header\n";

foreach my $line (@uselater)
    {   chomp($line);
        my @data = split('\t', $line);
        print "$data[0]\t";
        if ($TAX{"k_".$data[1]}{'display'} eq "TRUE" || $TAX{"k_".$data[1]}{'display'} eq "NA")
            {   print "$data[1]\t";
            }
        if ($TAX{"k_".$data[1]}{'display'} eq "FALSE")
            {   print "zzOther\t"; 
            }
        if ($TAX{"p_".$data[2]}{'display'} eq "TRUE" || $TAX{"p_".$data[2]}{'display'} eq "NA")
            {   print "$data[2]\t";
            }
        if ($TAX{"p_".$data[2]}{'display'} eq "FALSE")
            {   print "zzOther\t"; 
            }
        if ($TAX{"c_".$data[3]}{'display'} eq "TRUE" || $TAX{"c_".$data[3]}{'display'} eq "NA")
            {   print "$data[3]\t";
            }
        if ($TAX{"c_".$data[3]}{'display'} eq "FALSE")
            {   print "zzOther\t"; 
            }
        if ($TAX{"o_".$data[4]}{'display'} eq "TRUE" || $TAX{"o_".$data[4]}{'display'} eq "NA")
            {   print "$data[4]\t";
            }
        if ($TAX{"o_".$data[4]}{'display'} eq "FALSE")
            {   print "zzOther\t"; 
            }
        if ($TAX{"f_".$data[5]}{'display'} eq "TRUE" || $TAX{"f_".$data[5]}{'display'} eq "NA")
            {   print "$data[5]\t";
            }
        if ($TAX{"f_".$data[5]}{'display'} eq "FALSE")
            {   print "zzOther\t"; 
            }
        if ($TAX{"g_".$data[6]}{'display'} eq "TRUE" || $TAX{"g_".$data[6]}{'display'} eq "NA")
            {   print "$data[6]\t";
            }
        if ($TAX{"g_".$data[6]}{'display'} eq "FALSE")
            {   print "zzOther\t"; 
            }
        if ($TAX{"s_".$data[7]}{'display'} eq "TRUE" || $TAX{"s_".$data[7]}{'display'} eq "NA")
            {   print "$data[7]";
            }
        if ($TAX{"s_".$data[7]}{'display'} eq "FALSE")
            {   print "zzOther"; 
            }
        print "\n";
    }
    
##END OF RUN##