#!/usr/bin/perl -w
use strict;
use Getopt::Std;

# - - - - - H E A D E R - - - - - - - - - - - - - - - - -
#Goals of script:
#Convert ASV count data to presence absence table by replicate.
#Format is for box plots in R.
#
#Inputs:
# sample_metadata_forR.txt
# ASV counts table: one of:
#   ASV counts with unknowns, but with controls/low qual samples removed
#   ASVs_counts_NOUNKNOWNS_controlsRemoved.tsv or with low qual samples filtered
#   ASVs_counts_NOUNKNOWNS_collapsedOnTaxonomy_controlsRemoved or with low qual samples filtered

# - - - - - C O M M A N D    L I N E    O P T I O N S - - - - - - - -
my %options=();
getopts("i:m:h", \%options);

if ($options{h})
    {   print "\n\nHelp called:\nOptions:\n";
        print "-i = Input ASV count file\n";
        print "-m = Sample metadata file (for R)\n";
        print "-h = This help message\n\n";
        die;
    }

# - - - - - G L O B A L  V A R I A B L E S  - - - - - -
my %ASVs;
my %META;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - M A I N - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
open(IN, "<$options{i}") or die "\n\nNADA $options{i} file!!\n\n";
my @data = <IN>; close(IN);
my $sampleheaders = shift(@data); chomp($sampleheaders);
my @sampleheads = split('\t', $sampleheaders);

foreach my $line (@data)
    {   chomp($line);
        my @split_line = split('\t', $line);
        my $asv = $split_line[0];
        foreach my $i (1..$#split_line)
            {   my $sample = $sampleheads[$i];
                $ASVs{$asv}{$sample} = $split_line[$i];
            }
    }

open(IN2, "<$options{m}") or die "\n\nNADA $options{m} file!!\n\n";
my @sampdata = <IN2>; close(IN2);
my $metaheada = shift(@sampdata); chomp($metaheada);
my @metaheaders = split('\t', $metaheada);

my $sample_column;
my $replicate_column;

foreach my $i (0..$#metaheaders)
    {   my $header = $metaheaders[$i];
        if ($header eq 'Sample')
            {   $sample_column = $i; 
            }
        if ($header eq 'replicates')
            {   $replicate_column = $i;   
            }
    }

foreach my $line (@sampdata)
    {   chomp($line);
        my @split_line = split('\t', $line);
        my $sample = $split_line[$sample_column];
        my $replicatename = $split_line[$replicate_column];
        $META{$replicatename}{'samples'} .= $sample.";";
    }
    
my @replicates;
my $max_rep_number = 0;
foreach my $i (sort keys %META)
    {   my $samples = $META{$i}{'samples'};
        my @split_samples = split(';', $samples);
        my @new_samples;
        my $nohit = 0;
        foreach my $j (@split_samples)
            {   foreach my $g (@sampleheads)
                    {   unless($g eq "x")
                        {   if ($j eq $g)
                               {    push(@new_samples, $j);
                                    $nohit = 1;
                               }
                        }
                    }
            }
        if ($nohit == 0)
            {   $META{$i}{'count'} = 0;
            }
        if ($nohit == 1)
            {   my $new_samp = join(';', @new_samples);
                $META{$i}{'samples'} = $new_samp;
                $META{$i}{'count'} = $#new_samples + 1;
                if ($#new_samples + 1 > $max_rep_number)
                    {   $max_rep_number = $#new_samples + 1;
                    }
                push(@replicates, $i);
            }
    }

print "ASV\tsite\tnum_replicates\tnum_present\tsum_counts\tavg_counts_wherefound\tavg_counts_allrep\n";

foreach my $i (sort keys %ASVs)
    {   foreach my $j (@replicates)
            {   print "$i\t$j\t";
                my $samples = $META{$j}{'samples'};
                my @rep_samples = split(';', $samples);
                my $c = $#rep_samples + 1;
                print "$c\t";
                my $num_present = 0;
                foreach my $k (0..$#rep_samples)
                    {   my $final_samp = $rep_samples[$k];
                        my $value = $ASVs{$i}{$final_samp};
                        if ($value > 0)
                            {   $num_present += 1;
                            }
                    }
                print "$num_present\t";
                my $sum_count = 0;
                foreach my $k (0..$#rep_samples)
                {   my $final_samp = $rep_samples[$k];
                    my $value = $ASVs{$i}{$final_samp};
                    $sum_count += $value;
                }
                print "$sum_count\t";
                my $avg_wherefound;
                my $avg_allrep;
                if ($num_present == 0)
                    {   $avg_wherefound = 0;
                        $avg_allrep = 0;
                    }
                else
                    {   $avg_wherefound = $sum_count / $num_present;
                        $avg_allrep = $sum_count / $c;
                    }
                print "$avg_wherefound\t$avg_allrep\n";
            }
    }

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - S U B R O U T I N E S - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -
