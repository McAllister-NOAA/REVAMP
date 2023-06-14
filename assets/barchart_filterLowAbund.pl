#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use List::Util qw( min max );
use List::MoreUtils qw(uniq);

# - - - - - H E A D E R - - - - - - - - - - - - - - - - -
#Goals of script:
#Collapse hits for low abundance taxa to zzOther for barchart R input

# - - - - - C O M M A N D    L I N E    O P T I O N S - - - - - - - -
my %options=();
getopts("i:f:h", \%options);

if ($options{h})
    {   print "\n\nHelp called:\nOptions:\n";
        print "-i = Input barchart for R file\n";
        print "-f = Filter percent (low abund cutoff)\n";
        print "-h = This help message\n\n";
        die;
    }

# - - - - - G L O B A L  V A R I A B L E S  - - - - - -
my %BARCHART;
my @terminaltax;
my @samples;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - M A I N - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
open(IN,     "<$options{i}") or die "\n\nNADA $options{i} file!!\n\n";
my @data = <IN>; close(IN); shift(@data);
foreach my $line (@data)
    {   chomp($line);
        my @split_line = split('\t', $line);
        my $value = $split_line[0];
        my $sample = $split_line[1];
        push(@samples, $sample);
        my $termtax = $split_line[2];
        unless ($termtax eq "Unknown" || $termtax eq "Environmental Unknown")
        {   push(@terminaltax, $termtax);
            $BARCHART{$termtax}{$sample} = $value;
            $BARCHART{$termtax}{'keep'} = "FALSE";
        }
    }

my @unique_samples = uniq @samples;

foreach my $j (@unique_samples)
    {   my $sampletotal = 0;
        foreach my $i (sort keys %BARCHART)
            {   if (exists $BARCHART{$i}{$j}) {
                $sampletotal += $BARCHART{$i}{$j};}
            }
        foreach my $i (sort keys %BARCHART)
            {   if (exists $BARCHART{$i}{$j}) {
                my $value = $BARCHART{$i}{$j};
                my $newperc = (100 * $value)/$sampletotal;
                $BARCHART{$i}{$j} = $newperc;
                }
            }
    }

foreach my $i (sort keys %BARCHART)
    {   my $keep = "FALSE";
        foreach my $j (@unique_samples)
            {   if (exists $BARCHART{$i}{$j}) {
                if ($BARCHART{$i}{$j} >= $options{f})
                    {   $keep = "TRUE";
                    }
            }
            }
        if ($keep eq "TRUE")
            {   $BARCHART{$i}{'keep'} = "TRUE"
            }
    }

print "Percent\tSample\tTerminalTaxa\n";
foreach my $i (sort keys %BARCHART)
    {   foreach my $j (@unique_samples)
            {   if (exists $BARCHART{$i}{$j}) {
                if ($BARCHART{$i}{'keep'} eq "TRUE")
                    {   print "$BARCHART{$i}{$j}\t$j\t$i\n";
                    }
                else
                    {   print "$BARCHART{$i}{$j}\t$j\tzzOther\n";
                    }
            }
            }
    }




# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - S U B R O U T I N E S - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -
