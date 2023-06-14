#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use List::MoreUtils qw(uniq);
use List::Util qw( min max );

# - - - - - H E A D E R - - - - - - - - - - - - - - - - -
#Goals of script:
#Count the number of best hits and print out per query.

# - - - - - C O M M A N D    L I N E    O P T I O N S - - - - - - - -
my %options=();
getopts("i:c:h", \%options);

if ($options{h})
    {   print "\n\nHelp called:\nOptions:\n";
        print "-i = Input btab BLASTn file\n";
        print "-c = Count of ASVs\n";
        print "-h = This help message\n\n";
        die;
    }

# - - - - - G L O B A L  V A R I A B L E S  - - - - - -
my %ASV;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - M A I N - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
open(IN, "<$options{i}") or die "\n\nNADA $options{i} file!!\n\n";
my @data = <IN>; close(IN);
foreach my $line (@data)
    {   chomp($line);
        my @split_line = split('\t', $line);
        my $asv = $split_line[0];
        my $pid = $split_line[1];
        $ASV{$asv}{$pid} += 1;
        $ASV{$asv}{'pidoptions'} .= $pid.";";
    }

my $asv_count = $options{c};
chomp($asv_count);

foreach my $i (1..$asv_count)
    {   my $calledasv = "ASV_".$i;
        if (exists $ASV{$calledasv})
            {   my $pidoptions = $ASV{$calledasv}{'pidoptions'};
                my @split_pidoptions = split(';', $pidoptions);
                my @uniq_pidoptions = uniq @split_pidoptions;
                my $bestpid = max @uniq_pidoptions;
                print "$calledasv\t$ASV{$calledasv}{$bestpid}\t$bestpid\n";
            }
        else
            {   print "$calledasv\t0\tNA\n"
            }
    }
    
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - S U B R O U T I N E S - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -
