#!/usr/bin/perl -w
use strict;
use Getopt::Std;

# - - - - - H E A D E R - - - - - - - - - - - - - - - - -
#Goals of script:
#Create reproducible, trustworthy, metadata file for R from user input.

# - - - - - C O M M A N D    L I N E    O P T I O N S - - - - - - - -
my %options=();
getopts("i:h", \%options);

if ($options{h})
    {   print "\n\nHelp called:\nOptions:\n";
        print "-i = sample metadata file input\n";
        print "-h = This help message\n\n";
        die;
    }

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - M A I N - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
open(IN, "<$options{i}") or die "\n\nThere is no $options{i} file!!\n\n";
my @data = <IN>; close(IN);
my $heads = shift(@data);
$heads =~ s/\n//g;
$heads =~ s/\r//g;

my $num_tabs = () = $heads =~ /\t/g;

my $terminalWhitespace = "FALSE";
if ($heads =~ m/\t$/)
    {   $terminalWhitespace = "TRUE"; #GOING TO ASSUME IT IS NOT MEANT TO BE THERE and THAT THERE IS ONLY ONE
        $heads =~ s/\t$//;
    }

my $scalarHeaders;
if ($terminalWhitespace eq "FALSE")
    {   $scalarHeaders = $num_tabs + 1;
    }
else
    {   $scalarHeaders = $num_tabs;
    }

my @splitheads = split('\t', $heads);
if ($scalarHeaders != scalar @splitheads)
    {   die "\n\nThere is more than one trailing tab in the headers. Make sure you have a header for all columns.\n\n";
    }
my $numHeaders = scalar @splitheads;
my $lastEntry = $numHeaders - 1;

#Attempt to match controlled header vocabulary/case:
my $groupstartNUM = 1;
my $numBlankHeaders = 1;
foreach my $i (0..$#splitheads)
    {   my $thisHeader = $splitheads[$i];
        my $match = "FALSE";
        if ($thisHeader =~ m/^[sS][aA][mM][pP][lL][eE]$/ || $thisHeader =~ m/^[sS][aA][mM][pP][lL][eE][sS]$/ || $thisHeader =~ m/^[sS][aA][mM][pP]$/)
            {$splitheads[$i] = "Sample"; $match = "TRUE";}
        if ($thisHeader =~ m/^[lL][aA][tT]$/ || $thisHeader =~ m/^[lL][aA][tT][iI][tT][uU][dD][eE]$/)
            {$splitheads[$i] = "lat"; $match = "TRUE";}
        if ($thisHeader =~ m/^[lL][oO][nN][gG]$/ || $thisHeader =~ m/^[lL][oO][nN][gG][iI][tT][uU][dD][eE]$/)
            {$splitheads[$i] = "long"; $match = "TRUE";}
        if ($thisHeader =~ m/^[rR][eE][pP][lL][iI][cC][aA][tT][eE][sS]$/ || $thisHeader =~ m/^[rR][eE][pP][lL][iI][cC][aA][tT][eE]$/ || $thisHeader =~ m/^[rR][eE][pP]$/ || $thisHeader =~ m/^[rR][eE][pP][sS]$/)
            {$splitheads[$i] = "replicates"; $match = "TRUE";}
        if ($thisHeader =~ m/^[sS][iI][tT][eE][sS]$/ || $thisHeader =~ m/^[sS][iI][tT][eE]$/)
            {$splitheads[$i] = "sites"; $match = "TRUE";}
        if ($thisHeader =~ m/^[cC][oO][nN][tT][rR][oO][lL][sS]$/ || $thisHeader =~ m/^[cC][oO][nN][tT][rR][oO][lL]$/)
            {$splitheads[$i] = "controls"; $match = "TRUE";}
        if ($thisHeader =~ m/^[gG][rR][oO][uU][pP][a-zA-Z]+$/)
            {   $splitheads[$i] = "group".$groupstartNUM;
                $groupstartNUM += 1;
                $match = "TRUE";
            }
        if ($thisHeader =~ m/^[gG][rR][oO][uU][pP][0-9]+$/)
            {   $splitheads[$i] = "group".$groupstartNUM;
                $groupstartNUM += 1;
                $match = "TRUE";
            }
        if ($match eq "FALSE")
            {   $thisHeader =~ s/[^a-zA-Z0-9]/_/g;
                $splitheads[$i] = $thisHeader;
            }
        if ($thisHeader eq "")
            {   $splitheads[$i] = "BLANK_".$numBlankHeaders;
                $numBlankHeaders += 1;
            }
    }
    
my $samp_col;
foreach my $i (0..$#splitheads)
    {   if ($splitheads[$i] eq "Sample")
            {   $samp_col = $i;
            }
    }
    
if ($samp_col != 0)
    {   die "\n\nPlease make your Sample name column the first column in the metadata file!\n\n";
    }

foreach my $i (0..$#splitheads)
    {   print "$splitheads[$i]";
        if ($i == $#splitheads)
            {   print "\n";}
        else {  print "\t";}
    }

foreach my $line (@data)
    {   $line =~ s/\n//g;
        $line =~ s/\r//g;
        my @vals = split('\t', $line);
        foreach my $i (0..$lastEntry)
             {  if (exists $vals[$i])
                    {   if ($i == $samp_col)
                             {   my $name = $vals[$i];
                                 if ($name eq "")
                                     {   die "\n\nYour samples need to have names!\n\n";
                                     }
                                 $name =~ s/[^a-zA-Z0-9]/_/g;
                                 unless ($name =~ m/^MP_/)
                                     {   $name = "MP_".$name;
                                     }
                                 if ($i == $lastEntry)
                                     {   print "$name";
                                     }
                                 else
                                     {   print "$name\t";
                                     }
                             }
                         else
                             {   my $entry = $vals[$i];
                                 if ($entry eq "")
                                     {   $entry = "NA";
                                     }
                                 if ($i == $lastEntry)
                                     {   print "$entry";
                                     }
                                 else
                                     {   print "$entry\t";
                                     }
                             }
                    }
                else
                    {   if ($i == $lastEntry)
                            {   print "NA";
                            }
                        else
                            {   print "NA\t";
                            }
                    }
             }
        print "\n";
    }

# - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -