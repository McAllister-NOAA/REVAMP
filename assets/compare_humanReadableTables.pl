#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use List::Util qw( min max );
use List::MoreUtils qw(uniq);

# - - - - - H E A D E R - - - - - - - - - - - - - - - - -
#Goals of script:
#Merge taxonomy2PercentAbundance_humanReadable.txt from each marker into one file and set reference table (i.e. Morphology Carbon or the like).

# - - - - - C O M M A N D    L I N E    O P T I O N S - - - - - - - -
my %options=();
getopts("i:r:s:h", \%options);

if ($options{h})
    {   print "\n\nHelp called:\nOptions:\n";
        print "-i = Input directory with taxonomy2PercentAbundance_humanReadable.txt files renamed as Marker.txt\n";
        print "-r = Reference humanReadable file name (Marker.txt) within input directory for comparison to all other markers\n";
        print '-s = Tab delimited file: "Marker(no .txt)\tMarkerSample\tReferenceSample" one per line to establish equivalences. No header.';
        print "\n";
        print "-h = This help message\n\n";
        die;
    }

# - - - - - G L O B A L  V A R I A B L E S  - - - - - -
my %REF;
my %SAMPEQUIV;
my $maxsampnum;
my @taxaHierarchy = qw | Kingdom Phylum Class Order Family Genus Species |;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - M A I N - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
open(IN, "<$options{i}/$options{r}") or die "\n\nThere is no $options{i}/$options{r}!!\n\n";
my $ref_marker = $options{r};
$ref_marker =~ s/\.txt$//;
my @data = <IN>; close(IN);
shift(@data);
my $header = shift(@data); chomp($header);
my @head_array = split('\t', $header);
shift(@head_array); shift(@head_array);
my @samples;
foreach my $i (0..$#head_array)
    {   my $sample = $head_array[$i];
        chomp($sample);
        push(@samples, $sample);
    }
shift(@data); shift(@data); shift(@data);
my $ref_linenumb = 0;
$maxsampnum = scalar(@samples);
foreach my $line (@data)
    {   $ref_linenumb += 1;
        chomp($line);
        unless ($line =~ m/<<<.+>>>/ || $line !~ m/;/)
            {   my @line_array = split('\t', $line);
                my $taxaString = shift(@line_array); chomp($taxaString);
                my $terminalTaxa = shift(@line_array); chomp($terminalTaxa);
                $REF{$taxaString}{'terminal_taxa'} = $terminalTaxa;
                foreach my $i (0..$#samples)
                    {   my $sample = $samples[$i]; chomp($sample);
                        $REF{$taxaString}{$sample} = $line_array[$i];
                    }
            }
    }

open(IN, "<$options{s}") or die "\n\nThere is no $options{s}!!\n\n";
my @data2 = <IN>; close(IN);
foreach my $line (@data2)
    {   chomp($line);
        my @lineArray = split('\t', $line);
        my $marker = $lineArray[0]; chomp($marker);
        my $markersample = $lineArray[1]; chomp($markersample);
        my $referencesample = "NA";
        if (exists $lineArray[2])
            {   $referencesample = $lineArray[2]; chomp($referencesample);}
        my $key_text = $marker."---".$markersample;
        if ($referencesample ne "NA")
            {   $SAMPEQUIV{$key_text}{'refsample'} = $referencesample;}
    }


my $total_marker_linenumb = 0;
my $total_marker_num = 0;
opendir(DIR2, $options{i}) or die "\n\nInput directory $options{i} not found\n\n";
while (my $file = readdir(DIR2))
    {   if ($file =~ m/\.txt$/ && $file ne $options{r})
            {   open(IN, "<$options{i}/$file") or die "\n\nThere is no $file!!\n\n";
                my @data = <IN>; close(IN);
                $total_marker_num += 1;
                shift(@data);
                my $header = shift(@data); chomp($header);
                my @head_array = split('\t', $header);
                shift(@head_array); shift(@head_array);
                my @marker_samples;
                foreach my $i (0..$#head_array)
                    {   my $sample = $head_array[$i];
                        chomp($sample);
                        push(@marker_samples, $sample);
                    }
                my $marker_samp_scalar = scalar(@marker_samples);
                if ($marker_samp_scalar > $maxsampnum)
                    { $maxsampnum = $marker_samp_scalar;
                    }
                foreach my $i (@data)
                    {$total_marker_linenumb += 1;}
            }
    }

my $tot_lines_output = $total_marker_linenumb + ($total_marker_num * $ref_linenumb);
my $line_log;
if ($tot_lines_output < 1000)
    {$line_log = "1001";}
if ($tot_lines_output >= 1000 && $tot_lines_output < 10000)
    {$line_log = "10001";}
if ($tot_lines_output >= 10000 && $tot_lines_output < 100000)
    {$line_log = "100001";}
if ($tot_lines_output >= 100000 && $tot_lines_output < 1000000)
    {$line_log = "1000001";}
if ($tot_lines_output >= 1000000)
    {$line_log = "1000000001";}

print "LineNumb\t(M)arker/(R)eference\tTaxaString\tTerminalTaxa";
foreach my $i (1..$maxsampnum)
    {   print "\tSample$i";
    }
print "\n";
print "$line_log\t\tAll percentages range 0-100% of the total count of each marker for either Eukaryotes or Prokaryotes, with exception of the first three data rows (% Known, % Eukaryote, % Prokaryote), which represent %s of total reads per marker.\n";
$line_log += 1;
print "$line_log\n";
$line_log += 1;
opendir(DIR, $options{i}) or die "\n\nInput directory $options{i} not found\n\n";
while (my $file = readdir(DIR))
    {   if ($file =~ m/\.txt$/ && $file ne $options{r})
            {   my $markername = $file;
                $markername =~ s/\.txt$//;
                open(IN, "<$options{i}/$file") or die "\n\nThere is no $file!!\n\n";
                my @data = <IN>; close(IN);
                shift(@data);
                print "$line_log\t\tMarker: $markername against reference $ref_marker.\n";
                $line_log += 1;
                my $header = shift(@data); chomp($header);
                print "$line_log\t\t$header\n";
                $line_log += 1;
                my @head_array = split('\t', $header);
                shift(@head_array); shift(@head_array);
                my @marker_samples;
                foreach my $i (0..$#head_array)
                    {   my $sample = $head_array[$i];
                        chomp($sample);
                        push(@marker_samples, $sample);
                    }
                my %REF_currentSamples; #Reference database with only those taxa that are in samples found in the compared marker file.
                foreach my $refkey (sort keys %REF)
                    {   my $totalperc = 0;
                        foreach my $i (@marker_samples)
                            {   my $sampequivkey = $markername."---".$i;
                                my $refsamp = "NULL";
                                if (exists $SAMPEQUIV{$sampequivkey})
                                    {   $refsamp = $SAMPEQUIV{$sampequivkey}{'refsample'};
                                    }
                                my $value = 0;
                                if ($refsamp ne "NULL")
                                    {   $value = $REF{$refkey}{$refsamp};
                                    }
                                if ($value eq "D")
                                    {$totalperc += 1;}
                                else
                                    {$totalperc += $value;}
                            }
                        if ($totalperc > 0)
                            {   $REF_currentSamples{$refkey}{'printed'} = "FALSE";
                                foreach my $i (@marker_samples)
                                    {   my $sampequivkey = $markername."---".$i;
                                        if (exists $SAMPEQUIV{$sampequivkey})
                                            {   my $refsamp = $SAMPEQUIV{$sampequivkey}{'refsample'};
                                                my $value = $REF{$refkey}{$refsamp};
                                                my $terminaltaxa = $REF{$refkey}{'terminal_taxa'};
                                                $REF_currentSamples{$refkey}{$refsamp} = $value;
                                                $REF_currentSamples{$refkey}{'terminal_taxa'} = $terminaltaxa;
                                            }
                                    }
                            }
                    }
                my $line1 = shift(@data); chomp($line1);
                my $line2 = shift(@data); chomp($line2);
                my $line3 = shift(@data); chomp($line3);
                print "$line_log\tM\t$line1\n";
                $line_log += 1;
                print "$line_log\tM\t$line2\n";
                $line_log += 1;
                print "$line_log\tM\t$line3\n";
                $line_log += 1;
                my $currentPhylum = "";
                my $currentClass = "";
                foreach my $line (@data)
                    {   chomp($line);
                        if ($line =~ m/<<<.+>>>/)
                            {   print "$line_log\t\t$line\n";
                                $line_log += 1;
                            }
                        else
                            {   print "$line_log\tM\t$line\n";
                                $line_log += 1;
                                my @line_array = split('\t', $line);
                                my $taxaString = shift(@line_array); chomp($taxaString);
                                my $terminalTaxa = shift(@line_array); chomp($terminalTaxa);
                                my @taxa_array = split(';', $taxaString);
                                $currentPhylum = $taxa_array[1];
                                $currentClass = $taxa_array[2];
                                if (exists $REF_currentSamples{$taxaString})
                                    {   print "$line_log\tR\t\t";
                                        foreach my $i (0..$#marker_samples)
                                            {   my $sampequivkey = $markername."---".$marker_samples[$i];
                                                my $refsamp = "NULL";
                                                if (exists $SAMPEQUIV{$sampequivkey}{'refsample'})
                                                    {   $refsamp = $SAMPEQUIV{$sampequivkey}{'refsample'};
                                                    }
                                                if ($refsamp eq "NULL")
                                                    {   print "\tNA";
                                                    }
                                                else
                                                    {   if (exists $REF_currentSamples{$taxaString}{$refsamp})
                                                            {   print "\t$REF_currentSamples{$taxaString}{$refsamp}";
                                                            }
                                                        else
                                                            {   print "\tNA";}
                                                    }
                                            }
                                        print "\n";
                                        $line_log += 1;
                                        $REF_currentSamples{$taxaString}{'printed'} = "TRUE";
                                    }
                            }
                    }
                my $printnonmatchswitch = "FALSE";
                foreach my $i (sort keys %REF_currentSamples)
                    {   if ($REF_currentSamples{$i}{'printed'} eq "FALSE")
                            {   $printnonmatchswitch = "TRUE";
                            }
                    }
                if ($printnonmatchswitch eq "TRUE")
                    {   print "$line_log\t\tThe following are taxa strings found in reference that were not identified in this marker:\n";
                        $line_log += 1;
                        foreach my $i (sort keys %REF_currentSamples)
                            {   if ($REF_currentSamples{$i}{'printed'} eq "FALSE")
                                    {   print "$line_log\tR\t$i\t$REF_currentSamples{$i}{'terminal_taxa'}";
                                        foreach my $j (0..$#marker_samples)
                                            {   my $sampequivkey = $markername."---".$marker_samples[$j];
                                                my $refsamp = "NULL";
                                                if (exists $SAMPEQUIV{$sampequivkey}{'refsample'})
                                                    {   $refsamp = $SAMPEQUIV{$sampequivkey}{'refsample'};
                                                    }
                                                if ($refsamp eq "NULL")
                                                    {   print "\tNA";
                                                    }
                                                else
                                                    {   if (exists $REF_currentSamples{$i}{$refsamp})
                                                            {   print "\t$REF_currentSamples{$i}{$refsamp}";
                                                            }
                                                        else
                                                            {   print "\tNA";}
                                                    }
                                            }
                                        print "\n";
                                        $line_log += 1;
                                    }
                            }
                    }
                print "$line_log\n";
                $line_log += 1;
                print "$line_log\n";
                $line_log += 1;
                print "$line_log\n";
                $line_log += 1;
                print "$line_log\n";
                $line_log += 1;
            }
    }
    
    
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - S U B R O U T I N E S - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -
