#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use List::MoreUtils qw(uniq);


# - - - - - H E A D E R - - - - - - - - - - - - - - - - -
#Goals of script:
#Import ASVs_counts.tsv and x_asvTaxonomyTable.txt
#Output ASV table merged on identical taxonomy

# - - - - - C O M M A N D    L I N E    O P T I O N S - - - - - - - -
my %options=();
getopts("a:t:o:h", \%options);

if ($options{h})
    {   print "\n\nHelp called:\nOptions:\n";
        print "-a = DADA2 ASVs_counts.tsv file\n";
        print "-t = ASV taxonomy table from processing script (x_asvTaxonomyTable.txt)\n";
        print "-o = Output directory\n";
        print "-h = This help message\n\n";
        die;
    }

# - - - - - G L O B A L  V A R I A B L E S  - - - - - -
my %ASV;
my %TAX;
my %NEW_ASV;
my @taxHier = qw | Kingdom Phylum Class Order Family Genus Species |;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - M A I N - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

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
    
#IMPORT ASV Taxonomy table
open(IN2, "<$options{t}") or die "\n\nThere is no $options{t} file!!\n\n";
my @tax_dat = <IN2>; close(IN2); shift(@tax_dat);
foreach my $line (@tax_dat)
	{	chomp($line);
		my @data = split('\t', $line);
		my $asv = $data[0]; chomp($asv);
        my $taxstring = $data[1].";".$data[2].";".$data[3].";".$data[4].";".$data[5].";".$data[6].";".$data[7];
		$TAX{$taxstring}{'asvs'} .= $asv.";";
	}
    
#Collapse ASV counts on identical taxonomy
foreach my $i (sort keys %TAX)
    {   my $common_asvs = $TAX{$i}{'asvs'};
        my @commons = split(';', $common_asvs);
        my $asv_to_use = $commons[0];
        $NEW_ASV{$asv_to_use}{'taxonomy'} = $i;
        foreach my $j (1..$#file_sample_headers)
            {   my $sample = $file_sample_headers[$j]; chomp($sample);
                $NEW_ASV{$asv_to_use}{$sample} = $ASV{$asv_to_use}{$sample};
                foreach my $k (1..$#commons)
                    {   $NEW_ASV{$asv_to_use}{$sample} += $ASV{$commons[$k]}{$sample};
                    }
            }
    }
    
#Print out new ASV table
print "x";
foreach my $i (1..$#file_sample_headers)
    {   my $sample = $file_sample_headers[$i];
        chomp($sample);
        print "\t$sample";
    }
print "\n";
foreach my $i (sort keys %NEW_ASV)
    {   print "$i";
        foreach my $k (1..$#file_sample_headers)
            {   my $sample = $file_sample_headers[$k];
                chomp($sample);
                print "\t$NEW_ASV{$i}{$sample}";
            }
        print "\n";
    }
    
#Print out ASV table combining asv2Taxonomy with perc abundance for human viewing
open(OUT, ">".$options{o}."/taxonomy2PercentAbundance_humanReadable.txt");
open(OUT_NOROUND, ">".$options{o}."/taxonomy2PercentAbundance_humanReadable_NoRounding.txt");

my @sample_totals;
my @sample_totals_sansUnknowns;
my @sample_totals_Eukaryota;
my @sample_totals_Prokaryotes;
push(@sample_totals, "NULL");
push(@sample_totals_sansUnknowns, "NULL");
push(@sample_totals_Eukaryota, "NULL");
push(@sample_totals_Prokaryotes, "NULL");
foreach my $i (1..$#file_sample_headers)
    {   my $sample = $file_sample_headers[$i];
        chomp($sample);
        my $thisSampleTotal = 0;
        my $thisSampleNoUNKNOWN = 0;
        my $thisSampleEuk = 0;
        my $thisSampleProk = 0;
        foreach my $j (sort keys %NEW_ASV)
            {   $thisSampleTotal += $NEW_ASV{$j}{$sample};
                my $taxa = $NEW_ASV{$j}{'taxonomy'};
                chomp($taxa);
                unless ($taxa =~ m/^Unknown/)
                    {   $thisSampleNoUNKNOWN += $NEW_ASV{$j}{$sample};}
                if ($taxa =~ m/^Eukaryota/)
                    {   $thisSampleEuk += $NEW_ASV{$j}{$sample};}
                if ($taxa =~ m/^Bacteria/ || $taxa =~ m/^Archaea/)
                    {   $thisSampleProk += $NEW_ASV{$j}{$sample};}
            }
        push(@sample_totals, $thisSampleTotal);
        push(@sample_totals_sansUnknowns, $thisSampleNoUNKNOWN);
        push(@sample_totals_Eukaryota, $thisSampleEuk);
        push(@sample_totals_Prokaryotes, $thisSampleProk);
    }

my %orderedDB;
foreach my $i (sort keys %NEW_ASV)
    {   my $taxa = $NEW_ASV{$i}{'taxonomy'};
        $orderedDB{$taxa}{'asv_key'} = $i;
        foreach my $j (1..$#file_sample_headers)
            {   my $sample = $file_sample_headers[$j];
                chomp($sample);
                my $value = $NEW_ASV{$i}{$sample};
                $orderedDB{$taxa}{$sample} = $value;
            }
    }

print OUT "All percentages range 0-100% of the total count for either Eukaryotes or Prokaryotes, with exception of the first three data rows (% Known, % Eukaryote, % Prokaryote), which represent %s of total reads.\n";
print OUT "Taxonomy String (K->S)\t";
print OUT "Terminal Taxa";
print OUT_NOROUND "All percentages range 0-100% of the total count for either Eukaryotes or Prokaryotes, with exception of the first three data rows (% Known, % Eukaryote, % Prokaryote), which represent %s of total reads.\n";
print OUT_NOROUND "Taxonomy String (K->S)\t";
print OUT_NOROUND "Terminal Taxa";
foreach my $i (1..$#file_sample_headers)
    {   my $sample = $file_sample_headers[$i];
        chomp($sample);
        print OUT "\t$sample";
        print OUT_NOROUND "\t$sample";
    }
print OUT "\n";
print OUT_NOROUND "\n";

print OUT "\t% Known";
print OUT_NOROUND "\t% Known";
foreach my $i (1..$#file_sample_headers)
    {   my $sample = $file_sample_headers[$i];
        chomp($sample);
        my $number = 100 * ($sample_totals_sansUnknowns[$i] / $sample_totals[$i]);
        my $roundnum = &ROUND($number,2);
        if ($roundnum == 0 && $number > 0)
            {   print OUT "\tD";
            }
        else
            {   print OUT "\t$roundnum";   
            }
        print OUT_NOROUND "\t$number";
    }
print OUT "\n";
print OUT_NOROUND "\n";

print OUT "\t% Eukaryote";
print OUT_NOROUND "\t% Eukaryote";
foreach my $i (1..$#file_sample_headers)
    {   my $sample = $file_sample_headers[$i];
        chomp($sample);
        my $number = 100 * ($sample_totals_Eukaryota[$i] / $sample_totals[$i]);
        my $roundnum = &ROUND($number,2);
        if ($roundnum == 0 && $number > 0)
            {   print OUT "\tD";
            }
        else
            {   print OUT "\t$roundnum";   
            }
        print OUT_NOROUND "\t$number";
    }
print OUT "\n";
print OUT_NOROUND "\n";

print OUT "\t% Prokaryote";
print OUT_NOROUND "\t% Prokaryote";
foreach my $i (1..$#file_sample_headers)
    {   my $sample = $file_sample_headers[$i];
        chomp($sample);
        my $number = 100 * ($sample_totals_Prokaryotes[$i] / $sample_totals[$i]);
        my $roundnum = &ROUND($number,2);
        if ($roundnum == 0 && $number > 0)
            {   print OUT "\tD";
            }
        else
            {   print OUT "\t$roundnum";   
            }
        print OUT_NOROUND "\t$number";
    }
print OUT "\n";
print OUT_NOROUND "\n";

my $previousPhyla = "";
my $previousClass = "";

foreach my $i (sort keys %orderedDB)
    {   my @taxa_array = split(';', $i);
        my $taxastring = $i;
        chomp($taxastring);
        my $phylum = $taxa_array[1];
        my $class = $taxa_array[2];
        my $order = $taxa_array[3];
        my $family = $taxa_array[4];
        my $lastReal;
        foreach my $j (reverse 0..$#taxa_array)
            {   if ($taxa_array[$j] ne "NA")
                    {   $lastReal = $j;
                        last;
                    }
            }
        my $termtax = $taxa_array[$lastReal];
        if ($taxa_array[0] eq "Eukaryota")
            {   if ($phylum eq $previousPhyla)
                    {   if ($class eq $previousClass)
                            {   print OUT "$taxastring\t";
                                print OUT "$termtax";
                                print OUT_NOROUND "$taxastring\t";
                                print OUT_NOROUND "$termtax";
                                foreach my $j (1..$#file_sample_headers)
                                    {   my $sample = $file_sample_headers[$j];
                                        chomp($sample);
                                        my $number;
                                        if ($sample_totals_Eukaryota[$j] == 0)
                                            {   $number = 0;}
                                        else
                                            {   $number = 100 * ($orderedDB{$i}{$sample} / $sample_totals_Eukaryota[$j]);}
                                        my $roundnum = &ROUND($number,2);
                                        if ($roundnum == 0 && $number > 0)
                                            {   print OUT "\tD";
                                            }
                                        else
                                            {   print OUT "\t$roundnum";   
                                            }
                                        print OUT_NOROUND "\t$number";
                                    }
                                print OUT "\n";
                                print OUT_NOROUND "\n";
                            }
                        else
                            {   $previousClass = $class;
                                if ($phylum eq "Arthropoda" || $phylum eq "Chordata")
                                    {   print OUT "\t<<<$phylum, $class>>>\n";
                                        print OUT_NOROUND "\t<<<$phylum, $class>>>\n";
                                    }
                                print OUT "$taxastring\t";
                                print OUT "$termtax";
                                print OUT_NOROUND "$taxastring\t";
                                print OUT_NOROUND "$termtax";
                                foreach my $j (1..$#file_sample_headers)
                                    {   my $sample = $file_sample_headers[$j];
                                        chomp($sample);
                                        my $number;
                                        if ($sample_totals_Eukaryota[$j] == 0)
                                            {   $number = 0;}
                                        else
                                            {   $number = 100 * ($orderedDB{$i}{$sample} / $sample_totals_Eukaryota[$j]);
                                            }
                                        my $roundnum = &ROUND($number,2);
                                        if ($roundnum == 0 && $number > 0)
                                            {   print OUT "\tD";
                                            }
                                        else
                                            {   print OUT "\t$roundnum";   
                                            }
                                        print OUT_NOROUND "\t$number";
                                    }
                                print OUT "\n";
                                print OUT_NOROUND "\n";
                            }
                    }
                else
                    {   $previousPhyla = $phylum;
                        $previousClass = $class;
                        if ($phylum eq "Arthropoda" || $phylum eq "Chordata")
                            {   print OUT "\t<<<$phylum, $class>>>\n";
                                print OUT_NOROUND "\t<<<$phylum, $class>>>\n";
                            }
                        elsif ($phylum eq "NA")
                            {   print OUT "\t<<<Unknown Eukaryote>>>\n";
                                print OUT_NOROUND "\t<<<Unknown Eukaryote>>>\n";
                            }
                        else
                            {   print OUT "\t<<<$phylum>>>\n";
                                print OUT_NOROUND "\t<<<$phylum>>>\n";
                            }
                        print OUT "$taxastring\t";
                        print OUT "$termtax";
                        print OUT_NOROUND "$taxastring\t";
                        print OUT_NOROUND "$termtax";
                        foreach my $j (1..$#file_sample_headers)
                            {   my $sample = $file_sample_headers[$j];
                                chomp($sample);
                                my $number;
                                if ($sample_totals_Eukaryota[$j] == 0)
                                    {   $number = 0;}
                                else
                                    {   $number = 100 * ($orderedDB{$i}{$sample} / $sample_totals_Eukaryota[$j]);
                                    }
                                my $roundnum = &ROUND($number,2);
                                if ($roundnum == 0 && $number > 0)
                                    {   print OUT "\tD";
                                    }
                                else
                                    {   print OUT "\t$roundnum";   
                                    }
                                print OUT_NOROUND "\t$number";
                            }
                        print OUT "\n";
                        print OUT_NOROUND "\n";
                    }
            }
    }
    
$previousPhyla = "";
foreach my $i (sort keys %orderedDB)
    {   my @taxa_array = split(';', $i);
        my $taxastring = $i;
        chomp($taxastring);
        my $phylum = $taxa_array[1]; chomp($phylum);
        my $class = $taxa_array[2]; chomp($class);
        my $order = $taxa_array[3]; chomp($order);
        my $family = $taxa_array[4]; chomp($family);
        my $lastReal;
        foreach my $j (reverse 0..$#taxa_array)
            {   if ($taxa_array[$j] ne "NA")
                    {   $lastReal = $j;
                        last;
                    }
            }
        my $termtax = $taxa_array[$lastReal]; chomp($termtax);
        if ($taxa_array[0] eq "Bacteria" || $taxa_array[0] eq "Archaea")
            {   if ($phylum eq $previousPhyla)
                    {   print OUT "$taxastring\t";
                        print OUT "$termtax";
                        print OUT_NOROUND "$taxastring\t";
                        print OUT_NOROUND "$termtax";
                        foreach my $j (1..$#file_sample_headers)
                            {   my $sample = $file_sample_headers[$j];
                                chomp($sample);
                                my $number;
                                if ($sample_totals_Prokaryotes[$j] == 0)
                                    {   $number = 0;}
                                else
                                    {   $number = 100 * ($orderedDB{$i}{$sample} / $sample_totals_Prokaryotes[$j]);
                                    }
                                my $roundnum = &ROUND($number,2);
                                if ($roundnum == 0 && $number > 0)
                                    {   print OUT "\tD";
                                    }
                                else
                                    {   print OUT "\t$roundnum";   
                                    }
                                print OUT_NOROUND "\t$number"; 
                            }
                        print OUT "\n";
                        print OUT_NOROUND "\n";
                    }
                else
                    {   $previousPhyla = $phylum;
                        if ($phylum eq "NA")
                            {   print OUT "\t<<<Unknown $taxa_array[0]>>>\n";
                                print OUT_NOROUND "\t<<<Unknown $taxa_array[0]>>>\n";
                            }
                        else
                            {   print OUT "\t<<<$phylum>>>\n";
                                print OUT_NOROUND "\t<<<$phylum>>>\n";
                            }
                        print OUT "$taxastring\t";
                        print OUT "$termtax";
                        print OUT_NOROUND "$taxastring\t";
                        print OUT_NOROUND "$termtax";
                        foreach my $j (1..$#file_sample_headers)
                            {   my $sample = $file_sample_headers[$j];
                                chomp($sample);
                                my $number;
                                if ($sample_totals_Prokaryotes[$j] == 0)
                                    {   $number = 0;}
                                else
                                    {   $number = 100 * ($orderedDB{$i}{$sample} / $sample_totals_Prokaryotes[$j]);
                                    }
                                my $roundnum = &ROUND($number,2);
                                if ($roundnum == 0 && $number > 0)
                                    {   print OUT "\tD";
                                    }
                                else
                                    {   print OUT "\t$roundnum";   
                                    }
                                print OUT_NOROUND "\t$number";
                            }
                        print OUT "\n";
                        print OUT_NOROUND "\n";
                    }
            }
    }

$previousPhyla = "";
foreach my $i (sort keys %orderedDB)
    {   my @taxa_array = split(';', $i);
        my $taxastring = $i;
        chomp($taxastring);
        my $phylum = $taxa_array[1]; chomp($phylum);
        my $class = $taxa_array[2]; chomp($class);
        my $order = $taxa_array[3]; chomp($order);
        my $family = $taxa_array[4]; chomp($family);
        my $lastReal;
        foreach my $j (reverse 0..$#taxa_array)
            {   if ($taxa_array[$j] ne "NA")
                    {   $lastReal = $j;
                        last;
                    }
            }
        my $termtax = $taxa_array[$lastReal]; chomp($termtax);
        if ($taxa_array[0] ne "Bacteria" && $taxa_array[0] ne "Archaea" && $taxa_array[0] ne "Eukaryota" && $taxa_array[0] ne "Unknown")
            {   if ($phylum eq $previousPhyla)
                    {   print OUT "$taxastring\t";
                        print OUT "$termtax";
                        print OUT_NOROUND "$taxastring\t";
                        print OUT_NOROUND "$termtax";
                        foreach my $j (1..$#file_sample_headers)
                            {   my $sample = $file_sample_headers[$j];
                                chomp($sample);
                                my $number = 100 * ($orderedDB{$i}{$sample} / $sample_totals[$j]);
                                my $roundnum = &ROUND($number,2);
                                if ($roundnum == 0 && $number > 0)
                                    {   print OUT "\tD";
                                    }
                                else
                                    {   print OUT "\t$roundnum";   
                                    }
                                print OUT_NOROUND "\t$number";  
                            }
                        print OUT "\n";
                        print OUT_NOROUND "\n";
                    }
                else
                    {   $previousPhyla = $phylum;
                        if ($phylum eq "NA")
                            {   print OUT "\t<<<$taxa_array[0], Unknown Phylum (% out of total)>>>\n";
                                print OUT_NOROUND "\t<<<$taxa_array[0], Unknown Phylum (% out of total)>>>\n";
                            }
                        else
                            {   print OUT "\t<<<$taxa_array[0], $phylum (% out of total)>>>\n";
                                print OUT_NOROUND "\t<<<$taxa_array[0], $phylum (% out of total)>>>\n";
                            }
                        print OUT "$taxastring\t";
                        print OUT "$termtax";
                        print OUT_NOROUND "$taxastring\t";
                        print OUT_NOROUND "$termtax";
                        foreach my $j (1..$#file_sample_headers)
                            {   my $sample = $file_sample_headers[$j];
                                chomp($sample);
                                my $number = 100 * ($orderedDB{$i}{$sample} / $sample_totals[$j]);
                                my $roundnum = &ROUND($number,2);
                                if ($roundnum == 0 && $number > 0)
                                    {   print OUT "\tD";
                                    }
                                else
                                    {   print OUT "\t$roundnum";   
                                    }
                                print OUT_NOROUND "\t$number";
                            }
                        print OUT "\n";
                        print OUT_NOROUND "\n";
                    }
            }
    }


close(OUT);
close(OUT_NOROUND);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - S U B R O U T I N E S - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

sub ROUND
{	my $num = shift(@_);
	my $decimals = shift(@_);
	my $roundnum = int(($num * 10**$decimals) + 0.5)/(10**$decimals);
	return $roundnum;
}

