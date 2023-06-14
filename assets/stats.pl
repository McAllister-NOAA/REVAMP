#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use List::MoreUtils qw(uniq);

# - - - - - H E A D E R - - - - - - - - - - - - - - - - -
#Goals of script:
#Calculate some basic taxonomy/ASV stats per sample
#The calculated numbers per taxa represent the number of unique taxa at that level.
#   NOT the number of the taxa only resolved to that level.

# - - - - - C O M M A N D    L I N E    O P T I O N S - - - - - - - -
my %options=();
getopts("a:t:i:h", \%options);

if ($options{h})
    {   print "\n\nHelp called:\nOptions:\n";
        print "-a = ASV count file\n";
        print "-t = ASV collapsed on taxonomy count file\n";
        print "-i = ASV taxonomy table\n";
        print "-h = This help message\n\n";
        die;
    }

# - - - - - G L O B A L  V A R I A B L E S  - - - - - -
my %ASVs;
my %ASV_collapsed;
my $total_asv_count = 0;
my $total_taxaASV_count = 0;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - M A I N - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
open(IN, "<$options{a}") or die "\n\nNADA $options{a} file!!\n\n";
my @data = <IN>; close(IN);
my $sampleheaders = shift(@data); chomp($sampleheaders);
my @sampleheads = split('\t', $sampleheaders);

foreach my $line (@data)
    {   $total_asv_count += 1;
        chomp($line);
        my @split_line = split('\t', $line);
        my $asv = $split_line[0];
        foreach my $i (1..$#split_line)
            {   my $sample = $sampleheads[$i];
                $ASVs{$asv}{$sample} = $split_line[$i];
            }
    }

open(IN2, "<$options{t}") or die "\n\nNADA $options{t} file!!\n\n";
@data = <IN2>; close(IN2);
my $sampleheaders_coll = shift(@data); chomp($sampleheaders_coll);
my @sampleheads_coll = split('\t', $sampleheaders_coll);

foreach my $line (@data)
    {   $total_taxaASV_count += 1;
        chomp($line);
        my @split_line = split('\t', $line);
        my $asv = $split_line[0];
        foreach my $i (1..$#split_line)
            {   my $sample = $sampleheads_coll[$i];
                $ASV_collapsed{$asv}{$sample} = $split_line[$i];
            }
    }

foreach my $i (1..$#sampleheads)
    {   if ($sampleheads[$i] ne $sampleheads_coll[$i])
            {   die "\n\nSample header order does not match between ASV and collapsed ASV file\n\n";
            }
    }

open(IN3, "<$options{i}") or die "\n\nNADA $options{i} file!!\n\n";
@data = <IN3>; close(IN3); shift(@data);

foreach my $line (@data)
    {   chomp($line);
        my @split_line = split('\t', $line);
        my $asv = $split_line[0];
        $ASVs{$asv}{'kingdom'} = $split_line[1];
        $ASVs{$asv}{'phylum'} = $split_line[2];
        $ASVs{$asv}{'class'} = $split_line[3];
        $ASVs{$asv}{'order'} = $split_line[4];
        $ASVs{$asv}{'family'} = $split_line[5];
        $ASVs{$asv}{'genus'} = $split_line[6];
        $ASVs{$asv}{'species'} = $split_line[7];
        
        if (exists $ASV_collapsed{$asv})
        {   $ASV_collapsed{$asv}{'kingdom'} = $split_line[1];
            $ASV_collapsed{$asv}{'phylum'} = $split_line[2];
            $ASV_collapsed{$asv}{'class'} = $split_line[3];
            $ASV_collapsed{$asv}{'order'} = $split_line[4];
            $ASV_collapsed{$asv}{'family'} = $split_line[5];
            $ASV_collapsed{$asv}{'genus'} = $split_line[6];
            $ASV_collapsed{$asv}{'species'} = $split_line[7];
        }
    }

print "Summary stat";

foreach my $i (1..$#sampleheads)
    {   print "\t$sampleheads[$i]";
    }
print "\tAll_Samples\n";

my @total_asv_by_sample; push(@total_asv_by_sample, "N");
my @total_taxaASV_by_sample; push(@total_taxaASV_by_sample, "N");
foreach my $i (1..$#sampleheads)
    {   my $sampletotal = 0;
        my $sample = $sampleheads[$i];
        foreach my $j (sort keys %ASVs)
            {   if ($ASVs{$j}{$sample} > 0)
                    {   $sampletotal += 1;   
                    }
            }
        push(@total_asv_by_sample, $sampletotal);
        $sampletotal = 0;
        foreach my $j (sort keys %ASV_collapsed)
            {   if ($ASV_collapsed{$j}{$sample} > 0)
                    {   $sampletotal += 1;  
                    }
            }
        push(@total_taxaASV_by_sample, $sampletotal);
    }

print "Total ASVs";
foreach my $i (1..$#sampleheads)
    {   print "\t$total_asv_by_sample[$i]";
    }
print "\t$total_asv_count\n";

print "Total ASVs, unique taxonomy";
foreach my $i (1..$#sampleheads)
    {   print "\t$total_taxaASV_by_sample[$i]";
    }
print "\t$total_taxaASV_count\n";


my $count_known = 0;
my $count_unknown = 0;
foreach my $i (sort keys %ASVs)
    {   if ($ASVs{$i}{'kingdom'} eq "Unknown" || $ASVs{$i}{'kingdom'} eq "Environmental Unknown")
            {$count_unknown += 1;}
        else
            {$count_known += 1;}
    }
my $perc_unrounded = (100 * $count_known)/($count_unknown + $count_known);
my $perc_total_rounded = &ROUND($perc_unrounded,2);

my @perc_known_bysample; push(@perc_known_bysample, "N");
foreach my $i (1..$#sampleheads)
    {   my $sample_known = 0;
        my $sample_unknown = 0;
        my $sample = $sampleheads[$i];
        foreach my $j (sort keys %ASVs)
            {   if ($ASVs{$j}{$sample} > 0)
                    {   if ($ASVs{$j}{'kingdom'} eq "Unknown" || $ASVs{$j}{'kingdom'} eq "Environmental Unknown")
                            {$sample_unknown += 1;}
                        else
                            {$sample_known += 1;}
                    }
            }
        my $percknown_samp = (100 * $sample_known)/($sample_unknown + $sample_known);
        my $rounded = &ROUND($percknown_samp,2);
        push(@perc_known_bysample, $rounded)
    }

print "Percent ASVs with Taxonomic Classification";
foreach my $i (1..$#sampleheads)
    {   print "\t$perc_known_bysample[$i]";        
    }
print "\t$perc_total_rounded\n";



my @total_phylum_bysample; push(@total_phylum_bysample, "N");
foreach my $i (1..$#sampleheads)
    {   my @sampletotal;
        my $sample = $sampleheads[$i];
        foreach my $j (sort keys %ASV_collapsed)
            {   if ($ASV_collapsed{$j}{$sample} > 0)
                    {   if ( $ASV_collapsed{$j}{'phylum'} ne "NA")
                            {   push(@sampletotal, $ASV_collapsed{$j}{'phylum'});
                            }
                    }
            }
        my @uniq_sample = uniq @sampletotal;
        my $total = $#uniq_sample + 1;
        push(@total_phylum_bysample, $total);
    }

my @phylum_sum;
foreach my $i (sort keys %ASV_collapsed)
    {   if ( $ASV_collapsed{$i}{'phylum'} ne "NA")
            {   push(@phylum_sum, $ASV_collapsed{$i}{'phylum'});
            }
    }
my @uniq_phylum_sum = uniq @phylum_sum;
my $total_phylum_count = $#uniq_phylum_sum + 1;

print "Phylum Hits";
foreach my $i (1..$#sampleheads)
    {   print "\t$total_phylum_bysample[$i]";        
    }
print "\t$total_phylum_count\n";

my @total_class_bysample; push(@total_class_bysample, "N");
foreach my $i (1..$#sampleheads)
    {   my @sampletotal;
        my $sample = $sampleheads[$i];
        foreach my $j (sort keys %ASV_collapsed)
            {   if ($ASV_collapsed{$j}{$sample} > 0)
                    {   if ( $ASV_collapsed{$j}{'class'} ne "NA")
                            {   push(@sampletotal, $ASV_collapsed{$j}{'class'});
                            }
                    }
            }
        my @uniq_sample = uniq @sampletotal;
        my $total = $#uniq_sample + 1;
        push(@total_class_bysample, $total);
    }

my @class_sum;
foreach my $i (sort keys %ASV_collapsed)
    {   if ( $ASV_collapsed{$i}{'class'} ne "NA")
            {   push(@class_sum, $ASV_collapsed{$i}{'class'});
            }
    }
my @uniq_class_sum = uniq @class_sum;
my $total_class_count = $#uniq_class_sum + 1;

print "Class Hits";
foreach my $i (1..$#sampleheads)
    {   print "\t$total_class_bysample[$i]";        
    }
print "\t$total_class_count\n";

my @total_order_bysample; push(@total_order_bysample, "N");
foreach my $i (1..$#sampleheads)
    {   my @sampletotal;
        my $sample = $sampleheads[$i];
        foreach my $j (sort keys %ASV_collapsed)
            {   if ($ASV_collapsed{$j}{$sample} > 0)
                    {   if ( $ASV_collapsed{$j}{'order'} ne "NA")
                            {   push(@sampletotal, $ASV_collapsed{$j}{'order'});
                            }
                    }
            }
        my @uniq_sample = uniq @sampletotal;
        my $total = $#uniq_sample + 1;
        push(@total_order_bysample, $total);
    }

my @order_sum;
foreach my $i (sort keys %ASV_collapsed)
    {   if ( $ASV_collapsed{$i}{'order'} ne "NA")
            {   push(@order_sum, $ASV_collapsed{$i}{'order'});
            }
    }
my @uniq_order_sum = uniq @order_sum;
my $total_order_count = $#uniq_order_sum + 1;

print "Order Hits";
foreach my $i (1..$#sampleheads)
    {   print "\t$total_order_bysample[$i]";        
    }
print "\t$total_order_count\n";

my @total_family_bysample; push(@total_family_bysample, "N");
foreach my $i (1..$#sampleheads)
    {   my @sampletotal;
        my $sample = $sampleheads[$i];
        foreach my $j (sort keys %ASV_collapsed)
            {   if ($ASV_collapsed{$j}{$sample} > 0)
                    {   if ( $ASV_collapsed{$j}{'family'} ne "NA")
                            {   push(@sampletotal, $ASV_collapsed{$j}{'family'});
                            }
                    }
            }
        my @uniq_sample = uniq @sampletotal;
        my $total = $#uniq_sample + 1;
        push(@total_family_bysample, $total);
    }

my @family_sum;
foreach my $i (sort keys %ASV_collapsed)
    {   if ( $ASV_collapsed{$i}{'family'} ne "NA")
            {   push(@family_sum, $ASV_collapsed{$i}{'family'});
            }
    }
my @uniq_family_sum = uniq @family_sum;
my $total_family_count = $#uniq_family_sum + 1;

print "Family Hits";
foreach my $i (1..$#sampleheads)
    {   print "\t$total_family_bysample[$i]";        
    }
print "\t$total_family_count\n";

my @total_genus_bysample; push(@total_genus_bysample, "N");
foreach my $i (1..$#sampleheads)
    {   my @sampletotal;
        my $sample = $sampleheads[$i];
        foreach my $j (sort keys %ASV_collapsed)
            {   if ($ASV_collapsed{$j}{$sample} > 0)
                    {   if ( $ASV_collapsed{$j}{'genus'} ne "NA")
                            {   push(@sampletotal, $ASV_collapsed{$j}{'genus'});
                            }
                    }
            }
        my @uniq_sample = uniq @sampletotal;
        my $total = $#uniq_sample + 1;
        push(@total_genus_bysample, $total);
    }

my @genus_sum;
foreach my $i (sort keys %ASV_collapsed)
    {   if ( $ASV_collapsed{$i}{'genus'} ne "NA")
            {   push(@genus_sum, $ASV_collapsed{$i}{'genus'});
            }
    }
my @uniq_genus_sum = uniq @genus_sum;
my $total_genus_count = $#uniq_genus_sum + 1;

print "Genus Hits";
foreach my $i (1..$#sampleheads)
    {   print "\t$total_genus_bysample[$i]";        
    }
print "\t$total_genus_count\n";

my @total_species_bysample; push(@total_species_bysample, "N");
foreach my $i (1..$#sampleheads)
    {   my @sampletotal;
        my $sample = $sampleheads[$i];
        foreach my $j (sort keys %ASV_collapsed)
            {   if ($ASV_collapsed{$j}{$sample} > 0)
                    {   if ( $ASV_collapsed{$j}{'species'} ne "NA")
                            {   push(@sampletotal, $ASV_collapsed{$j}{'species'});
                            }
                    }
            }
        my @uniq_sample = uniq @sampletotal;
        my $total = $#uniq_sample + 1;
        push(@total_species_bysample, $total);
    }

my @species_sum;
foreach my $i (sort keys %ASV_collapsed)
    {   if ( $ASV_collapsed{$i}{'species'} ne "NA")
            {   push(@species_sum, $ASV_collapsed{$i}{'species'});
            }
    }
my @uniq_species_sum = uniq @species_sum;
my $total_species_count = $#uniq_species_sum + 1;

print "Species Hits";
foreach my $i (1..$#sampleheads)
    {   print "\t$total_species_bysample[$i]";        
    }
print "\t$total_species_count\n";


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - S U B R O U T I N E S - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub ROUND
{   my $num = shift(@_);
    my $decimals = shift(@_);
    my $roundnum = int(($num * 10**$decimals) + 0.5)/(10**$decimals);
    return $roundnum;
}

# - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -
