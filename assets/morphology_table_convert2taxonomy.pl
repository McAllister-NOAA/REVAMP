#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use List::MoreUtils qw(uniq);

$|=1;

# - - - - - H E A D E R - - - - - - - - - - - - - - - - -
#Goals of script:
#Convert morphology data so it can be compared with markers.


# - - - - - C O M M A N D    L I N E    O P T I O N S - - - - - - - -
my %options=();
getopts("i:o:m:f:ah", \%options);

if ($options{h})
    {   print "\n\nHelp called:\nOptions:\n";
        print "-i = Input tab delimited morphology spreadsheet\n";
        print "-a = Set to automatically fill in taxonkit output (as done in REVAMP)\n";
        print "-m = Location of REVAMP directory\n";
        print "-f = Filter cutoff for zzOther\n";
        print "-o = Output directory\n";
        print "-h = This help message\n\n";
        die;
    }

# - - - - - G L O B A L  V A R I A B L E S  - - - - - -
my %ASV;
my @sample_headers;
my @taxa_list;
my $asv_count = 1;
my %TAXLIST;
my %COLLAPSE;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - M A I N - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
system ("mkdir -p ".$options{o});

#Cleanup morphology table
system("cat ".$options{i}.' | awk '."'".'FNR==NR {gsub("[^a-zA-Z0-9]", "_", $1)} 1'."'".' OFS="\t" FS="\t" | sed -e '."'".'2,$ s/^/MP_/'."'".' | grep -v "^MP_$" | sed "s/$(printf '."'".'\r'."'".')\$//" > '.$options{o}.'/morphology_data_cleaned.txt');

#IMPORT morphology table info
open(IN1, "<$options{o}"."/morphology_data_cleaned.txt") or die "\n\nThere is no $options{o}"."/morphology_data_cleaned.txt file!!\n\n";
my @data = <IN1>; close(IN1);
my $header = shift(@data); chomp($header);
my @headers = split('\t', $header);

my @temp_headers;
foreach my $i (0..$#headers)
    {   my $h = $headers[$i];
        $h =~ s/[^a-zA-Z0-9]/_/g;
        push(@temp_headers, $h);
    }
@headers = @temp_headers;


foreach my $line (@data)
	{	chomp($line);
		my @info = split('\t', $line);
        my $asv_current = "ASV_".$asv_count;
        $ASV{$asv_current}{"sample"} = $info[0];
        if (exists $info[1])
            {   if (lc($info[1]) eq "unknown")
                    {   $ASV{$asv_current}{"taxonomy"} = "UNKNOWN";
                    }
                else
                    {   my $taxo = $info[1];
                        if ($taxo =~ m/ Sm$/ || $taxo =~ m/ Lg$/)
                            {   $taxo =~ m/(.+) [LS][gm]$/;
                                my $newtaxa = $1;
                                $taxo = $newtaxa;
                            }
                        $ASV{$asv_current}{"taxonomy"} = uc($taxo);
                    }
            }
        else
            {   $ASV{$asv_current}{"taxonomy"} = "UNKNOWN";
            }
		foreach my $i (2..$#headers)
            {   if (exists $info[$i])
                    {   if ($info[$i] =~ m/[0-9]+/)
                            {$ASV{$asv_current}{$headers[$i]} = $info[$i];}
                        else {$ASV{$asv_current}{$headers[$i]} = 0;}
                    }
                else {$ASV{$asv_current}{$headers[$i]} = 0;}
            }
        $asv_count += 1;
	}

foreach my $i (sort keys %ASV)
    {   push(@sample_headers, $ASV{$i}{"sample"});
        unless ($ASV{$i}{"taxonomy"} eq "UNKNOWN")
            {push(@taxa_list, $ASV{$i}{"taxonomy"});}
    }
    
my @unique_sampleHeaders = uniq @sample_headers;
my @unique_taxonomy = uniq @taxa_list;

open(TAXAOUT, ">".$options{o}."/morphology_unique_terminaltaxa.txt");
foreach my $i (0..$#unique_taxonomy-1)
    {   print TAXAOUT "$unique_taxonomy[$i]\n";
    }
print TAXAOUT "$unique_taxonomy[$#unique_taxonomy]";
close(TAXAOUT);
    
system("taxonkit name2taxid ".$options{o}."/morphology_unique_terminaltaxa.txt > ".$options{o}."/morphology_unique_terminaltaxa_w_taxids_round1.txt");

open(IN2, "<$options{o}"."/morphology_unique_terminaltaxa_w_taxids_round1.txt") or die "\n\nThere is no $options{o}"."/morphology_unique_terminaltaxa_w_taxids_round1.txt file!!\n\n";
my @data2 = <IN2>; close(IN2);
my @multi_taxID;
foreach my $line (@data2)
    {   chomp($line);
		my @info = split('\t', $line);
        unless ($info[0] eq "" || $info[0] eq "\n" || $info[0] eq "\t")
        {   if (exists $TAXLIST{$info[0]})
                {   push(@multi_taxID, $info[0]);
                }
            else
                {   if (exists $info[1])
                        {   $TAXLIST{$info[0]}{"taxID"} = $info[1];
                            $TAXLIST{$info[0]}{"run"} = "primary";
                        }
                    else
                        {   $TAXLIST{$info[0]}{"taxID"} = "UNMATCHED";
                            $TAXLIST{$info[0]}{"run"} = "unmatched";
                        }
                }
        }
    }
    
open(TAXAOUT5, ">".$options{o}."/morphology_finalTaxaNameMatch2TaxonomyAssignments.txt");

if (scalar @multi_taxID > 0)
    {   my $printnum = $#multi_taxID + 1;
        print "\n\nTaxa name(s) exist(s) as multiple potential taxids (n = $printnum)\n";
        print TAXAOUT5 "Manual choice decisions:";
        print TAXAOUT5 "\n\nTaxa name(s) exist(s) as multiple potential taxids (n = $printnum)\n";
        foreach my $i (@multi_taxID)
            {   print "\nReview Taxa <<$i>> and choose the appropriate TAXID or enter UNMATCHED:\n";
                print TAXAOUT5 "\nReview Taxa <<$i>> and choose the appropriate TAXID or enter UNMATCHED:\n";
                foreach my $line (@data2)
                    {   chomp($line);
                        my @info = split('\t', $line);
                        if ($info[0] eq $i)
                            {   print "$i\t$info[1]\n";
                                print TAXAOUT5 "$i\t$info[1]\n";
                            }
                    }
                print "Choose:\n";
                print TAXAOUT5 "Chose:\n";
                my $chosen_taxid = <STDIN>;
                chomp($chosen_taxid);
                print TAXAOUT5 "$chosen_taxid\n";
                if ($chosen_taxid eq "UNMATCHED")
                    {   $TAXLIST{$i}{"taxID"} = "UNMATCHED";
                        $TAXLIST{$i}{"run"} = "manual_choice";
                    }
                else
                    {   $TAXLIST{$i}{"taxID"} = $chosen_taxid;
                        $TAXLIST{$i}{"run"} = "primary_manual_choice";
                    }                
            }
    }

my @unmatched_taxa;
foreach my $i (sort keys %TAXLIST)
    {   if ($TAXLIST{$i}{"taxID"} eq "UNMATCHED")
            {   unless ($TAXLIST{$i}{"run"} eq "manual_choice")
                    {push(@unmatched_taxa, $i);}
            }
    }

my @rerun_taxa;
my @rerun_taxa_keys;
foreach my $i (@unmatched_taxa)
    {   if ($i =~ m/[a-zA-Z0-9]+ [a-zA-Z0-9]+/)
            {   $i =~ m/([a-zA-Z0-9]+) [a-zA-Z0-9]+/;
                my $newtaxa = $1;
                push(@rerun_taxa, $newtaxa);
                push(@rerun_taxa_keys, $i);
            }
    }

my @unique_rerun_taxa = uniq @rerun_taxa;
    
open(TAXAOUT2, ">".$options{o}."/morphology_rerun_terminaltaxa.txt");
foreach my $i (0..$#unique_rerun_taxa-1)
    {   print TAXAOUT2 "$unique_rerun_taxa[$i]\n";
    }
print TAXAOUT2 "$unique_rerun_taxa[$#unique_rerun_taxa]";
close(TAXAOUT2);

system("taxonkit name2taxid ".$options{o}."/morphology_rerun_terminaltaxa.txt > ".$options{o}."/morphology_rerun_terminaltaxa_w_taxids_round2.txt");

open(IN3, "<$options{o}"."/morphology_rerun_terminaltaxa_w_taxids_round2.txt") or die "\n\nThere is no $options{o}"."/morphology_rerun_terminaltaxa_w_taxids_round2.txt file!!\n\n";
my @data3 = <IN3>; close(IN3);
my @multi_taxID_round2;
foreach my $line (@data3)
    {   chomp($line);
		my @info = split('\t', $line);
        unless ($info[0] eq "" || $info[0] eq "\n" || $info[0] eq "\t")
        {   if (exists $info[1])
                {   foreach my $j (@rerun_taxa_keys)
                        {   $j =~ m/([a-zA-Z0-9]+) [a-zA-Z0-9]+/;
                            my $testname = $1;
                            if ($testname eq $info[0])
                                {   if ($TAXLIST{$j}{"run"} eq "genus")
                                        {   push(@multi_taxID_round2, $j);
                                        }
                                    else
                                        {   $TAXLIST{$j}{"taxID"} = $info[1];
                                            $TAXLIST{$j}{"run"} = "genus";   
                                        }
                                }   
                        }
                }
        }
    }
    
if (scalar @multi_taxID_round2 > 0)
    {   my $printnum = $#multi_taxID_round2 + 1;
        print "\n\nTaxa name(s) exist(s) as multiple potential taxids in Round2 (n = $printnum)\n";
        print TAXAOUT5 "\n\nTaxa name(s) exist(s) as multiple potential taxids in Round2 (n = $printnum)\n";
        foreach my $i (@multi_taxID_round2)
            {   print "\nReview Taxa <<$i>> and choose the appropriate TAXID or enter UNMATCHED:\n";
                print TAXAOUT5 "\nReview Taxa <<$i>> and choose the appropriate TAXID or enter UNMATCHED:\n";
                foreach my $line (@data3)
                    {   chomp($line);
                        my @info = split('\t', $line);
                        if ($info[0] eq $i)
                            {   print "$i\t$info[1]\n";
                                print TAXAOUT5 "$i\t$info[1]\n";
                            }
                    }
                print "Choose:\n";
                print TAXAOUT5 "Chose:\n";
                my $chosen_taxid = <STDIN>;
                chomp($chosen_taxid);
                print TAXAOUT5 "$chosen_taxid\n";
                if ($chosen_taxid eq "UNMATCHED")
                    {   $TAXLIST{$i}{"taxID"} = "UNMATCHED";
                        $TAXLIST{$i}{"run"} = "manual_choice";
                    }
                else
                    {   $TAXLIST{$i}{"taxID"} = $chosen_taxid;
                        $TAXLIST{$i}{"run"} = "genus_manual_choice";
                    }                
            }
    }
    
open(TAXAOUT3, ">".$options{o}."/morphology_taxIDs.txt");
foreach my $i (sort keys %TAXLIST)
    {   unless ($TAXLIST{$i}{"taxID"} eq "UNMATCHED")
            {   print TAXAOUT3 $TAXLIST{$i}{"taxID"}."\n";
            }
    }
close(TAXAOUT3);

system("taxonkit lineage ".$options{o}."/morphology_taxIDs.txt | awk '".'$2!=""'."' > ".$options{o}."/morphology_taxonkit_out.txt");
system("taxonkit reformat ".$options{o}."/morphology_taxonkit_out.txt > ".$options{o}."/morphology_reformatted_taxonkit_out.txt");

if ($options{a})
    {   system('perl '.$options{m}.'/assets/fillIn_taxonkit_mod.pl -i '.$options{o}.'/morphology_reformatted_taxonkit_out.txt > '.$options{o}.'/morphology_reformatted_taxonkit_out.txt_temp');
        system('mv '.$options{o}.'/morphology_reformatted_taxonkit_out.txt '.$options{o}.'/morphology_reformatted_taxonkit_out_ORIGINAL.txt');
        system('mv '.$options{o}.'/morphology_reformatted_taxonkit_out.txt_temp '.$options{o}.'/morphology_reformatted_taxonkit_out.txt');
    }
else
    {   print "\n\nMake your own adjustments to the morphology_reformatted_taxonkit_out.txt file and press enter (or rerun with automated option)\n";
        my $pause = <STDIN>;
    }

#foreach my $i (sort keys %TAXLIST)
#    {print "$i\t".$TAXLIST{$i}{"taxID"}."\t".$TAXLIST{$i}{"run"}."\n";}
open(TAXAOUT4, ">".$options{o}."/morphology_taxaWarning_synonym_or_mismatch.txt");

open(IN4, "<$options{o}"."/morphology_reformatted_taxonkit_out.txt") or die "\n\nThere is no $options{o}"."/morphology_reformatted_taxonkit_out.txt file!!\n\n";
my @data4 = <IN4>; close(IN4);
foreach my $line (@data4)
    {   chomp($line);
		my @info = split('\t', $line);
        my $taxid = $info[0];
        my $taxstring = $info[2];
        my $comparatorstring = $info[1];
        my $kingdom = "NA";
        my $phylumn = "NA";
        my $class = "NA";
        my $order = "NA";
        my $family = "NA";
        my $genus = "NA";
        my $species = "NA";
        my @taxstringarray = split(';',$taxstring);
        my $finalhit;
        if (exists $taxstringarray[0]) {$kingdom = $taxstringarray[0];
                                        $finalhit = 0;}
        if (exists $taxstringarray[1]) {$phylumn = $taxstringarray[1];
                                        $finalhit = 1;}
        if (exists $taxstringarray[2]) {$class = $taxstringarray[2];
                                        $finalhit = 2;}
        if (exists $taxstringarray[3]) {$order = $taxstringarray[3];
                                        $finalhit = 3;}
        if (exists $taxstringarray[4]) {$family = $taxstringarray[4];
                                        $finalhit = 4;}
        if (exists $taxstringarray[5]) {$genus = $taxstringarray[5];
                                        $finalhit = 5;}
        if (exists $taxstringarray[6]) {$species = $taxstringarray[6];
                                        $finalhit = 6;}
        
        my @compareString = split(';', $comparatorstring);
        my @additional_tax;
        my $store_tax = "FALSE";
        foreach my $k (0..$#compareString)
            {   my $compare = $compareString[$k];
                if ($store_tax eq "TRUE")
                    {   push(@additional_tax, $compare);
                    }
                if ($store_tax eq "FALSE")
                    {   if ($compare eq $taxstringarray[$finalhit])
                            {   $store_tax = "TRUE";
                            }
                    }
            }
        if (scalar @additional_tax > 0)
            {   my $final_add = $additional_tax[$#additional_tax];
                if ($finalhit == 0)
                    {   $phylumn = $final_add."__p";}
                if ($finalhit == 1)
                    {   $class = $final_add."__c";}
                if ($finalhit == 2)
                    {   $order = $final_add."__o";}
                if ($finalhit == 3)
                    {   $family = $final_add."__f";}
                if ($finalhit == 4)
                    {   $genus = $final_add."__g";}
                if ($finalhit == 5)
                    {   $species = $final_add."__s";}
            }
       
        foreach my $i (sort keys %TAXLIST)
            {   if ($TAXLIST{$i}{"taxID"} eq $taxid)
                    {   if ($TAXLIST{$i}{"run"} eq "primary" || $TAXLIST{$i}{"run"} eq "primary_manual_choice")
                            {   $TAXLIST{$i}{"cleantaxa"} = $kingdom.";".$phylumn.";".$class.";".$order.";".$family.";".$genus.";".$species;
                            }
                        if ($TAXLIST{$i}{"run"} eq "genus" || $TAXLIST{$i}{"run"} eq "genus_manual_choice")
                            {   my $speciesname = lc($i);
                                $speciesname = ucfirst($speciesname);
                                $TAXLIST{$i}{"cleantaxa"} = $kingdom.";".$phylumn.";".$class.";".$order.";".$family.";".$genus.";".$speciesname;
                            }   
                    }
                if ($TAXLIST{$i}{"taxID"} eq $taxid && uc($taxstringarray[$finalhit]) ne uc($i) && $TAXLIST{$i}{"cleantaxa"} !~ m/__/ && $TAXLIST{$i}{"run"} ne "genus" && $TAXLIST{$i}{"run"} ne "genus_manual_choice")
                    {   print TAXAOUT4 "$i\t".$TAXLIST{$i}{"taxID"}."\t".$TAXLIST{$i}{"run"}."\t".$TAXLIST{$i}{"cleantaxa"}."\n";
                    }
                if ($TAXLIST{$i}{"taxID"} eq $taxid && $TAXLIST{$i}{"cleantaxa"} !~ m/__/)
                    { if ($TAXLIST{$i}{"run"} eq "genus" || $TAXLIST{$i}{"run"} eq "genus_manual_choice")
                        {   my $taxname = $i;
                            $i =~ m/([a-zA-Z0-9]+) [a-zA-Z0-9]+/;
                            my $testname = $1;
                            if (uc($taxstringarray[$finalhit]) ne uc($testname))
                                {   print TAXAOUT4 "$i\t".$TAXLIST{$i}{"taxID"}."\t".$TAXLIST{$i}{"run"}."\t".$TAXLIST{$i}{"cleantaxa"}."\n";
                                }
                        }
                    }
            }
    }
close(TAXAOUT4);

foreach my $i (sort keys %TAXLIST)
    {   unless (exists $TAXLIST{$i}{"cleantaxa"})
            {   my $name = lc($i);
                $name = ucfirst($name);
                $TAXLIST{$i}{"cleantaxa"} = $name.";NA;NA;NA;NA;NA;NA";
            }
    }

print TAXAOUT5 "\n\nTaxa names to taxonomy results:\n";
print TAXAOUT5 "Taxa_name\tTaxID\tRun_Matched\tClean_Taxonomy\n";
foreach my $i (sort keys %TAXLIST)
    { print TAXAOUT5 "$i\t".$TAXLIST{$i}{"taxID"}."\t".$TAXLIST{$i}{"run"}."\t".$TAXLIST{$i}{"cleantaxa"}."\n";
    }
close(TAXAOUT5);

foreach my $i (sort keys %ASV)
    {   my $test_asv = $ASV{$i}{"taxonomy"};
        if (exists $TAXLIST{$test_asv})
            {   my $clean_taxa = $TAXLIST{$test_asv}{"cleantaxa"};
                $ASV{$i}{"tax_hierarchy"} = $clean_taxa;
            }
        elsif (uc($test_asv) eq "UNKNOWN")
            {   $ASV{$i}{"tax_hierarchy"} = "Unknown;NA;NA;NA;NA;NA;NA";
            }
        else
            {   die "\n\nUnexpected taxa name $test_asv not sent to be examined for taxonomic hierarchy\n\n";
            }
    }       
            
foreach my $i (sort keys %ASV)
    {   my $clean_taxa = $ASV{$i}{"tax_hierarchy"};
        my $sample = $ASV{$i}{"sample"};
        $COLLAPSE{$clean_taxa}{"ASVs"} .= $i.";";
        foreach my $j (2..$#headers)
            {   $COLLAPSE{$clean_taxa}{$sample."---".$headers[$j]} += $ASV{$i}{$headers[$j]};
            }   
    }
    
foreach my $i (sort keys %COLLAPSE)
    {   my $asv_list = $COLLAPSE{$i}{"ASVs"};
        my @asv_array = split(';', $asv_list);
        my $pick = $asv_array[0];
        foreach my $j (@asv_array)
            {   $pick =~ m/ASV_([0-9]+)/;
                my $orig = $1;
                $j =~ m/ASV_([0-9]+)/;
                my $comp = $1;
                if ($comp < $orig)
                    {   $pick = $j;
                    }
            }
        $COLLAPSE{$i}{"chosen_ASV"} = $pick;
    }

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - M A I N - O U T P U T - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#ASVs_counts.tsv = Plain counts for all "ASVs". In this case, each ASV simply represents a unique taxa.
foreach my $unique_measurement (2..$#headers)
    {   system("mkdir -p ".$options{o}."/morphology_REVAMPtables_$headers[$unique_measurement]");
        open(ASV_count, ">".$options{o}."/morphology_REVAMPtables_$headers[$unique_measurement]/ASVs_counts.tsv");
        print ASV_count "x";
        foreach my $i (0..$#unique_sampleHeaders)
            {   print ASV_count "\t$unique_sampleHeaders[$i]";
            }
        print ASV_count "\n";
        foreach my $i (sort keys %COLLAPSE)
            {   my $chosen_asv = $COLLAPSE{$i}{"chosen_ASV"};
                print ASV_count "$chosen_asv";
                foreach my $j (0..$#unique_sampleHeaders)
                    {   if (exists $COLLAPSE{$i}{$unique_sampleHeaders[$j]."---".$headers[$unique_measurement]})
                            {   my $print_value = $COLLAPSE{$i}{$unique_sampleHeaders[$j]."---".$headers[$unique_measurement]};
                                print ASV_count "\t$print_value";
                            }
                        else
                            {   print ASV_count "\t0";
                            }
                    }
                print ASV_count "\n";
            }
        close(ASV_count);
    }

#ASVs_counts_NOUNKNOWNS.tsv = Plain counts for all "ASVs", except those with unknown taxonomy.
foreach my $unique_measurement (2..$#headers)
    {   system("mkdir -p ".$options{o}."/morphology_REVAMPtables_$headers[$unique_measurement]");
        open(ASV_count, ">".$options{o}."/morphology_REVAMPtables_$headers[$unique_measurement]/ASVs_counts_NOUNKNOWNS.tsv");
        print ASV_count "x";
        foreach my $i (0..$#unique_sampleHeaders)
            {   print ASV_count "\t$unique_sampleHeaders[$i]";
            }
        print ASV_count "\n";
        foreach my $i (sort keys %COLLAPSE)
            {   unless ($i eq "Unknown;NA;NA;NA;NA;NA;NA")
                    {   my $chosen_asv = $COLLAPSE{$i}{"chosen_ASV"};
                        print ASV_count "$chosen_asv";
                        foreach my $j (0..$#unique_sampleHeaders)
                            {   if (exists $COLLAPSE{$i}{$unique_sampleHeaders[$j]."---".$headers[$unique_measurement]})
                                    {   my $print_value = $COLLAPSE{$i}{$unique_sampleHeaders[$j]."---".$headers[$unique_measurement]};
                                        print ASV_count "\t$print_value";
                                    }
                                else
                                    {   print ASV_count "\t0";
                                    }
                            }
                        print ASV_count "\n";
                    }
            }
        close(ASV_count);
    }

#morphology_asvTaxonomyTable.txt = ASV to taxonomy.
system("mkdir -p ".$options{o}."/morphology_REVAMPtables_ASV2Taxonomy");
open(ASV_taxa, ">".$options{o}."/morphology_REVAMPtables_ASV2Taxonomy/morphology_asvTaxonomyTable.txt");
print ASV_taxa "ASV\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n";
foreach my $i (sort keys %COLLAPSE)
    {   print ASV_taxa $COLLAPSE{$i}{"chosen_ASV"};
        my $clean_taxa = $i;
        my @clean_taxa_array = split(';', $clean_taxa);
        foreach my $j (0..$#clean_taxa_array)
            {   print ASV_taxa "\t".$clean_taxa_array[$j];
            }
        print ASV_taxa "\n";
    }
close(ASV_taxa);

#morphology_asvTaxonomyTable_NOUNKNOWNS.txt = ASV to taxonomy (without UNKNOWN).
open(ASV_taxa, ">".$options{o}."/morphology_REVAMPtables_ASV2Taxonomy/morphology_asvTaxonomyTable_NOUNKNOWNS.txt");
print ASV_taxa "ASV\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n";
foreach my $i (sort keys %COLLAPSE)
    {   unless ($i eq "Unknown;NA;NA;NA;NA;NA;NA")
            {   print ASV_taxa $COLLAPSE{$i}{"chosen_ASV"};
                my $clean_taxa = $i;
                my @clean_taxa_array = split(';', $clean_taxa);
                foreach my $j (0..$#clean_taxa_array)
                    {   print ASV_taxa "\t".$clean_taxa_array[$j];
                    }
                print ASV_taxa "\n";
            }
    }
close(ASV_taxa);

#KRONA Plots
foreach my $unique_measurement (2..$#headers)
    {   system("mkdir -p ".$options{o}."/morphology_REVAMPtables_$headers[$unique_measurement]/KRONA_plots/KRONA_inputs");
        foreach my $uniquesample (0..$#unique_sampleHeaders)
            {   open(ASV_taxa, ">".$options{o}."/morphology_REVAMPtables_$headers[$unique_measurement]/KRONA_plots/KRONA_inputs/".$unique_sampleHeaders[$uniquesample]."_KRONA.txt");
                foreach my $i (sort keys %COLLAPSE)
                    {   if (exists $COLLAPSE{$i}{$unique_sampleHeaders[$uniquesample]."---".$headers[$unique_measurement]})
                            {   my $print_value = $COLLAPSE{$i}{$unique_sampleHeaders[$uniquesample]."---".$headers[$unique_measurement]};
                                print ASV_taxa "$print_value";
                                my $tax = $i;
                                my @tax_array = split(';', $tax);
                                foreach my $entry (0..$#tax_array)
                                    {   if ($tax_array[$entry] ne "NA")
                                            {   print ASV_taxa "\t$tax_array[$entry]";
                                            }
                                        else
                                            {   last;
                                            }
                                    }
                                print ASV_taxa "\n";
                            }
                    }
                close(ASV_taxa);
            }
            
            open(ASV_taxa, ">".$options{o}."/morphology_REVAMPtables_$headers[$unique_measurement]/KRONA_plots/KRONA_inputs/Morphology_samplesSummedKRONA.txt");
            foreach my $uniquesample (0..$#unique_sampleHeaders)
            {   foreach my $i (sort keys %COLLAPSE)
                    {   if (exists $COLLAPSE{$i}{$unique_sampleHeaders[$uniquesample]."---".$headers[$unique_measurement]})
                            {   my $print_value = $COLLAPSE{$i}{$unique_sampleHeaders[$uniquesample]."---".$headers[$unique_measurement]};
                                print ASV_taxa "$print_value";
                                my $tax = $i;
                                my @tax_array = split(';', $tax);
                                foreach my $entry (0..$#tax_array)
                                    {   if ($tax_array[$entry] ne "NA")
                                            {   print ASV_taxa "\t$tax_array[$entry]";
                                            }
                                        else
                                            {   last;
                                            }
                                    }
                                print ASV_taxa "\n";
                            }
                    }
            }
            close(ASV_taxa);
    }

foreach my $unique_measurement (2..$#headers)
    {   my @locationKRONA;
        foreach my $i (@unique_sampleHeaders)
            {   push(@locationKRONA, $options{o}."/morphology_REVAMPtables_$headers[$unique_measurement]/KRONA_plots/KRONA_inputs/".$i."_KRONA.txt");
            }
        my $printkronasamples = join(' ', @locationKRONA);
        system("ImportText.pl -o ".$options{o}."/morphology_REVAMPtables_".$headers[$unique_measurement]."/KRONA_plots/Morphology_master_krona.html $printkronasamples");
        
        system("ImportText.pl -o ".$options{o}."/morphology_REVAMPtables_".$headers[$unique_measurement]."/KRONA_plots/Morphology_samplesSummedKRONA.html ".$options{o}."/morphology_REVAMPtables_".$headers[$unique_measurement]."/KRONA_plots/KRONA_inputs/Morphology_samplesSummedKRONA.txt");
    }

#ASVs_counts_NOUNKNOWNS_collapsedOnTaxonomy_percentabund.tsv = relative abundance ASV biom table
foreach my $unique_measurement (2..$#headers)
    {   system("mkdir -p ".$options{o}."/morphology_REVAMPtables_$headers[$unique_measurement]");
        open(ASV_count, ">".$options{o}."/morphology_REVAMPtables_$headers[$unique_measurement]/ASVs_counts_NOUNKNOWNS_collapsedOnTaxonomy_percentabund.tsv");
        print ASV_count "x";
        foreach my $i (0..$#unique_sampleHeaders)
            {   print ASV_count "\t$unique_sampleHeaders[$i]";
            }
        print ASV_count "\n";
        
        my @sample_sums;
        foreach my $j (0..$#unique_sampleHeaders)
            {   my $samp_sum = 0;
                foreach my $i (sort keys %COLLAPSE)
                    {   unless ($i eq "Unknown;NA;NA;NA;NA;NA;NA")
                            {   if (exists $COLLAPSE{$i}{$unique_sampleHeaders[$j]."---".$headers[$unique_measurement]})
                                            {   my $value = $COLLAPSE{$i}{$unique_sampleHeaders[$j]."---".$headers[$unique_measurement]};
                                                $samp_sum += $value;
                                            }
                            }
                    }
                push(@sample_sums, $samp_sum);
            }
        
        foreach my $i (sort keys %COLLAPSE)
            {   unless ($i eq "Unknown;NA;NA;NA;NA;NA;NA")
                    {   my $chosen_asv = $COLLAPSE{$i}{"chosen_ASV"};
                        print ASV_count "$chosen_asv";
                        foreach my $j (0..$#unique_sampleHeaders)
                            {   if (exists $COLLAPSE{$i}{$unique_sampleHeaders[$j]."---".$headers[$unique_measurement]})
                                    {   my $value = $COLLAPSE{$i}{$unique_sampleHeaders[$j]."---".$headers[$unique_measurement]};
                                        my $relabund = 100 * $value / $sample_sums[$j];
                                        print ASV_count "\t$relabund";
                                    }
                                else
                                    {   print ASV_count "\t0";
                                    }
                            }
                        print ASV_count "\n";
                    }
            }
        close(ASV_count);
    }
    
#open(SAMP, ">".$options{o}."/sample_order.txt");
#foreach my $i (sort @unique_sampleHeaders)
#    {   print SAMP "$i\n";
#    }
#close(SAMP);

#Run filter_lowabundance_taxa.pl on outfiles to create zzOther taxonomy file
foreach my $unique_measurement (2..$#headers)
    {   system('perl '.$options{m}.'/assets/filter_lowabundance_taxa.pl -a '.$options{o}.'/morphology_REVAMPtables_'.$headers[$unique_measurement].'/ASVs_counts_NOUNKNOWNS_collapsedOnTaxonomy_percentabund.tsv -t '.$options{o}.'/morphology_REVAMPtables_ASV2Taxonomy/morphology_asvTaxonomyTable_NOUNKNOWNS.txt -p '.$options{f}.' > '.$options{o}.'/morphology_REVAMPtables_'.$headers[$unique_measurement].'/ASVTaxonomyTable_NOUNKNOWNS_replaceLowAbund2zzOther.txt');
    }

#Print file with the unique header identifiers
open(HEADER, ">".$options{o}."/unique_measurements_list.txt");
foreach my $unique_measurement (2..$#headers)
    {   print HEADER "$headers[$unique_measurement]\n";
    }
close(HEADER);


#Move KRONA plots to Figures directory
foreach my $unique_measurement (2..$#headers)
    {   system("mkdir -p ".$options{o}."/morphology_REVAMPtables_".$headers[$unique_measurement]."/Figures/00_KRONA_plots");
        system("cp ".$options{o}."/morphology_REVAMPtables_".$headers[$unique_measurement]."/KRONA_plots/Morphology_master_krona.html ".$options{o}."/morphology_REVAMPtables_".$headers[$unique_measurement]."/Figures/00_KRONA_plots/Morphology_".$headers[$unique_measurement]."_master_krona.html");
        system("cp ".$options{o}."/morphology_REVAMPtables_".$headers[$unique_measurement]."/KRONA_plots/Morphology_samplesSummedKRONA.html ".$options{o}."/morphology_REVAMPtables_".$headers[$unique_measurement]."/Figures/00_KRONA_plots/Morphology_".$headers[$unique_measurement]."_samplesSummedKRONA.html");
    }

#Print out ASV table combining asv2Taxonomy with perc abundance for human viewing
foreach my $unique_measurement (2..$#headers)
    {   open(OUT, ">".$options{o}."/morphology_REVAMPtables_".$headers[$unique_measurement]."/taxonomy2PercentAbundance_humanReadable.txt");
        open(OUT_NOROUND, ">".$options{o}."/morphology_REVAMPtables_".$headers[$unique_measurement]."/taxonomy2PercentAbundance_humanReadable_NoRounding.txt");
        
        my @sample_totals;
        my @sample_totals_sansUnknowns;
        my @sample_totals_Eukaryota;
        my @sample_totals_Prokaryotes;
        foreach my $i (0..$#unique_sampleHeaders)
            {   my $sample = $unique_sampleHeaders[$i];
                chomp($sample);
                my $thisSampleTotal = 0;
                my $thisSampleNoUNKNOWN = 0;
                my $thisSampleEuk = 0;
                my $thisSampleProk = 0;
                
                foreach my $j (sort keys %COLLAPSE)
                    {   if (exists $COLLAPSE{$j}{$unique_sampleHeaders[$i]."---".$headers[$unique_measurement]})
                                            {   $thisSampleTotal += $COLLAPSE{$j}{$unique_sampleHeaders[$i]."---".$headers[$unique_measurement]};
                                            }
                        unless ($j =~ m/^Unknown/)
                            {   if (exists $COLLAPSE{$j}{$unique_sampleHeaders[$i]."---".$headers[$unique_measurement]})
                                            {   $thisSampleNoUNKNOWN += $COLLAPSE{$j}{$unique_sampleHeaders[$i]."---".$headers[$unique_measurement]};
                                            }
                            }
                        if ($j =~ m/^Eukaryota/)
                            {   if (exists $COLLAPSE{$j}{$unique_sampleHeaders[$i]."---".$headers[$unique_measurement]})
                                            {   $thisSampleEuk += $COLLAPSE{$j}{$unique_sampleHeaders[$i]."---".$headers[$unique_measurement]};
                                            }
                            }
                        if ($j =~ m/^Bacteria/ || $j =~ m/^Archaea/)
                            {   if (exists $COLLAPSE{$j}{$unique_sampleHeaders[$i]."---".$headers[$unique_measurement]})
                                            {   $thisSampleProk += $COLLAPSE{$j}{$unique_sampleHeaders[$i]."---".$headers[$unique_measurement]};
                                            }
                            }
                    }
                push(@sample_totals, $thisSampleTotal);
                push(@sample_totals_sansUnknowns, $thisSampleNoUNKNOWN);
                push(@sample_totals_Eukaryota, $thisSampleEuk);
                push(@sample_totals_Prokaryotes, $thisSampleProk);
            }
        
        my %orderedDB;
        foreach my $i (sort keys %COLLAPSE)
            {   my $taxa = $i; chomp($i);
                $orderedDB{$taxa}{'asv_key'} = $COLLAPSE{$i}{"chosen_ASV"};
                foreach my $j (0..$#unique_sampleHeaders)
                    {   my $sample = $unique_sampleHeaders[$j];
                        chomp($sample);
                        my $value;
                        if (exists $COLLAPSE{$i}{$unique_sampleHeaders[$j]."---".$headers[$unique_measurement]})
                            {   $value = $COLLAPSE{$i}{$unique_sampleHeaders[$j]."---".$headers[$unique_measurement]};
                            }
                        else
                            {   $value = 0;
                            }
                        $orderedDB{$taxa}{$sample} = $value;
                    }
            }
        
        print OUT "All percentages range 0-100% of the total count for either Eukaryotes or Prokaryotes, with exception of the first three data rows (% Known, % Eukaryote, % Prokaryote), which represent %s of total reads.\n";
        print OUT "Taxonomy String (K->S)\t";
        print OUT "Terminal Taxa";
        print OUT_NOROUND "All percentages range 0-100% of the total count for either Eukaryotes or Prokaryotes, with exception of the first three data rows (% Known, % Eukaryote, % Prokaryote), which represent %s of total reads.\n";
        print OUT_NOROUND "Taxonomy String (K->S)\t";
        print OUT_NOROUND "Terminal Taxa";
        foreach my $i (0..$#unique_sampleHeaders)
            {   my $sample = $unique_sampleHeaders[$i];
                chomp($sample);
                print OUT "\t$sample";
                print OUT_NOROUND "\t$sample";
            }
        print OUT "\n";
        print OUT_NOROUND "\n";
        
        print OUT "\t% Known";
        print OUT_NOROUND "\t% Known";
        foreach my $i (0..$#unique_sampleHeaders)
            {   my $sample = $unique_sampleHeaders[$i];
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
        foreach my $i (0..$#unique_sampleHeaders)
            {   my $sample = $unique_sampleHeaders[$i];
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
        foreach my $i (0..$#unique_sampleHeaders)
            {   my $sample = $unique_sampleHeaders[$i];
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
                                        foreach my $j (0..$#unique_sampleHeaders)
                                            {   my $sample = $unique_sampleHeaders[$j];
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
                                        foreach my $j (0..$#unique_sampleHeaders)
                                            {   my $sample = $unique_sampleHeaders[$j];
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
                                foreach my $j (0..$#unique_sampleHeaders)
                                    {   my $sample = $unique_sampleHeaders[$j];
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
                                foreach my $j (0..$#unique_sampleHeaders)
                                    {   my $sample = $unique_sampleHeaders[$j];
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
                                foreach my $j (0..$#unique_sampleHeaders)
                                    {   my $sample = $unique_sampleHeaders[$j];
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
                                foreach my $j (0..$#unique_sampleHeaders)
                                    {   my $sample = $unique_sampleHeaders[$j];
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
                                foreach my $j (0..$#unique_sampleHeaders)
                                    {   my $sample = $unique_sampleHeaders[$j];
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
    }


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - S U B R O U T I N E S - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub ROUND
{	my $num = shift(@_);
	my $decimals = shift(@_);
	my $roundnum = int(($num * 10**$decimals) + 0.5)/(10**$decimals);
	return $roundnum;
}
# - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -
