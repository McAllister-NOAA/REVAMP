#!/usr/bin/perl -w
#use strict;

use Getopt::Std;
use List::MoreUtils qw(uniq);

# - - - - - H E A D E R - - - - - - - - - - - - - - - - -
#Goals of script:
#Take ASV count table, x_asvTaxonomyTable.txt from BLAST and SILVAngs independent runs (on same ASVs),
#Required settings: output name.
#Optional inputs: list of samples in the order desired for output.
#If you want to remove particular ASVs, that should have been done in the output for each blast and silva run.
#Outputs: ASV2Taxonomy folder contents.
#
#Instructions:
#   1) Modify progress.txt in the main out directory to delete everything including "taxonomyscriptFinished=TRUE" and below.
#   2) Enter the ASV2Taxonomy folder from one of your runs. Move contents to "ORIGINAL_RUN".
#   3) Run this script from the ASV2Taxonomy folder.
#   4) Confirm that "taxonomyscriptFinished=TRUE" in the progress.txt file.
#   5) Rerun metapipe to complete tables/figures.


# - - - - - C O M M A N D    L I N E    O P T I O N S - - - - - - - -
my %options=();
getopts("a:b:s:n:o:m:h", \%options);

if ($options{h})
    {   print "\n\nHelp called:\nOptions:\n";
        print "-a = ASV counts table (make sure there is text in the upperleft)\n";
        print "-b = ASV taxonomy table from BLAST run\n";
        print "-s = ASV taxonomy table from SILVAngs run\n";
        print "-n = Allin Output basename\n";
        print "-o = List of samples (one per line) in the order you want them exported. Must be exact matches to ASV counts table.\n";
        print "       Does not have to include all samples. (optional)\n";
        print "-m = Metapipe directory\n";
        print "-h = This help message\n\n";
        die;
    }

# - - - - - G L O B A L  V A R I A B L E S  - - - - - -
my %ASV;
my @sample_headers;
my @file_sample_headers;


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - M A I N - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#IMPORT ASV count table info
open(IN1, "<$options{a}") or die "\n\nThere is no $options{a} file!!\n\n";
my @asv_dat = <IN1>; close(IN1);
my $sample_header = shift(@asv_dat);
@file_sample_headers = split('\t', $sample_header);
my @new_array_headers;
foreach my $userheaders (@file_sample_headers)
	{	my $value = $userheaders;
		chomp($userheaders);
		push(@new_array_headers, $userheaders);	
	}
@file_sample_headers = @new_array_headers;
foreach my $line (@asv_dat)
	{	chomp($line);
		my @data = split('\t', $line);
		my $asv = $data[0]; chomp($asv);
		foreach my $i (1..$#data)
			{	my $number = $data[$i];
				chomp($number);
				$ASV{$asv}{$file_sample_headers[$i]} = $number; #Stores ASV count by sample
			}
	}

#IMPORT SAMPLE ORDER (if applicable)
if ($options{o})
	{	open(SAMP, "<$options{o}") or die "\n\nThere is no $options{o} file!!\n\n";
		my @usersampleorder = <SAMP>; close(SAMP);
		foreach my $i (@usersampleorder)
			{	chomp($i);
				push(@sample_headers, $i);
			}
	}
else
	{	shift(@file_sample_headers);
		@sample_headers = @file_sample_headers;
	}

foreach my $i (sort keys %ASV)
	{	foreach my $j (0..$#sample_headers)
			{	if (exists $ASV{$i}{$sample_headers[$j]})
					{}
				else {die "\n\nImported sample headers don't match the ASV count table. Try again.\n\n";}
			}
	}

#Import SILVA ASV Taxonomy file
open(SILVAIN, "<$options{s}") or die "\n\nThere is no $options{s} file!!\n\n";
my @silvadata = <SILVAIN>; close(SILVAIN); shift(@silvadata);
foreach my $line (@silvadata)
    {   chomp($line);
        my @split_line = split('\t', $line);
        my $asv = shift(@split_line);
        chomp($asv);
        foreach my $i (0..6)
            {   if ($split_line[$#split_line] eq "NA")
                    {   pop(@split_line);
                    }
                else
                    {   last;
                    }
            }
        my $taxonomyString = join(';', @split_line);
        if (exists $ASV{$asv})
            {   $ASV{$asv}{'silvaChoice'} = $taxonomyString;   
            }
    }

#Import BLAST ASV Taxonomy file
open(BLASTIN, "<$options{b}") or die "\n\nThere is no $options{b} file!!\n\n";
my @blastdata = <BLASTIN>; close(BLASTIN); shift(@blastdata);
foreach my $line (@blastdata)
    {   chomp($line);
        my @split_line = split('\t', $line);
        my $asv = shift(@split_line);
        chomp($asv);
        foreach my $i (0..6)
            {   if ($split_line[$#split_line] eq "NA")
                    {   pop(@split_line);
                    }
                else
                    {   last;
                    }
            }
        my $taxonomyString = join(';', @split_line);
        if (exists $ASV{$asv})
            {   $ASV{$asv}{'blastChoice'} = $taxonomyString;   
            }
    }

#Make decisions for finale merged BLAST/SILVA taxonomy assignments
open(BLSIMERGE, ">merged_BLAST_SILVA_taxonomyChoices.txt");
print BLSIMERGE "ASV\tBLASTtaxa\tSILVAtaxa\tChoice\n";
foreach my $i (sort keys %ASV)
    {   if (exists $ASV{$i}{'silvaChoice'} && exists $ASV{$i}{'blastChoice'})
            {   my $silvataxa = $ASV{$i}{'silvaChoice'};
                my @silvatax_array = split(';', $silvataxa);
                my $blasttaxa = $ASV{$i}{'blastChoice'};
                my @blasttax_array = split(';', $blasttaxa);
                my $chosen = "FALSE";
                if ($silvatax_array[0] eq "Bacteria" || $silvatax_array[0] eq "Archaea")        #SILVA Bacteria and Archaea
                    {   $ASV{$i}{'finaltaxachoice'} = $ASV{$i}{'silvaChoice'};
                        $chosen = "TRUE";
                    }
                if ($blasttax_array[0] eq "Eukaryota")                                          #BLAST Eukaryota
                    {   unless ($chosen eq "TRUE")
                            {   $ASV{$i}{'finaltaxachoice'} = $ASV{$i}{'blastChoice'};
                                $chosen = "TRUE";
                            }
                    }
                if ($silvatax_array[0] eq "Eukaryota" && $blasttax_array[0] ne "Eukaryota")     #SILVA Eukaryota
                    {   unless ($chosen eq "TRUE")
                            {   $ASV{$i}{'finaltaxachoice'} = $ASV{$i}{'silvaChoice'};
                                $chosen = "TRUE";
                            }
                    }
                if ($blasttax_array[0] eq "Archaea" || $blasttax_array[0] eq "Bacteria")
                    {   if ($silvatax_array[0] ne "Bacteria" || $silvatax_array[0] ne "Archaea")
                            {   unless ($chosen eq "TRUE")
                                    {   $ASV{$i}{'finaltaxachoice'} = $ASV{$i}{'silvaChoice'};
                                        $chosen = "TRUE";
                                    }
                            }
                    }
                if ($blasttax_array[0] =~ m/Unknown/ && $silvatax_array[0] ne "Unknown")        #Unknowns for both datasets
                    {   unless ($chosen eq "TRUE")
                            {   $ASV{$i}{'finaltaxachoice'} = $ASV{$i}{'silvaChoice'};
                            }
                    }
                elsif ($silvatax_array[0] eq "Unknown" && $blasttax_array[0] !~ m/Unknown/)
                    {   unless ($chosen eq "TRUE")
                            {   $ASV{$i}{'finaltaxachoice'} = $ASV{$i}{'blastChoice'};
                            }
                    }
                elsif ($blasttax_array[0] =~ m/Unknown/ && $silvatax_array[0] eq "Unknown")
                    {   unless ($chosen eq "TRUE")
                            {   $ASV{$i}{'finaltaxachoice'} = $ASV{$i}{'blastChoice'};
                            }
                    }
                unless (exists $ASV{$i}{'finaltaxachoice'})
                    {   $ASV{$i}{'finaltaxachoice'} = "Unknown";
                    }
                print BLSIMERGE "$i\t$ASV{$i}{'blastChoice'}\t$ASV{$i}{'silvaChoice'}\t$ASV{$i}{'finaltaxachoice'}\n";
                
            }
        elsif (exists $ASV{$i}{'silvaChoice'})
            {   $ASV{$i}{'finaltaxachoice'} = $ASV{$i}{'silvaChoice'};
                print BLSIMERGE "$i\tNA\t$ASV{$i}{'finaltaxachoice'}\t$ASV{$i}{'finaltaxachoice'}\n";
            }
        elsif (exists $ASV{$i}{'blastChoice'})
            {   $ASV{$i}{'finaltaxachoice'} = $ASV{$i}{'blastChoice'};
                print BLSIMERGE "$i\t$ASV{$i}{'finaltaxachoice'}\tNA\t$ASV{$i}{'finaltaxachoice'}\n";
            }
        else
            {   $ASV{$i}{'finaltaxachoice'} = "DELETE"; 
            }
    }
close(BLSIMERGE);

foreach my $i (sort keys %ASV)
    {   if ($ASV{$i}{'finaltaxachoice'} eq "DELETE")
            {   delete $ASV{$i};
            }
    }

##Outputs
foreach my $i (0..$#sample_headers)
	{	my $sample = $sample_headers[$i];
		chomp($sample);
		open($sample, ">".$sample."_KRONA.txt");
		#print $sample "count\ttaxa\n";
	}

open(OUT, ">".$options{n}."_allin_KRONA.txt");
open(WHOLEKRONA, ">".$options{n}."_wholeKRONA.txt");
open(ASVTAX, ">".$options{n}."_asvTaxonomyTable.txt");
print ASVTAX "ASV\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n";
print OUT "Sample\tASV\tcount\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\t(common name)\n";
foreach my $i (sort keys %ASV)
	{   my $taxput = $ASV{$i}{'finaltaxachoice'};
		my @taxputs = split(';', $taxput);
        print ASVTAX "$i";
        foreach my $l (@taxputs) {
			print ASVTAX "\t$l";}
		if ($#taxputs < 6)
			{	my $actualtaxdepth = $#taxputs + 1;
				foreach my $printtab ($actualtaxdepth..6)
					{	print ASVTAX "\tNA";
					}
			}
        print ASVTAX "\n";
        foreach my $j (0..$#sample_headers)
			{	my $hello = $sample_headers[$j];
				chomp($hello);
				unless ($ASV{$i}{$sample_headers[$j]} == 0)
				{
				print OUT "$sample_headers[$j]\t";
				print OUT "$i\t";
				print OUT "$ASV{$i}{$sample_headers[$j]}";
				foreach my $l (@taxputs) {
				print OUT "\t$l";}
				if ($#taxputs < 6)
					{	my $actualtaxdepth = $#taxputs + 1;
						foreach my $printtab ($actualtaxdepth..6)
							{	print OUT "\t";
							}
					}
				if (exists $ASV{$i}{'common_name'})
					{	print OUT "\t(".$ASV{$i}{'common_name'}.")";
					}
				else {print OUT "\t";}
				print OUT "\n";
				print $hello "$ASV{$i}{$sample_headers[$j]}";
				print WHOLEKRONA "$ASV{$i}{$sample_headers[$j]}";
				foreach my $l (@taxputs) {print $hello "\t$l"; print WHOLEKRONA "\t$l";}
				if (exists $ASV{$i}{'common_name'})
					{	print $hello " (".$ASV{$i}{'common_name'}.")";
						print WHOLEKRONA " (".$ASV{$i}{'common_name'}.")";
					}
				print $hello "\n";
				print WHOLEKRONA "\n";
				}
			}
	}
close(OUT);
close(WHOLEKRONA);
close(ASVTAX);

foreach my $i (0..$#sample_headers)
	{	my $sample = $sample_headers[$i];
		chomp($sample);
		close($sample);
	}

#BARCHART OUT
my @bartaxa;
my @commonnamebartaxa;
foreach my $j (0..$#sample_headers)
	{	foreach my $i (sort keys %ASV)
			{	unless ($ASV{$i}{$sample_headers[$j]} == 0)
					{	my $pushchoice = $ASV{$i}{'finaltaxachoice'};
						chomp($pushchoice);
						if ($pushchoice =~ m/;/)
							{	my @splittingthis = split(';', $pushchoice);
								my $lastassignment = $splittingthis[$#splittingthis];
								chomp($lastassignment);
								if (exists $ASV{$i}{'common_name'})
									{	my $newlastassignment = $lastassignment." (".$ASV{$i}{'common_name'}.")";
										push(@commonnamebartaxa, $newlastassignment);
									}
								push(@bartaxa, $lastassignment);
							}
						else {push(@bartaxa, $pushchoice);}
					}
			}
	}
my @uniq_bartaxa = uniq @bartaxa;

open(TAXAOUT, ">".$options{n}."_unique_terminaltaxa.txt");
foreach my $i (@uniq_bartaxa)
	{	unless($i eq "Unknown" || $i eq "Environmental Unknown" || $i =~ m/__/){
		print TAXAOUT "$i\n";}
	}
close(TAXAOUT);
my $namme = $options{n};
system("taxonkit name2taxid ".$namme."_unique_terminaltaxa.txt > ".$namme."_unique_terminaltaxa_w_taxids.txt");

open(TID, "<".$options{n}."_unique_terminaltaxa_w_taxids.txt") or die "\n\nSomething wrong with taxonkit name2taxid output\n\n";
my @commonname_last_dat = <TID>; close(TID);
my %CommonName_Term;
foreach my $line (@commonname_last_dat)
	{	chomp($line);
		my @line_split = split('\t', $line);
		if (exists $line_split[1])
            {my $taxid_clean = $line_split[1];
            chomp($taxid_clean);
            $CommonName_Term{$taxid_clean}{'name'} = $line_split[0];
            }
	}

unless (scalar(@commonnamebartaxa) == 0)
	{	my @uniq_commonnamebartaxa = uniq @commonnamebartaxa;
		open(COMMONBAR, ">".$options{n}."_SPECIESlevel_commonNames_for_barchart.txt");
		foreach my $i (@uniq_commonnamebartaxa)
			{	print COMMONBAR "$i\n";
			}
		close(COMMONBAR);
	}

open(BARCHART, ">".$options{n}."_barchart.txt");
open(BARCHART_forR, ">".$options{n}."_barchart_forR.txt");
print BARCHART_forR "Value\tSample\tTerminalTaxa\n";
print BARCHART "Sample\t";
foreach my $i (@uniq_bartaxa)
	{	print BARCHART "$i\t";
	}
print BARCHART "\n";
foreach my $j (0..$#sample_headers)
	{	my $barsampleheader = $sample_headers[$j];
		chomp($barsampleheader);
		print BARCHART "$barsampleheader\t";
		my %BARCHART;
		foreach my $i (sort keys %ASV)
			{	unless ($ASV{$i}{$sample_headers[$j]} == 0)
					{	my $pushchoice = $ASV{$i}{'finaltaxachoice'};
						chomp($pushchoice);
						if ($pushchoice =~ m/;/)
							{	my @splittingthis = split(';', $pushchoice);
								my $lastassignment = $splittingthis[$#splittingthis];
								chomp($lastassignment);
								$BARCHART{$lastassignment} += $ASV{$i}{$sample_headers[$j]};
							}
						else {$BARCHART{$pushchoice} += $ASV{$i}{$sample_headers[$j]};}
					}
			}
		foreach my $uniq_taxa (@uniq_bartaxa)
			{	my $hit = 0;
				foreach my $k (sort keys %BARCHART)
					{	if ($uniq_taxa eq $k)
							{	$hit = 1;
								print BARCHART "$BARCHART{$k}\t";
                                print BARCHART_forR "$BARCHART{$k}\t$sample_headers[$j]\t$uniq_taxa\n";
							}
					}
				if ($hit == 0)
					{	print BARCHART "0\t";
					}
			}
		print BARCHART "\n";
	}
close(BARCHART);
close(BARCHART_forR);

#run KRONA plot
my @kronasamples;
foreach my $i (0..$#sample_headers)
	{	my $hebs = $sample_headers[$i];
		chomp($hebs);
		my $newhebs = $hebs."_KRONA.txt";
		push(@kronasamples, $newhebs);
	}
my $printkronasamples = join(' ', @kronasamples);
my $naming = $options{n};
system("ImportText.pl -o ".$naming."_master_krona.html $printkronasamples");
my $wholekronaout = $options{n}."_wholeKRONA.txt";
system("ImportText.pl -o ".$naming."_wholeKRONA.html $wholekronaout");

##Barchart without unknowns
open(NOUNKNOWN, ">".$options{n}."_NO_UNKNOWNS_barchart.txt");
print NOUNKNOWN "Sample\t";
foreach my $i (@uniq_bartaxa)
	{	unless ($i eq "Unknown" || $i eq "Environmental Unknown")
			{	print NOUNKNOWN "$i\t";
			}
	}
print NOUNKNOWN "\n";
foreach my $j (0..$#sample_headers)
	{	my $barsampleheader = $sample_headers[$j];
		chomp($barsampleheader);
		print NOUNKNOWN "$barsampleheader\t";
		my %BARCHART;
		foreach my $i (sort keys %ASV)
			{	unless ($ASV{$i}{$sample_headers[$j]} == 0)
					{	my $pushchoice = $ASV{$i}{'finaltaxachoice'};
						chomp($pushchoice);
						if ($pushchoice =~ m/;/)
							{	my @splittingthis = split(';', $pushchoice);
								my $lastassignment = $splittingthis[$#splittingthis];
								chomp($lastassignment);
								$BARCHART{$lastassignment} += $ASV{$i}{$sample_headers[$j]};
							}
						else {$BARCHART{$pushchoice} += $ASV{$i}{$sample_headers[$j]};}
					}
			}
		foreach my $uniq_taxa (@uniq_bartaxa)
			{	unless ($uniq_taxa eq "Unknown" || $uniq_taxa eq "Environmental Unknown")
					{	my $hit = 0;
						foreach my $k (sort keys %BARCHART)
							{	if ($uniq_taxa eq $k)
									{	$hit = 1;
										print NOUNKNOWN "$BARCHART{$k}\t";
									}
							}
						if ($hit == 0)
							{	print NOUNKNOWN "0\t";
							}
					}
			}
		print NOUNKNOWN "\n";
	}
close(NOUNKNOWN);

##Taxa with shared ASV w/ heatmap
my %TAXHEAT;

foreach my $i (sort keys %ASV)
	{	my $tax = $ASV{$i}{'finaltaxachoice'};
		chomp($tax);
		if ($tax =~ m/;/)
			{	my @split_tax = split(';', $tax);
				foreach my $j (0..$#split_tax)
					{	if (exists $TAXHEAT{$split_tax[$j]})
							{	$TAXHEAT{$split_tax[$j]}{'asvs'} .= ";".$i;
							}
						else
							{	$TAXHEAT{$split_tax[$j]}{'asvs'} = $i;
								$TAXHEAT{$split_tax[$j]}{'depth'} = $j + 1;
							}
					}
			}
		else
			{	if (exists $TAXHEAT{$tax})
					{	$TAXHEAT{$tax}{'asvs'} .= ";".$i;
					}
				else
					{	$TAXHEAT{$tax}{'asvs'} = $i;
						$TAXHEAT{$tax}{'depth'} = 1;
					}
			}	
	}

open(OUTDEEP1, ">".$options{n}."_heatmap_multiASV.txt");
foreach my $i (sort keys %TAXHEAT)
	{	if ($TAXHEAT{$i}{'asvs'} =~ m/;/)
			{	print OUTDEEP1 "$i <<Depth ".$TAXHEAT{$i}{'depth'}.">>\tSample\t";
				my $multiasv = $TAXHEAT{$i}{'asvs'};
				my @multiasv_split = split(';', $multiasv);
				foreach my $j (@multiasv_split)
					{	print OUTDEEP1 "$j\t";
					}
				print OUTDEEP1 "\n";
				foreach my $k (0..$#sample_headers)
					{	print OUTDEEP1 "$i <<Depth ".$TAXHEAT{$i}{'depth'}.">>\t";
						print OUTDEEP1 "$sample_headers[$k]\t";
						foreach my $j (@multiasv_split)
							{	print OUTDEEP1 "$ASV{$j}{$sample_headers[$k]}\t";
							}
						print OUTDEEP1 "\n";
					}
			}
	}
close(OUTDEEP1);

open(UNKNOWNS, ">".$options{n}."_unknown_asvids.txt");
foreach my $i (sort keys %ASV)
    {   my $taxonunknown = $ASV{$i}{'finaltaxachoice'};
        chomp($taxonunknown);
        if ($taxonunknown eq "Unknown" || $taxonunknown eq "Environmental Unknown")
            {   print UNKNOWNS "$i\t\n";
            }
    }
close(UNKNOWNS);

#Complete lines from shell to complete the ASV2Taxonomy folder
system('cat '.$options{n}.'_asvTaxonomyTable.txt | grep -v "Unknown" > '.$options{n}.'_asvTaxonomyTable_NOUNKNOWNS.txt');
system('cat ../dada2/ASVs_counts.tsv | grep -v -f '.$options{n}.'_unknown_asvids.txt > ASVs_counts_NOUNKNOWNS.tsv');
system('perl '.$options{m}.'/assets/merge_on_taxonomy.pl -a ../dada2/ASVs_counts.tsv -t '.$options{n}.'_asvTaxonomyTable.txt -o ./ > ASVs_counts_mergedOnTaxonomy.tsv');
system('cat ASVs_counts_mergedOnTaxonomy.tsv | grep -v -f '.$options{n}.'_unknown_asvids.txt > ASVs_counts_mergedOnTaxonomy_NOUNKNOWNS.tsv');
system('mkdir KRONA_plots');
system('mkdir KRONA_plots/KRONA_inputs');
system('mv MP* KRONA_plots/KRONA_inputs/');
system('mv '.$options{n}.'_wholeKRONA.txt KRONA_plots/KRONA_inputs/'.$options{n}.'_samplesSummedKRONA.txt');
system('mv '.$options{n}.'_master_krona.html KRONA_plots/');
system('mv '.$options{n}.'_wholeKRONA.html KRONA_plots/'.$options{n}.'_samplesSummedKRONA.html');
system('mv '.$options{n}.'_allin_KRONA.txt '.$options{n}.'_allin_TaxonomyASVSampleCount_byline.txt');
system('perl '.$options{m}.'/assets/stats.pl -a ../dada2/ASVs_counts.tsv -t ASVs_counts_mergedOnTaxonomy.tsv -i '.$options{n}.'_asvTaxonomyTable.txt > basic_ASV_taxonomy_stats.txt');
system('echo "taxonomyscriptFinished=TRUE" >> ../progress.txt');

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - S U B R O U T I N E S - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


# - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -
