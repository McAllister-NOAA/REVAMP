#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use List::Util qw( min max );
use List::MoreUtils qw(uniq);

# - - - - - H E A D E R - - - - - - - - - - - - - - - - -
#Goals of script:
#Compare # of shared/uniq taxa/terminal taxa at each level (P,C,O,F,G,S). Output the ASVs that are in each pie.
#Output counts for each marker alone as well.

# - - - - - C O M M A N D    L I N E    O P T I O N S - - - - - - - -
my %options=();
getopts("i:o:f:l:h", \%options);

if ($options{h})
    {   print "\n\nHelp called:\nOptions:\n";
        print "-i = Input directory with x_asvTaxonomyTable_NOUNKNOWNS.txt file renamed as Marker.txt\n";
        print "-o = Output directory\n";
        print "-f = File of taxa of interest (one per line)\n";
        print "-l = Taxonomic level (Kingdom/Phylum/Class/Order/Family/Genus/Species) of taxa of interest file\n";
        print "-h = This help message\n\n";
        die;
    }

# - - - - - G L O B A L  V A R I A B L E S  - - - - - -
my %TAXA;
my @markers;
my @taxaHierarchy = qw | Kingdom Phylum Class Order Family Genus Species |;
my @taxaOfInterest;
my $levelTaxaOfInterest;
my $levelTOI_loc;

system("mkdir -p ".$options{o});

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - M A I N - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if ($options{f})
    {   open(INTEREST, "<".$options{f}) or die "\n\nNo infile $options{f} given\n\n";
        my @dataint = <INTEREST>; close(INTEREST);
        foreach my $f (@dataint)
            {   chomp($f);
                $f =~ s/\ /\+/g;
                $f =~ s/\s//g;
                $f =~ s/\+/\ /g;
                push(@taxaOfInterest, $f);
            }
        $levelTaxaOfInterest = uc($options{l});
        foreach my $f (0..$#taxaHierarchy)
            {   if ($levelTaxaOfInterest eq uc($taxaHierarchy[$f]))
                    {   $levelTOI_loc = $f;
                    }
            }
    }


opendir(DIR, $options{i}) or die "\n\nInput directory $options{i} not found\n\n";
while (my $file = readdir(DIR))
    {   if ($file =~ m/\.txt$/)
            {   my $markername = $file;
                $markername =~ s/\.txt$//;
                push(@markers, $markername);
                open(IN, "<$options{i}/$file") or die "\n\nThere is no $file!!\n\n";
                my @data = <IN>; close(IN);
                my $header = shift(@data); chomp($header);
                foreach my $line (@data)
                    {   chomp($line);
                        my @linearray = split('\t', $line);
                        my $asv = shift(@linearray);
                        my $taxonomy = join(';', @linearray);
                        $TAXA{$taxonomy}{$markername}{"asvs"} .= $asv.";";
                    }
            }
    }
    
if ($options{f})
    {   foreach my $f (sort keys %TAXA)
            {   my @splitf = split(';', $f);
                my $passfilt = "FALSE";
                foreach my $j (@taxaOfInterest)
                    {   if ($j eq $splitf[$levelTOI_loc])
                            {   $passfilt = "TRUE";
                            }
                    }
                if ($passfilt eq "FALSE")
                    {   delete($TAXA{$f});
                    }
            }
    }
    
foreach my $i (sort keys %TAXA)
    {   foreach my $j (@markers)
            {   if (exists $TAXA{$i}{$j})
                {   my $asv = $TAXA{$i}{$j}{"asvs"};
                    my @asv_array = split(';', $asv);
                    my $new_asv = join(';', @asv_array);
                
                    $TAXA{$i}{$j}{"asvs"} = $new_asv;
                }
            }
    }

my @unique_markers = uniq @markers;
if (scalar @unique_markers != scalar @markers)
    {die "\n\nSomething not working with duplicate marker files\n\n";}
my $total_num_markers = scalar @markers;


#Print table showing #unique taxa at any level and the #terminal taxa at each level for each marker
open(MARKER, ">".$options{o}."/taxacounts_individualMarkersOnly.txt");
print MARKER "TaxaLevel_counts";
foreach my $i (@markers)
    {   print MARKER "\t$i";
    }
print MARKER "\n";

foreach my $i (0..$#taxaHierarchy)
    {   print MARKER "uniqueTotal_at_$taxaHierarchy[$i]";
        foreach my $j (@markers)
            {   my @taxonomyPresent;
                foreach my $k (sort keys %TAXA)
                    {   if (exists $TAXA{$k}{$j})
                            {   push(@taxonomyPresent, $k);
                            }
                    }
                my @chosenlevel;
                foreach my $k (@taxonomyPresent)
                    {   my @taxarray = split(';', $k);
                        unless ($taxarray[$i] eq "NA")
                            {   push(@chosenlevel, $taxarray[$i]);
                            }
                    }
                my @uniq_chosen = uniq @chosenlevel;
                my $printval = scalar @uniq_chosen;
                print MARKER "\t$printval";
            }
        print MARKER "\n";
    }
    
foreach my $i (0..$#taxaHierarchy)
    {   print MARKER "uniqueTerminal_at_$taxaHierarchy[$i]";
        foreach my $j (@markers)
            {   my @taxonomyPresent;
                foreach my $k (sort keys %TAXA)
                    {   if (exists $TAXA{$k}{$j})
                            {   push(@taxonomyPresent, $k);
                            }
                    }
                my @terminalTaxa;
                foreach my $k (@taxonomyPresent)
                    {   my @taxarray = split(';', $k);
                        my $lasthit;
                        foreach my $entry (0..$#taxarray)
                            {   if ($taxarray[$entry] ne "NA")
                                    {   $lasthit = $entry;
                                    }
                            }
                        my @newtaxarray;
                        foreach my $entry (0..$lasthit)
                            {   push(@newtaxarray, $taxarray[$entry])
                            }
                        my $goodtax = join(';', @newtaxarray);
                        push(@terminalTaxa, $goodtax);
                    }
                my @unique_terminalTaxa = uniq @terminalTaxa;
                
                my $count = 0;
                foreach my $k (@unique_terminalTaxa)
                    {   my @taxarray = split(';', $k);
                        if (scalar @taxarray == $i + 1)
                            {   $count += 1;
                            }
                    }
                print MARKER "\t$count";
            }
        print MARKER "\n";
    }
    
close(MARKER);


#Compare # of shared/uniq taxa at each level (P,C,O,F,G,S). Output the ASVs that are in each pie.

my @comparisonArray;
foreach my $entry (@markers)
    {   push(@comparisonArray, map [@$_, $entry], @comparisonArray, []);
    }
#print "@$_\n" for @comparisonArray;

my @newarray;
foreach my $i (@comparisonArray)
    {   my @comp = @{ $i };
        my $comp_list = join(';', @comp);
        push(@newarray, $comp_list);
    }

my @temp_array;
foreach my $i (1..$total_num_markers)
    {   foreach my $j (sort @newarray)
            {   my @splitting = split(';', $j);
                my $num = scalar @splitting;
                if ($num == $i)
                    {   push(@temp_array, $j);
                    }
            }
    }
@newarray = @temp_array;
    

open(COMP, ">".$options{o}."/taxacounts_comparisonsUniqueTotal.txt");
print COMP "Comparison_NumUniqueTotalSharedByXOnly";
foreach my $i (0..$#taxaHierarchy)
    {   print COMP "\t$taxaHierarchy[$i]";
    }
print COMP "\n";

open(TLIST, ">".$options{o}."/taxa_comparisonsUniqueTotal.txt");
print TLIST "Comparison_UniqueSharedByXOnly";
foreach my $i (0..$#taxaHierarchy)
    {   print TLIST "\t$taxaHierarchy[$i]";
    }
print TLIST "\n";

open(ASVlist, ">".$options{o}."/ASVs_comparisonsUniqueTotal.txt");
print ASVlist "Comparison_ASVsWithUniqueTaxaSharedByXOnly";
foreach my $i (0..$#taxaHierarchy)
    {   print ASVlist "\t$taxaHierarchy[$i]";
    }
print ASVlist "\n";

foreach my $i (@newarray)
    {   print COMP "$i";
        print TLIST "$i";
        print ASVlist "$i";
        my @markerArray = split(';', $i); #list of markers we are comparing
        my @anti; #list of markers that are not in the comparison (matches need to be exclusive)
        my $present;
        foreach my $j (0..$#markerArray)
            {   foreach my $k (0..$#markers)
                    {   if ($markerArray[$j] eq $markers[$k])
                            {   $present .= $k.";";
                            }
                    }
            }
        my @presentArr = split(';', $present);
        foreach my $k (0..$#markers)
            {   my $hit = "FALSE";
                foreach my $j (@presentArr)
                    {   if ($j == $k)
                            {   $hit = "TRUE";
                            }
                    }
                if ($hit eq "FALSE")
                    {   push(@anti, $markers[$k]);
                    }
            }
        
        #print "\n$i:\nMarkers:<<<@markerArray>>>\n";
        #print "Anti:<<<@anti>>>\n";
        
        foreach my $j (0..$#taxaHierarchy)
            {   my @sharedtaxa;
                foreach my $k (sort keys %TAXA)
                    {   my @hierarch = split(';', $k);
                        my $chosentax = $hierarch[$j];
                        if ($chosentax ne "NA")
                        {push(@sharedtaxa, $chosentax);}
                    }
                my @uniq_sharedtaxa = uniq @sharedtaxa;
                #print "$i\t<<<@uniq_sharedtaxa>>>\n";
                my $count = 0;
                my $listedtaxa;
                foreach my $share (@uniq_sharedtaxa)
                    {   my @hit_good;
                        my @hit_bad;
                        foreach my $b (@markerArray)
                            {   push(@hit_good, "FALSE");}
                        if (scalar @anti > 0)
                            {   foreach my $b (@anti)
                                    {   push(@hit_bad, "FALSE");
                                    }
                            }
                        foreach my $k (sort keys %TAXA)
                            {   my @hierarch = split(';', $k);
                                if ($hierarch[$j] eq $share)
                                    {   foreach my $h (0..$#markerArray)
                                            {   if (exists $TAXA{$k}{$markerArray[$h]})
                                                    {   $hit_good[$h] = "TRUE";
                                                    }
                                            }
                                        if (scalar @anti > 0)
                                        {   foreach my $h (0..$#anti)
                                                {   if (exists $TAXA{$k}{$anti[$h]})
                                                        {   $hit_bad[$h] = "TRUE";
                                                        }
                                                }
                                        }
                                    }
                            }
                        my $good_pass = "TRUE";
                        foreach my $h (@hit_good)
                            {   if ($h eq "FALSE")
                                    {   $good_pass = "FALSE";
                                    }
                            }
                        my $bad_fail = "FALSE";
                        if (scalar @anti > 0)
                            {   foreach my $h (@hit_bad)
                                    {   if ($h eq "TRUE")
                                            {   $bad_fail = "TRUE";
                                            }
                                    }
                            }
                        if ($good_pass eq "TRUE" && $bad_fail eq "FALSE")
                            {   $count += 1;
                                $listedtaxa .= $share.";";
                            }
                    }
                    
                my $asvs_keep;
                    if (defined $listedtaxa)
                        {   my @listedArray = split(';', $listedtaxa);
                            foreach my $bones (@listedArray)
                                {   foreach my $w (sort keys %TAXA)
                                        {   my @hierarch = split(';', $w);
                                            if ($hierarch[$j] eq $bones)
                                                {   foreach my $bp (@markerArray)
                                                        {   if (exists $TAXA{$w}{$bp})
                                                                {$asvs_keep .= $bp."<".$TAXA{$w}{$bp}{"asvs"}.">;";}
                                                        }
                                                }
                                        }
                                }
                        }
                    
                print COMP "\t$count";
                if (defined $listedtaxa)
                    {   print TLIST "\t$listedtaxa";
                    }
                else
                    {   print TLIST "\t";
                    }
                if (defined $asvs_keep)
                    {   print ASVlist "\t$asvs_keep";
                    }
                else
                    {   print ASVlist "\t";
                    }
            }
        print COMP "\n";
        print TLIST "\n";
        print ASVlist "\n";
    }

close(COMP);
close(TLIST);
close(ASVlist);

#Version of COMP with rows with all 0s removed
open(INCOMP, "<".$options{o}."/taxacounts_comparisonsUniqueTotal.txt");
my @data = <INCOMP>; close(INCOMP);
my $header = shift(@data); chomp($header);

my @keepdata;
foreach my $line (@data)
    {   chomp($line);
        my @split_line = split('\t', $line);
        my $pass = "FALSE";
        foreach my $i (1..$#split_line)
            {   if ($split_line[$i] != 0)
                    {   $pass = "TRUE";
                    }
            }
        if ($pass eq "TRUE")
            {   push(@keepdata, $line);
            }
    }
unshift(@keepdata, $header);
unshift(@data, $header);

my $comparedata = join('', @data);
my $comparekeep = join('', @keepdata);
if ($comparedata ne $comparekeep)
    {   open(OUTCOMP2, ">".$options{o}."/taxacounts_comparisonsUniqueTotal_allZeroRowsRemoved.txt");
        foreach my $i (@keepdata)
            {   chomp($i);
                print OUTCOMP2 "$i\n";
            }
        close(OUTCOMP2);
    }

#Compare the # of shared terminal taxa between markers.

open(TERMTAX1, ">".$options{o}."/taxacounts_comparisonsUniqueTerminal.txt");
print TERMTAX1 "Comparison_NumUniqueTerminalSharedByXOnly";
foreach my $i (0..$#taxaHierarchy)
    {   print TERMTAX1 "\t$taxaHierarchy[$i]";
    }
print TERMTAX1 "\n";

open(TERMTAX2, ">".$options{o}."/taxa_comparisonsUniqueTerminal.txt");
print TERMTAX2 "Comparison_UniqueTerminalSharedByXOnly";
foreach my $i (0..$#taxaHierarchy)
    {   print TERMTAX2 "\t$taxaHierarchy[$i]";
    }
print TERMTAX2 "\n";

open(TERMTAX3, ">".$options{o}."/ASVs_comparisonsUniqueTerminal.txt");
print TERMTAX3 "Comparison_ASVsWithinUniqueTerminalSharedByXOnly";
foreach my $i (0..$#taxaHierarchy)
    {   print TERMTAX3 "\t$taxaHierarchy[$i]";
    }
print TERMTAX3 "\n";

foreach my $i (@newarray)
    {   print TERMTAX1 "$i";
        print TERMTAX2 "$i";
        print TERMTAX3 "$i";
        my @markerArray = split(';', $i); #list of markers we are comparing
        my @anti; #list of markers that are not in the comparison (matches need to be exclusive)
        my $present;
        foreach my $j (0..$#markerArray)
            {   foreach my $k (0..$#markers)
                    {   if ($markerArray[$j] eq $markers[$k])
                            {   $present .= $k.";";
                            }
                    }
            }
        my @presentArr = split(';', $present);
        foreach my $k (0..$#markers)
            {   my $hit = "FALSE";
                foreach my $j (@presentArr)
                    {   if ($j == $k)
                            {   $hit = "TRUE";
                            }
                    }
                if ($hit eq "FALSE")
                    {   push(@anti, $markers[$k]);
                    }
            }
        foreach my $j (0..$#taxaHierarchy)
            {   my @terminalTaxa;
                foreach my $k (@markerArray)
                    {   my @taxonomyPresent;
                        foreach my $l (sort keys %TAXA)
                            {   if (exists $TAXA{$l}{$k})
                                    {   push(@taxonomyPresent, $l);
                                    }
                            }
                        foreach my $l (@taxonomyPresent)
                            {   my @taxarray = split(';', $l);
                                my $lasthit;
                                foreach my $entry (0..$#taxarray)
                                    {   if ($taxarray[$entry] ne "NA")
                                            {   $lasthit = $entry;
                                            }
                                    }
                                my @newtaxarray;
                                foreach my $entry (0..$lasthit)
                                    {   push(@newtaxarray, $taxarray[$entry])
                                    }
                                my $goodtax = join(';', @newtaxarray);
                                push(@terminalTaxa, $goodtax);
                            }
                    }
                my @unique_terminalTaxa = uniq @terminalTaxa;
                
                #Make sure all the terminal taxa are shared by all markers in @markerArray
                my @really_keep;
                foreach my $share (@unique_terminalTaxa)
                    {   my @hit_good;
                        foreach my $b (@markerArray)
                                    {   push(@hit_good, "FALSE");
                                    }
                        foreach my $key (sort keys %TAXA)
                                    {   my @hierarch = split(';', $key);
                                        my $lasthit;
                                        foreach my $entry (0..$#hierarch)
                                            {   if ($hierarch[$entry] ne "NA")
                                                    {   $lasthit = $entry;
                                                    }
                                            }
                                        my @newtaxarray;
                                        foreach my $entry (0..$lasthit)
                                            {   push(@newtaxarray, $hierarch[$entry])
                                            }
                                        my $goodtermtax = join(';', @newtaxarray);
                                        if ($goodtermtax eq $share)
                                            {   foreach my $a (0..$#markerArray)
                                                    {   if (exists $TAXA{$key}{$markerArray[$a]})
                                                            {   $hit_good[$a] = "TRUE";
                                                            }
                                                    }
                                            }
                                    }
                        my $good_pass = "TRUE";
                        foreach my $h (@hit_good)
                            {   if ($h eq "FALSE")
                                    {   $good_pass = "FALSE";
                                    }
                            }
                        if ($good_pass eq "TRUE")
                            {   push(@really_keep, $share);
                            }
                    }
                
                @unique_terminalTaxa = @really_keep;
                
                my @shared_terminalTaxa_atHierarch;    
                foreach my $k (@unique_terminalTaxa)
                     {   my @taxarray = split(';', $k);
                         if (scalar @taxarray == $j + 1)
                             {   push(@shared_terminalTaxa_atHierarch, $k);
                             }
                     }
                
                my @shared_terminalTaxa_notInAnti; #list of shared terminal taxa at the hierarchy tested
                if (scalar @anti > 0)
                    {   foreach my $share (@shared_terminalTaxa_atHierarch)
                            {   my @hit_bad;
                                foreach my $b (@anti)
                                    {   push(@hit_bad, "FALSE");
                                    }
                                foreach my $key (sort keys %TAXA)
                                    {   my @hierarch = split(';', $key);
                                        my $lasthit;
                                        foreach my $entry (0..$#hierarch)
                                            {   if ($hierarch[$entry] ne "NA")
                                                    {   $lasthit = $entry;
                                                    }
                                            }
                                        my @newtaxarray;
                                        foreach my $entry (0..$lasthit)
                                            {   push(@newtaxarray, $hierarch[$entry])
                                            }
                                        my $goodtermtax = join(';', @newtaxarray);
                                        if ($goodtermtax eq $share)
                                            {   foreach my $a (0..$#anti)
                                                    {   if (exists $TAXA{$key}{$anti[$a]})
                                                            {   $hit_bad[$a] = "TRUE";
                                                            }
                                                    }
                                            }
                                    }
                                my $bad_fail = "FALSE";
                                foreach my $h (@hit_bad)
                                    {   if ($h eq "TRUE")
                                            {   $bad_fail = "TRUE";
                                            }
                                    }
                                if ($bad_fail eq "FALSE")
                                    {   push(@shared_terminalTaxa_notInAnti, $share);
                                    }
                            }
                    }
                else
                    {   @shared_terminalTaxa_notInAnti = @shared_terminalTaxa_atHierarch;
                    }
                
                my $count = scalar @shared_terminalTaxa_notInAnti;
                print TERMTAX1 "\t$count";
                print TERMTAX2 "\t";
                if (scalar @shared_terminalTaxa_notInAnti > 0)
                    {   foreach my $p (@shared_terminalTaxa_notInAnti)
                            {   my @splittx = split(';', $p);
                                print TERMTAX2 "$splittx[$#splittx];";
                            }
                    }
                
                my $asvs_keep;
                foreach my $p (@shared_terminalTaxa_notInAnti)
                    {   my @splitpea = split(';', $p);
                        while (scalar @splitpea < 7)
                            {   push(@splitpea, "NA");
                            }
                        my $newcheck = join(';', @splitpea);
                        foreach my $mark (@markers)
                            {   if (exists $TAXA{$newcheck}{$mark})
                                    {   $asvs_keep .= $mark."<".$TAXA{$newcheck}{$mark}{"asvs"}.">;";
                                    }
                            }
                    }
                if (defined $asvs_keep)
                    {print TERMTAX3 "\t$asvs_keep";}
                else {print TERMTAX3 "\t";}
            }
        print TERMTAX1 "\n";
        print TERMTAX2 "\n";
        print TERMTAX3 "\n";
    }
close(TERMTAX1);
close(TERMTAX2);
close(TERMTAX3);

#Version of TERMTAX1 with rows with all 0s removed
open(INTERM, "<".$options{o}."/taxacounts_comparisonsUniqueTerminal.txt");
my @datas = <INTERM>; close(INTERM);
my $headerTERM = shift(@datas); chomp($headerTERM);

my @keepTERMdata;
foreach my $line (@datas)
    {   chomp($line);
        my @split_line = split('\t', $line);
        my $pass = "FALSE";
        foreach my $i (1..$#split_line)
            {   if ($split_line[$i] != 0)
                    {   $pass = "TRUE";
                    }
            }
        if ($pass eq "TRUE")
            {   push(@keepTERMdata, $line);
            }
    }
unshift(@keepTERMdata, $headerTERM);
unshift(@datas, $headerTERM);

my $compareTERMdata = join(' ', @datas);
my $compareTERMkeep = join(' ', @keepTERMdata);
if ($compareTERMdata ne $compareTERMkeep)
    {   open(OUTTERM2, ">".$options{o}."/taxacounts_comparisonsUniqueTerminal_allZeroRowsRemoved.txt");
        foreach my $i (@keepTERMdata)
            {   chomp($i);
                print OUTTERM2 "$i\n";
            }
        close(OUTTERM2);
    }


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - S U B R O U T I N E S - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -
