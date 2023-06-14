#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use List::MoreUtils qw(uniq);

# - - - - - H E A D E R - - - - - - - - - - - - - - - - -
#Goals of script:
#Take x number marker files for ASVs_counts_NOUNKNOWNS_collapsedOnTaxonomy_percentabund.txt, x_asvTaxonomyTable_NOUNKNOWNS.txt, sample_metadata_forR.txt
#Match taxonomy and counts where overlap. Simplify to one file each. Prep for phyloseq ordination/network R file.

# - - - - - C O M M A N D    L I N E    O P T I O N S - - - - - - - -
my %options=();
getopts("i:o:h", \%options);

if ($options{h})
    {   print "\n\nHelp called:\nOptions:\n";
        print "-i = Input folder with internal folders: ASV_relabund, ASV_taxonomy, and Sample_metadata\n";
        print "     Each folder should have their corresponding files renamed to Marker.txt (i.e. COI.txt)\n";
        print "-o = Output directory\n";
        print "-h = This help message\n\n";
        die;
    }

# - - - - - G L O B A L  V A R I A B L E S  - - - - - -
my %ASVs;
my %TAXA;
my %META;
my @markers;
my @samples;

my $taxonomyfileheader;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - M A I N - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
opendir(DIR, $options{i}."/ASV_taxonomy") or die "\n\nInput directory $options{i}/ASV_taxonomy not found\n\n";
while (my $file = readdir(DIR))
    {   if ($file =~ m/\.txt$/)
            {   my $markername = $file;
                $markername =~ s/\.txt$//;
                push(@markers, $markername);
                open(IN2, "<$options{i}/ASV_relabund/$file") or die "\n\nThere is no $file, which means the marker names in your input folders don't match!!\n\n";
                my @data2 = <IN2>; close(IN2);
                my $samp_header = shift(@data2); chomp($samp_header);
                my @split_samp = split('\t', $samp_header);
                shift(@split_samp);
                foreach my $i (@split_samp)
                    {   push(@samples, $markername."_".$i)}
                foreach my $line (@data2)
                    {   chomp($line);
                        my @linearray = split('\t', $line);
                        my $asv = shift(@linearray);
                        foreach my $i (0..$#linearray)
                            {   $ASVs{$markername."_".$asv}{$markername."_".$split_samp[$i]} = $linearray[$i];
                            }
                    }
                
                open(IN, "<$options{i}/ASV_taxonomy/$file") or die "\n\nThere is no $file!!\n\n";
                my @data = <IN>; close(IN);
                my $header = shift(@data); chomp($header);
                $taxonomyfileheader = $header;
                foreach my $line (@data)
                    {   chomp($line);
                        my @linearray = split('\t', $line);
                        my $asv = shift(@linearray);
                        my $taxonomy = join(';', @linearray);
                        foreach my $i (sort keys %ASVs)
                            {   if ($markername."_".$asv eq $i)
                                    {   $ASVs{$markername."_".$asv}{"taxonomy"} = $taxonomy;
                                    }
                            }
                    }
                
                open(MET, "<$options{i}/Sample_metadata/$file") or die "\n\nThere is no $file, which means the marker names in your input folders don't match!!\n\n";
                my @metadata = <MET>; close(MET);
                my $meta_header = shift(@metadata); chomp($meta_header);
                my @split_metaheader = split('\t', $meta_header);
                my $sample_column_num;
                foreach my $i (0..$#split_metaheader)
                    {   if ($split_metaheader[$i] eq "Sample")
                            {   $sample_column_num = $i;
                            }
                    }
                foreach my $line (@metadata)
                    {   chomp($line);
                        my @linearray = split('\t', $line);
                        my $sample_name = $markername."_".$linearray[$sample_column_num];
                        if ($linearray[$sample_column_num] !~ m/^MP_/)
                            {die "\n\nPlease provide sample metadata files that were prepared for R. They should have 'MP_' at the front of the sample names.\n\n";}
                        foreach my $i (0..$#linearray)
                            {   unless ($i == $sample_column_num)
                                    {   $META{$sample_name}{$split_metaheader[$i]} = $linearray[$i];
                                    }
                            }
                    }
            }
    }
    
foreach my $i (sort keys %ASVs)
    {   unless (exists $ASVs{$i}{"taxonomy"})
            {   die "\n\nNo match to taxonomy for $i\n\n";
            }
    }

#Merge identical taxonomy across markers (using %TAXA)
foreach my $i (sort keys %ASVs)
    {   my $taxonomy = $ASVs{$i}{"taxonomy"};
        foreach my $j (@samples)
            {   if (exists $ASVs{$i}{$j})
                    {   $TAXA{$taxonomy}{'asvs'} .= $i.";";
                        $TAXA{$taxonomy}{$j} += $ASVs{$i}{$j};
                    }
                else
                    {   $TAXA{$taxonomy}{$j} += 0;
                    }
            }
    }
    
#Add new ASV numbering
my $num = 1;
foreach my $i (sort keys %TAXA)
    {   $TAXA{$i}{'new_asvID'} = "ASV_".$num;
        $num += 1;
    }
    
#Remove redundant ASVs from original lists
foreach my $i (sort keys %TAXA)
    {   my $asvs = $TAXA{$i}{'asvs'};
        my @asv_array = split(';', $asvs);
        my @uniq_asv = uniq @asv_array;
        my $new_asvs = join(';', @uniq_asv);
        $TAXA{$i}{'asvs'} = $new_asvs;
    }

#Print out new taxonomy table, ASV_relabund table, and new ASV name 2 old ASV names mapping file
open(TAXA, ">".$options{o}."/MergedMarkers_asvTaxonomyTable_NOUNKNOWNS.txt");
open(NAMES, ">".$options{o}."/MergedMarkers_newASVnames_mappingFile.txt");
open(RELABUND, ">".$options{o}."/MergedMarkers_ASVs_counts_NOUNKNOWNS_collapsedOnTaxonomy_percentabund.txt");
print TAXA "$taxonomyfileheader\n";
print NAMES "NEW_ASV_Name\tList_OriginalASVs_w_MatchingTaxonomy\n";
print RELABUND "x";
foreach my $i (0..$#samples)
    {   print RELABUND "\t$samples[$i]";
    }
print RELABUND "\n";
foreach my $i (sort keys %TAXA)
    {   print NAMES "$TAXA{$i}{'new_asvID'}\t$TAXA{$i}{'asvs'}\n";
        print TAXA "$TAXA{$i}{'new_asvID'}";
        my @taxonomyArray = split(';', $i);
        foreach my $j (@taxonomyArray)
            {   print TAXA "\t$j";
            }
        print TAXA "\n";
        print RELABUND "$TAXA{$i}{'new_asvID'}";
        foreach my $j (0..$#samples)
            {   print RELABUND "\t$TAXA{$i}{$samples[$j]}";
            }
        print RELABUND "\n";
    }
close(TAXA);
close(NAMES);
close(RELABUND);

#Create new sample metadata file with: Sample name, Marker, group1-X, and sites/replicates allowed.
#Merge and rename as needed.
my $groupTF = "FALSE";
my $groupNum;
my @groupNames;
my $sitesTF = "FALSE";
my $replicatesTF = "FALSE";
foreach my $i (sort keys %META)
    {   if (exists $META{$i}{"sites"})
            {   $sitesTF = "TRUE";}
        if (exists $META{$i}{"replicates"})
            {   $replicatesTF = "TRUE";}
        if (exists $META{$i}{"group1"})
            {   $groupTF = "TRUE";
                my $lastnum;
                foreach my $j (1..1000000000)
                    {   my $test = "group".$j;
                        unless (exists $META{$i}{$test})
                            {   $lastnum = $j - 1;
                                last;
                            }
                    }
                $groupNum = $lastnum;
            }
    }

open(METOUT, ">".$options{o}."/MergedMarkers_sample_metadata_forR.txt");
print METOUT "Sample\tMarker";
if ($sitesTF eq "TRUE")
    {print METOUT "\tsites";}
if ($replicatesTF eq "TRUE")
    {print METOUT "\treplicates";}
if ($groupTF eq "TRUE")
    {   foreach my $i (1..$groupNum)
            {   print METOUT "\tgroup".$i;
            }
    }
print METOUT "\n";
foreach my $i (sort keys %META)
    {   my $samp_name = $i;
        my $markername = $samp_name;
        $markername =~ s/^(.+)_MP_.+$/$1/;
        my $base_samp = $samp_name;
        $base_samp =~ s/^.+_(MP_.+)$/$1/;
        print METOUT "$i\t$markername";
        if ($sitesTF eq "TRUE")
            {   if (exists $META{$i}{"sites"})
                    {   print METOUT "\t".$META{$i}{"sites"};
                    }
                else
                    {   print METOUT "\t".$base_samp;
                    }
            }
        if ($replicatesTF eq "TRUE")
            {   if (exists $META{$i}{"replicates"})
                    {   print METOUT "\t".$META{$i}{"replicates"};
                    }
                else
                    {   print METOUT "\t".$base_samp;
                    }
            }
        if ($groupTF eq "TRUE")
            {   foreach my $j (1..$groupNum)
                    {   my $groupname = "group".$j;
                        if (exists $META{$i}{$groupname})
                            {   print METOUT "\t".$META{$i}{$groupname};
                            }
                        else
                            {   print METOUT "\tNA";
                            }
                    }
            }
        print METOUT "\n";
    }
close(METOUT);

#Create file showing only the ASVs shared by more than one marker
open(MULTI, ">".$options{o}."/MergedMarkers_listASVsSharedBy_multipleMarkers.txt");
print MULTI "New_ASV_Name\tMarkers\n";
foreach my $i (sort keys %TAXA)
    {   my @markercount;
        foreach my $j (@samples)
            {   if ($TAXA{$i}{$j} > 0)
                    {   my $samplename = $j;
                        $samplename =~ s/^(.+)_MP_.+$/$1/;
                        push(@markercount, $samplename);
                    }
            }
        my @uniq_marker = uniq @markercount;
        if (scalar(@uniq_marker) > 1)
            {   my $markerstring = join(';', @uniq_marker);
                print MULTI "$TAXA{$i}{'new_asvID'}\t$markerstring\n";
            }
    }
close(MULTI);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - S U B R O U T I N E S - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -
