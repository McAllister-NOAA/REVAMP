#!/usr/bin/perl -w
#use strict;
use Getopt::Std;
use List::MoreUtils qw(uniq);

$|=1; #Autoflush so pauses within this script are pushed up to the shell script during run.

# - - - - - H E A D E R - - - - - - - - - - - - - - - - -
#Goals of script:
#Take ASV count table, taxAssignX.txt, reformatted_taxonkit_out, genbank common names from NCBI.
#Required settings: filtering options for species, genus, etc cutoffs, output name.
#Optional inputs: list of ASVs to ignore, list of samples in the order desired for output.
#Outputs: Reformatted out with one ASV per sample per line, Krona plots, Bar charts (with and without unknowns),
#         multi-ASV taxa heatmap, common names outputs, list of unknown ASVs. All with and without filtering.


# - - - - - C O M M A N D    L I N E    O P T I O N S - - - - - - - -
my %options=();
getopts("a:s:t:m:n:f:c:d:o:r:y:z:eh", \%options);

if ($options{h})
    {   print "\n\nHelp called:\nOptions:\n";
        print "-a = ASV counts table (make sure there is text in the upperleft)\n";
        print "-s = reformatted blast output (ASV_blastn_nt_formatted.txt; w/ headers ASV percent length taxid accession correction)\n";
        print "-t = reformatted taxonkit output (TaxID\tsimplified taxonomy)\n";
        print "-f = filtering options Species,Genus,Family,Order,Class,Phylum (e.g. 97,95,90,80,70,60)\n";
        print "-n = Allin Output basename\n";
        print '-c = Location of common names file (grep "genbank common name" from names.dmp NCBI taxonomy file)'; print "\n";
        print "-d = List of ASVs to ignore (one per line) for outputs ignoring contaminants and/or unknowns\n";
        print "-o = List of samples (one per line) in the order you want them exported. Must be exact matches to ASV counts table.\n";
        print "       Does not have to include all samples.\n";
        print "-e = Toggle use of SILVAngs taxonomy assignments by ASV (optional)\n";
        print "-y = SILVAngs results/[ls]su/exports/*---otus.csv File\n";
        print "-z = SILVAngs results/[ls]su/stats/sequence_cluster_map/data/*.fa.clstr File\n";
        print "-r = Reference taxonomy map for current SILVA database: i.e. tax_slv_ssu_138.1.txt\n";
        print "-m = Location of the REVAMP directory\n";
        print "-h = This help message\n\n";
        die;
    }

# - - - - - G L O B A L  V A R I A B L E S  - - - - - -
my %ASV;
my %TAXON;
my @uniq_bartaxa_ig;
my @sample_headers;
my @file_sample_headers;
my @commondat;
my %SILVAREF;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - M A I N - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if ($options{e})
    {   open(REFIN, "<$options{r}") or die "\n\nThere is no $options{r} file!!\n\n";
        my @refdat = <REFIN>; close(REFIN);
        foreach my $line (@refdat)
            {   chomp($line);
                my @split_line = split('\t', $line);
                my $header = $split_line[0];
                $header =~ s/[^A-Za-z0-9; ]/_/g;
                my $silvataxid = $split_line[1];
                my $level = $split_line[2];
                $SILVAREF{$header}{'silvataxid'} = $silvataxid;
                $SILVAREF{$header}{'level'} = $level;
            }
    }
    
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

#IMPORT List of ASVs to ingore in output
my @ignore_array;
my $output_ignore_string;
my $ignore_ASV_string;
if ($options{d})
    {	open(IGNORE, "<$options{d}") or die "\n\nThere is no $options{d} file!!\n\n";
        my @igdat = <IGNORE>; close(IGNORE);
        foreach my $linei (@igdat)
            {	chomp($linei);
                push(@ignore_array, $linei);
            }
        foreach my $i (@ignore_array)
            {	$ignore_ASV_string .= "_".$i."_";
                $output_ignore_string .= $i.";";
            }
        my @ignoreASVarray = split(';', $output_ignore_string);
        open(ASVsToIG, ">".$options{n}."_ASVs_to_IGNORE.txt");
        foreach my $igasv (@ignoreASVarray)
            {   chomp($igasv);
                print ASVsToIG "$igasv\n";
            }
        close(ASVsToIG);
    }

if ($options{e})
    {   #Import SILVAngs export file
        #Expected format: [2] = cluster ACC (ASV_#); [3] = # sequences; [8] = ncbi taxonomic classification; [9] = silva taxonomic classification
        my %SILVATAX;
        
        open(SILVAIN1, "<$options{y}") or die "\n\nThere is no $options{y} file!!\n\n";
        my @silvadata1 = <SILVAIN1>; close(SILVAIN1); shift(@silvadata1);
        my $heads = shift(@silvadata1);
        if ($heads !~ m/^sample name/)
            {   die "\n\nSILVA export header location has changed and will break this program. Ask for modification.\n\n";
            }
        foreach my $line (@silvadata1)
            {	chomp($line);
                my @split_line = split('\t', $line);
                my $asv = $split_line[2];
                chomp($asv);
                my $ncbi_tax;
                my $ncbi_taxid;
                if ($split_line[8] ne "")
                    {   $ncbi_tax = $split_line[8];
                        $ncbi_tax =~ m/^ncbi\|.+\|(.+)\|/;
                        $ncbi_taxid = $1;
                        $ncbi_tax =~ s/^ncbi\|.+\|.+\|//;
                        $ncbi_tax =~ s/[^A-Za-z0-9; ]/_/g;    
                    }
                else {  $ncbi_tax = "Unknown";
                        $ncbi_taxid = "NA";
                     }
                my $silva_tax;
                if ($split_line[9] ne "")
                    {   $silva_tax = $split_line[9];
                        $silva_tax =~ s/^silva\|.+\|.+\|//;
                        $silva_tax =~ s/[^A-Za-z0-9; ]/_/g;
                    }
                else {$silva_tax = "Unknown";}
                $SILVATAX{$asv}{'silva'} = $silva_tax;
                $SILVATAX{$asv}{'ncbi'} = $ncbi_tax;
                $SILVATAX{$asv}{'ncbi_taxid'} = $ncbi_taxid;
                if ($split_line[3] > 1)
                    {   $SILVATAX{$asv}{'multi'} = "TRUE";
                    }
                else
                    {   $SILVATAX{$asv}{'multi'} = "FALSE";
                    }
            }
        
        #Clean up taxonomy
        my @taxIDS;
        foreach my $i (sort keys %SILVATAX)
            {   unless ($SILVATAX{$i}{'ncbi_taxid'} eq "NA")
                    {   push(@taxIDS, $SILVATAX{$i}{'ncbi_taxid'});
                    }
            }
        my @uniq_taxids = uniq @taxIDS;
        
        open(NCBI, ">ncbiTaxonomyIDs.txt");
        foreach my $i (@uniq_taxids)
            {   print NCBI "$i\n";
            }
        close(NCBI);
        
        system("taxonkit lineage ncbiTaxonomyIDs.txt | awk '".'$2!=""'."' > taxonkit_out.txt");
        system("taxonkit reformat taxonkit_out.txt | cut -f1,3 > reformatted_taxonkit_out.txt");
        system("perl ".$options{m}."/assets/fillIn_taxonkit.pl -i reformatted_taxonkit_out.txt > reformatted_taxonkit_out.txt_temp");
        system("mv reformatted_taxonkit_out.txt reformatted_taxonkit_out_ORIGINAL.txt");
        system("cat reformatted_taxonkit_out.txt_temp | sed -E 's/[^A-Za-z0-9;[:blank:]]/_/g' > reformatted_taxonkit_out.txt");
        system("rm reformatted_taxonkit_out.txt_temp");

        open(NCBITAXIN, "<reformatted_taxonkit_out.txt") or die "\n\nThere is no reformatted_taxonkit_out.txt in the outdirectory!!\n\n";
        my @taxiddat = <NCBITAXIN>; close(NCBITAXIN);
        my %TAXONKIT;
        foreach my $line (@taxiddat)
            {   chomp($line);
                my @splitline = split('\t', $line);
                my $ncbi_taxid = $splitline[0];
                my $ncbi_taxonkit_format;
                if (exists $splitline[1])
                    {   $ncbi_taxonkit_format = $splitline[1];
                    }
                else
                    {   $ncbi_taxonkit_format = "Unknown";
                    }
                $TAXONKIT{$ncbi_taxid}{'taxonomystring'} = $ncbi_taxonkit_format;
            }
      
        foreach my $i (sort keys %SILVATAX)
            {   if (exists $TAXONKIT{$SILVATAX{$i}{'ncbi_taxid'}}{'taxonomystring'})
                    {   my $newNCBItax = $TAXONKIT{$SILVATAX{$i}{'ncbi_taxid'}}{'taxonomystring'};
                        my @splittax = split(';', $newNCBItax);
                        my $lastentry = $#splittax;
                        my $secondlastentry = $#splittax - 1;
                        if ($splittax[$lastentry] eq $splittax[$secondlastentry]." sp_" || $splittax[$lastentry] =~ m/uncultured/ || $splittax[$lastentry] =~ m/unidentified/ || $splittax[$lastentry] =~ m/metagenome/)
                            {   pop(@splittax);
                            }
                        if ($splittax[0] eq "NA" || $splittax[0] eq "")
                            {   $splittax[0] = "Unknown";
                            }
                        my $newnewNCBItax = join(';', @splittax);
                        $SILVATAX{$i}{'cleaned_ncbi'} = $newnewNCBItax;
                    }
                else
                    {   $SILVATAX{$i}{'cleaned_ncbi'} = "Unknown";
                    }
            }
        
        foreach my $i (sort keys %SILVATAX)
            {   unless (exists $SILVATAX{$i}{'cleaned_ncbi'})
                    {   die "\n\nSome of the taxa (like $i) are missing a cleaned NCBI taxonomy\n\n";
                    }
            }
        
        foreach my $i (sort keys %SILVATAX)
            {   my $inputtaxonomy = $SILVATAX{$i}{'silva'};
                my @silva_taxa = split(';', $inputtaxonomy);
                my @silva_taxa_copy = @silva_taxa;
                my $countSilva = scalar(@silva_taxa);
                my $kingdom = "NA";
                my $phylum = "NA";
                my $class = "NA";
                my $order = "NA";
                my $family = "NA";
                my $genus = "NA";
                my $species = "NA";
                foreach my $count (1..$countSilva)
                    {   my $comparison = join(';', @silva_taxa).";";
                        if (exists $SILVAREF{$comparison})
                            {   my $assignedLevel = $SILVAREF{$comparison}{'level'};
                                unless ($comparison =~ m/;uncultured;$/ || $comparison =~ m/;Unknown .+$/ || $comparison =~ m/;Incertae Sedis;$/)
                                    {   if ($assignedLevel eq "domain")
                                            {   $kingdom = $silva_taxa[$#silva_taxa];
                                            }
                                        elsif ($assignedLevel eq "phylum")
                                            {   $phylum = $silva_taxa[$#silva_taxa];
                                            }
                                        elsif ($assignedLevel eq "class")
                                            {   $class = $silva_taxa[$#silva_taxa];
                                            }
                                        elsif ($assignedLevel eq "order")
                                            {   $order = $silva_taxa[$#silva_taxa];
                                            }
                                        elsif ($assignedLevel eq "family")
                                            {   $family = $silva_taxa[$#silva_taxa];
                                            }
                                        elsif ($assignedLevel eq "genus")
                                            {   $genus = $silva_taxa[$#silva_taxa];
                                            }
                                        elsif ($assignedLevel eq "species")
                                            {   $species = $silva_taxa[$#silva_taxa];
                                            }
                                    }
                            }
                        else
                            {   unless ($comparison =~ m/^Unknown/)
                                    {   die "\n\nThe SILVA tax $comparison is missing from the reference file. Check that you have a current file.\n\n";
                                    }
                            }
                        pop(@silva_taxa);
                    }
                    
                my $comparison = join(';', @silva_taxa_copy).";"; #i.e. terminal subclass or other
                if (exists $SILVAREF{$comparison})
                    {   my $assignedLevel = $SILVAREF{$comparison}{'level'};
                        unless ($assignedLevel eq "domain" || $assignedLevel eq "phylum" || $assignedLevel eq "class" || $assignedLevel eq "order" || $assignedLevel eq "family" || $assignedLevel eq "genus" || $assignedLevel eq "species")
                            {   if ($species eq "NA" && $genus ne "NA")
                                    {   $species = $silva_taxa_copy[$#silva_taxa_copy]."__s";
                                    }
                                elsif ($species eq "NA" && $genus eq "NA" && $family ne "NA")
                                    {   $genus = $silva_taxa_copy[$#silva_taxa_copy]."__g";
                                    }
                                elsif ($species eq "NA" && $genus eq "NA" && $family eq "NA" && $order ne "NA")
                                    {   $family = $silva_taxa_copy[$#silva_taxa_copy]."__f";
                                    }
                                elsif ($species eq "NA" && $genus eq "NA" && $family eq "NA" && $order eq "NA" && $class ne "NA")
                                    {   $order = $silva_taxa_copy[$#silva_taxa_copy]."__o";
                                    }
                                elsif ($species eq "NA" && $genus eq "NA" && $family eq "NA" && $order eq "NA" && $class eq "NA" && $phylum ne "NA")
                                    {   $class = $silva_taxa_copy[$#silva_taxa_copy]."__c";
                                    }
                                elsif ($species eq "NA" && $genus eq "NA" && $family eq "NA" && $order eq "NA" && $class eq "NA" && $phylum eq "NA" && $kingdom ne "NA")
                                    {   $phylum = $silva_taxa_copy[$#silva_taxa_copy]."__p";
                                    }
                            }
                    }
                
                if ($genus ne "NA") #CHECK to FILL GAPS
                    {   if ($family eq "NA")
                            {   $family = $genus."__f";}
                        if ($order eq "NA")
                            {   $order = $family."__o";}
                        if ($class eq "NA")
                            {   $class = $order."__c";} 
                        if ($phylum eq "NA")
                            {   $phylum = $class."__p";}   
                    }
                if ($genus eq "NA")
                    {   if ($family ne "NA")
                            {   if ($order eq "NA")
                                    {   $order = $family."__o";}
                                if ($class eq "NA")
                                    {   $class = $order."__c";}
                                if ($phylum eq "NA")
                                    {   $phylum = $class."__p";}
                            }
                        if ($order ne "NA")
                            {   if ($class eq "NA")
                                    {   $class = $order."__c";}
                                if ($phylum eq "NA")
                                    {   $phylum = $class."__p";}
                            }
                        if ($class ne "NA")
                            {   if ($phylum eq "NA")
                                    {   $phylum = $class."__p";}
                            }
                    }
                
                my $final_silva_taxa = $kingdom.";".$phylum.";".$class.";".$order.";".$family.";".$genus.";".$species;
                $SILVATAX{$i}{'cleaned_silva'} = $final_silva_taxa;
                
                if ($comparison =~ m/^Unknown/)
                    {   $SILVATAX{$i}{'cleaned_silva'} = "Unknown;NA;NA;NA;NA;NA;NA";
                    }
            }
        
        foreach my $i (sort keys %SILVATAX)
            {   unless (exists $SILVATAX{$i}{'cleaned_silva'})
                    {   die "\n\nSome of the taxa (like $i) are missing a cleaned SILVA taxonomy\n\n";
                    }
            }
        
        #Populate $TAXA{$i}{'cleaned_merged'} filling in Euk assignments from NCBI
        open(EUKMERGE, ">merged_NCBI_SILVA_Eukaryotes_info.txt");
        print EUKMERGE "ASV\tOriginalSILVA\tNewMergedAssignmet\n";
        foreach my $i (sort keys %SILVATAX)
            {   my $ncbi_clean = $SILVATAX{$i}{'cleaned_ncbi'};
                my $silva_clean = $SILVATAX{$i}{'cleaned_silva'};
                if ($silva_clean =~ m/Bacteria\;Cyanobacteria\;Cyanobacteriia\;Chloroplast/)
                    {   print EUKMERGE "$i\t";
                        my @split_ncbi = split(';', $ncbi_clean);
                        if ($split_ncbi[0] eq "Eukaryota")
                            {   $SILVATAX{$i}{'cleaned_merged'} = $ncbi_clean;
                                print EUKMERGE "$silva_clean\t$ncbi_clean\n";
                            }
                        else
                            {   $SILVATAX{$i}{'cleaned_merged'} = "Eukaryota";
                                print EUKMERGE "$silva_clean\tEukaryota;Chloroplast\n";
                            }
                    }
                elsif ($silva_clean =~ m/Bacteria\;Proteobacteria\;Alphaproteobacteria\;Rickettsiales\;Mitochondria/)
                    {   print EUKMERGE "$i\t";
                        my @split_ncbi = split(';', $ncbi_clean);
                        if ($split_ncbi[0] eq "Eukaryota")
                            {   $SILVATAX{$i}{'cleaned_merged'} = $ncbi_clean;
                                print EUKMERGE "$silva_clean\t$ncbi_clean\n";
                            }
                        else
                            {   $SILVATAX{$i}{'cleaned_merged'} = "Eukaryota";
                                print EUKMERGE "$silva_clean\tEukaryota;Mitochondria\n";
                            }
                    }
                else
                    {   $SILVATAX{$i}{'cleaned_merged'} = $silva_clean;
                    }
            }
        close(EUKMERGE);
        
        open(SILVAIN2, "<$options{z}") or die "\n\nThere is no $options{z} file!!\n\n";
        my @silvadata2 = <SILVAIN2>; close(SILVAIN2);
        my %SILVACLUST;
        my $repseq;
        foreach my $line (@silvadata2)
            {   if ($line =~ m/^>ASV/)
                    {   $line =~ m/^>(ASV_[0-9]+)/;
                        $repseq = $1;
                    }
                else
                    {   $line =~ m/.+>(ASV_[0-9]+).+/;
                        my $clustermember = $1;
                        unless ($clustermember eq $repseq)
                            {   $SILVACLUST{$repseq}{'other'} .= $clustermember.";";
                            }
                    }
            }
        foreach my $i (sort keys %SILVACLUST)
            {   my $otherasv = $SILVACLUST{$i}{'other'};
                my @otherASVs = split(';', $otherasv);
                foreach my $j (@otherASVs)
                    {   if (exists $SILVATAX{$i}{'cleaned_merged'})
                            {   $SILVATAX{$j}{'cleaned_merged'} = $SILVATAX{$i}{'cleaned_merged'};
                            }
                    }
            }
        
        foreach my $i (sort keys %ASV)
            {   if (exists $SILVATAX{$i}{'cleaned_merged'})
                    {   $ASV{$i}{'finaltaxachoice'} = $SILVATAX{$i}{'cleaned_merged'};
                    }
                else
                    {   $ASV{$i}{'finaltaxachoice'} = "Unknown";
                    }
            }
    }
    
else
    {   #IMPORT blast reformatted taxonomy output
        open(IN2, "<$options{s}") or die "\n\nThere is no $options{s} file!!\n\n";
        my @dada_dat = <IN2>; close(IN2); shift(@dada_dat);
        foreach my $line (@dada_dat)
            {	chomp($line);
                my @data = split('\t', $line);
                my $blah = $data[0]; chomp($blah);
                $ASV{$blah}{'perchit'} = $data[1];                                                                         
                $ASV{$blah}{'lenhit'} = $data[2];
                $ASV{$blah}{'taxid'} = $data[3];
                $ASV{$blah}{'correction'} = $data[5];
                my $filterinfo = $options{f};
                my @splitfilter = split(',', $filterinfo);
                if ($data[1] >= $splitfilter[0])
                    {	$ASV{$blah}{'confidence'} = "SPECIES";	
                    }
                if ($data[1] < $splitfilter[0] && $data[1] >= $splitfilter[1])
                    {	$ASV{$blah}{'confidence'} = "GENUS";
                    }
                if ($data[1] < $splitfilter[1] && $data[1] >= $splitfilter[2])
                    {	$ASV{$blah}{'confidence'} = "FAMILY";
                    }
                if ($data[1] < $splitfilter[2] && $data[1] >= $splitfilter[3])
                    {	$ASV{$blah}{'confidence'} = "ORDER";
                    }
                if ($data[1] < $splitfilter[3] && $data[1] >= $splitfilter[4])
                    {	$ASV{$blah}{'confidence'} = "CLASS";
                    }
                if ($data[1] < $splitfilter[4] && $data[1] >= $splitfilter[5])
                    {	$ASV{$blah}{'confidence'} = "PHYLUM";
                    }
                if ($data[1] < $splitfilter[5])
                    {	$ASV{$blah}{'confidence'} = "NOCONFIDENCE";
                    }
            }
        
        #IMPORT reformatted taxonkit taxonomy out
        open(IN3, "<$options{t}") or die "\n\nThere is no $options{t} file!!\n\n";
        my @taxon_dat = <IN3>; close(IN3);
        foreach my $lines (@taxon_dat)
            {	chomp($lines);
                my @taxon_split = split('\t', $lines);
                my $taxonclean = $taxon_split[0];
                chomp($taxonclean);
                my $cleantaxastring = $taxon_split[1];
                chomp($cleantaxastring);
                $TAXON{$taxonclean}{'taxastring'} = $cleantaxastring;
            }
        
        #IMPORT common names
        open(INT, "<$options{c}") or die "\n\nThere is no $options{c} file!!\n\n";
        @commondat = <INT>; close(INT);
        foreach my $linet (@commondat)
            {	chomp($linet);
                my @common_split = split('\t', $linet);
                my $taxid_clean = $common_split[0];
                chomp($taxid_clean);
                my $common_name_clean = $common_split[2];
                chomp($common_name_clean);
                if (exists $TAXON{$taxid_clean})
                    {	$TAXON{$taxid_clean}{'common_name'} = $common_name_clean;
                    }
            }
            
        #Basic filtering of taxonomy to:
        #1) replace ;;;;;;;string with Environmental unknown
        #2) replace and double ;; with ;@@; for later trimming and detection
        #3) replace Eukarya;; with Environmental unknown
        
        foreach my $i (sort keys %TAXON)
            {	my $taxastring = $TAXON{$i}{'taxastring'};
                if ($taxastring =~ m/^;/)
                    {	$TAXON{$i}{'taxastring'} = "Environmental Unknown;@@;@@;@@;@@;@@;@@";
                        $taxastring = $TAXON{$i}{'taxastring'};
                    }
                if ($taxastring =~ m/;;/)
                    {	my @splitstring = split(';', $taxastring);
                        if ($#splitstring == 6)
                            {	foreach my $j (0..$#splitstring)
                                    {	if ($splitstring[$j] eq '')
                                            {	$splitstring[$j] = "@@";
                                                my $next = $j + 1;
                                                foreach my $k ($next..$#splitstring)
                                                    {	$splitstring[$k] = "@@";
                                                    }
                                                last;	
                                            }
                                    }
                                my $newgoodats = join(';', @splitstring);
                                if ($newgoodats !~ m/^Eukaryota;@/)
                                    {	$TAXON{$i}{'taxastring'} = $newgoodats;
                                    }
                                else {$TAXON{$i}{'taxastring'} = "Environmental Unknown;@@;@@;@@;@@;@@;@@";}
                            }
                        elsif ($#splitstring < 6) #TODO Think about this rare case more
                            {   foreach my $j (0..6)
                                    {	if (exists $splitstring[$j] && $splitstring[$j] eq '')
                                            {	$splitstring[$j] = "@@";
                                                my $next = $j + 1;
                                                foreach my $k ($next..6)
                                                    {	$splitstring[$k] = "@@";
                                                    }
                                                last;	
                                            }
                                        elsif (exists $splitstring[$j] && $splitstring[$j] ne '')
                                            {}
                                        else
                                            {   $splitstring[$j] = "@@";
                                                my $next = $j + 1;
                                                foreach my $k ($next..6)
                                                    {	$splitstring[$k] = "@@";
                                                    }
                                                last;	
                                            }
                                    }
                                my $newgoodats = join(';', @splitstring);
                                if ($newgoodats !~ m/^Eukaryota;@/)
                                    {	$TAXON{$i}{'taxastring'} = $newgoodats;
                                    }
                                else {$TAXON{$i}{'taxastring'} = "Environmental Unknown;@@;@@;@@;@@;@@;@@";}
                            }
                        else {die "\nYou forgot to load the reformatted taxonkit file with K/P/C/O/F/G/S only!\n\n";}
                    }
            }
        
        #CHOOSE TAXONOMY for output
        open(OUT_MULTI, ">".$options{n}."_singleBlastHits_with_MULTItaxid.txt");
        print OUT_MULTI "ASV\tMulti_taxID_entry\tUSED\n";
        foreach my $i (sort keys %ASV)
            {   my $passTaxTest = "FALSE";
                if (exists $ASV{$i}{'taxid'})	#TEST exists taxonomic assignment for taxid
                {	my $asv_taxid = $ASV{$i}{'taxid'};
                    my @tax_list;
                    if ($asv_taxid =~ m/\;/ && $asv_taxid =~ m/\,/)
                        {   my @multi = split(',', $asv_taxid);
                            my $totalMulti = scalar(@multi);
                            my $totalSemi = 0;
                            foreach my $j (@multi)
                                {	$j =~ s/\ //;
                                    if ($j =~ m/\;/)
                                        {   $totalSemi += 1;
                                        }
                                }
                            if ($totalSemi < $totalMulti) #purge semicolon hits to remove sequences with multiple potential taxonomies
                                {   foreach my $j (@multi)
                                        {	$j =~ s/\ //;
                                            unless ($j =~ m/\;/)
                                                {	push(@tax_list, $j);
                                                }
                                        }
                                }
                            elsif ($totalSemi == $totalMulti)
                                {   foreach my $j (@multi)
                                        {	$j =~ s/\ //;
                                            if ($j =~ m/\;/)
                                                {	my @multipli = split(';', $j);
                                                    foreach my $k (@multipli)
                                                        {	push(@tax_list, $k);	
                                                        }
                                                }
                                            else {push(@tax_list, $j);}
                                        }
                                }
                        }
                    elsif ($asv_taxid =~ m/\;/ && $asv_taxid !~ m/\,/)
                        {   my @multi = split(';', $asv_taxid);
                            foreach my $j (@multi)
                                {   push(@tax_list, $j);   
                                }
                        }
                    elsif ($asv_taxid =~ m/\,/ && $asv_taxid !~ m/\;/)
                        {   my @multi = split(',', $asv_taxid);
                            foreach my $j (@multi)
                                {	$j =~ s/\ //;
                                    push(@tax_list, $j);
                                }
                        }
                    else {push(@tax_list, $asv_taxid);}
                    my @new_tax_list;
                    foreach my $keepTaxHits (@tax_list)
                        {   if (exists $TAXON{$keepTaxHits})
                                {   push(@new_tax_list, $keepTaxHits);
                                }
                        }
                    @tax_list = @new_tax_list;
                    if (scalar @tax_list > 0)
                        {   $passTaxTest = "TRUE";
                        }
                }
                if (exists $ASV{$i}{'taxid'} && $passTaxTest eq "TRUE")	#CREATE list of taxids to consider
                {	my $asv_taxid = $ASV{$i}{'taxid'};
                    my @tax_list;
                    if ($asv_taxid =~ m/\;/ && $asv_taxid =~ m/\,/)
                        {   my @multi = split(',', $asv_taxid);
                            my $totalMulti = scalar(@multi);
                            my $totalSemi = 0;
                            foreach my $j (@multi)
                                {	$j =~ s/\ //;
                                    if ($j =~ m/\;/)
                                        {   $totalSemi += 1;
                                        }
                                }
                            if ($totalSemi < $totalMulti) #purge semicolon hits to remove sequences with multiple potential taxonomies
                                {   foreach my $j (@multi)
                                        {	$j =~ s/\ //;
                                            unless ($j =~ m/\;/)
                                                {	push(@tax_list, $j);
                                                }
                                        }
                                    print OUT_MULTI "$i\t$asv_taxid\tFALSE\n";
                                }
                            elsif ($totalSemi == $totalMulti)
                                {   foreach my $j (@multi)
                                        {	$j =~ s/\ //;
                                            if ($j =~ m/\;/)
                                                {	my @multipli = split(';', $j);
                                                    foreach my $k (@multipli)
                                                        {	push(@tax_list, $k);	
                                                        }
                                                }
                                            else {push(@tax_list, $j);}
                                        }
                                    print OUT_MULTI "$i\t$asv_taxid\tTRUE\n";
                                }
                        }
                    elsif ($asv_taxid =~ m/\;/ && $asv_taxid !~ m/\,/)
                        {   my @multi = split(';', $asv_taxid);
                            foreach my $j (@multi)
                                {   push(@tax_list, $j);   
                                }
                            print OUT_MULTI "$i\t$asv_taxid\tTRUE\n";
                        }
                    elsif ($asv_taxid =~ m/\,/ && $asv_taxid !~ m/\;/)
                        {   my @multi = split(',', $asv_taxid);
                            foreach my $j (@multi)
                                {	$j =~ s/\ //;
                                    push(@tax_list, $j);
                                }
                        }
                    else {push(@tax_list, $asv_taxid);}
                    my @new_tax_list;
                    foreach my $keepTaxHits (@tax_list)
                        {   if (exists $TAXON{$keepTaxHits})
                                {   push(@new_tax_list, $keepTaxHits);
                                }
                        }
                    @tax_list = @new_tax_list;
                    my $choiceindicator = 0;
                    my $commonnameYES = 0;
                    my $desiredcommonname;
                    my $choice; # SET the first taxid in the list as the first choice
                    if ($ASV{$i}{'confidence'} eq "SPECIES") # SET depth of first choice based on confidence (%)
                        {   my $tap = $TAXON{$tax_list[0]}{'taxastring'};
                            if ($tap !~ m/@@/)
                                {	$choice = $TAXON{$tax_list[0]}{'taxastring'};
                                    $choiceindicator = 1;
                                    if (exists $TAXON{$tax_list[0]}{'common_name'})
                                        {	$commonnameYES = 1;
                                            $desiredcommonname = $TAXON{$tax_list[0]}{'common_name'};
                                        }
                                }
                            else
                                {	my $newbie = $TAXON{$tax_list[0]}{'taxastring'};
                                    my @split_newbie = split(';', $newbie);
                                    my @pleasework;
                                    foreach my $position (0..$#split_newbie)
                                        {	if ($split_newbie[$position] eq "@@")
                                                {	last;
                                                }
                                            else {
                                                push(@pleasework, $split_newbie[$position]);
                                            }	
                                        }
                                    my $newwithout = join(';', @pleasework);
                                    $choice = $newwithout;
                                }
                        }
                    elsif ($ASV{$i}{'confidence'} eq "NOCONFIDENCE")
                        {	$choice = "Unknown";
                        }
                    else
                        {   my $firstchoice = $TAXON{$tax_list[0]}{'taxastring'};
                            chomp($firstchoice);
                            my @splitmodtax = split(';', $firstchoice);
                            my @conflist = qw | GENUS FAMILY ORDER CLASS PHYLUM |;
                            my $numpops;
                            foreach my $conf (0..$#conflist)
                                {   if ($ASV{$i}{'confidence'} eq $conflist[$conf])
                                        {   $numpops = $conf;
                                        }
                                }
                            foreach my $popnum (0..$numpops)
                                {   pop(@splitmodtax);
                                }
                            if ($splitmodtax[$#splitmodtax] ne "@@")
                                {	my $newmodtax = join(';', @splitmodtax);
                                    $choice = $newmodtax;
                                    $choiceindicator = 1;
                                }
                            else
                                {	my @pleasework;
                                    foreach my $position (0..$#splitmodtax)
                                        {	if ($splitmodtax[$position] eq "@@")
                                                {	last;
                                                }
                                            else {
                                                push(@pleasework, $splitmodtax[$position]);
                                            }	
                                        }
                                    my $newwithout = join(';', @pleasework);
                                    $choice = $newwithout;
                                }
                        }
                    
                    if ($#tax_list > 0) # SET COMPARISON taxonomy if taxid list has more than one
                        {	foreach my $heck (1..$#tax_list)
                                {	my $compare;
                                    my $commonnameCOMP = 0;
                                    my $compindicator = 0;
                                    if ($ASV{$i}{'confidence'} eq "SPECIES")
                                        {	my $tap = $TAXON{$tax_list[$heck]}{'taxastring'};
                                            if ($tap !~ m/@@/)
                                                {	$compare = $TAXON{$tax_list[$heck]}{'taxastring'};
                                                    $compindicator = 1;
                                                    if (exists $TAXON{$tax_list[$heck]}{'common_name'})
                                                        {	$commonnameCOMP = 1;
                                                        }
                                                }
                                            else
                                                {	if ($choiceindicator == 0)	# any ;; would end up here
                                                        {	my $newbie = $TAXON{$tax_list[$heck]}{'taxastring'};
                                                            my @split_newbie = split(';', $newbie);
                                                            my @pleasework;
                                                            foreach my $position (0..$#split_newbie)
                                                                {	if ($split_newbie[$position] eq "@@")
                                                                        {	last;
                                                                        }
                                                                    else {
                                                                        push(@pleasework, $split_newbie[$position]);
                                                                    }	
                                                                }
                                                            my $newwithout = join(';', @pleasework);
                                                            $compare = $newwithout;
                                                        }	
                                                    if ($choiceindicator == 1)	# only if no wonky tax assignment
                                                        {	$compare = $choice;
                                                        }
                                                    
                                                }
                                        }
                                    elsif ($ASV{$i}{'confidence'} eq "NOCONFIDENCE")
                                        {	$compare = "Unknown";
                                        }
                                    else
                                        {   my $firstchoice = $TAXON{$tax_list[$heck]}{'taxastring'};
                                            chomp($firstchoice);
                                            my @splitmodtax = split(';', $firstchoice);
                                            my @conflist = qw | GENUS FAMILY ORDER CLASS PHYLUM |;
                                            my $numpops;
                                            foreach my $conf (0..$#conflist)
                                                {   if ($ASV{$i}{'confidence'} eq $conflist[$conf])
                                                        {   $numpops = $conf;
                                                        }
                                                }
                                            foreach my $popnum (0..$numpops)
                                                {   pop(@splitmodtax);
                                                }
                                            if ($splitmodtax[$#splitmodtax] ne "@@")
                                                {	my $newmodtax = join(';', @splitmodtax);
                                                    $compare = $newmodtax;
                                                    $compindicator = 1;
                                                }
                                            else
                                                {	if ($choiceindicator == 0)
                                                        {	my @pleasework;
                                                            foreach my $position (0..$#splitmodtax)
                                                                {	if ($splitmodtax[$position] eq "@@")
                                                                        {	last;
                                                                        }
                                                                    else {
                                                                        push(@pleasework, $splitmodtax[$position]);
                                                                    }	
                                                                }
                                                            my $newwithout = join(';', @pleasework);
                                                            $compare = $newwithout;
                                                        }
                                                    if ($choiceindicator == 1)
                                                        {	$compare = $choice;
                                                        }
                                                }
                                        }
                                        
                                    my @new_choice;
                                    #print "$choice\t$compare\n";
                                    unless ($compindicator == 1 && $choiceindicator == 0)
                                        {
                                            if ($choice =~ m/;/ && $compare =~ m/;/)
                                                {	my @choice_array = split(';', $choice);
                                                    my @compare_array = split(';', $compare);
                                                    foreach my $s (0..$#choice_array)
                                                        {	if (exists $compare_array[$s])
                                                                {	if ($choice_array[$s] eq $compare_array[$s])
                                                                        {push(@new_choice, $choice_array[$s]);
                                                                        }
                                                                    else {last;}
                                                                }
                                                        }
                                                }
                                            if ($choice =~ m/;/ && $compare !~ m/;/)
                                                {	my @choice_array = split(';', $choice);
                                                    my @compare_array;
                                                    push(@compare_array,$compare);
                                                    foreach my $s (0)
                                                        {	if ($choice_array[$s] eq $compare_array[$s])
                                                                {push(@new_choice, $choice_array[$s]);
                                                                }
                                                            else {last;}
                                                        }
                                                }
                                            if ($choice !~ m/;/ && $compare =~ m/;/)
                                                {	my @choice_array;
                                                    push(@choice_array, $choice);
                                                    my @compare_array = split(';', $compare);
                                                    foreach my $s (0)
                                                        {	if ($choice_array[$s] eq $compare_array[$s])
                                                                {push(@new_choice, $choice_array[$s]);
                                                                }
                                                            else {last;}
                                                        }
                                                }
                                            if ($choice !~ m/;/ && $compare !~ m/;/)
                                                {	if ($choice eq $compare)
                                                        {	push(@new_choice, $choice);
                                                        }
                                                
                                                }
                                            if (scalar(@new_choice) == 0)
                                                {$choice = "Unknown";}
                                            else {	$choice = join(';', @new_choice);
                                                if ($#new_choice != 6)
                                                    {	$commonnameYES = 0;
                                                    }
                                                if ($#new_choice == 6)
                                                    {	if ($commonnameCOMP == 1 && $commonnameYES == 1)
                                                            {	my $choice_common = $desiredcommonname;
                                                                my $compare_common = $TAXON{$tax_list[$heck]}{'common_name'};
                                                                if ($choice_common ne $compare_common)
                                                                    {	$commonnameYES = 0;
                                                                    }
                                                                
                                                            }
                                                        if ($commonnameCOMP == 1 && $commonnameYES == 0)
                                                            {	$desiredcommonname = $TAXON{$tax_list[$heck]}{'common_name'};
                                                                $commonnameYES = 1;
                                                            }
                                                    }
                                                  }
                                        }
                                    if ($compindicator == 1 && $choiceindicator == 0)
                                        {	$choice = $compare;
                                            $choiceindicator = 1;
                                            if ($commonnameCOMP == 1)
                                                {	$desiredcommonname = $TAXON{$tax_list[$heck]}{'common_name'};
                                                    $commonnameYES = 1;
                                                }
                                        }
                                    
                                }
                        }
                    $ASV{$i}{'finaltaxachoice'} = $choice;
                    if ($commonnameYES == 1)
                        {	$ASV{$i}{'common_name'} = $desiredcommonname;
                        }
                }
                else {$ASV{$i}{'finaltaxachoice'} = "Unknown";}
            }
        close(OUT_MULTI);
    }

##Outputs
foreach my $i (0..$#sample_headers)
	{	my $sample = $sample_headers[$i];
		chomp($sample);
		open($sample, ">".$sample."_KRONA.txt");
		#print $sample "count\ttaxa\n";
	}

if ($options{d}) {
foreach my $i (0..$#sample_headers)
	{	my $sample = $sample_headers[$i];
		chomp($sample);
		my $sample_ig = $sample."_IGNORE";
		open($sample_ig, ">".$sample_ig."_KRONA.txt");
		#print $sample "count\ttaxa\n";
	}
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

if ($options{d}) {
open(OUT_IG, ">".$options{n}."_IGNORE_allin_KRONA.txt");
print OUT_IG "Sample\tASV\tcount\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\t(common name)\n";
open(ASVTAX_IG, ">".$options{n}."_IGNORE_asvTaxonomyTable.txt");
print ASVTAX_IG "ASV\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n";
foreach my $i (sort keys %ASV)
	{	unless ($ignore_ASV_string =~ m/_${i}_/)
		{
        
        my $taxput = $ASV{$i}{'finaltaxachoice'};
		my @taxputs = split(';', $taxput);
        print ASVTAX_IG "$i";
        foreach my $l (@taxputs) {
			print ASVTAX_IG "\t$l";}
		if ($#taxputs < 6)
			{	my $actualtaxdepth = $#taxputs + 1;
				foreach my $printtab ($actualtaxdepth..6)
					{	print ASVTAX_IG "\tNA";
					}
			}
        print ASVTAX_IG "\n";
		foreach my $j (0..$#sample_headers)
			{	my $hello = $sample_headers[$j];
				chomp($hello);
				my $hello_ig = $hello."_IGNORE";
				unless ($ASV{$i}{$sample_headers[$j]} == 0)
				{
				print OUT_IG "$sample_headers[$j]\t";
				print OUT_IG "$i\t";
				print OUT_IG "$ASV{$i}{$sample_headers[$j]}";
				foreach my $l (@taxputs) {
				print OUT_IG "\t$l";}
				if ($#taxputs < 6)
					{	my $actualtaxdepth = $#taxputs + 1;
						foreach my $printtab ($actualtaxdepth..6)
							{	print OUT_IG "\t";
							}
					}
				if (exists $ASV{$i}{'common_name'})
					{	print OUT_IG "\t(".$ASV{$i}{'common_name'}.")";
					}
				else {print OUT_IG "\t";}
				print OUT_IG "\n";
				print $hello_ig "$ASV{$i}{$sample_headers[$j]}";
				foreach my $l (@taxputs) {print $hello_ig "\t$l";}
				if (exists $ASV{$i}{'common_name'})
					{	print $hello_ig " (".$ASV{$i}{'common_name'}.")";
					}
				print $hello_ig "\n";
				}
			}
		}
	}
close(OUT_IG);
close(ASVTAX_IG);

foreach my $i (0..$#sample_headers)
	{	my $sample = $sample_headers[$i];
		chomp($sample);
		my $sample_ig = $sample."_IGNORE";
		close($sample_ig);
	}
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

if ($options{d})
	{ my @bartaxa_ig;
	  foreach my $j (0..$#sample_headers)
		{	foreach my $i (sort keys %ASV)
				{	unless ($ASV{$i}{$sample_headers[$j]} == 0 || $ignore_ASV_string =~ m/_${i}_/)
						{	my $pushchoice = $ASV{$i}{'finaltaxachoice'};
							chomp($pushchoice);
							if ($pushchoice =~ m/;/)
								{	my @splittingthis = split(';', $pushchoice);
									my $lastassignment = $splittingthis[$#splittingthis];
									chomp($lastassignment);
									push(@bartaxa_ig, $lastassignment);
								}
							else {push(@bartaxa_ig, $pushchoice);}
						}
				}
		}
	@uniq_bartaxa_ig = uniq @bartaxa_ig;
	}


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

unless ($options{e})
    {   open(COT, ">".$options{n}."_taxid_to_commonname_ALL.txt");
        foreach my $linet (@commondat)
            {	chomp($linet);
                my @common_split = split('\t', $linet);
                my $taxid_clean = $common_split[0];
                chomp($taxid_clean);
                my $common_name_clean = $common_split[2];
                chomp($common_name_clean);
                if (exists $CommonName_Term{$taxid_clean})
                    {	print COT "$CommonName_Term{$taxid_clean}{'name'} (".$common_name_clean.")\t$taxid_clean\n";
                    }
            }
        close(COT);
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


if ($options{d})
{
open(BARCHART_ig, ">".$options{n}."_IGNORE_barchart.txt");
open(BARCHART_forR_ig, ">".$options{n}."_IGNORE_barchart_forR.txt");
print BARCHART_forR_ig "Value\tSample\tTerminalTaxa\n";
print BARCHART_ig "Sample\t";
foreach my $i (@uniq_bartaxa_ig)
	{	print BARCHART_ig "$i\t";
	}
print BARCHART_ig "\n";
foreach my $j (0..$#sample_headers)
	{	my $barsampleheader = $sample_headers[$j];
		chomp($barsampleheader);
		print BARCHART_ig "$barsampleheader\t";
		my %BARCHART;
		foreach my $i (sort keys %ASV)
			{	unless ($ASV{$i}{$sample_headers[$j]} == 0 || $ignore_ASV_string =~ m/_${i}_/)
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
		foreach my $uniq_taxa (@uniq_bartaxa_ig)
			{	my $hit = 0;
				foreach my $k (sort keys %BARCHART)
					{	if ($uniq_taxa eq $k)
							{	$hit = 1;
								print BARCHART_ig "$BARCHART{$k}\t";
                                print BARCHART_forR_ig "$BARCHART{$k}\t$sample_headers[$j]\t$uniq_taxa\n";
							}
					}
				if ($hit == 0)
					{	print BARCHART_ig "0\t";
					}
			}
		print BARCHART_ig "\n";
	}
close(BARCHART_ig);
close(BARCHART_forR_ig);
}


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

if ($options{d})
	{	my @kronasamples_ig;
		foreach my $i (0..$#sample_headers)
			{	my $hebs = $sample_headers[$i];
				chomp($hebs);
				my $newhebs = $hebs."_IGNORE_KRONA.txt";
				push(@kronasamples_ig, $newhebs);
			}
		my $printkronasamples_ig = join(' ', @kronasamples_ig);
		my $naming = $options{n};
		system("ImportText.pl -o ".$naming."_IGNORE_master_krona.html $printkronasamples_ig");
	}

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

if ($options{d})
{
open(NOUNKNOWN_ig, ">".$options{n}."_IGNORE_NO_UNKNOWNS_barchart.txt");
print NOUNKNOWN_ig "Sample\t";
foreach my $i (@uniq_bartaxa_ig)
	{	unless ($i eq "Unknown" || $i eq "Environmental Unknown")
			{	print NOUNKNOWN_ig "$i\t";
			}
	}
print NOUNKNOWN_ig "\n";
foreach my $j (0..$#sample_headers)
	{	my $barsampleheader = $sample_headers[$j];
		chomp($barsampleheader);
		print NOUNKNOWN_ig "$barsampleheader\t";
		my %BARCHART;
		foreach my $i (sort keys %ASV)
			{	unless ($ASV{$i}{$sample_headers[$j]} == 0 || $ignore_ASV_string =~ m/_${i}_/)
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
		foreach my $uniq_taxa (@uniq_bartaxa_ig)
			{	unless ($uniq_taxa eq "Unknown" || $uniq_taxa eq "Environmental Unknown")
					{	my $hit = 0;
						foreach my $k (sort keys %BARCHART)
							{	if ($uniq_taxa eq $k)
									{	$hit = 1;
										print NOUNKNOWN_ig "$BARCHART{$k}\t";
									}
							}
						if ($hit == 0)
							{	print NOUNKNOWN_ig "0\t";
							}
					}
			}
		print NOUNKNOWN_ig "\n";
	}
close(NOUNKNOWN_ig);

}

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





if ($options{d})
{
my %TAXHEAT_ig;

foreach my $i (sort keys %ASV)
	{   unless ($ignore_ASV_string =~ m/_${i}_/)
        {
        my $tax = $ASV{$i}{'finaltaxachoice'};
		chomp($tax);
		if ($tax =~ m/;/)
			{	my @split_tax = split(';', $tax);
				foreach my $j (0..$#split_tax)
					{	if (exists $TAXHEAT_ig{$split_tax[$j]})
							{	$TAXHEAT_ig{$split_tax[$j]}{'asvs'} .= ";".$i;
							}
						else
							{	$TAXHEAT_ig{$split_tax[$j]}{'asvs'} = $i;
								$TAXHEAT_ig{$split_tax[$j]}{'depth'} = $j + 1;
							}
					}
			}
		else
			{	if (exists $TAXHEAT_ig{$tax})
					{	$TAXHEAT_ig{$tax}{'asvs'} .= ";".$i;
					}
				else
					{	$TAXHEAT_ig{$tax}{'asvs'} = $i;
						$TAXHEAT_ig{$tax}{'depth'} = 1;
					}
			}
        }
	}

open(OUTDEEP1_ig, ">".$options{n}."_IGNORE_heatmap_multiASV.txt");
foreach my $i (sort keys %TAXHEAT_ig)
	{	if ($TAXHEAT_ig{$i}{'asvs'} =~ m/;/)
			{	print OUTDEEP1_ig "$i <<Depth ".$TAXHEAT_ig{$i}{'depth'}.">>\tSample\t";
				my $multiasv = $TAXHEAT_ig{$i}{'asvs'};
				my @multiasv_split = split(';', $multiasv);
				foreach my $j (@multiasv_split)
					{	print OUTDEEP1_ig "$j\t";
					}
				print OUTDEEP1_ig "\n";
				foreach my $k (0..$#sample_headers)
					{	print OUTDEEP1_ig "$i <<Depth ".$TAXHEAT_ig{$i}{'depth'}.">>\t";
						print OUTDEEP1_ig "$sample_headers[$k]\t";
						foreach my $j (@multiasv_split)
							{	print OUTDEEP1_ig "$ASV{$j}{$sample_headers[$k]}\t";
							}
						print OUTDEEP1_ig "\n";
					}
			}
	}
close(OUTDEEP1_ig);

}

if ($options{e})
    {   open(UNKNOWNS, ">".$options{n}."_unknown_asvids.txt");
        foreach my $i (sort keys %ASV)
            {   my $taxonunknown = $ASV{$i}{'finaltaxachoice'};
                chomp($taxonunknown);
                if ($taxonunknown eq "Unknown" || $taxonunknown eq "Environmental Unknown")
                    {   print UNKNOWNS "$i\t\n";
                    }
            }
        close(UNKNOWNS);
    }
else
    {   open(UNKNOWNS, ">".$options{n}."_unknown_asvids.txt");
        foreach my $i (sort keys %ASV)
            {	my $taxonunknown = $ASV{$i}{'finaltaxachoice'};
                chomp($taxonunknown);
                if ($taxonunknown eq "Unknown" || $taxonunknown eq "Environmental Unknown")
                    {	print UNKNOWNS "$i\t";
                        if (exists $ASV{$i}{'confidence'})
                            {	print UNKNOWNS "$ASV{$i}{'confidence'}\t";
                            }
                        if (exists $ASV{$i}{'taxid'})
                            {	print UNKNOWNS "$ASV{$i}{'taxid'}\t";
                                my $string = $ASV{$i}{'taxid'};
                                my @unknown_tax_list;
                                if ($string =~ m/\;/ && $string =~ m/\,/)
                                    {   my @multi = split(',', $string);
                                        my $totalMulti = scalar(@multi);
                                        my $totalSemi = 0;
                                        foreach my $j (@multi)
                                            {	$j =~ s/\ //;
                                                if ($j =~ m/\;/)
                                                    {   $totalSemi += 1;
                                                    }
                                            }
                                        if ($totalSemi < $totalMulti) #purge semicolon hits to remove sequences with multiple potential taxonomies
                                            {   foreach my $j (@multi)
                                                    {	$j =~ s/\ //;
                                                        unless ($j =~ m/\;/)
                                                            {	push(@unknown_tax_list, $j);
                                                            }
                                                    }
                                            }
                                        elsif ($totalSemi == $totalMulti)
                                            {   foreach my $j (@multi)
                                                    {	$j =~ s/\ //;
                                                        if ($j =~ m/\;/)
                                                            {	my @multipli = split(';', $j);
                                                                foreach my $k (@multipli)
                                                                    {	push(@unknown_tax_list, $k);	
                                                                    }
                                                            }
                                                        else {push(@unknown_tax_list, $j);}
                                                    }
                                            }
                                    }
                                elsif ($string =~ m/\;/ && $string !~ m/\,/)
                                    {   my @multi = split(';', $string);
                                        foreach my $j (@multi)
                                            {   push(@unknown_tax_list, $j);   
                                            }
                                    }
                                elsif ($string =~ m/\,/ && $string !~ m/\;/)
                                    {   my @multi = split(',', $string);
                                        foreach my $j (@multi)
                                            {	$j =~ s/\ //;
                                                push(@unknown_tax_list, $j);
                                            }
                                    }
                                if (scalar(@unknown_tax_list) > 1)
                                    {   foreach my $entry (@unknown_tax_list)
                                            {   if (exists $TAXON{$entry}{'taxastring'})
                                                    {print UNKNOWNS "$TAXON{$entry}{'taxastring'}\t";}
                                                else
                                                    {print UNKNOWNS "TaxaStringDeleted_or_DoesNotExist\t";}
                                            }
                                        print UNKNOWNS "\n";
                                    }
                                else
                                {   if (exists $TAXON{$string}{'taxastring'})
                                        {   print UNKNOWNS "$TAXON{$string}{'taxastring'}\n";
                                        }
                                    else
                                        {   print UNKNOWNS "TaxaStringDeleted_or_DoesNotExist\n";
                                        }
                                }
                            }
                        else
                            {	print UNKNOWNS "NO TAXID ASSIGNMENT\n";
                            }
                    }
            }
        close(UNKNOWNS);
    }

if ($options{d})
	{	system("mkdir IGNORING_ASVs");
		system("mv *IGNORE* IGNORING_ASVs/");
	}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - S U B R O U T I N E S - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -
