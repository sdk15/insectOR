#!/usr/bin/perl

use strict;
use warnings FATAL => 'uninitialized';
use Getopt::Long qw(GetOptions);
use Tie::IxHash;
use HTML::Template;
use File::Basename;
use File::Spec;
use IO::File;
use POSIX qw(_exit);
use GD::Graph::bars;
use GD::Graph::Data;
use IO::Handle;
STDOUT->autoflush(1);
use lib File::Spec->catdir(
                File::Basename::dirname(File::Spec->rel2abs($0)),
                '..',
                'lib');
use My::utils qw(processAlignment processEveryAlignment findMet processFinalNew processSimilarity parseTMHMM parseHMMTOP parsePhobius makeConsensus assignColors joinAlignments three2one barplot);
use Cwd;
my $dir = getcwd;
my $accid= basename($dir);

#######################################################################################################################################################################################
##### INSECTOR tools path assignment ##################################################################################################################################################
##### Please provide correct $basepath as per the location of insectOR on your system #################################################################################################
#######################################################################################################################################################################################

my $basepath="/media/snehal/Data1/my_tools/insectOR/";
my $wiseconfigdir_path="$basepath/tools/wise2.4.1/wisecfg/";
my $hmmsearchpath="$basepath/tools/hmmer-3.1b2/src";
my $tmhmmpath="$basepath/tools/tmhmm-2.0c/bin";
my $genewisepath="$basepath/tools/wise2.4.1/bin";
my $hmmtoppath="$basepath/tools/hmmtop_2.1";
my $phobiuspath="$basepath/tools/phobius";
my $mastpath="$basepath/tools/meme_4.10.2/meme/bin";
my $mast_xslt_path="$basepath/tools/meme_4.10.2/etc/mast-to-html.xsl";
my $motifsfile="$basepath/tools/motifs/AfOr_motifs";
my $wiseconfigdir_path="$basepath/tools/wise2.4.1/wisecfg/";
my $gff_to_bed_loc="$basepath/tools/GFFtools-GX-master/";

my $usage=<<"END_MSG";

Usage: $0 -i <exonerate alignment output file> -s <genome sequence file> -q <query sequence file>
Please make sure that -
1) The exonerate is run using following parameters - 
   "exonerate --model protein2genome --maxintron <max intron size dependent on the organism> <protein query sequence fasta file> <genome fasta file> --showtargetgff TRUE"
2) The exonerate alignment contains 'Command line' and 'exonerate' in first line.

Example command - 
perl (path_to_insectOR/)insectOR/bin/scoreGenesOnScaffold.pl -i exonerate.txt -s seq.fasta -q query.fasta
perl (path_to_insectOR/)insectOR/bin/scoreGenesOnScaffold.pl -i exonerate.txt -s seq.fasta -q query.fasta -g ncbi.gff -tmh1 -tmh2 -tmh3 -p -m -mf motif.txt 


*** Input ***
-------------

* Mandatory
-exonerate_file|-i - Exonerate alignment file
-seq|-s            - Genome sequence (FASTA) file which was used to generate above exonerate file
-queryseq|-q       - OR/query sequence (FASTA) file which was used for generating above exonerate alignment

* Optional
-gff_file|-g       - User provided gene annotations (GFF format) with which InsectOR output will be compared
-cutoff|-c         - Alignment clusters are identified based on this cutoff. 
                     This is the minimum number of alignments needed at a nucleotide position for its inclusion into an alignment cluster. (Default: 1)
-lengthCutoff|-l   - Predicted proteins can be classified as complete or partial based on this cutoff. (Default: 300 amino acids)
-hmmsearch|-p      - Perform HMMSEARCH for insect olfactory receptor signature (PFAM 7tm_6 family) as additional validation.
-tmhmm|-tmh1       - Search for Transmembrane Helices (TMH) using TMHMM2 as additional validation.
-hmmtop|-tmh2      - Search for Transmembrane Helices (TMH) using HMMTOP2 as additional validation.
-phobius|-tmh3     - Search for Transmembrane Helices (TMH) using Phobius as additional validation.
                     (If all three TMH predictiors are selected, consensus TMH prediction will be performed.)
-mast|-m           - Search for known OR motifs in the predicted proteins using MAST motif tool
-motif_file|-mf    - Users can provide their own motifs (PSPM format of MEME) for MAST motif search into the predicted OR proteins. 
                     (If MAST option is checked without any file mentioned with this parameter, then default AfOR motif file will be used)
-helpmessage|-h    - Print this help message. Users can read README file to know more.


*** Output ***
--------------

...predictedORs.summary.txt    - Summary file of insectOR results
...OR_final_insector_table.txt - Table providing detailed information on each gene/fragment predicted by insectOR
...final_proteins.pep          - Protein sequence file of insectOR predicted ORs/fragments in FASTA format
...ORs.starRemoved.pep         - Protein sequence file of insectOR predicted ORs/fragments without any pseudogenizing elements in FASTA format
...ORs.cds                     - Nucleotide sequence file of insectOR predicted OR CDSs in FASTA format
...final_gff_file.gff          - Detailed gene annotations by insectOR in GFF format
...ORs_sorted.bed12            - Detailed gene annotations by insectOR in BED format

If gene annotation file from another resource is provided, following files are also generated - 
...gffcomparison                                 - Comparison of gene annotations from insectOR along with overlapping user provided gene annotations. 
...ORrelated_genesFromUserProvidedAnnotation.gff - A shorter version of user provided GFF file containing only overalapping annotations with those from this tool.
...gff.ORs_sorted.bed12                          - sorted BED formatted version of transcripts in the user-provided GFF file.

If HMMSEARCH option against 7tm_6 is selected, following files are also generated -
...ORs.starRemoved.pep.hmmsearchout              - HMMSEARCH output against 7tm_6 HMM file

If either of TMH prediction tools are used, following files are also generated -
...ORs.starRemoved.pep.tmhmmout if '-tmh1' is selected   - Output of TMHMM2 tool
...ORs.starRemoved.pep.hmmtopout if '-tmh2' is selected  - Output of HMMTOP2 tool
...ORs.starRemoved.pep.phobiusout if '-tmh3' is selected - Output of Phobius tool
...ORs.starRemoved.pep.consensusTMHpredout and ...ORs.starRemoved.pep.consensusTMHpred.html if all the three '-tmh1 -tmh2 -tmh3' are selected - Output of consensus TMH prediction tool.

If MAST motif search is selected, following files are also generated -
MAST.html - Motif search output file in HTML format
MAST.xml  - Motif search output file in HTML format


*** Contact ***
Please contact karpesnehal\@gmail.com or vikast\@ncbs.res.in to report any bugs.


END_MSG

#######################################################################################################################################################################################
##### INSECTOR input parameter assignment #############################################################################################################################################
##### Please do not change default settings hereafter. Specific parameters can be provided externally while executing insectOR ########################################################
#######################################################################################################################################################################################

my $filename=0;
my $cutoff = 1;
my $genRawFile=0;
my $gfffile;
my $seqfile=0;
my $queryFile=0;
my $lengthCutoff=300;

my $hmmsearch;
my $tmhmm;
my $hmmtop;
my $phobius;
my $mast;
my $helpmessage;

GetOptions(
	'base_path|b=s' => \$basepath,
    'exonerate_file|i=s' => \$filename,    # input exonerate file
    'cutoff|c=i' => \$cutoff,              # cutoff for editing gene clusters(How many minimum alignments are needed to declare a cluster)
    'genRawFile|f' => \$genRawFile,        # whether to generate raw files or not (0 or 1) - not very informative for users
    'gff_file|g=s' => \$gfffile,           # input GFF file containing known/predicted genome annotations from another resource
    'hmmsearch|p' => \$hmmsearch,          # select to perform hmmsearch against 7tm_6 protein family
    'tmhmm|tmh1' => \$tmhmm,			   # select to perform transmembrane helix prediction on predicted protein sequences using TMHMM
    'hmmtop|tmh2' => \$hmmtop,			   # select to perform transmembrane helix prediction on predicted protein sequences using HMMTOP
    'phobius|tmh3' => \$phobius,		   # select to perform transmembrane helix prediction on predicted protein sequences using Phobius
    'mast|m' => \$mast,					   # select to perform motif search using MAST tool from MEME suit
    'motif_file|mf=s' => \$motifsfile,	   # input motif file in PSPM format of MEME
    'seq|s=s' => \$seqfile,				   # input genome sequence file
    'queryseq|q=s' => \$queryFile,		   # input query protein / OR file
    'lengthCutoff|l=i' => \$lengthCutoff,  # protein length cut off to declare gene prediction to be complete
    'helpmessage|h' => \$helpmessage,
) or die $usage;

if(defined $helpmessage)
{
	die $usage;
}

if((! -e $filename) || (-z $filename))
{
	print "Exonerate alignment file does not exist or it is empty.\n"; die $usage;
}

if((! -e $seqfile) || (-z $seqfile))
{
	print "Genome sequence file does not exist or it is empty.\n"; die $usage;
}

if((! -e $queryFile) || (-z $queryFile))
{
	print "Query protein sequence file does not exist or it is empty.\n"; die $usage;
}

$cutoff=$cutoff-1;

my $filefolder;
if($filename=~/^(.*)\/.*$/)
{
	$filefolder=$1;
}
else
{
	$filefolder=".";
}


my $filename2=$filename;
$filename2=~s/\.//;

# Filecheck - if all these are files are of correct format
open(FH1,"<$filename");
$/="# --- END OF GFF DUMP ---
#

";

my @file=<FH1>;
close FH1;
$/="\n";

if(!($file[0]=~/Command line.*exonerate/))
{
	print "The input is not an exonerate alignment file.\n"; die $usage;
}

my %genomeSeqHash;
open(FH5,"<$seqfile");

my $ke="";
my $valu="";

while(<FH5>)
{
	chomp $_;
	my $lineInFile=$_;
	if($lineInFile=~/^>/)
	{
		if($ke ne "")
		{
			$genomeSeqHash{$ke}=$valu;
		}

		if($lineInFile=~/^(.*?)\s+/)
		{
			$ke=$1;
		}
		else
		{
			$ke=$lineInFile;
		}
		$valu="";
	}
	else
	{
		$valu.=$lineInFile;
	}
}
close FH5;
$genomeSeqHash{$ke}=$valu;

open(OUTN50,">$filename.modified_seqfile");
foreach my $keyys(keys %genomeSeqHash)
{
	if ($genomeSeqHash{$keyys} =~ m/^[ATGCNURYWSMKBHDV]*$/i)
	{
		$genomeSeqHash{$keyys}=~tr/[Uu]/T/;
		$genomeSeqHash{$keyys}=~tr/[RYWSMKBHDVrywsmkbhdv]/N/;
	}
	else
	{
		print "$keyys\nError in the sequence. Not a DNA sequence. It contains alphabets other than A,T,G,C,N,U,R,Y,W,S,M,K,B,H,D,V.\n"; die $usage;
	}

	print OUTN50 "$keyys\n$genomeSeqHash{$keyys}\n";
}
close OUTN50;
print "Degenerate characters like R,Y,W,S,M,K,B,H,D,V were replaced with N in the genomic DNA sequence.\n";


my %querySeqHash;
open(FH6,"<$queryFile");
my $qke="";
my $qvalu="";

while(<FH6>)
{
	chomp $_;
	my $qlineInFile=$_;
	if($qlineInFile=~/^>/)
	{
		if($qke ne "")
		{
			$querySeqHash{$ke}=$valu;
		}

		if($qlineInFile=~/^(.*?)\s+/)
		{
			$qke=$1;
		}
		else
		{
			$qke=$qlineInFile;
		}
		$qvalu="";
	}
	else
	{
		$qvalu.=$qlineInFile;
	}
}
close FH6;
$querySeqHash{$qke}=$qvalu;

foreach my $qkeyys(keys %querySeqHash)
{
	if ($querySeqHash{$qkeyys} =~ m/^[ACDEFGHIKLMNPQRSTVWYXZ]*$/i)
	{
		#print "$keyys\nFine\n";
	}
	else
	{
		print "$qkeyys\nError in the sequence. Not a protein sequence. It contains alphabets other than A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,X,Z.\n"; die $usage;
	}
}


#######################################################################################################################################################################################
##### INSECTOR - STEP 1 - Initialization of virtual scaffold with 0 at each position ##################################################################################################
#######################################################################################################################################################################################

##### Processing alignment file to generate a hash (%scaf_genes) of info on all alignments. Key of hash is a DNA scaffold and its values are all alignments on that scaffold. ####
##### Values in the array are start of target, end of the target, score, hit length, and full alignment entry ####

my %scaf_genes;
my $max_end=0;
my $min_start=2**53;
my %scaf_length;
my %scaf_offset;
my $vector_length=9;


print "\nProcessing alignment file...\n\n";
foreach my $i(@file)
{
#	if($i=~/Target: (.*?)(\:\[revcomp\])*\n.*?Raw score: (.*?)\n.*?Target range: (.*?) -> (.*?)\n.*?vulgar: (.*?)\n/ms)

	if($i=~/Query\: (.*?)\n.*?Target: (.*?)(\:\[revcomp\])*\n.*?Raw score: (.*?)\n.*?Query range: (.*?) -> (.*?)\n.*?vulgar: (.*?)\n.*?exonerate:protein2genome:local\s+gene\s+(.*?)\s+(.*?)\s+.*?\s+(.*?)\s+/ms)
	{
		my $quer=$1;
		my $target=$2;
		my $extra=$3;
		my $score=$4;
		my $query_start=$5+1;
		my $query_end=$6;
		my $vulgar=$7;
		my $start=$8;
		my $end=$9;
		my $direct=$10;
		my $hit_length=0;
		$target=~s/(.*?) .*$/$1/;
		$target=~s/ /~/g;
		$target=~s/,//g;
		#$target=~s/_//g;
		if(exists $scaf_length{$target}){$max_end=$scaf_length{$target}}
		else{$max_end=0}
		if(exists $scaf_offset{$target}){$min_start=$scaf_offset{$target}}
		else{$min_start=2**53}

		if($start<$end)
		{
			push(@{$scaf_genes{$target}},$start);
			push(@{$scaf_genes{$target}},$end);
			if ($end>$max_end){$max_end=$end}
			if ($start<$min_start){$min_start=$start}
		}
		else
		{
        	push(@{$scaf_genes{$target}},$end);
        	push(@{$scaf_genes{$target}},$start);
			if ($start>$max_end){$max_end=$start}
			if ($end<$min_start){$min_start=$end}
		}
		$scaf_length{$target}=$max_end;
		$scaf_offset{$target}=$min_start;
		push(@{$scaf_genes{$target}},$score);
		my @temp_vulgar=($vulgar=~/\sM (.*?) /g);
		foreach my $u(@temp_vulgar)
		{
			$hit_length+=$u;
		}
		push(@{$scaf_genes{$target}},$hit_length);
		push(@{$scaf_genes{$target}},$i);
		push(@{$scaf_genes{$target}},$query_start);
		push(@{$scaf_genes{$target}},$query_end);
		push(@{$scaf_genes{$target}},$direct);
		push(@{$scaf_genes{$target}},$quer);
	}
	else
	{
		print "No target found\n";  die $usage;
	}
}

print "Alignment file processed.\n\n";

####Initializing scaffold with zero score####

print "Initializing each DNA scaffold vector with score 0 for each nucleotide position...\n\n";

my %gene_loci;
foreach my $k(keys %scaf_length)
{
	#print "Initializing for $k with length $scaf_length{$k}\n";
	#print "Initializing for $k with offset $scaf_offset{$k}\n";
	@{$gene_loci{$k}}=(0) x ($scaf_length{$k} - $scaf_offset{$k} + 1);
}

my %scaf_length_bb=%scaf_length;


my %final_hits;
my %final_hits_suppl;

my %gene_clusters;
my $key_counter=1;

#######################################################################################################################################################################################
##### INSECTOR - STEP 2 and 3 - Increment the score at each alignment position, identify alignment clusters and best scoring query per cluster ########################################
#######################################################################################################################################################################################

#print "Scoring each nucleotide position per scaffold...\n";
print "Incrementing scaffold vector by score 1 for each alignment at respective nucleotide position...\n\n";

####Incrementing score by one for each alignment, identifying clusters and choosing the best query for each cluster with the subroutine processEveryAlignment()####
my @temporary=processEveryAlignment($filename, $cutoff, $vector_length, $genRawFile, \%scaf_genes, \%gene_loci, \%scaf_offset, \%scaf_length);
print "Scoring of nucleotide positions finished.\n\n";

#######################################################################################################################################################################################
##### INSECTOR - STEP 4 - Indentification of all alignments of the best queries with their repsective DNA scaffold sequences ##########################################################
#######################################################################################################################################################################################

my %scaf_scorebestquery_hits=%{$temporary[0]};
my %scaf_scorebestquery_all_hits;
my %scaf_scorebestquery_all_hits_suppl;
my $cou=1;

#####Retrieving all hits for best queries on a scaffold and their alignment details#####

#Again taking all alignments for scaffold $s and their related values in scaf_genes
foreach my $s (sort keys %scaf_scorebestquery_hits)
{	
	my $abc_ref=\@{$scaf_genes{$s}};
	my $out_stop = @{$scaf_genes{$s}};
	#Taking all the alignments for scaffold $s and with query $u and sorting them based on location - target start
	foreach my $u (sort keys %{$scaf_scorebestquery_hits{$s}})
	{

		my %forSortingHash;
		for (my $t=0;$t<$out_stop;$t=$t+$vector_length)
		{
			#$t should contain target start
			if(${$abc_ref}[$t+8] eq $u)
			{
				$forSortingHash{${$abc_ref}[$t]}=$t;
			}
		}

		
		#####Retrieving all hits for best queries on a scaffold sorted by target location#####
		foreach my $v(sort {$a <=> $b} keys %forSortingHash)
		{
			push(@{$scaf_scorebestquery_all_hits{$cou}},$s,0,1000000000,${$abc_ref}[$forSortingHash{$v}],${$abc_ref}[$forSortingHash{$v}+1],${$abc_ref}[$forSortingHash{$v}+2],${$abc_ref}[$forSortingHash{$v}+3],${$abc_ref}[$forSortingHash{$v}+5],${$abc_ref}[$forSortingHash{$v}+6],${$abc_ref}[$forSortingHash{$v}+7],${$abc_ref}[$forSortingHash{$v}+8],${$abc_ref}[$forSortingHash{$v}+4]);
			##### Process selected alignments and get their alignment details #####
			my @procAlignScoreTemp = processAlignment( ${$abc_ref}[ $forSortingHash{$v} + 4] );
			$scaf_scorebestquery_all_hits_suppl{$cou}{"QueryProt"}=$procAlignScoreTemp[0]; $scaf_scorebestquery_all_hits_suppl{$cou}{"TargetNuc"}=uc($procAlignScoreTemp[1]); $scaf_scorebestquery_all_hits_suppl{$cou}{"Similarity"}=$procAlignScoreTemp[2]; $scaf_scorebestquery_all_hits_suppl{$cou}{"TargetProt"}=$procAlignScoreTemp[3]; $scaf_scorebestquery_all_hits_suppl{$cou}{"CDS"}=$procAlignScoreTemp[4];
			$cou++;
		}


	}

}

#######################################################################################################################################################################################
##### INSECTOR - STEP 5 - Stitch the aligned and overlapping regions arising from the same query to decide plausible exon boundaries of the same gene #################################
##### INSECTOR - STEP 6 - Then, select the best stiched aligned query (and its alignement) out of multiple queries that might align at the same location ##############################
##### INSECTOR - STEP 7 - Remove other alignments at the same location and identify putative gene regions (P1) ########################################################################
#######################################################################################################################################################################################


#####Stitching alignments if same scaffold same query and overlapping hits. Also out of multiple best hits on scaffold from multiple protein queries - select one best scoring hit out of overlapping hits using processFinalNew() subroutine#####
(my $temphash,my $temphashsuppl)=processFinalNew(\%scaf_scorebestquery_all_hits,\%scaf_scorebestquery_all_hits_suppl, $filename, $lengthCutoff);
#####Core algorithm to select the best hits is over. #####

my %hash=%$temphash;
my %hash_suppl=%$temphashsuppl;
my %scafhash;

my $gene_count;

my $complete_count=0; my $partial_count=0;
my $start_present=0; my $start_absent=0;

foreach my $w(sort {$a <=> $b} keys %hash_suppl)
{
	$gene_count++;
	$scafhash{$hash_suppl{$w}{"Scaf"}}=1;
	if($hash_suppl{$w}{"Name"}=~/Complete/)
	{
		$complete_count++;
	}
	else
	{
		$partial_count++;
	}

	if($hash_suppl{$w}{"Name"}=~/StartCodonPresent/)
	{
		$start_present++;
	}
	else
	{
		$start_absent++;
	}
}

open(OUT6,">$filename.ORs.cds");
open(OUT7,">$filename.ORs.pep");

foreach my $w(sort {$a <=> $b} keys %hash_suppl)
{
	print OUT6 ">".$hash_suppl{$w}{"Name"}."\n".$hash_suppl{$w}{"TargetNuc"}."\n";
	print OUT7 ">".$hash_suppl{$w}{"Name"}."\n".$hash_suppl{$w}{"TargetProt"}."\n";
}
close OUT6; close OUT7;

open(OUT9,">$filename.OR_table.txt");


print OUT9 "Name\tScaffold\tStart\tStop\tStrand\tProt Len\tPfam\tExons\tTMHMM\tHMMTOP\tPhobius\tConsensus\tGene Model\tBest Query\n";


my $pseudo_count=0; my $normal_count=0;

my $pfam_count=0; my $notpfam_count=0;

my $start_chr; my $start_chr_start; my $start_chr_end;

my $cnt=0;
my @out9str_arr;
foreach my $w(sort {$a <=> $b} keys %hash_suppl)
{
	
	if($cnt==0)
	{
		$start_chr=$hash_suppl{$w}{"Scaf"};
		$start_chr_start=$hash_suppl{$w}{"GeneStart"};
		$start_chr_end=$hash_suppl{$w}{"GeneEnd"};
	}
	$cnt++;
	
	my $bedstring=$hash_suppl{$w}{"Scaf"}."\t".$hash_suppl{$w}{"GeneStart"}."\t".$hash_suppl{$w}{"GeneEnd"}."\t".$hash_suppl{$w}{"Scaf"}."_".$hash_suppl{$w}{"GeneStart"}."-".$hash_suppl{$w}{"GeneEnd"}."\t0\t".$hash_suppl{$w}{"QueryDir"}."\n";
	my @cdsboundary=split("-",$hash_suppl{$w}{"CDS"});
	my $out9str= $hash_suppl{$w}{"Name"}."\t"."${$hash{$w}}[0]\t$cdsboundary[0]\t$cdsboundary[-1]\t${$hash{$w}}[9]\t";
	
	$out9str.= length($hash_suppl{$w}{"TargetProt"})."\t";
	$out9str.=(@cdsboundary/2)."\t";

	if($hash_suppl{$w}{"Name"}=~/Normal/)
	{
		$normal_count++;
	}
	else
	{
		$pseudo_count++;
	}
	
	my $cdsstring=$cdsboundary[0]."..".$cdsboundary[1];
		
	$bedstring.=$hash_suppl{$w}{"Scaf"}."\t".$cdsboundary[0]."\t".$cdsboundary[1]."\tcds\t0\t".$hash_suppl{$w}{"QueryDir"}."\n";
	#print "Finding CDS...\n";
	for(my $y=2; $y<@cdsboundary; $y=$y+2)
	{
		$cdsstring.=",".$cdsboundary[$y]."..".$cdsboundary[$y+1];
		$bedstring.=$hash_suppl{$w}{"Scaf"}."\t".$cdsboundary[$y]."\t".$cdsboundary[$y+1]."\tcds\t0\t".$hash_suppl{$w}{"QueryDir"}."\n";
	}

	$out9str.= "$cdsstring\t".$hash_suppl{$w}{"Query"}."\n";
	print OUT9 $out9str;
	push(@out9str_arr,$out9str);
	$hash_suppl{$w}{"CDS"}=$cdsstring;

}
close OUT9; 

#######################################################################################################################################################################################
##### INSECTOR - STEP 8 - Extend boundaries of the previously predicted gene regions (P1) #############################################################################################
##### INSECTOR - STEP 9 - Perform GeneWise gene search for the extended boundaries with the correspondinf identified best query from P1 hits -> Generate P2 hits ######################
#######################################################################################################################################################################################

print "Running GeneWise for hit refinement. It may take a while...\n\n";


my %genomeSeqHash;
open(FH5,"<$filename.modified_seqfile");

my $ke="";
my $valu="";

while(<FH5>)
{
	chomp $_;
	my $lineInFile=$_;
	if($lineInFile=~/^>/)
	{
		if($ke ne "")
		{
			$genomeSeqHash{$ke}=$valu;
		}

		if($lineInFile=~/^(.*?)\s+/)
		{
			$ke=$1;
		}
		else
		{
			$ke=$lineInFile;
		}

		$valu="";
	}
	else
	{
		$valu.=$lineInFile;
	}
}
close FH5;
$genomeSeqHash{$ke}=$valu;

# foreach my $keyys(keys %genomeSeqHash)
# {
# 	if ($genomeSeqHash{$keyys} =~ m/^[ATGCNURYWSMKBHDV]*$/i)
# 	{
# 		#print "$keyys\nFine\n";
#                 $genomeSeqHash{$keyys}=~tr/[Uu]/T/;
#                 $genomeSeqHash{$keyys}=~tr/[RYWSMKBHDVrywsmkbhdv]/N/;
#                 #print "Degenerate characters like R,Y,W,S,M,K,B,H,D,V were replaced with N.\n";
# 	}
# 	else
# 	{
# 		print "$keyys\nError in the sequence. Not a DNA sequence. It contains alphabets other than A,T,G,C,N,U,R,Y,W,S,M,K,B,H,D,V.\n"; die;
# 	}
# }

my %querySeqHash;
open(FH6,"<$queryFile");
#open(FH6,"");
#$/=">";
my $qke="";
my $qvalu="";

while(<FH6>)
{
	chomp $_;
	my $qlineInFile=$_;
	if($qlineInFile=~/^>/)
	{
		if($qke ne "")
		{
			$querySeqHash{$qke}=$qvalu;
		}

		if($qlineInFile=~/^(.*?)\s+/)
		{
			$qke=$1;
		}
		else
		{
			$qke=$qlineInFile;
		}
		$qvalu="";
	}
	else
	{
		$qvalu.=$qlineInFile;
	}
}
close FH6;
$querySeqHash{$qke}=$qvalu;

# foreach my $qkeyys(keys %querySeqHash)
# {
# 	$querySeqHash{$qkeyys}=~s/ //g;
# 	if ($querySeqHash{$qkeyys} =~ m/^[ACDEFGHIKLMNPQRSTVWYXZ\*]*$/i)
# 	{
# 		#print "$qkeyys\nFine\n";
# 		#print "$qkeyys\n";
# 	}
# 	else
# 	{
# 		print "$qkeyys\nError in the sequence. Not a protein sequence. It contains alphabets other than A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,X,Z.\n"; die;
# 	}
# }

#sort table#
`tail -n +2 $filename.OR_table.txt | sort -o $filename.OR_table_sorted.txt -k 2,2 -k 3,3n`;

open(FH7,"<$filename.OR_table_sorted.txt");
#my $genewisepath="$basepath/tools/wise2.4.1/src/bin";
#$genewisepath 1,2,3,4,13
$/="\n";
my %hash_ali_file_OR_table_sorted_txt;
my $counter_alifile=0;
while(<FH7>)
{
	chomp $_;
	@{$hash_ali_file_OR_table_sorted_txt{$counter_alifile}}=split("\t",$_);
	$counter_alifile++;
}

my $size_of_hash= keys %hash_ali_file_OR_table_sorted_txt;
truncate '$filename.combinedGeneWiseRes_7tm6_nonfiltered', 0;

my $genewiserun_count=0;
my $max_dna_u;
my $max_dna_v;
my $first=0;
my $last=0;

#Extending boundaries
foreach my $tmline (sort {$a <=> $b} keys %hash_ali_file_OR_table_sorted_txt)
{
	my @tm7_entries=@{$hash_ali_file_OR_table_sorted_txt{$tmline}};
	my $dna_u; 
	my $dna_v;

	if (int($tm7_entries[5])>=350)
	{
		$dna_u=$tm7_entries[2]-500;
		$dna_v=$tm7_entries[3]+500;
	}	
	elsif (int($tm7_entries[5])>=300)
	{
		$dna_u=$tm7_entries[2]-1000;
		$dna_v=$tm7_entries[3]+1000;
	}
	elsif (int($tm7_entries[5])>=200)
	{
		$dna_u=$tm7_entries[2]-3000;
		$dna_v=$tm7_entries[3]+3000;
	}
	elsif (int($tm7_entries[5])>=100)
	{
		$dna_u=$tm7_entries[2]-4500;
		$dna_v=$tm7_entries[3]+4500;
	}
	elsif (int($tm7_entries[5])<100)
	{
		$dna_u=$tm7_entries[2]-6000;
		$dna_v=$tm7_entries[3]+6000;		
	}



	if((defined $hash_ali_file_OR_table_sorted_txt{$tmline-1}) && (${$hash_ali_file_OR_table_sorted_txt{$tmline-1}}[1] ne ${$hash_ali_file_OR_table_sorted_txt{$tmline}}[1]))
	{
		$first=0;
	}

	if ($first == 0)
	{
		$max_dna_u=$dna_u;
	}
	$first=1;



	if((defined $hash_ali_file_OR_table_sorted_txt{$tmline-1}) && (${$hash_ali_file_OR_table_sorted_txt{$tmline-1}}[1] eq ${$hash_ali_file_OR_table_sorted_txt{$tmline}}[1]))
	{
		
		if(${$hash_ali_file_OR_table_sorted_txt{$tmline-1}}[0]=~/Complete/)
		{
			$max_dna_u=int(${$hash_ali_file_OR_table_sorted_txt{$tmline-1}}[3]) + 500;
		}

		if(($max_dna_u - $dna_u) >= 0)
		{
			$dna_u=$max_dna_u;
		}		
		# elsif(${$hash_ali_file_OR_table_sorted_txt{$tmline-1}}[0]=~/Partial/)
		# {
		# 	$max_dna_u=${$hash_ali_file_OR_table_sorted_txt{$tmline-1}}[2];
		# }
		# else
		# {
		# 	print "Error!";
		# }
	}

	for (my $i=$tmline; $i<=$size_of_hash; $i++)
	{

		if((defined $hash_ali_file_OR_table_sorted_txt{$i+1}) && (${$hash_ali_file_OR_table_sorted_txt{$i+1}}[1] eq ${$hash_ali_file_OR_table_sorted_txt{$i}}[1]))
		{
			if(${$hash_ali_file_OR_table_sorted_txt{$i+1}}[0]=~/Complete/)
			{
				$max_dna_v=int(${$hash_ali_file_OR_table_sorted_txt{$i+1}}[2]) -500;
				$last=1;
				last;
			}		
		}
		else
		{
			$last=0;
			last;
		}

	}

	if ($last == 0)
	{
		$max_dna_v=$dna_v;
	}

	if(($dna_v - $max_dna_v) >= 0)
	{
		$dna_v=$max_dna_v;
	}

	$tm7_entries[8]=">".$tm7_entries[8];

	my @qry=split(/\s+/,"$tm7_entries[8]");
	my $qfil=$qry[0];
	$qfil=~s/[^a-zA-Z0-9]//g;
	#print "$qfil";
	if((! -e "$filename.$qfil.genewise.fasta") || (-z "$filename.$qfil.genewise.fasta"))
	{
		open(OUT19,">$filename.$qfil.genewise.fasta") or die "Couldn't open protein: $!";		
		print "$qry[0]\n";
		print "$querySeqHash{$qry[0]}\n";
		print OUT19 "$qry[0]\n$querySeqHash{$qry[0]}";
		close OUT19;
	}
	#print "z",$genewiserun_count+1,"\n$dna_u\t$dna_v\t$max_dna_u\t$max_dna_v\n";

	$tm7_entries[1]=">".$tm7_entries[1];
	my @trgt=split(/\s+/,"$tm7_entries[1]");
	my $gfil=$trgt[0];
	$gfil=~s/[^a-zA-Z0-9]//g;
	if((! -e "$filename.$gfil.genewise.fasta") || (-z "$filename.$gfil.genewise.fasta"))
	{
		open(OUT20,">$filename.$gfil.genewise.fasta") or die "Couldn't open genome: $!";		
		print OUT20 "$trgt[0]\n$genomeSeqHash{$trgt[0]}";
		close OUT20;
	}


	if($dna_u <= 0)
		{$dna_u=1}

	if($dna_v > length($genomeSeqHash{$trgt[0]}))
	{
		$dna_v = length($genomeSeqHash{$trgt[0]});
	}

	
	#perform GeneWise run
	if($tm7_entries[4] eq '+')
	{
		#`$genewisepath/genewise $qfil.fasta $gfil.fasta -u $dna_u -v $dna_v -tfor -para -pretty -genes -trans -cdna -embl -ace -gff -diana -init local -nosplice_gtag -null syn -alg 623 >>combinedGeneWiseRes_vikas`;
		#`$genewisepath/genewise $tm7_entries[8].fasta $tm7_entries[1].fasta`;
		#`export WISECONFIGDIR=$basepath/tools/wise2.4.1/wisecfg/;$genewisepath/genewise $filename.$qfil.genewise.fasta $filename.$gfil.genewise.fasta -u $dna_u -v $dna_v -tfor -para -pretty -genes -trans -cdna -embl -ace -gff -diana -nosplice_gtag -null syn -alg 623 -alln 0.5 -quiet -silent >> $filename.combinedGeneWiseRes_7tm6_nonfiltered`;
		`export WISECONFIGDIR=$wiseconfigdir_path;$genewisepath/genewise $filename.$qfil.genewise.fasta $filename.$gfil.genewise.fasta -u $dna_u -v $dna_v -tfor -para -pretty -genes -trans -cdna -embl -ace -gff -diana -nosplice_gtag -null syn -alg 623 -alln 0.5 -quiet -silent >> $filename.combinedGeneWiseRes_7tm6_nonfiltered`;
		$genewiserun_count++;
	}
	elsif($tm7_entries[4] eq '-')
	{
		#`$genewisepath/genewise $qfil.fasta $gfil.fasta -u $dna_u -v $dna_v -trev -para -pretty -genes -trans -cdna -embl -ace -gff -diana -init local -nosplice_gtag -null syn -alg 623 >>combinedGeneWiseRes_vikas`;			
		#`$genewisepath/genewise $tm7_entries[8].fasta $tm7_entries[1].fasta`;
		#`export WISECONFIGDIR=$basepath/tools/wise2.4.1/wisecfg/;$genewisepath/genewise $filename.$qfil.genewise.fasta $filename.$gfil.genewise.fasta -u $dna_u -v $dna_v -trev -para -pretty -genes -trans -cdna -embl -ace -gff -diana -nosplice_gtag -null syn -alg 623 -alln 0.5 -quiet -silent >> $filename.combinedGeneWiseRes_7tm6_nonfiltered`;			
		`export WISECONFIGDIR=$wiseconfigdir_path;$genewisepath/genewise $filename.$qfil.genewise.fasta $filename.$gfil.genewise.fasta -u $dna_u -v $dna_v -trev -para -pretty -genes -trans -cdna -embl -ace -gff -diana -nosplice_gtag -null syn -alg 623 -alln 0.5 -quiet -silent >> $filename.combinedGeneWiseRes_7tm6_nonfiltered`;			
		$genewiserun_count++;
	}
	else
	{
		print "vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvError";
	}
	
}

close FH7;
`rm $filename.*.genewise.fasta`;
print "Genewise run complete : $genewiserun_count\n";

open(FH8,"$filename.combinedGeneWiseRes_7tm6_nonfiltered") or die "Couldn't open genewise output: $!";
$/='genewise $Name: wise2-4-1 $ (unreleased release)';

my %genewise_hash;

open(OUT20,">$filename.genewise.table");
open(OUT22,">$filename.genewise_protein_seq.txt");

while(<FH8>)
{
	my $gene=$_;
	chomp $gene;
	if($gene=~/Query protein:\s+(.*?)\n.*?Target Sequence\s+(.*?)\nStrand:\s+(.*?)\n.*?\/\/\n.*?\/\/\n.*?\/\/\n.*?\/\/\n>(.*?.sp.tr)\n(.*?)\/\/\n>(.*?.sp)\n(.*?)\/\/\n.*?\/\/\n(.*?)\/\//ms)
	{
		my $quer=$1;
		my $target=$2;
		my $direct=$3;
		my $protHead=$4;
		my $protSeq=$5;
		my $nucHead=$6;
		my $nucSeq=$7;
		my $gffpart=$8;
		$protSeq=~s/\n//g;
		$nucSeq=~s/\n//g;
		my @gfftemp=($gffpart=~/cds\s+(.*?)\s+(.*?)\s+.*?\n/gms);
		my $firstExonLength=abs($gfftemp[1]-$gfftemp[0])+1;
		my $lastExonLength=abs($gfftemp[-1]-$gfftemp[-2])+1;
		my @phasetemp=($gffpart=~/cds\s+.*?\s+.*?\s+.*?\n/gms);

		if (scalar(@gfftemp)>2)
		{		
			if($firstExonLength < 45 && abs($gfftemp[2]-$gfftemp[1]) > 1000 ) 
			{
				#print "First exon is too small.\t@gfftemp\t".length($protSeq)."\t".length($nucSeq)."\t";
				shift(@gfftemp); shift(@gfftemp);
				if($firstExonLength%3==0)
				{
					$protSeq=substr($protSeq,int($firstExonLength/3));
					$nucSeq=substr($nucSeq,$firstExonLength);
				}
				elsif($firstExonLength%3==1)
				{
					$protSeq=substr($protSeq,int($firstExonLength/3)+1);
					$nucSeq=substr($nucSeq,$firstExonLength+2);
					$gfftemp[0]=$gfftemp[0]+2;			
				}
				else
				{
					$protSeq=substr($protSeq,int($firstExonLength/3)+1);
					$nucSeq=substr($nucSeq,$firstExonLength+1);
					$gfftemp[0]=$gfftemp[0]+1;		
				}

			}

			if($lastExonLength < 15 && abs($gfftemp[-2]-$gfftemp[-3]) > 1000 )
			{
				#print "Last exon is too small.\t@gfftemp\t".length($protSeq)."\t".length($nucSeq)."\t";
				pop(@gfftemp);pop(@gfftemp);
				if($lastExonLength%3==0)
				{
					$protSeq=substr($protSeq,0,length($protSeq)-int($lastExonLength/3));
					$nucSeq=substr($nucSeq,0,length($nucSeq)-$lastExonLength);			
				}
				elsif($lastExonLength%3==1)
				{
					$protSeq=substr($protSeq,0,length($protSeq)-int($lastExonLength/3)-1);
					$nucSeq=substr($nucSeq,0,(length($nucSeq)-$lastExonLength)-2);
					$gfftemp[-1]=$gfftemp[-1]-2;	
				}
				else
				{
					$protSeq=substr($protSeq,0,length($protSeq)-int($lastExonLength/3)-1);
					$nucSeq=substr($nucSeq,0,(length($nucSeq)-$lastExonLength)-1);
					$gfftemp[-1]=$gfftemp[-1]-1;	
				}
			}
		}
		
		$genewise_hash{$target}{$gfftemp[0]}{$gfftemp[-1]}{"protHead"}=$protHead;
		$genewise_hash{$target}{$gfftemp[0]}{$gfftemp[-1]}{"protSeq"}=$protSeq;
		$genewise_hash{$target}{$gfftemp[0]}{$gfftemp[-1]}{"nucSeq"}=$nucSeq;
		$genewise_hash{$target}{$gfftemp[0]}{$gfftemp[-1]}{"nucHead"}=$nucHead;
		$genewise_hash{$target}{$gfftemp[0]}{$gfftemp[-1]}{"strand"}=$direct;
		@{$genewise_hash{$target}{$gfftemp[0]}{$gfftemp[-1]}{"cds"}}=@gfftemp;
		my $protLength=length($protSeq);
		$protSeq=~s/\>.*?\.sp\.tr/\*/msg;

		if($direct eq "forward")
		{
			print OUT22"\>".$target."_".$gfftemp[0]."-".$gfftemp[-1]."_OR_";	
		}
		else
		{
			print OUT22"\>".$target."_".$gfftemp[-1]."-".$gfftemp[0]."_OR_";
		}
		#print $target."\t"."#@{$genewise_hash{$target}{$gfftemp[0]}{$gfftemp[-1]}{'cds'}}#\t".length($protSeq)."\t"."$protSeq\t".length($nucSeq)."\t$nucSeq\n";
		print OUT20"$target\t$gfftemp[0]\t$gfftemp[-1]\t$direct\t$protLength\t";
		if($protLength<$lengthCutoff)
		{
			print OUT20"Partial\t";
			print OUT22"Partial_";
		}
		else
		{
			print OUT20"Complete\t";
			print OUT22"Complete_";
		}

		if($protSeq=~/\*/)
		{
			print OUT20"Pseudo\t";
			print OUT22"Pseudo_";
		}
		else
		{
			print OUT20"Normal\t";
			print OUT22"Normal_";
		}

		my $flg=findMet($protSeq);
		print OUT22 "$flg\n$protSeq\n";

		print OUT20"@gfftemp\t$quer\n";
	}

}
close FH8;
close OUT22;
close OUT20;

#######################################################################################################################################################################################
##### INSECTOR - STEP 10 - Compare P1 and P2 predictions and retain the best unique hits ##############################################################################################
#######################################################################################################################################################################################

open(FH9,"$filename.OR_table_sorted.txt");
$/="\n";
my @insectOR_arr;
my $count_tbl=0;
while(<FH9>)
{
	chomp $_;
	my @insectORline=split("\t",$_);
	my @insectOR_arr_zero=split("_",$insectORline[0]);
	$insectORline[4]=~s/\+/forward/;
	$insectORline[4]=~s/\-/reverse/;
	$insectORline[7]=~s/\.\./ /g;
	$insectORline[7]=~s/\,/ /g;
	

	my @genemodel=split(/\s/,$insectORline[7]);
	
	$insectORline[7]=join(" ",@genemodel);
	#@{$insectOR_arr[$count_tbl]}=($insectORline[1],$insectORline[2],$insectORline[3],$insectORline[4],$insectORline[5],$insectOR_arr_zero[4],$insectOR_arr_zero[5],$insectORline[7],$insectORline[8]);
	@{$insectOR_arr[$count_tbl]}=($insectORline[1],$insectORline[2],$insectORline[3],$insectORline[4],$insectORline[5],$insectOR_arr_zero[-3],$insectOR_arr_zero[-2],$insectORline[7],$insectORline[8]);
	$count_tbl++;
}
close FH9;

##########################################   making hash of insectOR and Genewise proteins.......################################

my %insector_proteins_hash;
open(FH11, "$filename.ORs.pep");
$/="\n";
my $insector_proteins_key="";
my $insector_proteins_value="";

while(<FH11>)
{
	chomp $_;
	my $line=$_;

	if($line=~/^>/)
	{
		if($insector_proteins_key ne "")
		{
			$insector_proteins_hash{$insector_proteins_key}=$insector_proteins_value;
		}
		$insector_proteins_key=$line;
		$insector_proteins_value="";
	}
	else
	{
		$insector_proteins_value.=$line;
	}
}	

close FH11;
$insector_proteins_hash{$insector_proteins_key}=$insector_proteins_value;

my %genewise_proteins_hash;
#open(FH12, "genewise_7tm_6_protein_seq.txt");
open(FH12, "$filename.genewise_protein_seq.txt");
$/="\n";
my $genewise_proteins_key="";
my $genewise_proteins_value="";

while(<FH12>)
{
	
	chomp $_;
	my $line=$_;

	if($line=~/^>/)
	{
		if($genewise_proteins_key ne "")
		{
			$genewise_proteins_hash{$genewise_proteins_key}=$genewise_proteins_value;

		}
		$genewise_proteins_key=$line;
		$genewise_proteins_value="";
		
	}
	else
	{
		$genewise_proteins_value.=$line;
	}
}	

$genewise_proteins_hash{$genewise_proteins_key}=$genewise_proteins_value;
my $num=0;
foreach my $ke(keys %genewise_proteins_hash)
{
	$num++;
}

close FH12;
my $size_of_genewisehash=keys %genewise_proteins_hash;
my $size_of_insectorhash=keys %insector_proteins_hash;

open(FH10, "$filename.genewise.table");
$/="\n";

my @genewise_array=<FH10>;
close FH10;
$count_tbl=0;
my $count_tbl1=0;
my @bestgenes;

my %consolidatedproteins_hash;

truncate '$filename.mergedTable_1_withoverlaps_7tm6notfiltered.txt', 0;
truncate '$filename.consolidatedproteins_7tm6notfiltered.txt', 0;
truncate '$filename.consolidatedTable_1_7tm6notfiltered.txt', 0;
truncate '$filename.latest_consolidatedTable_1_7tm6notfiltered.txt', 0;
truncate '$filename.latest_consolidatedproteins_7tm6notfiltered.txt', 0;


open(OUT23,">$filename.mergedTable_1_withoverlaps_7tm6notfiltered.txt");

foreach(@genewise_array)
{
	my @genewiseline=split("\t",$genewise_array[$count_tbl]);
	chomp $genewiseline[8];

	my $temp_loc;
	if($genewiseline[3] eq "reverse")
	{
		$temp_loc=$genewiseline[1];
		$genewiseline[1]=$genewiseline[2];
		$genewiseline[2]=$temp_loc;
	}
	
	my $flag=1;
	if( ($insectOR_arr[$count_tbl][1] >= $genewiseline[1] && $insectOR_arr[$count_tbl][1] <= $genewiseline[2] ) || ($insectOR_arr[$count_tbl][2] >= $genewiseline[1] && $insectOR_arr[$count_tbl][2] <= $genewiseline[2]) ||
		($genewiseline[1] >= $insectOR_arr[$count_tbl][1] && $genewiseline[1] <= $insectOR_arr[$count_tbl][2]) || ($genewiseline[2] >= $insectOR_arr[$count_tbl][1] && $genewiseline[2] <= $insectOR_arr[$count_tbl][2]) )
	{
		$flag=0;
	}
	my $a="";
	my $b="";
	if(${$insectOR_arr[$count_tbl]}[5] eq "Partial" && $genewiseline[5] eq "Complete")
	{
		push(@{$bestgenes[$count_tbl1]},@genewiseline);
		foreach my $key_in_genewise_proteins (keys %genewise_proteins_hash)
		{
			$b=$genewiseline[0]."_";
			if (($key_in_genewise_proteins=~m/$b/) && (($key_in_genewise_proteins=~m/$genewiseline[1]/) || ($key_in_genewise_proteins=~m/$genewiseline[2]/)))
			{
				$consolidatedproteins_hash{$key_in_genewise_proteins}=$genewise_proteins_hash{$key_in_genewise_proteins};
				#print OUT23"$key_in_genewise_proteins\n$genewise_proteins_hash{$key_in_genewise_proteins}\n";
				foreach my $key_in_insector_proteins (keys %insector_proteins_hash)
				{
					$a=${$insectOR_arr[$count_tbl]}[0]."_";
					if (($key_in_insector_proteins=~m/$a/) && ($key_in_insector_proteins=~m/${$insectOR_arr[$count_tbl]}[1]/))	
					{
						if ($flag==1)
						{
							$count_tbl1++;
							$consolidatedproteins_hash{$key_in_insector_proteins}=$insector_proteins_hash{$key_in_insector_proteins};
							push(@{$bestgenes[$count_tbl1]},@{$insectOR_arr[$count_tbl]});
						}
					}
				}		
			}
		}
	}
	elsif(${$insectOR_arr[$count_tbl]}[5] eq "Complete" && $genewiseline[5] eq "Partial")
	{
		push(@{$bestgenes[$count_tbl1]},@{$insectOR_arr[$count_tbl]});
		foreach my $key_in_insector_proteins (keys %insector_proteins_hash)
		{
			$a=${$insectOR_arr[$count_tbl]}[0]."_";
			if (($key_in_insector_proteins=~m/$a/) && ($key_in_insector_proteins=~m/${$insectOR_arr[$count_tbl]}[1]/))
			{
				$consolidatedproteins_hash{$key_in_insector_proteins}=$insector_proteins_hash{$key_in_insector_proteins};
				#print OUT23"$key_in_insector_proteins\n$insector_proteins_hash{$key_in_insector_proteins}\n";
				foreach my $key_in_genewise_proteins (keys %genewise_proteins_hash)
				{
					$b=$genewiseline[0]."_";
					if (($key_in_genewise_proteins=~m/$b/) && (($key_in_genewise_proteins=~m/$genewiseline[1]/) || ($key_in_genewise_proteins=~m/$genewiseline[2]/)))
					{
						if ($flag==1)
						{
							$count_tbl1++;
							$consolidatedproteins_hash{$key_in_genewise_proteins}=$genewise_proteins_hash{$key_in_genewise_proteins};
							push(@{$bestgenes[$count_tbl1]},@genewiseline);
							#print "h222222222222\n";
						}
					}
				}		
			}
		}		
	}
	else
	{
		if(${$insectOR_arr[$count_tbl]}[6] eq "Normal" && $genewiseline[6] eq "Pseudo")
		{
			push(@{$bestgenes[$count_tbl1]},@{$insectOR_arr[$count_tbl]});
			foreach my $key_in_insector_proteins (keys %insector_proteins_hash)
			{
				$a=${$insectOR_arr[$count_tbl]}[0]."_";
				if (($key_in_insector_proteins=~m/$a/) && ($key_in_insector_proteins=~m/${$insectOR_arr[$count_tbl]}[1]/))
				{
					$consolidatedproteins_hash{$key_in_insector_proteins}=$insector_proteins_hash{$key_in_insector_proteins};
					#print OUT23"$key_in_insector_proteins\n$insector_proteins_hash{$key_in_insector_proteins}\n";
					foreach my $key_in_genewise_proteins (keys %genewise_proteins_hash)
					{
						$b=$genewiseline[0]."_";
						if (($key_in_genewise_proteins=~m/$b/) && (($key_in_genewise_proteins=~m/$genewiseline[1]/) || ($key_in_genewise_proteins=~m/$genewiseline[2]/)))
						{
							if ($flag==1)
							{
								$count_tbl1++;
								$consolidatedproteins_hash{$key_in_genewise_proteins}=$genewise_proteins_hash{$key_in_genewise_proteins};
								push(@{$bestgenes[$count_tbl1]},@genewiseline);
								#print "h222222222222\n";
							}
						}
					}		
				}
			}					
		}
		elsif(${$insectOR_arr[$count_tbl]}[6] eq "Pseudo" && $genewiseline[6] eq "Normal")
		{
			push(@{$bestgenes[$count_tbl1]},@genewiseline);
			foreach my $key_in_genewise_proteins (keys %genewise_proteins_hash)
			{
				$b=$genewiseline[0]."_";
				if (($key_in_genewise_proteins=~m/$b/) && (($key_in_genewise_proteins=~m/$genewiseline[1]/) || ($key_in_genewise_proteins=~m/$genewiseline[2]/)))
				{
					$consolidatedproteins_hash{$key_in_genewise_proteins}=$genewise_proteins_hash{$key_in_genewise_proteins};					
					#print OUT23"$key_in_genewise_proteins\n$genewise_proteins_hash{$key_in_genewise_proteins}\n";
					foreach my $key_in_insector_proteins (keys %insector_proteins_hash)
					{
						$a=${$insectOR_arr[$count_tbl]}[0]."_";
						if (($key_in_insector_proteins=~m/$a/) && ($key_in_insector_proteins=~m/${$insectOR_arr[$count_tbl]}[1]/))	
						{
							if ($flag==1)
							{
								$count_tbl1++;
								$consolidatedproteins_hash{$key_in_insector_proteins}=$insector_proteins_hash{$key_in_insector_proteins};
								push(@{$bestgenes[$count_tbl1]},@{$insectOR_arr[$count_tbl]});
								#print "h111111111\n";
								#print "9999999999999999999", $key_in_insector_proteins, $insector_proteins_hash{$key_in_insector_proteins} ,"\n";
							}
						}
					}							
				}
			}			
		}
		else
		{
			if(${$insectOR_arr[$count_tbl]}[4] < $genewiseline[4])
			{
				push(@{$bestgenes[$count_tbl1]},@genewiseline);
				foreach my $key_in_genewise_proteins (keys %genewise_proteins_hash)
				{
					#print "tttttttttttttttttt".$key_in_genewise_proteins;					
					$b=$genewiseline[0]."_";
					if (($key_in_genewise_proteins=~/$b/) && (($key_in_genewise_proteins=~m/$genewiseline[1]/) || ($key_in_genewise_proteins=~m/$genewiseline[2]/)))
					{
						$consolidatedproteins_hash{$key_in_genewise_proteins}=$genewise_proteins_hash{$key_in_genewise_proteins};						
						#print OUT23"$key_in_genewise_proteins\n$genewise_proteins_hash{$key_in_genewise_proteins}\n";
						foreach my $key_in_insector_proteins (keys %insector_proteins_hash)
						{
							$a=${$insectOR_arr[$count_tbl]}[0]."_";
							if (($key_in_insector_proteins=~m/$a/) && ($key_in_insector_proteins=~m/${$insectOR_arr[$count_tbl]}[1]/))	
							{
								if ($flag==1)
								{
									$count_tbl1++;
									$consolidatedproteins_hash{$key_in_insector_proteins}=$insector_proteins_hash{$key_in_insector_proteins};
									push(@{$bestgenes[$count_tbl1]},@{$insectOR_arr[$count_tbl]});
									#print "h111111111\n";
								}
							}
						}								
					}
				}
			}
			elsif(${$insectOR_arr[$count_tbl]}[4] > $genewiseline[4])
			{
				push(@{$bestgenes[$count_tbl1]},@{$insectOR_arr[$count_tbl]});
				foreach my $key_in_insector_proteins (keys %insector_proteins_hash)
				{
					$a=${$insectOR_arr[$count_tbl]}[0]."_";
					if (($key_in_insector_proteins=~m/$a/) && ($key_in_insector_proteins=~m/${$insectOR_arr[$count_tbl]}[1]/))
					{
						$consolidatedproteins_hash{$key_in_insector_proteins}=$insector_proteins_hash{$key_in_insector_proteins};
						#print OUT23"$key_in_insector_proteins\n$insector_proteins_hash{$key_in_insector_proteins}\n";
						foreach my $key_in_genewise_proteins (keys %genewise_proteins_hash)
						{
							$b=$genewiseline[0]."_";
							if (($key_in_genewise_proteins=~m/$b/) && (($key_in_genewise_proteins=~m/$genewiseline[1]/) || ($key_in_genewise_proteins=~m/$genewiseline[2]/)))
							{
								if ($flag==1)
								{
									$count_tbl1++;
									$consolidatedproteins_hash{$key_in_genewise_proteins}=$genewise_proteins_hash{$key_in_genewise_proteins};
									push(@{$bestgenes[$count_tbl1]},@genewiseline);
									#print "h222222222222\n";
								}
							}
						}								
					}
				}						
			}
			else
			{
				push(@{$bestgenes[$count_tbl1]},@{$insectOR_arr[$count_tbl]});
				foreach my $key_in_insector_proteins (keys %insector_proteins_hash)
				{
					$a=${$insectOR_arr[$count_tbl]}[0]."_";
					if (($key_in_insector_proteins=~m/$a/) && ($key_in_insector_proteins=~m/${$insectOR_arr[$count_tbl]}[1]/))
					{
						$consolidatedproteins_hash{$key_in_insector_proteins}=$insector_proteins_hash{$key_in_insector_proteins};
						#print OUT23"$key_in_insector_proteins\n$insector_proteins_hash{$key_in_insector_proteins}\n";

						foreach my $key_in_genewise_proteins (keys %genewise_proteins_hash)
						{
							$b=$genewiseline[0]."_";
							if (($key_in_genewise_proteins=~m/$b/) && (($key_in_genewise_proteins=~m/$genewiseline[1]/) || ($key_in_genewise_proteins=~m/$genewiseline[2]/)))
							{
								if ($flag==1)
								{
									$count_tbl1++;
									$consolidatedproteins_hash{$key_in_genewise_proteins}=$genewise_proteins_hash{$key_in_genewise_proteins};
									push(@{$bestgenes[$count_tbl1]},@genewiseline);
									#print "h222222222222\n";
								}
							}
						}								
					}
				}					
			}
		}		
	}

	$count_tbl++;
	$count_tbl1++;
}

@bestgenes=  sort { $a->[0] cmp $b->[0]  ||  $a->[1] <=> $b->[1] ||  $a->[2] <=> $b->[2] } @bestgenes;

foreach my $ke1 (0..(scalar(@bestgenes)-1))
{
	foreach my $ke2 (0..8)
	{
		print OUT23 "${$bestgenes[$ke1]}[$ke2]\t";
	}
	print OUT23"\n";
}

close OUT23;
my @ids;
my $newid=0;
push(@{$ids[$newid]},0);
my $large=0;
for ($count_tbl=1; $count_tbl < @bestgenes; $count_tbl++)
{
	my $i=${$bestgenes[$count_tbl]}[1];
	my $j=${$bestgenes[$count_tbl-1]}[2];
	if (${$bestgenes[$count_tbl]}[0] eq ${$bestgenes[$count_tbl-1]}[0])
	{
		if ($j>=$large)
		{
		$large=$j;

		}
	}
	else {$large=0;}		
	if(!((${$bestgenes[$count_tbl]}[0] eq ${$bestgenes[$count_tbl-1]}[0]) && ($i < $large)))
	{	
		$newid++;
	}
	push(@{$ids[$newid]},$count_tbl);
}

open(OUT21,">$filename.consolidatedTable_1_7tm6notfiltered.txt");
open(OUT22,">$filename.consolidatedproteins_7tm6notfiltered.txt");
my $table=0;
my $seq=0;
foreach (@ids)
{
	my @tmp_arr=@{$_};
	my %tmp_hash;
	tie %tmp_hash, "Tie::IxHash";
	foreach my $tmp_arr_ele (@tmp_arr)
	{
		$tmp_hash{$tmp_arr_ele}=${$bestgenes[$tmp_arr_ele]}[4];
	}
	
	my $counter=0;
	my $length=0;
	my @selected_keys;
	my $sel_flag=0;
	my $prev_start=0;
	my $prev_end=0;
	foreach my $tmpkey (sort { $tmp_hash{$b} <=> $tmp_hash{$a} } keys %tmp_hash)
	{

		$counter++;
		if ($counter==1)
		{
			$length=$tmp_hash{$tmpkey};
			$prev_start=$bestgenes[$tmpkey][1];
			$prev_end= $bestgenes[$tmpkey][2];	
			#print $bestgenes[$tmpkey][1], "\t", $bestgenes[$tmpkey][2], "\t", $length, "\n";		
			push(@selected_keys, $tmpkey);
			$sel_flag=1;
			next;		
		}

		my $flg=1;
		my $sel_key_1= $selected_keys[-1];
		#####################################   Compare Overlap with last added element in the array    #######################
		# foreach my $sel_key(@selected_keys) 
		# {
		if( ($bestgenes[$tmpkey][1] >= $bestgenes[$sel_key_1][1] && $bestgenes[$tmpkey][1] <= $bestgenes[$sel_key_1][2]) || ($bestgenes[$tmpkey][2] >= $bestgenes[$sel_key_1][1] && $bestgenes[$tmpkey][2] <= $bestgenes[$sel_key_1][2]) ||
			($bestgenes[$sel_key_1][1] >= $bestgenes[$tmpkey][1] && $bestgenes[$sel_key_1][1] <= $bestgenes[$tmpkey][2]) || ($bestgenes[$sel_key_1][2] >= $bestgenes[$tmpkey][1] && $bestgenes[$sel_key_1][2] <= $bestgenes[$tmpkey][2]) )
		{
			$flg=0;
		}
		# 	else
		# 	{}
		# }

		if (($tmp_hash{$tmpkey}==$length) && (($sel_flag==1) && ($flg==0)))
		{
			if (($bestgenes[$tmpkey][2]-$bestgenes[$tmpkey][1]) < ($prev_end-$prev_start))
			{
				pop @selected_keys;
				$sel_flag=0;

			}
		}
		$sel_flag=0;
		$flg=1;
		foreach my $sel_key(@selected_keys)
		{
			if( ($bestgenes[$tmpkey][1] >= $bestgenes[$sel_key][1] && $bestgenes[$tmpkey][1] <= $bestgenes[$sel_key][2]) || ($bestgenes[$tmpkey][2] >= $bestgenes[$sel_key][1] && $bestgenes[$tmpkey][2] <= $bestgenes[$sel_key][2]) ||
				($bestgenes[$sel_key][1] >= $bestgenes[$tmpkey][1] && $bestgenes[$sel_key][1] <= $bestgenes[$tmpkey][2]) || ($bestgenes[$sel_key][2] >= $bestgenes[$tmpkey][1] && $bestgenes[$sel_key][2] <= $bestgenes[$tmpkey][2]) )
			{
				$flg=0;
			}
			else
			{}
		}


		$length=$tmp_hash{$tmpkey};
		$prev_start=$bestgenes[$tmpkey][1];
		$prev_end= $bestgenes[$tmpkey][2];
		#print $bestgenes[$tmpkey][1], "\t", $bestgenes[$tmpkey][2], "\t", $length, "\n";

		if($flg==1)
		{
			push(@selected_keys, $tmpkey);
			$sel_flag=1;
			#print join(", ", @selected_keys), "\thiiiiiii\n";
			next;
		}
		

		#@{$bestgenes[$tmpkey]}[7]=~s/\D/,/g;
		# $counter++;
		# if($counter==1)
		# {
		# 	$selected_keys{$tmpkey}=$tmp_hash{$tmpkey};
		# 	next;
		# }

		# foreach my $tmpkey2 (sort { $selected_keys{$b} <=> $selected_keys{$a} } keys %selected_keys)
		# {
		# 	if( ${$bestgenes[$tmpkey]}[1] )
		# }
	#last;
	}

	foreach my $tmpkey(@selected_keys)
	{
		my @exon_struture=();	
		push(@exon_struture, sort { $a <=> $b } (split(" ", @{$bestgenes[$tmpkey]}[7])));

		print OUT21"$bestgenes[$tmpkey][0]\t$bestgenes[$tmpkey][1]\t$bestgenes[$tmpkey][2]\t$bestgenes[$tmpkey][3]\t$bestgenes[$tmpkey][4]\t$bestgenes[$tmpkey][5]\t$bestgenes[$tmpkey][6]\t@exon_struture\t$bestgenes[$tmpkey][8]\n";
		$table++;
		my $try=0;
		foreach my $ke(keys %consolidatedproteins_hash)
		{
			$try=0;
			my $x=${$bestgenes[$tmpkey]}[0]."_".${$bestgenes[$tmpkey]}[1]."-".${$bestgenes[$tmpkey]}[2];
			if($ke=~/$x/)
			{
				print OUT22 "$ke\n$consolidatedproteins_hash{$ke}\n";
				#$selected_keys{$tmpkey}=$tmp_hash{$tmpkey};
				$seq++;
				$try=1;
				last;
			}
		}
		if ($try==0)
		{
			#print ${bestgenes[$tmpkey]}[1];
		}
	}

}
#print "hiiiiiiiii\t", $table, "\t  jjjjjj", $seq, "\n";

close OUT21;
close OUT22;

`cat $filename.consolidatedTable_1_7tm6notfiltered.txt | sort -o $filename.consolidatedTable_1_7tm6notfiltered.txt -k 1,1 -k 2,2n`;

open(FH13, "$filename.consolidatedTable_1_7tm6notfiltered.txt");
my @final_genewise_array=<FH13>;
close FH13;

my $final_count_tbl=0;
my @remainingbestgenes;
my %remaining_proteins_hash;
my $c=0;
foreach (@insectOR_arr)
{

	my $final_flag=1;
	my $final_count_tbl1=0;
	foreach(@final_genewise_array)
	{
		#chomp $_;
		#my @insectORline=split("\t",$_);
		my @final_genewiseline=split("\t",$final_genewise_array[$final_count_tbl1]);
		chomp $final_genewiseline[8];
		#print "ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc", $insectOR_arr[$final_count_tbl][1], $final_genewiseline[0],"\n";
		if ($final_genewiseline[0] eq $insectOR_arr[$final_count_tbl][0])
		{
			if( ($insectOR_arr[$final_count_tbl][1] >= $final_genewiseline[1] && $insectOR_arr[$final_count_tbl][1] <= $final_genewiseline[2] ) || ($insectOR_arr[$final_count_tbl][2] >= $final_genewiseline[1] && $insectOR_arr[$final_count_tbl][2] <= $final_genewiseline[2]) ||
				($final_genewiseline[1] >= $insectOR_arr[$final_count_tbl][1] && $final_genewiseline[1] <= $insectOR_arr[$final_count_tbl][2]) || ($final_genewiseline[2] >= $insectOR_arr[$final_count_tbl][1] && $final_genewiseline[2] <= $insectOR_arr[$final_count_tbl][2]) )
			{
				$final_flag=0;
			}
		}

		$final_count_tbl1++;

	}

	if ($final_flag==1)
	{
		push(@remainingbestgenes, $insectOR_arr[$final_count_tbl]);
		foreach my $final_key_in_insector_proteins (keys %insector_proteins_hash)
		{
			if (($final_key_in_insector_proteins=~m/${$insectOR_arr[$final_count_tbl]}[0]/) && ($final_key_in_insector_proteins=~m/${$insectOR_arr[$final_count_tbl]}[1]/))
			{
				$remaining_proteins_hash{$final_key_in_insector_proteins}=$insector_proteins_hash{$final_key_in_insector_proteins};
			}	
		}	

		$c++;
	}
	$final_count_tbl++;
}
#print "remainin hits=", $c, "\n";

open(FH14,"$filename.consolidatedTable_1_7tm6notfiltered.txt");


my $number=0;
my @last_final_genewise_array;
while (<FH14>)
{
	chomp $_;
	my @last_final_genewise_array_line=split("\t", $_);
	@{$last_final_genewise_array[$number]}=($last_final_genewise_array_line[0], $last_final_genewise_array_line[1], $last_final_genewise_array_line[2], $last_final_genewise_array_line[3], $last_final_genewise_array_line[4], $last_final_genewise_array_line[5], $last_final_genewise_array_line[6], $last_final_genewise_array_line[7], $last_final_genewise_array_line[8]);
	$number++;
}



# my @last_final_genewise_array=<FH14>;

close FH14;
#print "checkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk", scalar(@remainingbestgenes);
# foreach my $ke3 (0..(scalar(@remainingbestgenes)-1))
# {
# 	foreach my $ke4 (0..8)
# 	{
# 		print "${$remainingbestgenes[$ke3]}[$ke4]\t";
# 	}
# 	print "\n";
# }


open(OUT24,">$filename.latest_consolidatedTable_1_7tm6notfiltered.txt");


my @latest_consolidatedTable_array=(@last_final_genewise_array, @remainingbestgenes);

@latest_consolidatedTable_array=  sort { $a->[0] cmp $b->[0]  ||  $a->[1] <=> $b->[1] ||  $a->[2] <=> $b->[2] } @latest_consolidatedTable_array;

#print "saksasaanccccmccnmcnmcnkjcHAkcncncjnac\n", scalar(@latest_consolidatedTable_array), "\n", ${$latest_consolidatedTable_array[0]}[8];

foreach my $ke5 (0..(scalar(@latest_consolidatedTable_array)-1))
{

	foreach my $ke6 (0..8)
	{
		if ($ke6<=7)
		{
			print OUT24"${$latest_consolidatedTable_array[$ke5]}[$ke6]\t";
		}
		else
		{
			print OUT24"${$latest_consolidatedTable_array[$ke5]}[$ke6]";
		}
	}	

	print OUT24"\n";
}

close OUT24;


my %final_consolidatedproteins_hash;
open(FH15, "$filename.consolidatedproteins_7tm6notfiltered.txt");
$/="\n";
my $consolidated_proteins_key="";
my $consolidated_proteins_value="";

while(<FH15>)
{
	
	chomp $_;
	my $line=$_;

	if($line=~/^>/)
	{
		if($consolidated_proteins_key ne "")
		{
			$final_consolidatedproteins_hash{$consolidated_proteins_key}=$consolidated_proteins_value;

		}
		$consolidated_proteins_key=$line;
		$consolidated_proteins_value="";
		
	}
	else
	{
		$consolidated_proteins_value.=$line;
	}
}	

$final_consolidatedproteins_hash{$consolidated_proteins_key}=$consolidated_proteins_value;

close FH15;

open(OUT25,">$filename.latest_consolidatedproteins_7tm6notfiltered.txt");

my %latest_consolidated_protein_hash;
%latest_consolidated_protein_hash=(%final_consolidatedproteins_hash, %remaining_proteins_hash);
foreach my $t (sort { $a cmp $b  ||  $a <=> $b } keys %latest_consolidated_protein_hash)
{
	print OUT25 "$t\n$latest_consolidated_protein_hash{$t}\n";
}

close OUT25;

#######################################################################################################################################################################################
##### INSECTOR - STEP 11 - Modify the ends for start and stop codons ##################################################################################################################
#######################################################################################################################################################################################


`$basepath/bin/new_finalprocessing.py $filename.modified_seqfile $filename.latest_consolidatedproteins_7tm6notfiltered.txt $filename.latest_consolidatedTable_1_7tm6notfiltered.txt $filename $lengthCutoff`;


#######################################################################################################################################################################################
##### INSECTOR - STEP 12 - Perform validations and other analyses on insectOR predicted proteins ######################################################################################
#######################################################################################################################################################################################


print "Putative protein sequences for pseudogenes contain '*' in place of pseudogenizing elements\nlike stop codons and frame-shifts. These were removed for further analysis.\n\n";
#`sed s/\*//g $filename.final_proteins.txt > $filename.ORs.starRemoved.pep`;
`wc -l $filename.ORs.starRemoved.pep`;


my $tmhmmpred_ref; my $hmmtoppred_ref; my $phobiuspred_ref; my %tmhmmpred_TMHs_hash; my %hmmtoppred_TMHs_hash; my %phobiuspred_TMHs_hash;
my $tmhmmpred_ref2; my $hmmtoppred_ref2; my $phobiuspred_ref2; my $consensuspred_ref; my $consensusscore_ref; my $consensusTMHs_ref; my $seq_ref;
my %tmhmmpred_str_hash; my %hmmtoppred_str_hash; my %phobiuspred_str_hash; my %consensuspred_str_hash; my %consensusscore_str_hash; my %consensus_TMHs_hash; my %seqhash;
my %hmmsearch_hash;

my $hmmsearch_selected=0;
if((!defined $hmmsearch) || (!defined $hmmsearchpath) || (! -e "$hmmsearchpath\/hmmsearch") || (-z "$hmmsearchpath\/hmmsearch") || (! -x "$hmmsearchpath\/hmmsearch"))
{
	print "Hmmsearch option not chosen OR file path not provided OR file does not exist OR it is empty OR is not an executable. Skipping the search for 7tm_6 domain in the protein sequences predicted by this program. This will not affect the core analysis.\n\n";
}
else
{
	$hmmsearch_selected=1;
	print "Searching for 7tm_6 domain...\n\n";
	my $pid = fork();
	if (!defined $pid)
	{
		warn "Fork failed: $!\n";
	}
	elsif($pid==0)
	{
	    `$hmmsearchpath/hmmsearch --domtblout $filename.ORs.starRemoved.pep.hmmsearchout $basepath/hmm/7tm_6.hmm $filename.ORs.starRemoved.pep`;
	    print "7tm_6 domain search finished.\n\n";
	    POSIX::_exit(0);
	}
	else
	{}
}


#my $mast_selected=0;
if((!defined $mast) || (!defined $mastpath) || (! -e "$mastpath\/mast") || (-z "$mastpath\/mast") || (! -x "$mastpath\/mast") || (!defined $motifsfile) || (! -e "$motifsfile") || (-z "$motifsfile"))
{
	print "For Mast option/motifs file -  not chosen OR file path not provided OR file does not exist OR it is empty OR Mast is not an executable. Skipping the search for given motifs using Mast in the protein sequences predicted by this program. This will not affect the core analysis.\n\n";
}
else
{
	#$mast_selected=1;
	print "Searching for motifs...\n\n";
	my $pid = fork();	
	if (!defined $pid)
	{
		warn "Fork failed: $!\n";
	}
	elsif($pid==0)
	{
		#print "$mastpath/mast $motifsfile $filename.ORs.starRemoved.pep -oc . -nostatus\n";
		`$mastpath/mast $motifsfile $filename.ORs.starRemoved.pep -oc . -nostatus`;
		if(!-z "mast.xml")
		{
			`$mastpath/xsltproc_lite $mast_xslt_path mast.xml mast.html`;
			print "Motif search finished.\n\n";
		}
		else
		{
			print "Motif search using MAST could not be executed, probably due to wrong motif input file\n";
		}
		POSIX::_exit(0);
	}
	else
	{}
}


my $tmhmmm_selected = 0;
if((!defined $tmhmm) || (!defined $tmhmmpath) || (! -e "$tmhmmpath\/tmhmm") || (-z "$tmhmmpath\/tmhmm") || (! -x "$tmhmmpath\/tmhmm"))
{
	print "Tmhmm option not chosen OR file path not provided OR file does not exist OR it is empty OR is not an executable. Skipping the search for transmembrane helices using TMHMM in the protein sequences predicted by this program. This will not affect the core analysis.\n\n";
}
else
{
	$tmhmmm_selected = 1;
	print "TMHMM search started...\n\n";
	my $pid = fork();
	if (!defined $pid)
	{
		warn "Fork failed: $!\n";
	}
	elsif($pid==0)
	{
		`$tmhmmpath/tmhmm $filename.ORs.starRemoved.pep -short > $filename.ORs.starRemoved.pep.tmhmmout`;
		print "TMHMM search finished.\n\n";
		POSIX::_exit(0);
	}
	else
	{}
    
}

my $hmmtop_selected=0;
if((!defined $hmmtop) || (!defined $hmmtoppath) || (! -e "$hmmtoppath\/hmmtop") || (-z "$hmmtoppath\/hmmtop") || (! -x "$hmmtoppath\/hmmtop"))
{
	print "Hmmtop option not chosen OR file path not provided OR file does not exist OR it is empty OR is not an executable. Skipping the search for transmembrane helices using HMMTOP in the protein sequences predicted by this program. This will not affect the core analysis.\n\n";
}
else
{
	$hmmtop_selected=1;
	print "HMMTOP search started...\n\n";
	my $pid = fork();
	if (!defined $pid)
	{
		warn "Fork failed: $!\n";
	}
	elsif($pid==0)
	{
		`export HMMTOP_ARCH=$hmmtoppath/hmmtop.arch; export HMMTOP_PSV=$hmmtoppath/hmmtop.psv; $hmmtoppath/hmmtop -if=$filename.ORs.starRemoved.pep -of=$filename.ORs.starRemoved.pep.hmmtopout;`;
		# open(OUTP,">$filename.ORs.starRemoved.pep.hmmtopout");
		# open(FHP,"<$filename.ORs.starRemoved.pep");
		# $/=">";
		# my @prot_file=<FHP>;
		# close FHP;
		# shift(@prot_file);
		# my $prot_header;
		# foreach my $protseq(@prot_file)
		# {
		# 	$protseq=~s/>//;
		# 	my @splitprot=split("\n",$protseq);
		# 	$prot_header=shift(@splitprot);
		# 	$protseq=join("",@splitprot);
		# 	#my $tmphmmtop=`wget -q -O - 'http://www.enzim.hu/hmmtop/server/hmmtop.cgi?if=$protseq&ol=1' -e use_proxy=yes -e http_proxy=proxy.ncbs.res.in:3128 | grep '^>HP: '`;
		# 	`export HMMTOP_ARCH=`
		# 	my $tmphmmtop=`-if=$protseq | grep '^>HP: '`;
		# 	$tmphmmtop=~s/noname/$prot_header/;
		# 	print OUTP "$tmphmmtop";
		# }
		# close OUTP;

		print "HMMTOP search finished.\n\n";
		POSIX::_exit(0);
	}
	else
	{}
}

my $phobius_selected = 0;
if((!defined $phobius) || (!defined $phobiuspath) || (! -e "$phobiuspath\/phobius.pl") || (-z "$phobiuspath\/phobius.pl"))
{
	print "Phobius option not chosen OR file path not provided OR file does not exist OR it is empty OR is not an executable. Skipping the search for transmembrane helices using Phobius in the protein sequences predicted by this program. This will not affect the core analysis.\n\n";
}
else
{
	print "Phobius search started...\n\n";
	$phobius_selected = 1;
	my $pid = fork();
	if (!defined $pid)
	{
		warn "Fork failed: $!\n";
	}
	elsif($pid==0)
	{
		`perl $phobiuspath/phobius.pl -short $filename.ORs.starRemoved.pep > $filename.ORs.starRemoved.pep.phobiusout`;
		print "Phobius search finished.\n\n";
		POSIX::_exit(0);
	}
	else
	{}
}

while (wait() != -1) {};

my %TMHMMTMHcount; my %HMMTOPTMHcount; my %PhobiusTMHcount;
if($tmhmmm_selected == 1 )
{
	$tmhmmpred_ref=parseTMHMM("$filename.ORs.starRemoved.pep.tmhmmout");
	%tmhmmpred_TMHs_hash=%$tmhmmpred_ref;
	foreach my $i(keys %tmhmmpred_TMHs_hash)
	{
		$TMHMMTMHcount{(@{$tmhmmpred_TMHs_hash{$i}}/2)}++;
	}
}
if($hmmtop_selected == 1)
{
	$hmmtoppred_ref=parseHMMTOP("$filename.ORs.starRemoved.pep.hmmtopout");
	%hmmtoppred_TMHs_hash=%$hmmtoppred_ref;
	foreach my $i(keys %hmmtoppred_TMHs_hash)
	{
		$HMMTOPTMHcount{(@{$hmmtoppred_TMHs_hash{$i}}/2)}++;
	}
}
if ($phobius_selected == 1) 
{
	$phobiuspred_ref=parsePhobius("$filename.ORs.starRemoved.pep.phobiusout");
	%phobiuspred_TMHs_hash=%$phobiuspred_ref;
	foreach my $i(keys %phobiuspred_TMHs_hash)
	{
		$PhobiusTMHcount{(@{$phobiuspred_TMHs_hash{$i}}/2)}++;
	}	
}

if($hmmsearch_selected == 1)
{
	open(FH3,"<$filename.ORs.starRemoved.pep.hmmsearchout");
    while(<FH3>)
    {
    	my $i=$_; chomp $i;
    	if($i=~/^\#/){}
    	else
    	{
    		if($i=~/^(.*?)\s+.*7tm_6/)
    		{
    			$hmmsearch_hash{$1}=$i;
    		}
    	}
    }
    close FH3;
}

my %consTMHcount; 
###calculating consensus###
if( ((!defined $phobius)) || (!defined $hmmtop) || (!defined $tmhmm) || (! -e "$filename.ORs.starRemoved.pep.tmhmmout") || (! -e "$filename.ORs.starRemoved.pep.hmmtopout") || (! -e "$filename.ORs.starRemoved.pep.phobiusout"))
{
	print "Either of the three transmembrane helix prediction programs could not run. Skipping the search for transmembrane helices using consensus in the protein sequences predicted by this program. This will not affect the core analysis.\n\n";
}
else
{
	print "Consensus TMH search started...\n\n";
	($tmhmmpred_ref2, $hmmtoppred_ref2, $phobiuspred_ref2, $consensuspred_ref, $consensusscore_ref, $consensusTMHs_ref, $seq_ref)=makeConsensus(\%$tmhmmpred_ref, \%$hmmtoppred_ref, \%$phobiuspred_ref, "$filename.ORs.starRemoved.pep");

	%tmhmmpred_str_hash=%$tmhmmpred_ref2; 
	%hmmtoppred_str_hash=%$hmmtoppred_ref2; 
	%phobiuspred_str_hash=%$phobiuspred_ref2; 
	%consensuspred_str_hash=%$consensuspred_ref; 
	%consensusscore_str_hash=%$consensusscore_ref; 
	%consensus_TMHs_hash=%$consensusTMHs_ref;
	%seqhash=%$seq_ref;

	open(OUT10,">$filename.ORs.starRemoved.pep.consensusTMHpredout");
	open(OUT11,">$filename.ORs.starRemoved.pep.consensusTMHpred.html");

	#include sequence later!!!
	print OUT10 "SequenceName\tNoOfHelices\tTMHs\n";
	print OUT11 "<div class=\"panel-group\" id=\"accordion\">";
	my $countnew=2;
	my $in;

	foreach my $i(sort keys %seqhash)
	{

		print OUT10 "$i\t".(@{$consensus_TMHs_hash{$i}}/2)."\t@{$consensus_TMHs_hash{$i}}\n";
		$consTMHcount{(@{$consensus_TMHs_hash{$i}}/2)}++;
		if($countnew==1){$in="in"}
		else{$in=""}
		print OUT11 "<div class=\"panel panel-default\"><div class=\"panel-heading\"><h4 class=\"panel-title\"><a data-toggle=\"collapse\" href=\"#collapse$countnew\">$i</a></h4></div><div id=\"collapse$countnew\" class=\"panel-collapse collapse $in\"><div class=\"panel-body\">\n";
		print OUT11 "<br/>";
		$countnew++;
		my $tmhmm_temp; my $hmmtop_temp; my $phobius_temp; my $consensus_temp; my $consscore_temp; my $seq_temp;
		for (my $j=0;$j<length($tmhmmpred_str_hash{$i});$j=$j+70)
		{
			$tmhmm_temp=substr($tmhmmpred_str_hash{$i},$j,70);
			$hmmtop_temp=substr($hmmtoppred_str_hash{$i},$j,70);
			$phobius_temp=substr($phobiuspred_str_hash{$i},$j,70);
			$consensus_temp=substr($consensuspred_str_hash{$i},$j,70);
			$consscore_temp=substr($consensusscore_str_hash{$i},$j,70);
			$seq_temp=substr($seqhash{$i},$j,70);

			$tmhmm_temp=assignColors($tmhmm_temp);			
			$hmmtop_temp=assignColors($hmmtop_temp);
			$phobius_temp=assignColors($phobius_temp);
			$consensus_temp=assignColors($consensus_temp);

			my $template = HTML::Template->new(filename => "$basepath/html_templates/table.tmpl");
			$template->param(
				SEQUENCEStart => $j+1,
				SEQUENCEend => $j+70,
				SEQUENCESTR => $seq_temp,
				HMMTOPSTR => $hmmtop_temp,
				TMHMMSTR => $tmhmm_temp,
				PHOBIUSSTR => $phobius_temp,
				CONSENSUSSTR => $consensus_temp,
				CONSENSUSSCORE => $consscore_temp
				);

			print OUT11 $template->output;
			print OUT11 "<br/>"
		}
		print OUT11 "</div></div></div>";
	}

	#print OUT11 "</body>\n</html>";
	print OUT11 "</div>";

	close OUT10;
	close OUT11;
	print "Consensus TMH search finished.\n\n";
}

if((%TMHMMTMHcount) || (%HMMTOPTMHcount) || (%PhobiusTMHcount) || (%consTMHcount))
{
	my $graph=barplot(\%TMHMMTMHcount,\%HMMTOPTMHcount,\%PhobiusTMHcount,\%consTMHcount);
	my $file = "$filename.bars.png";
	open(my $out, '>', $file) or die "Cannot open '$file' for write: $!";
	binmode $out;
	print $out $graph->gd->png;
	close $out;
}



open(OUTN9,">$filename.OR_final_insector_table.txt");


print OUTN9 "Name\tScaffold\tStart\tStop\tStrand\tProt Len\tPfam\tExons\tTMHMM\tHMMTOP\tPhobius\tConsensus\tGene Model\tBest Query\n";
my %ig_final_protein_hash;
open(FHN1,"<$filename.final_proteins.pep");

$/="\>";
my @ig_prot_seq_file=<FHN1>;
close FHN1;
shift @ig_prot_seq_file;
foreach(@ig_prot_seq_file)
{
	if($_=~/((.*?)_OR_.*?)\n/ms)
	{
		$ig_final_protein_hash{$2}=$1;
	}
}


open(FHN2,"<$filename.final_table.txt");
$/="\n";
while(<FHN2>)
{
	my $t_line=$_;
	chomp $t_line;
	my @t_line_arr=split("\t", $t_line);
	my $name=$t_line_arr[0]."_".$t_line_arr[1]."-".$t_line_arr[2];
	print OUTN9 "$ig_final_protein_hash{$name}\t$t_line_arr[0]\t$t_line_arr[1]\t$t_line_arr[2]\t$t_line_arr[3]\t$t_line_arr[4]\t";

	 if(exists $hmmsearch_hash{$ig_final_protein_hash{$name}})
	 {
	 	#$out8str.="\;Pfam=7tm_6";
	 	#$tempout8.="\;Pfam=7tm_6";
	 	print OUTN9 "7tm_6\t" ;

	 	$pfam_count++;
	 }
	 else
	 {
	 	#$out8str.="\;Pfam=-";
	 	#$tempout8.="\;Pfam=-";
	  	print OUTN9 "-\t" ;
	  	$notpfam_count++;
	 }
	my $exo_length=split(/\s/, $t_line_arr[7]); 
	$exo_length=$exo_length/2;
	print OUTN9 "$exo_length\t";
	#$out9str.=(@cdsboundary/2)."\t";

	 if(%tmhmmpred_TMHs_hash)
	 {
	 	print OUTN9 (@{$tmhmmpred_TMHs_hash{$ig_final_protein_hash{$name}}}/2)."\t";
	 }else{print OUTN9 "-\t"}

	 if(%hmmtoppred_TMHs_hash && (exists $hmmtoppred_TMHs_hash{$ig_final_protein_hash{$name}}))
	 {
	 	print OUTN9 (@{$hmmtoppred_TMHs_hash{$ig_final_protein_hash{$name}}}/2)."\t";
	 }else{print OUTN9 "-\t"}

	 if(%phobiuspred_TMHs_hash)
	 {
	 	print OUTN9 (@{$phobiuspred_TMHs_hash{$ig_final_protein_hash{$name}}}/2)."\t";
	 }else{print OUTN9 "-\t"}

	 if(%consensus_TMHs_hash)
	 {
	 	print OUTN9 (@{$consensus_TMHs_hash{$ig_final_protein_hash{$name}}}/2)."\t";
	 }else{print OUTN9 "-\t"}	

	 my @t_line_exons=split(/\s+/,$t_line_arr[7]);
	 my $exon_temp="";
	 for(my $i=0; $i<@t_line_exons; $i+=2)
	 {
		$exon_temp.= "$t_line_exons[$i]..$t_line_exons[$i+1],";
	 }
	 chop $exon_temp;
	 print OUTN9 "$exon_temp\t$t_line_arr[8]\n";
}
close FHN2;
close OUTN9;

#######################################################################################################################################################################################
#############################################################GFF Comparison##########################################################################################################################

my $gfffile_exists=0;
my $gffbedfile_exists=0;
my %overlappingGFFgenes;
my $overlapping_genes_no=0;
#my $abc=1;

if((!defined $gfffile) || (! -e $gfffile) || (-z $gfffile))
{
	print "GFF annotation file not provided OR does not exist OR it is empty. Skipping the comparison of OR annotation results from this program with overlapping annotations from another source. This will not affect the core analysis.\n";
}
else
{
	#`awk -F'\t' -v OFS='\t' '{sub(/_/, "", $1)} 1' $gfffile > $gfffile.new`;
	#`mv $gfffile.new $gfffile`;
	$gfffile_exists=1;
	print "Comparing with user-provided gff file...\n\n";
	EXIT_IF:{
		my %gffhash;
		my %genegffhash;
		my %gffhash_suppl;
		my @gene_suppl;
		open(FH2,"<$gfffile");
		#open(FH2,"<ncbi.gff");

		$/="\n";
		while(<FH2>)
		{
			chomp $_;
			my $tempvar=$_;
			if($tempvar=~/^(.*?)\s+.*?\s+.*?\s+.*?\s+.*?\s+/)
			{
				my $scaftemp=$1;
				if(exists $scaf_length_bb{$scaftemp})
				{}
				else
				{
					$scaf_length_bb{$scaftemp}=1000000000;
				}
			}

			if($tempvar=~/^(.*?)\s+.*?\s+exon\s+(.*?)\s+(.*?)\s+.*?\s+(.*?)\s+.*?\s+.*?(gene|Name)\=(.*?)\;.*?product\=(.*?)\;.*?protein_id\=(.*?)(;|$)/)
			{
				push(@{$gffhash{$1."~".$6."~".$8."~".$7."~".$4}},$2);
				push(@{$gffhash{$1."~".$6."~".$8."~".$7."~".$4}},$3);
				$gffhash_suppl{$6}.=$tempvar."\n";
			}
			elsif($tempvar=~/^(.*?)\s+.*?\s+exon\s+(.*?)\s+(.*?)\s+.*?\s+(.*?)\s+.*?\s+.*?Parent\=(.*?)(;|$)/)
			{
				push(@{$gffhash{$1."~".$5."~".$4}},$2);
				push(@{$gffhash{$1."~".$5."~".$4}},$3);
				$gffhash_suppl{$5}.=$tempvar."\n";
			}
			elsif($tempvar=~/^(.*?)\s+.*?\s+gene\s+(.*?)\s+(.*?)\s+/)
			{
				$genegffhash{"$1.$2.$3"}=$tempvar."\n";
			}
			elsif($tempvar=~/^(.*?)\s+.*?\s+pseudogene\s+(.*?)\s+(.*?)\s+/)
			{
                                $genegffhash{"$1.$2.$3"}=$tempvar."\n";
                        }
			elsif($tempvar=~/^(.*?)\s+.*?\s+region\s+(.*?)\s+(.*?)\s+/)
			{
				my $tmplength=$3-$2+1;
				$scaf_length_bb{$1}=$tmplength;
			}
			else
			{}
		}

		#print "%gffhash created.\n";
		if((keys %gffhash)==0)
		{
			print "Either 1) The 'exon' tag is missing in the third column of the input gff file or 2) None of the usual tags \(gene OR name OR parent\) were found in the input gff file. Skipping the comparison of OR annotation results from this program with overlapping annotations from another source.\n";
			last EXIT_IF;
		}
		
		`$gff_to_bed_loc/gff_to_bed.py $gfffile | awk 'tolower(\$4) ~ /rna/' | sort -k1,1 -k2,2n | sed 's/\$/\tprotein_coding/' | sed 's/\t\\.\t/\t0\t/' > $gfffile.ORs_sorted.bed12`;
		$gffbedfile_exists=1;

		open (OUT4,">$filename.gffcomparison");
		open (OUT5,">$filename.ORrelated_genesFromUserProvidedAnnotation.gff");

		foreach my $i(sort keys %gffhash)
		{
			@{$gffhash{$i}}=sort{$a <=> $b}@{$gffhash{$i}};
		}

		#print "Comparing GWS OR boundaries with preexisting GFF gene boundaries...\n";
		my %unique_hash;

		open (FHN4,"<$filename.OR_final_insector_table.txt");

		my @or_table_final=<FHN4>;
		close FHN4;

		shift @or_table_final;
		my $w=0;
		foreach(@or_table_final)
		{
			$w++;
			my @line=split("\t","$_");
		 	print OUT4 $line[0]."\n";
			print OUT4 "InsectOR annotation:	    ";
			#my @cdstemp=split(/\,/,$line[-2]);
			#my $cdstempstr;
			#for(my $cdsno=0; $cdsno<@cdstemp; $cdsno+=2)
			#{
			#	$cdstempstr.=$cdstemp[$cdsno]."..".$cdstemp[$cdsno+1].",";
			#}

			#$cdstempstr=~s/\,$//;
			#print OUT4 "$cdstempstr\n";

			print OUT4 "$line[-2]\n";
			my $scafname=$line[1];
			foreach my $j(keys %gffhash)
			{
				#print "gffhash key is $j\n";
				if($j=~/^$scafname~(.*?)~/)
				{
					my $loc=$1;
					#print "loc is $loc\n";
					if(${$gffhash{$j}}[0]==$line[2] || ${$gffhash{$j}}[-1]==$line[3] || (${$gffhash{$j}}[0]>$line[2] && ${$gffhash{$j}}[0]<$line[3]) || (${$gffhash{$j}}[-1]>$line[2] && ${$gffhash{$j}}[-1]<$line[3]) || (${$gffhash{$j}}[0]<$line[2] && ${$gffhash{$j}}[-1]>$line[2]) || (${$gffhash{$j}}[0]<$line[3] && ${$gffhash{$j}}[-1]>$line[3]) )
					{
						my $temp1;
						#print "$scafname.${$gffhash{$j}}[0].${$gffhash{$j}}[-1]\n";
						if (exists $genegffhash{"$scafname.${$gffhash{$j}}[0].${$gffhash{$j}}[-1]"})
						{
							$temp1 = $genegffhash{"$scafname.${$gffhash{$j}}[0].${$gffhash{$j}}[-1]"};
							$overlappingGFFgenes{$w}=1;
						}
						else
						{
							#my $assignedflag=0;
							my @tps;
							foreach my $i(keys %genegffhash)
							{
								@tps=split(/\s+/,$genegffhash{$i});
								if($tps[0]=~/$scafname/)
								{
									if($tps[3]<=${$gffhash{$j}}[0] && $tps[4]>=${$gffhash{$j}}[-1])
									{
										$temp1 = $genegffhash{$i};
										$overlappingGFFgenes{$w}=1;
										#$assignedflag=1;
									}
								}
								else
								{}
							}
							# if ($assignedflag==0)
							# {
							# 	print "$tps[0], $tps[3], $tps[4], ${$gffhash{$j}}[0], ${$gffhash{$j}}[-1]\n"; 
							# }
						}
						#$temp1=~s/\tgene\t/\tprotein_match\t/g;
						#if($temp1 eq ""){print $grepstr}
						#my $temp2=`grep $loc $gfffile`;
						my $temp2=$gffhash_suppl{$loc};
						$unique_hash{$temp1}=0;
						$unique_hash{$temp2}=0;
						print OUT4"User-provided annotation:   ";
						my $cdsstring=${$gffhash{$j}}[0]."..".@{$gffhash{$j}}[1];
						for (my $i=2;$i<@{$gffhash{$j}};$i=$i+2)
						{
							$cdsstring.=",".${$gffhash{$j}}[$i]."..".@{$gffhash{$j}}[$i+1];
						}
						print OUT4 "$cdsstring\n$j\n";
					}
				}
			}
			print OUT4 "\n";			
		}

		foreach my $z(sort keys %unique_hash)
		{
			print OUT5 $z;
		}

		close FH2; close OUT4; close OUT5;

		print "Comparison with user-provided gff completed.\n\n";

		my $filename3=$filename2."ORrelated_genesFromAutomatedAnnotation";
		$overlapping_genes_no = (keys %overlappingGFFgenes);
		# my $non_overlapping_genes= $gene_count - $overlapping_genes_no;
		# $summary.="OR gene regions non-overlapping with any user provided genes (Completely novel genes) - $non_overlapping_genes\n";
		# $summary.="OR gene regions overlapping with any user provided genes - $overlapping_genes_no\n(The gene models may not match perfectly. The user provided gene might be non-OR gene as well.)\n\n";
	}
}

#######################################################################################################################################################################################
#############################################################GFF Comparison over##########################################################################################################################


`$gff_to_bed_loc/gff_to_bed.py $filename.final_gff_file.gff | sort -k1,1 -k2,2n | sed 's/\$/\tprotein_coding/' | sed 's/\t\\.\t/\t0\t/' > $filename.ORs_sorted.bed12`;


open(FHN102,"$filename.OR_final_insector_table.txt");
$/="\n";
my @insector_final_table=<FHN102>;
close FHN102;
shift @insector_final_table;

open(OUTN102,">$filename.final_proteins.cds");
foreach (@insector_final_table)
{
	my @line=split(/\t/,$_);
	my @exons=split(/\,/,$line[-2]);
	my $cds="";
	for(my $i=0; $i<@exons; $i+=1)
	{
		my @ex=split(/\.\./,$exons[$i]);
		$cds.=substr($genomeSeqHash{">".$line[1]},$ex[0]-1,($ex[1]+1-$ex[0]))
	}
	$cds=uc $cds;
	if($line[4] eq "reverse")
	{
		$cds=reverse $cds;
		$cds=~tr/ATGC/1234/;
		$cds=~tr/1234/TACG/;
	}
	print OUTN102 ">$line[0]\n$cds\n";
}
close OUTN102;

open(OUTN100,">$filename.predictedORs.summary.txt");

print OUTN100 "Total number of genes: ";
my $t=`wc -l <$filename.OR_final_insector_table.txt`;
$t=$t-1;
print OUTN100 $t,"\n";

print OUTN100 "Total number of genes overlapping with user-provided genes: ", $overlapping_genes_no,"\n";

print OUTN100 "Total number of genes non-overlapping with user-provided genes: ", $t-$overlapping_genes_no,"\n";

print OUTN100 "Total number of Complete genes: ";

print OUTN100 `grep -c "Complete" $filename.OR_final_insector_table.txt`;

print OUTN100 "Total number of Partial genes: ";

print OUTN100 `grep -c "Partial" $filename.OR_final_insector_table.txt`;

print OUTN100 "Total number of Normal genes: ";

print OUTN100 `grep -c "Normal" $filename.OR_final_insector_table.txt`;

print OUTN100 "Total number of Pseudo genes: ";

print OUTN100 `grep -c "Pseudo" $filename.OR_final_insector_table.txt`;

print OUTN100 "Total number of genes with start codon present: ";

print OUTN100 `grep -c "StartCodonPresent" $filename.OR_final_insector_table.txt`;

print OUTN100 "Total number of genes without start codon: ";

print OUTN100 `grep -c "StartCodonAbsent" $filename.OR_final_insector_table.txt`;

print OUTN100 "Total number of genes with 7tm_6 domain: ";

#print OUTN100 `grep -c "7tm_6" $filename.OR_final_insector_table.txt`;

my $b=`grep -c "7tm_6" $filename.OR_final_insector_table.txt`;
print OUTN100 $b;
print OUTN100 "Total number of genes without 7tm_6: ";


my $a= int($t) - int($b);

print OUTN100 $a,"\n\n";

print OUTN100 "**********************************Within 7tm_6 containg hits*********************************\n";
print OUTN100 "Total number of Complete genes: ";
print OUTN100 `grep "7tm_6" $filename.OR_final_insector_table.txt|grep -c "Complete" `;

print OUTN100 "Total number of Partial genes: ";
print OUTN100 `grep "7tm_6" $filename.OR_final_insector_table.txt|grep -c "Partial" `;

print OUTN100 "Total number of Normal genes: ";
print OUTN100 `grep "7tm_6" $filename.OR_final_insector_table.txt|grep -c "Normal" `;

print OUTN100 "Total number of Pseudo genes: ";
print OUTN100 `grep "7tm_6" $filename.OR_final_insector_table.txt|grep -c "Pseudo" `;

print OUTN100 "Total number of genes with start codon present: ";
print OUTN100 `grep "7tm_6" $filename.OR_final_insector_table.txt|grep -c "StartCodonPresent" `;

print OUTN100 "Total number of without start codon: ";
print OUTN100 `grep "7tm_6" $filename.OR_final_insector_table.txt|grep -c "StartCodonAbsent" `;


print OUTN100 "\n";

close OUTN100;
print "\n\nAdditional user-selected analyses over.\n\n";
print "Pipeline completed.\n\nCleaning up.\n";

`rm $filename.modified_seqfile`;
###
#$filename.geneclusters.0.txt

`rm $filename.combinedGeneWiseRes_7tm6_nonfiltered`;
`rm $filename.consolidatedproteins_7tm6notfiltered.txt`;
`rm $filename.consolidatedTable_1_7tm6notfiltered.txt`;
`rm $filename.final_table.txt`;
`rm $filename.genewise_protein_seq.txt`;
`rm $filename.genewise.table`;
`rm $filename.latest_consolidatedproteins_7tm6notfiltered.txt`;
`rm $filename.latest_consolidatedTable_1_7tm6notfiltered.txt`;
`rm $filename.mergedTable_1_withoverlaps_7tm6notfiltered.txt`;
`rm $filename.OR_table_sorted.txt`;
`rm $filename.OR_table.txt`;
`rm $filename.ORs.pep`;
`rm $filename.ORs.cds`;

exit;

#######################################################################################################################################################################################
##### INSECTOR - pipeline over ########################################################################################################################################################
#######################################################################################################################################################################################
