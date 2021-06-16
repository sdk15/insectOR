InsectOR package
================
*Authors: Snehal Dilip Karpe and Vikas Tiwari*


This package helps to annotate genes coding for any specific family of proteins from a medium sized genome (currently tested on insect olfactory receptors / ORs). 

The web version of the package can be accessed at http://caps.ncbs.res.in/insectOR/ 

More details about the package can be found at - http://caps.ncbs.res.in/cgi-bin/gws_ors/about.py?help=about&help_t=About%20insectOR

This package is under development.


Operating system requirements:
------------------------------
This package is made for use on unix systems and requires python2.7.


Installation:
-------------
Please follow these instructions for installing and using this package -

1) Extract the contents of this folder.

2) Make sure that Perl programming language (tested on v5.26.1) and the following perl modules are installed -
```	
	Getopt::Long qw(GetOptions)	
	Tie::IxHash
	HTML::Template
	File::Basename
	File::Spec
	IO::File
	GD::Graph::bars
	GD::Graph::Data
	IO::Handle
	Cwd
```
3) Install following packages and keep their respective folders/binaries inside ```insectOR/tools``` or ```insectOR-main/tools```folder-

	A) GeneWise ( e.g. wise2.4.1 - https://www.ebi.ac.uk/~birney/wise2/) 
	  
	  - Copy the ```wiseX.X.X``` folder inside ```insectOR/tools```
	  - Change the ```$wiseconfigdir_path``` in the ```bin/scoreGenesOnScaffold.pl``` accordingly (If the GeneWise package wise2.4.1 is installed at the correct location as given above there is no need to edit ```$wiseconfigdir_path```).
		
	B) GFFtools-GX (e.g. GFFtools-GX-master - https://github.com/vipints/GFFtools-GX)
		
	  - copy the ```GFFtools-GX-master``` folder inside ```insectOR/tools```
		
	C) TMHMM2 (e.g. tmhmm-2.0c - https://services.healthtech.dtu.dk/software.php) 
		
	  - mandatory if using option ```-tmhmm``` or ```-tmh1```
	  - change the ```$tmhmmpath``` in the ```bin/scoreGenesOnScaffold.pl``` accordingly (If the TMHMM2 package tmhmm-2.0c is installed at the correct location as given above there is no need to edit ```$tmhmmpath```).
	  - make sure to give executable permission
		
	D) HMMTOP2 (e.g. hmmtop_2.1 - http://www.enzim.hu/hmmtop/html/download.html) 
		
	  - mandatory if using option ```-hmmtop``` or ```-tmh2```
	  - change the ```$hmmtoppath``` in the ```bin/scoreGenesOnScaffold.pl``` accordingly (If the HMMTOP2 package hmmtop_2.1 is installed at the correct location as given above there is no need to edit ```$hmmtoppath```).
	  - make sure to give executable permission
		
	E) Phobius (http://phobius.sbc.su.se/data.html)
	
	  - mandatory if using option ```-phobius``` or ```-tmh3```
	  - check of the ```$phobiupath``` is set correctly in ```bin/scoreGenesOnScaffold.pl```
	  - make sure to give executable permission
		
	F) HMMER (e.g. hmmer-3.1b2 - http://hmmer.org/download.html)
	
	  - mandatory if using option ```-hmmsearch``` or ```-p```
	  - change the ```$hmmsearchpath``` in the ```bin/scoreGenesOnScaffold.pl``` accordingly (If the hmmerpackage hmmer-3.1b2 is installed at the correct location as given above there is no need to edit ```$hmmsearchpath```).
		
	G) MEME (e.g. meme_4.10.2 - http://meme-suite.org/doc/download.html) 
	
	  - mandatory if using option ```-mast``` or ```-m```
	  - change the ```$mastpath``` and ```$mast_xslt_path``` in the ```bin/scoreGenesOnScaffold.pl``` accordingly (If the MEME suit meme_4.10.2 is installed at the correct location as given above there is no need to edit ```$hmmsearchpath```).

4) Change the ```$basepath``` in ```bin/scoreGenesOnScaffold.pl``` as per the location of ```insectOR``` folder on your system.

5) Download and keep '7tm_6.hmm' inside ```insectOR/hmm``` folder (https://pfam.xfam.org/family/PF02949/hmm)
	
	  - mandatory if using option ```-hmmsearch``` or ```-p```

You are ready to use insectOR!


Notes
-----
1) The exonerate is run using following parameters - 
```   
exonerate --model protein2genome --maxintron <max intron size dependent on the organism> <protein query sequence fasta file> <genome fasta file> --showtargetgff TRUE
```
2) The exonerate alignment contains 'Command line' and 'exonerate' in first line.


Example command
---------------
```
perl (path_to_insectOR/)insectOR/bin/scoreGenesOnScaffold.pl -i exonerate.txt -s seq.fasta -q query.fasta 
perl (path_to_insectOR/)insectOR/bin/scoreGenesOnScaffold.pl -i exonerate.txt -s seq.fasta -q query.fasta -g ncbi.gff -tmh1 -tmh2 -tmh3 -p -m -mf motif.txt 
```

Input
-----
* Mandatory
```
-exonerate_file|-i - Exonerate alignment file 
```
(Sample exonerate file is provided in the test folder - ```exonerate.txt```)
```
-seq|-s            - Genome sequence (FASTA) file which was used to generate above exonerate file 
```
(Sample genome sequence file is provided in the test folder - ```seq.fasta```. Please note that this is a *Habropoda laboriosa* genome scaffold sequence NCBI Genbank accession LHQN01028732.1.)
```
-queryseq|-q       - OR/query sequence (FASTA) file which was used for generating above exonerate alignment 
````
(Sample OR sequence file is provided in the test folder - ```query.fasta```)

* Optional
```
-gff_file|-g       - User provided gene annotations (GFF format) with which InsectOR output will be compared 
```

(Sample genome annotation file is provided in the test folder - ```ncbi.gff```. Please note that this is a *Habropoda laboriosa* scaffold gene annotation file in GFF format for NCBI Genbank sequence LHQN01028732.1.)

```
-cutoff|-c         - Alignment clusters are identified based on this cutoff. 
                     This is the minimum number of alignments needed at a nucleotide position for its inclusion into an alignment cluster. (Default: 1)
-lengthCutoff|-l   - Predicted proteins can be classified as complete or partial based on this cutoff. (Default: 300 amino acids)
-hmmsearch|-p      - Perform HMMSEARCH for insect olfactory receptor signature (PFAM 7tm_6 family) as additional validation.
-tmhmm|-tmh1       - Search for Transmembrane Helices (TMH) using TMHMM2 as additional validation.
-hmmtop|-tmh2      - Search for Transmembrane Helices (TMH) using HMMTOP2 as additional validation.
-phobius|-tmh3     - Search for Transmembrane Helices (TMH) using Phobius as additional validation.
```
(If all three TMH predictiors are selected, consensus TMH prediction will be performed.)
```
-mast|-m           - Search for known OR motifs in the predicted proteins using MAST motif tool
```
```
-motif_file|-mf    - Users can provide their own motifs (PSPM format of MEME) for MAST motif search into the predicted OR proteins. 
```
(If MAST option is checked without any file mentioned with this parameter, then default AfOR motif file will be used)                    
(Sample motif file is provided in the test folder - motif.txt)
```
-helpmessage|-h    - Print this help message
```

Output
------
```
...predictedORs.summary.txt    - Summary file of insectOR results
...OR_final_insector_table.txt - Table providing detailed information on each gene/fragment predicted by insectOR
...final_proteins.pep          - Protein sequence file of insectOR predicted ORs/fragments in FASTA format
...ORs.starRemoved.pep         - Protein sequence file of insectOR predicted ORs/fragments without any pseudogenizing elements in FASTA format
...ORs.cds                     - Nucleotide sequence file of insectOR predicted OR CDSs in FASTA format
...final_gff_file.gff          - Detailed gene annotations by insectOR in GFF format
...ORs_sorted.bed12            - Detailed gene annotations by insectOR in BED format
```
If gene annotation file from another resource is provided, following files are also generated - 
```
...gffcomparison                                 - Comparison of gene annotations from insectOR along with overlapping user provided gene annotations. 
...ORrelated_genesFromUserProvidedAnnotation.gff - A shorter version of user provided GFF file containing only overalapping annotations with those from this tool.
...gff.ORs_sorted.bed12                          - sorted BED formatted version of transcripts in the user-provided GFF file.
```
If HMMSEARCH option against 7tm_6 is selected, following files are also generated -
```
...ORs.starRemoved.pep.hmmsearchout              - HMMSEARCH output against 7tm_6 HMM file
```
If either of TMH prediction tools are used, following files are also generated -
```
...ORs.starRemoved.pep.tmhmmout if '-tmh1' is selected   - Output of TMHMM2 tool
...ORs.starRemoved.pep.hmmtopout if '-tmh2' is selected  - Output of HMMTOP2 tool
...ORs.starRemoved.pep.phobiusout if '-tmh3' is selected - Output of Phobius tool
...ORs.starRemoved.pep.consensusTMHpredout and ...ORs.starRemoved.pep.consensusTMHpred.html if all the three '-tmh1 -tmh2 -tmh3' are selected 
                                                         - Output of consensus TMH prediction tool.
```
If MAST motif search is selected, following files are also generated -
```
MAST.html - Motif search output file in HTML format
MAST.xml  - Motif search output file in HTML format
```


Suggestions and known bugs
--------------------------
* Please use simple names for your genome sequences. Issues might occur if longer sequence lengths and special characters including _ or - are part of these names. 
* The sequence names should match across the Exonerate and the sequence files for genome and protein queries.
* InsectOR performs several parallel executions (upto 6 threads) of the validation packages near the end of the execution. 
* Please be warned that this might need computational resources.
* Run the InsectOR in a folder with just input files and avoid having multiple runs of InsectOR within the same folder.

Please cite
-----------
* Karpe SD, Tiwari V, Ramanathan S (2021) InsectOR—Webserver for sensitive identification of insect olfactory receptor genes from non-model genomes. PLoS ONE 16(1): e0245324. https://doi.org/10.1371/journal.pone.0245324 (Updated from : Karpe SD, Tiwari V & Sowdhamini R. InsectOR - webserver for sensitive identification of insect olfactory receptor genes from non-model genomes. bioRxiv doi: https://doi.org/10.1101/2020.04.29.067470)
* Please also cite respective papers for tools used herewith - e.g. GeneWise, GFFtools, TMH prediction methods, 7tm_6 hmmsearch, MAST motif search tool, etc.


Contact
-------
Please contact karpesnehal@gmail.com or vikast@ncbs.res.in to report any bugs.
