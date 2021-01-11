package My::utils;
use strict;
#use warnings;
use warnings FATAL => 'uninitialized';

use Exporter qw(import);
 
our @EXPORT_OK = qw(processAlignment processEveryAlignment findMet processFinalNew processSimilarity parseTMHMM parseHMMTOP parsePhobius makeConsensus assignColors joinAlignments three2one barplot);

##### Process selected alignments and get their alignment details #####
sub processAlignment{
	my $s=$_[0];
	my $queryProt; my $similarity; my $targetProt; my $targetNuc; my $spaces; my $cds;
	if ($s=~/Target range: .*? -> .*?\n\n(.*?)vulgar.*?\# seqname source feature start end score strand frame attributes\n\#\n(.*?)\# --- END OF GFF DUMP ---/ms)
	{
		my $intron_string=0;
		my @ali=split("\n",$1);
		my @cdss=($2=~/\tcds\t(.*?)\t(.*?)\t/msg);
		@cdss=sort {$a <=> $b} @cdss;
		$cds=join("-",@cdss);
		#print "CDS are @cdss\n";
		#print "$ali[0]\n";
		for (my $t=0;$t<@ali-3;$t=$t+5)
		{
			#print "$ali[$t]\n$ali[$t+1]\n$ali[$t+2]\n$ali[$t+3]\n"; 
			if ($ali[$t]=~/^(\s+\d+\s+\:\s)(.*)\s\:\s+\d+$/)
			{$queryProt.=$2; $spaces=length($1)}
			else{print"ErrorA!\n"}
			if ($ali[$t+1]=~/^\s{$spaces}(.*)$/)
			{$similarity.=$1}
			else{print"ErrorB!\n"}
			if ($ali[$t+2]=~/^\s{$spaces}(.*)$/)
			{$targetProt.=$1}
			else{print"ErrorC!\n"}
			if ($ali[$t+3]=~/^\s+\d+\s+\:\s(.*)\s\:\s+\d+$/)
			{$targetNuc.=$1}
			else{print"ErrorD!\n"}
		}

		$queryProt=~s/\{//g;
		$queryProt=~s/\}//g;
		while($queryProt=~/(\s\s>>>> Target Intron \d+ >>>>\s\s)/)
		{
			my $sstr_intron=$1;
			#elsif(length($1) == 2 ){$sstr="            "}
			if($similarity=~/(\d+)\sbp/)
			{
				my $digits=$1;
				#print "The intron is a $digits number\n";
				my $intro_length=length($1)+3;
				#print "Into length is $intro_length\n";
				my $tot_length=length($sstr_intron);
				#print "Tot length is $tot_length\n";
				my $ca=(($tot_length - $intro_length)/2);
				$ca= int($ca);
				#print "ca is $ca\n";
				my $lsstr=" " x int(($tot_length - $intro_length)/2);
				#print "lsstr is *$lsstr*\n";
				my $sstr=" " x ($tot_length-$intro_length-length($lsstr));
				#print "sstr is *$sstr*\n";

				$similarity=~s/($sstr\d+\sbp$lsstr)//;
			}
	        $queryProt=~s/  >>>> Target Intron \d+ >>>>  //;
		}
		$queryProt=~s/  >>>> Target Intron \d+ >>>>  //g;
		$queryProt=~s/ //g;
		#print "$queryProt\n";
		$targetNuc=~s/ //g;
		$targetNuc=~s/\{//g;
		$targetNuc=~s/\}//g;
		$targetNuc=~s/(.)(.)\.{1,}(.)(.)//g;
		#print "$targetNuc\n";
		$similarity=~s/\{//g;
		$similarity=~s/\}//g;
		#print "$similarity\n";#currently wrong - may remove some white spaced entries in alignment;
		$targetProt=~s/\}//g;
		$targetProt=~s/\{//g;
		$targetProt=~s/..\s+..//g;
		#print "$targetProt\n";
		if(length($queryProt)==length($targetNuc) && length($targetNuc)==length($similarity) && length($similarity)==length($targetProt))
		{}
		else
		{
			print "Error! Lengths are not similar!\n$queryProt\n$similarity\n$targetProt\n$targetNuc\n$s\n";
			exit;
		}
	}
	else
	{
	print "$s\nError!\n";
	}

	return(($queryProt,$targetNuc,$similarity,$targetProt,$cds));
}

sub processSimilarity{
	(my $sim_prev_over_reg, my $sim_over_reg)=@_;
	#$sim_prev_over_reg=~tr/\|/4/;$sim_prev_over_reg=~tr/\!/3/;$sim_prev_over_reg=~tr/\:/2/;$sim_prev_over_reg=~tr/\./1/;$sim_prev_over_reg=~tr/ /0/;
	$sim_prev_over_reg=~tr/\|/5/;$sim_prev_over_reg=~tr/\!/4/;$sim_prev_over_reg=~tr/\:/3/;$sim_prev_over_reg=~tr/\./2/;$sim_prev_over_reg=~tr/ /1/;$sim_prev_over_reg=~tr/#/0/;

    #$sim_over_reg=~tr/\|/4/;$sim_over_reg=~tr/\!/3/;$sim_over_reg=~tr/\:/2/;$sim_over_reg=~tr/\./1/;$sim_over_reg=~tr/ /0/;
   	$sim_over_reg=~tr/\|/5/;$sim_over_reg=~tr/\!/4/;$sim_over_reg=~tr/\:/3/;$sim_over_reg=~tr/\./2/;$sim_over_reg=~tr/ /1/;$sim_over_reg=~tr/#/0/;


    my @sim_over_reg_count=split("",$sim_over_reg);
    my @sim_prev_over_reg_count=split("",$sim_prev_over_reg);
    my $prev_count= eval join '+', @sim_prev_over_reg_count;
    my $curr_count= eval join '+', @sim_over_reg_count;
    #print "The sim_prev_over_reg is $sim_prev_over_reg and sim_over_reg is $sim_over_reg\n";
    #print "The sim_prev_over_reg_count is $prev_count and sim_over_reg_count is $curr_count\n";
	return(($prev_count,$curr_count));
}

sub translateNucToProt{
	my $str=$_[0];
	if((length($str) % 3) != 0){print "The nucelotide length is not in multiples of three\n$str\n"; exit;}
	#what if 2 nucleotides left at the end? check algo
	my %codon2aa=("TTT" => "F", "TTC" => "F", "TTA" => "L", "TTG" => "L", "TCT" => "S", "TCC" => "S", "TCA" => "S", "TCG" => "S", "TAT" => "Y", "TAC" => "Y", "TAA" => "*", "TAG" => "*", "TGT" => "C", "TGC" => "C", "TGA" => "*", "TGG" => "W", "CTT" => "L", "CTC" => "L", "CTA" => "L", "CTG" => "L", "CCT" => "P", "CCC" => "P", "CCA" => "P", "CCG" => "P", "CAT" => "H", "CAC" => "H", "CAA" => "Q", "CAG" => "Q", "CGT" => "R", "CGC" => "R", "CGA" => "R", "CGG" => "R", "ATT" => "I", "ATC" => "I", "ATA" => "I", "ATG" => "M", "ACT" => "T", "ACC" => "T", "ACA" => "T", "ACG" => "T", "AAT" => "N", "AAC" => "N", "AAA" => "K", "AAG" => "K", "AGT" => "S", "AGC" => "S", "AGA" => "R", "AGG" => "R", "GTT" => "V", "GTC" => "V", "GTA" => "V", "GTG" => "V", "GCT" => "A", "GCC" => "A", "GCA" => "A", "GCG" => "A", "GAT" => "D", "GAC" => "D", "GAA" => "E", "GAG" => "E", "GGT" => "G", "GGC" => "G", "GGA" => "G", "GGG" => "G", "NNN" => "X", "---" => "");
	$str=uc($str);
	if($str=~/N/)
	{
		my @tp=($str=~/(.{1,3})/g);
		my $temp_str;
		foreach my $x(@tp)
		{
			if($x=~/N/)
			{
				$x="X";
			}
			else
			{
				$x=$codon2aa{$x};
			}
			$temp_str.=$x;
		}
		$str=$temp_str;
	}
	else
	{
		$str=~s/(...)/$codon2aa{$1}/g;
	}
	#print "Middle seq is $str\n";
	return($str);
}

sub three2one{
	my $str=$_[0];
	#print "*$str*\n";
	my %three2oneAA=( "Ala"=>'A',  "Arg"=>'R',  "Asn"=>'N',  "Asp"=>'D',  "Cys"=>'C',  "Glu"=>'E',  "Gln"=>'Q',  "Gly"=>'G',  "His"=>'H',  "Ile"=>'I',  "Leu"=>'L',  "Lys"=>'K',  "Met"=>'M',  "Phe"=>'F',  "Pro"=>'P',  "Ser"=>'S',  "Thr"=>'T',  "Trp"=>'W',  "Tyr"=>'Y',  "Val"=>'V',  "Sec"=>'U',  "Pyl"=>'O', "Unk"=>'X', "---"=>'', "***"=>'*');
	#$str=~s/(#)#/$three2oneAA{$1}/g;
	#$str=~s/(#)/$three2oneAA{$1}/g;
	$str=~s/#{1,}/***/g;
	$str=~s/(...)/$three2oneAA{$1}/g;

	#print "$str\n";
	#print "3-2-1 protein is $str\n";
	return($str);
}

sub one2three{
	my $str=$_[0];
	my %one2threeAA=('A' => "Ala", 'R' => "Arg", 'N' => "Asn", 'D' => "Asp", 'C' => "Cys", 'E' => "Glu", 'Q' => "Gln", 'G' => "Gly", 'H' => "His", 'I' => "Ile", 'L' => "Leu", 'K' => "Lys", 'M' => "Met", 'F' => "Phe", 'P' => "Pro", 'S' => "Ser", 'T' => "Thr", 'W' => "Trp", 'Y' => "Tyr", 'V' => "Val", 'U' => "Sec", 'O' => "Pyl", 'X' => "Unk", '*' => "***");
	$str=~s/(.)/$one2threeAA{$1}/g;
	return($str);
}

sub startCodonPresence{
	my $str=$_[0];
	my $decision;
	if($str=~/^M/)
	{
		$decision=1;
	}
	else
	{
		$decision=0;
	}
	return($decision)
}

sub pseudoPresence{
	my $str=$_[0];
	my $decision;
	if($str=~/\*/)
	{
		$decision=1;
	}
	else
	{
		$decision=0;
	}
	return($decision);
}

sub joinAlignments{
	(my $prevOverlap, my $prevOverlapSimilarity, my $newOverlap, my $newOverlapSimilarity)=@_;
	my @prevGTsPos;
	my %scores;
	(my $prevscore, my $newscore)=processSimilarity($prevOverlapSimilarity,$newOverlapSimilarity);
	@{$scores{$prevOverlap."~"}}=($prevscore,$prevOverlap,"",$prevOverlapSimilarity,"");
	# print "$prevOverlap.'~'\t*@{$scores{$prevOverlap.'~'}}*\n";

	@{$scores{"~".$newOverlap}}=($newscore,"",$newOverlap,"",$newOverlapSimilarity);
	# print "$newOverlap.'~'\t*@{$scores{'~'.$newOverlap}}*\n";

	my $prevOverlapCopy=$prevOverlap;
	while($prevOverlapCopy=~/(.*)GT/i)
	{
    	push(@prevGTsPos,$+[-1]);
        $prevOverlapCopy=$1;
	}
	#push(@prevGTsPos,(length($prevOverlap)-1));#if full GT can not be seen at the end of alignment #not sure if -1 should be there
	# print "@prevGTsPos\n";

	foreach my $i(sort {$a <=> $b} @prevGTsPos)
	{
    	#$i=$prevPos+$i;
        #print "***$i***\n";
    	my $overlapLength=$i-2;
        if($newOverlap=~/^.({$overlapLength,$overlapLength})AG(.*)$/i)
    	{
        	my $ag=$2;
            my $sim=substr($prevOverlapSimilarity,0,$i).substr($newOverlapSimilarity,0-length($ag));
	        ($prevscore, $newscore)=processSimilarity(substr($prevOverlapSimilarity,0,$i),substr($newOverlapSimilarity,0-length($ag)));
        	$prevscore+=5;
            $newscore+=5;
	        my $simScore=$prevscore+$newscore;
        	@{$scores{substr($prevOverlap,0,$i)."~".$ag}}=($simScore,substr($prevOverlap,0,$i),$ag,substr($prevOverlapSimilarity,0,$i),substr($newOverlapSimilarity,length($newOverlapSimilarity)-length($ag)));
    	}
    	else{}
    	($prevscore, $newscore)=(0,0);
        #$prevPos=$i;
	}

	foreach my $j(sort { ${$scores{$b}}[0] <=> ${$scores{$a}}[0] } keys %scores)
	{
        #print "$j\t@{$scores{$j}}\n";
		return(${$scores{$j}}[1],${$scores{$j}}[2],${$scores{$j}}[3],${$scores{$j}}[4]);
		exit;
	}

}

sub joinCDS{
	(my $prevCDS, my $newCDS, my $dir, my $prevJoint, my $newJoint, my $scf)=@_;
	#print "The prevCDS is $prevCDS, new CDS is $newCDS, direction is $dir, prevJoint is $prevJoint, newJoint is $newJoint\n";
    my $finalcds;
    my @oldcds=split("-",$prevCDS);
    my @newcds=split("-",$newCDS);

    if($dir eq "+")
	{
		# if($prevJoint < $oldcds_lengths[-1] && $newJoint < $newcds_lengths[0])
		# {
			my $last=pop @oldcds;
		    $last=$last-$prevJoint;
			#print "Last of prev is $last\n";
			my $first=shift @newcds;
	        $first=$first+$newJoint;
			#print "First of new is $first\n";
	        my $newcd=$first."-".join("-",@newcds);
			my $oldcd=join("-",@oldcds)."-".$last;
			$finalcds=$oldcd."-".$newcd;
		# }
		# else
		# {
		# 	#####what to do???
		# }
	}
    elsif($dir eq "-")
	{
		my $last=pop @newcds;
        $last=$last-$newJoint;
        #print "Last of new is $last\n";
        my $first=shift @oldcds;
        $first=$first+$prevJoint;
        #print "First of old is $first\n";
        my $oldcd=$first."-".join("-",@oldcds);
        my $newcd=join("-",@newcds)."-".$last;
        $finalcds=$newcd."-".$oldcd;
	}
    #print "The final CDS is $finalcds\n";

    my @finalcds_split=split("-",$finalcds);
    my @finalcds_final;

    for(my $i=0; $i < @finalcds_split; $i=$i+2)
    {
    	if($finalcds_split[$i+1]-$finalcds_split[$i] > 0)
    	{
    		push(@finalcds_final,$finalcds_split[$i]);
    		push(@finalcds_final,$finalcds_split[$i+1]);  		
    	}
    	else
    	{
    		print "\nWarning! Gene model overlapping with the following location might be wrong!\nHence the CDS/protein sequence provided for this location may not match with the nucleotides at this location.\n";
    		print "$scf\t$finalcds_split[0]\t$finalcds_split[-1]\n\n";
    	}
    }
   	
   	my $finalcds_string=join("-",@finalcds_final);

    return ($finalcds_string);
}


sub parseHMMTOP{
	my $filename=shift;
	my %hash;
	open(FH1,"<$filename");
	$/="\n";
	while(<FH1>)
	{
		my $k=$_;
		chomp $k;
		my @line="";
		if ($k=~/OUT/){@line=split(/OUT/,$k)}
		elsif ($k=~/IN/){@line=split(/IN/,$k)};
		my @key=split(/\s+/,$line[0]);
		@{$hash{$key[2]}}=split(/\s+/,$line[1]);
		shift @{$hash{$key[2]}};
		shift @{$hash{$key[2]}};
		#print "$key[2]\t@{$hash{$key[2]}}\n";
	}
	close FH1;
	return(\%hash);
}

sub parseTMHMM{
	my $filename=shift;
	my %TMHMM;
	open(FH1,"<$filename");
	$/="\n";
	while(<FH1>)
	{
		my $k=$_;
		chomp $k;
		my @line=split(/\t/,$k);
		@{$TMHMM{$line[0]}}=($line[5]=~/(\d+)/g);
		#print "$line[0]\t@{$hash{$line[0]}}\n";
	}
	close FH1;
	return (\%TMHMM);
}

sub parsePhobius{
	my $filename=shift;
	my %hash;
	open(FH1,"<$filename");
	$/="\n";
	while(<FH1>)
	{
		my $k=$_;
		chomp $k;
		my @tp;
		my @line=split(/\s+/,$k);if($line[2] eq 'TM'){next;}
		if($line[2] eq 'Y')
			{
			if($line[3]=~/(.*?)\/(\d+?)[o|i](.*)/)
				{
				@tp=($3=~/(\d+)/g);
				}
			}
		elsif($line[2]==0)
			{
			@tp=($line[3]=~/(\d+)/g);
			}
		else{print "Error\n"}
		@{$hash{$line[0]}}=@tp;
		#print "$line[0]\t@{$hash{$line[0]}}\n";
	}
	close FH1;
	return(\%hash);
}

sub makeConsensus{
	my $tmhmm_handler=shift;
	my $hmmtop_handler=shift;
	my $phobius_handler=shift;
	my $seqfile=shift;

	my %TMHMM=%$tmhmm_handler;
	my %HMMTOP=%$hmmtop_handler;
	my %PHOBIUS=%$phobius_handler;
	my %CONSENSUS;

	#print "KQ414573.1_628482-653065_OR_Complete_Normal_StartCodonPresent\n".@{$TMHMM{"KQ414573.1_628482-653065_OR_Complete_Normal_StartCodonPresent"}}."\n";

	my %tmhmm_str_hash;
	my %hmmtop_str_hash;
	my %phobius_str_hash;
	my %consensus_str_hash;
	my %consscore_str_hash;

	open(FH1,"<$seqfile");
	$/=">";
	my @sequences=<FH1>;
	close FH1;
	shift @sequences;
	#print "$sequences[0]\n";

	my %seq;
	foreach my $y(@sequences)
	{
		if($y=~/>$/){chop $y;}
		my @sm=split(/\n/,$y);
		my $header=shift @sm;
		#my @important=split(/\s+/,$imp);
		$seq{$header}=join("",@sm);
		#print "$header\n$seq{$header}\n";
	}

	my %hash;
	my %hashstrings;
	my %HMMTOP_temp; my %TMHMM_temp; my %PHOBIUS_temp; my %CONSENSUS_temp;
	#my $raw_scores;

	foreach my $x(sort keys %seq)
	{
		@{$HMMTOP_temp{$x}}=(0) x length($seq{$x});
		@{$TMHMM_temp{$x}}=(0) x length($seq{$x});
		@{$PHOBIUS_temp{$x}}=(0) x length($seq{$x});
		@{$CONSENSUS_temp{$x}}=(0) x length($seq{$x});
		#print "$x\n$seq{$x}\n$HMMTOP_temp{$x}\n";
	}

	foreach my $i(sort keys %seq)
	{
		#$raw_scores=">$i\n$seq{$i}\nHMMTOP\t";
		#$raw_scores.="@{$HMMTOP{$i}}\n";

		my $helix_no=0;
		#$raw_scores.="\nTMHMM\t";
		#$raw_scores.="@{$TMHMM{$i}}\n";
		for(my $j=0; $j<@{$TMHMM{$i}}; $j=$j+2)
		{
			if($helix_no<7){$helix_no++;}
			else{$helix_no="x"}
			for(my $k=${$TMHMM{$i}}[$j]-1; $k<=${$TMHMM{$i}}[$j+1]-1; $k++)
			{
				${$TMHMM_temp{$i}}[$k]=$helix_no;
				${$CONSENSUS_temp{$i}}[$k]+=1;
			}
		}

		#print "$i\n@{$TMHMM{$i}}\n";
		#print "@{$TMHMM_temp{$i}}\n";
		my $tmhmm_str=join("",@{$TMHMM_temp{$i}});
		$tmhmm_str=~tr/0/-/;
		$tmhmm_str_hash{$i}=$tmhmm_str;
		#$raw_scores.=$tmhmm_str;

		$helix_no=0;
		@{$CONSENSUS{$i}}=();
		if(!exists $HMMTOP{$i})
		{
			@{$HMMTOP{$i}}=();
		}
		for(my $j=0; $j<@{$HMMTOP{$i}}; $j=$j+2)
		{
			if($helix_no<7){$helix_no++;}
			else{$helix_no="x"}
			for(my $k=${$HMMTOP{$i}}[$j]-1; $k<=${$HMMTOP{$i}}[$j+1]-1; $k++)
			{
				${$HMMTOP_temp{$i}}[$k]=$helix_no;
				${$CONSENSUS_temp{$i}}[$k]+=1;
			}
		}
		my $hmmtop_str=join("",@{$HMMTOP_temp{$i}});
		$hmmtop_str=~tr/0/-/;
		$hmmtop_str_hash{$i}=$hmmtop_str;
		#$raw_scores.=$hmmtop_str;
		
		$helix_no=0;
		#$raw_scores.="\nPHOBIUS\t";
		#$raw_scores.="@{$PHOBIUS{$i}}\n";
		for(my $j=0; $j<@{$PHOBIUS{$i}}; $j=$j+2)
		{
			if($helix_no<7){$helix_no++;}
			else{$helix_no="x"}
			for(my $k=${$PHOBIUS{$i}}[$j]-1; $k<=${$PHOBIUS{$i}}[$j+1]-1; $k++)
			{
				${$PHOBIUS_temp{$i}}[$k]=$helix_no;
				${$CONSENSUS_temp{$i}}[$k]+=1;
			}
		}

		my $phobius_str=join("",@{$PHOBIUS_temp{$i}});
		$phobius_str=~tr/0/-/;
		$phobius_str_hash{$i}=$phobius_str;
		#$raw_scores.=$phobius_str;
		#my $consensus_raw_scores.="\nConsensusRawScore\n";
		#$raw_scores.=join("",@{$CONSENSUS_temp{$i}});
		#$raw_scores.="\n";
		my $consensus_raw_scores.=join("",@{$CONSENSUS_temp{$i}});
		$consscore_str_hash{$i}=$consensus_raw_scores;
		
		my $start=1; my $stop=0;
		for(my $j=$start; $j<length($seq{$i}); $j++)
		{
			if((${$CONSENSUS_temp{$i}}[$j]>1) && (${$CONSENSUS_temp{$i}}[$j-1]<2))
			{
				$start=$j+1; push(@{$CONSENSUS{$i}},$start);
			}
			elsif((${$CONSENSUS_temp{$i}}[$j-1]>1) && (${$CONSENSUS_temp{$i}}[$j]<2))
			{
				$stop=$j; push(@{$CONSENSUS{$i}},$stop);
			}
			else
			{}
		}

		#$raw_scores.=@{$CONSENSUS{$i}};
		#$raw_scores.="\nConsensusFinalPrediction\t";
		#$raw_scores.="@{$CONSENSUS{$i}}\n";
		
		$helix_no=0;
		@{$CONSENSUS_temp{$i}}=(0) x length($seq{$i});
		for(my $j=0; $j<@{$CONSENSUS{$i}}; $j=$j+2)
		{
			if($helix_no<7){$helix_no++;}
			else{$helix_no="x"}
			for(my $k=${$CONSENSUS{$i}}[$j]-1; $k<=${$CONSENSUS{$i}}[$j+1]-1; $k++)
			{
				${$CONSENSUS_temp{$i}}[$k]=$helix_no;
			}
		}

		my $consensus_str=join("",@{$CONSENSUS_temp{$i}});
		$consensus_str=~tr/0/-/;
		$consensus_str_hash{$i}=$consensus_str;
		#raw_scores.=$consensus_str;
		#$raw_scores.=$consensus_raw_scores;
		#print "$raw_scores\n\n";	
	}

	return(\%tmhmm_str_hash, \%hmmtop_str_hash, \%phobius_str_hash, \%consensus_str_hash, \%consscore_str_hash, \%CONSENSUS, \%seq);
}

sub assignColors{
	my $tmhstr=shift;
	my @colors=("violet","indigo","blue","green","yellow","orange","red","grey");
	$tmhstr=~s/(1{1,})/<font color=$colors[0]>$1<\/font>/g;
	$tmhstr=~s/(2{1,})/<font color=$colors[1]>$1<\/font>/g;
	$tmhstr=~s/(3{1,})/<font color=$colors[2]>$1<\/font>/g;
	$tmhstr=~s/(4{1,})/<font color=$colors[3]>$1<\/font>/g;
	$tmhstr=~s/(5{1,})/<font color=$colors[4]>$1<\/font>/g;
	$tmhstr=~s/(6{1,})/<font color=$colors[5]>$1<\/font>/g;
	$tmhstr=~s/(7{1,})/<font color=$colors[6]>$1<\/font>/g;
	$tmhstr=~s/(x{1,})/<font color=$colors[7]>$1<\/font>/g;
	#$tmhstr=~s/([8|9|10|11|12|13|14|15]{1,})/<font color=$colors[7]>$1<\/font>/;
	return($tmhstr);
}

sub barplot{
	#my $xaxis_handler=shift;
	my $tmhmm_handler=shift;
	my $hmmtop_handler=shift;
	my $phobius_handler=shift;
	my $consensus_handler=shift;
	#my %xaxis=%$xaxis_handler;
	my %tmhmm=%$tmhmm_handler;
	my %hmmtop=%$hmmtop_handler;
	my %phobius=%$phobius_handler;
	my %consensus=%$consensus_handler;

	my %combined;

	foreach(sort keys %tmhmm)
	{
		${$combined{$_}}[0]=$tmhmm{$_};
	}

	foreach(sort keys %hmmtop)
	{
		${$combined{$_}}[1]=$hmmtop{$_};
	}

	foreach(sort keys %phobius)
	{
		${$combined{$_}}[2]=$phobius{$_};
	}

	foreach(sort keys %consensus)
	{
		${$combined{$_}}[3]=$consensus{$_};
	}

	foreach my $j(sort keys %combined)
	{
		foreach my $i(0..3)
		{
			if( !defined ${$combined{$j}}[$i] )
			{
				${$combined{$j}}[$i]=0;
			}
		}
	}

	my @tmhmm_arr; my @hmmtop_arr; my @phobius_arr; my @consensus_arr;
	foreach my $i(sort keys %combined)
	{
		push(@tmhmm_arr,${$combined{$i}}[0]);
		push(@hmmtop_arr,${$combined{$i}}[1]);
		push(@phobius_arr,${$combined{$i}}[2]);
		push(@consensus_arr,${$combined{$i}}[3]);
	}

	my $data = GD::Graph::Data->new([
	    #["1st","2nd","3rd","4th","5th","6th","7th", "8th", "9th"],
	    [(sort keys %combined)],
	    [@tmhmm_arr],
	    [@hmmtop_arr],
	    [@phobius_arr],
	    [@consensus_arr],
	]) or die GD::Graph::Data->error;
	 
	 
	my $graph = GD::Graph::bars->new(800, 500 );
	 
	$graph->set( 
	    x_label         => 'No. of transmembrane helices per sequence',
	    y_label         => 'No. of ORs',
	    title           => 'TMH predictions',
	 	x_label_position => 0.5,
	    #y_max_value     => 7,
	    #y_tick_number   => 8,
	    #y_label_skip    => 3,
	 
	    #x_labels_vertical => 1,
	 	#bar_width => 10,
	    bar_spacing     => 7,
	    #shadow_depth    => 4,
	    #shadowclr       => 'dred',
	    transparent     => 0,
	) or die $graph->error;
	$graph->set_legend( qw(TMHMM HMMTOP Phobius Consensus));
	#$graph->set_legend_font('/fonts/arial.ttf', 26);
	#$graph->set_x_label_font('/fonts/arial.ttf', 26);
	#$graph->set_y_label_font('/fonts/arial.ttf', 26);
	$graph->plot($data) or die $graph->error;
	
	return($graph);
	
}

#0:Scaffold_name	1:Cluster_start	2:Cluster_end	3:ScoreBestStart	4:ScoreBestEnd	5:ScoreBestScore	6:ScoreBestLength	7:ScoreBestQstart	8:ScoreBestQend	
#9:ScoreBasedDirection	10:ScoreBasedQuery	11:FullEntry	12:CoverageBestStart	13:CoverageBestEnd	14:CoverageBestScore	15:CoverageBestLength	
#16:CoverageBestQstart	17:CoverageBestQend	18:CoverageBasedDirection	19:CoverageBasedQuery	20:FullEntry	21:ClusterCoverage	22:HitOverlapBtwnMethods	
#23:Decision	24:FinalMaxLength	25:FullORPartial

sub processFinalNew{
	#####Stitching alignments if same scaffold same query and overlapping hits#####
	print "Retrieving all alignments for the best query and scaffold pair...\n\n";
	print "Stitching the alignments with overlapping query region in the correct direction...\n\n";
    my $has=shift;
	my %hash=%$has;
	my $hash_sup=shift;
    my $filename=shift;
    my $lengthCutoff=shift;
	my %hash_suppl=%$hash_sup;
    #my %new_hash;
    #my %new_hash_suppl;
    my $prev_scaf=0;
    my $prev_key;
    my $prev_SQstart=0; my $prev_CQstart=0;
    my $prev_SQend=0; my $prev_CQend=0;
    my $prev_Sdir=0; my $prev_Cdir=0;
    my $prev_clust;
    my $prev_Squery=0; my $prev_Cquery=0;
    ###my $prev_score=0;###
	my @ali_details;
	my @prev_ali_details;

	
    # make sure we run through the hash in a sorted fashion
    foreach my $v(sort {$a <=> $b} keys %hash)
    {
    	#print "Current key is $v\t\n";
        #print "${$hash{$v}}[0]\t${$hash{$v}}[1]\t${$hash{$v}}[2]\n";
        #Initialize prev start and end of protein.
        my $curr_clust="@{$hash{$v}}[0]~@{$hash{$v}}[1]~@{$hash{$v}}[2]";
        my $curr_scaf=@{$hash{$v}}[0];
        my $curr_SQstart=@{$hash{$v}}[7];
        my $curr_SQend=@{$hash{$v}}[8];
        my $curr_Sdir=@{$hash{$v}}[9];
        my $curr_Squery=@{$hash{$v}}[10];
        my $curr_Slength=@{$hash{$v}}[6];
        ###my $curr_score=@{$hash{$v}}[5];###
        #my $curr_CQstart=@{$hash{$v}}[16];
        #my $curr_CQend=@{$hash{$v}}[17];
        #my $curr_Cdir=@{$hash{$v}}[18];
        #my $curr_Cquery=@{$hash{$v}}[19];
        #my $curr_Clength=@{$hash{$v}}[15];
		my %alignscores;
		my $finalCDS;
		my $finalProt;
		my $finalTranscript;
		my $finalSimilarity;
		my $flag=0;
		#print "prevScaf-$prev_scaf currScaf-$curr_scaf prevSDir-$prev_Sdir currSDir-$curr_Sdir currSLength-$curr_Slength prevCDir-$prev_Cdir currCDir-$curr_Cdir currCLength-$curr_Clength\n";
		#in following section only score based hits are processed
        if(($prev_scaf eq $curr_scaf) && ($prev_Squery eq $curr_Squery) && ($prev_Sdir eq $curr_Sdir) && ($curr_Slength<420))
        {
			my $temp_curr_SQstart; my $temp_curr_SQend ;my $temp_prev_SQend; my $temp_prev_SQstart;
            if($curr_Sdir eq "+" )
			{
				#print "The hits are in the positive direction\n";
				$temp_curr_SQstart=$curr_SQstart; $temp_curr_SQend=$curr_SQend; $temp_prev_SQend=$prev_SQend ; $temp_prev_SQstart=$prev_SQstart;
				#@ali_details=processAlignment(@{$hash->{$v}}[11]);
	            #@prev_ali_details=processAlignment(@{$hash->{$prev_key}}[11]);
				#$queryProt,$targetNuc,$similarity,$targetProt,$cds
				@ali_details=($hash_suppl{$v}{"QueryProt"},$hash_suppl{$v}{"TargetNuc"},$hash_suppl{$v}{"Similarity"},$hash_suppl{$v}{"TargetProt"},$hash_suppl{$v}{"CDS"});
				@prev_ali_details=($hash_suppl{$prev_key}{"QueryProt"},$hash_suppl{$prev_key}{"TargetNuc"},$hash_suppl{$prev_key}{"Similarity"},$hash_suppl{$prev_key}{"TargetProt"},$hash_suppl{$prev_key}{"CDS"});
			}
			elsif($curr_Sdir eq "-")
			{
				#print "The hits are in the negative direction\n";
				$temp_curr_SQstart=$prev_SQstart; $temp_curr_SQend=$prev_SQend; $temp_prev_SQend=$curr_SQend ; $temp_prev_SQstart=$curr_SQstart;
				#@prev_ali_details=processAlignment(@{$hash->{$v}}[11]);
	            #@ali_details=processAlignment(@{$hash->{$prev_key}}[11]);
				@prev_ali_details=($hash_suppl{$v}{"QueryProt"},$hash_suppl{$v}{"TargetNuc"},$hash_suppl{$v}{"Similarity"},$hash_suppl{$v}{"TargetProt"},$hash_suppl{$v}{"CDS"});
	            @ali_details=($hash_suppl{$prev_key}{"QueryProt"},$hash_suppl{$prev_key}{"TargetNuc"},$hash_suppl{$prev_key}{"Similarity"},$hash_suppl{$prev_key}{"TargetProt"},$hash_suppl{$prev_key}{"CDS"});
			}
			else
			{
				print "The direction is ambiguous!!!!\n"; 
			}

			

			my @prevcds_tmp=split("-",$prev_ali_details[4]);
			my @currcds_tmp=split("-",$ali_details[4]);



			if ($temp_curr_SQstart<=$temp_prev_SQend && $temp_curr_SQstart>=$temp_prev_SQstart && $temp_curr_SQend>$temp_prev_SQend)
            {
            	#handle complete overlap case
            	#print "";
				####Need to change things in here
                #print "For @{$hash->{$v}}[0], @{$hash->{$prev_key}}[1]~@{$hash->{$prev_key}}[2] and @{$hash->{$v}}[1]~@{$hash->{$v}}[2] can be combined in that order according to best scoring hit\n";
				my $overlap=3*($temp_prev_SQend-$temp_curr_SQstart+1);
				#my $joiningSeqLength=length(substr($prev_ali_details[1],0,0-$overlap).$ali_details[1]);
				#my $joiningSeqLength=length($prev_ali_details[1])+length($ali_details[1])-$overlap;
				# print "Overlap is $overlap\nJoining seq length is $joiningSeqLength\n";
				#include overlap criterian here to stop joining two separate partial genes
			################# User provided input to distinguish between two different isoforms??? Instead pf 20 some other percentage of overlap could be given#############################################				
				if( (($overlap / length($prev_ali_details[1]))*100) < 20 && (($overlap / length($ali_details[1]))*100) <20 )
				#if((($overlap/$joiningSeqLength)*100) < 20)
				{
					#print "$prev_ali_details[1]\n";
					#print "Overlap is $overlap\n";
					#print "$prev_ali_details[2]\n";
					#print "$ali_details[1]\n";
					#print "$ali_details[2]\n\n";
					
					if( ($curr_Sdir eq "+" && int($prevcds_tmp[-1]) >= int($currcds_tmp[0])) || ( $curr_Sdir eq "-" && int($prevcds_tmp[0]) <= int($currcds_tmp[-1])) )
					{
						#print "Direction is $curr_Sdir\n";
						#print "Prev ali CDS $prev_ali_details[4]\n";
						#print "Current ali CDS $ali_details[4]\n";
						#print "$prevcds_tmp[-1]\t$currcds_tmp[0]\n";
						#print "$prevcds_tmp[0]\t$currcds_tmp[-1]\nLengths of prev and curr are - ";
						#print length($prev_ali_details[1])."\t".length($ali_details[1])."\n";
						
						if(length($prev_ali_details[1]) > length($ali_details[1]))
						{
							# print "Deleting current key\n";
							#$v=$prev_key;
							#delete $hash{$v};
							#delete $hash_suppl{$v};

							$hash{$v}=$hash{$prev_key};
							$hash_suppl{$v}=$hash_suppl{$prev_key};
							delete $hash{$prev_key};
							delete $hash_suppl{$prev_key};
							#print "Current length is \t";
							#print length($hash_suppl{$v}{"TargetNuc"});
							#print "\n";
						}
						else
						{
							# print "Deleting prev key\n";
							delete $hash{$prev_key};
							delete $hash_suppl{$prev_key};
						}

						#TODO remove either current or prev depending on score
						#delete $hash{$prev_key};
						#delete $hash_suppl{$prev_key};
					}
					else
					{
						my $prev_substr=substr($prev_ali_details[2],0-$overlap);

						#my $prev_substr_sim=substr($prev_ali_details[2],0-$overlap);
						my $substr=substr($ali_details[2],0,$overlap);
						#my $substr_sim=substr($ali_details[2],0,$overlap);
						#print "~$prev_substr~\t~$substr~\n";

						my @prev_hashes_array=($prev_substr=~/(#{1,})/g);
						my $prev_hashes=length(join("",@prev_hashes_array));

						my @hashes_array=($substr=~/(#{1,})/g);
						my $hashes=length(join("",@hashes_array));

						#print "# $prev_hashes #\t# $hashes #\n";

						my @overlapTargetNuc=joinAlignments(substr($prev_ali_details[1],(0-$overlap-$prev_hashes)),substr($prev_ali_details[2],(0-$overlap-$prev_hashes)),substr($ali_details[1],0,($overlap+$hashes)),substr($ali_details[2],0,($overlap+$hashes)));
						#print "Array overlapTargetNuc is @overlapTargetNuc\n";
						#print "PrevQuerySequence and newQuerysequence is\n$prev_ali_details[0] $ali_details[0]\n";
						#print "PrevCDSsequences and newCDSseq is\n$prev_ali_details[1] $ali_details[1]\n";
						$finalTranscript=substr($prev_ali_details[1],0,(0-$overlap-$prev_hashes)).$overlapTargetNuc[0].$overlapTargetNuc[1].substr($ali_details[1],($overlap+$hashes));
						#print "Final CDS Nucleotide sequence is\n$finalTranscript\n";
						#print "Final protein sequence is\n";
						my $str1=substr($prev_ali_details[3],0,(0-$overlap-$prev_hashes));
						my $str3=substr($ali_details[3],($overlap+$hashes));
						#print "$overlapTargetNuc[0].$overlapTargetNuc[1]\n";
						my $str2=translateNucToProt($overlapTargetNuc[0].$overlapTargetNuc[1]);
						$str2=one2three($str2);
						$finalProt=$str1.$str2.$str3;
						$finalSimilarity=substr($prev_ali_details[2],0,(0-$overlap-$prev_hashes)).$overlapTargetNuc[2].$overlapTargetNuc[3].substr($ali_details[2],($overlap+$hashes));
						# print "$finalProt\n";
						#translateNucToProt("ATGAGTAGTTTCACCACCGATGACATATCAATCAGCCTGACGTCCGTTTTCATGAAG");
						#print "PrevCDS is $prev_ali_details[4]\n"; print "CDS is $ali_details[4]\n";
						my @prevCDSs=split("-",$prev_ali_details[4]); my @newCDSs=split("-",$ali_details[4]);						
						$finalCDS=joinCDS($prev_ali_details[4],$ali_details[4],$curr_Sdir,($overlap-length($overlapTargetNuc[0])),($overlap-length($overlapTargetNuc[1])),$curr_scaf);
						${$hash{$v}}[1]=${$hash{$prev_key}}[1];
						${$hash{$v}}[3]=${$hash{$prev_key}}[3];
						${$hash{$v}}[5]=${$hash{$v}}[5]+${$hash{$prev_key}}[5];
						${$hash{$v}}[6]=length($finalProt);
						if(${$hash{$v}}[9] eq "+"){${$hash{$v}}[7]=${$hash{$prev_key}}[7];}
						else{${$hash{$v}}[8]=${$hash{$prev_key}}[8];}
						${$hash{$v}}[11]=${$hash{$prev_key}}[11].${$hash{$v}}[11];
						#@{$new_hash{$v}}=@{$hash{$v}};
						$hash_suppl{$v}{"TargetNuc"}=$finalTranscript;
						$hash_suppl{$v}{"Similarity"}=$finalSimilarity;
						$hash_suppl{$v}{"TargetProt"}=$finalProt;
						$hash_suppl{$v}{"CDS"}=$finalCDS;
						#$new_hash_suppl{$v}=$hash_suppl{$v};
						#print "Final CDS is $finalCDS\n";
						delete $hash{$prev_key};
						delete $hash_suppl{$prev_key};
						#print"\n\n";
					}
				}
				else
				{
					#print "The fragments should be kept as separate as overlap between them is too high\n";
                }
			}
			elsif($temp_curr_SQstart==$temp_prev_SQend+1)
			{
				if( ($curr_Sdir eq "+" && int($prevcds_tmp[-1]) >= int($currcds_tmp[0])) || ( $curr_Sdir eq "-" && int($prevcds_tmp[0]) <= int($currcds_tmp[-1])) )
				{
					# print "Direction is $curr_Sdir\n";
					# print "Prev ali CDS $prev_ali_details[4]\n";
					# print "Current ali CDS $ali_details[4]\n";
					# print "$prevcds_tmp[-1]\t$currcds_tmp[0]\n";
					# print "$prevcds_tmp[0]\t$currcds_tmp[-1]\n";
					
					if(length($prev_ali_details[1]) > length($ali_details[1]))
					{
						# print "Deleting current key when equal\n";
						$hash{$v}=$hash{$prev_key};
						$hash_suppl{$v}=$hash_suppl{$prev_key};
						delete $hash{$prev_key};
						delete $hash_suppl{$prev_key};
					}
					else
					{
						# print "Deleting prev key\n";
						delete $hash{$prev_key};
						delete $hash_suppl{$prev_key};
					}
					
					#TODO remove either current or prev depending on score
					#delete $hash{$prev_key};
					#delete $hash_suppl{$prev_key};
				}
				else
				{
					#handle negative direction here
					#print "PrevQuerySequence and newQuerysequence is\n$prev_ali_details[0] $ali_details[0]\n";
					#print "PrevCDSsequences and newCDSseq is\n$prev_ali_details[1] $ali_details[1]\n";
					$finalTranscript=$prev_ali_details[1].$ali_details[1];
					#print "Final CDS Nucleotide sequence is\n$finalTranscript\n";
					$finalProt=$prev_ali_details[3].$ali_details[3];
					#print "Final protein sequence is\n$finalProt\n";
					$finalSimilarity=$prev_ali_details[2].$ali_details[2];
					${$hash{$v}}[1]=${$hash{$prev_key}}[1];
					${$hash{$v}}[3]=${$hash{$prev_key}}[3];
					${$hash{$v}}[5]=${$hash{$v}}[5]+${$hash{$prev_key}}[5];
					${$hash{$v}}[6]=length($finalProt);
					if(${$hash{$v}}[9] eq "+")
					{
						${$hash{$v}}[7]=${$hash{$prev_key}}[7];
						$finalCDS=$prev_ali_details[4]."-".$ali_details[4];
					}
					else
					{
						${$hash{$v}}[8]=${$hash{$prev_key}}[8];
						$finalCDS=$ali_details[4]."-".$prev_ali_details[4];
					}
					#print "PrevCDS is $prev_ali_details[4]\n"; print "CDS is $ali_details[4]\n";
					#print "The final CDS is $finalCDS\n";
					${$hash{$v}}[11]=${$hash{$prev_key}}[11].${$hash{$v}}[11];
					$hash_suppl{$v}{"TargetNuc"}=$finalTranscript;
					$hash_suppl{$v}{"Similarity"}=$finalSimilarity;
					$hash_suppl{$v}{"TargetProt"}=$finalProt;
					$hash_suppl{$v}{"CDS"}=$finalCDS;
					delete $hash{$prev_key};
					delete $hash_suppl{$prev_key};
					#print"\n\n";
				}
			}
			else
			{

			}
        }
        else
        {
            #print "prevScaf-$prev_scaf currScaf-$curr_scaf prevDir-$prev_dir currDir-$curr_dir currLength-$curr_length Hence skipped\n";
			#$flag=1;
        }


	    $prev_SQstart=$curr_SQstart; $prev_SQend=$curr_SQend; $prev_Sdir=$curr_Sdir; $prev_Squery=$curr_Squery;
	    #in following section only coverage based hits are processed
	    #$prev_CQstart=$curr_CQstart; $prev_CQend=$curr_CQend; $prev_Cdir=$curr_Cdir; $prev_Cquery=$curr_Cquery;
		$prev_key=$v; $prev_scaf=$curr_scaf; $prev_clust=$curr_clust;
	}



	#print "\n\n***************************************************************************************************************************************\n\n";
	print "Stitching finished.\n\n";


	#####Out of overlapping hits - select best scoring one#####
	foreach my $v(sort {$a <=> $b} keys %hash_suppl)
	{
		my $ORname="${$hash{$v}}[0]_${$hash{$v}}[3]-${$hash{$v}}[4]_OR";
		my @cdsboundary=split("-",$hash_suppl{$v}{"CDS"});

		#call three2one before declring final targetProt
		$hash_suppl{$v}{"TargetProt"}=three2one($hash_suppl{$v}{"TargetProt"});
		$hash_suppl{$v}{"TargetNuc"}=~s/---//g;
		$hash_suppl{$v}{"TargetProt"}=~s/\*\*/\*/g;
		$hash_suppl{$v}{"Scaf"}=${$hash{$v}}[0]; $hash_suppl{$v}{"ClustStart"}=${$hash{$v}}[1]; $hash_suppl{$v}{"ClustEnd"}=${$hash{$v}}[2];
		$hash_suppl{$v}{"GeneStart"}=$cdsboundary[0]; $hash_suppl{$v}{"GeneEnd"}=$cdsboundary[-1];
		$hash_suppl{$v}{"Query"}=${$hash{$v}}[10]; $hash_suppl{$v}{"QueryStart"}=${$hash{$v}}[7]; $hash_suppl{$v}{"QueryEnd"}=${$hash{$v}}[8]; $hash_suppl{$v}{"QueryDir"}=${$hash{$v}}[9];
		#if(length($hash_suppl{$v}{"TargetNuc"}) != 3*length($hash_suppl{$v}{"TargetProt"}))
		#{
		#	print "$hash{$v}}[0]_${$hash{$v}}[3]-${$hash{$v}}[4]_OR has some problem\n";
		#}
		

		# if (length($hash_suppl{$v}{"TargetProt"}) > 370 )
		# {
		# 	$ORname.="_Complete";
		# 	#print OUT4 "Complete\t";
		# 	$hash_suppl{$v}{"Completeness"}=3;
		# }
		# elsif(length($hash_suppl{$v}{"TargetProt"}) > 350)
		# {
		# 	$ORname.="_CompleteAlmost";
		# 	#print OUT4 "CompleteAlmost\t";
		# 	$hash_suppl{$v}{"Completeness"}=2;
		# }
		if (length($hash_suppl{$v}{"TargetProt"}) >= $lengthCutoff )
		{
			$ORname.="_Complete";
			#print OUT4 "Complete\t";
			$hash_suppl{$v}{"Completeness"}=2;
		}
		else
		{
			$ORname.="_Partial";
			#print OUT4 "Partial\t";
			$hash_suppl{$v}{"Completeness"}=1;
		}

		my $pseudoDecision=pseudoPresence($hash_suppl{$v}{"TargetProt"});
		if($pseudoDecision==1)
		{
			$ORname.="_Pseudo";
			#print OUT4 "Pseudo\t";
			$hash_suppl{$v}{"PseudoNature"}=1;
		}
		else
		{
			$ORname.="_Normal";
			#print OUT4 "Normal\t";
			$hash_suppl{$v}{"PseudoNature"}=2;
		}

		my $decision=startCodonPresence($hash_suppl{$v}{"TargetProt"});
		if($decision==1)
		{
			$ORname.="_StartCodonPresent";
			#print OUT4 "StartCodonPresent\n";
			$hash_suppl{$v}{"StartCodon"}=2;

		}
		else
		{
			$ORname.="_StartCodonAbsent";
			#print OUT4 "StartCodonAbsent\n";
			$hash_suppl{$v}{"StartCodon"}=1;
		}

		$hash_suppl{$v}{"Name"}=$ORname;
	}
	
	print "Retaining best hits only...\n\n";

	my %refined_hash; my %refined_hash_suppl; my $prev_w;
	foreach my $w(sort {$a <=> $b} keys %hash_suppl)
	{
		if(exists $hash_suppl{$w}){}
		else{next;}
		my $wdel=0;
		foreach my $prev_w(sort {$a <=> $b} keys %hash_suppl)
		{
			if(exists $hash_suppl{$prev_w}){}
			else{next;}

			if($w == $prev_w){next;}
			my $dec;
			if($hash_suppl{$w}{"Scaf"} eq $hash_suppl{$prev_w}{"Scaf"})
			{
				if(($hash_suppl{$w}{"GeneStart"} >= $hash_suppl{$prev_w}{"GeneStart"}) && ($hash_suppl{$w}{"GeneStart"} < $hash_suppl{$prev_w}{"GeneEnd"}))
				{
					$dec=1
				}
				elsif(($hash_suppl{$w}{"GeneEnd"} <= $hash_suppl{$prev_w}{"GeneEnd"}) && ($hash_suppl{$w}{"GeneEnd"} > $hash_suppl{$prev_w}{"GeneStart"}))
				{
					$dec=1;
				}
				elsif(($hash_suppl{$w}{"GeneStart"} == $hash_suppl{$prev_w}{"GeneStart"}) && ($hash_suppl{$w}{"GeneEnd"} == $hash_suppl{$prev_w}{"GeneEnd"}))
				{
					$dec=1;
				}
				else
				{
					$dec=0;
				}
			}
			else{next;}
	

			if($dec==1)
			{
			#can include a check of pfam signature
				if($hash_suppl{$w}{"Completeness"} == $hash_suppl{$prev_w}{"Completeness"})
				{
					if(length($hash_suppl{$w}{"TargetProt"}) == length($hash_suppl{$prev_w}{"TargetProt"}))
					{
						if($hash_suppl{$w}{"PseudoNature"} == $hash_suppl{$prev_w}{"PseudoNature"})
						{
							if($hash_suppl{$w}{"StartCodon"} == $hash_suppl{$prev_w}{"StartCodon"})
							{
								delete $hash_suppl{$prev_w}; delete $hash{$prev_w};
							}
							elsif($hash_suppl{$w}{"StartCodon"} > $hash_suppl{$prev_w}{"StartCodon"})
							{
								delete $hash_suppl{$prev_w}; delete $hash{$prev_w};
							}
							else
							{
								$wdel=1; #delete $hash_suppl{$w}; delete $hash{$w}; last;
							}
						}
						elsif( $hash_suppl{$w}{"PseudoNature"} > $hash_suppl{$prev_w}{"PseudoNature"} )
						{
							delete $hash_suppl{$prev_w}; delete $hash{$prev_w};
						}
						else
						{
							$wdel=1; #delete $hash_suppl{$w}; delete $hash{$w}; last;
						}
					}
					elsif(length($hash_suppl{$w}{"TargetProt"}) > length($hash_suppl{$prev_w}{"TargetProt"}))
					{
						delete $hash_suppl{$prev_w}; delete $hash{$prev_w};
					}
					else
					{
						$wdel=1; #delete $hash_suppl{$w}; delete $hash{$w}; last;
					}
				}
				elsif($hash_suppl{$w}{"Completeness"} > $hash_suppl{$prev_w}{"Completeness"})
				{
					delete $hash_suppl{$prev_w}; delete $hash{$prev_w};
				}
				else
				{
					$wdel=1; #delete $hash_suppl{$w}; delete $hash{$w}; last;
				}
			}
		}

		if($wdel==1)
		{
			delete $hash_suppl{$w}; delete $hash{$w};
		}

	}

	#Currently the phase of each each gene is not available - but can be calculated. So as to compare within a cluster like older papers used to do. 
	#For this need to calculate length of CDS covered by the previous exons - divide by 3 - see what is the remainder. 3-remainder should be the phase.
	print "Best hits decided.\n\n*********************************************************************************************************\nCore algorithm to find out best possible OR gene regions from the genome is over.\n*********************************************************************************************************\n\nGeneWise gene prediction will be performed next...\n\n";
	return(\%hash,\%hash_suppl)
}

sub processEveryAlignment{
	print "Choosing best queries for each alignment cluster...\n\n";
	my $filename=shift;
	my $cutoff=shift;
	my $vector_length=shift;
	my $genRawFile=shift;
	my $has=shift;
	my $geneloc=shift;
	my $scafoff=shift;
	my $scaflen=shift;
	my %scaf_genes=%$has;
	my %gene_loci=%$geneloc;
	my %scaf_offset=%$scafoff;
	my %scaf_length=%$scaflen;
	my $key_counter=1;
	my %final_hits;

	open(OUT2,">$filename.geneclusters.$cutoff.txt");
	print OUT2 "Scaffold_name\tCluster_start\tCluster_end\tScoreBestStart\tScoreBestEnd\tScoreBestScore\tScoreBestLength\tScoreBestQstart\tScoreBestQend\tScoreBasedDirection\tScoreBasedQuery\tFullEntry\n";
	my %scaf_scorebestquery_hits;
	foreach my $m(sort keys %scaf_genes)
	{
		#print "scoring for scaffold $m ...\n";
		#print "scoring : ".time()."\n";
		my $out_stop = @{$scaf_genes{$m}};
		my $xyz_ref=\@{$gene_loci{$m}};
		my $abc_ref=\@{$scaf_genes{$m}};
		####Incrementing score by one for each alignment in following for loop for one scaffold at a time####
		for(my $n=0;$n<$out_stop;$n=$n+$vector_length)
		{
			my $start=${$scaf_genes{$m}}[$n]-1-$scaf_offset{$m}+1;
			my $stop = ${$scaf_genes{$m}}[$n+1]-$scaf_offset{$m}+1;
			#print "Incrementing score from $start to $stop\n";
			@{$gene_loci{$m}}[$start..($stop-1)]=map {$_+=1}@{$gene_loci{$m}}[$start..($stop-1)];
		}

		#print scalar(@{$gene_loci{$m}});
		#print OUT2"$m\n";
	    my $r=0;
	    
	    ####Identifying alignment clusters based on cutoff value defined by user $cutoff####
		while(1)
		{
			#print "cluster boundary : ".time()."\n";
			my $max_scaf_length=$scaf_length{$m}-$scaf_offset{$m}+1;
			if ($r >= $max_scaf_length - 1)
			{last;}
			my $start_boundary;
			if ($r == 0 && ${$xyz_ref}[0] > $cutoff)
			{
				$start_boundary = 0
			}
			else 
			{
				while(${$xyz_ref}[$r] <= $cutoff ) 
				{
					
					$r++;
					if ($r >= $max_scaf_length - 1) {last;}
				}
				$start_boundary = $r;
		    }

	        #accumulate/run through regions belonging to cluster boundaries(these are above cutoff)
			while(${$xyz_ref}[$r] > $cutoff ) 
			{
				$r++;
				if ($r >= $max_scaf_length -1 ) {last;}
			}
			#loop exits when cutoff is broken, either a cluster boundary is deducted or no cluster region

			$r = $r+$scaf_offset{$m}-1;
			$start_boundary += $scaf_offset{$m}-1;
			#print "Cluster start : $start_boundary, end : $r\n";
			if ($r > $start_boundary)
			{
			# a cluster is found where  $start_boundary+1 is start of cluster and $r is end of cluster
	            #print $start_boundary+1,"\t",$r,"\t";
				#print "\n";
				print OUT2"$m\t";
				print OUT2 $start_boundary+1,"\t",$r,"\t";
				#my $final_hits_key=sprintf("%s_%d-%d",$m,$start_boundary+1,$r);
				my $final_hits_key=$key_counter;
				
				push(@{$final_hits{$final_hits_key}},$m,$start_boundary+1,$r);	
				my $confidence = 0;
	            my $high_conf_hit_pos = 0;
				
				#looking for the best scoring hit overlapping the start boundary and end ($r)

				#print "Best hit : ".time()."\n";
				for(my $n=0;$n<$out_stop;$n=$n+$vector_length)
				{
					my $hit_start = ${$abc_ref}[$n];
					my $hit_end = ${$abc_ref}[$n+1];
					my $overlapFound = 0;
					
					if($hit_start==($start_boundary+1) || $hit_end==$r || ($hit_start>($start_boundary+1) && $hit_start<$r) || ($hit_end>($start_boundary+1) && $hit_end<$r) || ($hit_start<($start_boundary+1) && $hit_end>($start_boundary+1)) || ($hit_start<$r && $hit_end>$r) )
					{
						$overlapFound = 1;
						my $hit_confidence = ${$abc_ref}[$n+2];
						if ($hit_confidence > $confidence)
						{
							$confidence = $hit_confidence;
							$high_conf_hit_pos = $n;
						}

					}
				}
	            print OUT2 "${$abc_ref}[$high_conf_hit_pos]\t${$abc_ref}[$high_conf_hit_pos+1]\t${$abc_ref}[$high_conf_hit_pos+2]\t${$abc_ref}[$high_conf_hit_pos+3]\t${$abc_ref}[$high_conf_hit_pos+5]\t${$abc_ref}[$high_conf_hit_pos+6]\t${$abc_ref}[$high_conf_hit_pos+7]\t${$abc_ref}[$high_conf_hit_pos+8]\t-\n";
				##Hash of two keys - scaffold and best query -
				$scaf_scorebestquery_hits{$m}{${$abc_ref}[$high_conf_hit_pos+8]}=1;
				$key_counter++;
			}
			$r = $r - $scaf_offset{$m} + 1;
			$r++;
			#print "Cluster end : ".time()."\n";
		}	

		if($genRawFile==1)
		{
			open(OUT1,">$filename.$m.genes_new.txt");
			my $counter=1;
			foreach my $q(@{$gene_loci{$m}})
			{
				print OUT1"$counter	$q\n";
				$counter++;
			}
			close OUT1;
		}
	#last;
	}
	close OUT2;
	#Returning only high scoring hit and related values ##Hash of two keys - scaffold and best query - is returned
	return (\%scaf_scorebestquery_hits);
}

sub findMet{
	my @temp=@_;
	my $flg;
	if($_[0]=~/^M/)
	{
		$flg="StartCodonPresent";
	}
	else
	{
		$flg="StartCodonAbsent";
	}
	return ($flg);
}
