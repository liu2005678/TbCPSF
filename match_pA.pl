#!/usr/bin/perl
use String::Approx 'amatch';
$fil= shift @ARGV;
$fil=~s/.sam//;
$file= $fil.".sam";
$file2=$fil."_perl.sam";
$file3=$fil."_notAT.sam";
open File, $file;
@sequences=<File>;
close(File);
for($i=0;$i<=$#sequences;$i++){
	if($sequences[$i]=~/^@/) {push @sam, $sequences[$i];} 
	else{
		my @temp=split("\t", $sequences[$i]);
		my $readname=$temp[0];
		my $cigar=$temp[5];
		#print $cigar,"\n";
		my $seq=$temp[9];
		if($cigar=~/^(\d+)S/ and $1 >= 8 and $1 <= 150){
			$soft=substr $seq, 0, $1+2;
			if(amatch("TTTTTTTTTT", $soft)){ 
				#print $readname,"\tT", "\n",$soft,"\n"; 
				push @sam, $sequences[$i];
			}
			elsif($cigar=~/(\d+)S$/ and $1 >= 8 and $1 <= 150){
				$soft=substr $seq, -($1+2);
				if(amatch("AAAAAAAAAA", $soft)){ #print $readname,"\t","A","\n",$soft, "\n";
												 push @sam, $sequences[$i];} 
				else{push @sam2, $sequences[$i];}
			}
			else{push @sam2, $sequences[$i];}
		}elsif($cigar=~/(\d+)S$/ and $1 >= 8 and $1 <= 150){
				$soft=substr $seq, -($1+2);
				if(amatch("AAAAAAAAAA", $soft)){ #print $readname,"\t","A","\n",$soft, "\n";
												push @sam, $sequences[$i];}
				else {push @sam2,$sequences[$i];}
		}else  {push @sam2,$sequences[$i];}
		 
	}
}

open File, ">$file2";
print File  @sam;
close(File);
open File, ">$file3";
print File  @sam2;
close(File);
