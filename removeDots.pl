#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
# makeSequenceACTG.pl
# Holly Arnold
# PURPOSE:	To remove any characters in an alignment file that are "." and replace with a "-".
# INPUTS: 	-i: This is an input alignment file where each sequence identifier starts with a >
# OUTPUTS:	-o: An alignment file where all "." letters have been replaced with a gap 


my $in_seqs;
my $out_seqs;

GetOptions(
	"i:s" => \$in_seqs,
	"o:s" => \$out_seqs,
);

open INSEQS, $in_seqs or die ("I cant find your sequences file. Please specify another -i option.");

open OUT, ">$out_seqs" or die "I cant open outfile, please check -o option. \n";

print("Reading in $in_seqs file.\n");

my $lineNumber = 0;
my $charCount;
while(<INSEQS>){
	my $line = $_;
	$lineNumber = $lineNumber + 1;
	if ($lineNumber%10000000 == 0){
		print "Processing line number: $lineNumber\n";
	}
	if($_ =~/^>/){ # if we are at the head of a sequence, print ID
		#print("Processing sequence: $_\n");
		chomp $line;
		my @array = split (/\t/, $line);
		print OUT "\n", $array[0], "\n";
		#print OUT $line;
		$charCount = 0;
	}else{ # we are at a sequence line
		
		for my $c (split //, $line){
			#print $c, "\n";
			#print $charCount, "\n";
			if ($charCount == 60){
				$charCount = 0;
				print OUT "\n";
				#print $charCount, "\n";
			}
			if ($c =~ /[.]/){ #the character matches an odd character
				chomp;
				#print ("Found $c character on line $lineNumber\n");
				$c = "-";
				$charCount = $charCount + 1;
			}elsif($c =~ / /){
				$c = "";
			}elsif( $c =~/\n/){
				#print("Found a newline char");
				$c = "";
			}else{
				chomp;
				#print $c;
				$charCount = $charCount + 1;
			}
			print OUT $c;
		}
		
	}
	#print OUT "\n";
	
	
}

