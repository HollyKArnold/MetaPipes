#!/usr/bin/perl
# getAccessions.pl
# PURPOSE: Get a list of file accessions using fastq-dump.
#  1. Check data is paired or not. Download paired data or forward only

use warnings;
use strict;
use Getopt::Long;

my $inFile; # accession file
my $outDir; # file to dump data to

GetOptions(

    "i:s" => \$inFile,
    "o:s" => \$outDir,
);


open ACC, $inFile or die ("getAccessions.pl: cannot find infile at ${inFile}. Please check the -i option");

my $count = 0;
my $paired;

while(<ACC>){
    $count = $count + 1;
    my $line = $_;
    chomp $line;

    if($line eq ""){
	print "Reached end of accessions file. \n";
	return;
    }

    if($count == 1){
	my $cmd = "bash sraPaired.sh $line";
	$cmd = `$cmd`;
	
	if($cmd =~ m/TRUE/){
	    $paired = 1;
	    print "The data are paired. Downloading in paired end mode now...\n";

	}elsif($cmd =~ m/FALSE/){
	    $paired = 0;
	    print "The data are NOT paired. Downloading in single end mode now...\n";
	}else{
	    die( "\n ERROR: Cannot determine if reads are paired or not\n");
	}

     
    }


    if($paired == 1){

	my $cmd = "fastq-dump -O ${outDir} -I --split-files $line\n";
	`$cmd`;

    }elsif($paired == 0){

	my $cmd = "fastq-dump -O ${outDir} $line";
	`$cmd`;

    }
    
}
