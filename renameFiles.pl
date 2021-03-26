
#!/usr/bin/perl
# renameFiles.pl
# PURPOSE: To rename files that have been downloaded in fastq-dump
#          format for dada2 pipeline.
# INPUTS: 
#  -i: inDir - an input directory containing files with which to rename.
#  -p: 1 for paired, 0 for single reads.
# OUTPUTS:
#  1. Will ask if the file renaming scheme is appropriate (sent to STDOUT)
#  2. If no, will stop. If yes, will rename files

use warnings;
use strict;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);


my $inDir; # Directory containing files to rename (raw fastq files)
my $paired; # Do files contain paired data or not
my $name; # study name
my $forwardFileTail = "_1.fastq"; # forward file tail pattern of old file name
my $reverseFileTail = "_2.fastq"; # reverse file tail pattern of new file name
my $forwardReplace = ".R1.fastq"; # new file tail pattern forward
my $reverseReplace = ".R2.fastq"; # new file tail pattern reverse 

GetOptions(

    "i:s" => \$inDir,
    "p:i" => \$paired,
    "n:s" => \$name,
    "f:s" => \$forwardFileTail,
    "r:s" => \$reverseFileTail,
    "g:s" => \$forwardReplace,
    "s:s" => \$reverseReplace,

);

opendir(DH, $inDir) or die ("renameFiles.pl: Can't open input directory at ${inDir}. Please check the -i option\n");
my @files = readdir(DH);
closedir(DH);

my @cmds;

foreach my $file (@files){

    # For paired reads 
    if($paired == 1){

	my $newName = getNewFileHandlePaired($name, $file, $forwardFileTail, $reverseFileTail, $forwardReplace, $reverseReplace);
	my $cmd = "mv ${inDir}${file} ${inDir}${name}_${newName}";

	# Don't push any empty file strings that don't match our pattern
	if ($newName eq ""){}

	# Push an Array of commands so that we can make sure we are renaming the right thing
	else{

	    push @cmds, $cmd;
    
	}
    }

    # For non paired end reads
    elsif($paired == 0){
	
	print "Branch still under development for non-paired files\n";
    
    }
    
    # Paired value is something not allowed.
    else{

	print "renameFiles.pl: Sorry, that is not a valid option. Please check the -p option and input 1 for paired files and 0 for single files\n";

    }

}

# Print user message
my $message = "The following previous commands are the basic format for how files will be renamed. Choose a number for how to proceed.\n  1 - Rename with these commands.\n  2 - Cancel, these names aren't correct.\n  3 - Print out all commands that will be executed\n";



# Get user input
my $toPrint = 5;
getUserInput($message, $toPrint, @cmds);


sub getNewFileHandlePaired{

    my ($na, $fi, $of, $or, $nf, $nr) = @_;
    my $newFH = $fi;

    if($fi =~ m/${of}/){

	$newFH =~ s/${of}/${nf}/g;


    }elsif($fi =~ m/${or}/){
	
	$newFH =~ s/${or}/${nr}/g;

    }else{

	$newFH = "";

    }

    return($newFH);

}

sub getUserInput{

    my ($m, $p, @a) = @_;

    my $Ncmds = scalar @a;

    # Print some of the first few commands
    if($Ncmds < $p){print $m . printArray(@a);}
    else{print $m . printArray(@a[1.."${p}"]);}

    my $choice = <STDIN>;


    if(!looks_like_number($choice)){
	getUserInput($m, $p, @a);
    }
    elsif($choice == 1){
	
	print "Proceeding with commands\n";
	foreach my $moveCmd (@a){
	    print $moveCmd . "\n";
	    `$moveCmd`;
	}
    }

    elsif($choice == 2){
	
	die("renameFiles.pl: Cancelled. you have killed the run.\n");

    }

    elsif($choice == 3){
	print("Printing out all commands that will be performed\n");
	getUserInput($m, $Ncmds-1, @a);

    }else{
	
	print "Please choose a valid option for what to do next\n";
	getUserInput($m, $p, @a);
    }
}

sub printArray{
    my (@array) = @_;
    foreach my $e (@array){
	print $e, "\n";
    }
}


