# AUTHOR: ARNOLD
# DAY: Feb 5, 2021
# DATE: 20210205
# PURPOSE: The purpose of this file is to provide a set of helper functions for metaanalysis integration
# Eventually I will make this into an R package


# primerHits
# INPUTS
#   primer: a list from AllOrients of all possible orientations of your primer
#   fn: a path to the file that you wish to read
# OUTPUTS
#   Prints to screen the number of primers it finds in each file
# USE 
#   Uses max.mismatch of 0.
primerHits <- function(primer, fn) {
    
	#Count the number of reads that the primer sequence is found in
	nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = F, max.mismatch = 0)
	return(sum(nhits > 0))
}

# allOrients
# PURPOSE: Create a vector of all possible oreintations of primer sequence
# INPUTS:
#   primer: a string of letters representing a primer
# OUTPUTS
#   a vector containing all possible orientations of the provided primer sequence
allOrients <- function(primer) {
    # Create all orientations of the input sequence
    require(Biostrings)
    dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
    orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
        RevComp = reverseComplement(dna))
    return(sapply(orients, toString))  # Convert back to character vector
}

# getN
# INPUTS: a vector
# OUTPUTS: returns a count of unique items in the vector
getN <- function(x) sum(getUniques(x))

# my.cat
# INPUTS: Something to print to the screen
# OUTPUTS: Prints out print to screen in a nice way
my.cat <- function(x) {cat(paste0("\n### ", x, "\n"), sep = "\n")}

# my.bash 
# INPUTS: Something to print into a bash file
# OUTPUTS: prints out something for bash file. You still need to start sink before this.
my.bash = function(x){
	cat(paste0(x, "\n"))
}

# dada2.envPrep
# PURPOSE: Loads files and paths that should be the same from dada2 run to run. 
# INPUTS: NONE
# OUPUTS: 
#  run.env: returns the run environment
#  if it doesn't exist, makes a dada 2 output directory
dada2.envPrep = function(){
	
	# Load Libraries
	load.libs()

	# Project name
	proj.name <- tail(strsplit(getwd(), "/")[[1]], 1)

	# Dada2 Version Number	
	dada.pkg.ver <- paste0("dada2_", packageVersion("dada2"))
	
	# Date of run
	process.date <- Sys.Date()
	
	# Process ID
	process.id <- paste(dada.pkg.ver, process.date, sep = "_")
	
	# Output for dada2 folder
	output <- paste0(c(proj.name, "_", process.id, "_output"), sep = "", collapse = "")

	#Create a dada2 directory
	if (!dir.exists(output)) {
 		dir.create(output)
	}
	
	# Create a new dada2 environment to do these runs in
	run.env <- new.env()

	# Set run environment variables
	run.env$output = output
	run.env$process.id = process.id
	run.env$process.date = process.date
	run.env$proj.name = proj.name
	run.env$dada.pkg.ver = dada.pkg.ver

	# Get session info
	run.env$sessionInfo = sessionInfo()
	
	# Write session info
	my.cat("Writing out session info to sessionInfo.txt")
	writeLines(capture.output(sessionInfo()), paste0(output, "/sessionInfo.txt"))	

	return(run.env)

}

# load.libs
# PURPOSE: Load libraries necessary for the run
# INPUTS: None
# OUTPUTS: String with success of library load
load.libs = function(){
	 # List of required packages
	 packages = c("dada2", "R.utils", "ShortRead", "Biostrings", "stringr", "utils", "ggplot2", "reshape2", "phyloseq", "xtable")
	  
	 package.check(packages)

	  
}

# package.check()
# PURPOSE: checks for appropriate packages installed for the run.
# INPUTS: 
#   packages: A list of packages to load for the run
# OUTPUTS: Attempts to install packages if not present.
package.check =  function(packages){
	      
	      # Make a variable to decide if all packages were loaded correctly
	      success = 1
	      
	      packageListString = paste0(packages, collapse = " ")
	      lapply(packages, FUN = function(x) {if (!require(x, character.only = TRUE)) {
	      		       	     		success = 0
	      		       	     		my.cat("The following required packages did not load properly:")
						my.cat(x)
						my.cat("Attempting to install them now.")
                                                install.packages(x, dependencies = TRUE)
                                                library(x, character.only = TRUE)
                                                }else{
							#my.cat(paste0(c("Loading: ", x), collapse = " ", sep = " "))
						}
									                }
                )
		if(success == 1){
			   my.cat(paste0(c("Successfully loaded all packages in required package list:", packageListString), collapse = " "))
		}else{
			my.cat("Warning: all of the required packages were not loaded.")
		}

}

# dada2.userDefinedParams
# INPUTS:
#   file: A link to the user defined parameter file with the following attributes
#   	  col1: the user parameters that need to be defined
#	  col2: the user provided variables that must be filled out.
#	  Tab separated file
#   env: An environment object from dada2.envPrep() which contains the output directory
# OUTPUTS:
#	env: Returns 
dada2.userDefinedParams = function(file, env){

        # Read in options file
	my.cat(paste0("Setting user defined params at", file))
	opts = read.table(file, sep = "\t", header = TRUE)

	# Set silva file
	env$silvaTax = as.vector(opts[which(opts[,1] == "silvaTrainingSet"), 2])
	my.cat(paste0("Setting silva training set to: ", env$silvaTax))

	# Set silva species file
	env$silvaSpecies = as.vector(opts[which(opts[,1] == "silvaSpeciesSet"), 2])
	my.cat(paste0("Setting silva species training set to: ", env$silvaSpecies))

	# Set cut adapt file path
	env$cutAdapt = as.vector(opts[which(opts[,1] == "cutadapt"), 2])
	my.cat(paste0("Setting cut adapt file path to: ", env$cutAdapt))

        # set silva alignment file
        env$silvaAln = as.vector(opts[which(opts[,1] == "silvaAlign"), 2])
        my.cat(paste0("Setting silva alignment file to: ", env$silvaAln))

        # set silva one line alignment file
        env$silvaAlnOneLine = as.vector(opts[which(opts[,1] == "silvaAlignOneLine"), 2])
        my.cat(paste0("Setting silva alignment one line file to: ", env$silvaAlnOneLine))

        # Set the silva Info file
        env$silvaInfo = as.vector(opts[which(opts[,1] == "info"), 2])
        my.cat(paste0("Setting the silva info file to: ", env$silvaInfo))

        # Set the silva tree file
        env$silvaTree = as.vector(opts[which(opts[,1] == "silvaTree"), 2])
        my.cat(paste0("Setting the silva tree file to: ", env$silvaTree))
    
        # Ask what primer you want to use. 	
  	message = "Choose option for what primer to remove: \n1. EMP 16S sequencing protocol 515F(Parada)-806R(Apprill) \n2. EMP 16S original 515F(Caporaso)-806RCaporaso)\n3. Enter your own custom primer sequences.\n4. Don't remove primer sequences.\n\n"	

	primers = readintegerPrimer(max = 4, min = 1, message = message)
	env$FWD = as.vector(primers["FWD"])
	env$REV = as.vector(primers["REV"])
	env$DIST = as.vector(primers["DIST"])
	
	# Setting primer sequences
	#env$FWD = as.vector(opts[which(opts[,1] == "FWD"), 2])
	#env$REV = as.vector(opts[which(opts[,1] == "REV"), 2])
	#my.cat(paste0("Setting forward primer to: ", env$FWD))
	#my.cat(paste0("Setting reverse primer to: ", env$REV))
	
	# Path to python3
	env$python3 = as.vector(opts[which(opts[,1] == "python3"), 2])
	my.cat(paste0("Setting python3 path to: ", env$python3))

	# Path to figaro
	env$figaro = as.vector(opts[which(opts[,1] == "figaro"), 2])
	my.cat(paste0("Setting figaro path to: ", env$figaro))

	# Return environment
	return(env)
}

# FUNCTION: readintegerPrimer
# INPUTS:
#  max: the maximum integer that is allowed for input
#  min: the minimum integer that is allowed for input
# OUTPUT: Prints to screen until recieves valid option and then returns appropriate primer
readintegerPrimer <- function(max, min, message)
{



  # Get user input
  n <- readline(prompt=message)

  # Here are some commonly used primers
  empPrimerNewForward = "GTGYCAGCMGCCGCGGTAA"
  empPrimerNewReverse = "GGACTACNVGGGTWTCTAAT"
  empNewDist = 806 - 515
  empPrimerOldForward = "GTGCCAGCMGCCGCGGTAA"
  empPrimerOldReverse = "GGACTACHVGGGTWTCTAAT"
  empOldDist = 806 - 515

  if(!grepl("^[0-9]+$",n))
  {
    return(readintegerPrimer(max = max, min = min, message = message))
  
  }else if(as.integer(n) > max){

    return(readintegerPrimer(max = max, min = min, message = message))

  }else if(as.integer(as.integer(n)<min)){
    
    return(readintegerPrimer(max = max, min = min, message = message))
  
  }else if(as.integer(n == 1)){

    # User selects new EMP primers
    my.cat("You have selected option 1 for the new Earth Microbiome Project Primers!")
    my.cat(paste0(c("Setting forward primer to: ", empPrimerNewForward, ". Setting reverse primer to: ", empPrimerNewReverse, 
    		". Setting expected amplifed bp distance to 806 - 515 = ", empNewDist), sep = "", collapse = ""))
    return(list = c(CUT = TRUE, FWD = empPrimerNewForward, REV = empPrimerNewReverse, DIST = 806-515))

  }else if (as.integer(n == 2)){

    # User selects original EMP primers
    my.cat("You have selected option 2 for the original Earth Microbiome Project Primers!")
    my.cat(paste0(c("Setting forward primer to: ", empPrimerOldForward, ". Setting reverse primer to: ", empPrimerOldReverse,
    		". Setting distance expected amplified distance to 806 - 515 = ", empOldDist), sep = "", collapse = ""))
    return(list = c(CUT = TRUE, FWD = empPrimerOldForward, REV = empPrimerOldReverse, DIST = 806 - 515))
    
  }else if(as.integer(n == 3)){

    # User has decided to enter their own custom primers
    my.cat("You have selected option 3 to input your own primer sequences.")
    FWD = dada2.getPrimer(message = "Please enter sequence in for forward primer:")
    REV = dada2.getPrimer(message = "Please enter sequence in for reverse primer:")
    primers = dada2.getPrimerStartStop(messageStart = "Please enter integer of start primer site (i.e. 515):\n", 
    	    			     messageStop =  "Please enter integer of reverse primer start site (i.e. 806):\n")
    DIST = as.vector(primers["stop"]) - as.vector(primers["start"])
    return(list = c(CUT = TRUE, FWD = FWD, REV = REV, DIST = DIST))
    
  }else if(as.integer(n == 4)){
    my.cat("You have selected option 4 to not remove any primers. Development of this branch of code is still under process.")
    return(list = c(CUT = FALSE))
  }

  return(as.integer(n))
}

# FUNCTION: dada2.getPrimer
# INPUTS:
#   message: message to print to screen
# OUTPUTS:
#   returns user entered primer string. 
dada2.getPrimer = function(message){

     # Get user input		       
     n <- readline(prompt=message)
     x = substring(n, 1:nchar(n), 1:nchar(n))
     
     possibleChars = c("A", "C", "T", "G", "U", "N", "Y", "R", "W", "S", "M", "K", "B", "D", "H", "V")

     validString = TRUE
     for(i in 1:length(x)){
     
     # Make sure there is valid input
     if(! x[i] %in% possibleChars)
     	  validString = FALSE
     }
     
     # If not valid input, ask for user to re-enter code
     if(! validString){
     	  
	  my.cat("Please recheck your primer sequence. Valid characters are: ")
	  print(possibleChars)
	  dada2.getForwardPrimer(message = message)
    
     }

     # Print warning if primer string is shorter than 20 characters
     if(length(x) < 20){

          warning("Primer sequences are usually at least 20 characters long")

     }


     return(n)
}

# FUNCTION
# INPUTS
#   messageStart: Message to print to user to get first primer
#   messageStop: Message to print to user to get second primer
# OUTPUTS
#   primer: returns a list with "start" and "stop" for the primer
dada2.getPrimerStartStop = function(messageStart, messageStop){
  
     # Get primer start
     n <- readline(prompt=messageStart)

     # If not an integer, redirect
     if(!grepl("^[0-9]+$",n)){
	
	dada2.getPrimerStartStop(messageStart = messageStart, messageStop = messageStop)
     
     }
     
     # If an invalid integer, redirect
     if(as.integer(n) < 1 | as.integer(n) > 1550){
     	  my.cat("Valid integers for start sites on the 16S rRNA gene are between 1 and 1550 basepairs.")
     	  dada2.getPrimerStartStop(messageStart = messageStart, messageStop = messageStop)
     }
     start = as.integer(n)

     n <- readline(prompt = messageStop)
     
     # If not an integer, redirect
     if(!grepl("^[0-9]+$",n)){
        dada2.getPrimerStartStop(messageStart = messageStart, messageStop = messageStop)
     }
     
     # If an invalid integer, redirect
     if(as.integer(n) < 1 |  as.integer(n) > 1550){
          my.cat("Valid integers for stop sites on the 16S rRNA gene are between 1 and 1550 basepairs.")
     	  dada2.getPrimerStartStop(messageStart = messageStart, messageStop = messageStop)
     }
     stop = as.integer(n)

     # Check that stop is greater than start
     if(start > stop){
     	 my.cat("Stop site should be greater than the start site. ")
	 dada2.getPrimerStartStop(messageStart = messageStart, messageStop = messageStop)
     }

     return(list = c(start = start, stop = stop))

}
dada2.getPrimerEnd = function(message){
  n <- readline(prompt=message)
    my.cat(message)
    return(100)
}

# dada2.checkMetadataMatchFastQ
# PURPOSE: Read in metadata file, see if it matches the fastq files
# INPUTS:
#   fastq.path: The path to fastq files.
#   metadata: The metadata file.
# OUTPUTS:
dada2.checkMetadataMatchFastQ = function(fastq.path, metadataFile, env){

	print(head(list.files(fastq.path), 20))
  	proceed <- readline(prompt = "Above are the first 20 (or fewer) files in the provided path, do you want to proceed? [y/n]: ")
  	while (!(proceed %in% c("y", "n"))) {
    		proceed <- readline(prompt="Please answer y or n: ")
  	}
  	if (proceed == "n") {
    		my.cat("Terminated")
  	} else {
	  names = basename(list.files(fastq.path))
	}

	#metadata = read.table(metadataFile, sep = "\t", header = T)
        metadata = read.csv(metadataFile, sep = "\t")
        print(head(metadata))
	proceed = readline(prompt = "Above is the head of the metadatafile in the provided file, do you want to proceed? [y/n]: ")
	
        while (!(proceed %in% c("y", "n"))) {
                proceed <- readline(prompt="Please answer y or n: ")
	}
	if(proceed == "n"){
		   my.cat("Terminated")
	}else{
		metadatanames = rownames(metadata)

	}


	my.cat("Checking if metadatanames and fastqfiles are identical...")
	fastqNameBase = unique(gsub(pattern=".R[12].fastq.gz", replacement="", x = names))
	if("filtN" %in% fastqNameBase){
	    	idx = which(fastqNameBase == "filtN")
		fastqNameBase = fastqNameBase[-idx]
	}
	if("cut" %in% fastqNameBase){
                idx = which(fastqNameBase == "cut")
                fastqNameBase = fastqNameBase[-idx]
        }

	if(sum(fastqNameBase %in% metadatanames) != length(fastqNameBase)){
		my.cat("Some fastQ files were not found in the metadata file.")
		my.cat("Please delete these fastQ files and rerun.")
		print(fastqNameBase[!fastqNameBase %in% metadatanames])
		toRemove = fastqNameBase[!fastqNameBase %in% metadatanames]
		sink(paste0(env$output, "/fastQFilesToRemove.sh"))
		my.bash("#!/bin/bash")
		for (i in 1:length(toRemove)){
		    my.bash(paste0(c("rm ", toRemove[i], ".R1.fastq.gz"), sep = "", collapse = ""))
		    my.bash(paste0(c("rm ", toRemove[i], ".R2.fastq.gz"), sep = "", collapse = ""))
		}
		sink()
		my.cat("You will have to rerun your R script, but I printed a bash file for you at:")
		my.cat(paste0(env$output, "/fastQFilesToRemove.sh"))
		my.cat("Use with caution!")
		stop("Please make sure metadata rownames and fastq filenames match before proceeding...")
	}else{
		my.cat("Perfect! The metadata file rownames and the fastq files match! Proceeding...")
	}
	

	# Assign values
	env$fastqNames = names
	env$metadataNames = metadatanames
	env$fastq.path = fastq.path
	env$metadataFile = metadataFile
	env$metadata = metadata

	# Return environemnt file
	return(env)
}

# dada2.preFilter
# INPUTS: 
#   env = an enviornment obtained from dada2.checkMetadataMatchFastQ
# OUTPUT:
#   prefilters sequences before cutting into fastq directory under filtN/ with the number of cores for processors
dada2.preFilter = function(env, cores){
	
	# Get forward and reverse files
	fnFs = sort(list.files(env$fastq.path, pattern = ".R1.fastq.gz", full.names = T))
	fnRs = sort(list.files(env$fastq.path, pattern = ".R2.fastq.gz", full.names = T))

	# Get all orientations of the primer
	FWD.orients = allOrients(env$FWD)
	REV.orients = allOrients(env$REV)

	# Create file names for all the filtered seqs	
	my.cat("Pre-filtering sequences before cutting primers:")
	fnFs.filtN = file.path(env$fastq.path, "filtN", basename(fnFs))
	fnRs.filtN = file.path(env$fastq.path, "filtN", basename(fnRs))

	# Now filter the sequences
	out.fn.filtN = filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = cores, compress = T)

	# Add our steps for this function to the environment varaible
	env$fnFs = fnFs
	env$fnRs = fnRs
	env$FWD.orients = FWD.orients
	env$REV.orients = REV.orients
	env$fnFs.filtN = fnFs.filtN
	env$fnRs.filtN = fnRs.filtN
	env$out.fn.filtN = out.fn.filtN

	return(env)
}

# dada2.checkForPrimersInitial
# INPUTS
#   env: an environment obtained from dada2.preFilter
# OUTPUTS
#   Prints out a list of primers found pre cut
dada2.checkForPrimersInitial = function(env){
	
	# Number of files	
	N = length(env$fnFs.filtN)

	# Make a message for user to see how many primers to look for before cutting
	message = paste0(c("Please enter number of files to initially scan for primer presence. If you would like to print all files then enter in ", N, "\n"), sep = "", collapse = "")	
	proceed = readinteger(max = N, message = message)
	if(proceed == 0){
	}else{
	
		# Print out initial primer finds
		my.cat(paste0(c("Printing out initial primers found for first ", proceed, " files:"), collapse = "", sep = ""))
		sink(paste0(env$output, "/initialPrimers.txt"))
		for(i in 1:proceed){ 
	      	      print(i)
	      	      if (file.exists(env$fnFs.filtN[[i]])){print(rbind(
	      	      	 FWD.ForwardReads = sapply(env$FWD.orients, primerHits, fn = env$fnFs.filtN[[i]]), 
		 	 FWD.ReverseReads = sapply(env$FWD.orients, primerHits, fn = env$fnRs.filtN[[i]]), 
		 	 REV.ForwardReads = sapply(env$REV.orients, primerHits, fn = env$fnFs.filtN[[i]]), 
		 	 REV.ReverseReads = sapply(env$REV.orients, primerHits, fn = env$fnRs.filtN[[i]])))
	      	      } 
		}
		sink()
	}

}

# dada2.cut
# INPUTS
#  env: an environment obtained from dada2.preFilter
#  cores: the number of cores for cut adapt
#  errors: the number of errors allowed for cut adapt. Default is 0.1
# OUTPUTS
#  a set of cut sequences and a file with theoutput from cutadapt
#  modifies environment by adding paths for cut variables
dada2.cut = function(env, cores, errors = 0.1){
	# Make a path for cut files
	env$path.cut = file.path(env$fastq.path, "cut")

	# Make the directory for cut files
	if(!dir.exists(env$path.cut)) dir.create(env$path.cut)
	env$fnFs.cut = file.path(env$path.cut, basename(env$fnFs))
	env$fnRs.cut = file.path(env$path.cut, basename(env$fnRs))

	# Make primer reverse compliment
	env$FWD.RC = dada2:::rc(env$FWD)
	env$REV.RC = dada2:::rc(env$REV)
	
	# Trim FWD and reverse compliment off of foward reads
	env$R1.flags = paste("-g", env$FWD, "-a", env$REV.RC)

	# Trim reverse and foward reverse compliment off of reverse reads
	env$R2.flags = paste("-G", env$REV, "-A", env$FWD.RC)

	# Now, run cut adapt!
	# -n is required to remove FWD and REV from reads
	# -m required to not cut sequence to length 0 and kill the quality plot later
	# fnFs.cut and fnFs.rev are output files
	# fnFs.filtN and fnRs.filtN are the input files.
	my.cat("Running cut adapt:")
	sink(paste0(env$output, "/cutAdaptOutput.txt"))
	
	
	for(i in seq_along(env$fnFs)){
	      if(file.exists(env$fnFs[i])){
		x = system2(env$cutAdapt,  args = c(env$R1.flags, env$R2.flags, 
				  "-n", 2, "-m", 1, "-j", cores, "-e", errors, "-o", env$fnFs.cut[i], "-p", env$fnRs.cut[i], 
				  env$fnFs.filtN[i], env$fnRs.filtN[i]), stdout = TRUE, wait = TRUE)
	      	my.cat(x)
	      }
	}
	sink()

	# Name cut forward and reverse files
	env$cutFs = sort(list.files(env$path.cut, pattern = ".R1.fastq.gz", full.names = T))
	env$cutRs = sort(list.files(env$path.cut, pattern = ".R2.fastq.gz", full.names = T))
	get.sample.names = function(fname) strsplit(basename(fname), ".R1.fastq.gz")[[1]][1]
	env$sample.names = unname(sapply(env$cutFs, get.sample.names))


	# Return environment
	return(env)
}

# readinteger
# INPUTS:
#  max: the maximum integer that is allowed for input
# OUTPUT: Function to make sure an integer is read that is below max
readinteger <- function(max, message)
{ 

  n <- readline(prompt=message)


  if(!grepl("^[0-9]+$",n))
  {
    return(readinteger(max = max, message = message))
  }else if(as.integer(n) > max){
    return(readinteger(max = max, message = message))
  }
  
  return(as.integer(n))
}

# dada2.checkForPrimersPostCut
# PURPOSE: TO check for primers in the files post cut
# INPUTS: an environment returned from dada2.cut
# OUTPUTS: 
#   makes no changes to environment variable
#   outputs a file with primers found in the cut sequences. 
dada2.checkForPrimersPostCut = function(env){
        
	# Get the total number of files to check for primer post cut
	N = length(env$fnFs.filtN)

	# Print message
	message = paste0(c("Please enter number of files to scan for primer presence post cut. If you would like to print all files then enter in ", N, "\n"), sep = "", collapse = "")

	# Get the number of files to print
	proceed = readinteger(max = N, message = message)

	# User doesn't want to check any post files
        if(proceed == 0){
        }else{

		# Print out primer finds post cut
                my.cat(paste0(c("Printing out primers found in cut sequences for first ", proceed, " files:"), collapse = "", sep = ""))
                sink(paste0(env$output, "/cutPrimers.txt"))
		for(i in 1:proceed){
                      print(i)
                      if (file.exists(env$fnFs.cut[i])){print(rbind(
                         FWD.ForwardReads = sapply(env$FWD.orients, primerHits, fn = env$fnFs.cut[i]),
                         FWD.ReverseReads = sapply(env$FWD.orients, primerHits, fn = env$fnRs.cut[i]),
                         REV.ForwardReads = sapply(env$REV.orients, primerHits, fn =env$ fnFs.cut[i]),
                         REV.ReverseReads = sapply(env$REV.orients, primerHits, fn = env$fnRs.cut[i])))
                      }
                }
		sink()
	}

}

# dada2.printCutQualityPlots
# INPUTS:
#   env: An enviornment function returned from dada2.cut
#   thresh: A number to not print a quality plot (i.e. less than 1000 seqs will get thrown out)
#   space: The space between quality plot figures
#   figsize: size of quality plot figures for latex
#   nsubfigure: The number of subfigures per latex plot 
# OUTPUTS:
#   A set of quality plots in "qualityPlotsCut"
dada2.printCutQualityPlots = function(env, thresh = 1000, space = 0.5, figsize = 50, nsubfigure = 8){

     # Number of files
     N = length(env$sample.names)

     # Make a message for user to see how many primers to look for before cutting
     message = paste0(c("Please enter the number of quality plots to print out. If you want to print them all, then enter: ", N, "\n"), sep = "", collapse = "")

     proceed = readinteger(max = N, message = message)
     if(proceed == 0){
        my.cat("You have chosen to not print any quality plots.")
	return(1)
     }else{
	
	
	# Make a path for quality plots
        env$path.qualCut = file.path(env$output, "qualityPlotsCut")

        # Make the directory for cut files
        if(!dir.exists(env$path.qualCut)) dir.create(env$path.qualCut)

	# Make a progress bar for forward plots
	my.cat("Printing forward quality plot scores...")
	pb = txtProgressBar(min = 0, max = proceed)
	
	# Make the forward quality plots
	for(i in 1:proceed){
	      setTxtProgressBar(pb, i)
      	      if(file.exists(env$cutFs[i]) & env$out.fn.filtN[i,"reads.out"] > thresh){

	           # Name the plot
	           plotName = paste0(c("/", env$sample.names[i], "ForwardCutQualityScorePlot.pdf"), sep = "", collapse = "")
 		   plotName = gsub(x = plotName, pattern = "_", replace = "")
		   
		   # Print the plot
		   pdf(file = paste0(env$path.qualCut, plotName))
  	     	   print(plotQualityProfile(env$cutFs[i]))
	      	   dev.off()
	      }
  	}
	close(pb)

	# Make a progress bar for reverse plots
	my.cat("Printing reverse quality plot scores...")
	pb = txtProgressBar(min = 0, max = proceed)
	
	# Make the reverse quality plots
	for(i in 1:proceed){
	      setTxtProgressBar(pb, i)
              if(file.exists(env$cutRs[i]) & env$out.fn.filtN[i,"reads.out"] > thresh){

	           # Name the plot
                   plotName = paste0(c("/", env$sample.names[i], "ReverseCutQualityScorePlot.pdf"), sep = "", collapse = "")
		   plotName = gsub(x = plotName, pattern = "_", replace = "")                   

		   # Print the plot
		   pdf(file = paste0(env$path.qualCut, plotName))
                   print(plotQualityProfile(env$cutRs[i]))
                   dev.off()
              }
        }
	close(pb)

	# Now make latex output file so that we can easily import these into latex
	latexwritefile = paste0(c(env$path.qualCut, "/latexQualityCutPlotsSubfigures.txt"), sep = "", collapse = "") 
	latexwritefileFigs = paste0(c(env$path.qualCut, "latexFigs.txt"), sep = "", collapse = "")

	qualityScoreNames = gsub(x=list.files(env$path.qualCut, pattern = ".pdf", full.names =F), pattern = ".pdf", replacement = "")
	writeLatexQualityPlotFiles(writefile = latexwritefile, append = F, space= space, figsize= figsize, filelist=qualityScoreNames, nsubfigure=nsubfigure)
   }
   return(env)

}



# dada2.printFilteredQualityPlots
# INPUTS:
#   env: An environment function returned from dada2.cut
#   thresh: A number below which to not print a quality for
# OUTPUTS:
#   A set of quality plots in "qualityPlotsFiltN"
dada2.printFiltNQualityPlots = function(env, thresh = 1000, space = 0.5, figsize = 50, nsubfigure = 8){

			     
     # Number of files
     N = length(env$sample.names)

     # Make a message for user to see how many primers to look for before cutting
     message = paste0(c("Please enter the number of quality plots to print out. If you want to print them all, then enter: ", N, "\n"), sep = "", collapse = "")

     proceed = readinteger(max = N, message = message)
     if(proceed == 0){
        my.cat("You have chosen to not print any quality plots.")
        return(1)
     }else{


       # Make a path for quality plots
       env$path.qualFiltN = file.path(env$output, "qualityPlotsFiltN")

       # Make the directory for quality plots of filtered files
       if(!dir.exists(env$path.qualFiltN)) dir.create(env$path.qualFiltN)

       # Make a progress bar for forward plots
       my.cat("Printing forward quality plot scores...")
       pb = txtProgressBar(min = 0, max = proceed)

       # Make the forward quality plots
       for(i in 1:proceed){
              setTxtProgressBar(pb, i)
              if(file.exists(env$fnFs.filtN[i]) & env$out.fn.filtN[i,"reads.out"] > thresh){

                   # Name the plot
                   plotName = paste0(c("/", env$sample.names[i], "ForwardFilteredQualityScorePlot.pdf"), sep = "", collapse = "")
		   plotName = gsub(x = plotName, pattern = "_", replace = "")

                   # Print the plot
                   pdf(file = paste0(env$path.qualFiltN, plotName))
                   print(plotQualityProfile(env$fnFs.filtN[i]))
                   dev.off()
              }
        }
        close(pb)

        # Make a progress bar for reverse plots
        my.cat("Printing reverse filtered quality plot scores...")
        pb = txtProgressBar(min = 0, max = proceed)

        # Make the reverse quality plots
        for(i in 1:proceed){
              setTxtProgressBar(pb, i)
              if(file.exists(env$fnRs.filtN[i]) & env$out.fn.filtN[i,"reads.out"] > thresh){

                   # Name the plot
                   plotName = paste0(c("/", env$sample.names[i], "ReverseFilteredQualityScorePlot.pdf"), sep = "", collapse = "")
		   plotName = gsub(x = plotName, pattern = "_", replace = "")		   

                   # Print the plot
                   pdf(file = paste0(env$path.qualFiltN, plotName))
                   print(plotQualityProfile(env$fnRs.filtN[i]))
                   dev.off()
              }
        }
        close(pb)

        # Now make latex output file so that we can easily import these into latex
        latexwritefile = paste0(c(env$path.qualFiltN, "/latexQualityPlotsFilteredSubfigures.txt"), sep = "", collapse = "")
        latexwritefileFigs = paste0(c(env$path.qualFiltN, "latexFigs.txt"), sep = "", collapse = "")

        qualityScoreNames = gsub(x=list.files(env$path.qualFiltN, pattern = ".pdf", full.names =F), pattern = ".pdf", replacement = "")
        writeLatexQualityPlotFiles(writefile = latexwritefile, append = F, space = space, figsize= figsize, filelist=qualityScoreNames, nsubfigure=nsubfigure)
   }
   return(env)
}


# FUNCTION: writeLatexQualityPlotFiles
# PURPOSE:  Given a list of subfigure names outputs a latex documentation
#           for each file to be added as a subfigure.
# INPUTS:
#           fileList - list of file names 
#           space - space between figures
#           figSize - figure Size
#           writefile - the file to write to
#           nSubfigure - the number of plots to output to a single subfigure
# OUTPUTS:
#           prints out the latex subfigures.

writeLatexQualityPlotFiles = function(writefile, space, figsize, filelist, nsubfigure, append){
  
  sink(writefile, append = append)
  
  N = length(filelist)
  cur = 1
  curFig = 0
  
  while(cur <= N ){ # Make a new figure
    curFig = curFig + 1
    cat("\\begin{figure}[H]\n")
    
    while(cur <= nsubfigure*curFig & cur <= N){
      
      cat(paste0(c("\\begin{subfigure}[b]{", space, "\\textwidth}\n"), collapse = "", sep = ""))
      cat(paste0(c("\\includegraphics[width=", figsize, "mm]{", filelist[cur], "}\n"), collapse = "", sep = ""))
      cat(paste0(c("\\caption{", filelist[cur], "}\n"), collapse = "", sep = ""))
      cat(paste0(c("\\label{fig:", filelist[cur], "}\n"), collapse = "", sep = ""))
      cat("\\end{subfigure}\n")
      
      cur = cur + 1 #update current subfigure
      
    }
    
    cat("\\end{figure}\n\n")
  }
  
  sink()
}

# FUNCTION: dada2.filterCutSeqs
# INPUTS:
#    cores: The number of cores to run on
#    fwdTrimLength: The length to trim the forward reads
#    revTrimLength: The length to trim the reverse reads to
#    env: An environment output from dada2.printCutQualityPlots
#    thresh: The number of reads required to further consider sampels (i.e. throw out reads under this thresh)
# OUTPUTS
#   filters cut seqneuces
dada2.filterCutSeqs = function(env, fwdTrimLength, revTrimLength, cores, thresh){

    # Print message
    my.cat("Filtering cut sequences...")

    # Create file paths to filtered cut sequences		   
    env$filtFs = file.path(env$path.cut, "filtered", basename(env$cutFs))
    env$filtRs = file.path(env$path.cut, "filtered", basename(env$cutRs))
    
    # Filter sequences
    env$out.cut.filt = filterAndTrim(env$cutFs, env$filtFs, env$cutRs, env$filtRs, truncLen = c(fwdTrimLength, revTrimLength), maxN = 0, maxEE = c(1, 1), truncQ = 2, rm.phix = TRUE, compress = TRUE, multithread = cores)

    # Return env
    return(env)
}

# FUNCTION: dada2.filterCutSeqsForwardOnly
# INPUTS:
#    cores: The number of cores to run on
#    fwdTrimLength: The length to trim the forward reads
#    env: An environment output from dada2.printCutQualityPlots
#    thresh: The number of reads required to further consider sampels (i.e. throw out reads under this thresh)
# OUTPUTS
#   filters cut seqneuces for forward reads only
dada2.filterCutSeqsForwardOnly = function(env, fwdTrimLength, cores, thresh){

    # Print message
    my.cat("Filtering cut sequences...")

    # Create file paths to filtered cut sequences		   
    env$filtFs = file.path(env$path.cut, "filtered", basename(env$cutFs))
    
    # Filter sequences
    env$out.cut.filt = filterAndTrim(env$cutFs, env$filtFs, truncLen = c(fwdTrimLength), maxN = 0, maxEE = c(1), truncQ = 2, rm.phix = TRUE, compress = TRUE, multithread = cores)

    # Return env
    return(env)
}

# FUNCTION dada2.learnErrors
# INPUTS: 
#    env: an environment returned from dada2.filterCutSeqs
#    cores: the number of cores to use for error learning
# OUTPUTS
#    learns errors for dada2 and plots graphs. 
#    errF: Forward errors added to env
#    errR: Reverse errors
dada2.learnErrors = function(env, cores, thresh){

    # Print message
    my.cat("Learning errors...")

    # Get the files that exist
    env$exists = file.exists(env$filtFs) & as.vector(env$out.cut.filt[,"reads.out"] > thresh) & file.exists(env$filtRs)
    
    # Learn errors
    env$errF = learnErrors(env$filtFs[env$exists], multithread = cores)
    env$errR = learnErrors(env$filtRs[env$exists], multithread = cores)

    #Plot Forward Errors
    errF.plot = plotErrors(env$errF, nominalQ = TRUE)
    ggsave(errF.plot, file = file.path(env$output, paste0(env$proj.name, "errFplot.pdf")))

    #Plot Reverse Errors
    errR.plot <- plotErrors(env$errR, nominalQ = TRUE)
    ggsave(errR.plot, file = file.path(env$output, paste0(env$proj.name, "errRplot.pdf")))

    
    return(env)
}  

# FUNCTION dada2.learnErrorsForwardOnly
# INPUTS: 
#    env: an environment returned from dada2.filterCutSeqs
#    cores: the number of cores to use for error learning
# OUTPUTS
#    learns errors for dada2 and plots graphs. 
#    errF: Forward errors added to env
dada2.learnErrorsForwardOnly = function(env, cores, thresh){

    # Print message
    my.cat("Learning errors...")

    # Get the files that exist
    env$exists = file.exists(env$filtFs) & as.vector(env$out.cut.filt[,"reads.out"] > thresh) 
    
    # Learn errors
    env$errF = learnErrors(env$filtFs[env$exists], multithread = cores)

    #Plot Forward Errors
    errF.plot = plotErrors(env$errF, nominalQ = TRUE)
    ggsave(errF.plot, file = file.path(env$output, paste0(env$proj.name, "errFplot.pdf")))
    
    return(env)
}  

# FUNCTION dada2.dada2
# PURPOSE: Run dada2
# INPUTS:
#    env: an env from dada2.learnErrors
#    pool: should sequences be pooled
#    cores: number of cores to use to run
# OUTPUTS:
#    dadaF: forward dada2 objects
#    dadaR: reverse dada2 object
dada2.dada2 = function(env, pool, cores){
     
     # Print
     my.cat("Running dada2...")
     
     # Run dada2 
     env$dadaFs = dada(env$filtFs[env$exists], err = env$errF, multithread = cores, pool = pool)
     env$dadaRs = dada(env$filtRs[env$exists], err = env$errR, multithread = cores, pool = pool)
     
     # Return env
     return(env)

}

# FUNCTION dada2.dada2ForwardOnly
# PURPOSE: Run dada2 for forward reads only
# INPUTS:
#    env: an env from dada2.learnErrors
#    pool: should sequences be pooled
#    cores: number of cores to use to run
# OUTPUTS:
#    dadaF: forward dada2 objects
dada2.dada2ForwardOnly = function(env, pool, cores){
     
     # Print
     my.cat("Running dada2...")
     
     # Run dada2 
     env$dadaFs = dada(env$filtFs[env$exists], err = env$errF, multithread = cores, pool = pool)
     
     # Return env
     return(env)

}

# FUNCTION dada2.mergers
# PURPOSE: To merge dadaF and dadaR objects
# INPUTS: env from dada2.dada2
# OUTPUTS:
#   mergers: the merged pairs
dada2.mergers = function(env){
    
    # Starting dada2.mergers
    my.cat("Merging dada2 forward and reverse... ")

    # Merge Pairs
    env$mergers = mergePairs(env$dadaFs, env$filtFs[env$exists], env$dadaRs, env$filtRs[env$exists], verbose = TRUE)
    
    # Print mergers
    my.cat(head(env$mergers[[1]]))

    # Return
    return(env)        	      
}

# FUNCTION: dada2.sequenceTableForwardOnly
# PURPOSE: Make sequence table from merged dada2 objects
# INPUTS:
#    env from dada2.mergers
# OUTPUTS:
#    seqtab: a sequence table from the mergers object
dada2.sequenceTableForwardOnly = function(env){
     
     # Make a sequence table
     env$seqtab = makeSequenceTable(env$dadaFs)
     
     # Get the dimensions of the sequence table
     my.cat("The original sequence table table dimensions:")
     my.cat(dim(env$seqtab))
     
     # Get the range of sequence lengths
     my.cat(table(nchar(getSequences(env$seqtab))))

     # Return the sequence table
     return(env)
      
}


# FUNCTION: dada2.sequenceTable
# PURPOSE: Make sequence table from merged dada2 objects
# INPUTS:
#    env from dada2.mergers
# OUTPUTS:
#    seqtab: a sequence table from the mergers object
dada2.sequenceTable = function(env){
     
     # Make a sequence table
     env$seqtab = makeSequenceTable(env$mergers)
     
     # Get the dimensions of the sequence table
     my.cat("The original sequence table table dimensions:")
     my.cat(dim(env$seqtab))
     
     # Get the range of sequence lengths
     my.cat(table(nchar(getSequences(env$seqtab))))

     # Return the sequence table
     return(env)
      
}

# FUNCTION: dada2.removeChimeras
# INPUTS
#     env: an env from dada2.sequenceTable
# OUTPUTS
#     seqtab.nochim: a sequence table witout chimeras
dada2.removeChimeras = function(env, cores){

     # Print out function heading
     my.cat("Removing chimeras from the sequence table...")
     
     # Remove chimeras
     env$seqtab.nochim = removeBimeraDenovo(env$seqtab, method = "consensus", multithread = cores, verbose = TRUE)
     
     # Print out the dimensions of the ASV table with removed chimeras
     my.cat("The dimensions of the sequence table with chimeras removed")
     my.cat(dim(env$seqtab.nochim))

     # Print out the percent of reads retained when chimeras are removed
     my.cat("The percentage of sequences retained after chimeras are removed")
     my.cat(sum(env$seqtab.nochim)/sum(env$seqtab))

     #Return environment
     return(env)
}

# FUNCTION: dada2.makeASVTable
# INPUTS:
#    env: an env from dada2.removeChimeras
#    fileTailPattern: A string that may be at the end of the rownames of the ASV file matching forward file names
# OUTPUTS
#    asv: makes an ASV table from the removed chimera table
dada2.makeASVTable = function(env, fileTailPattern = ".R1.fastq.gz"){
     
     # Run dada2.makeASVTable
     my.cat("Making ASV table...")

     # Make ASV Table
     env$asv = env$seqtab.nochim
     
     #Rename rownames
     rownames(env$asv) = sub(x = rownames(env$asv), pattern = fileTailPattern, replace = "")

     
     #Return env
     return(env)
}

# FUNCTION: dada2.makeCountTable
# INPUTS:
#     env: from dada2.makeASVTable
# OUTPUTS:
#     counts: a table of counts for each step in the filtering process
dada2.makeCountTable = function(env, fileTailPatternForward = ".R1.fastq.gz", fileTailPatternReverse = ".R2.fastq.gz"){
     
     
     # Count forward original reads
     my.cat("Counting the reads in original files...")

     # Count the number of sequences in the forward and reverse reads
     env$count.fnFs = sapply(env$fnFs, getN)
     env$count.fnRs = sapply(env$fnRs, getN)     

     # Rename
     names(env$count.fnFs) = sub(x = basename(path = names(env$count.fnFs)), pattern = fileTailPatternForward, replace = "")
     names(env$count.fnRs) = sub(x = basename(path = names(env$count.fnRs)), pattern = fileTailPatternReverse, replace = "")

     # Count the sequences filtered for N's
     my.cat("Counting the reads filtered for Ns...")

     # Count the number of filtered N sequences in the forward and reverse reads
     env$count.fnFs.filtN = sapply(env$fnFs.filtN[file.exists(env$fnFs.filtN)], getN)
     env$count.fnRs.filtN = sapply(env$fnRs.filtN[file.exists(env$fnRs.filtN)], getN)

     # Rename the filtered N sequences
     names(env$count.fnFs.filtN) = sub(x = basename(path = names(env$count.fnFs.filtN)), pattern = fileTailPatternForward, replace = "")
     names(env$count.fnRs.filtN) = sub(x = basename(path = names(env$count.fnRs.filtN)), pattern = fileTailPatternReverse, replace = "")

     # Message to count filtered reads
     my.cat("Counting the numbers of filtered cut reads...")
     
     # Count the filtered cut sequences
     env$count.filtFs = sapply(env$filtFs[file.exists(env$filtFs)], getN)
     env$count.filtRs = sapply(env$filtRs[file.exists(env$filtRs)], getN)

     # Rename the filtered cut sequences
     names(env$count.filtFs) = sub(x = basename(path = names(env$count.filtFs)), pattern = fileTailPatternForward, replace = "")
     names(env$count.filtRs) = sub(x = basename(path = names(env$count.filtRs)), pattern = fileTailPatternReverse, replace = "")

     # Count dadaFs and dadaRs
     my.cat("Counting the number of reads after dada2 algorithm...")
     
     # Count dada reads
     env$count.dadaFs = sapply(env$dadaFs, getN)
     env$count.dadaRs = sapply(env$dadaRs, getN)     

     # Rename dada reads
     names(env$count.dadaFs) = sub(x = basename(path = names(env$count.dadaFs)), pattern = fileTailPatternForward, replace = "") 
     names(env$count.dadaRs) = sub(x = basename(path = names(env$count.dadaRs)), pattern = fileTailPatternReverse, replace = "")

     # Count merged
     my.cat("Counting the number of mergers...")
     env$count.merged = sapply(env$mergers, getN)
     names(env$count.merged) = sub(x = basename(path = names(env$count.merged)), pattern = fileTailPatternForward, replace = "")

     # Count Nonchim
     my.cat("Counting the number of sequences after chimeras are removed...")
     env$count.nonchim = rowSums(env$seqtab.nochim)
     names(env$count.nonchim) = sub(x = basename(path = names(env$count.nonchim)), pattern = fileTailPatternForward, replace = "")

     # Make track
     env$track = data.frame(matrix(nrow = length(env$count.fnFs), ncol = 10, data = NA))
     colnames(env$track) = c("rawF", "rawR", "filteredNF", "filteredNR", "cutThreshF", "cutThreshR", "denoisedF", "denoisedR", "merged", "nonchim")

     rownames(env$track) = names(env$count.fnFs)
     env$track[names(env$count.fnFs), "rawF"] = env$count.fnFs 
     env$track[names(env$count.fnRs), "rawR"] = env$count.fnRs 
     env$track[names(env$count.fnFs.filtN), "filteredNF"] = env$count.fnFs.filtN 
     env$track[names(env$count.fnRs.filtN), "filteredNR"] = env$count.fnRs.filtN 
     env$track[names(env$count.filtFs), "cutThreshF"] = env$count.filtFs 
     env$track[names(env$count.filtRs), "cutThreshR"] = env$count.filtRs 
     env$track[names(env$count.dadaFs), "denoisedF"] = env$count.dadaFs 
     env$track[names(env$count.dadaRs), "denoisedR"] = env$count.dadaRs 
     env$track[names(env$count.merged), "merged"] = env$count.merged 
     env$track[names(env$count.nonchim), "nonchim"] = env$count.nonchim    

     # Convert each of the tracks into percents
     env$track[,"filteredNF"] = env$track[,"filteredNF"]/env$track[,"rawF"]
     env$track[,"filteredNR"] = env$track[,"filteredNR"]/env$track[,"rawR"]
     env$track[,"cutThreshF"] = env$track[,"cutThreshF"]/env$track[,"rawF"]
     env$track[,"cutThreshR"] =	env$track[,"cutThreshR"]/env$track[,"rawR"]
     env$track[,"denoisedF"] = env$track[,"denoisedF"]/env$track[,"rawF"]
     env$track[,"denoisedR"] = env$track[,"denoisedR"]/env$track[,"rawR"]
     env$track[,"merged"] = env$track[,"merged"]/env$track[,"rawF"]
     env$track[,"nonchim"] = env$track[,"nonchim"]/env$track[,"rawF"]
     env$track[,"rawF"] = env$track[,"rawF"]/env$track[,"rawF"]
     env$track[,"rawR"] = env$track[,"rawR"]/env$track[,"rawR"]

     ## F. Visualize read input / output
     env$track.m = melt(as.matrix(env$track))
     colnames(env$track.m) = c("sample", "dataClass", "reads")
     env$track.m$dataClass = factor(env$track.m$dataClass, levels = c("rawF", "rawR", "filteredNF", "filteredNR",
     			   "cutThreshF", "cutThreshR", "denoisedF", "denoisedR", "merged", "nonchim"))
     
     p = ggplot(data = env$track.m, aes(x = dataClass, y = sample, fill = reads)) + 
       	  geom_tile(aes(fill = reads)) +
     	  scale_fill_gradient(low = "white", high = "steelblue") +
  	  xlab("Pipeline Phase") + 
	  ylab("Sample") + 
	  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  	  guides(fill=guide_legend(title="%ReadsRetained"))
     ggsave(p, file = file.path(env$output, paste0(env$proj.name, "heatMapReads.pdf")))
 
    
     # Return env
     return(env)
     
}

# FUNCTION: dada2.makeCountTableForwardOnly
# INPUTS:
#     env: from dada2.makeASVTable
# OUTPUTS:
#     counts: a table of counts for each step in the filtering process
dada2.makeCountTableForwardOnly = function(env, fileTailPatternForward = ".R1.fastq.gz", fileTailPatternReverse = ".R2.fastq.gz"){
     
     
     # Count forward original reads
     my.cat("Counting the reads in original files...")

     # Count the number of sequences in the forward 
     env$count.fnFs = sapply(env$fnFs, getN)

     # Rename
     names(env$count.fnFs) = sub(x = basename(path = names(env$count.fnFs)), pattern = fileTailPatternForward, replace = "")

     # Count the sequences filtered for N's
     my.cat("Counting the reads filtered for Ns...")

     # Count the number of filtered N sequences in the forward and reverse reads
     env$count.fnFs.filtN = sapply(env$fnFs.filtN[file.exists(env$fnFs.filtN)], getN)

     # Rename the filtered N sequences
     names(env$count.fnFs.filtN) = sub(x = basename(path = names(env$count.fnFs.filtN)), pattern = fileTailPatternForward, replace = "")

     # Message to count filtered reads
     my.cat("Counting the numbers of filtered cut reads...")
     
     # Count the filtered cut sequences
     env$count.filtFs = sapply(env$filtFs[file.exists(env$filtFs)], getN)

     # Rename the filtered cut sequences
     names(env$count.filtFs) = sub(x = basename(path = names(env$count.filtFs)), pattern = fileTailPatternForward, replace = "")

     # Count dadaFs 
     my.cat("Counting the number of reads after dada2 algorithm...")
     
     # Count dada reads
     env$count.dadaFs = sapply(env$dadaFs, getN)

     # Rename dada reads
     names(env$count.dadaFs) = sub(x = basename(path = names(env$count.dadaFs)), pattern = fileTailPatternForward, replace = "") 

     # Count Nonchim
     my.cat("Counting the number of sequences after chimeras are removed...")
     env$count.nonchim = rowSums(env$seqtab.nochim)
     names(env$count.nonchim) = sub(x = basename(path = names(env$count.nonchim)), pattern = fileTailPatternForward, replace = "")

     # Make track
     env$track = data.frame(matrix(nrow = length(env$count.fnFs), ncol = 5, data = NA))
     colnames(env$track) = c("rawF", "filteredNF", "cutThreshF",  "denoisedF",  "nonchim")

     rownames(env$track) = names(env$count.fnFs)
     env$track[names(env$count.fnFs), "rawF"] = env$count.fnFs 
     env$track[names(env$count.fnFs.filtN), "filteredNF"] = env$count.fnFs.filtN 
     env$track[names(env$count.filtFs), "cutThreshF"] = env$count.filtFs 
     env$track[names(env$count.dadaFs), "denoisedF"] = env$count.dadaFs 
     env$track[names(env$count.nonchim), "nonchim"] = env$count.nonchim    

     # Convert each of the tracks into percents
     env$track[,"filteredNF"] = env$track[,"filteredNF"]/env$track[,"rawF"]
     env$track[,"cutThreshF"] = env$track[,"cutThreshF"]/env$track[,"rawF"]
     env$track[,"denoisedF"] = env$track[,"denoisedF"]/env$track[,"rawF"]
     env$track[,"nonchim"] = env$track[,"nonchim"]/env$track[,"rawF"]
     env$track[,"rawF"] = env$track[,"rawF"]/env$track[,"rawF"]

     ## F. Visualize read input / output
     env$track.m = melt(as.matrix(env$track))
     colnames(env$track.m) = c("sample", "dataClass", "reads")
     env$track.m$dataClass = factor(env$track.m$dataClass, levels = c("rawF", "filteredNF",
     			     				      "cutThreshF", "denoisedF", "nonchim"))
     
     p = ggplot(data = env$track.m, aes(x = dataClass, y = sample, fill = reads)) + 
       	   geom_tile(aes(fill = reads)) +
     	     scale_fill_gradient(low = "white", high = "steelblue") +
  	       xlab("Pipeline Phase") + 
	         ylab("Sample") + 
		   theme(axis.text.x=element_text(angle=90,hjust=1)) +
  		     guides(fill=guide_legend(title="%ReadsRetained"))
     ggsave(p, file = file.path(env$output, paste0(env$proj.name, "heatMapReads.pdf")))
 
    
     # Return env
     return(env)
     
}

# FUNCTION: assignTaxonomy
# INPUT: 
#     env: an environment from dada2.removeChimeras
# OUTUT:
#     taxa: the assigned taxonomy table
dada2.assignTaxonomy = function(env, cores){
      
     # Assign taxonomy table
     my.cat("Assigning taxonomy table...")
     env$taxa = assignTaxonomy(env$seqtab.nochim, env$silvaTax, multithread = cores)

     # Assign species table
     my.cat("Assigning species table...")
     env$taxa = addSpecies(env$taxa, env$silvaSpecies)

     # Print out tax table
     my.cat("Taxa assigned! ...")
     taxa.print = env$taxa
     rownames(taxa.print) = NULL
     print(head(taxa.print))
     
     # Return environment
     return(env)
}

# FUNCTION: dada2.makePhyloseq
# INPUT:
#    env: an environment variable from dada2.assignTaxonomy
# OUTPUT:
#    phylo: a phyloseq object with taxa table, metadata, and asvs
dada2.makePhyloseq = function(env){
     

     my.cat("Removing the following rows of the metadata table which no longer have samples:")
     print(rownames(env$metadata)[which(!rownames(env$metadata) %in% rownames(env$asv))])
     
     my.cat("Creating phyloseq object...")
     env$dada2Metadata = env$metadata[rownames(env$asv),]
     
     # Make components of the phyloseq object
     ASV = otu_table(as.matrix(env$asv), taxa_are_rows = FALSE)
     TAX = tax_table(as.matrix(env$taxa))
     META = sample_data(env$dada2Metadata)
     
     # Make the phyloseq object
     env$phyloseq = phyloseq(ASV, TAX, META)  

     # Print phyloseq object
     print(env$phyloseq)
     
     # Return env
     return(env)
}

# FUNCTION: dada2.save
# INPUT: 
#   env: an env object from dada2.makePhyloseq
# OUTPUT: 
#   Saves image of run and saves phyloseq object as a .RDS object
dada2.save = function(env){
     
     # Name
     baseName = paste0(c(env$proj.name, "_", env$process.id), sep = "", collapse = "")

     # Save Image
     fileImageName = file.path(env$output, paste0(baseName, "_dada2Run.RData"))
     my.cat(paste0("Saving image to... ", fileImageName))
     save.image(file = fileImageName, compress = TRUE) 
     
     # Save phyloseq
     phyloseqImageName = file.path(env$output, paste0(baseName, "_phyloseq.rds"))
     my.cat(paste0("Saving phyloseq object to... ", phyloseqImageName))
     saveRDS(file = phyloseqImageName, env$phyloseq)

     # Save ASV Table
     asvImageName = file.path(env$output, paste0(baseName, "asv.rds"))
     my.cat(paste0("Saving ASV table to ... ", asvImageName))
     saveRDS(file = asvImageName, env$asv)

     # Save TAX Table
     taxImageName = file.path(env$output, paste0(baseName, "tax.rds"))
     my.cat(paste0("Saving TAX table to ... ", taxImageName))
     saveRDS(file = taxImageName, env$tax)

}

# FUNCTION: dada2.begin
# PURPOSE: Run dada2 up to producting quality plots
# INPUTS:
#    cores: number of cores to run dada2
#    cutAdapt: number of cores to run cutAdapt
#    pool: should sequences be pooled for dada2
#    thresh: the number of reads for which to throw out a sample 
#    params: the user defined params
#    metadataPath: the path to the metadatafile
#    fastq.path: the path to the fastq files
#    error: the amount of mismatches allowed for cut adapt. Default 0.1
# OUTPUTS: Runs the dada2 script for the following functions
#    dada2.envPrep()
#    dada2.userDefinedParms(params, env)
#    dada2.preFilter(env)
#    dada2.checkMetadataMatchFasQ(metadataFile, fastq.path, env)
#    dada2.preFilter(env, cores)
#    dada2.checkForPrimersInitial(env)
#    dada2.cut(env, cores, error)
#    dada2.checkForPrimersPostCut(env)
#    dada2.printCutQualityPlots(env, thresh)
#    dada2.printFiltNQualityPlots(env, thresh)
dada2.begin = function(cores, coresAdapt, params, metadataPath, fastq.path, thresh = 1000, error = 0.1){
     

     # Prep details about the environment
     master = dada2.envPrep()

     # Set some user defined parameters
     master = dada2.userDefinedParams(file = params, env = master)

     # Check that the fastq filenames and the metadata files match
     master = dada2.checkMetadataMatchFastQ(metadataFile = metadataPath, fastq.path = fastq.path, env = master)

     # Prefilter the reads - remove Ns
     master = dada2.preFilter(env = master, cores = cores)

     # Check for primers up front
     dada2.checkForPrimersInitial(env = master)

     # Cut off primers
     master = dada2.cut(env = master, cores = coresAdapt, error = 0)

     # Check for primers after cut
     dada2.checkForPrimersPostCut(env = master) 

     # Print out the cut sequence quality plots
     dada2.printCutQualityPlots(env = master, thresh = thresh)

     # Print out the filtered quality plots
     master = dada2.printFiltNQualityPlots(env = master, thresh = thresh)
     
     # Return environment
     return(master)
}

# FUNCTION: dada2.finish
# PURPOSE: Run dada2 algorithm from quality plots to end
# INPUTS: 
#    env: an env from dada2.start
#    pool: should dada2 sequences be pooled?
#    thresh: the threshold to remove samples if below that threshold
#    fwdTrimLength: Forward read trim length
#    revTrimLength: The reverse read trim length
#    cores: The numbers of cores to run dada2
#    merge: set merge to false for consideration of forward reads only
# OUTPUTS:
#    dada2.filterCutSeqs: Filter the cut sequences
#    dada2.learnErrors: Learn errors from sequences
#    dada2.mergers: merge forward and reverse reads
#    dada2.sequenceTable: make sequence table from mergers
#    dada2.removeChimeras: Remove chimeras from sequence table
#    dada2.makeASVTable: ASV table from removeChiemras
#    dada2.makeCountTable: makes a figure of percent reads lost at each step
#    dada2.assignTaxonomy: assign taxonomy tables
#    dada2.makePhyloseq: make phyloseq environment
#    dada2.saveImages: save key parts of the dada2 algorithm
dada2.finish = function(env, pool = F, thresh = 1000, fwdTrimLength, revTrimLength = NA, cores, merge = T){

# Check if merge = T and no revTrimLength
if(merge & is.na(revTrimLength)){
     stop("Error: please provide revTrimLength to merge forward and reverse reads\n")
}


if(merge){     
     # Filter cut sequences
     master = dada2.filterCutSeqs(env = env, fwdTrimLength = fwdTrimLength, revTrimLength = revTrimLength, cores = cores, thresh = thresh)

     # Learn errors
     master = dada2.learnErrors(env = master, cores = cores, thresh = thresh) 
    
     # Run dada2
     master = dada2.dada2(env = master, pool = pool, cores = cores) 

     # Merge forward and reverse reads
     master = dada2.mergers(env = master)

     # Make sequence table
     master = dada2.sequenceTable(env = master)
     
     # Remove Chimeras
     master = dada2.removeChimeras(env = master, cores = cores)

     # Make ASV Table
     master = dada2.makeASVTable(env = master)

     # Make track table
     master = dada2.makeCountTable(env = master)
     
     # Assign Taxonomy
     master = dada2.assignTaxonomy(env = master, cores = cores)

     # Make Phyloseq object
     master = dada2.makePhyloseq(env = master)

     # Save
     dada2.save(env = master)

     # Get stats
     dada2.getStats(env = master)

     # Return
     return(master)
     }else{ # Run dada2 on forward seqs only

     # Filter and cut reads for forward sequences only
     master = dada2.filterCutSeqsForwardOnly(env = env, fwdTrimLength = fwdTrimLength, cores = cores, thresh = thresh)

     # Learn errors for forward reads only
     master = dada2.learnErrorsForwardOnly(env = master, cores = cores, thresh = thresh)

     # Run Dada2 on forward reads only
     master = dada2.dada2ForwardOnly(env = master, cores = cores, pool = pool)

     # Make sequence table from forward reads only
     master = dada2.sequenceTableForwardOnly(env = master)

     # Remove Chimeras
     master = dada2.removeChimeras(env = master, cores = cores)

     # Make ASV Table
     master = dada2.makeASVTable(env = master)

     # Make counts table
     master = dada2.makeCountTableForwardOnly(env = master)

     # Assign Taxonomy
     master = dada2.assignTaxonomy(env = master, cores = cores)

     # Make Phyloseq object
     master = dada2.makePhyloseq(env = master)

     # Save
     dada2.save(env = master)

     # Get Stats
     dada2.getStats(env = master)
}
}

# FUNCTION: runFigaro
# INPUTS: 
#    fastqDir: a directory with a set of files named with the convention
#         study_sampleName.R[12].fastq.gz
#    dir: the output directory where to make the soft links
# OUTPUTS:
#    outDir: an out directory with links to the original files
#         renamed in the zymo format of samName_16s_R[12].fastq.gz
runFigaro = function(fastqDir, env){

     # Make a path for links
     env$path.figaro = file.path(env$output, "figaroInDir")

     # Make a new directory if it doesn't exist already
     if(!dir.exists(env$path.figaro)) dir.create(env$path.figaro)

     # Make a list of files to loop through
     files = list.files(fastqDir, pattern = "fastq.gz")

     for(i in 1:length(env$fnFs.cut)){
          targetName = env$fnFs.cut[i]
	  

	  figaroFileName = basename(gsub(env$fnFs.cut[i], pattern = "_", replacement = ""))
	  figaroFileName = gsub(figaroFileName, pattern = ".R", replacement = "_16s_R")
	  figaroFileName = file.path(env$path.figaro, figaroFileName)

	  my.cat(paste0("Target Name: ", targetName))
	  my.cat(paste0("Figaro Name: ", figaroFileName))
          createLink(link = figaroFileName, target = targetName)     

}

  for(i in 1:length(env$fnRs.cut)){
          targetName = env$fnRs.cut[i]

          figaroFileName = basename(gsub(env$fnRs.cut[i], pattern = "_", replacement = ""))
          figaroFileName = gsub(figaroFileName, pattern = ".R", replacement = "_16s_R")
          figaroFileName = file.path(env$path.figaro, figaroFileName)

          my.cat(paste0("Target Name: ", targetName))
          my.cat(paste0("Figaro Name: ", figaroFileName))
          createLink(link = figaroFileName, target = targetName)

}
     system2(env$python3, args = c(env$figaro, "-i", env$path.figaro, "-o", env$output, "-f", 1, "-r", 1,
	     	     "-a", env$DIST, "-F zymo"), stdout = TRUE, wait = TRUE)
     

}

# FUNCTION: dada2.getTrackStats
# INPUTS: env from dada2 with env$track
# OUTPUTS:
#   trackStats.txt file with a set of track stats
dada2.getTrackStats = function(env){

    # Make a dataframe
    d = as.data.frame(matrix(ncol = 2, nrow = 4))
    colnames(d) = c("StatDescription", "Stat")
    d[1, 1] = "Ave Reads Retained"
    d[1, 2] = mean(env$track[,"nonchim"], na.rm = TRUE)
    d[2, 1] = "Min Reads Retained"
    d[2, 2] = min(env$track[,"nonchim"], na.rm = TRUE)
    d[3, 1] = "Max Reads Retained"
    d[3, 2] = max(env$track[,"nonchim"], na.rm = TRUE)
    d[4, 1] = "SD Reads Retained"
    d[4, 2] = sd(env$track[, "nonchim"], na.rm = TRUE)
    d[5, 1] = "Samples Removed"
    d[5, 2] = sum(is.na(env$track[, "nonchim"]))
    d[6, 1] = "Percent Samples Removed"
    d[6, 2] = sum(is.na(env$track[, "nonchim"])) / length(env$track[, "nonchim"])

    x = xtable(d)
    print(x, file = file.path(env$output, paste0(env$proj.name, "TrackStatsLatex.txt")))
    
}

# FUNCTION: getHistogramReadLengths
# INPUTS: env
# OUTPUTS: plots read length histograms
dada2.getHistogramReadLengths = function(env, row_names_are_samples = T){
  
  if(!row_names_are_samples){
    asv = t(env$asv)
  }else{
    asv = env$asv
  }
  
  # Get read lengths
  d = as.data.frame(nchar(colnames(asv)))
  colnames(d) = c("ReadLength")
  
  # Plot read lengths
  h = ggplot(d, aes(x = ReadLength)) + geom_histogram() + labs(x = "Read Length", y = "Count") + ggtitle("Read Length Distribution")
  
  ggsave(h, file = file.path(env$output, paste0(env$proj.name, "readLengthDistribution.pdf")))

}

# FUNCTION: dada2.getReadDepthHistogram
# INPUTS: env after running dada2
# 
dada2.getReadDepthHistogram = function(env, row_names_are_samples = T){
  
  if(!row_names_are_samples){
    asv = t(env$asv)
  }else{
    asv = env$asv
  }

  d = as.data.frame(apply(asv, 1, sum))
  colnames(d) = "ReadDepth"
  
  statsReadDepth = as.data.frame(matrix(ncol = 2, nrow = 3, data = NA))
  colnames(statsReadDepth) = c("StatDescription", "Stat")
  statsReadDepth[1, 1] = "AveReadDepth"
  statsReadDepth[1, 2] = ave(d$ReadDepth)[1]
  statsReadDepth[2, 1] = "MinReadDepth"
  statsReadDepth[2, 2] = min(d$ReadDepth)
  statsReadDepth[3, 1] = "MaxReadDepth"
  statsReadDepth[3, 2] = max(d$ReadDepth)
  write.table(file = file.path(env$output, "statsReadDepth.txt"), x = statsReadDepth, quote = F, sep = "\t")

  env$statsReadDepth = statsReadDepth
  

  xtableStats = xtable(env$statsReadDepth)
  print(xtableStats, file = file.path(env$output, paste0(env$proj.name, "statsReadDepthXtable.txt")))
  
  h = ggplot(d, aes(x = ReadDepth)) + geom_histogram() + labs(x = "Read Depth", y = "Count") + ggtitle("Read Depth Distribution")
  ggsave(h, file = file.path(env$output, paste0(env$proj.name, "readDepthDistribution.pdf")))
  
}



# FUNCTION: dada2.getCollectorsCurveBySample
# INPUTS: An environment variable output from dada2 pipeline
# OUTPUTS: A pdf containing collectors curves
dada2.getCollectorsCurveBySample = function(env, row_names_are_samples = T){
 
  if(!row_names_are_samples){
    asv = t(env$asv)
  }else{
    asv = env$asv
  }
  
  maxReadDepth = max(apply(asv, 1, sum))
  rare = seq(from = 100, to = maxReadDepth, length.out = 10)
  rare = floor(rare)  

  d = data.frame(matrix(nrow = nrow(asv), ncol = length(rare)))
  rownames(d) = rownames(asv)
  colnames(d) = paste0("r", rare)
  
  for(i in 1:length(rare)){
    #print(i/length(rare))
    curASVRare = rarefy_even_depth(env$phyloseq, sample.size = rare[i], rngseed = 3, replace = F, trimOTUs = F, verbose = F)
    curASVRare = otu_table(curASVRare)
    
    for(j in 1:nrow(curASVRare)){
      
      curSamp = rownames(curASVRare)[j]
      d[curSamp, i] = length(which(curASVRare[curSamp, ] >0))
      
    }
   
  }
  
  d$id = rownames(d)
  d.melt = melt(d)
  
  p = ggplot(d.melt, aes(x = variable, y = value, group = id)) +
    geom_line()+
    labs(title = "Collector's Curve", x = "Sequencing Depth", y = "# Unique ASV") + 
    theme(axis.text = element_text(size = 18), axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 22), axis.title = element_text(size = 22))
  ggsave(p, file = file.path(env$output, paste0(env$proj.name, "CollectorsCurves.pdf")))
  
  return(d)
}


# FUNCTION dada2.getStats
# INPUTS: env environment from dada2.finish
# OUTPUTS: calls the following
#   dada2.getCollectorsCurveBySample
#   dada2.getTrackStats
#   dada2.getReadDepthHistogram
#   dada2.getHistogramReadLengths
dada2.getStats = function(env){

     
     # Get Collectors curve
     my.cat("Making collector's curves...")
     dada2.getCollectorsCurveBySample(env = env)

     # Print track stats
     my.cat("Getting readloss statistics...")
     dada2.getTrackStats(env = env)     

     # Get Read depth histogram
     my.cat("Making a histogram of read depth...")
     dada2.getReadDepthHistogram(env = env)

     # Get Read Length Histogram
     my.cat("Making a histogram of read lengths...")
     dada2.getHistogramReadLengths(env = env)

     # Get Run Variables
     my.cat("Printing output of run variables...")
     dada2.printRunVars(env = env)

}

# FUNCTION dada2.printRunVars
# INPUTS: env environment from dada2.finish
# OUTPUTS: a text file with stats on the dada2 run
dada2.printRunVars = function(env){

     d = as.data.frame(matrix(ncol = 2, nrow = 17))
     colnames(d) = c("RunVariable", "Value")

     d[1, 1] = "CutAdaptPath"
     d[1, 2] = env$cutAdapt
     d[2, 1] = "SilvaTax"
     d[2, 2] = env$silvaTax
     d[3, 1] = "metadataFilePath"
     d[3, 2] = env$metadataFile
     d[4, 1] = "ForwardPrimer"
     d[4, 2] = env$FWD
     d[5, 1] = "ReversePrimer"
     d[5, 2] = env$REV
     d[6, 1] = "DIST"
     d[6, 2] = env$DIST
     d[7, 1] = "fastq.path"
     d[7, 2] = env$fastq.path
     d[8, 1] = "outputDir"
     d[8, 2] = env$output
     d[9, 1] = "path.cut"
     d[9, 2] = env$path.cut
     d[10, 1] = "pathQualityCut"
     d[10, 2] = env$path.qualCut
     d[11, 1] = "path.qualFiltN"
     d[11, 2] = env$path.qualFiltN
     d[12, 1] = "processDate"
     d[12, 2] = env$process.date
     d[13, 1] = "processID"
     d[13, 2] = env$process.id
     d[14, 1] = "projectName"
     d[14, 2] = env$proj.name
     d[15, 1] = "R1 Flags"
     d[15, 2] = env$R1.flags
     d[16, 1] = "R2 Flags"
     d[16, 2] = env$R2.flags
     d[17, 1] = "SilvaSpecies"
     d[17, 2] = env$silvaSpecies


     write.table(d, file = file.path(env$output, paste0(env$proj.name, "runSessionStats.txt")), quote = F, row.names = F)

     dLatex = xtable(d)
     
     print(dLatex, file = file.path(env$output, paste0(env$proj.name, "runSessionStatsLatex.txt")))
     
     return(d)
}

# FUNCTION: metapipes.humanReadable
# INPUTS:
#   env: an environment passed from dada2.getStats, dada2.finish, dada2.begin
# OUTPUTS:
#   asvHR.rds # asv table which is in human readable format
#   taxHR.rds # tax table which is in human readable dada2 format
#   taxOTUHR.rds # tax table which is in human readable QIIME format
#   seqKeyHR.rds # sequence key which is in human readable format
#   asvs.fastq $ fastq file containing all asvs.
metapipes.humanReadable = function(env){

    # Rename ASV and Tax table
    env$phyloseqHR = metapipes.renameASVs(physeq = env$phyloseq)
    env$sequenceKey = data.frame("seqID" = env$phyloseqHR$sequenceKey.seqID, "ASV" = env$phyloseqHR$sequenceKey.ASV)
    env$psHR = phyloseq(otu_table(env$phyloseqHR$phyloseq), tax_table(env$phyloseqHR$phyloseq), sample_data(env$metadata))
    env$path.sequenceKeyfasta = file = file.path(env$output, "seqKeyHR.fasta")
    # Write out the fastq file
    metapipes.writeFasta(key = env$sequenceKey, file = file.path(env$output, "seqKeyHR.fasta"))

    # Get the human readible QIIME tax table format
    env$taxHRQiimeFormat = metapipes.getTaxString(taxTable = tax_table(env$phyloseqHR$phyloseq))

    # Base name for files
    baseName = paste0(c(env$proj.name, "_", env$process.id), sep = "", collapse = "")

    #Save phyloseq human readable
    phyloseqImageName = file.path(env$output, paste0(baseName, "_phyloseqHumanReadable.rds"))
    my.cat(paste0("Saving human readable phyloseq object to...", phyloseqImageName ))
    saveRDS(file = phyloseqImageName, env$psHR)

    # Save ASV table
    asvImageName = file.path(env$output, paste0(baseName, "asvHumanReadable.rds"))
    my.cat(paste0("Saving human readable ASV table to ... ", asvImageName))
    saveRDS(file = asvImageName, otu_table(env$phyloseqHR$phyloseq))

    # Save TAX table dada2 format
    taxImageName = file.path(env$output, paste0(baseName, "taxHumanReadable.rds"))
    my.cat(paste0("Saving tax table to ...", taxImageName))
    saveRDS(file = taxImageName, tax_table(env$phyloseqHR$phyloseq))

    # Save TAX table in QIIME format
    taxQIIMEImageName = file.path(env$output, paste0(baseName, "taxQIIMEHumanReadable.rds"))
    my.cat(paste0("Saving QIIME formatted tax table to ...", taxQIIMEImageName))
    saveRDS(file = taxQIIMEImageName, env$taxHRQiimeFormat)

    # Save sequence key
    sequenceKeyName = file.path(env$output, paste0(baseName, "sequenceKeyHumanReadable.rds"))
    my.cat(paste0("Saving sequence key object to...", sequenceKeyName))
    saveRDS(file = sequenceKeyName, env$sequenceKey)
    
    return(env)
}



# FUNCTION: metapipes.renameASVs
# INPUT:
#   physeq: a phyloseq object with tax and ASV table minimum
# OUTPUT:
#   returns a list with renamed phyloseq object and the sequence key
#      sequenceKey.seqID
#      sequenceKey.ASV
#      phyloseq
metapipes.renameASVs = function(physeq){

         N = nrow(tax_table(physeq))

         # Make Sequence Key
         sequenceKey = as.data.frame(matrix(nrow = N, ncol = 2))
         colnames(sequenceKey) = c("seqID", "ASV")
         rownames(sequenceKey) = rownames(tax_table(physeq))
         sequenceKey$seqID = paste0("seq", seq(from = 1, to = N, by = 1))
         sequenceKey$ASV = rownames(sequenceKey)

         # Get ASV table where taxa are rows
         asv = meta.getOTUTable(physeq)
         tax = tax_table(physeq)

         # For every row
         for(i in 1:N){

                # Get the current ASV to rename
                curASV = rownames(sequenceKey)[i]

                # Rename the ASV rows
                curASVIdx = which(rownames(asv) == curASV)
                rownames(asv)[curASVIdx] = as.vector(sequenceKey[curASV, 1])

                # Rename the tax rows
                curTaxIdx = which(rownames(tax) == curASV)
                rownames(tax)[curTaxIdx] = as.vector(sequenceKey[curASV, 1])

                  }

         # Rename rownmaes of sequenceKey to match tax table and asv Table
         rownames(sequenceKey) = sequenceKey$seqID
             
         # Return relabeled phyloseq object
         ps = phyloseq(OTU = otu_table(asv, taxa_are_rows = T), TAX = tax_table(tax))

    
         return(list = c("phyloseq" = ps, "sequenceKey" = as.data.frame(sequenceKey)))

}

## FUNCTION: metapipes.writeFasta
## PURPOSE: To convert a sequence table to a fastQ format
## INPUTS:
##   key: a sequence key file with the following requirements
##     Col 1: seqID the sequence ID (short version)
##     Col 2: sequence - The sequence itself (dada2 output)
##   file: the file name that you want to write the fastQ to
## OUTPUTS:
##   fastQ file written within your current working directory with the name provided to file
metapipes.writeFasta = function(key, file){

    # Create a sink to filename
    sink(file)

    # For each sequence, write a fastq entry
    for(i in 1:nrow(key)){
        cur = paste0(c(">", as.vector(key[i,"seqID"]), "\n", as.character(key[i,"ASV"]), "\n"), sep= "", collapse = "")
        cat(cur)
    }

    #Close connection
    closeAllConnections()
}

##FUNCTION: metapipes.getTaxString
## PURPOSE: To convert a dada2 formatted taxonomic table to a qiime formated
##          taxonomic table
## INPUTS   taxTable:
##          columnNames must be titled Kindom, Phylum, Class, Order, Family, Genus
##          rows correspond to each sequence.
##          An entry taxTable[i,j] corresponds to the jth taxonomic classifier
##          of sequence i.
## OUTPUT   taxReturn
##          A matrix with two colums
##          column1: the sequence ID from the rows in TaxTable
##          column2: the collapsed taxString in QIIME format
##          if the taxonomic entry in taxTable is NA, then the output string will
##          contain "Unclassified".
metapipes.getTaxString = function(taxTable){

    taxReturn = as.data.frame(matrix(nrow = nrow(taxTable), ncol= 2))
    for (i in 1:nrow(taxTable)){
        string = ""

        #Get Kingdom
        if(!is.na(taxTable[i, "Kingdom"])){
            string = paste(c(string, "k__", toString(taxTable[i, "Kingdom"])), collapse = "")
        }else{
            #string = paste(c(string, " Unassigned;"), collapse = "")
        }

        #Get Phylum
        if(!is.na(taxTable[i, "Phylum"])){
            string = paste(c(string, "; p__", toString(taxTable[i, "Phylum"])), collapse = "")
        }else{
                                        #string = paste(c(string, " Unassigned;"), collapse = "")
        }

        #Get Class
        if(!is.na(taxTable[i, "Class"])){
            string = paste(c(string, "; c__", toString(taxTable[i, "Class"])), collapse = "")
        }else{
                                          #string = paste(c(string, " Unassigned;"), collapse = "")
        }

        #Get Order
        if(!is.na(taxTable[i, "Order"])){
            string = paste(c(string, "; o__", toString(taxTable[i, "Order"])), collapse = "")
        }else{
                                        #string = paste(c(string, " Unassigned;"), collapse = "")
        }
 
        #Get Family
        if(!is.na(taxTable[i, "Family"])){
            string = paste(c(string, "; f__", toString(taxTable[i, "Family"])), collapse = "")
        }else{
                                        #string = paste(c(string, " Unassigned;"), collapse = "")
        }

        #Get Genus
        if(!is.na(taxTable[i, "Genus"])){
            string = paste(c(string, "; g__", toString(taxTable[i, "Genus"])), collapse = "")
        }else{
                                        #string = paste(c(string, " Unassigned;"), collapse = "")
        }

        if(!is.na(taxTable[i, "Kingdom"])){
            if(taxTable[i, "Kingdom"] == "No blast hit"){
                print("Unassigned")
                string = "Unassigned"
                taxReturn[i,1] = rownames(taxTable)[i]
                taxReturn[i,2] = string
 
            }
        }
        taxReturn[i,1] = rownames(taxTable)[i]
        taxReturn[i,2] = string
    }

    # Return tax table
    colnames(taxReturn) = c("sequence", "tax")
    return(taxReturn)
}

# FUNCTION: metapipes.mothurAlign
# INPUTS
#   env from metapipes.humanReadable
# OUTPUTS
#   returns an env which has paths to produced alignment files. 
metapipes.mothurAlign = function(env){

    cmd = paste(c("mothur \"#align.seqs(candidate=", env), collapse = "", sep = "")
    system2()
    
}

#FUNCTION clipASVSeqsTable
# PURPOSE: Clips ASVs to specified length, then merges now no longer
#   ASVs into a single row / column in the matrix as well as tax labels
# INPUTS:
#  ps: A phyloseq object
#  length: The length to clip sequences to
# OUTPUTS:
#  a phyloseq object

meta.clipASVSeqsTable = function(ps, length){
  
  # Get ASV and tax Table to work with
  asv = otu_table(ps)
  if(taxa_are_rows(ps)){asv = t(asv)}else{}
  tax = tax_table(ps)
  
  # Determine clip length is appropriate for the min seq length
  minSeqLength = as.numeric(min(names(table(nchar(colnames(otu_table(asv)))))))
  print(minSeqLength)
  if(minSeqLength < length){
    
    my.cat("Warning: some of the sequences are less than desired clip length")
    
  }
  
  # Now get the clipped sequences
  my.cat("Making clipped ASV table...")
  asvNames = unique(str_sub(colnames(asv), start = 1, end = length))
  asvNew = as.data.frame(matrix(nrow = nrow(asv), ncol = length(asvNames), data = 0))
  rownames(asvNew) = rownames(asv)
  colnames(asvNew) = asvNames
  
  
  for(i in 1:nrow(asv)){
    
    curSamp = rownames(asv)[i]
    
    for(j in 1:ncol(asv)){
      
      curASV = str_sub(colnames(asv)[j], start = 1, end = length)	   
      asvNew[curSamp, curASV] = asvNew[curSamp, curASV] + asv[i, j]
      
    }
  }
  
  my.cat("Making clipped TAX table...")
  taxNames = asvNames
  taxNew = as.data.frame(matrix(nrow = length(taxNames), ncol = ncol(tax)))
  rownames(taxNew) = taxNames 
  colnames(taxNew) = c("Kindom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  
  for(i in 1:nrow(tax)){
    
    curASV = rownames(tax)[i]
    curASV = str_sub(curASV, start = 1, end = length)
    
    if(is.na(taxNew[curASV, "Kindom"])){
      taxNew[curASV, ] = as.vector(tax[i,])
    }else{
      taxNew[curASV, ] = getTaxString(t1 = as.vector(tax[i,]), t2 = as.vector(taxNew[curASV,]))
    }
    
  }
  
  my.cat("Making phyloseq object...")
  
  clippedPS = phyloseq(ASV = as.matrix(otu_table(asvNew, taxa_are_rows = F)), TAX = tax_table(as.matrix(taxNew)))
  return(clippedPS)
}

getTaxString = function(t1, t2){
  
  N = length(t1)
  
  if(N != length(t2)){
    stop("Cannot compare tax strings of unequal length")
  }
  
  matched = vector()
  for(i in 1:N){
    
    # Case 1: t1[i] is not present
    if(is.na(t1[i])){ 
      
      for(j in i:N){matched = c(matched, NA)}
      return(matched)
      
      
    }
    
    # Case 2: t2[i] is not present
    else if(is.na(t2[i])){
      
      for(j in i:N){matched = c(matched, NA)}
      return(matched)
      
    }
    
    # Case 3: Both t1[i] and t2[i] are present and match
    else if(t1[i] == t2[i]){
      
      matched = c(matched, t1[i])
      
    }
    
    # Case 4: Both t1[i] and t2[i] are present but don't match
    else{
      
      for(j in i:N){matched = c(matched, NA)}		
      return(matched)
      
    }
  }
  
  return(matched)
}


# LIBRARIES
#library(phyloseq)
#library(stringr)
#library(ggplot2)
#setwd(directory.data)
#psOld = readRDS("Kundu_2021_dada2_1.16.0_2021-02-15_phyloseq.rds")

# SOURCE
#setwd(directory.scripts)
#source("functionDebug.R")

#asvOld = otu_table(psOld)
#taxOld = tax_table(psOld)

# GLOBAL
#LENGTH = 240

#psNew = meta.clipASVSeqsTable(ps = psOld, length = LENGTH)
#test.meta.clipASVSeqsTable.asv(nPS = psNew, oPS = psOld, length = LENGTH)
#debug(test.meta.clipASVSeqsTable.asv)
#psFAIL = psNew
#otu_table(psFAIL)[1,1]
#otu_table(psFAIL)[1, 1] = 100000000
#otu_table(psFAIL)[1,1]
#test.meta.clipASVSeqsTable.asv(nPS = psFAIL, oPS = psOld, length = LENGTH)

# FUNCTION: test.meta.clipASVSeqsTable.asv
# INPUTS: 
#   nPS: a merged phyloseq object from meta.clipASVSeqsTable
#   oPS: an original phyloseq object provided to input to meta.clipASVSeqsTable
# OUTPUTS: 
#   will print nothing if two entries are equivalent
#   will print error if two values are not identical.
test.meta.clipASVSeqsTable.asv = function(nPS, oPS, length){
  
  # Get OTU tables with taxa as rows
  asvNew = meta.getOTUTable(ps = nPS)
  asvOld = meta.getOTUTable(ps = oPS)
    
  # Get duplicated and single names
 
  for(i in 1:ncol(asvOld)){
    print(i / ncol(asvOld))  
    curSampName = colnames(asvOld)[i]
    
    for(j in 1:nrow(asvOld)){
      
      curASVold = rownames(asvOld)[j]
      curASVnew = str_sub(curASVold, start = 1, end = length)
      
      entryOld = rownames(asvOld)[str_which(rownames(asvOld), pattern = paste0("^",curASVnew))]
      entryOld = sum(asvOld[entryOld, curSampName])
      entryNew = asvNew[curASVnew, curSampName]
      
      if(entryOld != entryNew){
       print(paste0("Testing entry i=", i, " j= ", j, ". Entry Old =", entryOld, ". Entry New = ", entryNew, "."), sep = "", collapse = "")
        stop ("ERROR")
      }else{

        
      }
      
    }
    
  }
  
}

# FUNCTION: test.meta.clipASVSeqsTable.tax
# INPUTS: 
#   nPS: a merged phyloseq object from meta.clipASVSeqsTable
#   oPS: an original phyloseq object provided to input to meta.clipASVSeqsTable
# OUTPUTS: 
#   will print nothing if two entries are equivalent
#   will print error if two values are not identical.
#   will print duplicated tax entries for user to check 
test.meta.clipASVSeqsTable.tax = function(nPS, oPS, length){

  newTax = tax_table(nPS)
  oldTax = tax_table(oPS)
  
  for(i in 1:nrow(oldTax)){
    #print(i/nrow(oldTax))
    curASVold = rownames(oldTax)[i]
    curASVnew = str_sub(curASVold, start = 1, end = length)
    
    entryOld = rownames(oldTax)[str_which(rownames(oldTax), pattern = paste0("^",curASVnew))]
    entryOld = oldTax[entryOld, ]
    entryNew = newTax[curASVnew, ]
    
    
    if(nrow(entryOld) == 1){
      
      if(identical(as.vector(entryOld), as.vector(entryNew))){
        
      }else{
        print(paste0("Testing entry i=", i, ". Entry Old =", entryOld, ". Entry New = ", entryNew, "."), sep = "", collapse = "")
        stop("ERROR")
      }
      
    }else{
      print("Old entry: ")
      print(entryOld)
      print("New entry: ")
      print(entryNew)
    }
  }
}
#debug(test.meta.clipASVSeqsTable.tax)
#test.meta.clipASVSeqsTable.tax(nPS = psNew, oPS = psOld, length = LENGTH)
#psFAIL = psNew
#tax_table(psFAIL)[1, 1]
#tax_table(psFAIL)[1, 1] = "ALIEN"
#tax_table(psFAIL)[1, 1]
#test.meta.clipASVSeqsTable.tax(nPS = psFAIL, oPS = psOld, length = LENGTH)


# FUNCTION: meta.getOTUTable
# INPUTS: Phyloseq object
# OUTPUTS: return OTU table where the taxa are rows and samples are columsn
meta.getOTUTable = function(ps){
  
  if(taxa_are_rows(ps)){
    return(otu_table(ps))
  }else{
    return(t(otu_table(ps)))
  }
}

# FUNCTION: meta.getCollectorsCurveBySample
# INPUTS: A phyloseq object
# OUTPUTS: A pdf containing collectors curves
meta.getCollectorsCurveBySample = function(ps, title){

  asv = meta.getOTUTable(ps = ps)

  maxReadDepth = max(apply(asv, 2, sum))
  rare = seq(from = 100, to = maxReadDepth, length.out = 20)
  rare = floor(rare)

  d = data.frame(matrix(nrow = ncol(asv), ncol = length(rare)))
  rownames(d) = colnames(asv)
  colnames(d) = paste0("r", rare)

  for(i in 1:length(rare)){
    #print(i/length(rare))
    curASVRare = rarefy_even_depth(ps, sample.size = rare[i], rngseed = 3, replace = F, trimOTUs = F, verbose = F)
    curASVRare = otu_table(curASVRare)

    for(j in 1:nrow(curASVRare)){

      curSamp = rownames(curASVRare)[j]
      d[curSamp, i] = length(which(curASVRare[curSamp, ] >0))

    }

  }

  d$id = rownames(d)
  d.melt = melt(d)

  p = ggplot(d.melt, aes(x = variable, y = value, group = id)) +
    geom_line()+
    labs(title = paste0(title, " Collector's Curve"), x = "Sequencing Depth", y = "# Unique ASV") +
    theme(axis.text = element_text(size = 18), axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 22), axis.title = element_text(size = 22))
  ggsave(p, file = paste0(title, "CollectorsCurves.pdf"))

  return(d)
}




# FUNCTION: meta.rarefyEvenDepth
# INPUTS: 
#   physeq: the phyloseq object with ASV and tax table minimum
#   sample.size: Defaults to the minimum sequencing depth
#   rngseed: random integer so results are reproducible. Record this. Defaults to 3
#   replace: Should ASVs be sampled with replacement? Default is FALSE
#   trimOTUs: Should OTUs be trimmed from the ASV table if they are not present after rarefaction in any sample? Default TRUE
#   verbose: Print out things? Default TRUE
# OUTPUTS
#   returns a phyloseq object that has been rarefied.
meta.rarefyEvenDepth = function(physeq, sample.size = min(sample_sums(physeq)), rngseed = 3, replace = FALSE, trimOTUs = TRUE, verbose = TRUE){

    #Run rarefy  in phyloseq
    ps.rare = rarefy_even_depth(physeq = physeq, sample.size = sample.size, rngseed = rngseed, replace = replace, trimOTUs = trimOTUs, verbose = verbose)

    return(ps.rare)
}

# FUNCTION: meta.renameASVs
# INPUT:
#   physeq: a phyloseq object with tax and ASV table minimum
#   name: the name of the sequence key table
# OUTPUT:
#   physeq: with ASVs relabeled for ease of use seq1, ..., seqN
#   name: name to write out sequence key
meta.renameASVs = function(physeq, name){
     
     N = nrow(tax_table(physeq))

     # Make Sequence Key
     sequenceKey = as.data.frame(matrix(nrow = N, ncol = 2))
     colnames(sequenceKey) = c("seqID", "ASV")
     rownames(sequenceKey) = rownames(tax_table(physeq))
     sequenceKey$seqID = paste0("seq", seq(from = 1, to = N, by = 1))
     sequenceKey$ASV = rownames(sequenceKey)
     
     # Get ASV table where taxa are rows
     asv = meta.getOTUTable(physeq)     
     tax = tax_table(physeq)

     # For every row
     for(i in 1:N){
     	   
	   # Get the current ASV to rename
	   curASV = rownames(sequenceKey)[i]

	   # Rename the ASV rows
	   curASVIdx = which(rownames(asv) == curASV)
	   rownames(asv)[curASVIdx] = as.vector(sequenceKey[curASV, 1])

	   # Rename the tax rows
	   curTaxIdx = which(rownames(tax) == curASV)
	   rownames(tax)[curTaxIdx] = as.vector(sequenceKey[curASV, 1])
     
     }

     # Rename rownmaes of sequenceKey to match tax table and asv Table
     rownames(sequenceKey) = sequenceKey$seqID


     # Write out sequence Key
     write.table(file = paste0(name, "sequenceKey.txt"), quote = F, sep = "\t", x = sequenceKey)
         
     # Return relabeled phyloseq object
     ps = phyloseq(OTU = otu_table(asv, taxa_are_rows = T), TAX = tax_table(tax))
     return(ps)
     
}

## FUNCTION: writeFastQ
## PURPOSE: To convert a sequence table to a fastQ format
## INPUTS: 
##   t: a sequence key file with the following requirements
##     Col 1: seqID the sequence ID (short version)
##     Col 2: sequence - The sequence itself (dada2 output)
##   file: the file name that you want to write the fastQ to
## OUTPUTS:
##   fastQ file written within your current working directory with the name provided to file
meta.writeFastQ = function(t, file){
  sink(file)
  for(i in 1:nrow(t)){
    cur = paste0(c(">", as.vector(t[i,"seqID"]), "\n",
                   as.character(t[i,"ASV"]), "\n"), sep= "", collapse = "")
    cat(cur)
  }
  closeAllConnections()
}



##FUNCTION: getTaxString
## PURPOSE: To convert a dada2 formatted taxonomic table to a qiime formated 
##          taxonomic table
## INPUTS   taxTable:
##          columnNames must be titled Kindom, Phylum, Class, Order, Family, Genus
##          rows correspond to each sequence. 
##          An entry taxTable[i,j] corresponds to the jth taxonomic classifier
##          of sequence i. 
## OUTPUT   taxReturn
##          A matrix with two colums
##          column1: the sequence ID from the rows in TaxTable
##          column2: the collapsed taxString in QIIME format 
##          if the taxonomic entry in taxTable is NA, then the output string will
##          contain "Unclassified".
getTaxString = function(taxTable){
  
  taxReturn = as.data.frame(matrix(nrow = nrow(taxTable), ncol= 2))
  
  for (i in 1:nrow(taxTable)){
    string = ""
    
    
    
    #get kingdom
    if(!is.na(taxTable[i, "Kingdom"])){
      string = paste(c(string, "k__", toString(taxTable[i, "Kingdom"])), collapse = "")
    }else{
      #string = paste(c(string, " Unassigned;"), collapse = "")
    }
    
    #get Phylum
    if(!is.na(taxTable[i, "Phylum"])){
      string = paste(c(string, "; p__", toString(taxTable[i, "Phylum"])), collapse = "")
    }else{
      #string = paste(c(string, " Unassigned;"), collapse = "")
    }
    
    #Get Class
    if(!is.na(taxTable[i, "Class"])){
      string = paste(c(string, "; c__", toString(taxTable[i, "Class"])), collapse = "")
    }else{
      #string = paste(c(string, " Unassigned;"), collapse = "")
    }
    
    #Get Order
    if(!is.na(taxTable[i, "Order"])){
      string = paste(c(string, "; o__", toString(taxTable[i, "Order"])), collapse = "")
    }else{
      #string = paste(c(string, " Unassigned;"), collapse = "")
    }
    
    #Get Family
    if(!is.na(taxTable[i, "Family"])){
      string = paste(c(string, "; f__", toString(taxTable[i, "Family"])), collapse = "")
    }else{
      #string = paste(c(string, " Unassigned;"), collapse = "")
    }
    
    #Get Genus
    if(!is.na(taxTable[i, "Genus"])){
      string = paste(c(string, "; g__", toString(taxTable[i, "Genus"])), collapse = "")
    }else{
      #string = paste(c(string, " Unassigned;"), collapse = "")
    }
    
    
    if(!is.na(taxTable[i, "Kingdom"])){
      if(taxTable[i, "Kingdom"] == "No blast hit"){
        print("Unassigned")
        string = "Unassigned"
        taxReturn[i,1] = rownames(taxTable)[i]
        taxReturn[i,2] = string
        
      }
    }
    taxReturn[i,1] = rownames(taxTable)[i]
    taxReturn[i,2] = string
  }
  
  
  colnames(taxReturn) = c("sequence", "tax")
  return(taxReturn)
}

