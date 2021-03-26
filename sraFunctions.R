
# AUTHOR: ARNOLD
# DAY: Wednesday, March 24, 2021
# DATE: 20210324
# PURPOSE: This file contains functions for prepping SRA data before running dada2




# FUNCTION: sra.prepMetadataTable
# INPUTS:
#  metadata: the SRA metadata table unmodified, downoladed from the site
# OUTPUTS:
#  sra.metadata.txt: A tab separated file where the rownmaes match the sample ERR or SRR numbers from fastq-dump
sra.prepMetadataTable = function(metadata, name){

    m = read.csv(metadata)
    name = paste0(name, "_")
    rownames(m) = paste0(name, as.vector(m$Run))
    write.table(file = paste(c(name, "_", "metadata.txt"), sep = "", collapse = ""), m, sep = "\t", quote = F)
    
    
}


