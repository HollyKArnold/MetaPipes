
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
     if(taxa_are_rows(ps)){}else{asv = t(asv)}
     tax = tax_table(ps)     		    
     
     # Determine clip length is appropriate for the min seq length
     minSeqLength = as.numeric(min(names(table(nchar(colnames(otu_table(asv)))))))
     if(minSeqLength < length){

          my.cat("Warning: some of the sequences are less than desired clip length")

     }


     return(clippedPS)
}