#' Title: Load Single End reads
#'
#' This function takes a fastq files with  PAIRED-END READS mapping to HXB2 genome with Kraken2 and loads it into the workspace as a dataframe
#'
#' @param metadata A table with sample meta data requiered to run virKra
#' @param library_se A string indicating data type (e.g.: "Chromatin Accessibility")
#' @param sample Sample id (same as orig.ident from Seurat)
#' @return A data.table with reads information
#' @examples
#' loadPairedEnd(metadata, library_pe = "Chromatin Accessibility", sample = "sample_id")
#' @importFrom data.table fread
#' @export

loadPairedEnd <- function(metadata, library_pe = c("Gene Expression", "Chromatin Accessibility"), sample = NULL){
  
  library_pe <- match.arg(library_pe, choices = c("Gene Expression", "Chromatin Accessibility"))
  
  # load metadata
  metadata <- metadata[library_type == library_pe, ]
  
  if (!is.null(sample)) {
    sample <- unique(metadata$sample)
  }
  
  files_r1 <- if (!is.null(metadata$read1_path)) metadata$read1_path else stop("\nERROR: missing R1 FASTQ file. ")
  files_r2 <- if (!is.null(metadata$read2_path)) metadata$read2_path else stop("\nERROR: missing R2 FASTQ file. ")
  
  # Transform read1_path and read2_path into a single column pe_reads_path
  metadata_long <- melt(metadata, id.vars = setdiff(names(metadata), c("read1_path", "read2_path")),
                        measure.vars = c("read1_path", "read2_path"), value.name = "pe_reads_path")
  data.table::setnames(metadata_long, "variable", "read_pair")
  metadata_long$read_pair <- gsub("read", "R", metadata_long$read_pair) %>% gsub("_path", "", .)
  
  # extract names and sequences of reads assigned to HXB2 genomes
  reads_r1 <- lapply(files_r1, .seq_info) # here we generate a list of tables
  names(reads_r1) <- files_r1 #basename(files)
  reads_r2 <- lapply(files_r2, .seq_info) # here we generate a list of tables
  names(reads_r2) <- files_r2 #basename(files)
    
  # remove empty fastq files_r1 from the list
  reads_r1 <- reads_r1[!unlist(lapply(reads_r1, is.null))]
  reads_r2 <- reads_r2[!unlist(lapply(reads_r2, is.null))]
  
  reads_pe <- c(reads_r1, reads_r2)
  #remove(reads_r1, reads_r2); gc()
  
  # no hiv+ reads found
  if(length(reads_pe) == 0){
    cat("\nWARNING!: any read has been reported as infected\n")
    return(NULL)
  } 
  
  # extract hiv+ identified reads
  reads_table <- data.table::rbindlist(
    lapply(names(reads_pe), function(file_name) 
      # extract fastq content
      sapply(reads_pe[file_name], function(record){
        data.frame(
          sample_file = file_name,
          read_name = record["read_name"],
          seq = record["seq"],
          quality = record["quality"])}
        )))
  
  data.table::setnames(reads_table, names(reads_table), c("sample_file", "read_name", "seq", "quality"))
  
  reads_table <- merge(reads_table, metadata_long[,.(library_type, sample, pe_reads_path, read_pair, index_path, cellranger_output, whitelist)], by.x = "sample_file", by.y = "pe_reads_path", all.x = TRUE)
  
  return(reads_table)
}
