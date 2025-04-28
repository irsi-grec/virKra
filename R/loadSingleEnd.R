#' Title: Load Single End reads
#'
#' This function takes a fastq files with  SINGLE-END READS mapping to HXB2 genome with Kraken2 and loads it into the workspace as a dataframe
#'
#' @param metadata A table with sample meta data requiered to run virKra
#' @param library_se A string indicating data type (e.g.: "Gene Expression")
#' @param sample Sample id (same as orig.ident from Seurat)
#' @return A data.table with reads information
#' @examples
#' loadSingleEnd(metadata, library_se = "Gene Expression", sample = "sample_id")
#' @importFrom data.table fread
#' @export

loadSingleEnd <- function(metadata, library_se = c("Gene Expression", "Chromatin Accessibility"), sample = NULL){

  library_se <- match.arg(library_se, choices = c("Gene Expression", "Chromatin Accessibility"))
    
  # load metadata
  metadata <- metadata[library_type == library_se, ]
  
  if (!is.null(sample)) {
    sample <- unique(metadata$sample)
  }

  files_r1 <- if (!is.null(metadata$read1_path)) metadata$read1_path else stop("\nERROR: missing R1 FASTQ file. ")
    
  # extract names and sequences of reads assigned to HXB2 genomes
  reads_se <- lapply(files_r1, .seq_info) # here we generate a list of tables
  names(reads_se) <- files_r1 #basename(files_r1) 
  
  # remove empty fastq files from the list
  reads_se <- reads_se[!unlist(lapply(reads_se, is.null))]
  
  # no hiv+ reads found
  if(length(reads_se) == 0){
    cat("\nWARNING!: any read has been reported as infected... ")
    return(NULL)
  } 
  
  # extract hiv+ identified reads
  reads_table <- data.table::rbindlist(
    lapply(names(reads_se), function(file_name) 
      # extract fastq content
      sapply(reads_se[file_name], function(record){
        data.frame(
          sample_file = file_name,
          read_name = record["read_name"],
          seq = record["seq"],
          quality = record["quality"])}
        )))
  
  data.table::setnames(reads_table, names(reads_table), c("sample_file", "read_name", "seq", "quality"))
  
  reads_table <- merge(reads_table, metadata[,.(library_type, sample, read1_path, index_path, cellranger_output, whitelist)], by.x = "sample_file", by.y = "read1_path", all.x = TRUE)
  
  return(reads_table)
}
