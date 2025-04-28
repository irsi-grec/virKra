#' Title: Load Metadata
#'
#' This function takes a csv file with metadata information and loads it into the workspace as a data.table.
#' The csv must have the following columns: library_type, sample, read1, read2, index, whitelist.
#'
#' @param x A csv file.
#' @return A data.table with the sample metadata.
#' @examples
#' loadMetadata("mymetadata.csv")
#' @importFrom data.table fread
#' @export

loadMetadata <- function(file){
  metadata <- data.table::fread(file)
  
  # check table content (minimum one row of data)
  if (nrow(metadata) < 1) {
    stop("The metadata file must contain at least one row of data.")
  }
  
  # check required columns
  required_columns <- c("library_type", "sample", "read1", "read2", "index", "whitelist")
  if (!all(required_columns %in% colnames(metadata))) {
    stop("The metadata file must contain the following columns: ",
         paste(required_columns, collapse = ", "))
  }
  
  # set names
  data.table::setnames(metadata, c("read1", "read2", "index"), c("read1_path", "read2_path", "index_path"))
  
  # check accepted values in library_type
  valid_libraries <- c("Gene Expression", "Chromatin Accessibility")
  if (!all(metadata$library_type %in% valid_libraries)) {
    stop('Invalid library_type. Must be "Gene Expression" or "Chromatin Accessibility".')
  }
  
  # check sample column for non-empty values
  if (any(is.na(metadata$sample) | metadata$sample == "" | length(unique(metadata$sample)) != 1)) {
    stop("The sample column contains empty or missing values.")
  } else {
    cat("\nREMEMBER! sample column MUST contain the same value as orig.ident within the Seurat object metadata")
  }
  
  # check read1 column
  if (!all(.is_valid_fastq(metadata$read1_path))) {
    stop("The read1 column must contain paths to files with .fastq or .fastq.gz extensions.")
  }
  
  # check read2 column if present
  if (!all(is.na(metadata$read2_path) | metadata$read2_path == "" | .is_valid_fastq(metadata$read2_path))) {
    stop("The read2 column, if present, must contain paths to files with .fastq or .fastq.gz extensions.")
  }
  
  # check index column
  if (!all(.is_valid_fastq(metadata$index_path))) {
    stop("The index column must contain paths to files with .fastq or .fastq.gz extensions.")
  }
  
  # check whitelist column for 10x
  if (!all(grepl("\\.txt\\.gz$", metadata$whitelist, ignore.case = TRUE))) {
    stop("The whitelist column must contain paths to files with .txt.gz extensions.")
  }
  
  # return metadata if all checks pass
  return(metadata)

}
