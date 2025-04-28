#' @import Biostrings
#' @import data.table
#' @import stringr
#' @import parallel
#' @importFrom ShortRead readFastq
#' @importFrom tidyr %>%

# check if a file path has valid FASTQ extension
# @param path A character vector of file paths.
# @return Logical vector indicating if each path is valid.

.is_valid_fastq <- function(path) {
  grepl("\\.fastq$|\\.fastq\\.gz$", path, ignore.case = TRUE)
}

# recover read's name, nucleotide sequence and quality of hiv+ reads mapped with Kraken2 to HXB2 genome.
# notice that seq names are shared between read sequences and indexes (R1, R2 and R3)
# @param path A character vector of fastq file.
# @return Dataframe indicating each read name, sequence and quality.

.seq_info <- function(fastq){
  if (!file.size(fastq) == 0L){
    sr <- ShortRead::readFastq(fastq)
    rseq <- sr@sread %>% as.data.frame()
    rname <- data.table::tstrsplit(sr@id, " ")[[1]]
    names(rseq) <- 'seq'
    rseq$read_name <- rname
    rquality <- as(Biostrings::quality(sr), "PhredQuality") %>% as.data.frame()
    rseq$quality <- rquality$x
    return(rseq)
  }
}

.setNumberThreads <- function(num_threads){
  if (isFALSE(num_threads)){
    num_threads <- 1
  } else if (isTRUE(num_threads)) {
    num_threads <- parallel::detectCores()
  } else {
    # detect number of available threads or select a number if provided
    # if the number of provided threads is higher than the number of available threads, select the minimum
    avail_threads <- parallel::detectCores()
    num_threads <- min(num_threads, avail_threads)
  }
  return(num_threads)
}

.atac_barcodes <- function(read_name, index_file){ 
  # Find target read within the index file
  # -A3 shows 3 lines after the match
  # -m1 shows only the firs match
  grep_command <- ifelse(grepl("\\.fastq\\.gz$", index_file), "zgrep", "grep")
  
  output <- tryCatch({
    suppressWarnings(system2(grep_command, args = c("-A3", "-m1", read_name, index_file), stdout = TRUE, stderr = FALSE))
  }, warning = function(w) NULL)  
  # Extrae las líneas correctas del formato FASTQ
  
  # Extrae las líneas correctas del formato FASTQ
  barcode_atac <- if (length(output) >= 2) output[2] else NA # second line: nucleotide sequences (barcode)
  quality_seq <- if (length(output) >= 4) output[4] else NA # forth line: sequencing quality
  
  # Extract atac barcode
  # atac barcode is located at the last 16 characters
  index_length <- nchar(barcode_atac)
  
  barcode_atac <- substr(barcode_atac, index_length-15, index_length) 
  barcode_atac_revcomp <- Biostrings::DNAString(barcode_atac) %>% 
    Biostrings::reverseComplement() %>% 
    as.character() # Reverse complement
  
  # Extract atac barcode's quality
  quality_seq <- substr(quality_seq, index_length-15, index_length)
  quality_seq_rev <- paste0(rev(strsplit(quality_seq, NULL)[[1]]), collapse = "") # reverse the quality
  
  return(list(barcode_atac = barcode_atac_revcomp, barcode_atac_quality = quality_seq_rev))
}

# uses bash command to grep a particular read index using the read name
# @param string with a read name and path to index file
# @return String character with the barcode + UMI sequences.

.gex_barcodes <- function(read_name, index_file){ 

  # Find target read within the index file
  # -A3 shows 3 lines after the match
  # -m1 shows only the firs match
  grep_command <- ifelse(grepl("\\.fastq\\.gz$", index_file), "zgrep", "grep")
  
  output <- tryCatch({
    suppressWarnings(system2(grep_command, args = c("-A3", "-m1", read_name, index_file), stdout = TRUE, stderr = FALSE))
  }, warning = function(w) NULL)  
  # Extrae las líneas correctas del formato FASTQ
  barcode_umi <- if (length(output) >= 2) output[2] else NA # second line: nucleotide sequences (barcode + UMI)
  quality_seq <- if (length(output) >= 4) output[4] else NA # forth line: sequencing quality
  
  return(list(barcode_umi = barcode_umi, barcode_umi_quality = quality_seq))
}

# checks if input object is a Seurat Object, if TRUE, extracts meta.data
# @params obj seurat object or seurat's meta.data
# @return data.table with meta.data information

.export_seurat_metadata <- function(obj){
  # checks if obj is a Seurat object, if TRUE, export Seurat's metadata
  if (inherits(obj, "Seurat")) {
    obj_metadata <- obj@meta.data %>% data.table::data.table()
    obj_metadata[, barcode := row.names(obj@meta.data)]
  } else {
    cat("\nWARNING: make sure you have a 'barcode' column with individual cell's barcodes (e.g.: AAACCCAAGACATGCG-1_1)\n")
    obj_metadata <- obj %>% data.table::data.table()
    
  }
  return(obj_metadata)
}

# checks if input object is a Seurat Object, if TRUE, adds new meta.data
# @params obj seurat object or seurat's meta.data
# @params obj_metadata new meta.data to add if Seurat object is provided
# @return Seurat object or data.table with meta.data information

.import_seurat_metadata <- function(obj, obj_metadata){
  # checks if obj is a Seurat object, if TRUE, export Seurat's metadata
  if (inherits(obj, "Seurat")) {
    df_obj_metadata <- obj_metadata %>% as.data.frame()
    obj@meta.data <- df_obj_metadata
    row.names(obj@meta.data) <- df_obj_metadata$barcode
    return(obj)
  }
  return(obj_metadata)
}

# computes summary table of GEX transcripts count or ATAC features count per sample
# @param res A table with recoverBarcode() output from one or several samples
# @return A table with sample, barcode and N columns

summary_counts <- function(res, library_type = NULL){
  # process gex output
  if (library_type == "Gene Expression"){
    gex_umi_count <- res[,.N, by = .(sample, barcode, umi)]
    gex_umi_count$N <- NULL # remove umi count
    summary_count <- gex_umi_count[,.N, by = .(sample, barcode)]
    
  # process atac output
  } else if (library_type == "Chromatin Accessibility"){
    atac_read_count <- res[,.N, by = .(sample, barcode, read_pair)]
    if ( nrow(atac_read_count[read_pair == "R1"]) != nrow(atac_read_count[read_pair == "R2"]) ){
      stop("ERROR: discorcordant paired-end reads assignment (both reads of the pair must be assigned")
    }
    summary_count <- res[,.N, by = .(sample, barcode)]
    summary_count$N <- summary_count$N/2 # report paired-end reads
   
  # no library indicated  
  } else {
    stop("ERROR: missing library_type")
  }
  
  return(summary_count)
}

# computes summary table of infected cells per sample
# @param res A table with recoverBarcode() output from one or several samples
# @return A table with sample, barcode and N columns

summary_cell_counts <- function(res, library_type = NULL){
  
  # process gex and atac output
  if (!is.null(library_type)){
    
    res_counts <- summary_counts(res, library_type)
    #res_counts$N <- NULL # remove count of dna/rna per cell
    summary_count <- res_counts[,.N, by = .(sample)]

  # no library indicated  
  } else {
    stop("ERROR: missing library_type")
  }
  
  return(summary_count)
}

