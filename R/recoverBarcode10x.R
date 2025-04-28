#' Title: Recover cells barcodes from indexes
#'
#' This function takes a a table with hiv+ reads and its fastq index files to search for the 10x barcodes of the reads
#'
#' @param reads_table A table with hiv reads information
#' @param metadata A table with hiv+ reads metadata.
#' @param nthreads A boolean (T/F) or an integer with the number of threads to use
#' @return A table containing the input table along with the barcode information
#' @examples
#' recoverBarcode10x()
#' @importFrom data.table fread
#' @export

recoverBarcode10x <- function(reads_table, metadata = NULL, cellrang_output_path = NULL, nthreads = FALSE){
  
  # recover library_type from reads_table
  library_sc <- unique(reads_table$library_type)
  
  # read white list file
  if (nrow(unique(reads_table[library_type == library_sc, .(whitelist)])) != 1) {
    stop("\nERROR: More than one white lists paths provided within the meta data file")
  }
  
  # set up number of threads
  nthreads <- .setNumberThreads(nthreads)
  cat("\nRemember that the searching of barcodes can be speeded up by increasing the number of cores. :)\n")
  
  if (library_sc == "Chromatin Accessibility"){

    # recovering ATAC barcodes
    barcode_results <- mclapply(1:nrow(reads_table), function(x) {
      .atac_barcodes(reads_table$read_name[x], reads_table$index_path[x])
    }, mc.cores = nthreads)
    
    barcode_dt <- rbindlist(barcode_results)
    reads_table[, c("barcode_atac", "barcode_atac_quality") := barcode_dt]
    #reads_table[, c("barcode_atac", "barcode_atac_quality") := .atac_barcodes(read_name, index_path, nthreads), by = read_name]
    
    # convert atac to gex barcodes
    whitelist_path <- unique(reads_table[library_type == library_sc, whitelist])
    whitelist_atac_barcodes <- readLines(whitelist_path,
                                         warn = FALSE)

    # recover GEX barcodes (used to associate reads to cells)
    whitelist_gex_path <- unique(metadata[library_type == "Gene Expression", whitelist])
    whitelist_gex_path <- if (length(whitelist_gex_path) == 1) whitelist_gex_path else gsub("/atac/barcodes/", "/cellranger/barcodes/", whitelist_path)    
    barcodes_atac2gex <- readLines(whitelist_gex_path,
                                   warn = FALSE)    
    names(barcodes_atac2gex) <- whitelist_atac_barcodes
    remove(whitelist_atac_barcodes)
    
    # add gex barcode as "barcode" column
    reads_table[, barcode := barcodes_atac2gex[barcode_atac]]
    reads_table[is.na(barcode), barcode := "WRONG BARCODE"]
    
  } else if (library_sc == "Gene Expression"){

    # look for the barcode within the index file and append the barcode + UMI sequence in the table
    barcode_results <- mclapply(1:nrow(reads_table), function(x) {
      .gex_barcodes(reads_table$read_name[x], reads_table$index_path[x])
    }, mc.cores = nthreads)
    
    barcode_dt <- rbindlist(barcode_results)
    reads_table[, c("barcode_umi", "barcode_umi_quality") := barcode_dt]

    # extract barcodes' sequences and qualities (notice that barcode is located at the 16 starting nucleotides)
    reads_table[, barcode := substr(barcode_umi, 1, 16), by = read_name]
    reads_table[, barcode_quality := substr(barcode_umi_quality, 1, 16), by = read_name]
    
    # extract UMI (notice that UMI is located at the last 12 nucleotides)
    reads_table[, umi := substr(barcode_umi, 17, nchar(barcode_umi)), by = read_name]
    
    reads_table$barcode_umi <- NULL
    reads_table$barcode_umi_quality <- NULL
    gc()

  } else {
    stop('\nERROR: Invalid data type, neither Chromatin Accessibility or Gene Expression provided')
  }
  
  return(reads_table)
}
