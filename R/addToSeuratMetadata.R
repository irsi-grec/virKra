#' Title: Export infected cell information into the Seurat Object's meta data
#'
#' This function takes a a table with hiv+ cell information: barcodes, number of transcripts/dna and sample
#' and uses  barcodes and sample information to include hiv+ cells data into the seurat object meta data.
#'
#' @param x A table with hiv+ cells metadata.
#' @return A table containing the seurat object metadata with infected cells information
#' @examples
#' addToSeuratMetadata()
#' @importFrom data.table fread
#' @export

addToSeuratMetadata <- function(reads_table, meta_data, library_type = NULL){
  # check if seurat_metada is a Seurat object or table
  # returns a data.table with barcode information
  seurat_metadata <- .export_seurat_metadata(meta_data) 
  
  # formatting seurat meta data
  seurat_metadata$barcode_cp <- seurat_metadata$barcode
  seurat_metadata$barcode <- NULL
  seurat_metadata[, barcode := sub("[-_].*", "", barcode_cp)]
  
  if (library_type == "Chromatin Accessibility"){
    # formatting results into cell counts
    count_table <- summary_counts(reads_table, library_type)
    data.table::setnames(count_table, "N", "N_infected_atac")
    
  } else if (library_type == "Gene Expression"){
    # formatting results into cell counts
    count_table <- summary_counts(reads_table, library_type)
    data.table::setnames(count_table, "N", "N_infected_gex")
    
  } else {
    stop('\nERROR: Invalid data type, neither Chromatin Accessibility or Gene Expression provided')
  }
  
  # ensure that "samples" are named as "orig.ident"
  if( !all(unique(count_table$sample) %in% unique(seurat_metadata$orig.ident)) ){
    stop("\nERROR: sample name cannot be fount within the Seurat Object meta data")
  }
  
  # include infection information
  seurat_new_metadata <- merge(seurat_metadata, count_table, by.x = c("orig.ident", "barcode"), by.y = c("sample", "barcode"), all.x = T)
  seurat_new_metadata[is.na(seurat_new_metadata)] <- 0  
  
  # remove auxilary barcode column
  seurat_new_metadata$barcode <- seurat_new_metadata$barcode_cp
  seurat_new_metadata$barcode_cp <- NULL
  
  # return metadata as a table or within Seurat's object meta.data
  out_meta_data <- .import_seurat_metadata(meta_data, seurat_new_metadata)
  return(out_meta_data)
}
