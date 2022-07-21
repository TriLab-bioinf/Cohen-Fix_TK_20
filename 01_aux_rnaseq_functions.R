# Auxiliary functions for RNAseq analysis

# Replace ensembl_ids by genbank gene_ids.
replace_gene_acc_by_gb_ids <- function(ensemble_ids, return_all = TRUE){
  entrezids <- mapIds(org.Sc.sgd.db, keys = ensemble_ids, column = c('ENTREZID'), keytype = 'ENSEMBL', multiVals = "first")
  entrezids <- entrezids[!is.na(entrezids)]
  to_name <- ensemble_ids %in% names(entrezids)
  ensemble_ids[to_name] <- as.vector(entrezids)
  if (return_all){
    return(ensemble_ids)
  }
  else {
    return(ensemble_ids[to_name])
  }
}


get_gene_names_from_gene_ids <- function(ensemble_ids, annotation_db, 
                                         look_for = 'ENSEMBL', 
                                         fetch = 'GENENAME'){
  # Reference organism: Saccharomyces cerevisiae => DATABASE = org.Sc.sgd.db
  symbols <- mapIds(annotation_db, keys = ensemble_ids, column = fetch, 
                    keytype = look_for)
  symbols <- symbols[!is.na(symbols)]
  to_name <- ensemble_ids %in% names(symbols)
  ensemble_ids[to_name] <- as.vector(symbols)
  return(ensemble_ids)
}

generate_volcano_plot <- function(res.tmp, my_file_name, log_scale = FALSE){
  dir.create("./Plots", showWarnings = FALSE)
  rownames(res.tmp) <- res.tmp$gene_names
  vp <- EnhancedVolcano(res.tmp, 
                        lab = rownames(res.tmp), 
                        x = 'log2FoldChange', 
                        y = 'padj',
                        pCutoff = 0.05,
                        FCcutoff = 1,
                        pointSize = 1,
                        colAlpha = 4/5,
                        labSize = 2,  # Controls labels size
                        title = res.tmp@elementMetadata$description[2],
                        titleLabSize = 14,
                        subtitle = '', # add subtitle here
                        subtitleLabSize = 12,
                        legendPosition = 'right',
                        legendLabSize = 12,
                        legendIconSize = 2.0,
                        axisLabSize = 10,
                        drawConnectors = TRUE,
                        widthConnectors = 0.25,
                        arrowheads = FALSE,
                        gridlines.major = FALSE,
                        gridlines.minor = FALSE
                        )
  
  vp <- vp + ggplot2::coord_cartesian(xlim=c(-6, 6)) +
              ggplot2::scale_x_continuous(
                  breaks=seq(-6,6, 1))
  
  ggsave2(filename = my_file_name, plot = vp, path = "./Plots", width = 8, height = 11)
  print(vp)
}