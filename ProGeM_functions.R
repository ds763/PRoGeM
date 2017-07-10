

GR_object_creator <- function(CHR, START, END, NAMES) {
  ranges_object <- GRanges(seqnames = CHR,
                           ranges = IRanges(start = START, end = END))
  names(ranges_object) <- NAMES
  return(ranges_object)
}


LD_region_range_finder <- function(sentinel_table, proxy_table, overhang) {
  unique_sentinels <- unique(sentinel_data$rsID)
  LD_ranges <- matrix(nrow = length(unique_sentinels), ncol = 3, 
                      dimnames = list(unique_sentinels, c("CHR", "START", "END")))
  for(i in 1:length(unique_sentinels)) {
    LD_ranges[i, 1] <- sentinel_table$CHR[which(sentinel_table$rsID == unique_sentinels[i])]
    LD_ranges[i, 2] <- min(proxy_table$PROXY_START[which(proxy_table$LEAD_rsID == unique_sentinels[i])],
                           sentinel_table$START[which(sentinel_table$rsID == unique_sentinels[i])])
    LD_ranges[i, 3] <- max(proxy_table$PROXY_END[which(proxy_table$LEAD_rsID == unique_sentinels[i])],
                           sentinel_table$END[which(sentinel_table$rsID == unique_sentinels[i])])
  }
  GR_object_creator(LD_ranges[,1], as.numeric(LD_ranges[,2]), as.numeric(LD_ranges[,3]), 
                    rownames(LD_ranges))
}


find_overlapping_genes <- function(ranges, gene_model, overhang) {
  hits <- findOverlaps(ranges, gene_model, maxgap = overhang * 1000) # convert kb to bp
  overlapping_genes <- gene_model[subjectHits(hits)]
  mcols(overlapping_genes) <- names(sentinel_ranges[queryHits(hits)])
  names(mcols(overlapping_genes)) <- "LEAD_rsID"
  return(overlapping_genes)
}


gene_annotator <- function(ids) {
  ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", 
                    dataset="hsapiens_gene_ensembl")
  gene_info <- data.table(getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype", 
                                               "chromosome_name"), filters = "ensembl_gene_id", 
                                values = unique(ids), mart = ensembl))
  # format the gene_info object slightly
  if(length(grep("CHR", gene_info$chromosome_name))) {
    gene_info <- gene_info[-grep("CHR", gene_info$chromosome_name),]
  }
  gene_info[is.na(gene_info)] <- "-"
  gene_info[gene_info == ""] <- "-"
  gene_info <- gene_info[,c(1:3)]
  return(gene_info)
}


gene_variant_distance_finder <- function(sentinel_variant_ranges, local_gene_ranges) {
  distance <- vector(mode = "numeric", length = length(local_gene_ranges))
  rank <- vector(mode = "numeric", length = length(local_gene_ranges))
  for(i in 1:length(sentinel_variant_ranges)) {
    indices <- which(mcols(local_gene_ranges)[[1]] == names(sentinel_variant_ranges)[i])
    if(length(indices) > 0) {
      distance[indices] <- distance(sentinel_variant_ranges[i], 
                                    local_genes[which(mcols(local_gene_ranges)[[1]] == names(sentinel_variant_ranges)[i])])
      temp <- vector(mode = "numeric", length = length(indices))
      temp[order(distance[indices])] <- 1:length(indices)
      rank[indices] <- temp
    }
  }
  mcols(local_gene_ranges)[2] <- distance
  names(mcols(local_gene_ranges))[2] <- "distance_to_sentinel"
  mcols(local_gene_ranges)[3] <- rank
  names(mcols(local_gene_ranges))[3] <- "distance_ranking"
  return(local_gene_ranges)
}


nearest_genes_selector <- function(local_gene_ranges, annotation, biotype, number) {
  nearest_genes <- NULL
  annotation <- annotation[which(annotation$gene_biotype %in% biotype),]
  unique_sentinels <- unique(mcols(local_gene_ranges)[[1]])
  for(i in 1:length(unique_sentinels)) {
    temp <- local_gene_ranges[which(mcols(local_gene_ranges)[[1]] == unique_sentinels[i])]
    temp <- temp[which(names(temp) %in% annotation$ensembl_gene_id)]
    temp_length <- length(temp)
    if(temp_length > 0) {
      for(x in 1:min(number, temp_length)) {
        index <- which(mcols(temp)[[3]] == min(mcols(temp)[[3]]))
        nearest_genes <- rbind(nearest_genes, c(as.character(mcols(temp[index])[[1]]), 
                                                as.character(mcols(temp[index])[[2]]), x, 
                                                as.character(annotation[which(annotation$ensembl_gene_id == names(temp)[index]),])))
        temp <- temp[-index]
      }
    }
  }
  colnames(nearest_genes) <- c("LEAD_rsID", "distance_to_sentinel", "distance_ranking", "ensembl_id",
                               "hgnc_symbol", "gene_biotype")
  return(data.table(nearest_genes))
}



cis_eQTL_target_finder <- function(dir, tissues, sentinels, proxies) {

  sentinel_eQTL_hits <- NULL
  proxy_eQTL_hits <- NULL
  for(i in tissues) {
    current_tissue <- unlist(strsplit(i[1], "_Analysis.nominal.filtered.txt"))
    temp_tissue_data <- data.table(read.table(file = file.path(dir, i), header = TRUE, quote = NULL,
                                              sep = "\t", stringsAsFactors = FALSE))
    temp_search_terms <- gsub("^([^_]*_[^_]*_).*$", "\\1", temp_tissue_data$variant_id)
    
    sentinel_indices <- which(temp_search_terms %in% sentinels$GTEx_search_term)
    tissue <- rep(current_tissue, length(sentinel_indices))
    search_term <- temp_search_terms[sentinel_indices]
    
    sentinel_eQTL_hits <- rbind(sentinel_eQTL_hits, 
                                cbind(search_term, 
                                      temp_tissue_data[sentinel_indices, c(1:4, 9)], 
                                      tissue))
    
    proxy_indices <- which(temp_search_terms %in% proxies$GTEx_search_term)
    tissue <- rep(current_tissue, length(proxy_indices))
    search_term <- temp_search_terms[proxy_indices]
    
    proxy_eQTL_hits <- rbind(proxy_eQTL_hits,
                             cbind(search_term,
                                   temp_tissue_data[proxy_indices, c(1:4, 9)], 
                                   tissue))
    
    names(sentinel_eQTL_hits) <- c("search_term", names(temp_tissue_data)[c(1:4, 9)], "tissue")
    names(proxy_eQTL_hits) <- c("search_term", names(temp_tissue_data)[c(1:4, 9)], "tissue")
    
    cat(paste("\t\tSearch of ", current_tissue, " complete!\n", sep = ""))
  }
  return(list(sentinel_eQTL_hits, proxy_eQTL_hits))
}
   

cis_eQTL_formatter <- function(sentinels, proxies, cis_eQTLs) {
  cis_eQTLs[[1]]$gene_id <- gsub("[.].*$", "", cis_eQTLs[[1]]$gene_id)
  cis_eQTLs[[2]]$gene_id <- gsub("[.].*$", "", cis_eQTLs[[2]]$gene_id)
  
  LEAD_rsID <- vector(mode = "character", length = length(cis_eQTLs[[1]]$search_term))
  for(i in 1:length(sentinels$rsID)) {
    indices <- which(cis_eQTLs[[1]]$search_term == sentinels$GTEx_search_term[i])
    LEAD_rsID[indices] <- sentinels$rsID[i]
  }
  
  cis_eQTLs[[1]] <- cbind(LEAD_rsID, cis_eQTLs[[1]])
  
  PROXY_rsID <- vector(mode = "character", length = length(cis_eQTLs[[2]]$search_term))
  for(x in 1:length(proxies$PROXY_rsID)) {
    indices <- which(cis_eQTLs[[2]]$search_term == proxies$GTEx_search_term[x])
    PROXY_rsID[indices] <- proxies$PROXY_rsID[x]
  }
  
  cis_eQTLs[[2]] <- cbind(PROXY_rsID, cis_eQTLs[[2]])

  return(cis_eQTLs)
}



cis_eQTL_gene_extractor <- function(sentinels, proxies, cis_eQTLs) {
  cis_eQTL_targets <- NULL
  for(i in 1:length(sentinels$rsID)) {
    
    sentinel_genes <- NULL
    if(sentinels$rsID[i] %in% cis_eQTLs[[1]]$LEAD_rsID) {
      sentinel_indices <- which(cis_eQTLs[[1]]$LEAD_rsID == sentinels$rsID[i])
      sentinel_genes <- unique(cis_eQTLs[[1]]$gene_id[sentinel_indices])
    }
    proxy_variants <- proxies$PROXY_rsID[which(proxies$LEAD_rsID == sentinels$rsID[i])]
    
    proxy_genes <- unique(cis_eQTLs[[2]]$gene_id[cis_eQTLs[[2]]$PROXY_rsID %in% proxy_variants])
    if(length(sentinel_genes) > 0 | length(proxy_genes) > 0) {
      all_genes <- unique(unlist(c(sentinel_genes, proxy_genes)))
      sentinel_proxy <- vector(mode = "character", length = length(all_genes))
      for(x in 1:length(all_genes)) {
        if(all_genes[x] %in% sentinel_genes & all_genes[x] %in% proxy_genes) {
          sentinel_proxy[x] <- "lead_and_proxy"
        }
        else if(all_genes[x] %in% proxy_genes) {
          sentinel_proxy[x] <- "proxy"
        }
        else { sentinel_proxy[x] <- "lead" }
      }
      cis_eQTL_targets <- rbind(cis_eQTL_targets, 
                                cbind(rep(sentinels$rsID[i], length(all_genes)), 
                                      all_genes, sentinel_proxy))
      colnames(cis_eQTL_targets) <- c("rsID", "ensembl_id", "sentinel_or_proxy")
    }
  }
  return(data.table(cis_eQTL_targets))
}


top_down_candidate_identifier <- function(SNP_ranges, metabolic_gene_ranges, overhang) {
  metabolic_overlapping_genes <- vector(mode = "list", length = length(metabolic_gene_ranges))
  names(metabolic_overlapping_genes) <- names(metabolic_gene_ranges)
  for(i in 1:length(metabolic_gene_ranges)) {
    hits <- findOverlaps(SNP_ranges, metabolic_gene_ranges[[i]], maxgap = overhang * 1000)
    overlapping_genes <- metabolic_gene_ranges[[i]][subjectHits(hits)]
    sentinels <- names(sentinel_ranges[queryHits(hits)])
    mcols(overlapping_genes)[length(mcols(overlapping_genes)) + 1] <- sentinels
    names(mcols(overlapping_genes))[length(mcols(overlapping_genes))] <- "lead_variant"
    metabolic_overlapping_genes[[i]] <- overlapping_genes
  }
  return(metabolic_overlapping_genes)
}


top_down_scorer <- function(sentinels, metabolic_overlapping_genes) {
  top_down_scores <- NULL
  for(i in 1:length(sentinels$rsID)) {
    lead <- sentinels$rsID[i]
    temp <- vector(mode = "list", length = length(metabolic_overlapping_genes))
    names(temp) <- names(metabolic_overlapping_genes)
    genes <- NULL
    for(x in 1:length(temp)) {
      temp[[x]] <- as.data.frame(mcols(metabolic_overlapping_genes[[x]][which(mcols(metabolic_overlapping_genes[[x]])$lead_variant == lead)]))
      genes <- unique(rbind(genes, temp[[x]][,c(1:5)]))
    }
    if(length(genes[,1]) > 0) {
      top_down_tally <- matrix(data = 0, 
                               nrow = length(genes[,1]), ncol = 6, dimnames = list(NULL, c("GO", "KEGG", "MGI",
                                                                                           "orphanet", "reactome",
                                                                                           "score")))
      for(y in 1:length(genes[,1])) {
        for(z in 1:length(temp)) {
          if(genes[y, 3] %in% temp[[z]][,3]) {
            top_down_tally[y, z] <- 1
          }
        }
        top_down_tally[y, 6] <- sum(top_down_tally[y,])
      }
      LEAD_rsID <- rep(lead, length(genes[,1]))
      top_down_scores <- rbind(top_down_scores, cbind(LEAD_rsID, genes[,c(3:5)], top_down_tally))
    }
  }
  colnames(top_down_scores)[c(2, 3)] <- c("ensembl_id", "hgnc_symbol")
  for(i in 1:length(sentinels$rsID)) {
    indices <- which(top_down_scores[,1] == sentinels$rsID[i])
    ordered <- order(top_down_scores[indices, 10], decreasing = TRUE)
    top_down_scores[indices,] <- top_down_scores[indices[ordered], ]
  }
  
  
  return(data.table(top_down_scores))
}


bottom_up_summariser <- function(sentinels, nearest, LD_overlapping, LD_annotation, cis_eQTLs, 
                                 cis_eQTL_annotation) {
  bottom_up_summary <- NULL
  for(i in 1:length(sentinels$rsID)) {
    lead <- sentinels$rsID[i]
    temp_nearest <- nearest$ensembl_id[which(nearest$LEAD_rsID == sentinels$rsID[i])]
    temp_LD_overlapping <- names(LD_overlapping[which(mcols(LD_overlapping)[[1]] == sentinels$rsID[i])])
    temp_cis_eQTLs <- cis_eQTLs$ensembl_id[which(cis_eQTLs$rsID == sentinels$rsID[i])]
    lead_or_proxy <- cis_eQTLs$sentinel_or_proxy[which(cis_eQTLs$rsID == sentinels$rsID[i])]
    
    all_genes <- unique(c(temp_nearest, temp_LD_overlapping, temp_cis_eQTLs))
    hgnc <- vector(mode = "character", length = length(all_genes))
    biotype <- vector(mode = "character", length = length(all_genes))
    for(z in 1:length(all_genes)) {
      if(all_genes[z] %in% nearest$ensembl_id) {
        hgnc[z] <- nearest$hgnc_symbol[match(all_genes[z], nearest$ensembl_id)]
        biotype[z] <- nearest$gene_biotype[match(all_genes[z], nearest$ensembl_id)]
      }
      else if(all_genes[z] %in% LD_annotation$ensembl_gene_id) {
        hgnc[z] <- LD_annotation$hgnc_symbol[match(all_genes[z], LD_annotation$ensembl_gene_id)]
        biotype[z] <- LD_annotation$gene_biotype[match(all_genes[z], LD_annotation$ensembl_gene_id)]
      }
      else if(all_genes[z] %in% cis_eQTL_annotation$ensembl_gene_id) {
        hgnc[z] <- cis_eQTL_annotation$hgnc_symbol[match(all_genes[z], cis_eQTL_annotation$ensembl_gene_id)]
        biotype[z] <- cis_eQTL_annotation$gene_biotype[match(all_genes[z], cis_eQTL_annotation$ensembl_gene_id)]
      }
    }
    temp_matrix <- matrix(data = 0, 
                          nrow = length(all_genes), ncol = 10, 
                          dimnames = list(NULL, c("LEAD_rsID", "ensembl_id",
                                                  "hgnc_symbol", "gene_biotype", "nearest", 
                                                  "second_nearest", "third_nearest", 
                                                  "LD_overlapping", "lead_eQTL", "proxy_eQTL")))
    for(x in 1:length(all_genes)) {
      if(all_genes[x] %in% temp_nearest) {
        index <- which(temp_nearest == all_genes[x])
        if(index == 1) {
          temp_matrix[x, 5] <- 1
        }
        else if(index == 2) {
          temp_matrix[x, 6] <- 1
        }
        else if(index == 3) {
          temp_matrix[x, 7] <- 1
        }
      }
      if(all_genes[x] %in% temp_LD_overlapping) {
        temp_matrix[x, 8] <- 1
      }
      if(all_genes[x] %in% temp_cis_eQTLs) {
        index <- which(temp_cis_eQTLs == all_genes[x])
        
        
        
        if(lead_or_proxy[index] == "lead") {
          temp_matrix[x, 9] <- 1          
        }
        else if(lead_or_proxy[index] == "proxy") {
          temp_matrix[x, 10] <- 1
        }
        else { temp_matrix[x, c(9, 10)] <- 1 }
      }
    }
    temp_matrix[,1] <- rep(lead, length(all_genes))
    temp_matrix[,2] <- all_genes
    temp_matrix[,3] <- hgnc
    temp_matrix[,4] <- biotype
    bottom_up_summary <- rbind(bottom_up_summary, temp_matrix)
  }
  return(data.table(bottom_up_summary))
}


VEP_integrator <- function(sentinels, bottom_up, VEP, sentinel_col, proxy_col, ensembl_col, IMPACT_col) {
  VEP_matrix <- matrix(data = 0, nrow = length(bottom_up$LEAD_rsID), ncol = 2, 
                       dimnames = list(NULL, c("lead_IMPACT", "proxy_IMPACT")))
  for(i in 1:length(sentinels$rsID)) {
    lead <- sentinels$rsID[i]
    if(lead %in% VEP[[sentinel_col]]) {
      temp <- VEP[which(VEP[[sentinel_col]] == lead),]
      for(x in unique(temp[[ensembl_col]])) {
        temp1 <- temp[which(temp[[ensembl_col]] == x),]
        if("MODERATE" %in% temp1[[IMPACT_col]]) {
          temp2 <- temp1[which(temp1[[IMPACT_col]] == "MODERATE"),]
          boolean <- temp2[[sentinel_col]] == temp2[[proxy_col]]
          if(TRUE %in% boolean) {
            VEP_matrix[which(bottom_up$LEAD_rsID == lead & bottom_up$ensembl_id == x), 1] <- 1
          }
          if(FALSE %in% boolean) {
            VEP_matrix[which(bottom_up$LEAD_rsID == lead & bottom_up$ensembl_id == x), 2] <- 1
          }
        }
        if("HIGH" %in% temp1[[IMPACT_col]]) {
          temp3 <- temp1[which(temp1[[IMPACT_col]] == "HIGH"),]
          boolean <- temp3[[sentinel_col]] == temp3[[proxy_col]]
          if(TRUE %in% boolean) { 
              VEP_matrix[which(bottom_up$LEAD_rsID == lead & bottom_up$ensembl_id == x), 1] <- 2
            }
          if(FALSE %in% boolean) {
            VEP_matrix[which(bottom_up$LEAD_rsID == lead & bottom_up$ensembl_id == x), 2] <- 2
            }
          }
        }
      }
    }
  return(data.table(VEP_matrix))
}


bottom_up_classifier <- function(bottom_up) {
  category <- NULL
  for(i in 1:length(unique(bottom_up$LEAD_rsID))) {
    indices <- which(bottom_up$LEAD_rsID == unique(bottom_up$LEAD_rsID)[i])
    temp <- vector(mode = "integer", length = length(indices))
    for(x in 1:length(temp)) {
      if(sum(bottom_up[indices[x], c(11, 12)]) > 0) {
        temp[x] <- 1
      }
      else if(sum(as.numeric(bottom_up[indices[x], c(5, 6, 7, 8)])) > 0) {
        temp[x] <- 2
      }
      else if(sum(as.numeric(bottom_up[indices[x], c(9, 10)])) > 0) {
        temp[x] <- 3
      }
    }
    ordered <- order(temp)
    bottom_up[indices,] <- bottom_up[indices[ordered],]
    category <- c(category, temp[ordered])
  }
  return(cbind(bottom_up, category))
}


converging_candidate_finder <- function(sentinels, bottom_up, top_down) {
  cooccurring_candidates <- NULL
  for(i in 1:length(sentinels$rsID)) {
    temp_bottom <- bottom_up[which(as.character(bottom_up$LEAD_rsID) == sentinels$rsID[i]),]
    temp_top <- top_down[which(as.character(top_down$LEAD_rsID) == sentinels$rsID[i]),]
    cooccurring_indices <- which(as.character(temp_bottom$ensembl_id) %in% as.character(temp_top$ensembl_id))
    cooccurring_candidates <- rbind(cooccurring_candidates, temp_bottom[cooccurring_indices, c(1:4)])
  }
  colnames(cooccurring_candidates)[c(1, 2)] <- c("LEAD_rsID", "ensembl_id")
  return(data.table(cooccurring_candidates))
}







