# common libraries
library(glue)
library(RColorBrewer)

#' Compute ssGSEA scores on bulk data given SC-derived signatures
#' 
#' Given path to cluster marker input files,
#' reformat the matrix, and compute the scores with the
#' bulk RNA-seq data
#'
#' @param marker_path character Path to marker incidence matrix
#' @param counts list Matrix of counts from the RNA-seq data
#' @param out_prefix character Prefix for path/file name where scores should be saved as tsv
#' will be appended with ".ssgsea_scores.txt"
computeScores <- function(marker_path, counts, out_prefix, test_mode = FALSE) {
  
  markers.mat <- read_tsv(marker_path)
  
  require(GSVA)
  
  markers <- markers.mat %>%
    .[, 6:ncol(.)] %>%
    sapply(function (i) markers.mat[["hsapiens_ensembl_gene_id"]][i], 
           simplify = FALSE)
  
  if (test_mode) { # Just run for two samples, to make sure it works
    
    scores <- GSVA::gsva(expr          = counts[,c(1, 2)],
                         gset.idx.list = markers,
                         mx.diff       = FALSE,
                         method        = "ssgsea",
                         rnaseq        = TRUE,
                         ssgsea.norm   = FALSE,
                         tau           = 0.75)
    
    return(scores)
    
  }
  
  scores <- GSVA::gsva(expr          = counts,
                       gset.idx.list = markers,
                       mx.diff       = FALSE,
                       method        = "ssgsea",
                       rnaseq        = TRUE,
                       ssgsea.norm   = FALSE,
                       tau           = 0.75)
  
  scores %>%
    data.frame %>% 
    magrittr::set_colnames(colnames(counts)) %>% 
    tibble::rownames_to_column(var = "Signature") %>%
    {write.table(file = glue("{out_prefix}.ssgsea_scores.txt"),
                 x = .,
                 sep = "\t",
                 quote = T,
                 row.names = FALSE)}
  
  invisible(scores)
  
} # this function doesn't seem to be useful in my case <- just to simplify the process

plotHallmark <- function(seurat, prefix, figout = NULL,
                         summary_type = "cluster_mean",
                         cluster_rows = TRUE,
                         cluster_cols = TRUE,
                         drop_sigs = NULL,
                         print = FALSE,
                         hallmark_sigs = NULL,
                         clusters = NULL) {
  
  require(pheatmap)
  
  scores <- suppressMessages(read_tsv(glue("{prefix}.hallmark.{summary_type}.ssgsea_scores.txt")))
  
  if (is.null(seurat@misc$colours)) seurat@misc$colours <- getClusterColours(seurat)
  
  anno_row <- data.frame(cluster = names(seurat@misc$colours))
  rownames(anno_row) <- anno_row$cluster
  side_colors <- list(cluster = seurat@misc$colours)
  
  if (!is.null(drop_sigs)) scores <- scores %>% 
    filter(!(Signature %in% drop_sigs))
  
  # Tidy names
  scores$Signature <- scores$Signature %>%
    lapply(str_split_fixed, "_", n = 2) %>% 
    sapply(getElement, 2)
  
  mat <- scores %>%
    as.data.frame %>%
    tibble::column_to_rownames(var = "Signature") %>%
    t()
  
  if (is.null(hallmark_sigs)) hallmark_sigs <- colnames(mat)
  
  fig <- figout# glue("{prefix}.hallmark")
  
  rdbu <- rev(colorRampPalette(brewer.pal(8, "RdBu"))(n = 100))
  
  if (is.null(clusters)) clusters <- levels(seurat@ident)
  
  hm_fun <- purrr::partial(pheatmap,
                           mat[clusters, hallmark_sigs],
                           scale = "none",
                           color = rdbu,
                           cluster_rows = cluster_rows,
                           cluster_cols = cluster_cols,
                           border_color = NA,
                           cellwidth = 13,
                           cellheight = 13,
                           main = paste0(seurat@project.name, " - Hallmark gene set"),
                           annotation_row = anno_row,
                           annotation_colors = side_colors,
                           annotation_legend = FALSE)
  
  if (print) invisible(hm_fun())
  else {
    
    hm_fun(filename = glue("{fig}.pdf"))
    hm_fun(filename = glue("{fig}.png"))
    knitr::include_graphics(glue("{fig}.png"))
    
  }  
  
}

#' integrate_cdfs
#' 
#' This code was adapted from GSVA and the GSEA code from the Broad Institute
#' 
#' @param gene_set Character vector of genes
#' @param gene_ranking Numeric vector specifying gene ranking
#' @param R Matrix of gene ranks
#' @param j Integer, specifying sample (column of expression matrix)
#' @param alpha Numeric, exponent giving the weight in the weighted ECDF $P^w_{in}$
integrate_cdfs <- function(gene_set, gene_ranking, expr_mat, R, j, alpha, raw_rank = TRUE) {
  
  # Binary 1/0 vector indicating whether each gene (in ranked order)
  # is in the gene set or not
  indicator_in_gene_set <- match(rownames(expr_mat)[gene_ranking], gene_set)
  indicator_in_gene_set[!is.na(indicator_in_gene_set)] <- 1
  indicator_in_gene_set[ is.na(indicator_in_gene_set)] <- 0
  indicator_out <- 1 - indicator_in_gene_set
  
  # P_in (the weighted ECDF of genes in the set)
  # i.e. it will take a weighted step whenever a gene is in the set
  P_in <- cumsum( (abs(R[gene_ranking, j]) * indicator_in_gene_set)^alpha ) /
    sum( (abs(R[gene_ranking, j]) * indicator_in_gene_set)^alpha )
  
  # P_out (the ECDF of genes not in the set)
  # i.e. it will be length N - Ng, with a step whenever the gene is not in the set
  P_out <- cumsum( !indicator_in_gene_set ) /
    sum( !indicator_in_gene_set )
  
  # The difference in CDFs
  cdf_diffs <- P_in - P_out
  es <- sum(cdf_diffs)
  
  # Calculate running sum statistic, and leading edge subset
  n_g <- length(gene_set)
  n <- nrow(R)
  steps_in <- (abs(R[gene_ranking, j]) * indicator_in_gene_set)^alpha
  step_norm_in <- 1/sum(steps_in)
  step_norm_out <- 1/(n - n_g)
  
  running_sum <- cumsum(indicator_in_gene_set * steps_in * step_norm_in -  # Increment the score by a weighted step if a gene is a hit
                          indicator_out * step_norm_out)                   # Decrement by an fixed size step if a gene is a miss
  
  # The leading edge consists of all the genes in the gene set, which come up
  # in the ranked list prior to the max enrichment score (diff. between ECDFs)
  leading_edge <- names(running_sum)[1:which.max(running_sum)]
  
  # We only care about the ones in the gene set
  if (raw_rank) { # whether output raw ranks
    le_raw_ranks <- data.frame(gene = leading_edge, rank = 1:which.max(running_sum)) %>% 
      add_row(gene = rownames(R)[which(!rownames(R) %in% leading_edge)], rank = 0) %>% 
      column_to_rownames("gene")
    leading_edge_indicator <- le_raw_ranks[rownames(R),]
    names(leading_edge_indicator) <- rownames(R)
  } else {
    leading_edge_indicator <- 1 * ((rownames(R) %in% leading_edge) & (rownames(R) %in% gene_set))
    names(leading_edge_indicator) <- rownames(R)
  }
  
  return(list(
    leading = leading_edge_indicator,
    es = sum(cdf_diffs),
    running = running_sum))
  
}

ssgsea_le <- function(expr_mat, gene_sets, alpha = 0.25, normalize = TRUE,
                      save_le = TRUE, raw_rank = TRUE, n_cores = 1) {
  
  require(pbapply)
  
  message("1. Converting gene expr to ranks...")
  
  # Convert gene expression data in X to ranks within each sample
  R <- apply(expr_mat, 2, function(x) as.integer(rank(x, na.last = TRUE)))
  rownames(R) <- rownames(expr_mat)
  
  # # Collect the differences between CDFs
  # diff_list <- vector("list", length = n)
  
  message(glue("2. Calculating enrichment, parallelizing over ", n_cores, " cores..."))
  
  # For each sample S_1 to S_n, calculate the score for each gene set
  # Parallelize over samples
  ssgsea_out <- pblapply(1:ncol(expr_mat), function(j) {
    
    # The order (by indices) by which to retrieve genes from the matrix
    # to get them from high -> low
    
    # BUG: I think that subsetting the dataframe in this way does not properly
    # deal with any NAs for that sample
    gene_ranking <- order(R[, j], decreasing = TRUE, na.last = TRUE)
    sample_out <- lapply(gene_sets, integrate_cdfs, gene_ranking, expr_mat, R, j, alpha, raw_rank)
    
    return(sample_out)
    
  })
  
  # TODO: this is super unwieldy and unneeded!
  # Would be better to simply return the genes in the set which are in the leading
  # edge as a char vector, that could be saved.
  # But, that may not be the solution, since the leading edge of interest,
  # in our case, is defined per-sample, and so we'll only really care about it
  # with respect to biological groups of samples
  filter_zeros <- function(df) {
    
    df[rowSums(df) != 0, ]
    
  }
  
  assign_nas <- function(df) {
    na_if(df, 0)
  }
  
  message("3. Tidying ouput...")
  
  if (save_le) {
    
    # Tidy the binary data indicating if genes are within the leading edge subset
    leading_mat <- lapply(ssgsea_out,
                          function(sample_out) purrr::transpose(sample_out)$leading) %>% 
      set_names(colnames(expr_mat)) %>%
      purrr::transpose() %>%
      lapply(bind_rows) %>% 
      lapply(as.data.frame) %>% 
      lapply(magrittr::set_rownames, rownames(expr_mat)) %>% 
      lapply(filter_zeros)
    
    if (raw_rank) {
      leading_mat <- leading_mat %>% 
        lapply(assign_nas)
      leading_mat <- lapply(1:length(gene_sets), function(i) {
        df <- leading_mat[[names(gene_sets)[i]]]
        return(df[rownames(df) %in% gene_sets[[i]], ])
      })
      names(leading_mat) <- names(gene_sets)
    }
    
    # Tidy the running sum data
    # Can't coerce to a dataframe because then order would be lost
    running_sum <- lapply(ssgsea_out,
                          function(sample_out) purrr::transpose(sample_out)$running) %>% 
      set_names(colnames(expr_mat)) %>%
      purrr::transpose()
    
  } else {
    
    leading_mat <- NULL
    
  }
  
  # Tidy the enrichment scores for each set for each sample
  es <- lapply(ssgsea_out, function(sample_out) unlist(map(sample_out, "es"))) %>%
    data.frame %>%
    as.matrix()
  
  # Normalize the whole dataset
  normES <- function(mat) apply(mat, 2, function(x, mat) x /
                                  (range(mat)[2] - range(mat)[1]), mat)
  
  if (normalize) es <- normES(es)
  
  rownames(es) <- names(gene_sets)
  colnames(es) <- colnames(expr_mat)
  
  return(list("enrichment_scores" = es,
              "leading_edge" = leading_mat))
  
  message("4. Done.")
  
}

#' @param genelists A named list, where each element is a character vector of
#' genes. Expects gene symbols.
#' @param summary_type Character specifying how to treat single cells for 
#' the enrichment. If "cluster_mean", computes cluster mean expression and
#' evaluates enrichment per cluster. If "metacells", assumes the seurat object
#' contains metacells (by \code{seurat2metacells}) and evaluates enrichment
#' per metacell.
#' @param gene_space Character vector with gene symbols, used to restrict
#' the space of genes considered in ssGSEA. Note, this is *not* the list of genes
#' for which enrichment will be assessed. Rather, the enrichment of gene set
#' will be computed among these genes.
compute_scores_sc_le <- function(seurat, genelists,
                                 out_prefix = "./",
                                 summary_type = "cluster_mean",
                                 n_cores = 1,
                                 save_le = FALSE,
                                 gene_space = NULL,
                                 return_counts = FALSE) {
  
  require(pbapply)
  
  # Get counts
  if (summary_type == "cluster_mean") counts <- seurat %>% cytokit::meanClusterExpression()
  else if (summary_type %in% c("metacells", "single_cells")) {
    
    counts <- as.matrix(seurat@data)
    
  }
  
  # Reduce counts to the space of genes to consider
  if (!is.null(gene_space)) {
    gene_space_detected <- unique(gene_space[gene_space %in% rownames(counts)])
    counts <- counts[gene_space_detected, ]
    message("0. Reducing the gene space to the ", nrow(counts), " genes detectd from the provided genes")
  }
  
  if (return_counts) return(counts)
  
  out <- ssgsea_le(expr_mat  = counts,
                   gene_sets = genelists,
                   alpha     = 0.75,
                   normalize = FALSE,
                   save_le   = save_le,
                   n_cores   = n_cores)
  
  # Tidy scores and write to file
  out$enrichment_scores %>%
    data.frame %>% 
    magrittr::set_rownames(NULL) %>% 
    magrittr::set_rownames(rownames(out$enrichment_scores)) %>% 
    magrittr::set_colnames(colnames(counts)) %>% 
    tibble::rownames_to_column(var = "Signature") %>%
    {write.table(file = glue("{out_prefix}.{summary_type}.ssgsea_scores.txt"),
                 x = .,
                 sep = "\t",
                 quote = TRUE,
                 row.names = FALSE)}
  
  if (isTRUE(save_le)) {
    
    leading_edge <- out$leading_edge
    save(leading_edge, file = glue("{out_prefix}.{summary_type}.ssgsea_le.Rda"))  
    
  }
  
  if (length(genelists) == 1) return(unlist(out$enrichment_scores)[1, ])
  
}

# From: /mnt/KLEINMAN_JBOD1/SCRATCH/projects/sjessa/from_hydra/side_projects/2019/2019-10_ssgsea/ssgsea.R
#' Rank the leading edge genes
#' 
#' For a given set of leading edge genes, rank them by their median rank expression
#' across a given sample group.
#' 
#' @param leading_edge_only Logical, if TRUE, compute the ranks only for genes in 
#' the leading edge subset. Otherwise, compute ranks for all genes in the signature.
#' Default: TRUE
#' @param return_symbols Logical, if TRUE, return the gene symbols, otherwise return
#' the Ensembl IDs. Default: TRUE
#' @param relative_rank Logical, if TRUE, return ranks starting from 1, 2, 3... etc, 
#' describing how highly expressed these genes are across the samples, relative to each other.
#' If FALSE, return the actual median ranks of the genes across all samples. Default: FALSE
rank_leading_edge <- function(leading_edge_path, counts, cluster,
                              proportion = 1,
                              samples = NULL,
                              leading_edge_only = TRUE,
                              return_symbols = TRUE,
                              relative_rank = FALSE,
                              return_raw_ranks = FALSE,
                              return_gene_set = FALSE,
                              debug = FALSE) {
  
  expr_mat <- counts
  
  message("1/3: Getting leading edge...")
  if (leading_edge_only) gene_set <- get_leading_edge(leading_edge_path,
                                                      samples = samples,
                                                      cluster = cluster,
                                                      proportion = proportion,
                                                      return_symbols = FALSE)
  else gene_set <- get_leading_edge(leading_edge_path,
                                    samples = samples,
                                    cluster = cluster,
                                    proportion = proportion,
                                    return_symbols = FALSE,
                                    leading_edge_only = FALSE)
  
  
  if (debug) return(gene_set)
  
  message("2/3: Computing ranks...")
  # Convert gene expression data in expr_mat to ranks within each sample
  # By negating the expression values, we will rank the genes from high expression
  # to low, i.e. highly-expressed genes will have small ranks
  R <- apply(expr_mat, 2, function(x) as.integer(rank(-x)))
  rownames(R) <- rownames(expr_mat)
  R_gene_set <- R[gene_set, ]
  
  if (!is.null(samples) & length(samples) == 1) {
    df <- data.frame(R_gene_set[, samples]) %>%
      setNames("rank") %>% 
      tibble::rownames_to_column(var = "gene") %>% 
      arrange(rank) %>% 
      as.data.frame %>% 
      tibble::column_to_rownames(var = "gene")
    return(df)
  }
  else if (!is.null(samples) & length(samples) > 1) R_gene_set <- R_gene_set[, samples]
  
  if (return_gene_set) return(R_gene_set)
  
  message("3/3: Computing medians...")
  # Calculate the median rank of each gene across the relevant samples
  
  message("...Computing median ranks across ", ncol(R_gene_set), " samples")
  
  if (!return_raw_ranks) {
    
    if (return_symbols) ranks <- data.frame(gene = ensembl2symbols_safe(rownames(R_gene_set)),
                                            median_rank = matrixStats::rowMedians(R_gene_set))
    else ranks <- data.frame(gene = rownames(R_gene_set),
                             median_rank = matrixStats::rowMedians(R_gene_set))
    
    # Make the ranks relative, to start from 1, 2, 3, etc.
    if (relative_rank) ranks <- ranks %>% mutate(median_rank = min_rank(median_rank))
    
  } else {
    
    if (return_symbols) {
      
      ranks <- as.data.frame(R_gene_set) %>% 
        tibble::add_column(gene = ensembl2symbols_safe(rownames(R_gene_set)), .before = 1) %>% 
        tibble::add_column(median_rank = matrixStats::rowMedians(R_gene_set), .after = "gene")
      
    } else {
      
      ranks <- as.data.frame(R_gene_set) %>% 
        tibble::add_column(gene = rownames(R_gene_set), .before = 1) %>% 
        tibble::add_column(median_rank = matrixStats::rowMedians(R_gene_set), .after = "gene")
      
    }
    
  }
  
  return(arrange(ranks, median_rank))
  
}

# From the same folder as above
#' Get leading edge genes
#'
#' For a given signature (cluster), return the genes which are
#' in the leading edge subset in some proportion of the samples. Since more, rather
#' than fewer, genes in the gene set tend to be in the leading edge, the default
#' proportion is 1, i.e. the most stringent.
#' 
#' @param leading_edge The leading edge list object produced by \code{compute_scores_le}
get_leading_edge <- function(leading_edge_path,
                             cluster,
                             samples = NULL,
                             proportion = 1,
                             return_mat = FALSE,
                             return_symbols = TRUE,
                             leading_edge_only = TRUE) {
  
  load(leading_edge_path)
  
  if (!(cluster %in% names(leading_edge))) stop("Cluster not found in signatures...")
  
  if (is.null(samples)) le_clust <- leading_edge[[cluster]]
  else if (length(samples) == 1) {
    
    # Simply return the leading edge genes, by definition, for one sample  
    return(rownames(leading_edge[[cluster]])[which(leading_edge[[cluster]][, samples] == 1)])
    
  } else {
    
    # Subset to the leading edge info for all relevant samples
    le_clust <- leading_edge[[cluster]][, samples]
    message("...Identifying leading edge genes for ", ncol(le_clust), " samples")
    
  }
  
  if (return_mat) return(le_clust)
  
  if (leading_edge_only) le_genes <- le_clust[rowSums(le_clust) >= proportion * ncol(le_clust), ] %>%
      rownames()
  else le_genes <- le_clust %>% rownames()
  
  if (return_symbols & grepl("^ENS", head(le_genes[[1]]))) return(ensembl2symbols(le_genes))
  else return(le_genes)
  
  
}


#' Changes from Selin's plotHallmark function:
#' * added long_name option
#' * changed text output name and figure main title 
#' @colrange if NULL, use the default color range, i.e. [min, max]; 
#' otherwise, input a vector of min and max to be used as color range.
plotGSEA_pathway <- function(seurat, prefix, figout = NULL, 
                             out_dir = NULL, 
                             summary_type = "cluster_mean",
                             cluster_rows = TRUE,
                             cluster_cols = TRUE,
                             zscale = FALSE, 
                             drop_sigs = NULL,
                             save_sigs = NULL, 
                             suffix = NULL, 
                             print = FALSE,
                             geneset_sigs = NULL,
                             clusters = NULL, 
                             long_name = TRUE, 
                             input_dir = NULL, 
                             colrange = c(-3, 3)) {
  
  require(pheatmap)
  
  if(is.null(input_dir)) input_dir <- out_dir
  
  scores <- suppressMessages(read_tsv(glue("{input_dir}", "{prefix}.{summary_type}.ssgsea_scores.txt")))
  num_col <- ncol(scores)
    
  if (is.null(seurat@misc$colours)) seurat@misc$colours <- getClusterColours(seurat)
  
  anno_row <- data.frame(cluster = names(seurat@misc$colours))
  rownames(anno_row) <- anno_row$cluster
  side_colors <- list(cluster = seurat@misc$colours)
  
  # Tidy names
  if(!long_name) {
    scores$Signature <- scores$Signature %>%
      lapply(str_split_fixed, "_", n = 2) %>% 
      sapply(getElement, 2)
  }
  
  mat <- scores %>%
    as.data.frame %>%
    tibble::column_to_rownames(var = "Signature") %>%
    t()
  
  if (zscale) {
    
    mat <- zscoring_df(mat, 2)
    
    mat %>% 
      t() %>% 
      data.frame %>% 
      magrittr::set_colnames(colnames(scores)[2:num_col]) %>%
      tibble::rownames_to_column(var = "Signature") %>% 
      {write.table(file = glue("{out_dir}", "{prefix}.{summary_type}.ssgsea_scores_zscale.txt"),
                   x = .,
                   sep = "\t",
                   quote = T,
                   row.names = FALSE)}
    
    figout <- glue(figout, "_zscale")
  }
  
  if (!is.null(drop_sigs)) mat <- mat[, -which(colnames(mat) %in% drop_sigs)]
  
  if (!is.null(save_sigs)) mat <- mat[, save_sigs]
  
  if (is.null(geneset_sigs)) geneset_sigs <- colnames(mat)
  
  if (!is.null(suffix)) figout <- glue(figout, "_", suffix)
  
  if (cluster_cols) {
    fig <- glue(figout, "_clucolum")
  } else {
    fig <- figout
  }
  
  rdbu <- rev(colorRampPalette(brewer.pal(8, "RdBu"))(n = 100))
  if (zscale && !is.null(colrange)) {
    breakList <- seq(colrange[1], colrange[2], length.out = 100)
  } else {
    breakList <- NULL
  }
  
  if (is.null(clusters)) clusters <- levels(seurat@ident)
  
  hm_fun <- purrr::partial(pheatmap,
                           mat[clusters, geneset_sigs],
                           scale = "none",
                           color = rdbu, 
                           breaks = breakList, 
                           cluster_rows = cluster_rows,
                           cluster_cols = cluster_cols,
                           border_color = NA,
                           cellwidth = 13,
                           cellheight = 13,
                           main = paste0(seurat@project.name, " - ", gsub("\\.", "", prefix), " gene set"),
                           annotation_row = anno_row,
                           annotation_colors = side_colors,
                           annotation_legend = FALSE)
  
  if (print) knitr::include_graphics(glue("{out_dir}", "{fig}.png")) # invisible(hm_fun()) #why needs to make it invisible when print?
  else {
    
    hm_fun(filename = glue("{out_dir}", "{fig}.pdf"))
    hm_fun(filename = glue("{out_dir}", "{fig}.png"))
    
  }
  
}

# use zscore table as the input
plotGSEA_pathway1 <- function(seurat, scores, 
                              prefix = NULL, 
                              figout = NULL, 
                              out_dir = NULL, 
                              summary_type = "cluster_mean",
                              cluster_rows = TRUE,
                              cluster_cols = TRUE,
                              zscale = FALSE, 
                              drop_sigs = NULL,
                              save_sigs = NULL, 
                              suffix = NULL, 
                              print = FALSE, 
                              output = TRUE, 
                              geneset_sigs = NULL,
                              clusters = NULL, 
                              long_name = TRUE, 
                              input_dir = NULL, 
                              colrange = c(-3, 3)) {
  
  require(pheatmap)
  
  #if(is.null(input_dir)) input_dir <- out_dir
  
  #if (zscale) {
  #  scores <- suppressMessages(read_tsv(glue("{input_dir}", "{prefix}.{summary_type}.ssgsea_zscores.txt")))
  #  figout <- glue(figout, "_zscale")
  #} else {
  #  scores <- suppressMessages(read_tsv(glue("{input_dir}", "{prefix}.{summary_type}.ssgsea_scores.txt")))
  #}
  
  print("Please make sure the scores are z-score normalized")
  
  num_col <- ncol(scores)
  
  if (is.null(seurat@misc$colours)) seurat@misc$colours <- getClusterColours(seurat)
  
  anno_row <- data.frame(cluster = names(seurat@misc$colours))
  rownames(anno_row) <- anno_row$cluster
  side_colors <- list(cluster = seurat@misc$colours)
  
  # Tidy names
  if(!long_name) {
    scores$Signature <- scores$Signature %>%
      lapply(str_split_fixed, "_", n = 2) %>% 
      sapply(getElement, 2)
  }
  
  mat <- scores %>%
    as.data.frame %>%
    tibble::column_to_rownames(var = "Signature") %>%
    t()
  
  if (!is.null(drop_sigs)) mat <- mat[, -which(colnames(mat) %in% drop_sigs)]
  
  if (!is.null(save_sigs)) mat <- mat[, save_sigs]
  
  if (is.null(geneset_sigs)) geneset_sigs <- colnames(mat)
  
  if (!is.null(suffix)) figout <- glue(figout, "_", suffix)
  
  if (cluster_cols) {
    fig <- glue(figout, "_clucolum")
  } else {
    fig <- figout
  }
  
  rdbu <- rev(colorRampPalette(brewer.pal(8, "RdBu"))(n = 100))
  if (zscale && !is.null(colrange)) {
    breakList <- seq(colrange[1], colrange[2], length.out = 100)
  } else {
    breakList <- NULL
  }
  
  if (is.null(clusters)) clusters <- levels(seurat@ident)
  
  hm_fun <- purrr::partial(pheatmap,
                           mat[clusters, geneset_sigs],
                           scale = "none",
                           color = rdbu, 
                           breaks = breakList, 
                           cluster_rows = cluster_rows,
                           cluster_cols = cluster_cols,
                           border_color = NA,
                           cellwidth = 13,
                           cellheight = 13,
                           main = paste0(seurat@project.name, " - ", gsub("\\.", "", prefix), " gene set"),
                           annotation_row = anno_row,
                           annotation_colors = side_colors,
                           annotation_legend = FALSE)
  
  if (print) knitr::include_graphics(glue("{out_dir}/{fig}.png")) # invisible(hm_fun()) #why needs to make it invisible when print?
  else {
    hm_fun(filename = glue("{out_dir}/{fig}.pdf"))
    hm_fun(filename = glue("{out_dir}/{fig}.png"))
  }
  
}

plotGSEA_pathway_allfigs <- function(seurat, prefix, figout = NULL, 
                                     out_dir = NULL, 
                                     summary_type = "cluster_mean",
                                     cluster_rows = TRUE,
                                     cluster_cols = "both",
                                     zscale = "both", 
                                     drop_sigs = NULL,
                                     save_sigs = NULL, 
                                     suffix = NULL, 
                                     print = FALSE,
                                     geneset_sigs = NULL,
                                     clusters = NULL, 
                                     long_name = TRUE, 
                                     input_dir = NULL, 
                                     colrange = c(-3, 3)) {
  
  # list the situations that needs to draw more than one figures
  
  if (cluster_cols == "both" & zscale == "both") {
    
    plotGSEA_pathway(seurat, prefix, figout, out_dir, summary_type, cluster_rows,
                     cluster_cols = TRUE, zscale = TRUE, 
                     drop_sigs, save_sigs, suffix, print, 
                     geneset_sigs, clusters, long_name, input_dir, colrange)
    plotGSEA_pathway(seurat, prefix, figout, out_dir, summary_type, cluster_rows,
                     cluster_cols = FALSE, zscale = TRUE, 
                     drop_sigs, save_sigs, suffix, print, 
                     geneset_sigs, clusters, long_name, input_dir, colrange)
    
    plotGSEA_pathway(seurat, prefix, figout, out_dir, summary_type, cluster_rows,
                     cluster_cols = TRUE, zscale = FALSE, 
                     drop_sigs, save_sigs, suffix, print, 
                     geneset_sigs, clusters, long_name, input_dir, colrange = NULL)
    plotGSEA_pathway(seurat, prefix, figout, out_dir, summary_type, cluster_rows,
                     cluster_cols = FALSE, zscale = FALSE, 
                     drop_sigs, save_sigs, suffix, print, 
                     geneset_sigs, clusters, long_name, input_dir, colrange = NULL)
  
    } else if (cluster_cols == "both" & zscale == TRUE) {
    
      plotGSEA_pathway(seurat, prefix, figout, out_dir, summary_type, cluster_rows,
                     cluster_cols = TRUE, zscale = TRUE, 
                     drop_sigs, save_sigs, suffix, print, 
                     geneset_sigs, clusters, long_name, input_dir, colrange)
      plotGSEA_pathway(seurat, prefix, figout, out_dir, summary_type, cluster_rows,
                     cluster_cols = FALSE, zscale = TRUE, 
                     drop_sigs, save_sigs, suffix, print, 
                     geneset_sigs, clusters, long_name, input_dir, colrange)
      } else if (cluster_cols == "both" & zscale == FALSE) {
        
        plotGSEA_pathway(seurat, prefix, figout, out_dir, summary_type, cluster_rows,
                         cluster_cols = TRUE, zscale = FALSE, 
                         drop_sigs, save_sigs, suffix, print, 
                         geneset_sigs, clusters, long_name, input_dir, colrange = NULL)
        plotGSEA_pathway(seurat, prefix, figout, out_dir, summary_type, cluster_rows,
                         cluster_cols = FALSE, zscale = FALSE, 
                         drop_sigs, save_sigs, suffix, print, 
                         geneset_sigs, clusters, long_name, input_dir, colrange = NULL)
        
      } else if (cluster_cols == TRUE & zscale == "both") {
        
        plotGSEA_pathway(seurat, prefix, figout, out_dir, summary_type, cluster_rows,
                         cluster_cols = TRUE, zscale = TRUE, 
                         drop_sigs, save_sigs, suffix, print, 
                         geneset_sigs, clusters, long_name, input_dir, colrange)
        plotGSEA_pathway(seurat, prefix, figout, out_dir, summary_type, cluster_rows,
                         cluster_cols = TRUE, zscale = FALSE, 
                         drop_sigs, save_sigs, suffix, print, 
                         geneset_sigs, clusters, long_name, input_dir, colrange = NULL)
        
      } else if (cluster_cols == FALSE & zscale == "both") {
        
        plotGSEA_pathway(seurat, prefix, figout, out_dir, summary_type, cluster_rows,
                         cluster_cols = FALSE, zscale = TRUE, 
                         drop_sigs, save_sigs, suffix, print, 
                         geneset_sigs, clusters, long_name, input_dir, colrange)
        plotGSEA_pathway(seurat, prefix, figout, out_dir, summary_type, cluster_rows,
                         cluster_cols = FALSE, zscale = FALSE, 
                         drop_sigs, save_sigs, suffix, print, 
                         geneset_sigs, clusters, long_name, input_dir, colrange = NULL)
        
        }
  
}

#' to output tables to text, and violin plots (to be deleted?)
get_geneset_sigs <- function(geneset, seurat = cca_joint_ct1915s, clu_cols = cca_joint_pris_cols) {
  geneset_sigs <- msig_df_wide[[geneset]] %>% as.character()
  write.table(geneset_sigs, file = paste0(geneset, "_gene_list.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  geneset_sigs_inData <- geneset_sigs[which(geneset_sigs %in% rownames(seurat@data))]
  if(length(geneset_sigs_inData)/16 > 1) {
    fig_n <- length(geneset_sigs_inData)/16
    lapply(1:fig_n, function(i) {
      p <- VlnPlot(seurat, features.plot = geneset_sigs_inData[(16*(i-1)+1):min(16*i, length(geneset_sigs_inData))], cols.use = clu_cols, point.size.use = 0, do.return = TRUE)
      show(p)
    })
  } else {
    p <- VlnPlot(seurat, features.plot = geneset_sigs_inData, cols.use = clu_cols, point.size.use = 0, do.return = TRUE)
    show(p)
  }
}

#' In this function, the output folder will be created if it's not existed.
#' outputs are .txt tables and figures in both pdf and png
ranking_sigs <- function(seurat, pathway, out_dir, summary_type = "cluster_mean", prefix = NULL, alpha = 0.75, save_le = TRUE) {
  
  n_clu <- length(unique(seurat@ident))
  
  counts <- seurat %>% cytokit::meanClusterExpression()
  expr_mat <- as.matrix(seurat@data)
  
  ssgsea_le <- ssgsea_le(counts, pathway, alpha = alpha, normalize = FALSE, save_le = save_le)
  
  geneset_name <- names(pathway)
  genes_of_interest <- intersect(pathway[[geneset_name]], rownames(counts))
  expr_mat <- expr_mat[genes_of_interest, ] %>% 
    t() %>% as.data.frame() %>% # need to add as.data.frame for proceeding to add_column
    add_column(cluster = seurat@ident)
  
  dir.create(paste0(out_dir, geneset_name, "_top100"), showWarnings = FALSE)
  out_dir <- paste0(out_dir, geneset_name, "_top100", "/")
  if(is.null(prefix)) {
    prefix <- seurat@project.name
  }
  
  lapply(1:n_clu, function(i) {
    # all pathway genes (that expressed in the data)
    counts_of_interest <- counts[genes_of_interest, ]
    gene_ranking <- order(counts_of_interest[, i], decreasing = TRUE)
    ranking_df <- data.frame(gene = genes_of_interest[gene_ranking], 
                             rank = 1:length(genes_of_interest), 
                             exp = counts_of_interest[gene_ranking, i])
    write.table(ranking_df, file = glue("{out_dir}", "{prefix}", "_cluster", i-1, "_", "{geneset_name}", "_ranked_gene_list.txt"), 
                row.names = FALSE, quote = FALSE, sep = "\t")
    p <- ggplot(ranking_df[1:100, ], aes(x = rank, y = exp)) +  
      geom_point() + 
      scale_x_continuous(breaks = ranking_df$rank, 
                         labels = ranking_df$gene) + 
      ggtitle(label = glue("Top 100 ranked genes of cluster ", i-1, " in ", "{geneset_name}", " pathway")) + 
      theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1))
    save_plot(filename = glue("{out_dir}", "{prefix}", "_cluster", i-1, "_", "{geneset_name}", "_top100_gene_list.pdf"), plot = p, base_height = 4, base_aspect_ratio = 10)
    save_plot(filename = glue("{out_dir}", "{prefix}", "_cluster", i-1, "_", "{geneset_name}", "_top100_gene_list.png"), plot = p, base_height = 4, base_aspect_ratio = 10)
    
    genes_of_interest <- as.character(genes_of_interest[gene_ranking])
    
    ## boxplot (TODO)
    clu_cell_exp <- expr_mat %>% filter(cluster == i-1)
    clu_cell_exp <- clu_cell_exp[, -nrow(expr_mat)] %>% 
      t() %>% .[genes_of_interest, ] %>% as.data.frame(stringsAsFactors = FALSE)
    colnames(clu_cell_exp) <- rownames(expr_mat[which(expr_mat$cluster == i-1), ])
    
    ranking_df1 <- data.frame(gene = genes_of_interest, 
                              rank = 1:length(genes_of_interest)) %>% cbind(clu_cell_exp)
    rownames(ranking_df1) <- NULL
    tmp <- ranking_df1 %>% tidyr::gather(rank, gene, 3:length(.)) # data factor issue with gather
    
    # le genes
    le_mat <- ssgsea_le$leading_edge %>% .[[geneset_name]]
    le_genes <- rownames(le_mat[which(le_mat[, i] == 1), ]) %>% intersect(., rownames(counts))
    counts_le <- counts[le_genes, ]
    le_gene_ranking <- order(counts_le[, i], decreasing = TRUE)
    le_gene_ranking_df <- data.frame(gene = le_genes[le_gene_ranking], rank = 1:length(le_genes), exp = counts_le[le_gene_ranking, i])
    write.table(ranking_df, file = glue("{out_dir}", "{prefix}", "_cluster", i-1, "_", "{geneset_name}", "_le_ranked_gene_list.txt"), 
                row.names = FALSE, quote = FALSE, sep = "\t")
    
    p <- ggplot(ranking_df[1:100, ], aes(x = rank, y = exp)) +  
      geom_point() + 
      scale_x_continuous(breaks = ranking_df$rank, 
                         labels = ranking_df$gene) + 
      ggtitle(label = paste0("Top 100 ranked genes of cluster ", i-1, " in ", "{geneset_name}", " leading edge gene set")) + 
      theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1))
    save_plot(filename = glue("{out_dir}", "{prefix}", "_cluster", i-1, "_", "{geneset_name}", "_top100_le_gene_list.pdf"), plot = p, base_height = 4, base_aspect_ratio = 10)
    save_plot(filename = glue("{out_dir}", "{prefix}", "_cluster", i-1, "_", "{geneset_name}", "_top100_le_gene_list.png"), plot = p, base_height = 4, base_aspect_ratio = 10)
  })
}

pathway_gene_exp <- function(seurat, pathway, out_dir) {
  
  clusters <- levels(seurat@ident)
  
  counts <- suppressMessages(cytokit::meanClusterExpression(seurat))
  expr_mat <- as.matrix(seurat@data) 
  
  lapply(1:length(pathway), function(i) {
    geneset <- pathway[[i]]
    genes <- intersect(geneset, rownames(expr_mat)) 
    ## some genes in the genesets not expressed in the data
    pathway_name <- names(pathway)[i]
    dir.create(glue("{out_dir}", "{pathway_name}"), showWarnings = F)
    
    write.table(geneset, file = glue("{out_dir}", "{pathway_name}/", "the_whole_gene_list.txt"), 
                row.names = F, col.names = F, quote = F, sep = "\t")
    write.table(genes, file = glue("{out_dir}", "{pathway_name}/", "the_genes_expressed_in_data.txt"), 
                row.names = F, col.names = F, quote = F, sep = "\t")
    
    expr_mat1 <- expr_mat[genes, ] %>% 
      t() %>% as.data.frame() %>% # need to add as.data.frame for proceeding to add_column
      add_column(cluster = seurat@ident)
    
    # all pathway genes (that expressed in the data)
    counts1 <- counts[genes, ]
    write.table(counts1, file = glue("{out_dir}", "{pathway_name}/", "gene_clu_mean_exp.txt"), 
                quote = F, sep = "\t", col.names = NA)
    
    lapply(clusters, function(j) {
      
      gene_ranking <- order(counts1[, j], decreasing = T)
      ranking_df <- data.frame(gene = genes[gene_ranking], 
                               exp = counts1[gene_ranking, j])
      write.table(ranking_df, file = glue("{out_dir}", "{pathway_name}/", "gene_exp_ordered_by_cluster_{j}_mean_exp.txt"), 
                  row.names = F, quote = F, sep = "\t")
    })
    
    message("Work on ", pathway_name, " is done")
    
  })
  
}

zscoring_df <- function(df, margin = 2) {
  mat <- as.matrix(df)
  z_scores <- apply(mat, margin, function(x) {
    x <- as.numeric(x)
    z_scores <- (x - mean(x))/sd(x)
    return(z_scores)
  })
  colnames(z_scores) <- colnames(df)
  rownames(z_scores) <- rownames(df)
  return(z_scores)
}

zscoring_table <- function(prefix, 
                           input_dir = NULL, 
                           out_dir = NULL, 
                           summary_type = "cluster_mean",
                           zscale = FALSE) {
  scores <- suppressMessages(read_tsv(glue("{input_dir}", "{prefix}.{summary_type}.ssgsea_scores.txt")))
  num_col <- ncol(scores)
  
  mat <- scores %>%
    as.data.frame %>%
    tibble::column_to_rownames(var = "Signature") %>%
    t()
  
  if (zscale) {
    
    mat <- zscoring_df(mat, 2)
    
    mat %>% 
      t() %>% 
      data.frame %>% 
      magrittr::set_colnames(colnames(scores)[2:num_col]) %>%
      tibble::rownames_to_column(var = "Signature") %>% 
      {write.table(file = glue("{out_dir}", "{prefix}.{summary_type}.ssgsea_scores_zscale.txt"),
                   x = .,
                   sep = "\t",
                   quote = T,
                   row.names = FALSE)}
  }
}
