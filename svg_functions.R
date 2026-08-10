###################################
# Create SVG detection functions
# Daniel Israel Kakou, updated Aug, 2026
###################################

svg.spline <- function(count_in, locus_in, x_in = NULL, df = 5, block_size = NULL){
  
  ## Packages
  require(Matrix)
  require(splines)
  
  ## Convert to sparse only if needed
  if (!inherits(count_in, "dgCMatrix"))
    count_in <- Matrix(count_in, sparse = TRUE)
  
  ## Remove empty spots
  keep_spot <- Matrix::colSums(count_in) > 0
  count_in <- count_in[, keep_spot, drop = FALSE]
  locus_in <- locus_in[keep_spot,,drop=FALSE]
  
  if(!is.null(x_in))
    x_in <- x_in[keep_spot,,drop=FALSE]
  
  ## Remove genes with zero counts
  keep_gene <- Matrix::rowSums(count_in) > 0
  count_in <- count_in[keep_gene,,drop=FALSE]
  gene_names <- rownames(count_in)
  
  ## Center coordinates
  locus_in[,1] <- locus_in[,1] - mean(locus_in[,1])
  locus_in[,2] <- locus_in[,2] - mean(locus_in[,2])
  
  ## Spline basis
  Bx <- bs(locus_in[,1], df=df)
  By <- bs(locus_in[,2], df=df)
  Z <- cbind(Bx, By)
  Z <- sweep(Z, 2, colMeans(Z), "-")
  
  if(!is.null(x_in))
    Z <- cbind(x_in, Z)
  
  ## Precompute
  cholXtX <- chol(crossprod(Z))
  n <- nrow(Z)
  p <- ncol(Z)
  ngenes <- nrow(count_in)

  ## -------------------------------
  ## Automatic block size rule
  ## -------------------------------
  if (is.null(block_size)) {
    if (ngenes > 20000 ) {
      block_size <- 5000   # larger block for large matrices
    } else {
      block_size <- 1000   # default
    }
    # message(sprintf("Auto-selected block_size = %d", block_size))
  }
  
  ## -------------------------------
  ## Main computation
  ## -------------------------------
  stat <- numeric(ngenes)
  starts <- seq(1, ngenes, by = block_size)
  
  for(k in seq_along(starts)) {
    
    s <- starts[k]
    e <- min(s + block_size - 1, ngenes)
    
    block <- count_in[s:e,,drop=FALSE]
    XtY <- t(block %*% Z)
    
    beta <- backsolve(
      cholXtX,
      forwardsolve(t(cholXtX), XtY)
    )
    
    SSR <- colSums(beta * XtY)
    sy <- Matrix::rowSums(block)
    sy2 <- Matrix::rowSums(block^2)
    SST <- sy2 - sy^2/n
    
    stat[s:e] <- n * SSR / pmax(SST - SSR, .Machine$double.eps)
    
    if(k %% 20 == 0)
      gc(FALSE)
  }
  
  pval <- pchisq(stat, df = p, lower.tail = FALSE)
  
  data.frame(
    gene = gene_names,
    stat = stat,
    pvalue = pval,
    p.adj = p.adjust(pval, "BY"),
    row.names = NULL
  )
}


## Spark-x Wald and Score equivalence

spark.x <- function(count_in, locus_in, x_in = NULL,
                    method = "wald") {
  
  method <- match.arg(method, c("wald", "score"))
  
  ## Keep genes x spots
  if (inherits(count_in, "sparseMatrix")) {
    
    count <- count_in
    
  } else if (is.matrix(count_in)) {
    
    count <- Matrix::Matrix(count_in, sparse = TRUE )
    
  } else {
    
    stop("count_in must be a matrix or sparseMatrix.")
    
  }
  
  locus <- as.matrix(locus_in)
  
  # if (ncol(locus) != 2)
  #   stop("locus_in must have two columns.")

  if (ncol(count) != nrow(locus))
    stop(
      "Number of spots in count_in must equal ",
      "number of rows in locus_in."
    )
  
  
  ## Remove empty spots
  keep_spot <- Matrix::colSums(count) > 0
  
  count <- count[, keep_spot, drop = FALSE]
  locus <- locus[keep_spot, , drop = FALSE]
  
  if (!is.null(x_in))
    x_in <- as.matrix(x_in)[keep_spot, , drop = FALSE]
  
  
  ## Remove genes with zero counts
  keep_gene <- Matrix::rowSums(count) > 0
  
  count <- count[keep_gene, , drop = FALSE]
  
  gene_names <- rownames(count)
  
  if (is.null(gene_names))
    gene_names <- paste0("gene_", seq_len(nrow(count)))
  
  
  ## Number of spots

  n <- ncol(count)
  
    ## Center spatial coordinates

  S <- sweep( locus, 2, colMeans(locus),"-"  )
  
  rm(locus)
  
  ## Gene-wise centered sum of squares
  gene_sum <- Matrix::rowSums(count)
  
  gene_sum_sq <- Matrix::rowSums( count^2 )
  
  yg2 <- gene_sum_sq / n -
    (gene_sum / n)^2
  
  rm(gene_sum, gene_sum_sq)
  
  
  ## Quantiles
  qq <- apply(abs(S), 2, quantile, probs = seq(0.2, 1, by = 0.2), names = FALSE)
  
    ## Kernel calculation
  compute_stat <- function(Smat) {
    
    ## Center kernel
    Smat <- sweep(Smat, 2, colMeans(Smat), "-")
    SY <- crossprod(Smat, t(count))

    ## Solve
    R <- chol(crossprod(Smat))
    beta <- backsolve(R, forwardsolve( t(R), SY))
    
    SSR <- colSums( SY * beta )
    
    ## MSE
    if (method == "score") {
      MSE <- yg2
    } else {
      MSE <- yg2 - SSR / n
    }
    
    MSE <- pmax(MSE, .Machine$double.eps)
    
    stat <- SSR / MSE
    
    pval <- pchisq( stat, df = ncol(Smat),lower.tail = FALSE)
    
    rbind(stat = stat, pval = pval)
  }
  
  ## Projection
  res_list <- vector("list", 11)
  
  res_list[[1]] <- compute_stat(S)
  
  ## Gaussian kernels
  S2 <- S^2
  for (i in seq_len(5)) {
    S.new <- exp( sweep(-S2, 2, 2 * qq[i, ]^2, "/"))
    res_list[[i + 1]] <- compute_stat(S.new)
    rm(S.new)
  }
  
  rm(S2)
  
  
  ## Cosine kernels
  two_pi_S <- 2 * pi * S
  for (i in seq_len(5)) {
    S.new <- cos( sweep( two_pi_S, 2, qq[i, ], "/"))
    
    res_list[[i + 6]] <- compute_stat(S.new)
    rm(S.new)
  }
  
  
  ## Combine
  res <- do.call( rbind, res_list )
  
  stats <- res[ seq(1, 21, by = 2),   , drop = FALSE]
  
  pval <- res[ seq(2, 22, by = 2),  ,  drop = FALSE ]
  
  
  out.labels <- c("projection",  paste0("gaus", 1:5), paste0("cos", 1:5))
  
  rownames(stats) <- rownames(pval) <- out.labels
  
  
  ## ACAT
  comb_pval <- apply( pval, 2, function(x) {
      
      x <- x[!is.na(x)]
      
      if (length(x) == 0)
        return(NA_real_)
      ACAT::ACAT(x)
    }
  )
  
  
  ## BY adjustment
  joint_pval <- data.frame(
    combinedPval = comb_pval,
    adjustedPval = p.adjust(
      comb_pval,
      method = "BY"
    )
  )
  
    ## Return genes x kernels
  stats <- t(stats)
  pval <- t(pval)
  
  rownames(stats) <- gene_names
  rownames(pval) <- gene_names
  rownames(joint_pval) <- gene_names
  
  
  list(
    stats = stats,
    res_stest = pval,
    res_mtest = joint_pval
  )
}







