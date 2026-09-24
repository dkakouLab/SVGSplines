###################################
# Create SVG detection functions
# Daniel Israel Kakou, updated Aug, 2026
###################################

# B-spline --------------------------------------------------
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
  
 
  if (!is.null(x_in)) {
    
    x_in <- as.matrix(x_in)
    
    x_in <- scale(
      x_in,
      center = TRUE,
      scale = FALSE
    )
    
    Z <- cbind( x_in, Z )
  }
  
  ## Precompute
  cholXtX <- chol(crossprod(Z))
  n <- nrow(Z)
  p <- ncol(Z)
  ngenes <- nrow(count_in)

  ## Automatic block size rule
  if (is.null(block_size)) {
    if (ngenes > 20000 ) {
      block_size <- 5000   # larger block for large matrices
    } else {
      block_size <- 1000   # default
    }
  }
  
  ## Main computation
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


# B-spline interaction -----------------------------------------

svg.spline.int <- function(count_in, locus_in, x_in = NULL, df = 5, block_size = NULL){
  
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
  

  ## Spline basis
  z1 <- bs(locus_in[,1], df=df)
  z2 <- bs(locus_in[,2], df=df)
  
  Z = model.matrix(~ -1 + z1 * z2)   # a new line added to spark.bs function
  Z <- sweep(Z, 2, colMeans(Z), FUN = "-") # another new line
  
  if (!is.null(x_in)) {
    
    x_in <- as.matrix(x_in)
    
    x_in <- scale(
      x_in,
      center = TRUE,
      scale = FALSE
    )
    
    Z <- cbind( x_in, Z )
  }
  
  ## Precompute
  cholXtX <- chol(crossprod(Z))
  n <- nrow(Z)
  p <- ncol(Z)
  ngenes <- nrow(count_in)
  
  ## Automatic block size rule
  if (is.null(block_size)) {
    if (ngenes > 20000 ) {
      block_size <- 5000   # larger block for large matrices
    } else {
      block_size <- 1000   # default
    }
  }
  
  ## Main computation
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


# B-spline cct -----------------------------------------

svg.spline.cct <- function(count_in, locus_in, df = 5) {
  
  gene_names <- rownames(count_in)
  locus_in <- as.matrix(locus_in)
  
  locus_in[, 1] <- locus_in[, 1] - mean(locus_in[, 1])
  locus_in[, 2] <- locus_in[, 2] - mean(locus_in[, 2])
  
  
  Bx <-splines::bs(locus_in[, 1], df = df)
  By <-splines::bs(locus_in[, 2], df = df)
  
  ## Center each spline basis column
  Bx <- sweep(Bx, 2, colMeans(Bx), "-")
  By <- sweep(By, 2, colMeans(By), "-")
  
  
  
  ## Prepare expression matrix
  ##
  ## count_in:  genes x locations
  ## Y:        locations x genes
  Y <- t(count_in)
  
  ## Center each gene
  Y <- sweep(Y, 2, colMeans(Y), "-")
  
  
  ## Number of observations
  n <- nrow(Y)
  
  ## Combine the 10 individual spline basis functions
  
  Z <- cbind(Bx, By)
  
  ## Perform 10 individual 1-df tests
  ##
  ## For each spline basis:
  ##
  ##     Y_g ~ Z_j
  ##
  ## The test statistic is
  ##
  ##     T_j = n * SSR_j / (SST - SSR_j)
  ##
  ## which is asymptotically chi-square(1).
  
  SST <- colSums(Y^2)
  
  stat_mat <- matrix(
    NA_real_,
    nrow = ncol(Z),
    ncol = ncol(Y)
  )
  
  p_mat <- matrix(
    NA_real_,
    nrow = ncol(Z),
    ncol = ncol(Y)
  )
  
  for (j in seq_len(ncol(Z))) {
    
    z <- Z[, j]
    
    ## z'z
    ztz <- sum(z^2)
    
    ## z'Y for every gene
    zty <- crossprod(z, Y)
    
    ## Regression sum of squares for each gene
    SSR <- as.numeric(zty^2 / ztz)
    
    ## 1-df test statistic
    stat <- n * SSR / (SST - SSR)
    
    ## Numerical protection
    stat[!is.finite(stat)] <- NA_real_
    
    ## Chi-square(1) p-value
    pval <- pchisq(
      stat,
      df = 1,
      lower.tail = FALSE
    )
    
    stat_mat[j, ] <- stat
    p_mat[j, ] <- pval
  }
  
  ##  Cauchy Combination Test
  
  ## Avoid exactly 0 or 1 because tan() can become
  ## infinite or numerically unstable.
  p_for_cct <- p_mat
  
  p_for_cct[p_for_cct <= 0] <- .Machine$double.xmin
  p_for_cct[p_for_cct >= 1] <- 1 - .Machine$double.eps
  
  ## Equal weights
  cct_stat <- colMeans(
    tan((0.5 - p_for_cct) * pi),
    na.rm = TRUE
  )
  
  ## Cauchy combined p-value
  cct_pvalue <- 0.5 - atan(cct_stat) / pi
  
  ## Numerical protection
  cct_pvalue <- pmax(
    pmin(cct_pvalue, 1),
    0
  )
  
  result <- data.frame(
    gene = gene_names,
    stat = cct_stat,
    pvalue = cct_pvalue
  )
  
  ## BY adjustment across genes
  result$p.adj <- p.adjust(
    result$pvalue,
    method = "BY"
  )
  
  ## Add the 10 individual p-values
  for (j in seq_len(ncol(Z))) {
    result[[paste0("p", j)]] <- p_mat[j, ]
  }
  
  ## Add the 10 individual statistics
  for (j in seq_len(ncol(Z))) {
    result[[paste0("stat", j)]] <- stat_mat[j, ]
  }
  
  result
}

# svg.spline.cct <- function(count_in,
#                             locus_in,
#                             x_in = NULL,
#                             df = 5,
#                             block_size = NULL) {
#   
#   require(Matrix)
#   require(splines)
#   require(ACAT)
#   
#   ##  Input checks
# 
#   if (!is.matrix(count_in) &&
#       !inherits(count_in, "sparseMatrix")) {
#     stop("count_in must be a matrix or sparseMatrix.")
#   }
#   
#   if (ncol(count_in) != nrow(locus_in)) {
#     stop(
#       "count_in must have the same number of spots ",
#       "as rows in locus_in."
#     )
#   }
#   
# 
#   ## Convert counts to sparse matrix
#   
#   if (!inherits(count_in, "dgCMatrix")) {
#     count_in <- Matrix::Matrix(
#       count_in,
#       sparse = TRUE
#     )
#   }
#   
# 
#   ##  Remove empty spots
#   
#   keep_spot <- Matrix::colSums(count_in) > 0
#   
#   count_in <- count_in[
#     ,
#     keep_spot,
#     drop = FALSE
#   ]
#   
#   locus_in <- locus_in[
#     keep_spot,
#     ,
#     drop = FALSE
#   ]
#   
#   if (!is.null(x_in)) {
#     x_in <- x_in[
#       keep_spot,
#       ,
#       drop = FALSE
#     ]
#   }
#   
#   ##  Remove zero-count genes
#   
#   keep_gene <- Matrix::rowSums(count_in) > 0
#   
#   count_in <- count_in[
#     keep_gene,
#     ,
#     drop = FALSE
#   ]
#   
#   gene_names <- rownames(count_in)
#   
#   if (is.null(gene_names)) {
#     gene_names <- paste0(
#       "gene_",
#       seq_len(nrow(count_in))
#     )
#   }
#   
# 
#   ##  Center spatial coordinates
#   
#   locus_in[, 1] <- locus_in[, 1] -
#     mean(locus_in[, 1])
#   
#   locus_in[, 2] <- locus_in[, 2] -
#     mean(locus_in[, 2])
#   
#   ##  Construct spline basis
#   
#   z1 <- splines::bs(
#     locus_in[, 1],
#     df = df
#   )
#   
#   z2 <- splines::bs(
#     locus_in[, 2],
#     df = df
#   )
#   
#   ## Center each basis
#   z1 <- sweep(
#     z1,
#     2,
#     colMeans(z1),
#     "-"
#   )
#   
#   z2 <- sweep(
#     z2,
#     2,
#     colMeans(z2),
#     "-"
#   )
#   
#   ## Combined spline basis
#   Z_spline <- cbind(z1,  z2 )
#   
#   rm(z1, z2)
#   
#   ## Number of spline basis functions
#   nspline <- ncol(Z_spline)
#   
# 
#   ##  Note about x_in
#   
#   ## This method is specifically:
#   ##
#   ## Y ~ spline_1
#   ## Y ~ spline_2
#   ## ...
#   ## Y ~ spline_10
#   ##
#   ## followed by ACAT.
#   ##
#   ## Therefore x_in is not incorporated into these individual
#   ## regressions.
#   ##
#   ## If adjustment for x_in is desired, the spline and response
#   ## need to be residualized with respect to x_in.
#   
#   if (!is.null(x_in)) {
#     warning(
#       "x_in is ignored in svg.spline.cct(). ",
#       "This method tests the spline basis functions individually ",
#       "without nuisance-covariate adjustment."
#     )
#   }
#   
#   ## Dimensions
# 
#   n <- ncol(count_in)
#   ngenes <- nrow(count_in)
#   
#   ## Automatic block size
# 
#   if (is.null(block_size)) {
#     
#     if (ngenes <= 5000) {
#       
#       block_size <- 1000
#       
#     } else if (ngenes <= 10000) {
#       
#       block_size <- 2000
#       
#     } else {
#       
#       block_size <- 5000
#     }
#   }
#   
#   ##  Gene-specific SST
# 
#   
#   gene_sum <- Matrix::rowSums(count_in  )
#   
#   gene_sum_sq <- Matrix::rowSums(count_in^2)
#   
#   ## SST = sum((Y - mean(Y))^2)
#   gene_ss <- gene_sum_sq -  gene_sum^2 / n
#   
#   ## Numerical protection
#   gene_ss <- pmax(gene_ss,0 )
#   
#   ## Genes with zero variance cannot be tested
#   valid_gene <- gene_ss > 0
#   
#   ## Precompute Z'Z for each individual spline
#   
#   ZSS <- colSums(
#     Z_spline^2
#   )
#   
#     ## 12. Output vectors
# 
#   pvalue <- rep(
#     NA_real_,
#     ngenes
#   )
#   
#   stat <- rep(
#     NA_real_,
#     ngenes
#   )
#   
# 
#   ##  Blockwise computation
# 
#   starts <- seq(
#     1,
#     ngenes,
#     by = block_size
#   )
#   
#   for (k in seq_along(starts)) {
#     
#     s <- starts[k]
#     
#     e <- min(
#       s + block_size - 1,
#       ngenes
#     )
#     
#     idx <- s:e
#     
#     ## Only genes with non-zero variance
#     idx_valid <- idx[
#       valid_gene[idx]
#     ]
#     
#     if (length(idx_valid) == 0) {
#       next
#     }
#     
#     ## Extract sparse count block
#     
#     count_block <- count_in[
#       idx_valid,
#       ,
#       drop = FALSE
#     ]
#     
#     ## Z'Y
#     ##
#     ## Z_spline:
#     ##       N × 10
#     ##
#     ## t(count_block):
#     ##       N × G
#     ##
#     ## result:
#     ##       10 × G
# 
#     
#     ZtY <- crossprod(
#       Z_spline,
#       t(count_block)
#     )
#     
#     ## Convert only the small 10 × G result to dense
#     ZtY <- as.matrix(ZtY)
#     
#     ## Gene SST
# 
#     SST <- as.numeric(
#       gene_ss[idx_valid]
#     )
#     
#     ## Calculate SSR for each individual spline
#     ##
#     ## For one predictor:
#     ##
#     ## SSR = (Z'Y)^2 / (Z'Z)
# 
#     SSR <- sweep(
#       ZtY^2,
#       1,
#       ZSS,
#       "/"
#     )
#     
# 
#     ## Chi-square statistic
# 
#     residual_ss <- sweep(
#       SSR,
#       2,
#       SST,
#       "-"
#     )
#     
#     stat_matrix <-
#       n * SSR /
#       residual_ss
#     
# 
#     ## Chi-square p-values
#     p_matrix <- pchisq(
#       stat_matrix,
#       df = 1,
#       lower.tail = FALSE
#     )
#     
# 
#     ## ACAT across 10 spline basis functions
#     
#     combined_p <- apply(
#       p_matrix,
#       2,
#       function(p) {
#         
#         p <- p[
#           is.finite(p)
#         ]
#         
#         if (length(p) == 0) {
#           return(NA_real_)
#         }
#         
#         SPARK::ACAT(p)
#       }
#     )
#     
#     ## Store results
#     
#     pvalue[idx_valid] <- combined_p
#     
#     ## Store -log10(p) as statistic
#     stat[idx_valid] <- -log10(
#       pmax(
#         combined_p,
#         .Machine$double.xmin
#       )
#     )
#     
# 
#     ## Clean temporary objects
# 
#     rm(
#       count_block,
#       ZtY,
#       SSR,
#       stat_matrix,
#       p_matrix,
#       combined_p
#     )
#     
#     if (k %% 10 == 0) {
#       gc(FALSE)
#     }
#   }
#   
#   ##  BY adjustment
#   p_adj <- rep(
#     NA_real_,
#     ngenes
#   )
#   
#   valid_p <- is.finite(
#     pvalue
#   )
#   
#   p_adj[valid_p] <- p.adjust(
#     pvalue[valid_p],
#     method = "BY"
#   )
#   
# 
#   ## Return
# 
#   result <- data.frame(
#     gene = gene_names,
#     stat = stat,
#     pvalue = pvalue,
#     p.adj = p_adj,
#     row.names = NULL
#   )
#   
#   result
# }

# B-spline kernel ------------------------------------------------------------------------

spline.kernel <- function(count_in, locus_in, x_in = NULL, df = 5) {
  
  ## count_in = genes x spots
  ## locus_in = spots x 2

  ## Packages
  require(Matrix)
  require(splines)
  
  if (!is.matrix(count_in) && !inherits(count_in, "sparseMatrix")) {
    stop(
      "count_in must be a matrix or sparseMatrix."
    )
  }
  
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
  
    ## Preserve gene names

  gene_names <- rownames(count_in)
  
  if (is.null(gene_names)) {
    gene_names <- paste0(
      "gene_",
      seq_len(nrow(count_in))
    )
  }
  
    ## Spline basis
  z1 <- splines::bs(locus_in[, 1], df = df)
  z2 <- splines::bs(locus_in[, 2], df = df)
  
  Z <- cbind(z1, z2)
  
  Z <- sweep(Z, 2,  colMeans(Z), "-")
  
  rm(z1, z2)
  
  if (!is.null(x_in)) {
    
    x_in <- as.matrix(x_in)
    
    x_in <- scale(
      x_in,
      center = TRUE,
      scale = FALSE
    )
    
    Z <- cbind( x_in, Z )
  }
  

  ## Eigenvalues
  lambda <- eigen(
    crossprod(Z),
    symmetric = TRUE,
    only.values = TRUE
  )$values
  
  
  ## Instead of:
  ## Y <- t(count_in)
  ## Y <- scale(Y)
  ## calculate the equivalent standardized cross-product
  ## without creating the huge dense Y matrix.
  n <- ncol(count_in)
  
  ## Gene means
  gene_sum <- Matrix::rowSums(count_in)
  gene_mean <- gene_sum / n

  ## Gene sum of squares
  gene_sum_sq <- Matrix::rowSums( count_in^2 )
  
  ## scale() uses:
  ## sqrt(sum((x - mean(x))^2) / (n - 1))
  gene_ss <- gene_sum_sq - n * gene_mean^2
  gene_sd <- sqrt(gene_ss / (n - 1) )
  
  ## Calculate Z'X'
  ZtX <- crossprod( Z, t(count_in))
  
  
  ## Because Z is centered:
  ## Z' (X - mean(X)) = Z'X
  ## Therefore we only need to divide by the gene SD.

  ZtY <- sweep( ZtX, 2, gene_sd, "/" )
  
  stat <- colSums(ZtY * ZtY)
  
  

  ## Davies p-values
  pvalue <- vapply(stat, function(s) {
      
      davies( s,lambda = lambda)$Qq
    },
    numeric(1)
  )
  
  ## Output 
  result <- data.frame(
    gene = gene_names,
    stat = stat,
    pvalue = pvalue,
    p.adj = p.adjust(
      pvalue,
      method = "BY"
    )
  )
  
  rownames(result) <- gene_names
  
    ## Cleanup

  rm( count_in, Z, ZtX, ZtY, gene_sum, gene_mean, gene_sum_sq, gene_ss, gene_sd)
  
  gc(FALSE)
  
  result
}



# B-spline + Kernel cct -------------------------------------------

svg.spline.kernel.cct <- function(
    count_in,
    locus_in,
    x_in = NULL,
    df = 5
) {
  
  
  ## Run SVG-Spline

  res_svg <- svg.spline(
    count_in = count_in,
    locus_in = locus_in,
    x_in = x_in,
    df = df
  )
  
  
  ## Run Spline Kernel
  res_kernel <- spline.kernel(
    count_in = count_in,
    locus_in = locus_in,
    x_in = x_in,
    df = df
  )
  
    ## Extract p-values
  p1 <- res_svg$pvalue
  p2 <- res_kernel$pvalue
  
    ## Combine p-values using ACAT package
  p_acat <- vapply(
    seq_along(p1),
    function(i) {
      
      p <- c(p1[i], p2[i])
      
      ## Remove missing values
      p <- p[is.finite(p)]
      
      ## No valid p-values
      if (length(p) == 0) {
        return(NA_real_)
      }
      
      ## Protect ACAT from exact 0 and 1
      p <- pmax(p, .Machine$double.xmin)
      p <- pmin(p, 1 - .Machine$double.eps)
      
      ## ACAT package
      ACAT::ACAT(p)
    },
    numeric(1)
  )
  
  
  ## BY adjustment
  p_adj <- rep(NA_real_, length(p_acat))
  
  valid <- is.finite(p_acat)
  
  p_adj[valid] <- p.adjust(
    p_acat[valid],
    method = "BY"
  )
  
    ## Output

  data.frame(
    gene = res_svg$gene,
    p_svg_spline = p1,
    p_spline_kernel = p2,
    pvalue = p_acat,
    p.adj = p_adj,
    row.names = NULL
  )
}

# svg.spline.kernel.cct <- function(
#     count_in,
#     locus_in,
#     x_in = NULL,
#     df = 5
# ) {
#     ## Run SVG-Spline
#   res_svg <- svg.spline(
#     count_in = count_in,
#     locus_in = locus_in,
#     x_in = x_in,
#     df = df
#   )
#   
#     ## Run spline kernel
#   res_kernel <- spline.kernel(
#     count_in = count_in,
#     locus_in = locus_in,
#     x_in = x_in,
#     df = df
#   )
#   
#   
#  
#   ## Combine p-values
#   p1 <- res_svg$pvalue
#   p2 <- res_kernel$pvalue
#   
#   pmat <- cbind(
#     p1,
#     p2
#   )
#   
#   p_acat <- apply(
#     pmat,
#     1,
#     acat_combine
#   )
#   
#     ## Output
#   
#   data.frame(
#     gene = res_svg$gene,
#     p_svg_spline = p1,
#     p_spline_kernel = p2,
#     pvalue = p_acat,
#     p.adj = p.adjust(
#       p_acat,
#       method = "BY"
#     ),
#     row.names = NULL
#   )
# }


# symmetric transformation -------------------------------

svg.sphere <- function(count_in,
                       locus_in,
                       x_in = NULL,
                       block_size = NULL) {
  
  require(Matrix)
  
 
  if (!is.matrix(count_in) &&
      !inherits(count_in, "sparseMatrix")) {
    stop("count_in must be a matrix or sparseMatrix.")
  }
  
  if (ncol(count_in) != nrow(locus_in)) {
    stop(
      "count_in must have the same number of spots ",
      "as rows in locus_in."
    )
  }
  
  if (ncol(locus_in) < 2) {
    stop(
      "locus_in must contain at least two spatial coordinates."
    )
  }
  

  ##  Convert counts to sparse

  
  if (!inherits(count_in, "dgCMatrix")) {
    count_in <- Matrix::Matrix(
      count_in,
      sparse = TRUE
    )
  }
  

  ##  Remove empty spots

  
  keep_spot <- Matrix::colSums(
    count_in
  ) > 0
  
  count_in <- count_in[
    ,
    keep_spot,
    drop = FALSE
  ]
  
  locus_in <- locus_in[
    keep_spot,
    ,
    drop = FALSE
  ]
  
  if (!is.null(x_in)) {
    
    x_in <- x_in[
      keep_spot,
      ,
      drop = FALSE
    ]
  }
  

  ## Remove zero-count genes

  
  keep_gene <- Matrix::rowSums(
    count_in
  ) > 0
  
  count_in <- count_in[
    keep_gene,
    ,
    drop = FALSE
  ]
  
  gene_names <- rownames(count_in)
  
  if (is.null(gene_names)) {
    
    gene_names <- paste0(
      "gene_",
      seq_len(nrow(count_in))
    )
  }
  
  ## 5. Normalize observed coordinates to [-1, 1]

  
  ## Coordinate 1
  x_coord <- locus_in[, 1]
  
  ## Coordinate 2
  y_coord <- locus_in[, 2]
  
  ## Scale each coordinate independently to [-1, 1]
  x_coord <- 2 * (
    x_coord - min(x_coord)
  ) / (
    max(x_coord) - min(x_coord)
  ) - 1
  
  y_coord <- 2 * (
    y_coord - min(y_coord)
  ) / (
    max(y_coord) - min(y_coord)
  ) - 1
  
  ## Convert coordinates to spherical coordinates
  
  theta <- x_coord * pi
  
  phi <- y_coord * (pi / 2)
  

  ## Construct spherical X, Y, Z

  X_matrix <- cos(phi) * cos(theta)
  
  Y_matrix <- cos(phi) * sin(theta)
  
  Z_matrix <- sin(phi)
  
  ## Combine
  S <- cbind(
    X_matrix,
    Y_matrix,
    Z_matrix
  )
  

  ## Center spatial predictors
  
  S <- sweep(
    S,
    2,
    colMeans(S),
    "-"
  )
  

  ## Add nuisance covariates if supplied
  
  ## If x_in is supplied, center it and include it in the
  ## regression model.
  ##
  ## The joint chi-square test below tests the spherical
  ## X/Y/Z components conditional on x_in.
  
  if (!is.null(x_in)) {
    
    x_in <- scale(
      x_in,
      center = TRUE,
      scale = FALSE
    )
    
    ## Full design matrix
    Z <- cbind(
      x_in,
      S
    )
    
    ## Number of nuisance covariates
    q_nuisance <- ncol(x_in)
    
  } else {
    
    Z <- S
    
    q_nuisance <- 0
  }
  
  ## Number of spherical predictors
  q_sphere <- 3
  
  ## Total number of predictors
  q_total <- ncol(Z)
  
  ## Dimensions
  n <- ncol(count_in)
  
  ngenes <- nrow(count_in)
  

  ##  Automatic block size
  
  if (is.null(block_size)) {
    
    if (ngenes <= 5000) {
      
      block_size <- 1000
      
    } else if (ngenes <= 10000) {
      
      block_size <- 2000
      
    } else {
      
      block_size <- 5000
    }
  }
  

  ## Gene SST
  
  gene_sum <- Matrix::rowSums(
    count_in
  )
  
  gene_sum_sq <- Matrix::rowSums(
    count_in^2
  )
  
  gene_ss <- gene_sum_sq -
    gene_sum^2 / n
  
  ## Numerical protection
  gene_ss <- pmax(
    gene_ss,
    0
  )
  
  ## Zero-variance genes
  valid_gene <- gene_ss > 0
  

  ##  Fit full model once through cross-products
  
  XtX <- crossprod(Z)
  
  cholXtX <- chol(
    XtX
  )
  

  ##  Output
  pvalue <- rep(
    NA_real_,
    ngenes
  )
  
  stat <- rep(
    NA_real_,
    ngenes
  )
  

  ## Blockwise regression
  
  starts <- seq(
    1,
    ngenes,
    by = block_size
  )
  
  for (k in seq_along(starts)) {
    
    s <- starts[k]
    
    e <- min(
      s + block_size - 1,
      ngenes
    )
    
    idx <- s:e
    
    idx_valid <- idx[
      valid_gene[idx]
    ]
    
    if (length(idx_valid) == 0) {
      next
    }
    

    ## Count block
    
    count_block <- count_in[
      idx_valid,
      ,
      drop = FALSE
    ]
    
    ## Z'Y
    ZtY <- crossprod(
      Z,
      t(count_block)
    )
    
    ZtY <- as.matrix(
      ZtY
    )
    

    ## Regression coefficients
    beta <- backsolve(
      cholXtX,
      forwardsolve(
        t(cholXtX),
        ZtY
      )
    )
    

    ## SSR for full model

    SSR_full <- colSums(
      ZtY * beta
    )
    

    ## If nuisance covariates are present, we need the
    ## incremental SSR for X/Y/Z beyond the nuisance model.

    if (q_nuisance > 0) {
      
      Z_nuis <- Z[, seq_len(q_nuisance), drop = FALSE]
      
      XtX_nuis <- crossprod(
        Z_nuis
      )
      
      chol_nuis <- chol(
        XtX_nuis
      )
      
      ZtY_nuis <- crossprod(
        Z_nuis,
        t(count_block)
      )
      
      beta_nuis <- backsolve(
        chol_nuis,
        forwardsolve(
          t(chol_nuis),
          ZtY_nuis
        )
      )
      
      SSR_nuis <- colSums(
        ZtY_nuis * beta_nuis
      )
      
      ## Additional explained variation from spherical terms
      SSR_sphere <- SSR_full -
        SSR_nuis
      
      rm(
        Z_nuis,
        XtX_nuis,
        chol_nuis,
        ZtY_nuis,
        beta_nuis,
        SSR_nuis
      )
      
    } else {
      
      ## No nuisance covariates
      SSR_sphere <- SSR_full
    }
    
    ## Chi-square statistic

    denominator <- pmax(
      gene_ss[idx_valid] -
        SSR_full,
      .Machine$double.eps
    )
    
    stat_block <- n *
      SSR_sphere /
      denominator
    

    ## 3-df chi-square p-value

    p_block <- pchisq(
      stat_block,
      df = q_sphere,
      lower.tail = FALSE
    )
    
    ## Store

    stat[idx_valid] <- stat_block
    
    pvalue[idx_valid] <- p_block
    

    ## Clean
    rm(
      count_block,
      ZtY,
      beta,
      SSR_full,
      SSR_sphere,
      denominator,
      stat_block,
      p_block
    )
    
    if (k %% 10 == 0) {
      gc(FALSE)
    }
  }
  

  ## BY adjustment
  
  p_adj <- rep(
    NA_real_,
    ngenes
  )
  
  valid_p <- is.finite(
    pvalue
  )
  
  p_adj[valid_p] <- p.adjust(
    pvalue[valid_p],
    method = "BY"
  )
  
  ## Return

  
  result <- data.frame(
    gene = gene_names,
    stat = stat,
    pvalue = pvalue,
    p.adj = p_adj,
    row.names = NULL
  )
  
  result
}

# Spark-x Wald and Score equivalence -----------------------------------

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







