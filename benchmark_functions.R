###################################
# Create benchmarking functions
# Daniel Israel Kakou, updated Aug, 2026
###################################


## Evaluate power and Type 1 based on groundtruth SVG and non-SVG
evaluate_svg <- function(results,
                         alpha = 0.05,
                         adjust = FALSE){
  
  out <- lapply(names(results), function(nm){
    
    dat <- results[[nm]]
    
    methods <- c( "p_splines",
                  "p_spark_wald",
                  "p_spark_score",
                  "p_SPARKX")
    
    do.call(rbind,
            lapply(methods, function(m){
              
              p <- dat[[m]]
              
              if(adjust)
                p <- p.adjust(p, method = "BY")
              
              reject <- p < alpha
              
              data.frame(
                
                Dataset = sub( "\\.RData$", "", nm),
                
                Method = sub("^p_", "", m),
                
                Power =
                  mean(reject[ dat$is.de ]),
                
                Type1 =
                  mean(reject[ !dat$is.de ]),
                
                stringsAsFactors = FALSE
                
              )
              
            })
            
    )
    
  })
  
  do.call(rbind, out)
  
}



## To Plot TPR vs FDR

compute_tpr_fdr <- function(dat,
                             method,
                             alpha = c(0.01, 0.05, 0.10)){
  
  padj <- p.adjust(dat[[method]], method = "BY")
  
  out <- data.frame()
  
  for(a in alpha){
    
    sig <- padj <= a
    
    TP <- sum(sig & dat$is.de)
    FP <- sum(sig & !dat$is.de)
    FN <- sum(!sig & dat$is.de)
    
    out <- rbind(out,
                 data.frame(
                   Cutoff = factor(a,
                                   levels = c(0.01,0.05,0.10)),
                   TPR = TP/(TP+FN),
                   FDR = ifelse(TP+FP==0,0,FP/(TP+FP))
                 ))
  }
  
  out
}



pvalue_summary_all <- function(results){
  
  method_names <- c(
    p_splines = "Splines",
    p_spark_wald  = "Wald",
    p_spark_score = "Score",
    p_SPARKX = "SPARK-X"
  )
  
  out <- list()
  
  for(ds in names(results)){
    
    dat <- results[[ds]]
    
    for(m in names(method_names)){
      
      padj <- p.adjust(dat[[m]], method = "BY")
      
      grp <- cut(
        padj,
        breaks = c(-Inf,0,1e-10,0.01,0.05,Inf),
        labels = c(
          "0",
          "(0,1e-10]",
          "(1e-10,0.01]",
          "(0.01,0.05]",
          ">0.05"
        )
      )
      
      tab <- prop.table(table(grp))*100
      
      out[[length(out)+1]] <-
        data.frame(
          Dataset = sub("\\.RData$","",ds),
          Method = method_names[m],
          Category = factor(
            names(tab),
            levels = c(
              "0",
              "(0,1e-10]",
              "(1e-10,0.01]",
              "(0.01,0.05]",
              ">0.05"
            )
          ),
          Percent = as.numeric(tab)
        )
      
    }
    
  }
  
  bind_rows(out)
  
}



plotUpSet <- function(dat, dataset_name = "", alpha = 0.05){
  
  svg <- data.frame(
    svg.splines    = p.adjust(dat$p_splines, "BY") < alpha, #& dat$is.de,
    spark.x.wald    = p.adjust(dat$p_spark_wald, "BY") < alpha, #& dat$is.de,
    spark.x.score   = p.adjust(dat$p_spark_score, "BY") < alpha, #& dat$is.de,
    SPARKX  = p.adjust(dat$p_SPARKX, "BY") < alpha #& dat$is.de
  )
  
  ComplexUpset::upset(
    svg,
    intersect = c("svg.splines", "spark.x.wald", "spark.x.score", "SPARKX"),
    name = "SVGs"
  ) +
    ggtitle(dataset_name)
}

plotVenn <- function(dat,
                     dataset_name = "",
                     alpha = 0.05){
  
  sets <- list(
    Splines = rownames(dat)[p.adjust(dat$p_splines, "BY") < alpha], 
    Wald    = rownames(dat)[p.adjust(dat$p_spark_wald, "BY") < alpha],
    Score   = rownames(dat)[p.adjust(dat$p_spark_score, "BY") < alpha],
    SPARKX =  rownames(dat)[p.adjust(dat$p_SPARKX, "BY") < alpha]
    
  )
  
  ggVennDiagram(sets, label_alpha = 0) +
    scale_fill_gradient(low = "#fee08b", high = "#d53e4f") +
    labs(title = dataset_name) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
}


plotJaccard <- function(dat,
                        dataset_name = "",
                        alpha = 0.05){
  
  ## Significant genes
  svg <- list(
    Splines = p.adjust(dat$p_splines, "BY") < alpha,
    Wald    = p.adjust(dat$p_spark_wald, "BY") < alpha,
    Score   = p.adjust(dat$p_spark_score, "BY") < alpha,
    SPARKX  = p.adjust(dat$p_SPARKX, "BY") < alpha
  )
  
  methods <- names(svg)
  J <- matrix(0,
              length(methods),
              length(methods),
              dimnames = list(methods, methods))
  
  for(i in seq_along(methods)){
    for(j in seq_along(methods)){
      
      A <- svg[[i]]
      B <- svg[[j]]
      
      J[i,j] <- sum(A & B) / sum(A | B)
    }
  }
  
  pheatmap(
    J,
    display_numbers = TRUE,
    number_format = "%.2f",
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    main = dataset_name
  )
}


plotCompare <- function(dat,
                        method1,
                        method2,
                        truth = NULL,
                        main = NULL){
  
  x <- -log10(dat[[method1]])
  y <- -log10(dat[[method2]])
  
  keep <- is.finite(x) & is.finite(y)
  
  x <- x[keep]
  y <- y[keep]
  
  ## Same scale on both axes
  lim <- range(c(x, y))
  
  if(!is.null(truth))
    truth <- truth[keep]
  
  plot(x, y,
       pch = 19,
       cex = 0.35,
       col = if(is.null(truth))
         "black"
       else
         ifelse(truth, "red", "grey70"),
       xlab = paste0("-log10(", method1, ")"),
       ylab = paste0("-log10(", method2, ")"),
       xlim = lim,
       ylim = lim,
       main = main)
  
  abline(0,1,lwd=1)
}


