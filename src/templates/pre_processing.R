# Pre-processing ----
# |||||||||||||||||||


## Normalize data with SCT transformation ----
## |||||||||||||||||||||||||||||||||||||||||||

# SCT transformation uses regularized negative binomial regression to 
# normalize UMI count data

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# A toolkit for quality control, analysis, and exploration of single cell RNA
# sequencing data. 'Seurat' aims to enable users to identify and interpret 
# sources of heterogeneity from single cell transcriptomic measurements, and to 
# integrate diverse types of single cell data.
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# - return.only.var.genes (default is TRUE)
#   If set to TRUE the scale.data matrices in output assay are subset to contain 
#   only the variable genes.

# - variable.features.n (default is 3000)
#   Use this many features as variable features (spots) after ranking 
#   by residual variance.
# - variable.features.rv.th (default is 1.3.)
#   Instead of setting a fixed number of variable features, use this residual
#   variance cutoff; this is only used when variable.features.n is set to NULL.

return.only.var.genes <- FALSE
variable.features.n <- NULL
variable.features.rv.th <- 1.3

Seu <- CreateSeuratObject(counts = counts,
                          min.cells = min.loctions, min.features = min.features)
Seu <- SCTransform(Seu,
                   return.only.var.genes = return.only.var.genes,
                   variable.features.n = variable.features.n, 
                   variable.features.rv.th = variable.features.rv.th, 
                   verbose = VERBOSE)

counts.normalized <- as(Seu@assays$SCT@scale.data, "dgCMatrix")
names.significant_locations <- colnames(counts.normalized)
names.significant_genes <- rownames(counts.normalized)
locations.significant <- locations[names.significant_locations,]
counts.significant <- counts[names.significant_genes, names.significant_locations]

## Select spatial genes with SPARK ----
## ||||||||||||||||||||||||||||||||||||

spark = function(rawcount, location, numCores, verbose){
  location <- as.data.frame(location)
  rownames(location) = colnames(rawcount)
  
  spark <- CreateSPARKObject(counts=rawcount,
                             location=location,
                             percentage = 0.1,
                             min_total_counts = 10)
  spark@lib_size <- apply(rawcount, 2, sum)
  spark <- spark.vc(spark, 
                    covariates = NULL, 
                    lib_size = spark@lib_size, 
                    num_core = numCores,
                    fit.model="gaussian",
                    verbose = verbose)
  spark <- spark.test(spark, 
                      check_positive = T, 
                      verbose = FALSE)
  return(spark)  
}

# SPARK is an efficient tool for identifying spatial expression patterns.
# SPARK directly models raw count data generated from various spatial resolved 
# transcriptomic techniques. With a new efficient penalized quasi-likelihood 
# based algorithm, SPARK is scalable to data sets with tens of thousands of 
# genes measured on thousands of samples. Build upon a non-parametric framework, 
# SPARK-X is scalable to large-scale data sets with tens of thousands of genes 
# measured on hundred thousands of samples.

# In particular:
# - spark small sample size data for higher detection power of spatial genes
# - sparkX for large sample size data for saving time and memory

if(sparkversion=="spark"){
  spark_result <- spark(counts.significant, locations.significant,
                        numCores = numCores_spark,
                        verbose = VERBOSE)
  number.significant_genes <- sum(spark_result@res_mtest$adjusted_pvalue <= 0.05)
  names.significant_genes <- rownames(spark_result@res_mtest[order(spark_result@res_mtest$adjusted_pvalue),])[1:number.significant_genes]
}else if(sparkversion=="sparkX"){
  # Spark results
  sparkX_result <- sparkx(counts.significant, as.matrix(locations.significant),
                          numCores=numCores_spark,
                          verbose = VERBOSE)
  number.significant_genes <- sum(sparkX_result$res_mtest$adjustedPval<=0.05)
  names.significant_genes <- rownames(sparkX_result$res_mtest[order(sparkX_result$res_mtest$adjustedPval),])[1:number.significant_genes]
}

# Subset normalized data with significant spatial genes
if(is.null(number.genes)){
  if(VERBOSE){
    cat(paste("- Gene number is not specified, all", number.significant_genes, "spatially variable genes have been considered\n")) 
  }
}else {
  if(length(number.significant_genes) < number.genes){
    if(VERBOSE){
      cat("- The  number of significant spatial genes is less than the specified number of spatial genes. \n")
      cat(paste("  (Using only the", number.significant_genes, "significant spatially variable genes.)\n"))
    }
  }else{
    if(VERBOSE){
      cat(paste("- Using top",number.genes,"significant spatially variable genes. \n"))
    }
    names.significant_genes <- names.significant_genes[1:number.genes]
  }
}

counts.normalized <- counts.normalized[names.significant_genes,]
counts.significant <- counts.significant[names.significant_genes, ]


## Save results ----
## |||||||||||||||||

# Results
# - counts.significant (counts of the significant genes in the significant locations)
# - counts.normalized (normalized counts of the significant genes in the significant locations)

save(locations.significant,
     counts.normalized,
     counts.significant,
     # Saving options
     file = paste(directory.results, name.dataset, "_processed.RData", sep = ""))

## Clean ----
## ||||||||||

rm(Seu, spark_result, 
   return.only.var.genes, variable.features.n, variable.features.rv.th,
   number.significant_genes,
   names.significant_locations, names.significant_genes, 
   locations.significant, counts.normalized, counts.significant)
