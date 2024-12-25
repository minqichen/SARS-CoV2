
if(F){
  # source activate cmq
  # export R_LIBS=/media/user/sdf/R.lib/library/wj/reg  
  # cd /media/user/sdh/cmq_wj
  # nohup Rscript gsea_R.R &
}
.libPaths("/home/wj/R/x86_64-pc-linux-gnu-library/reg")
library("utils")
library("tools")

# library("GSEA")
#library("dplyr")
#library("GSEA")

# load GSEA package
pkgs<-dir("/media/user/sda/R/GSEA/GSEA_R-master/R")
for (i in pkgs){
  source(paste0('/media/user/sda/R/GSEA/GSEA_R-master/R/',i))
  print(paste("LOADING:",i))
}

cat("\n")
cat("Starting...\n")
cat("\n")
setwd('/media/user/sdh/cmq_wj/GSEA')
path<-getwd()



# inputds_all<-dir(path)[grep(".gct.txt", dir(path))]  #前两列列名必须为NAME	Description
# inputcls<-dir(path)[grep(".cls", dir(path))]
gsdb<-dir('/media/user/sdh/cmq_wj/cellcycle_gmt')[grep(".gmt", dir('/media/user/sdh/cmq_wj/cellcycle_gmt'))]
active.ident<-c( "NK_cell" , "T_cell" , "B_cell" ,"Neutrophil","Monocyte_CCL3" , 
                 "Monocyte_Ly6c_low" , "Macrophage_Ki67" , "Monocyte_Ly6c_high" , "Macrophage_alveolar")

for (cell_group in active.ident[4:9]){
  for (gg in c('G6_G5','G7_G5','G8_G6','G8_G7')){
    inputds<-paste0(cell_group,'_',gg,'_expr.gct.txt')
    inputcls<-paste0(cell_group,'_',gg,'_pheno.cls')
    for (gmt_path in dir('/media/user/sdh/cmq_wj/cellcycle_gmt')[grep(".gmt", 
                          dir('/media/user/sdh/cmq_wj/cellcycle_gmt'))][5]){
      
      outdir <- "/media/user/sdh/cmq_wj/GSEA_out"
      
      # Enter a prefix to label output files:
      outname <-paste(cell_group,gg,gmt_path,sep='_')
      gsdb<-paste0('/media/user/sdh/cmq_wj/cellcycle_gmt/',gmt_path)
      ###
      # Collapse data set to Gene Symbols?
      collapsedataset <- FALSE
      inputchip <- "NOCHIP"
      collapsemode <- "NOCOLLAPSE"
      # Select GSEA Permutation Type (recommended: Phenotype)": c("Phenotype", "gene_set")
      reshuffetype <- "Phenotype"; permutation <- "sample.labels"
      # reshuffetype <- "gene_set"; permutation <- "gene.labels"
      
      # Use default number of permutations for significance testing? (default: 1000)
      nperms <- 1000
      # Use default maximum gene set size filter? (default: 500 genes)
      maxsize <- 2000
      # Use default minimum gene set size filter? (default: 15 genes)
      minsize <- 15
      # Use default Signal2Noise metric for ranking genes? "S2N" or "ttest"
      rankmetric <- "S2N"
      
      rankmethod <- "GSEA"
      #rankmethod <- "preranked"
      
      ###
      GSEA(
        # Input/Output Files :-------------------------------------------------------------------------------
        input.ds = inputds,                      # Input gene expression dataset file in GCT format
        input.cls = inputcls,                    # Input class vector (phenotype) file in CLS format
        gs.db = gsdb,                            # Gene set database in GMT format
        input.chip = inputchip,                  # CHIP File
        output.directory = outdir,               # Directory where to store output and results (default: "")
        # Program parameters :-------------------------------------------------------------------------------
        doc.string            = outname,         # Documentation string used as a prefix to name result files (default: "GSEA.analysis")
        reshuffling.type      = permutation,     # Type of permutation reshuffling: "sample.labels" or "gene.labels" (default: "sample.labels"
        nperm                 = as.integer(nperms),          # Number of random permutations (default: 1000)
        weighted.score.type   =  1,              # Enrichment correlation-based weighting: 0=no weight (KS), 1= weigthed, 2 = over-weigthed (default: 1)
        nom.p.val.threshold   = -1,              # Significance threshold for nominal p-vals for gene sets (default: -1, no thres)
        fwer.p.val.threshold  = -1,              # Significance threshold for FWER p-vals for gene sets (default: -1, no thres)
        fdr.q.val.threshold   = 0.25,            # Significance threshold for FDR q-vals for gene sets (default: 0.25)
        topgs                 = 20,              # Besides those passing test, number of top scoring gene sets used for detailed reports (default: 10)
        adjust.FDR.q.val      = F,               # Adjust the FDR q-vals (default: F)
        gs.size.threshold.min = as.integer(minsize),         # Minimum size (in genes) for database gene sets to be considered (default: 25)
        gs.size.threshold.max = as.integer(maxsize),         # Maximum size (in genes) for database gene sets to be considered (default: 500)
        reverse.sign          = F,               # Reverse direction of gene list (pos. enrichment becomes negative, etc.) (default: F)
        preproc.type          = 0,               # Preproc.normalization: 0=none, 1=col(z-score)., 2=col(rank) and row(z-score)., 3=col(rank). (def: 0)
        random.seed           = as.integer(as.POSIXct(Sys.time())),            # Random number generator seed. (default: 123456)
        perm.type             = 0,               # For experts only. Permutation type: 0 = unbalanced, 1 = balanced (default: 0)
        fraction              = 1.0,             # For experts only. Subsampling fraction. Set to 1.0 (no resampling) (default: 1.0)
        replace               = F,               # For experts only, Resampling mode (replacement or not replacement) (default: F)
        collapse.dataset      = collapsedataset, # Collapse dataset to gene symbols using a user provided chip file (default: F)
        collapse.mode         = collapsemode,
        save.intermediate.results = F,           # For experts only, save intermediate results (e.g. matrix of random perm. scores) (default: F)
        use.fast.enrichment.routine = T,         # Use faster routine to compute enrichment for random permutations (default: T)
        gsea.type = rankmethod,                  # Select Standard GSEA (default) or preranked
        rank.metric = rankmetric
      )
      
    }
  }
}
# inputds<-'B_cell_G6_G5_expr.gct.txt'
# inputcls<-'B_cell_G6_G5_pheno.cls'
# gsdb<-'/media/user/sdh/cmq_wj/cellcycle_gmt/cell_cycle.gmt'
# 
# 
# # Input path to experiment CLS file (or drop file into R window):
# 
# # Input path to GMT formatted gene set database (or drop file into R window):
# 
# # Drop a directory into R window to use as the output folder or enter directory path:
# outdir <- "/media/user/sdh/cmq_wj/GSEA_out"
# 
# # Enter a prefix to label output files:
# outname <- "B_cell_G6_G5_cell_cycle"
# 
# ###
# # Collapse data set to Gene Symbols?
# collapsedataset <- FALSE
# inputchip <- "NOCHIP"
# collapsemode <- "NOCOLLAPSE"
# # Select GSEA Permutation Type (recommended: Phenotype)": c("Phenotype", "gene_set")
# reshuffetype <- "Phenotype"; permutation <- "sample.labels"
# # reshuffetype <- "gene_set"; permutation <- "gene.labels"
# 
# # Use default number of permutations for significance testing? (default: 1000)
# nperms <- 1000
# # Use default maximum gene set size filter? (default: 500 genes)
# maxsize <- 2000
# # Use default minimum gene set size filter? (default: 15 genes)
# minsize <- 15
# # Use default Signal2Noise metric for ranking genes? "S2N" or "ttest"
# rankmetric <- "S2N"
# 
# rankmethod <- "GSEA"
# #rankmethod <- "preranked"
# 
# ###
# GSEA(
#   # Input/Output Files :-------------------------------------------------------------------------------
#   input.ds = inputds,                      # Input gene expression dataset file in GCT format
#   input.cls = inputcls,                    # Input class vector (phenotype) file in CLS format
#   gs.db = gsdb,                            # Gene set database in GMT format
#   input.chip = inputchip,                  # CHIP File
#   output.directory = outdir,               # Directory where to store output and results (default: "")
#   # Program parameters :-------------------------------------------------------------------------------
#   doc.string            = outname,         # Documentation string used as a prefix to name result files (default: "GSEA.analysis")
#   reshuffling.type      = permutation,     # Type of permutation reshuffling: "sample.labels" or "gene.labels" (default: "sample.labels"
#   nperm                 = as.integer(nperms),          # Number of random permutations (default: 1000)
#   weighted.score.type   =  1,              # Enrichment correlation-based weighting: 0=no weight (KS), 1= weigthed, 2 = over-weigthed (default: 1)
#   nom.p.val.threshold   = -1,              # Significance threshold for nominal p-vals for gene sets (default: -1, no thres)
#   fwer.p.val.threshold  = -1,              # Significance threshold for FWER p-vals for gene sets (default: -1, no thres)
#   fdr.q.val.threshold   = 0.25,            # Significance threshold for FDR q-vals for gene sets (default: 0.25)
#   topgs                 = 20,              # Besides those passing test, number of top scoring gene sets used for detailed reports (default: 10)
#   adjust.FDR.q.val      = F,               # Adjust the FDR q-vals (default: F)
#   gs.size.threshold.min = as.integer(minsize),         # Minimum size (in genes) for database gene sets to be considered (default: 25)
#   gs.size.threshold.max = as.integer(maxsize),         # Maximum size (in genes) for database gene sets to be considered (default: 500)
#   reverse.sign          = F,               # Reverse direction of gene list (pos. enrichment becomes negative, etc.) (default: F)
#   preproc.type          = 0,               # Preproc.normalization: 0=none, 1=col(z-score)., 2=col(rank) and row(z-score)., 3=col(rank). (def: 0)
#   random.seed           = as.integer(as.POSIXct(Sys.time())),            # Random number generator seed. (default: 123456)
#   perm.type             = 0,               # For experts only. Permutation type: 0 = unbalanced, 1 = balanced (default: 0)
#   fraction              = 1.0,             # For experts only. Subsampling fraction. Set to 1.0 (no resampling) (default: 1.0)
#   replace               = F,               # For experts only, Resampling mode (replacement or not replacement) (default: F)
#   collapse.dataset      = collapsedataset, # Collapse dataset to gene symbols using a user provided chip file (default: F)
#   collapse.mode         = collapsemode,
#   save.intermediate.results = F,           # For experts only, save intermediate results (e.g. matrix of random perm. scores) (default: F)
#   use.fast.enrichment.routine = T,         # Use faster routine to compute enrichment for random permutations (default: T)
#   gsea.type = rankmethod,                  # Select Standard GSEA (default) or preranked
#   rank.metric = rankmetric
# )
