# FROM RASflow RESULTS, GENERATES A TABLE AND A PCA PLOT WITH ALL SAMPLES
# USING DESEQ2 NORMALIZED VALUES
# Rscript scripts/merge_deseq2.R 
# ------------------------------------------------------------------------------
suppressMessages(library(yaml, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(tidyverse, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(FactoMineR, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(factoextra, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(ggpubr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(ggrepel, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(fpc, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

# PCA PLOTS FUNCTION
# ------------------------------------------------------------------------------
# FROM RESULTS OF A PCA ANALYSIS, CREATES PLOTS AND CLUSTERING EVALUATION ON THE 
# 2 FIRST COMPONANTS
#
# INPUT:
# res.pca: result object of the PCA() function for a table T (features and samples)
# plot.title: string to be used as title
# metadata: a data frame annotating the samples in table T: a column named sample (sample ids, chr) and a column named group (group annotation, chr)
#
# OUTPUT: a list containing 
#  a screeplot (output of fviz_eig() function), 
#  a barplot for contributing features of pca dimension 1 (output of fviz_contrib() function),
#  a barplot for contributing features of pca dimension 2 (output of fviz_contrib() function),
#  a pca plot (ggplot),
#  a list with clustering evaluation measures (ouput object from function cluster.stats())
#
#
PCA_plots <- function(res.pca, plot.title, metadata){
  var <- get_pca_var(res.pca)
  
  # Diagnostic plots
  scree.plot0    <- fviz_eig(res.pca, addlabels = TRUE, rotate=T, subtitle=plot.title)
  contrib1.plot0 <- fviz_contrib(res.pca, choice = "var", axes = 1, top = 20, rotate=T, subtitle=plot.title)
  contrib2.plot0 <- fviz_contrib(res.pca, choice = "var", axes = 2, top = 20, rotate=T, subtitle=plot.title)
  
  # Integration of the metadata table
  res.pca.enriched.table <-  as_tibble(res.pca$ind$coord, rownames="id") %>% 
    left_join(metadata, by=c("id"="sample")) %>% 
    column_to_rownames(var="id") %>%
    mutate(group=ifelse(is.na(group), "NA", group))
  res.pca.enriched.dist <- dist(res.pca.enriched.table[, c("Dim.1", "Dim.2")])
  clustering.labels <- as.numeric(as.factor(res.pca.enriched.table$group))
  
  # clustering evaluation
  pca.valid0 <- cluster.stats(res.pca.enriched.dist, clustering.labels)
  
  # pca plot
  pca.plot0.data <- as_tibble(res.pca$ind$coord, rownames="id") %>% 
    left_join(metadata, by=c("id"="sample")) %>% 
    mutate(labels=paste(group, id, sep="_"))
  pca.plot0 <- pca.plot0.data %>% 
    ggplot(aes(x = Dim.1, y = Dim.2))+
    geom_point(aes(color = group), size=5) +
    geom_text_repel(aes(label=labels), size = 2) +
    labs(title=plot.title, 
         subtitle = paste0("Norm. Gamma: ", round(pca.valid0$pearsongamma, 2), 
                           ", Dunn index: ", round(pca.valid0$dunn, 2)) ) + 
    guides(size=FALSE) +
    theme_minimal() + 
    theme(panel.background = element_rect(fill = "white", colour = "grey50"))
  
  return(list(scree.plot0, contrib1.plot0, contrib2.plot0, pca.plot0, pca.valid0))
}

# CONFIG FILE AND VARIABLES
# ------------------------------------------------------------------------------
yaml.file <- yaml.load_file('configs/config_main.yaml')
PROJECT <- yaml.file$PROJECT
FINALOUTPUT <- yaml.file$FINALOUTPUT
DEA <- yaml.file$DEA
METAFILE <- yaml.file$METAFILE
OUTPUTPATH <- yaml.file$OUTPUTPATH
DEATOOL <- yaml.file$DEATOOL
GENE_LEVEL <- yaml.file$GENE_LEVEL

# CONFIG FILE AND VARIABLES
# ------------------------------------------------------------------------------
metadata <- read_tsv(METAFILE)
dir.path <- file.path(FINALOUTPUT, PROJECT, "trans", "dea", "countGroup")

extensions <- c("_gene_norm", "_gene_abundance", "_trans_norm", "_trans_abundance")

for (ext in extensions) {
  
  file.names.vector <- paste0(unique(metadata$group), ext, ".tsv")
  file.paths.vector <- file.path(dir.path, file.names.vector)
  
  plot.title <- "PCA on Gene Abundance"
  if(grepl("trans", ext)){ plot.title <- gsub("Gene", "Transcript", plot.title) }
  if(grepl("norm", ext)){ plot.title <- gsub("Abundance", "Normalized Values", plot.title) }

  table.file.name <- paste0("all_samples", ext, ".tsv")
  pca.plot.file.name <- paste0("all_samples", ext, "_pca.png")
  pcadiag.plot.file.name <- paste0("all_samples", ext, "_pcadiag.png")

  # CONTROLS
  # ------------------------------------------------------------------------------
  if(!DEA | !GENE_LEVEL) {
    q("Nothing to do: options DEA and GENE_LEVEL not 'TRUE'")
  }
  if(length(file.names.vector)<2){
    q("Nothing to do: only 1 file to merge")
  }
  
  # MERGE TABLES AND SAVE
  # ------------------------------------------------------------------------------
  data <- read.delim(file.paths.vector[1], check.names=FALSE) %>% rownames_to_column("ensembl_id") %>% as_tibble()
  
  for(f in file.paths.vector[-1]){
    temp_table <- read.delim(f, check.names=FALSE) %>% rownames_to_column("ensembl_id") %>% as_tibble()
    data <- data %>% left_join(temp_table, by = "ensembl_id")
  }
  
  data %>% write_tsv(file.path(dir.path, table.file.name))

  # PCA PLOTTING
  # ------------------------------------------------------------------------------
  norm.table <- as.data.frame(data %>% column_to_rownames("ensembl_id"))
  pca.table <- log2(norm.table[which(apply(norm.table, 1, var) != 0), ]+0.00000001)
  pca.table <- t(pca.table)
  res.pca <- PCA(pca.table, scale.unit = TRUE, ncp = 5, graph = FALSE)
  plots_and_valid <- PCA_plots(res.pca, plot.title, metadata)
  scree.plot0    <- plots_and_valid[[1]]
  contrib1.plot0 <- plots_and_valid[[2]]
  contrib2.plot0 <- plots_and_valid[[3]]
  pca.plot0      <- plots_and_valid[[4]]
  pca.valid0     <- plots_and_valid[[5]]
  
  pca.plot0 %>% ggsave(filename=file.path(dir.path, pca.plot.file.name),
                       units="cm",
                       width=16,
                       height=12,
                       dpi=300)
  
  capture.output(suppressWarnings(suppressMessages(
    ggexport(plotlist = list(pca.plot0, scree.plot0, contrib1.plot0, contrib2.plot0), 
             nrow = 2, ncol = 2, width=3600, height=3600, res=300, verbose=F,
             filename = file.path(dir.path, pcadiag.plot.file.name)
    )
  )),  file='/dev/null')

}
