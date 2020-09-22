# FROM RASflow RESULTS, GENERATES A TABLE AND A PCA PLOT WITH ALL SAMPLES
# USING TPM VALUES
# Rscript scripts/merge_tpm.R 
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
#  a barplot for contributing features of pca dimension 3 (output of fviz_contrib() function),
#  a pca plot of dim 1 and 2 (ggplot),
#  a pca plot of dim 2 and 3 (ggplot),
#  a list with clustering evaluation measures (ouput object from function cluster.stats())
#
PCA_plots <- function(res.pca, plot.title, metadata){
  var <- get_pca_var(res.pca)
  percent_var_pc1 <- round(res.pca$eig["comp 1",2], 2)
  percent_var_pc2 <- round(res.pca$eig["comp 2",2], 2)
  percent_var_pc3 <- round(res.pca$eig["comp 3",2], 2)
  
  # Diagnostic plots
  scree.plot0    <- fviz_eig(res.pca, addlabels = TRUE, rotate=T, subtitle=plot.title)
  contrib1.plot0 <- fviz_contrib(res.pca, choice = "var", axes = 1, top = 20, rotate=T, subtitle=plot.title)
  contrib2.plot0 <- fviz_contrib(res.pca, choice = "var", axes = 2, top = 20, rotate=T, subtitle=plot.title)
  contrib3.plot0 <- fviz_contrib(res.pca, choice = "var", axes = 3, top = 20, rotate=T, subtitle=plot.title)
  
  # Integration of the metadata table
  res.pca.enriched.table <-  as_tibble(res.pca$ind$coord, rownames="id") %>% 
    left_join(metadata, by=c("id"="sample")) %>% 
    column_to_rownames(var="id")
  res.pca.enriched.dist.a <- dist(res.pca.enriched.table[, c("Dim.1", "Dim.2")])
  res.pca.enriched.dist.b <- dist(res.pca.enriched.table[, c("Dim.2", "Dim.3")])
  clustering.labels <- as.numeric(as.factor(res.pca.enriched.table$group))

  # clustering evaluation
  pca.valid0a <- cluster.stats(res.pca.enriched.dist.a, clustering.labels)
  pca.valid0b <- cluster.stats(res.pca.enriched.dist.b, clustering.labels)
  
  # pca plot
  pca.plot0.data <- as_tibble(res.pca$ind$coord, rownames="id") %>% 
    left_join(metadata, by=c("id"="sample")) %>% 
    mutate(labels=paste(group, id, sep="_"))
  pca.plot0a <- pca.plot0.data %>% 
    ggplot(aes(x = Dim.1, y = Dim.2))+
    geom_point(aes(color = group), size=4) +
    geom_text_repel(aes(label=labels), size = 2) +
    labs(title=plot.title, 
         subtitle = paste0("Norm. Gamma: ", round(pca.valid0a$pearsongamma, 2), 
                           ", Dunn index: ", round(pca.valid0a$dunn, 2)) ) + 
    xlab(paste0("Dim 1 (", percent_var_pc1, "% variance)")) +
    ylab(paste0("Dim 2 (", percent_var_pc2, "% variance)")) +
    guides(size=FALSE) +
    theme_minimal() + 
    theme(panel.background = element_rect(fill = "white", colour = "grey50"))
  pca.plot0b <- pca.plot0.data %>% 
    ggplot(aes(x = Dim.3, y = Dim.2))+
    geom_point(aes(color = group), size=4) +
    geom_text_repel(aes(label=labels), size = 2) +
    labs(title=plot.title, 
         subtitle = paste0("Norm. Gamma: ", round(pca.valid0b$pearsongamma, 2), 
                           ", Dunn index: ", round(pca.valid0b$dunn, 2)) ) + 
    xlab(paste0("Dim 3 (", percent_var_pc3, "% variance)")) +
    ylab(paste0("Dim 2 (", percent_var_pc2, "% variance)")) +
    guides(size=FALSE) +
    theme_minimal() + 
    theme(panel.background = element_rect(fill = "white", colour = "grey50"))

  return(list(scree.plot0, contrib1.plot0, contrib2.plot0, contrib3.plot0, pca.plot0a, pca.plot0b, pca.valid0a, pca.valid0b))
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
BIOMART_ENS_IDS <- yaml.file$BIOMART_ENS_IDS

# CONFIG FILE AND VARIABLES
# ------------------------------------------------------------------------------
metadata <- read_tsv(METAFILE)

dir.path <- file.path(FINALOUTPUT, PROJECT, "trans", "tpmFile")
names.vector <- unique(metadata$sample)
file.names.vector <- paste0(unique(metadata$sample), "_tpm.tsv")
file.paths.vector <- file.path(dir.path, file.names.vector)

plot.title <- "PCA on TPM values"
table.file.name <- "all_samples_tpm.tsv"
pca.plot1.file.name <- "all_samples_tpm_pca.dim1_2.png"
pca.plot2.file.name <- "all_samples_tpm_pca.dim2_3.png"
pcadiag.plot.file.name <- "all_samples_tpm_pcadiag.png"

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
data <- read.delim(file.paths.vector[1], header=F) %>% as_tibble()
colnames(data) <- c("trans_id_ver", names.vector[1])

for(i in 2:length(names.vector)){
  f <- file.paths.vector[i]
  temp_table <- read.delim(f, header=F) %>% as_tibble()
  colnames(temp_table) <- c("trans_id_ver", names.vector[i])
  data <- data %>% left_join(temp_table, by = "trans_id_ver")
}
data <- data %>% 
  mutate(trans_id_ver=as.character(trans_id_ver)) %>% 
#  mutate(total = select(., sample1:sample18) %>% rowSums) %>% 
  mutate(total = select(., -trans_id_ver) %>% rowSums) %>% 
  filter(total>0) %>% 
  select(-total)

biomart <- read_tsv(BIOMART_ENS_IDS) %>% 
  select("Transcript stable ID version", "Transcript name", "Gene stable ID version", "Gene name") %>% 
  rename(trans_id_ver=`Transcript stable ID version`, trans_name=`Transcript name`,
         gene_id_ver=`Gene stable ID version`, gene_name=`Gene name`)

data %>% 
  left_join(biomart, by="trans_id_ver") %>% 
  select(trans_id_ver, trans_name, gene_id_ver, gene_name, everything()) %>% 
  write_tsv(file.path(dir.path, table.file.name))

# PCA PLOTTING
# ------------------------------------------------------------------------------
norm.table <- as.data.frame(
  data %>% 
    left_join(biomart, by="trans_id_ver") %>% 
    select(trans_id_ver, trans_name, gene_id_ver, gene_name, everything()) %>% 
    mutate(rowlabel=paste0(trans_name, "_", trans_id_ver)) %>% 
    select(-trans_id_ver, -trans_name, -gene_id_ver, -gene_name) %>% 
    column_to_rownames("rowlabel")
  )
pca.table <- log2(norm.table[which(apply(norm.table, 1, var) != 0), ]+0.00000001)
pca.table <- t(pca.table)
res.pca <- PCA(pca.table, scale.unit = TRUE, ncp = 5, graph = FALSE)

#  return(list(scree.plot0, contrib1.plot0, contrib2.plot0, contrib3.plot0, pca.plot0a, pca.plot0b, pca.valid0a, pca.valid0b))

plots_and_valid <- PCA_plots(res.pca, plot.title, metadata)
scree.plot0    <- plots_and_valid[[1]]
contrib1.plot0 <- plots_and_valid[[2]]
contrib2.plot0 <- plots_and_valid[[3]]
contrib3.plot0 <- plots_and_valid[[4]]
pca.plot0a     <- plots_and_valid[[5]]
pca.plot0b     <- plots_and_valid[[6]]
pca.valid0a    <- plots_and_valid[[7]]
pca.valid0b    <- plots_and_valid[[8]]

pca.plot0a %>% ggsave(filename=file.path(dir.path, pca.plot1.file.name),
                     units="cm",
                     width=16,
                     height=12,
                     dpi=300)

pca.plot0b %>% ggsave(filename=file.path(dir.path, pca.plot2.file.name),
                     units="cm",
                     width=16,
                     height=12,
                     dpi=300)

capture.output(suppressWarnings(suppressMessages(
  ggexport(plotlist = list(pca.plot0a, scree.plot0, contrib1.plot0, contrib2.plot0,
                           pca.plot0b, scree.plot0, contrib2.plot0, contrib3.plot0), 
           nrow = 2, ncol = 4, width=7200, height=3600, res=300, verbose=F,
           filename = file.path(dir.path, pcadiag.plot.file.name)
  )
)),  file='/dev/null')

