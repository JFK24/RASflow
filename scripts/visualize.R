# 2020-06-30: Use local Biomart annotation file (remove EnsDb.Hsapiens.v86)
# 2020-03-17: Use log2 values for heatmap and pca
#             Fix cutoff lines on volcano plots for adjusted p values
# 2020-03-16: Add pca plots and bug fixes (free axis limits of volcano plot)
#             Remove dependency to online resources (ensembldb)
#             but it is now restricted to human data (EnsDb.Hsapiens.v86)!
# To do: fix restriction to human data by reading in a gene annotation file

# ==============================================================================
# load the libraries
# ==============================================================================
if (!require("plotscale")) install.packages('scripts/plotscale_0.1.6.tar.gz', repos = NULL, type="source")
suppressMessages(library(yaml, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(hash, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(mygene, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(EnhancedVolcano, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(plotscale, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(GenomicFeatures, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(FactoMineR, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(factoextra, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(ggpubr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(readr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(dplyr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(stringr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

# ==============================================================================
#  load parameters in config file
# ==============================================================================

# passing the params from command line
args <- commandArgs(TRUE)
norm.path <- args[1] # norm.path <- 'output/test/trans/dea/countGroup'
dea.path <- args[2]  # dea.path <- 'output/test/trans/dea/DEA/gene-level'
out.path <- args[3]  # out.path <- 'output/test/trans/dea/visualization'

# norm.path <- '../../RASflowResults/MZB_WT_KO_2/trans/dea/countGroup'
# dea.path <- '../../RASflowResults/MZB_WT_KO_2/trans/dea/DEA/gene-level'
# out.path <- '../../RASflowResults/MZB_WT_KO_2/trans/dea/visualization'
# name.control <- "WT"
# name.treat <- "KO"

# load the config file
if (length(args) > 3) {  # this script is used in  visualize_test.rules
    yaml.file <- yaml.load_file(args[4])
} else {  # this script is used in visualize.rules
    yaml.file <- yaml.load_file('configs/config_main.yaml')
}

# extract the information from the yaml file
controls <- yaml.file$CONTROL
treats <- yaml.file$TREAT
dea.tool <- yaml.file$DEATOOL
ensembl_annotation <- yaml.file$ANNOTATION
ensembl_ids  <- as_tibble(read.delim(ensembl_annotation, header=F,  comment.char = "#", stringsAsFactors=F)) %>% 
  filter(V3=="transcript") %>% 
  mutate(gene_id=str_extract(V9, regex("(?<=gene_id )(.+?)(?=;)", dotall = TRUE))) %>% 
  mutate(gene_name=str_extract(V9, regex("(?<=gene_name )(.+?)(?=;)", dotall = TRUE))) %>% 
  mutate(trans_id=str_extract(V9, regex("(?<=transcript_id )(.+?)(?=;)", dotall = TRUE))) %>% 
  mutate(trans_name=str_extract(V9, regex("(?<=transcript_name )(.+?)(?=;)", dotall = TRUE))) %>% 
  select(gene_id, trans_id, gene_name, trans_name)

# check the number of comparisons
num.control <- length(controls)  # number of comparisons that the user wants to do
num.treat <- length(treats)  # should equals to num.control

if (num.control != num.treat) {
  message("Error: Control groups don't mathch with treat groups!")
  message("Please check config_dea.yaml")
  quit(save = 'no')
}

num.comparison <- num.control

# ==============================================================================
# function to plot volcano plot and heatmap
# ==============================================================================
plot.volcano.heatmap <- function(name.control, name.treat, dea.path, norm.path, ensembl_ids) {
  file.dea.table <- paste(dea.path, "/dea_", name.control, "_", name.treat, ".tsv", sep = "")
  norm.control <- paste(norm.path, "/", name.control, "_gene_norm.tsv", sep = "")  # normalized table of control
  norm.treat <- paste(norm.path, "/", name.treat, "_gene_norm.tsv", sep = "")  # normalized table of treat

  dea.table <- read.table(file.dea.table, header = TRUE, row.names = 1)
  # sort the dea table: ascending of FDR then descending of absolute valued of logFC
  if (dea.tool == 'edgeR') {
    dea.table <- dea.table[order(dea.table$FDR, -abs(dea.table$logFC), decreasing = FALSE), ]  
  } else if (dea.tool == 'DESeq2') {
    dea.table <- dea.table[order(dea.table$padj, -abs(dea.table$log2FoldChange), decreasing = FALSE), ]
  }

  gene.id.dea <- row.names(dea.table)

  if (length(intersect(row.names(dea.table), ensembl_ids$gene_id))>0){
    gene.dea <- data.frame(id=gene.id.dea, stringsAsFactors=FALSE) %>% 
      left_join(ensembl_ids %>% group_by(gene_id) %>% top_n(1, trans_id), by=c("id" = "gene_id")) %>% 
      select(id, gene_name) %>% 
      rename(GENEID=id, SYMBOL=gene_name)
  } else {
    gene.dea <- data.frame(id=gene.id.dea, stringsAsFactors=FALSE) %>% 
      left_join(ensembl_ids %>% group_by(trans_id) %>% top_n(1, gene_id), by=c("id" = "trans_id")) %>% 
      select(id, trans_name) %>% 
      rename(GENEID=id, SYMBOL=trans_name)
  }

  # volcano plot
  #-----------------------------------------------------------------------------
  if (dea.tool == 'edgeR') {
    # fig.volcano <- EnhancedVolcano(dea.table, lab = gene.dea, xlab = bquote(~Log[2]~ "fold change"), x = 'logFC', y = 'FDR', pCutoff = 10e-5, col = c("grey30", "orange2", "royalblue", "red2"),
    #                                FCcutoff = 1, xlim = c(-5, 5), ylim = c(0, 10), transcriptPointSize = 1.5, title = NULL, subtitle = NULL)  
    fig.volcano <- EnhancedVolcano(dea.table, lab = gene.dea$SYMBOL, xlab = bquote(~Log[2]~ "fold change"), x = 'logFC', y = 'FDR', pCutoff = 1e-2, col = c("grey30", "orange2", "royalblue", "red2"),
                                   FCcutoff = 1, transcriptPointSize = 1.5, title = NULL, subtitle = NULL)  
  } else if (dea.tool == 'DESeq2') {
    # fig.volcano <- EnhancedVolcano(dea.table, lab = gene.dea, xlab = bquote(~Log[2]~ "fold change"), x = 'log2FoldChange', y = 'padj', pCutoff = 10e-5, col = c("grey30", "orange2", "royalblue", "red2"),
    #                                FCcutoff = 1, xlim = c(-5, 5), ylim = c(0, 10), transcriptPointSize = 1.5, title = NULL, subtitle = NULL)
    fig.volcano <- EnhancedVolcano(dea.table, lab = gene.dea$SYMBOL, xlab = bquote(~Log[2]~ "fold change"), x = 'log2FoldChange', y = 'padj', pCutoff = 1e-2, col = c("grey30", "orange2", "royalblue", "red2"),
                                   FCcutoff = 1, transcriptPointSize = 1.5, title = NULL, subtitle = NULL)
  }
  
  as.pdf(fig.volcano, width = 9, height = 6, scaled = TRUE, file = file.path(out.path, paste('volcano_plot_', name.control, '_', name.treat, '.pdf', sep = '')))
  png(file.path(out.path, paste('volcano_plot_', name.control, '_', name.treat, '.png', sep = '')), , width = 20, height = 16, units="cm", res=300)
  fig.volcano
  dev.off()
  
  # heatmap
  #-----------------------------------------------------------------------------
  norm.table.control <- read.table(norm.control, header = TRUE, row.names = 1)
  norm.table.treat <- read.table(norm.treat, header = TRUE, row.names = 1)

  num.control <- dim(norm.table.control)[2]
  num.treat <- dim(norm.table.treat)[2]

  norm.table <- cbind(norm.table.control, norm.table.treat)
  groups <- c(name.control, name.treat)

  # instead using all genes, only use the top 20 genes in dea.table
  index.deg <- which(row.names(norm.table) %in% gene.id.dea[1:20])
  norm.table.deg <- norm.table[index.deg,]

  gene.id.norm.table <- rownames(norm.table.deg)
  gene.symbol.norm.table <- data.frame(GENEID=gene.id.norm.table, stringsAsFactors=FALSE) %>% 
    left_join(gene.dea, by=c("GENEID"))
  gene.symbol.norm.table <- gene.symbol.norm.table$SYMBOL

  # if can't find a symbol for the id, then keep the id as it is
  gene.norm.table <- gene.symbol.norm.table
  for (i in c(1:length(gene.norm.table))) {
    if (is.na(gene.norm.table[i])) {
      gene.norm.table[i] <- gene.id.norm.table[i]
    }
  }

  palette <- c("#999999", "#377EB8")
  palette.group <- c(rep(palette[1], num.control), rep(palette[2], num.treat))
  palette.group.names <- c(rep(name.control, num.control), rep(name.treat, num.treat))

  ## draw heatmap
  pdf(file = file.path(out.path, paste('heatmap_', name.control, '_', name.treat, '.pdf', sep = '')), width = 15, height = 15, title = 'Heatmap using the top features')
#  heatmap(as.matrix(norm.table.deg), ColSideColors = palette.group, margins = c(15,15), labRow = gene.norm.table, cexRow = 1.9, cexCol = 1.9)
  heatmap(as.matrix(log2(norm.table.deg+0.00000001)), ColSideColors = palette.group, margins = c(15,15), labRow = gene.norm.table, cexRow = 1.9, cexCol = 1.9)
  legend("topleft", title = 'Group', legend=groups, text.font = 15,
         col = palette, fill = palette, cex=1.8)
  dev.off()

  png(file = file.path(out.path, paste('heatmap_', name.control, '_', name.treat, '.png', sep = '')), width = 1500, height = 1500, title = 'Heatmap using the top features')
#  heatmap(as.matrix(norm.table.deg), ColSideColors = palette.group, margins = c(15,15), labRow = gene.norm.table, cexRow = 1.9, cexCol = 1.9)
  heatmap(as.matrix(log2(norm.table.deg+0.00000001)), ColSideColors = palette.group, margins = c(15,15), labRow = gene.norm.table, cexRow = 1.9, cexCol = 1.9)
  legend("topleft", title = 'Group', legend=groups, text.font = 15,
         col = palette, fill = palette, cex=1.8)
  dev.off()

  # PCA
  #-----------------------------------------------------------------------------
#  http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/
#  pca.table <- norm.table[which(apply(norm.table, 1, var) != 0), ]
  pca.table <- log2(norm.table[which(apply(norm.table, 1, var) != 0), ]+0.00000001)
  pca.table.ids <- rownames(pca.table)
  pca.table.symbols <- data.frame(GENEID=pca.table.ids, stringsAsFactors=FALSE) %>% 
    left_join(gene.dea, by=c("GENEID"))
  pca.table.symbols <- pca.table.symbols$SYMBOL
  rownames(pca.table) <- paste(pca.table.symbols, pca.table.ids, sep="_")
  pca.table <- t(pca.table)
  
  res.pca <- PCA(pca.table, scale.unit = TRUE, ncp = 5, graph = FALSE)
  var <- get_pca_var(res.pca)

  # Screeplot
  scree.plot <- fviz_eig(res.pca, addlabels = TRUE, rotate=T)
  # Contributions of variables to PC1
  contrib1.plot <- fviz_contrib(res.pca, choice = "var", axes = 1, top = 20, rotate=T)
  # Contributions of variables to PC2
  contrib2.plot <- fviz_contrib(res.pca, choice = "var", axes = 2, top = 20, rotate=T)
  # pca plot
  pca.plot <- fviz_pca_ind(res.pca,
                col.ind = palette.group.names, # color by groups
                palette = palette,
                addEllipses = TRUE, # Concentration ellipses
                legend.title = "Groups",
                repel=TRUE
               )

  ggexport(plotlist = list(scree.plot, contrib1.plot, contrib2.plot), 
         nrow = 3, ncol = 1, width=10, height=12,
         filename = file.path(out.path, paste('pca_diag_', name.control, '_', name.treat, '.pdf', sep = ''))
         )
  ggexport(plotlist = list(pca.plot), 
         nrow = 3, ncol = 1,
         filename = file.path(out.path, paste('pca_', name.control, '_', name.treat, '.pdf', sep = ''))
         )
  ggexport(plotlist = list(scree.plot, contrib1.plot, contrib2.plot), 
         nrow = 3, ncol = 1, width=1000, height=1200,
         filename = file.path(out.path, paste('pca_diag_', name.control, '_', name.treat, '.png', sep = ''))
         )
  ggexport(plotlist = list(pca.plot), 
         nrow = 3, ncol = 1, width=1000, height=1000,
         filename = file.path(out.path, paste('pca_', name.control, '_', name.treat, '.png', sep = ''))
         )
}

# ==============================================================================
# MAIN FUNCTION
# ==============================================================================
for (ith.comparison in c(1:num.comparison)) {
  name.control <- controls[ith.comparison]
  name.treat <- treats[ith.comparison]
  plot.volcano.heatmap(name.control, name.treat, dea.path, norm.path, ensembl_ids)
}
