# Load packages
library(tidyverse) 
library(limma) 
library(edgeR)
library(gt)
library(DT)
library(plotly)

# Set up your design matrix 
group <- factor(studydesign$group)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# Model mean-variance trend and fit linear model to data
v.DEGList.filtered.norm <- voom(myDGEList.filtered.norm, design, plot = TRUE)

# fit a linear model to your data
fit <- lmFit(v.DEGList.filtered.norm, design)

# get Bayesian stats for your linear model fit
ebFit <- eBayes(fits)

# TopTable to view DEGs -----
myTopHits <- topTable(ebFit, adjust ="BH", coef=1, number=20000, P.value=0.01, lfc=2, sort.by="logFC")

# convert to a tibble
myTopHits.df <- myTopHits %>%
  as_tibble(rownames = "geneID")

# Volcano Plots
vplot <- ggplot(myTopHits.df) +
  aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", geneID)) +
  geom_point(size=2) +
  geom_hline(yintercept = -log10(0.01), linetype="longdash", colour="grey", linewidth=1) +
  geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", linewidth=1) +
  geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", linewidth=1) +
  #annotate("rect", xmin = 1, xmax = 12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="#BE684D") +
  #annotate("rect", xmin = -1, xmax = -12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="#2C467A") +
  labs(title="Volcano plot",
       subtitle = "resistant_vs_sensitive",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

# making the volcano plot above interactive with plotly
ggplotly(vplot)

# decideTests to pull out the DEGs and make Venn Diagram ----
results <- decideTests(ebFit, method="global", adjust.method="BH", p.value=0.01, lfc=2)

# take a look at what the results of decideTests looks like
head(results)
summary(results)
vennDiagram(results, include="both")

# retrieve expression data for the DEGs
head(v.DEGList.filtered.norm$E)
colnames(v.DEGList.filtered.norm$E) <- sampleLabels

diffGenes <- v.DEGList.filtered.norm$E[results[,1] !=0,]
head(diffGenes)
dim(diffGenes)
#convert the DEGs to a dataframe using as_tibble
diffGenes.df <- as_tibble(diffGenes, rownames = "geneID")

# create interactive tables to display your DEGs 
datatable(diffGenes.df,
          extensions = c('KeyTable', "FixedHeader"),
          caption = 'Table 1: DEGs in resistance and sensitive cell lines',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(2:7), digits=2)

#write the DEGs to a file
write_tsv(diffGenes.df,"DiffGenes.txt")

#filter for upregulated and downregulated genes
diffGenes.up <- myTopHits.df %>% filter(logFC > 2)
diffGenes.down <- myTopHits.df %>% filter(logFC < -2)
diffGenes.up.df <- diffGenes.up %>%
  as_tibble(rownames = "geneID")
diffGenes.down.df <- diffGenes.down %>%
  as_tibble(rownames = "geneID")
library(writexl)
write_xlsx(diffGenes.up.df, "DiffGenes_upregulated.xlsx")
write_xlsx(diffGenes.down.df, "DiffGenes_downregulated.xlsx")
gt(myTopHits.df)
