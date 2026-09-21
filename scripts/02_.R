# Load packages
library(tidyverse)
library(edgeR)
library(matrixStats)
library(cowplot)

# Examine and create summary stats for your data
myTPM <- Txi_gene$abundance
myCounts <- Txi_gene$counts
colSums(myTPM)
colSums(myCounts)
myTPM.stats <- transform(myTPM, 
                         SD=rowSds(myTPM), 
                         AVG=rowMeans(myTPM),
                         MED=rowMedians(myTPM))
head(myTPM.stats) # look at what you created  

# make a DGEList from the gene-level counts
myDGEList <- DGEList(Txi_gene$counts)
sf<-myDGEList$counts
view(sf)
save(myDGEList, file = "myDGEList")
load(file = "myDGEList")

# CPM to correct for sequencing depth and log2 transformed
cpm <- edgeR::cpm(myDGEList) 
colSums(cpm)
log2.cpm <- edgeR::cpm(myDGEList, log=TRUE)

# reshape unfiltered data from wide to long (tidy) format
log2.cpm.df <- as_tibble(log2.cpm, rownames = "geneID")
colnames(log2.cpm.df) <- c("geneID", sampleLabels) # add your sample names to this dataframe
log2.cpm.df.pivot <- pivot_longer(log2.cpm.df, # dataframe to be pivoted
                                  cols = -1, # column names to be stored as a SINGLE variable
                                  names_to = "samples", # name of that new variable (column)
                                  values_to = "expression") # name of new variable (column) storing all the values (data)

# plot 1: unfiltered, non-normalized
p1<-ggplot(log2.cpm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="unfiltered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

# filter low-expressed genes
table(rowSums(myDGEList$counts==0)==6) #tells you how many genes have zero counts in all 6 samples
cpm <- cpm(myDGEList) # counts per million (not logged) used for filtering
keepers <- rowSums(cpm>1) >= 3 # user-defined: keep genes with CPM > 1 in at least 3 samples; set 3 to the size of your smallest group
myDGEList.filtered <- myDGEList[keepers,] # subset the DGEList to the genes that passed

# reshape filtered data
log2.cpm.filtered <- cpm(myDGEList.filtered, log=TRUE)
log2.cpm.filtered.df <- as_tibble(log2.cpm.filtered, rownames = "geneID")
colnames(log2.cpm.filtered.df) <- c("geneID", sampleLabels)
log2.cpm.filtered.df.pivot <- pivot_longer(log2.cpm.filtered.df, # dataframe to be pivoted
                                           cols = -1, # column names to be stored as a SINGLE variable
                                           names_to = "samples", # name of that new variable (column)
                                           values_to = "expression") # name of new variable (column) storing all the values (data)

# normalize with TMM
myDGEList.filtered.norm <- calcNormFactors(myDGEList.filtered, method = "TMM")
log2.cpm.filtered.norm <- edgeR::cpm(myDGEList.filtered.norm, log=TRUE)
log2.cpm.filtered.norm.df <- as_tibble(log2.cpm.filtered.norm, rownames = "geneID")
colnames(log2.cpm.filtered.norm.df) <- c("geneID", sampleLabels)
log2.cpm.filtered.norm.df.pivot <- pivot_longer(log2.cpm.filtered.norm.df, # dataframe to be pivoted
                                                cols = -1, # column names to be stored as a SINGLE variable
                                                names_to = "samples", # name of that new variable (column)
                                                values_to = "expression") # name of new variable (column) storing all the values (data)

# plot 2: filtered, TMM normalized
p2<-ggplot(log2.cpm.filtered.norm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, TMM normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()
