# load packages----
library(rhdf5) #provides functions for handling hdf5 file formats (kallisto outputs bootstraps in this format)
library(tidyverse) # provides access to Hadley Wickham's collection of R packages for data science, which we will use throughout the course
library(tximport) # package for getting Kallisto results into R
library(ensembldb) #helps deal with ensembl
library(EnsDb.Hsapiens.v86) #replace with your organism-specific database package

# read in your study design using the readr package from tidyverse
studydesign <- read_tsv("studydesign.txt")

# Locate Kallisto abundance files
path <- file.path(studydesign$sample, "abundance.tsv")
# Now check to make sure this path is correct by seeing if the files exist
all(file.exists(path)) 

# Do transcript-to-gene annotations
Tx <- transcripts(EnsDb.Hsapiens.v86, columns=c("tx_id", "gene_name"))
Tx <- as_tibble(Tx)
# Need to change the first column name to 'target_id'
Tx <- dplyr::rename(Tx, target_id = tx_id)
# transcript ID needs to be the first column in the dataframe
Tx <- dplyr::select(Tx, "target_id", "gene_name")

# Import Kallisto results at gene level
Txi_gene <- tximport(path, 
                     type = "kallisto", 
                     tx2gene = Tx, 
                     txOut = FALSE, #How does the result change if this =FALSE vs =TRUE?
                     countsFromAbundance = "lengthScaledTPM",
                     ignoreTxVersion = TRUE)

# To look at the type of object you just created
class(Txi_gene)
names(Txi_gene)
