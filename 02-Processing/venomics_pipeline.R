# =============================================================================
# BIOINFORMATIC PIPELINE: Glycera alba Multi-Omics Venomics Project
# Input: BLASTp tabular output (Swiss-prot, e-value 1e-5, top hit per peptide)
# =============================================================================


# packages
library(dplyr)



# =============================================================================
# PROBOSCIS: Process proboscis replicates (P1+P2)
# =============================================================================


# Read BLASTp tabular output
p1_peptides_annotation <- read.csv("./01-Input/BLASTp_Proteomica_P1_Peptides_annotation.csv", sep = ";", header = FALSE)
p2_peptides_annotation <- read.csv("./01-Input/BLASTp_Proteomica_P2_Peptides_annotation.csv", sep = ";", header = FALSE)


# Standardize column names (BLAST tbular format)
blast_colnames <-c('ID','sacc','evalue','pident','qcovs','sstile')

colnames(p1_peptides_annotation) <- blast_colnames
colnames(p2_peptides_annotation) <- blast_colnames


# Merge technical replicares (all = TRUE)
proboscis <-merge(x = p1_peptides_annotation, y = p2_peptides_annotation,
                  by = c('ID','sacc','evalue','pident','qcovs','sstile'),
                  all = TRUE)
                  



