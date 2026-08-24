# =============================================================================
# Title: Multi-Omics for Drug Discovery: Finding Novel toxins from marine invertebrates
# Class: Multi-Omics project
# Research group: SeaTox Lab
# Author: Andrea Duque
# =============================================================================


# packages
library(dplyr)
library(UniprotR)


# =============================================================================
# PROBOSCIS TISSUE (P1 + P2 Replicates)
# TASK: Cleaning data
# =============================================================================


# Read BLASTp tabular output
p1_peptides_annotation <- read.csv("./01-Input/BLASTp_Proteomica_P1_Peptides_annotation.csv", sep = ";", header = FALSE)
p2_peptides_annotation <- read.csv("./01-Input/BLASTp_Proteomica_P2_Peptides_annotation.csv", sep = ";", header = FALSE)


# Standardize column names (BLAST tabular format)
blast_colnames <-c('ID','sacc','evalue','pident','qcovs','sstile')

colnames(p1_peptides_annotation) <- blast_colnames
colnames(p2_peptides_annotation) <- blast_colnames


# Merge technical replicates (all = TRUE)
P1P2 <-merge(x = p1_peptides_annotation, y = p2_peptides_annotation,
                  by = c('ID','sacc','evalue','pident','qcovs','sstile'),
                  all = TRUE)
                  

# Count occurrences per UniProt accession in merged data
# why? Multiple peptides may hit the SAME protein -> 9 hits ≠ 9 proteins
# Instead: 9 peptides matching 1 protein (higher confidence!).

occurences_P1P2 <- table(P1P2$sacc)
occurences_P1P2


# Aggregate by UniProt accession: mean metrics + peptide count

## mean metrics applied to: evalue, pident, qcovs,ID
proboscis <- aggregate(.~sacc, data = P1P2, FUN = function(x) if(is.numeric(x)) mean(x) else x[1])
proboscis <- proboscis[,-2] # remove ID column

## Add 'peptides' column with occurrence count per protein
## Why important? Single peptide + 100% coverage on long protein = low confidence
## Multiple peptides (e.g., 4+) = much higher identification reliability
proboscis$peptide <- occurences_P1P2 
View(proboscis)



# Export analysis-ready table
# write.csv(proboscis, './03-Output/proboscis_protein_hits.csv', row.names = FALSE)



# =============================================================================
# SKIN TISSUE (S1 + S2 Replicates)
# TASK: Cleaning data
# =============================================================================

# Read BLASTp tabular output
s1_peptides_annotation <- read.csv("./01-Input/BLASTp_Proteomica_S1_Peptides_annotation.csv", sep = ";", header = FALSE)
s2_peptides_annotation <- read.csv("./01-Input/BLASTp_Proteomica_S2_Peptides_annotation.csv", sep = ";", header = FALSE)

# Standardize column names (BLAST tabular format)
blast_colnames <-c('ID','sacc','evalue','pident','qcovs','sstile')
colnames(s1_peptides_annotation) <- blast_colnames
colnames(s2_peptides_annotation) <- blast_colnames

# Merge technical replicates (all = TRUE)
S1S2 <-merge(x = s1_peptides_annotation, y = s2_peptides_annotation,
             by = c('ID','sacc','evalue','pident','qcovs','sstile'),
             all = TRUE)

# Count occurrences per UniProt accession in merged data
occurences_S1S2 <- table(S1S2$sacc)
occurences_S1S2

# Aggregate by UniProt accession: mean metrics + peptide count
skin <- aggregate(.~sacc, data = S1S2, FUN = function(x) if(is.numeric(x)) mean(x) else x[1])
skin <- skin[,-2] 
skin$peptide <- occurences_S1S2 
View(skin)


# Export analysis-ready table
# write.csv(skin, './03-Output/skin_protein_hits.csv', row.names = FALSE)



# =============================================================================
# 
# =============================================================================


# Import analysis-ready data:

proboscis_protein <- read.csv('./03-Output/proboscis_protein_hits.csv', header = TRUE)
skin_protein <- read.csv('./03-Output/skin_protein_hits.csv', header = TRUE)



# Proboscis:

# Create a new column 'protein_name' by extracting name from ssitle
head(proboscis_protein$sstile,1) 
proboscis_protein <- proboscis_protein %>% 
  mutate(protein_name = gsub('OS.*', '', sstile))

# Q: Which proteins are of interest?
proboscis_protein$protein_name
View(proboscis_protein %>% select(protein_name, peptide))

# Q: Which protein has the highest number of peptides?
MaxPeptides_proboscis <- proboscis_protein %>%
  # Arrange by peptides in descending order and get the top protein
  arrange(desc(peptide)) %>%
  # Filter to keep all rows with the maximum peptide count
  filter(peptide == max(peptide)) %>%
  pull(protein_name)

cat('Protein with highest number of peptides:\n',
    paste(MaxPeptides_proboscis, collapse = ', '),
    '->',
    max(proboscis_protein$peptide), 'peptides')




## Skin:

# Create a new column 'protein_name'
skin_protein <- skin_protein %>% 
  mutate(protein_name = gsub('OS.*', '', sstile)) %>%
  arrange(desc(peptide))

# Q: Which proteins are of interest?
skin_protein$protein_name
View(skin_protein %>% select(protein_name, peptide))

# Q: Which protein has the highest number of peptides?
MaxPeptides_skin <- skin_protein %>%
  arrange(desc(peptide)) %>%
  filter(peptide == max(peptide)) %>%
  pull(protein_name)

cat('Protein with highest number of peptides:\n',
    paste(MaxPeptides_skin, collapse = ', '),
    '->',
    max(skin_protein$peptide), 'peptides')




## Proboscis vs. Skin

# Q: Are there common proteins between organs? How many and which ones?

# Find common proteins at protein name level
commonProteins <- intersect(proboscis_protein$protein_name,
                            skin_protein$protein_name)

cat('Common proteins between organs:\n',
    length(commonProteins), 'common proteins\n',
    paste(commonProteins, collapse = '\n'), '\n')


# Validation checks using original identifiers (optional)
commonSACC <- intersect(proboscis_protein$sacc,skin_protein$sacc)
commonSSTILES <- intersect(proboscis_protein$sstile,skin_protein$sstile)

cat('- Common protein name:',length(commonProteins),'\n',
    '- Common SACC IDs:', length(commonSACC),'\n',
    '- Common titles:', length(commonSSTILES))







## Analysis of Identified Proteins (UniProt data retrieval)

# Accession list for both organs
ProteinAccList_prob <-proboscis_protein$sacc
ProteinAccList_skin <-skin_protein$sacc





## Proboscis

# Q: What species were identified? (sample sequence matched to x,y,z specie sequence)

prob_NamesTaxa <- GetNamesTaxa(ProteinAccList_prob)
colnames(prob_NamesTaxa)
print(prob_NamesTaxa$Organism)

# Create organism summary (clean names + count)
proboscisOrganism <- prob_NamesTaxa %>% 
  mutate(Organism = gsub("\\s*\\(.*", "", Organism)) %>%
  count(Organism, sort = TRUE) %>%
  rename(Occurrences = n) 

View(proboscisOrganism)

cat('Number of species identified:\n',
    nrow(proboscisOrganism), 'which are: \n\n',
    paste(proboscisOrganism$Organism, collapse = '\n'), '\n')

           
# Q: What is the top species?
cat("\nTop species (Proboscis):\n", 
    "- ", proboscisOrganism$Organism[1], "[", proboscisOrganism$Occurrences[1], "]\n")







## Skin (same pattern)

# Q: What species were identified? 

skin_NamesTaxa <- GetNamesTaxa(ProteinAccList_skin)
print(skin_NamesTaxa$Organism)

# Create organism summary (clean names + count)
skinOrganism <- skin_NamesTaxa %>% 
  mutate(Organism = gsub("\\s*\\(.*", "", Organism)) %>%
  count(Organism, sort = TRUE) %>%
  rename(Occurrences = n) 

View(skinOrganism)

cat('Number of species identified:\n',
    nrow(skinOrganism), 'which are: \n\n',
    paste(skinOrganism$Organism, collapse = '\n'), '\n')


# Q: What is the top species?
cat("\nTop species (skin):\n", 
    "- ", skinOrganism$Organism[1], "[", skinOrganism$Occurrences[1], "]\n")





























