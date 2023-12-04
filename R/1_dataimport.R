############################################################
## Metabolomic & 16S integration for correlation analysis ##
## Dorothée L. Berthold, ETH Zürich                       ##
############################################################

#Import required packages and functions
source("gutPackages.R")
gutPackages()

#Import metabolite intensities and convert to long format
##############################################################################################################################
#metabolite_conc: df with columns as different metabolites & row different samples (identifier in column "barcode")
#meta_metabolties: df with meta_data containing "barcode" column to match and join with metabolite_conc

metabolite_conc <- read.csv("examples/intdata_all.csv")
meta_metabolite <- read.csv("examples/meta_metabolites.csv")

#generates new joined df "metabolites" containing both meta-data and measured data if rownames match

if (all(metabolite_conc$barcode == meta_metabolite$barcode)) {
  metabolites <- left_join(meta_metabolite, metabolite_conc, by = "barcode")
} else {
  print("Barcode columns don't match")
}

#convert "metabolites" into long format
metabolites <- metabolites[,-c(1,6)]
metabolites_long <- reshape2::melt(metabolites, id.vars = c("mouse.number", "site", "colonization", "origin"))
colnames(metabolites_long)[5:6] <- c("metabolite", "intensity")

#metabolites dataset have 15 sample sites whereas 16S only has five larger ones (without stomach)
#1: stomach
#2-4: duodenum
#5-8: jejunum
#9-11: ileum
#12: cecum
#13-15: colon
####################. Recode metabolites_long to match 16S sampling sites
metabolite_long2 <- metabolites_long %>% 
  mutate(site2 = recode(site, "1" = "sto", "2" = "duo", "3" = "duo", "4" = "duo", "5" = "jej", "6" = "jej", "7" = "jej", "8" = "jej", 
                        "9" = "ile", "10" = "ile", "11" = "ile", "12" = "cec", "13" = "col", "14" = "col", "15" = "col"))

#calculate mean per site recoded
metabolite_long2 <- metabolite_long2 %>% 
  group_by(site2, metabolite, origin, colonization) %>% 
  mutate(mean_intensity = mean(intensity, na.rm = T)) %>% 
  ungroup()

#select unique values for metabolites per site recoded, keeping the mean
metabolite_long2 <- metabolite_long2 %>% 
  group_by(mean_intensity) %>% 
  distinct(metabolite, site2, origin, colonization) %>% 
  ungroup()

#remove stomach as not included in 16S samples
keep <- c("duo","jej", "ile", "cec", "col")
metabolite_long2 <- metabolite_long2[metabolite_long2$site2 %in% keep,]


#Import 16S abundances and convert to long format
################################################################################################################################
#read in 16S and corresponding meta data
#family_16S: fraction table containing bacterial families in a column "rank" and all samples in columns
#meta_16S: long-table format with a column "barcode" matching column names of family_16S

#read in 16S sequencing data
family_16S <- as.data.frame(read_excel("examples/rarefied_OTU_table_RA_Family.xlsx"))

#read in metadata
meta_16S <- read.csv("examples/metadata_16S.csv")
meta_16S <- subset(meta_16S, microbiota == "OMM")

#multiply fractions * 100 in order to create relative abundances instead of fractions
#here column 1 are the bacteria, so they excluded from the calculations
family_16S[, 2:64] <- family_16S[, 2:64] * 100

#make bacteria names into rownames
rownames(family_16S) <- family_16S$rank
family_16S <- family_16S[,2:64]

#map to metadata
family_16S <- family_16S[,names(family_16S) %in% meta_16S$barcode]

#define bacteria families to keep for downstream analysis
OMM_family <- c("Bacteroidaceae", "Muribaculaceae", "Akkermansiaceae", "Bifidobacteriaceae", "Sutterellaceae", "Oscillospiraceae", 
                "Ruminococcaceae", "Lachnospiraceae", "Erysipelotrichaceae", "Enterococcaceae", "Lactobacillaceae")

#just keep OMM members in family plot
family_16S <- family_16S[rownames(family_16S) %in% OMM_family,]

#transpose dataframe to add metadata
family_16S <- as.data.frame(t(family_16S))
family_16S$barcode <- rownames(family_16S)

#join with metadata by barcode
microbes <- left_join(meta_16S, family_16S, by = "barcode")
microbes <- microbes[,-c(1,4)]

#generate long format
microbes_long <- reshape2::melt(microbes, id.vars = c("origin", "site", "mouse number"))
colnames(microbes_long) [4:5] <- c("Bacteria", "Abundance")

#calculate mean abundance across biological replicates
microbes_long_mean <- microbes_long  %>%  
  group_by(Bacteria, site, origin) %>% 
  mutate(mean_abundance = mean(Abundance, na.rm = T)) %>% 
  ungroup() %>%  
  group_by(mean_abundance) %>%  
  distinct(Bacteria, site, origin)

microbes_long_mean <- as.data.frame(microbes_long_mean)
