############################################################
## Metabolomic & 16S integration for correlation analysis ##
## Dorothée L. Berthold, ETH Zürich                       ##
############################################################

#################### Create wide-formats for correlation calculations #################################

#use wide formats, using sites as replicates for correlations
#here I used content + lumen to find more correlations
#use OMM_long created in "2_foldchange.R" - filtered and subset for colonization "sDMDM2"
#use microbes_long_mean created in "1_dataimport.R"

metabolite_input <- reshape2::acast(OMM_long, formula = metabolite ~ site2, 
                                        fill = 0, value.var = "log_OMM_concentration", 
                                        fun.aggregate = sum, drop = FALSE) 

metabolite_input <- apply(metabolite_input, 2, exp)
metabolite_input <- as.data.frame(metabolite_input)

OMM_input <- reshape2::acast(microbes_long_mean, formula = Bacteria ~ site, 
                                fill = 0, value.var = "mean_abundance", 
                                fun.aggregate = sum, drop = FALSE) 
OMM_input <- as.data.frame(OMM_input)


##################### Calculate correlations ##########################################################

# Create empty data frames for correlations and p-values
correlations <- data.frame(matrix(NA, nrow = nrow(metabolite_input), ncol = nrow(OMM_input)))
colnames(correlations) <- rownames(OMM_input)
rownames(correlations) <- rownames(metabolite_input)
pvals <- data.frame(matrix(NA, nrow = nrow(metabolite_input), ncol = nrow(OMM_input)))
colnames(pvals) <- rownames(OMM_input)
rownames(pvals) <- rownames(metabolite_input)


for (seq in rownames(OMM_input)) {
  seqLvl <- as.numeric(OMM_input[seq,])
  
  for (met in rownames(metabolite_input)) {
    metLvl <- as.numeric(metabolite_input[met,])
    
    # Exclude NA values from correlation calculation
    validIndices <- complete.cases(seqLvl, metLvl)
    validSeqLvl <- seqLvl[validIndices]
    validMetLvl <- metLvl[validIndices]
    
    if (length(validSeqLvl) >= 2 && length(validMetLvl) >= 2) {
      results <- cor.test(validSeqLvl, validMetLvl, method = "spearman", na.rm = TRUE, exact = F)
      correlations[met, seq] <- as.numeric(results$estimate)
      pvals[met, seq] <- as.numeric(results$p.value)
    } else {
      # Insufficient non-NA values for correlation calculation
      correlations[met, seq] <- NA
      pvals[met, seq] <- NA
    }
  }
}


#add metabolite as column for reshaping to long format
correlations$metabolite <- rownames(correlations)
pvals$metabolite <- rownames(pvals)


# Convert correlations dataframe to long format and rename columns
correlations_long <- reshape2::melt(correlations, id.vars = "metabolite", value.name = "correlation")
colnames(correlations_long) <- c("metabolite", "Bacteria", "correlation")

# Convert pvals dataframe to long format and rename column
pvals_long <- reshape2::melt(pvals, id.vars = "metabolite", value.name = "pval")
colnames(pvals_long) <- c("metabolite", "Bacteria", "pval")

# Concatenate dataframes along columns
corrDat <- cbind(correlations_long, pvals_long)
corrDat <- corrDat[,-c(4,5)]
corrDat <- corrDat %>% 
  filter(!is.na(correlation))

################# Filter correlation matrix ##########################################################

# Filter out metabolites from GF_OMM_complete df that are upregulated in GF (DOWN in OMM vs GF)
# GF_OMM_complete created in "foldchange.R"

GF_down <- GF_OMM_complete |> 
  filter(diffexpressed == "DOWN")

#create vector for GF metabolites
GF_metabolites <- unique(as.character(GF_down$metabolite))

#filter corrDat for those metabolites
filtered_corrDat <- corrDat[!(corrDat$metabolite) %in% GF_metabolites,]



#set treshold for correlations (based on Meier et al., 2023)
thres <- 0.76

# Plot histogram
hist(corrDat$correlation, breaks = 20, xlim = c(-1, 1), ylim = c(0, 200),
     col = "darkgrey", border = "black", lwd = 0.8, xlab = "Correlation coefficient", ylab = "Number of metabolite-microbe pairs", main = "Microbe-metabolite correlation")

# Add vertical lines
abline(v = 0, col = "black", lwd = 2, lty = 1)
abline(v = thres, col = "black", lwd = 1.2, lty = "dotted")
abline(v = -thres, col = "black", lwd = 1.2, lty = "dotted")

#filter correlations higher than thresh (+&-) & p-value
filtered_corrDat <- filtered_corrDat[abs(filtered_corrDat$correlation) >= thres, ] #307 higher correlations total 
filtered_corrDat <- filtered_corrDat[abs(filtered_corrDat$pval) < 0.05, ] #223 significant correlations

################################# Heatmap for correlation matrix #################################################################

#make filtered_corrDat into matrix for heatmap plotting
corrDat_wide <- reshape2::acast(filtered_corrDat, 
                                formula = metabolite ~ Bacteria, 
                                fill = 0, 
                                value.var = "correlation", 
                                fun.aggregate = sum, 
                                drop = FALSE)
#heatmap plotting
pheatmap(corrDat_wide, 
         cluster_rows = T, 
         cluster_cols = F, 
         border_color = "black",
         color = hcl.colors(50, "BluYl"))
