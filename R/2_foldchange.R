############################################################
## Metabolomic & 16S integration for correlation analysis ##
## Dorothée L. Berthold, ETH Zürich                       ##
############################################################


#use metabolites_long created in "dataimport.R"

#subset metabolites into GF and OMM only to calculate fold-changes
GF_OMM <- subset(metabolites_long, metabolites$colonization == c("GF", "sDMDM"))

#recode sites to match 16S data
GF_OMM <- GF_OMM %>%
  mutate(site2 = recode(site, "1" = "sto", "2" = "duo", "3" = "duo", 
                        "4" = "duo", "5" = "jej", "6" = "jej", 
                        "7" = "jej", "8" = "jej", "9" = "ile", 
                        "10" = "ile", "11" = "ile", "12" = "cec", 
                        "13" = "col", "14" = "col", "15" = "col"))

#get rid of stomach
keep <- c("duo","jej", "ile", "cec", "col")
GF_OMM <- GF_OMM[GF_OMM$site2 %in% keep,]

########################## Filter data for FC calculations ################################################################

#check how many observations per site etc for filtering
#here we group origins together in order to get a more thorough understanding of metabolite abundances
counts_metabolites <- GF_OMM %>%
  group_by(site2, metabolite, colonization) %>%
  summarise(count = n(),
            count_above_0 = sum(intensity > 0 & !is.na(intensity)))

#filter out metabolites that have less than 5 counts over 0 per site/colonization
counts_metabolites_low <- counts_metabolites |> 
  filter(count_above_0 < 5)

sum(counts_metabolites_low$count) #1342 observations should be filtered out --> 18882 remain

#create barcodes for filtering
counts_metabolites_low$barcode <- paste0(counts_metabolites_low$site2, "_", counts_metabolites_low$metabolite, "_", counts_metabolites_low$colonization)
GF_OMM$barcode <- paste0(GF_OMM$site2, "_", GF_OMM$metabolite, "_", GF_OMM$colonization)

#now onto filtering
# Find the set difference between the barcodes in GF_OMM and low metabolite counts
barcodes_to_keep <- setdiff(GF_OMM$barcode, counts_metabolites_low$barcode)

# Create a new dataframe by filtering GF_OMM based on the barcodes_to_keep
GF_OMM_filtered <- GF_OMM[GF_OMM$barcode %in% barcodes_to_keep, ]


########################### Normality check + ttest calculation ########################################

#check for normalcy with Q-Q plot
#here we use testdata to randomly sample 3000 from our df
testdata <- dplyr::sample_n(GF_OMM_filtered, 3000)
ggqqplot(testdata$intensity)
shapiro.test(testdata$intensity)
#data is not normally distributed --> nonparametric test

#try to log transform Abundances
testdata <- testdata %>%  
  mutate(logInt = log(intensity))
ggqqplot(testdata$logInt)
shapiro.test(testdata$logInt)
#log Abundance is normally distributed


#logtransform (log2) GF_OMM_filtered with pseudo +1 for all 0's
GF_OMM_filtered <- GF_OMM_filtered %>% 
  mutate(logInt = log(intensity +1))


#ttest on log values: parametric, unpaired, BH correction by metabolite & site
ttest <- compare_means(logInt ~ colonization, data = GF_OMM_filtered, paired = F, 
                       method = "t.test", p.adjust.method = "BH", 
                       group.by = c("metabolite", "site2"), na.rm = T)


#subset dataframe based on different microbiota colonization
GF_long <- subset(GF_OMM_filtered, GF_OMM_filtered$colonization == "GF")
OMM_long <- subset(GF_OMM_filtered, GF_OMM_filtered$colonization == "sDMDM")


#calculate average per metabolite per site recoded
OMM_long <- OMM_long %>% 
  group_by(metabolite, site2) %>% 
  mutate(mean = mean(logInt, na.rm = T))

GF_long <- GF_long %>% 
  group_by(metabolite, site2) %>% 
  mutate(mean = mean(logInt, na.rm = T))

#keep distinct values per site
OMM_long <- OMM_long %>% 
  group_by(mean) %>% 
  distinct(metabolite, site2)

GF_long <- GF_long %>% 
  group_by(mean) %>% 
  distinct(metabolite, site2)

colnames(OMM_long)[1] = "log_OMM_concentration"
colnames(GF_long)[1] = "log_GF_concentration"

#join the two datatables
GF_OMM_long <- left_join(GF_long, OMM_long, by = c("metabolite", "site2"), keep = F)

#for log transformed values (foldchange calculation):
GF_OMM_long <- GF_OMM_long %>% mutate(log2FC = (log_OMM_concentration-log_GF_concentration))

#add pvals
GF_OMM_complete <- left_join(GF_OMM_long, ttest, by = c("metabolite", "site2"))

#log10 for padj
GF_OMM_complete <- GF_OMM_complete %>%  
  mutate(log10_padj = -log(p.adj))

GF_OMM_complete <- as.data.frame(GF_OMM_complete)

############################# Plot & Labelling #####################################################

# Create scatter plot for correlations
ggplot(data=GF_OMM_complete, aes(x=log2FC, y=log10_padj)) + geom_point() +facet_wrap(~site2)

### based on log2FC 0.6 & p-value < 0.05, set diffexpressed as up or down-regulated
# add a column of NAs
GF_OMM_complete$diffexpressed <- "NO"
# if log2Foldchange > 1 and pvalue < 0.05, set as "UP" 
GF_OMM_complete$diffexpressed[GF_OMM_complete$log2FC > 1 & GF_OMM_complete$p.adj < 0.05] <- "UP"
# if log2Foldchange < -1 and pvalue < 0.05, set as "DOWN"
GF_OMM_complete$diffexpressed[GF_OMM_complete$log2FC < -1 & GF_OMM_complete$p.adj < 0.05] <- "DOWN"


#label correspondingly with the metabolites matching the up- and downregulated metabolites
GF_OMM_complete$delabel <- NA
GF_OMM_complete$delabel[GF_OMM_complete$diffexpressed != "NO"] <- as.character(GF_OMM_complete$metabolite[GF_OMM_complete$diffexpressed != "NO"])

# Re-plot but this time color the points with "diffexpressed" & textlabels

ggplot(data=GF_OMM_complete, aes(x=log2FC, y=-log10(p.adj), col=diffexpressed, label=delabel)) + 
  geom_point() + 
  theme_classic() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-1, 1), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") + 
  labs(title = "Fold change OMM vs GF mice", x = "log2 Fold change", y= "-log10 adj p-value") +
  facet_wrap(~site2, scales = "free")
