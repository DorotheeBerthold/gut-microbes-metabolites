



#read in metabolomic data
OMM_conc <- read.csv("tables/intdata_all.csv")
colnames(OMM_conc)[1] <- "barcode"

#extract metadata from barcode
OMM_conc$string <- str_right(OMM_conc$barcode, 6)
OMM_conc$mouse.number <- str_left(OMM_conc$string, 1)
OMM_conc$origin <- str_mid(OMM_conc$string, 3, 1)
OMM_conc$site <- str_mid(OMM_conc$string, 5, 6)
OMM_conc$colonization <- gsub(".*_(GF|sDMDM|SPF)_.*", "\\1", OMM_conc$barcode)

#convert into long format
int_long <- reshape2::melt(OMM_conc, id.vars = c("barcode", "mouse.number", "site", "colonization", "origin"))
int_long <- int_long |> 
  filter(variable != "string")

#write.csv(int_long, "results/int_all_long.csv")

#subset for only OMM
OMM_long <- subset(int_long, int_long$colonization == "sDMDM")

#write.csv(OMM_long, "results/OMM_intensity_long.csv")
