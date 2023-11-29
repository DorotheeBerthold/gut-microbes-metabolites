############################################################
## Metabolomic & 16S integration for correlation analysis ##
## Dorothée L. Berthold, ETH Zürich                       ##
############################################################



#read in metabolomic data
metabolite_conc <- read.csv("tables/intdata_all.csv")
colnames(metabolite_conc)[1] <- "barcode"

#extract metadata from barcode
metabolite_conc$string <- str_right(metabolite_conc$barcode, 6)
metabolite$mouse.number <- str_left(metabolite$string, 1)
metabolite$origin <- str_mid(metabolite$string, 3, 1)
metabolite$site <- str_mid(metabolite$string, 5, 6)
metabolite$colonization <- gsub(".*_(GF|sDMDM|SPF)_.*", "\\1", metabolite$barcode)

#convert into long format
int_long <- reshape2::melt(metabolite_conc, id.vars = c("barcode", "mouse.number", "site", "colonization", "origin"))
int_long <- int_long |> 
  filter(variable != "string")
