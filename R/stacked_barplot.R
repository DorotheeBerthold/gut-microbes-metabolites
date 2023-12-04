############################################################
## Metabolomic & 16S integration for correlation analysis ##
## Dorothée L. Berthold, ETH Zürich                       ##
############################################################

#use microbes_long_mean generated in 1_dataimport.R #######

#factor the categorical variables
microbes_long_mean$site <- factor(microbes_long_mean$site, levels = c("duo", "jej", "ile", "cec", "col"))

ggplot(microbes_long_mean, aes(y = mean_abundance, x = site, fill = Bacteria)) +
  geom_bar(position = "stack", stat = "identity") +
  theme_classic() +
  labs(x = "", y = "Relative abundance (%)") +
  ggtitle("Relative abundances in the lumen of 12w male mice") +
  facet_wrap(~origin, scales = "free", ncol = 2) +
  scale_fill_brewer(palette = "RdYlBu")
