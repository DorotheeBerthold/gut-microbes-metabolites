#calculate mean abundance across biological replicates
microbes_long_mean <- microbes_long  %>%  
  group_by(Bacteria, site, origin) %>% 
  mutate(mean_abundance = mean(Abundance, na.rm = T)) %>% 
  ungroup() %>%  
  group_by(mean_abundance) %>%  
  distinct(Bacteria, site, origin)

microbes_long_mean <- as.data.frame(microbes_long_mean)

#factor the categorical variables
microbes_long_mean$site <- factor(microbes_long_mean$site, levels = c("duo", "jej", "ile", "cec", "col"))

ggplot(microbes_long_mean, aes(y = mean_abundance, x = site, fill = Bacteria)) +
  geom_bar(position = "stack", stat = "identity") +
  theme_classic() +
  labs(x = "", y = "Relative abundance (%)") +
  ggtitle("Relative abundances in the lumen of 12w male mice") +
  facet_wrap(~origin, scales = "free", ncol = 2) +
  scale_fill_brewer(palette = "RdYlBu")
