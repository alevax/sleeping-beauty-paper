library(tidyverse)
getBeltranNEAnnotation <- function(){
  Beltran_metadata <- readRDS("experiments/oncomatch-analysis/processed_data/metadata/Beltran_metadata_forOMPlot.rds")
  Beltran_metadata_NE_feats <- Beltran_metadata$Neuroendocrine_Features
  names(Beltran_metadata_NE_feats) <- Beltran_metadata$Sample_Identifier
  return(Beltran_metadata_NE_feats)
}
getSU2CEastCoastNEAnnotation <- function(){
  SU2CEastCoast_metadata <- readRDS("experiments/oncomatch-analysis/processed_data/metadata/SU2C_EastCoast_metadata_forOMPlot.rds")
  SU2CEastCoast_metadata <- SU2CEastCoast_metadata %>%
    filter(!is.na(Neuroendocrine_Features)) #%>%
  SU2CEastCoast_metadata_NE_feats <- SU2CEastCoast_metadata$Neuroendocrine_Features
  names(SU2CEastCoast_metadata_NE_feats) <- SU2CEastCoast_metadata$Sample_Identifier
  return(SU2CEastCoast_metadata_NE_feats)
}

getSU2CWestCoastNEAnnotation <- function(){
  SU2CWestCoast_metadata <- readRDS("experiments/oncomatch-analysis/processed_data/metadata/SU2C_WestCoast_metadata_forOMPlot.rds")
  SU2CWestCoast_metadata <- SU2CWestCoast_metadata %>%
    filter(!is.na(Neuroendocrine_Features))# %>%
  SU2CWestCoast_metadata_NE_feats <- SU2CWestCoast_metadata$Neuroendocrine_Features
  names(SU2CWestCoast_metadata_NE_feats) <- SU2CWestCoast_metadata$Sample_Identifier
  return(SU2CWestCoast_metadata_NE_feats)
}


Beltran_om <- readRDS("experiments/oncomatch-analysis/processed_data/oncomatch/SleepingBeauty_Beltran_OncoMatch.rds")
Su2c_EC_om <- readRDS("experiments/oncomatch-analysis/processed_data/oncomatch/SleepingBeauty_SU2CEastCoast_OncoMatch.rds")
Su2c_WC_om <- readRDS("experiments/oncomatch-analysis/processed_data/oncomatch/SleepingBeauty_SU2CWestCoast_OncoMatch.rds")
gemm_metadata <- readRDS("experiments/oncomatch-analysis/processed_data/metadata/SleepingBeauty_metadata_forOMPlot.rds")


Beltran_metadata_NE_feats <- getBeltranNEAnnotation()
Su2c_EC_metadata_NE_feats <- getSU2CEastCoastNEAnnotation()
Su2c_WC_metadata_NE_feats <- getSU2CWestCoastNEAnnotation()
SB_metadata_NE_feats <- gemm_metadata$Phenotype

om_threshold <- 5

getPercMatches <- function(Human_om, SB_metadata_NE_feats, Human_metadata_NE_feats){
  round(sum(apply(Human_om[SB_metadata_NE_feats == "NEPC", Human_metadata_NE_feats == "Yes"] >= om_threshold, 2, any))/sum(Human_metadata_NE_feats == "Yes")*100,1)
}
SB_NE_Beltran_NE_perc_matches <- getPercMatches(Beltran_om, SB_metadata_NE_feats, Beltran_metadata_NE_feats)
SB_NE_Su2c_EC_NE_perc_matches <- getPercMatches(Su2c_EC_om, SB_metadata_NE_feats, Su2c_EC_metadata_NE_feats)
SB_NE_Su2c_WC_NE_perc_matches <- getPercMatches(Su2c_WC_om, SB_metadata_NE_feats, Su2c_WC_metadata_NE_feats)

plot_df <- data.frame(
  "Beltran" = SB_NE_Beltran_NE_perc_matches,
  "SU2C_East_Coast" = SB_NE_Su2c_EC_NE_perc_matches,
  "SU2C_West_Coast" = SB_NE_Su2c_WC_NE_perc_matches
) %>%
  pivot_longer(cols = c("Beltran", "SU2C_East_Coast", "SU2C_West_Coast"), names_to = "Dataset", values_to = "Percent_Matched") %>%
  mutate(perc_label = paste0(Percent_Matched, "%"))
bar_plot <- ggplot(plot_df, mapping = aes(x = Dataset, y = Percent_Matched, fill = Dataset, label = perc_label)) +
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(values = c("#E09F3E", "#003884", "#B71234")) +
  ylim(0, 100) + 
  theme_light() +
  # geom_text(hjust=+0.5, vjust=+5, size = 10, color = "white", fontface=2)
  geom_text(hjust=+0.5, vjust=-0.5, size = 2, color = "black", fontface=2) +
  theme(plot.title = element_text(size = 5),
        legend.title=element_text(size=5), 
        legend.text=element_text(size=5),
        axis.title = element_text(size=5),
        axis.text = element_text(size=5),
        axis.text.x = element_text(size=3))# +
  # theme(legend.position="none")
my_filename <- paste0("experiments/oncomatch-analysis/bar_plot.pdf")
ggsave(filename = my_filename,
       plot = bar_plot,
       device = pdf(),
       width = 3,
       height = 3,
       units = "in")
dev.off()

# "SU2C_East_Coast" = "#003884",
# "SU2C_West_Coast" = "#B71234",
# "Beltran" = "#E09F3E"

# Super small
# White background - theme light
# Add the number of matches to the plot


