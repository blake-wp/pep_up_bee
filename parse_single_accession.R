# Script for scraping 'Peptide Atlas' for tissue locations based on specific peptides or all peptides of a protein.
# 
# Rationale: can't seem to find on the website a way to export 'Experiment Peptide Map' data.
# Given a 'PeptideAtlas Build' and 'Protein Name', the underlying data is gathered from individual peptides
# in the 'Observed in Experiments' section.

library(tidyverse)
# Library to mask uncharacterised protein id.
library(dotenv)
# Library to scrape pages:
library(rvest)


#### ---- Get the list of peptides for the protein. ----

url_builder1 <- "https://db.systemsbiology.net/sbeams/cgi/PeptideAtlas/GetProtein?atlas_build_id=399&protein_name="
url_builder3 <- "&action=QUERY"

load_dot_env()
PA_protein_id <- Sys.getenv("UP_ACCESSION")

up_page <- read_html(paste0(url_builder1, PA_protein_id, url_builder3))
up_peptide_table <- data.frame(up_page %>% html_element(".PA_sort_table") %>% html_table())

pepatlas_pep_list <- 
  up_peptide_table %>% 
  select(Accession) %>% 
  unlist

#### ---- Get the data for each peptide. ----

url_builder4 <- "https://db.systemsbiology.net/sbeams/cgi/PeptideAtlas/GetPeptide?_tab=3&atlas_build_id=399&searchWithinThis=Peptide+Name&searchForThis="
url_builder6 <- "&action=QUERY"

pepatlas_data_list <- list()

# Loop over peptides, extract the table data as a separate list for each peptide.
for (i in 1:length(pepatlas_pep_list)) {
  # This loop can take a while if there are lots of peptides, so uncomment to print a message for fetching each peptide.
  #print(paste0("Fetching", i, " of ",  length(pepatlas_pep_list)))
  x <- read_html(paste(url_builder4, pepatlas_pep_list[i], url_builder6, sep = ""))
  y <- data.frame((x %>% html_elements(".PA_sort_table") %>% html_table())[[2]])
  y$Peptide <- rep(pepatlas_pep_list[i], nrow(y))
  pepatlas_data_list[i] <- list(y)
}

# Collapse the lists.
df <- do.call(rbind,pepatlas_data_list)

# Mask peptides with number.
df <-
  df %>% 
  mutate(pep_mask = as_factor(as.integer(as_factor(Peptide))))

#### ---- Plot. ----
ggplot(df, aes(pep_mask, Experiment.Tag)) +
  labs(x = "Peptide",
       y = "Peptide Atlas listed tissue",
       title = "Tissue localisation of peptides from the uncharacterised protein",
       subtitle = "Data collated from different experiments/instruments") +
  scale_fill_gradient(low = "grey60",
                      high = "green",
                      name = "Obs per\nmillion\nspectra") +
  geom_tile(aes(fill = NObs),
            colour = "black",
            width = 1,
            height = 1) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4, size = 8),
        axis.text.y = element_text(size = 8),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),
        panel.grid.major = element_blank(),
        plot.caption.position = "plot",
        plot.caption = element_text(hjust = 0),
        plot.title.position = "plot",
        panel.background = element_blank(),
        panel.border = element_rect(color = "grey60", fill = NA),
        legend.position = "right",
        strip.background = element_rect(color = "grey60", fill = "white"))
