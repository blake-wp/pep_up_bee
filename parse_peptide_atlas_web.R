library(rvest)
library(ggplot2)
library(tidyverse)
library(data.table)

gene_accesssions <- read.csv("gene-accession_list.txt", header = F)
gene_accesssions <- gene_accesssions[!(is.na(gene_accesssions$V2) | gene_accesssions$V2==""), ]

pepatlas_accessions <- gene_accesssions[,2][gene_accesssions[,2] != ""]


##### peptides may be similar across proteins
##### assign protein identifiers to the pepatlas_pep_list


url_builder1 <- "https://db.systemsbiology.net/sbeams/cgi/PeptideAtlas/GetProtein?atlas_build_id=399&protein_name="
url_builder3 <- "&action=QUERY"

pepatlas_pep_list <- list()

for (i in 1:length(pepatlas_accessions)) {
  print(paste(i, " of ",  length(pepatlas_accessions)))
  x <- read_html(paste(url_builder1, pepatlas_accessions[i], url_builder3, sep = ""))  
  y <- data.frame(x %>% html_element(".PA_sort_table") %>% html_table())
  y$pepAtlas_Acc <- rep(pepatlas_accessions[i], nrow(y))
  pepatlas_pep_list[i] <- list(y)
}
pepatlas_pep_list <- do.call(rbind,pepatlas_pep_list)
pepatlas_pep_list_check <- pepatlas_pep_list[pepatlas_pep_list$Accession %like% "PAp", ]

pepatlas_pep_list <- pepatlas_pep_list_check

url_builder4 <- "https://db.systemsbiology.net/sbeams/cgi/PeptideAtlas/GetPeptide?_tab=3&atlas_build_id=399&searchWithinThis=Peptide+Name&searchForThis="
url_builder6 <- "&action=QUERY"

pepatlas_data_list <- list()
for (i in 1:nrow(pepatlas_pep_list)) {
  print(paste(i, " of ",  nrow(pepatlas_pep_list)))
  x <- read_html(paste(url_builder4, pepatlas_pep_list[i,1], url_builder6, sep = ""))
  y <- data.frame((x %>% html_elements(".PA_sort_table") %>% html_table())[[2]])
  y$Peptide <- rep(pepatlas_pep_list[i,1], nrow(y))
  pepatlas_data_list[i] <- list(y)
}

pepatlas_data_list <- mapply(cbind, pepatlas_data_list, "pepAtlas_ACC"=pepatlas_pep_list$pepAtlas_Acc, SIMPLIFY=F)
df <- do.call(rbind,pepatlas_data_list)
df$Gene <- gene_accesssions$V1[match(df$pepAtlas_ACC, gene_accesssions$V2)]

gene_list <- gene_accesssions[,1]

# in df, sum NObs for each Experiment.Tag of each gene.
sum_obs <- df %>% group_by(across(c(Experiment.Tag, Gene))) %>% summarise(NObs = sum(NObs)) %>% ungroup
total_obs <- sum_obs %>% group_by(across(Gene)) %>% summarise(NObs = sum(NObs)) %>% ungroup
max_obs <- sum_obs %>% group_by(across(Gene)) %>% summarise(NObs = max(NObs)) %>% ungroup

#sum_obs$normalised <- apply(sum_obs, 1, function(i) { as.numeric(i['NObs']) / as.numeric(total_obs[match(i['Gene'],total_obs$Gene),2]) })

sum_obs$normalised <- apply(sum_obs, 1, function(i) { as.numeric(i['NObs']) / as.numeric(max_obs[match(i['Gene'],max_obs$Gene),2]) * 100 })

ggplot(sum_obs) +
  geom_bar(aes(NObs, Experiment.Tag, colour = Gene), stat = "identity", position = "dodge")

grid <- expand.grid(Experiment.Tag = unique(sum_obs$Experiment.Tag), Gene = unique(sum_obs$Gene))
grid$Gene <- factor(grid$Gene, levels = c("LOC408608", "A4", "*LOC100578744", "*LOC408961", "*LOC409051", "*LOC409278", "LOC410793", "LOC411955", "LOC552009", "LOC725217", "*LOC726446", "Vg"))
sum_obs$category <- ifelse(sum_obs$Gene == "LOC408608", "UP", "Apolipophorins\nor similar")
sum_obs$category <- factor(sum_obs$category, levels = c("UP", "Apolipophorins\nor similar"))
sum_obs$Gene <- factor(sum_obs$Gene, levels = c("LOC408608", "A4", "*LOC100578744", "*LOC408961", "*LOC409051", "*LOC409278", "LOC410793", "LOC411955", "LOC552009", "LOC725217", "*LOC726446", "Vg"))

ggplot(sum_obs, aes(Gene, Experiment.Tag)) +
  labs(y = "Peptide Atlas listed tissue") +
  scale_fill_gradient(low = "grey60", high = "green", name = "Normalised\nobservations") +
  #geom_tile(data = grid, color = "grey60", height = 1, width = 1, fill = NA) +
  geom_tile(aes(fill = normalised), colour = "black", width = 1, height = 1) +
  facet_grid(. ~ category, scales = "free", space = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4, size = 8),
        axis.text.y = element_text(size = 6),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),
        panel.grid.major = element_blank(),
        #axis.ticks.y = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "grey60", fill = NA),
        legend.position = "right",
        strip.background = element_rect(color = "grey60", fill = "white")
  )
  #scale_x_discrete(expand = c(0, 0.5))

#write.csv(unique(sum_obs$Experiment.Tag), "tissues.csv", quote = F)
