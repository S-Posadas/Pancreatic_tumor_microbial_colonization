
# Taxonomy analysis #

#### Package setup ####

Sys.setenv(language = "EN")
library(phyloseq)
library(microViz)
library(ggplot2)
library(dplyr)

theme_set(theme_bw())

#### End ####

#### Load RData with filtered and normalized phyloseq object ####

res.dir <- "Results_check"
R.dir <- file.path(res.dir, "RData")
load(file.path(R.dir, "3.phyloseq.filtered.RData"))

# Create directories for results
tax.dir <- file.path(res.dir,"6.Taxonomy")
bar.dir <- file.path(tax.dir,"Bar_plots")

dir.create(bar.dir, recursive = T)
dir.create(stat.dir, recursive = T)

#### End ####

#### Check that the variables have the right class ####

# Change column values
sample_data(physeq) <- data.frame(sample_data(physeq)) %>%
  mutate(Sepsis = ifelse(Sepsis == 1, "Yes", ifelse(Sepsis == 2, "No", Sepsis)))
sample_data(physeq) <- data.frame(sample_data(physeq)) %>%
  mutate(Culture_growth = ifelse(Culture_growth == 1, "Microbial growth", ifelse(Culture_growth == 0, "No growth", Culture_growth)))
sample_data(physeq) <- data.frame(sample_data(physeq)) %>%
  mutate(Stent = ifelse(Stent == "Yes", "Stent", ifelse(Stent == "No", "No stent", Stent)))

sapply(sample_data(physeq), class)

# Transform all variables to factors in phyloseq object with absolute abundance
df <- as.data.frame(lapply(sample_data(physeq),function (y) if(class(y)!="factor" ) as.factor(y) else y),stringsAsFactors=T)
row.names(df) <- sample_names(physeq)
sample_data(physeq) <- sample_data(df)

# Transform all variables to factors in phyloseq object with relative abundance
df <- as.data.frame(lapply(sample_data(physeq_re),function (y) if(class(y)!="factor" ) as.factor(y) else y),stringsAsFactors=T)
row.names(df) <- sample_names(physeq_re)
sample_data(physeq_re) <- sample_data(df)

sapply(sample_data(physeq), class)

#### End ####

#### Change Species names for visualization ####

# Absolute abundance
tax <- as.data.frame(tax_table(physeq))
tax[,"Species"] = tax[,"Species.1"] ; tax[,"Species.1"] = NULL
tax[,"Species"] <- tax[,"Species"] %>% gsub("\\(RS_.*", "", .) %>% gsub("\\(GB_.*", "", .) %>% gsub("_RS_GCF_.*", "", .)
tax[!is.na(tax$Species),"Species"] = paste(tax[!is.na(tax$Species), "Genus"], tax[!is.na(tax$Species),"Species"])
tax <- as.matrix.data.frame(tax)
tax <- phyloseq::tax_table(tax)
phyloseq::tax_table(physeq) <- tax; rm(tax)

# Assign upper taxonomic ranks to unassigned taxa

#physeq <- tax_fix(physeq)

# Relative abundance

tax <- as.data.frame(tax_table(physeq_re))
tax[,"Species"] = tax[,"Species.1"] ; tax[,"Species.1"] = NULL
tax[,"Species"] <- tax[,"Species"] %>% gsub("\\(RS_.*", "", .) %>% gsub("\\(GB_.*", "", .) %>% gsub("_RS_GCF_.*", "", .)
tax[!is.na(tax$Species),"Species"] = paste(tax[!is.na(tax$Species), "Genus"], tax[!is.na(tax$Species),"Species"])
tax <- as.matrix.data.frame(tax)
tax <- phyloseq::tax_table(tax)
phyloseq::tax_table(physeq_re) <- tax; rm(tax)

# Assign upper taxonomic ranks to unassigned taxa

#physeq_re <- tax_fix(physeq_re)

#### End ####

#### Bar plots ####

Variable = "Stent"

# Absolute abundance based on Variable

plot_bar(physeq, y =  "Abundance", title = "Abundance based on Phylum", fill="Phylum")+ facet_wrap(as.formula(paste("~", Variable)), scales="free_x") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat = "Identity", position = "stack")
plot_bar(physeq, y =  "Abundance", title = "Abundance based on Class", fill="Class")+ facet_wrap(as.formula(paste("~", Variable)), scales="free_x") + 
  geom_bar(aes(color=Class, fill=Class), stat = "Identity", position = "stack")

# Relative abundance

plot_bar(physeq_re, y =  "Abundance", title = "Abundance based on Phylum", fill="Phylum")+ facet_wrap(as.formula(paste("~", Variable)), scales="free_x") +
  geom_bar(aes(color=Phylum, fill=Phylum), stat = "Identity", position = "stack")
plot_bar(physeq_re, y =  "Abundance", title = "Abundance based on Class", fill="Class")+ facet_wrap(as.formula(paste("~", Variable)), scales="free_x") +
  geom_bar(aes(color=Class, fill=Class), stat = "Identity", position = "stack")

#### End ####

#### Plot abundance from top n ASVs ####

Variable = "Stent"

n = 100 # adjust number to wished top ones
top <- names(sort(taxa_sums(physeq), decreasing=TRUE))[1:n]  

ps.top <- prune_taxa(top, physeq)
ps.top_re <- prune_taxa(top, physeq_re)

# Absolute abundance

plot_bar(ps.top, y =  "Abundance", title = paste("Abundance based on top", n, "ASV"), fill="Genus")+ facet_wrap(as.formula(paste("~", Variable)), scales="free_x") +
  geom_bar(aes(color=Genus, fill=Genus), stat = "Identity", position = "stack")

# Relative abundance

plot_bar(ps.top_re, y =  "Abundance", title = paste("Abundance based on top", n, "ASV"), fill="Genus")+ facet_wrap(as.formula(paste("~", Variable)), scales="free_x") + 
  geom_bar(aes(color=Genus, fill=Genus), stat = "Identity", position = "stack")

#### End ####

#### Plot abundance per condition ####

## Stent relative
ps <- merge_samples(physeq, "Stent")  # Merging 
sample_data(ps)$Stent <- rownames(sample_data(ps))
ps_per <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU)) #2nd transformation to make it again in percentage

plot_bar(ps_per, y =  "Abundance", fill="Phylum") + facet_wrap(~Stent, scales="free_x") +
  geom_bar(aes(color=Phylum, fill=Phylum), stat = "Identity", position = "stack") +
  ylab("Absolute abundance") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size = 20),
        legend.key.height = unit(0.03, "npc")) 

## Culture growth absolute
ps.growth <- merge_samples(physeq, "Culture_growth")  # Merging 
sample_data(ps.growth)$Culture_growth <- rownames(sample_data(ps.growth))

plot_bar(ps.growth, y =  "Abundance", fill="Phylum") + facet_wrap(~Culture_growth, scales="free_x") +
  geom_bar(aes(color=Phylum, fill=Phylum), stat = "Identity", position = "stack") +
  ylab("Absolute abundance") + xlab(NULL) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size = 20),
        legend.key.height = unit(0.03, "npc")) 

# Top n Phyla
n = 10
sum = tapply(taxa_sums(physeq), tax_table(physeq)[, "Phylum"], sum, na.rm=TRUE)
top = names(sort(sum, TRUE))[1:n]
topphy = prune_taxa((tax_table(ps.growth)[, "Phylum"] %in% top), ps.growth)

plot_bar(topphy, y =  "Abundance", fill="Phylum") + facet_wrap(~Culture_growth, scales="free_x") +
  geom_bar(aes(color=Phylum, fill=Phylum), stat = "Identity", position = "stack") +
  ylab("Absolute abundance") + xlab(NULL) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size = 20),
        legend.key.height = unit(0.03, "npc")) 

ggsave(file.path(bar.dir,"top.10.phyla.culture_growth.tiff"), width = 8, height = 6, dpi = 300)
ggsave(file.path(bar.dir,"top.10.phyla.culture_growth.svg"), width = 8, height = 6, dpi = 300)

# Top n Genus
n = 10
sum = tapply(taxa_sums(physeq), tax_table(physeq)[, "Genus"], sum, na.rm=TRUE)
top = names(sort(sum, TRUE))[1:n]
topgen = prune_taxa((tax_table(ps.growth)[, "Genus"] %in% top), ps.growth)

plot_bar(topgen, y =  "Abundance", fill="Genus") + facet_wrap(~Culture_growth, scales="free_x") +
  geom_bar(aes(color=Genus, fill=Genus), stat = "Identity", position = "stack") +
  ylab("Absolute abundance") + xlab(NULL) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size = 20),
        legend.key.height = unit(0.03, "npc")) 

ggsave(file.path(bar.dir,"top.10.genus.culture_growth.tiff"), width = 8, height = 6, dpi = 300)
ggsave(file.path(bar.dir,"top.10.genus.culture_growth.svg"), width = 8, height = 6, dpi = 300)

## Family of cultured bacteria absolute

# Normalize counts to number of samples per group

# 1.Calculate the number of samples in each group
sample_counts <- sample_data(physeq) %>%
  data.frame() %>%
  group_by(cultured_family) %>%
  summarise(sample_count = n())

# 2.Normalize the abundance values
otu_table_normalized <- as.data.frame(t(otu_table(physeq))) %>%  # Extract abundance data
  merge(sample_data(physeq), by = 0) %>%                         # Merge with sample counts
  left_join(sample_counts, by = "cultured_family") %>%
  mutate(across(starts_with("ASV"), ~ . / sample_count))         # Normalize by dividing by sample_count
rownames(otu_table_normalized) <- otu_table_normalized$Row.names

# 3.Select only counts with normalized data and convert back to matrix
otu_table_normalized_df <- otu_table_normalized[,grep("ASV", colnames(otu_table_normalized))] %>%
  t(.) %>% as.matrix(.)

# 4.Create new phyloseq object
ps_family_normalized <- phyloseq(otu_table(otu_table_normalized_df, taxa_are_rows = TRUE),
                                 tax_table(tax_table(physeq)),
                                 sample_data(sample_data(physeq)))

# 5.Transform all variables to factors
df <- as.data.frame(lapply(sample_data(ps_family_normalized),function (y) if(class(y)!="factor" ) as.factor(y) else y),stringsAsFactors=T)
row.names(df) <- sample_names(ps_family_normalized)
sample_data(ps_family_normalized) <- sample_data(df)

# Merge samples with the same cultured family
ps_cultured_family <- merge_samples(ps_family_normalized, "cultured_family")  # Merging 
sample_data(ps_cultured_family)$cultured_family <- rownames(sample_data(ps_cultured_family))

plot_bar(ps_cultured_family, y =  "Abundance", fill="Phylum") + facet_wrap(~cultured_family, scales="free_x") +
  geom_bar(aes(color=Phylum, fill=Phylum), stat = "Identity", position = "stack") +
  ylab("Absolute abundance") + xlab(NULL) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size = 20),
        legend.key.height = unit(0.03, "npc")) 

# Top n Phyla
n = 10
sum = tapply(taxa_sums(physeq), tax_table(physeq)[, "Phylum"], sum, na.rm=TRUE)
top = names(sort(sum, TRUE))[1:n]
topphy = prune_taxa((tax_table(ps_cultured_family)[, "Phylum"] %in% top), ps_cultured_family)

plot_bar(topphy, y =  "Abundance", fill="Phylum") + facet_wrap(~cultured_family, scales="free_x") +
  geom_bar(aes(color=Phylum, fill=Phylum), stat = "Identity", position = "stack") +
  ylab("Absolute abundance") + xlab(NULL) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size = 20),
        legend.key.height = unit(0.03, "npc")) 

ggsave(file.path(bar.dir,"top.10.phyla.culture_family_normalized.tiff"), width = 8, height = 10, dpi = 300)
ggsave(file.path(bar.dir,"top.10.phyla.culture_family_normalized.svg"), width = 8, height = 10, dpi = 300)

# Top n Genus
n = 10
sum = tapply(taxa_sums(physeq), tax_table(physeq)[, "Genus"], sum, na.rm=TRUE)
top = names(sort(sum, TRUE))[1:n]
topgen = prune_taxa((tax_table(ps_cultured_family)[, "Genus"] %in% top), ps_cultured_family)

plot_bar(topgen, y =  "Abundance", fill="Genus") + facet_wrap(~cultured_family, scales="free_x") +
  geom_bar(aes(color=Genus, fill=Genus), stat = "Identity", position = "stack") +
  ylab("Absolute abundance") + xlab(NULL) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size = 20),
        legend.key.height = unit(0.03, "npc")) 

ggsave(file.path(bar.dir,"top.10.genus.culture_family_normalized.tiff"), width = 8, height = 10, dpi = 300)
ggsave(file.path(bar.dir,"top.10.genus.culture_family_normalized.svg"), width = 8, height = 10, dpi = 300)

# Cultured genus
gen = c("Enterococcus", "Enterococcus_A", "Enterococcus_B", "Staphylococcus", "Escherichia", "Citrobacter", "Proteus", "Klebsiella", "Enterobacter")
selgen = prune_taxa((tax_table(ps_cultured_family)[, "Genus"] %in% gen), ps_cultured_family)

plot_bar(selgen, y =  "Abundance", fill="Genus") + facet_wrap(~cultured_family, scales="free_x") +
  geom_bar(aes(color=Genus, fill=Genus), stat = "Identity", position = "stack") +
  ylab("Absolute abundance") + xlab(NULL) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size = 20),
        legend.key.height = unit(0.03, "npc")) 

ggsave(file.path(bar.dir,"selected.genus.culture_family_normalized.tiff"), width = 8, height = 10, dpi = 300)
ggsave(file.path(bar.dir,"selected.genus.culture_family_normalized.svg"), width = 8, height = 10, dpi = 300)

#### End ####

#### Box plots ####

# Box plot function
plot_abundance = function(physeq, ylabn = "Absolute abundance",
                          Facet = "Class",
                          Color = "Stent",
                          n = NULL,
                          size = 20){
  mphyseq = psmelt(physeq)
  mphyseq <- subset(mphyseq, Abundance > 0)
  ggplot(data = mphyseq,
         mapping = aes_string(x = Color, y = "Abundance",
                              color = Color, fill = Color)) +
    theme(text = element_text(size = size)) +
    geom_boxplot(fill = NA) +
    geom_point() +
    facet_wrap(facets = Facet, nrow = n) + ylab(ylabn) +
    scale_y_log10()
}

# Agglomerate taxa according to desired taxa 
p.genus <- tax_glom(physeq, taxrank = "Genus")
p.phylum <- tax_glom(physeq, taxrank = "Phylum")

# Top n Phyla
n = 10
sum = tapply(taxa_sums(physeq), tax_table(physeq)[, "Phylum"], sum, na.rm=TRUE)
top = names(sort(sum, TRUE))[1:n]
topphy = prune_taxa((tax_table(p.phylum)[, "Phylum"] %in% top), p.phylum)

# Option 1
plot_abundance(topphy, Facet = "Phylum", Color = "Stent", n= 2) 

ggsave(file.path(bar.dir, "top.phyla.stent.tiff"), width = 12, height = 7, dpi = 300)
ggsave(file.path(bar.dir, "top.phyla.stent.svg"), width = 12, height = 7, dpi = 300)

# Option 2
plot_abundance(topphy, Facet = "Stent", Color = "Phylum", n= 1) + theme(axis.text.x.bottom = element_text(angle = -30, hjust = 0))

ggsave(file.path(bar.dir, "top.phyla.stent2.tiff"), width = 12, height = 7, dpi = 300)
ggsave(file.path(bar.dir, "top.phyla.stent2.svg"), width = 12, height = 7, dpi = 300)

# Top n Genus
n = 10
sum = tapply(taxa_sums(physeq), tax_table(physeq)[, "Genus"], sum, na.rm=TRUE)
top = names(sort(sum, TRUE))[1:n]
topgen = prune_taxa((tax_table(p.genus)[, "Genus"] %in% top), p.genus)

# Option 1
plot_abundance(topgen, Facet = "Genus", Color = "Stent", n= 2) 

ggsave(file.path(bar.dir, "top.genus.stent.tiff"), width = 13, height = 7, dpi = 300)
ggsave(file.path(bar.dir, "top.genus.stent.svg"), width = 13, height = 7, dpi = 300)

# Option 2
plot_abundance(topgen, Facet = "Stent", Color = "Genus", n= 1) + theme(axis.text.x.bottom = element_text(angle = -30, hjust = 0))

ggsave(file.path(bar.dir, "top.genus.stent2.tiff"), width = 12, height = 7, dpi = 300)
ggsave(file.path(bar.dir, "top.genus.stent2.svg"), width = 12, height = 7, dpi = 300)

#### End ####