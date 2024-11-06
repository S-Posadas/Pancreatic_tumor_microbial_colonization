
# Alpha diversity #

#### Package setup ####

Sys.setenv(language = "EN")
library(phyloseq)
library(ggplot2)
library(ggpubr)
library(openxlsx)
library(dplyr)
library(cowplot)
theme_set(theme_bw())

#### End ####

#### Load output from DADA2 pipeline with original phyloseq object ####

res.dir <- "Results_check"
R.dir <- file.path(res.dir, "RData")
load(file.path(R.dir, "2.physeq.original.RData"))

# Create directories for results
a.stats <- file.path(res.dir,"2.Alpha_stats")
a.plots <- file.path(res.dir,"3.Alpha_plots")

dir.create(a.stats, recursive = T)
dir.create(a.plots, recursive = T)

#### End ####

#### Check that the variables have the right class ####

# Change column values
sample_data(physeq) <- data.frame(sample_data(physeq)) %>%
  mutate(Sepsis = ifelse(Sepsis == 1, "Yes", ifelse(Sepsis == 2, "No", Sepsis)))
sample_data(physeq) <- data.frame(sample_data(physeq)) %>%
  mutate(Culture_growth = ifelse(Culture_growth == 1, "Yes", ifelse(Culture_growth == 0, "No", Culture_growth)))
sample_data(physeq)$ClavienDindo_simple = gsub("a", "", sample_data(physeq)$ClavienDindo) %>% gsub("b", "", .)

sapply(sample_data(physeq), class)

# Transform all variables to factors in phyloseq object 
df <- as.data.frame(lapply(sample_data(physeq),function (y) if(class(y)!="factor" ) as.factor(y) else y),stringsAsFactors=T)
row.names(df) <- sample_names(physeq)
sample_data(physeq) <- sample_data(df)

sapply(sample_data(physeq), class)                           

#### End ####

#### Create richness matrix ####
# Has to be counts not relative abundances

rich = estimate_richness(physeq, measures = c("Observed", "Shannon","InvSimpson"))

# Merge with sample_data for easy plotting

rich <- merge(rich, sample_data(physeq), by = 0)

# Save table with alpha diversity

write.xlsx(rich, file.path(a.stats, "richness.xlsx"))

#### End ####

#### Alpha diversity statistics####

#### Mann-Whitney-U-Test ####

mwut <- list()

for(v in c("Stent", "Culture_growth", "Sepsis")){
for (n in c("Observed", "Shannon", "InvSimpson")) {
  U_test <- wilcox.test(rich[!is.na(rich[,v]), n] ~ rich[!is.na(rich[,v]), v], data = rich, exact = T)
  z <- abs(qnorm(U_test$p.value/2))
  r <- z/sqrt(nrow(rich))
  
  tab <- c(U_test$method, n, v, U_test$statistic, U_test$p.value, r)
  mwut[[paste0(n, "_", v)]] <- tab
  
}}

mwut <- as.data.frame(do.call(rbind, mwut))
colnames(mwut) <- c("Test", "Variable1", "Variable2", "Statistic", "p-value", "Effect size")

mwut$p.adj <- p.adjust(mwut$`p-value`, method = "BH")

write.xlsx(mwut, file.path(a.stats, "Mann_Whitney_u_test.xlsx"), rowNames = T)

#### Kruskal-Wallis ####

kw <- list()

for (v in c("ClavienDindo", "ClavienDindo_simple")) {
  for (n in c("Observed", "Shannon", "InvSimpson")) {
      kruskal <- kruskal.test(rich[,n] ~ rich[,v])
      eta_squared <- (kruskal$statistic - length(unique(rich[,v])) + 1)/(nrow(rich) - length(unique(rich[,v])))
      f <- sqrt(eta_squared/(1-eta_squared))
      
      tab <- c(kruskal$method,  n, v, kruskal$statistic, kruskal$p.value, f)
      kw[[paste0(n, "_", v)]] <- tab
      
    }}

kw <- as.data.frame(do.call(rbind, kw))
colnames(kw) <- c("Test", "Variable1", "Variable2", "Chi2", "p-value", "Effect size")

write.xlsx(kw, file.path(a.stats, "Kruskal.wallis_time.xlsx"))

#### End ####

#### Alpha diversity plots ####

# Visualize relevant variables

Variable = "ClavienDindo_simple"
V_name = "Clavien-Dindo"
V_levels <- dim(unique(sample_data(physeq)[,Variable]))[1]

P = plot_richness(physeq, x=Variable, color = Variable, title = "Alpha Diversity", measures=c("Observed", "Shannon", "InvSimpson"))

P <- P+
  geom_boxplot(alpha = 0.5) +
  theme_bw() +
  labs(x = NULL, color = V_name) +
  theme(text = element_text(size=20),
        axis.text.y = element_text(size = 13),
        legend.position = "bottom") 

  if(V_levels == 2){
    P.1 <- P + stat_compare_means(method = "wilcox.test", paired = F, label.y = 1.8)
  }else{
    P.1 <- P + stat_compare_means(method = "kruskal.test", paired = F, label.y = 1.8)
  }

P.1  

# Save plot
ggsave(file.path(a.plots, paste(Variable, "alpha.tiff", sep = "_")), width = 10, height = 5.5, dpi = 300)
ggsave(file.path(a.plots, paste(Variable, "alpha.svg", sep = "_")), width = 10, height = 5.5, dpi = 300)


# Create plots of several variables of interest and combine

plots <- list()

for (Variable in c("Stent", "Culture_growth", "Sepsis", "ClavienDindo_simple")) {
  
  V_name = gsub("_", " ", Variable)
  if(Variable == "ClavienDindo_simple"){V_name = "Clavien-Dindo"}
  V_levels <- dim(unique(sample_data(physeq)[,Variable]))[1]
  
  P = plot_richness(physeq, x=Variable, color = Variable, measures=c("Observed", "Shannon", "InvSimpson"))
  
  P <- P+
    geom_boxplot(alpha = 0.5) +
    theme_bw() +
    labs(x = NULL, color = V_name) +
    theme(text = element_text(size=20),
          axis.text.y = element_text(size = 13),
          legend.position = "bottom") 
  
  if(V_levels == 2){
    P.1 <- P + stat_compare_means(method = "wilcox.test", paired = F, label.y = 1.8)
  }else{
    P.1 <- P + stat_compare_means(method = "kruskal.test", paired = F, label.y = 1.8)
  }
  
  plots[[Variable]] <- P.1
  
  # ggsave(file.path(a.plots, paste(Variable, "alpha.tiff", sep = "_")), width = 10, height = 5.5, dpi = 300)
  # ggsave(file.path(a.plots, paste(Variable, "alpha.svg", sep = "_")), width = 10, height = 5.5, dpi = 300)
  
}

# Without labels
plot_grid(plots[["Stent"]] , plots[["Culture_growth"]] , plots[["Sepsis"]] , plots[["ClavienDindo_simple"]] )

ggsave(file.path(a.plots, paste("alpha.tiff", sep = "_")), width = 20, height = 11, dpi = 300)
ggsave(file.path(a.plots, paste("alpha.svg", sep = "_")), width = 20, height = 11, dpi = 300)

# With labels
plot_grid(plots[["Stent"]] , plots[["Culture_growth"]] , plots[["Sepsis"]] , plots[["ClavienDindo_simple"]], labels = "AUTO", label_size = 20 )

ggsave(file.path(a.plots, paste("alpha_labs.tiff", sep = "_")), width = 20, height = 11, dpi = 300)
ggsave(file.path(a.plots, paste("alpha_labs.svg", sep = "_")), width = 20, height = 11, dpi = 300)

#### End ####
