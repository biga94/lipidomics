#Caricamento librerie
library(readxl)
library(lipidr)
library(tidyverse)
library(FactoMineR)
library(factoextra)

#####
#Upload dati + aggiustamenti nomenclatura + nomenclatura
data <- read_excel("22QT37 Results_lipidomica rielab_26.04.2023.xlsx", 
                                                           sheet = "Results TOT (2)")
data$Ontology <- gsub("Ether DG","EtherDG", data$Ontology)
classi <- read_excel("22QT37 Results_lipidomica rielab_26.04.2023.xlsx")

swiss_lipids <- read_tsv(file = "lipids.tsv", col_names = T)

background <- merge(data,swiss_lipids,by.x = "Formula bruta", by.y = "Formula (pH7.3)")
background <- background |> filter(Level == "Species")
background$`Metabolite name` <- gsub("\\|.*", "", background$`Metabolite name`) #rimuovo tutto il supplementary
background$`Metabolite name` <-  sub(";(O\\d*)", "(\\1)", background$`Metabolite name`)
background$`Metabolite name` <-  sub(" O-", "O ", background$`Metabolite name`)
background$`Metabolite name` <-  sub(" P-", "P ", background$`Metabolite name`)
background$`Metabolite name` <-  sub("-FA", "/", background$`Metabolite name`)
background$`Metabolite name` <- gsub("\\|.*", "", background$`Metabolite name`) #rimuovo tutto il supplementary

data$`Metabolite name` <- gsub("\\|.*", "", data$`Metabolite name`) #rimuovo tutto il supplementary
data$`Metabolite name` <-  sub(";(O\\d*)", "(\\1)", data$`Metabolite name`)
data$`Metabolite name` <-  sub(" O-", "O ", data$`Metabolite name`)
data$`Metabolite name` <-  sub(" P-", "P ", data$`Metabolite name`)
data$`Metabolite name` <-  sub("-FA", "/", data$`Metabolite name`)
data$`Metabolite name` <- gsub("\\|.*", "", data$`Metabolite name`)
data[,9:24] <- log(data[,9:24]) #log transformation


#Creazione dataset informazione clinica
annotations <- names(data)[9:24]
data_clin <- data.frame(Sample = annotations)
data_clin$Group <- sub(".*_.*_(.*)", "\\1", data_clin$Sample)
data_clin$Group <- factor(data_clin$Group)
Dati_anonimizzati_OB_Lipidomica <- read_excel("Dati anonimizzati_OB_Lipidomica.xlsx")
data_clin <- cbind(data_clin, Dati_anonimizzati_OB_Lipidomica[,3:6])
data_clin$Età <- round(data_clin$Età)
names(data_clin)[6] <- "BMI"
rm(annotations)

#####
#Analisi multivariata
data_pca <- data |> select(-c(1:6))
pca0 <- PCA(data_pca, scale.unit = T, graph = T, quali.sup = 1:2, )
fviz_eig(pca0)

fviz_pca_biplot(pca0, geom = "point", col.var = "black", addEllipses = F,
                col.ind = "darkgoldenrod1",alpha.ind = 0.4, label = c("ind", "ind.sup", "var"),
                repel = T, title = "", axes = c(1,2))


#Transposta

rownames(data_new) <- data_new$`Metabolite name`
data_new <- t(data_new)
data_new <- data_new[-c(1:2),]

data_new <- as.data.frame(data_new)
data_new[] <- lapply(data_new, as.numeric)
data_new <- cbind(data_new, data_clin$Group)

pca1 <- PCA(data_new, scale.unit = T, graph = T, quali.sup = 823)



fviz_eig(pca1)

fviz_pca_biplot(pca1, axes = c(1,2), geom = "point",
                col.var = "black" , addEllipses = F,
                habillage = data_clin$Group,
                alpha.ind = 1, alpha.var = 0.05, label = c("ind"),
                repel = T, title = "PCs 1 and 2", geom.ind = c("text", "point"))

fviz_pca_biplot(pca1, axes = c(2,3), geom = "point",
                col.var = "black", addEllipses = F,
                habillage = data_clin$Group,
                alpha.ind = 1, alpha.var = 0.05, label = c("ind"),
                repel = T, title = "PCs 2 and 3", geom.ind = c("text", "point"))


#####
#Analisi univariata con LipidR
data_lipidr <- data |> select(-c(1:6,8))

d <- as_lipidomics_experiment(data_lipidr, logged = T, normalized = T)
d <- add_sample_annotation(d, data_clin)
rowData(d)$Class <- data$Ontology

plot_lipidclass(d, "boxplot") #distribuzione dei valori di area per ogni classe di lipidi

OWvsNW <- de_analysis(d, OW-NW)
OBvsNW <- de_analysis(d, OB-NW)
SVvsNW <- de_analysis(d, SV-NW)

tmp <- de_analysis(d, OB-NW,OW-NW, SV-NW)
plot_results_volcano(tmp)
tmp2 <- merge(tmp,classi,by.x="Class",by.y="Abbreviation")
tmp3 <- tmp2 |> select(-`Lipid subclass`,-Class) |> relocate(Categories, .after = Molecule) |> 
  rename(Class = Categories)

plot_results_volcano(tmp3)
#OW vs NW
plot_results_volcano(OWvsNW)
OWvsNW_merged <- merge(OWvsNW,background,
                       by.x="Molecule", by.y = "Metabolite name")
OWvsNW_merged <- OWvsNW_merged |> 
  filter(logFC > 1 & adj.P.Val <= 0.05) |> 
  select(`Abbreviation*`, adj.P.Val)

OWvsNW_merged <- merge(OWvsNW,background,
                       by.x="Molecule", by.y = "Metabolite name") #devo ricaricare il precedente dataset
OWvsNW_merged_down <- OWvsNW_merged |> 
  filter(logFC < -1 & adj.P.Val <= 0.05) |> 
  select(`Abbreviation*`, adj.P.Val)

#OB vs NW
plot_results_volcano(OBvsNW)
OBvsNW_merged <- merge(OBvsNW, background, 
                       by.x="Molecule", by.y="Metabolite name")
OBvsNW_merged <- OBvsNW_merged |> 
  filter(logFC > 1 & adj.P.Val <= 0.05) |> 
  select(`Abbreviation*`, adj.P.Val)

OBvsNW_merged <- merge(OBvsNW, background, 
                       by.x="Molecule", by.y="Metabolite name")
OBvsNW_merged_down <- OBvsNW_merged |> 
  filter(logFC < -1 & adj.P.Val <= 0.05) |> 
  select(`Abbreviation*`, adj.P.Val)

#SV vs NW
plot_results_volcano(SVvsNW)
SVvsNW_merged <- merge(SVvsNW, background,
                       by.x="Molecule", by.y="Metabolite name")
SVvsNW_merged <- SVvsNW_merged |> 
  filter(logFC > 1 & adj.P.Val <= 0.05) |> 
  select(`Abbreviation*`, adj.P.Val)

SVvsNW_merged <- merge(SVvsNW, background,
                       by.x="Molecule", by.y="Metabolite name")
SVvsNW_merged_down <- SVvsNW_merged |> 
  filter(logFC < -1 & adj.P.Val <= 0.05) |> 
  select(`Abbreviation*`, adj.P.Val)






######
#prova nuovo grafico enrichment
library(readr)

OWvsNWup_enrich <- read_csv("LION/OWvsNWup/LION-enrichment-job2.csv")
OWvsNWup_enrich$Ratio <- OWvsNWup_enrich$Significant / 336
OWvsNWup_enrich <- OWvsNWup_enrich |> rename(p.value = "p-value") |>
  mutate(Condition = "OW")

OWvsNWdown_enrich <- read_csv("LION/OWvsNWdown/LION-enrichment-job3.csv")
OWvsNWdown_enrich$Ratio <- OWvsNWdown_enrich$Significant / 377
OWvsNWdown_enrich <- OWvsNWdown_enrich |> rename(p.value = "p-value") |>
  mutate(Condition = "OW")
###
OBvsNWup_enrich <- read_csv("LION/OBvsNWup/LION-enrichment-job1.csv")
OBvsNWup_enrich$Ratio <- OBvsNWup_enrich$Significant / 345
OBvsNWup_enrich <- OBvsNWup_enrich |> rename(p.value = "p-value") |>
  mutate(Condition = "OB")

OBvsNWdown_enrich <- read_csv("LION/OBvsNWdown/LION-enrichment-job2.csv")
OBvsNWdown_enrich$Ratio <- OBvsNWdown_enrich$Significant / 384
OBvsNWdown_enrich <- OBvsNWdown_enrich |> rename(p.value = "p-value") |>
  mutate(Condition = "OB")
###
SVvsNWup_enrich <- read_csv("LION/SVvsNWup/LION-enrichment-job3.csv")
SVvsNWup_enrich$Ratio <- SVvsNWup_enrich$Significant / 295
SVvsNWup_enrich <- SVvsNWup_enrich |> rename(p.value = "p-value") |>
  mutate(Condition = "SV")

SVvsNWdown_enrich <- read_csv("LION/SVvsNWdown/LION-enrichment-job4.csv")
SVvsNWdown_enrich$Ratio <- SVvsNWdown_enrich$Significant / 407
SVvsNWdown_enrich <- SVvsNWdown_enrich |> rename(p.value = "p-value") |>
  mutate(Condition = "SV")

allUP <- rbind(OWvsNWup_enrich,OBvsNWup_enrich,SVvsNWup_enrich)
allUP <- allUP |> filter(p.value <= 0.05)
allUP$Condition <- factor(allUP$Condition, levels = c("OW", "OB", "SV"))

ggplot(data = allUP, aes(x = Condition, y = Discription, 
                                   color = Ratio, size = Significant)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red")  +
  geom_text(aes(label = round(Ratio, 3)), hjust = -0.4, vjust = 0.5, size = 3)+
  theme_bw() +
  ylab("") +
  xlab("") +
  ggtitle("LOEnA dotplot: Up-regulated")

allDOWN <- rbind(OWvsNWdown_enrich,OBvsNWdown_enrich,SVvsNWdown_enrich)
allDOWN <- allDOWN |> filter(p.value <= 0.05)
allDOWN$Condition <- factor(allDOWN$Condition, levels = c("OW", "OB", "SV"))

ggplot(data = allDOWN, aes(x = Condition, y = Discription, 
                         color = -log2(Ratio), size = Significant)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red")  +
  geom_text(aes(label = round(-log2(Ratio), 3)), hjust = -0.4, vjust = 0.5, size = 3) +
  theme_bw() +
  ylab("") +
  xlab("") +
  ggtitle("LOEnA dotplot: Down-regulated")

#questo serve solo per visualizzare un singolo listato di lipidi
ggplot(data = OWvsNWup_enrich, aes(x = Ratio, y = reorder(Discription, p.value), 
                                   color = p.value, size = Significant)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  geom_text(aes(label = round(p.value, 3)), hjust = -0.2, vjust = 0.5, size = 4) +
  theme_bw() +
  ylab("") +
  xlab("") +
  ggtitle("LOEn analysis")
