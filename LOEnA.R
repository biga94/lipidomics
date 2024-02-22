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
                           color = -log2(Ratio), size = Significant)) +
  geom_point() +
  scale_color_gradientn(colors = c("red","#2d00f7"))  +
  geom_text(aes(label = round(-log2(Ratio), 2)), hjust = -0.4, vjust = 0.5, size = 3) +
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
  scale_color_gradientn(colors = c("red","#2d00f7"))  +
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
