
##### read in nt_content data
## nt content path
tp53_path <- "/Volumes/feder-vol1/project/li_fraumeni/results/nt_content_sequencing_panels.csv"
nt_content <- read_delim(tp53_path, delim = ",")

## convert df to long form
nt_content_long <- nt_content %>%
  pivot_longer(cols = A:T, names_to = "Nt", values_to = "Count") %>%
  group_by(Panel) %>%
  mutate(Proportion = Count/sum(Count)) %>%
  mutate(Total = sum(Count))

## GC content comparison
GC_content <- nt_content_long %>%
  filter(Nt == "G" | Nt == "C") %>%
  mutate(GC_prop = sum(Count)/Total) %>%
  filter(Nt == "G") %>%
  dplyr::select(Panel, GC_prop, Nt) %>%
  print()


## counts
# ggplot(nt_content_long, aes(x = Panel, y = Count, fill = Nt)) +
#   geom_col() +
#   scale_fill_manual(values = c("A" = "blue", "G" = "green", "C" = "red", "T" = "purple")) +  # Customize colors
#   theme_minimal() +
#   labs(title = "Nucleotide Composition by Panel",
#        x = "Panel",
#        y = "Count",
#        fill = "Nucleotide") +
#   theme(text = element_text(size = 12))



## proportion
ggplot(nt_content_long, aes(x = Panel, y = Proportion, fill = Nt)) +
  geom_col() +
  scale_fill_manual(values = c("G" = "#377EB8", "C" = "#FF7F00", "A" = "#4DAF4A", "T" = "#E41A1C")) +  # Customize colors
  theme_minimal() +
  labs(title = "Nucleotide Composition by Panel",
       x = "Panel",
       y = "Proportion",
       fill = "Nucleotide") +
  theme(text = element_text(size = 12)) +
  geom_text(data = GC_content, aes(x = Panel, y = 0.4, label = sprintf("GC = %.2f", GC_prop))) 
            #vjust = -0.5, size = 5, color = "black", fontface = "bold")

CpG <- nt_content_long %>%
  filter(Nt == "A") %>%
  mutate(Total_DN = Total-1) %>%
  mutate(CpG_prop = CpG/Total_DN) %>%
  print()


ggplot(CpG, aes(x = Panel, y = CpG_prop)) +
  geom_col() + 
  labs(title = "CpG Composition by Panel")





