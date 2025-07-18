# Load necessary libraries
library(tidyverse)      # Collection of R packages for data science
library(readxl)         # For reading Excel files
library(writexl)        # For writing Excel files

#-------------------- Load protein-coding genes --------------------#

# Read a file containing genes with protein products
protein_coding_genes <- read_delim("C:/Users/HP-ssd/Desktop/Short term project/protein coding genes/gene_with_protein_product.txt", 
                                   delim = "\t", escape_double = FALSE, 
                                   trim_ws = TRUE)

# Extract the 'symbol' column as a list of protein-coding genes
protein_coding_genes_list <- protein_coding_genes$symbol

#-------------------- Load FUSIL data --------------------#

# Read in FUSIL dataset
fusil_m_gene <-  read_delim("C:/Users/HP-ssd/Desktop/Short term project2/fusil.csv")

# Check how many unique gene symbols are in FUSIL data
length(unique(fusil_m_gene$gene_symbol))

#-------------------- Load gene paralogues from Biomart --------------------#

# Read paralogues file
human_gene_paralogues <- read.csv("C:/Users/HP-ssd/Desktop/Short term project2/paralogues/human_gene_paralogues.csv")

# Drop unnecessary columns and rename for consistency
human_gene_paralogues <- human_gene_paralogues %>%
  dplyr::select(-1,-2,-4)%>%
  rename(gene_symbol = external_gene_name )

# Check number of unique genes with paralogues
length(unique(human_gene_paralogues$gene_symbol))

#-------------------- Merge FUSIL with paralogue data --------------------#

# Left join FUSIL info with paralogue info
paralogue_fusil <- human_gene_paralogues %>%
  left_join(fusil_m_gene)%>%
  dplyr::select(-6, -7) %>% # Drop unnecessary columns
  mutate(hsapiens_paralog_associated_gene_name = na_if(hsapiens_paralog_associated_gene_name,"")) %>%
  distinct()

# Check number of unique genes in merged data
length(unique(paralogue_fusil$gene_symbol))

#-------------------- Gene count per FUSIL bin --------------------#

# Count genes in each FUSIL bin and calculate percentages
gene_counts <- fusil_m_gene %>%
  count(fusil) %>%
  mutate(percentage = (n/sum(n)*100)) %>%
  mutate(total_gene_count =n)

# Reorder factor levels for plotting
gene_counts$fusil <- factor(gene_counts$fusil, 
                            levels = c("CL", "DL", "SV", "VP", "VnP"  ))

# Plot gene distribution by FUSIL
ggplot(gene_counts, aes(x= fusil, y=percentage, fill = fusil))+
  geom_bar(stat = "identity", position = "dodge")+
  labs(y = "Percentage of Genes", x = "FUSIL bin", title = "Gene Counts by FUSIL Category") +
  theme_minimal()+
  scale_fill_manual(
    values = c(
      "CL" = "#E41A1C",
      "DL" = "#377EB8",
      "SV" = "#4DAF4A",
      "VnP" = "#984EA3",
      "VP" = "#FF7F00"
    ))

#-------------------- Genes with vs without paralogues in each FUSIL --------------------#

# Filter genes with no paralogues (i.e., NA similarity values)
na_paralogue_count <- paralogue_fusil %>%
  group_by(gene_symbol) %>%
  filter(all(is.na(hsapiens_paralog_perc_id)| hsapiens_paralog_perc_id < 30 )) %>%
  ungroup() %>%
  select(1,6) %>%
  distinct() %>%
  count(fusil) %>%
  mutate(Gene_with_no_paralogues =n)

# Total count of such genes
sum(na_paralogue_count$n)

# Filter genes with at least one paralogue
Genes_with_paralogue <- paralogue_fusil %>%
  group_by(gene_symbol) %>%
  na.omit() %>%
  filter(hsapiens_paralog_perc_id > 30) %>%
  ungroup() %>%
  distinct() %>%
  select(1,6) %>%
  distinct()%>%
  count(fusil) %>%
  mutate(Genes_with_paralogue =n)


# Total count of these genes
sum(gene_with_paralogue_count$n)

# Combine paralogue info with total gene counts
Summary_gene_count <- gene_counts %>%
  left_join(na_paralogue_count, by = "fusil") %>%
  left_join(Genes_with_paralogue, by = "fusil") %>%
  select(1,4,6,8)

# Convert wide table to long format for plotting
Summary_gene_count_long <- Summary_gene_count %>%
  pivot_longer(cols = c(Genes_with_paralogue, Gene_with_no_paralogues),
               names_to = "Paralogue Status",
               values_to = "Count")

Summary_gene_count_long <- Summary_gene_count_long %>%
  group_by(fusil) %>%
  mutate(percentage = (Count/total_gene_count)*100)


# Plot stacked bar chart comparing genes with/without paralogues by FUSIL

Summary_gene_count_long$fusil <- factor(Summary_gene_count_long$fusil, 
                                        levels = c("CL", "DL", "SV", "VP", "VnP" ))

df_totals <- Summary_gene_count %>%
  select(fusil, total_gene_count) %>%
  distinct()

ggplot(Summary_gene_count_long, aes(x= fusil, y= Count, fill = `Paralogue Status`))+
  geom_bar(stat = "identity") +
  geom_text(data = df_totals, aes(x = fusil, y = total_gene_count, label = total_gene_count),
            vjust = -0.5, inherit.aes = FALSE) +
  labs(
    title = "Gene Counts by Paralogue Presence",
    x = "FUSIL bin", y = "Gene Count"
  )

# Create the bar plot with percentages

ggplot(Summary_gene_count_long, aes(x = fusil, y = percentage, fill = `Paralogue Status`)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_text(aes(label = round(percentage, 1)),   # Round for cleaner labels
            position = position_dodge(width = 0.9),
            vjust = -0.5, size = 3) +
  labs(
    title = "Paralogue Presence Across FUSIL Bins",
    x = "FUSIL Category",
    y = "Percentage of Genes",
    fill = "Paralogue Status"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +  # Prevent bars from touching top
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

############ Paralogue counts Per FUSIL Category while accounting for protein coding genes --------

# Remove NA entries and filter to only include paralogues that are protein-coding genes
paralogue_fusil_filtered <- paralogue_fusil %>%
  na.omit() %>%
  filter(hsapiens_paralog_associated_gene_name %in% protein_coding_genes_list)

# Check number of unique genes in the filtered dataset
length(unique(paralogue_fusil_filtered$gene_symbol))

# Set FUSIL levels to maintain consistent order in plots
paralogue_fusil_filtered$fusil <- factor(paralogue_fusil_filtered$fusil, 
                                         levels = c("VnP", "VP", "SV", "DL","CL"  ))

# Define thresholds for similarity percentages
thresholds <- c(30,50,70)

# For each threshold, count paralogues with similarity >= threshold per gene and FUSIL bin
fusil_cat_thresh <-  lapply(thresholds, function(thresh) {
  paralogue_fusil_filtered %>%
    filter(hsapiens_paralog_perc_id >= thresh) %>%
    group_by(gene_symbol, fusil) %>%
    tally(name = "Paralogues_above_threshold") %>%
    mutate(Threshold =thresh)
}) %>%
  bind_rows()

# Check number of unique genes across all thresholds
length(unique(fusil_cat_thresh$gene_symbol))

# Summarize average number of paralogues per FUSIL bin and threshold
summary_thresh <- fusil_cat_thresh %>%
  group_by(fusil, Threshold) %>%
  summarise(avg_paralogues = mean(Paralogues_above_threshold), .groups= "drop")

# Set FUSIL bin order for plot
summary_thresh$fusil <- factor(summary_thresh$fusil, 
                               levels = c("CL", "DL", "SV", "VP", "VnP"  ))

# Bar plot: average number of paralogues per FUSIL bin by threshold
ggplot(summary_thresh, aes(x= factor(Threshold), y=avg_paralogues, fill = fusil))+
  geom_bar(stat = "identity", position = "dodge")+
  facet_wrap(~Threshold)+
  labs(x = "Homology Threshold (%)", y = "Average Number of Paralogues", title = "Paralogue Count by FUSIL Category") +
  theme_minimal()+
  theme(axis.text.x = element_blank())+
  scale_fill_manual(
    values = c(
      "CL" = "#E41A1C",
      "DL" = "#377EB8",
      "SV" = "#4DAF4A",
      "VnP" = "#984EA3",
      "VP" = "#FF7F00"
    ))

### Plot density distribution of similarity scores across FUSIL categories

# Set FUSIL bin order again for consistency
paralogue_fusil_filtered$fusil <- factor(paralogue_fusil_filtered$fusil, 
                                         levels = c("CL", "DL", "SV", "VP", "VnP" ))

# Density plot of paralogue similarity scores by FUSIL bin
ggplot(paralogue_fusil_filtered, aes( x= hsapiens_paralog_perc_id, fill = fusil ))+
  geom_density(alpha = 0.5) +
  labs(title = "Distribution of Paralogue Similarity by FUSIL Category",
       x = "Similarity Score",
       y = "Density")+
  theme_minimal()

#-------------------- Binary classification of paralogues by similarity threshold --------------------#

library(stringr)  # For regex string operations

# Custom function to determine if a gene has a paralogue above a given threshold
has_paralogues_above_thresholds <- function(gene, paralogue_perc_id, pc_genes, threshold) {
  
  if (!(gene %in% pc_genes)) {
    return(0)
  }
  
  if(is.na(paralogue_perc_id) || paralogue_perc_id == "") {
    return(0)
  }
  
  perc <- as.numeric(str_extract(paralogue_perc_id, "\\d+\\.?\\d*")[[1]])
  
  if(!is.na(perc) && perc >= threshold ) {
    return(1)
  } else {
    return(0)
  }
}

# Define thresholds
thresholds <- c(30, 50, 70)

# Extract list of protein coding gene symbols
pc_genes <- protein_coding_genes$symbol

# Apply binary classification for each threshold
binary_thresh <- lapply(thresholds, function (thresh) {
  
  binary_column <- mapply(has_paralogues_above_thresholds,
                          gene = paralogue_fusil$gene_symbol,
                          paralogue_perc_id = paralogue_fusil$hsapiens_paralog_perc_id,
                          MoreArgs = list(pc_genes = pc_genes, threshold = thresh)
  )
  col_name <- paste0("Has_paralogue_over_", thresh)
  setNames(data.frame(binary_column), col_name)
})

# Combine binary classification results
binary_columns_df <- bind_cols(binary_thresh)

# Merge with original paralogue + FUSIL data
fusil_cat_thresh_2 <- bind_cols(paralogue_fusil, binary_columns_df)

# Count unique genes after binary transformation
length(unique(fusil_cat_thresh_2$gene_symbol))

# Summarize binary values per gene (if any paralogue meets threshold, keep 1)
binary_summary_per_gene <- fusil_cat_thresh_2 %>%
  dplyr::select(gene_symbol, fusil, starts_with("Has_paralogue_over_")) %>%
  group_by(gene_symbol, fusil) %>%
  summarise(across(starts_with("Has_paralogue_over_"), max, na.rm = TRUE), .groups = "drop")

# Count number of paralogues meeting thresholds
count_summary_per_gene <- fusil_cat_thresh_2 %>%
  dplyr::select(gene_symbol,fusil, starts_with("Has_paralogue_over_")) %>%
  group_by(gene_symbol, fusil) %>%
  summarise(across(starts_with("Has_paralogue_over_"), sum, na.rm = TRUE), .groups = "drop") %>%
  rename_with(~ gsub("Has", "Count", .x), starts_with("Has"))

#-------------------- Visualizations of binary paralogue data --------------------#

# Convert to long format for plotting
binary_summary_per_gene_long <- binary_summary_per_gene %>%
  pivot_longer( cols = starts_with("Has_paralogue_over_"),
                names_to = "Threshold",
                values_to = "Has_Paralogue"
  )

# Clean up threshold names
binary_summary_per_gene_long <- binary_summary_per_gene_long %>%
  mutate(Threshold = as.numeric(gsub("Has_paralogue_over_", "", Threshold)))

# Summarize proportions of 1s and 0s per FUSIL bin and threshold
plot_binary_summary_per_gene_long <- binary_summary_per_gene_long %>%
  group_by(fusil, Threshold, Has_Paralogue ) %>%
  summarise(Count=n(), , .groups = "drop") %>%
  group_by(fusil, Threshold) %>%
  mutate(Proportions = Count/sum(Count))

# Proportion bar plot by threshold (x-axis = threshold)
ggplot(plot_binary_summary_per_gene_long, aes( x = factor(Threshold), y = Proportions, fill = factor(Has_Paralogue)))+
  geom_col(position = "dodge") +
  facet_wrap(~ fusil)+
  labs(
    title = "Proportion of Genes With and Without Paralogue by FUSIL Class",
    x = "Similarity Threshold (%)",
    y = "Proportion of Genes",
    fill = "Has Paralogue"
  ) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("0" = "#FF9999", "1" = "#66CC99"),
                    labels = c("No", "Yes")) +
  theme_minimal(base_size = 14)

# Proportion bar plot by FUSIL (x-axis = FUSIL)
ggplot(plot_binary_summary_per_gene_long, aes( x = factor(fusil), y = Proportions, fill = factor(Has_Paralogue)))+
  geom_col(position = "dodge") +
  facet_wrap(~ Threshold)+
  labs(
    title = "Proportion of Genes With and Without Paralogue by FUSIL Class",
    x = "FUSIL bin",
    y = "Proportion of Genes",
    fill = "Has Paralogue"
  ) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("0" = "#FF9945", "1" = "#690C99"),
                    labels = c("No", "Yes")) +
  theme_minimal(base_size = 14)

# Absolute count version (instead of proportion)
ggplot(plot_binary_summary_per_gene_long, aes( x = factor(fusil), y = Count, fill = factor(Has_Paralogue)))+
  geom_col(position = "dodge") +
  facet_wrap(~ Threshold)+
  labs(
    title = "Proportion of Genes With and Without Paralogue by FUSIL Class",
    x = "FUSIL bin",
    y = "Proportion of Genes",
    fill = "Has Paralogue"
  ) +
  scale_fill_manual(values = c("0" = "#FF9945", "1" = "#690C99"),
                    labels = c("No", "Yes")) +
  theme_minimal(base_size = 14)

# Log-transformed count version
ggplot(plot_binary_summary_per_gene_long, aes(
  x = factor(fusil),
  y = Count,
  fill = factor(Has_Paralogue))
) +
  geom_col(position = "dodge") +
  facet_wrap(~ Threshold) +
  labs(
    title = "Proportion of Genes With and Without Paralogue by FUSIL Class",
    x = "FUSIL bin",
    y = "Log10(Count of Genes)",
    fill = "Has Paralogue"
  ) +
  scale_y_log10() +
  scale_fill_manual(
    values = c("0" = "#FF9945", "1" = "#690C99"),
    labels = c("No", "Yes")
  ) +
  theme_minimal(base_size = 14)


##################### FUSIL comparison for Genes and their Paralogues --------------------
# Are the genes and their paralogues in the same FUSIL Category?

# Create a table with FUSIL categories for both genes and their paralogues

# Select only gene symbol and FUSIL bin from the main gene list
Fusil_genes <- fusil_m_gene %>%
  dplyr::select(-1,-2)

# Join paralogue data to get FUSIL category of each paralogue gene
Fusil_genes_paralogues <- human_gene_paralogues %>%
  left_join(Fusil_genes, by = c("hsapiens_paralog_associated_gene_name" = "gene_symbol")) %>%
  filter(hsapiens_paralog_associated_gene_name %in% protein_coding_genes_list) %>%
  rename("fusil_paralogue" = "fusil")

# Join again to get FUSIL category of the original gene
Fusil_genes_paralogues <- Fusil_genes_paralogues %>%
  left_join(Fusil_genes, by = c("gene_symbol" = "gene_symbol"))

# Rearranging columns and removing any rows with missing data
Fusil_genes_paralogues <- Fusil_genes_paralogues %>%
  relocate(fusil, .after = "gene_symbol") %>%
  relocate(fusil_paralogue, .after = "hsapiens_paralog_associated_gene_name") %>%
  na.omit()

#---------------------- PLOTTING ----------------------------#

# Load additional libraries for alluvial and heatmap plots
library(ggalluvial) # For Sankey-like alluvial plots
library(reshape2)   # For generating heatmaps (optional)

# Count genes per FUSIL category
df_fusil <- Fusil_genes_paralogues %>%
  count(fusil)

# Count FUSIL → FUSIL_PARALOGUE transitions and compute percentages
fusil_table <- Fusil_genes_paralogues %>%
  count(fusil, fusil_paralogue) %>%
  left_join(df_fusil, by = c("fusil" = "fusil")) %>%
  mutate(percentage = n.x * 100 / n.y)

# Sankey-like alluvial plot showing FUSIL transitions between genes and their paralogues
ggplot(fusil_table, aes(axis1 = fusil, axis2 = fusil_paralogue, y = percentage / 2)) +
  geom_alluvium(aes(fill = fusil), width = 0.2, show.legend = FALSE) +
  geom_stratum(width = 0.2) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Gene FUSIL", "Paralog FUSIL"), expand = c(.01, .01)) +
  scale_fill_manual(
    values = c("CL" = "#E41A1C", "DL" = "#377EB8", "SV" = "#4DAF4A", "VnP" = "#984EA3", "VP" = "#FF7F00")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    plot.margin = margin(5, 5, 5, 5)
  ) +
  ggtitle("FUSIL Category Flow: Gene → Paralog") +
  geom_text(aes(label = paste0(round(percentage, 1), "%")),
            stat = "alluvium",
            nudge_x = 0.15,     # move to the left
            size = 3,
            fontface = "bold",  # make text bold
            color = "black")



#---------------------- Similarity bins and FUSIL match ----------------------#

# Check whether each gene and its paralogue are in the same FUSIL bin
fusil_match <- Fusil_genes_paralogues %>%
  mutate(FUSIL_match = ifelse(fusil == fusil_paralogue, "Match", "Mismatch")) %>%
  mutate(SIMILARITY_bin = case_when(
    hsapiens_paralog_perc_id >= 80 ~ "High >80%",
    hsapiens_paralog_perc_id >= 60 ~ "Medium-High 60-80%",
    hsapiens_paralog_perc_id >= 40 ~ "Medium 40-60%",
    hsapiens_paralog_perc_id >= 20 ~ "Medium-Low 20-50%",
    TRUE ~ "Low <20%"
  ))

#---------------------- Alluvial plots by similarity thresholds ----------------------#

thresholds <- c(30, 50, 70)

# Generate Sankey plot for each similarity threshold
for (bin in thresholds) {
  
  cat("\n===== SIMILARITY_bin:", bin, "=====\n")
  
  bin_data <- fusil_match %>%
    filter(hsapiens_paralog_perc_id >= bin)
  
  df_fusil2 <- bin_data %>%
    count(fusil)
  
  fusil_table2 <- bin_data %>%
    count(fusil, fusil_paralogue) %>%
    left_join(df_fusil2, by = c("fusil" = "fusil")) %>%
    mutate(percentage = n.x * 100 / n.y)
  
  p <- ggplot(fusil_table2, aes(axis1 = fusil, axis2 = fusil_paralogue, y = percentage)) +
    geom_alluvium(aes(fill = fusil), width = 0.1, show.legend = FALSE) +
    geom_stratum(width = 0.1) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
    scale_x_discrete(limits = c("Gene FUSIL", "Paralog FUSIL"), expand = c(.05, .05)) +
    scale_fill_manual(values = c("CL" = "#E41A1C", "DL" = "#377EB8", "SV" = "#4DAF4A", "VnP" = "#984EA3", "VP" = "#FF7F00")) +
    theme_minimal() +
    theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
    ggtitle(paste("Gene vs Paralog FUSIL -", bin, "Similarity")) +
    geom_text(aes(label = paste0(round(percentage, 1), "%")),
              stat = "alluvium",
              nudge_x = 0.2,
              size = 3,
              color = "black")
  
  print(p)
}

#---------------------- % Similarity bin distribution per FUSIL bin ----------------------#

# Total paralogue counts per FUSIL pair (gene, paralogue)
similarity_df <- fusil_match %>%
  group_by(fusil, fusil_paralogue) %>%
  tally()

# Count of SIMILARITY_bin types per gene/paralog FUSIL pair
similarity_paralogues <- fusil_match %>%
  group_by(fusil, fusil_paralogue, SIMILARITY_bin) %>%
  tally() %>%
  left_join(similarity_df, by = c("fusil", "fusil_paralogue")) %>%
  mutate(percentage = (n.x / n.y) * 100) %>%
  na.omit()

# Set FUSIL order for plot aesthetics
similarity_paralogues$fusil <- factor(similarity_paralogues$fusil,
                                      levels = c("CL", "DL", "SV", "VP", "VnP"))
similarity_paralogues$fusil_paralogue <- factor(similarity_paralogues$fusil_paralogue,
                                                levels = c("CL", "DL", "SV", "VP", "VnP"))

# Plot percentage of similarity bins between FUSIL gene-paralogue pairs
ggplot(similarity_paralogues, aes(x = SIMILARITY_bin, y = percentage, fill = fusil_paralogue)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~fusil) +
  labs(title = "Similarity Distribution Between Paralogues",
       x = "Similarity Bin", y = "Percentage") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("CL" = "#E41A1C", "DL" = "#377EB8", "SV" = "#4DAF4A", "VnP" = "#984EA3", "VP" = "#FF7F00"))


#-------------------- FUSIL gene with MCRA --------------------------------

# Get list of unique FUSIL gene symbols
fusil_m_gene_list <- unique(fusil_m_gene$gene_symbol)

# Load BioMart interface to access Ensembl
library(biomaRt)

# Connect to Ensembl BioMart, human dataset
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# View available attributes (optional; shows what data we can pull)
attributes_human <- listAttributes(human)

# Retrieve paralogue subtype (e.g., ancestral gene family) for FUSIL genes
fusil_gene_mcra <- getBM(
  attributes = c("external_gene_name", "hsapiens_paralog_subtype"),
  filters = "external_gene_name",
  values = fusil_m_gene_list,
  mart = human
)

# (This line replaces the earlier result?) — loads from local paralogues dataset again
fusil_gene_mcra <- human_gene_paralogues %>%
  rename("external_gene_name" = "gene_symbol")

# Filter out rows where paralog subtype is empty
fusil_gene_mcra_df <- fusil_gene_mcra %>%
  filter(!(hsapiens_paralog_subtype == ""))

# Merge MCRA info with FUSIL classification
fusil_gene_mcra_df2 <- fusil_m_gene %>%
  left_join(fusil_gene_mcra_df, by = c("gene_symbol" = "external_gene_name")) %>%
  dplyr::select(-1,-2) %>% # Drop unneeded columns
  na.omit()  # Remove incomplete rows

# Count genes by paralogue subtype
count_mcra_fusil <- fusil_gene_mcra_df2 %>%
  count(hsapiens_paralog_subtype)

# Focus on broad paralogue subtypes (higher-order evolutionary clades)
broad_subtypes <- c("Bilateria", "Chordata", "Gnathostomata", "Opisthokonta", "Vertebrata")

# Count how many genes per FUSIL bin have MCRA in each broad subtype
mcra_df_fusil <- fusil_gene_mcra_df2  %>%
  filter(hsapiens_paralog_subtype %in% broad_subtypes) %>%
  count(fusil)

# Calculate percentages of each MCRA subtype within FUSIL bins
mcra_fusil <- fusil_gene_mcra_df2 %>%
  filter(hsapiens_paralog_subtype %in% broad_subtypes) %>%
  group_by(fusil, hsapiens_paralog_subtype) %>%
  tally() %>%
  left_join(mcra_df_fusil, by =c ("fusil" ="fusil")) %>%
  mutate(percentage = (n.x/n.y)*100)

# Set factor levels for consistent plotting order
mcra_fusil$hsapiens_paralog_subtype <-factor(mcra_fusil$hsapiens_paralog_subtype, 
                                             levels = c ("Opisthokonta" ,"Bilateria" ,"Chordata","Vertebrata", "Gnathostomata"))

mcra_fusil$fusil <- factor(mcra_fusil$fusil,
                           levels = c("CL", "DL", "SV", "VP", "VnP" ))

#-------------------- Plot MCRA subtype distribution by FUSIL bin --------------------

ggplot(mcra_fusil, aes(x = factor(hsapiens_paralog_subtype), y = percentage, fill = fusil))+
  geom_col(position = "dodge") +
  labs(
    title = "FUSIL Gene's Paralogues Distribution",
    y = "Percentage of Paralogue Genes",
    x = "MCRA Subtype Group",
    fill = "Fusil"
  ) +
  theme_minimal(base_size = 14) +
  theme_minimal() +  # (duplicated but harmless)
  scale_fill_manual(
    values = c(
      "CL" = "#E41A1C",
      "DL" = "#377EB8",
      "SV" = "#4DAF4A",
      "VnP" = "#984EA3",
      "VP" = "#FF7F00"
    ))

# Plot MCRA subtypes grouped by FUSIL category instead
ggplot(mcra_fusil, aes(x = fusil, y = percentage, fill = hsapiens_paralog_subtype))+
  geom_col(position = "dodge") +
  labs(
    title = "FUSIL Gene Distribution by MCRA Subtype",
    y = "Percentage of Paralogues Genes",
    x = "FUSIL bin",
    fill = "MCRA Subtype"
  ) +
  theme_minimal(base_size = 14)
