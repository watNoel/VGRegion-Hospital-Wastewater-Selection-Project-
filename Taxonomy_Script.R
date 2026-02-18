# ==============================
# 0. Required libraries
# ==============================
library(tidyverse)   # Data wrangling + ggplot2
library(ggsci)       # Scientific palettes (optional)
library(randomcoloR) # Distinct colors
library(scales)      # Axis formatting

# ==============================
# 1. Working directory & parameters
# ==============================
setwd("path/to/file/")

taxa_file     <- "merged_taxonomy_abundance.tsv"
metadata_file <- "metadata.csv"

tax_level <- "p__"  # Taxonomic level of interest (Phylum)
top_n     <- 10     # Top taxa per Group

# ==============================
# 2. Load taxonomy table
# ==============================
taxa_df <- read_tsv(taxa_file, col_names = TRUE)

# Rename taxonomy column for clarity
colnames(taxa_df)[1] <- "clade_name"

# ==============================
# 3. Clean sample column names
# ==============================
sample_cols <- colnames(taxa_df)[-1]

# Remove directory paths if present
clean_names <- basename(sample_cols)
colnames(taxa_df)[-1] <- clean_names

# ==============================
# 4. Separate taxonomy levels
# ==============================
taxa_df <- taxa_df %>%
  separate(
    clade_name,
    into = paste0("Level", 1:10),
    sep = "\\|",
    fill = "right",
    extra = "drop"
  )

# ==============================
# 5. Identify requested taxonomic level
# ==============================
level_col <- taxa_df %>%
  select(starts_with("Level")) %>%
  summarise(across(everything(), ~ any(str_starts(.x, tax_level)))) %>%
  pivot_longer(everything(), names_to = "level", values_to = "has_tax") %>%
  filter(has_tax) %>%
  slice(1) %>%
  pull(level)

level_num <- as.integer(str_extract(level_col, "\\d+"))

# ==============================
# 6. Extract taxa at desired level
# ==============================
taxa_df <- taxa_df %>%
  filter(
    !is.na(.data[[level_col]]) &
      is.na(.data[[paste0("Level", level_num + 1)]])  # ensures terminal phylum
  ) %>%
  mutate(
    taxon = str_replace(.data[[level_col]], paste0("^", tax_level), "")
  ) %>%
  select(taxon, all_of(clean_names))

# ==============================
# 7. Convert to long format
# ==============================
taxa_long <- taxa_df %>%
  pivot_longer(
    -taxon,
    names_to  = "filename",
    values_to = "abundance"
  )

# ==============================
# 8. Load metadata
# ==============================
meta_df <- read.csv(
  metadata_file,
  stringsAsFactors = FALSE,
  fileEncoding = "ISO-8859-1"
) %>%
  mutate(
    filename = paste0(Sample, "_R1_001.fastq.gz"),
    label    = Site
  ) %>%
  select(filename, label, Group) %>%
  mutate(Group = factor(Group))

# ==============================
# 9. Join taxonomy with metadata
# ==============================
plot_df <- taxa_long %>%
  left_join(meta_df, by = "filename") %>%
  filter(!is.na(Group))

# ==============================
# 10. Determine Top N taxa per Group
# ==============================
taxa_top <- plot_df %>%
  group_by(Group, taxon) %>%
  summarise(total_abundance = sum(abundance), .groups = "drop") %>%
  group_by(Group) %>%
  arrange(desc(total_abundance)) %>%
  slice_head(n = top_n) %>%
  ungroup() %>%
  select(Group, taxon) %>%
  mutate(taxon_group = taxon)

# ==============================
# 11. Assign Top taxa / Other
# ==============================
complete_df <- plot_df %>%
  left_join(taxa_top, by = c("Group", "taxon")) %>%
  mutate(
    taxon_group = ifelse(is.na(taxon_group), "Other", taxon_group)
  ) %>%
  group_by(label, Group, taxon_group) %>%
  summarise(abundance = sum(abundance), .groups = "drop") %>%
  group_by(label, Group) %>%
  mutate(
    abundance = abundance / sum(abundance)  # Relative abundance per sample
  ) %>%
  ungroup()

# ==============================
# 12. Global taxa ordering (VERY IMPORTANT)
# ==============================
top_taxa_order <- complete_df %>%
  group_by(taxon_group) %>%
  summarise(total_abundance = sum(abundance), .groups = "drop") %>%
  arrange(desc(total_abundance)) %>%
  pull(taxon_group)

# Always force "Other" to bottom of stack
top_taxa_order <- c(setdiff(top_taxa_order, "Other"), "Other")

complete_df <- complete_df %>%
  mutate(taxon_group = factor(taxon_group, levels = top_taxa_order))

# ==============================
# 13. Summarise for group-level plot
# ==============================
complete_df_summarized <- complete_df %>%
  group_by(Group, taxon_group) %>%
  summarise(total_abundance = sum(abundance), .groups = "drop") %>%
  group_by(Group) %>%
  mutate(
    abundance = total_abundance / sum(total_abundance)
  ) %>%
  ungroup()

# ==============================
# 14. Factor ordering (CRITICAL FOR FIGURES)
# ==============================
desired_order <- c(
  "Hospital Wastewater",
  "Hospital Biofilm",
  "Municipal Influent Wastewater",
  "Municipal Influent Biofilm",
  "Municipal Effluent Biofilm"
)

complete_df_summarized$Group <- factor(
  complete_df_summarized$Group,
  levels = desired_order
)

# ==============================
# 15. Color palette
# ==============================
n_cols <- length(levels(complete_df_summarized$taxon_group))

# Distinct colors work best for stacked bars
palette_random <- distinctColorPalette(n_cols)

# ==============================
# 16. Plot
# ==============================
text_size <- 6

tax_plot <- ggplot(
  complete_df_summarized,
  aes(x = Group, y = abundance, fill = taxon_group)
) +
  
  geom_col(position = "stack", width = 0.9) +
  
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.2),
    minor_breaks = seq(0, 1, by = 0.05),
    labels = function(x) sprintf("%.0f", x * 100),
    expand = expansion(mult = c(0, 0.02))
  ) +
  
  scale_fill_manual(values = palette_random) +
  
  labs(
    y = "Relative Abundance (%)",
    x = "Sample Group"
  ) +
  
  guides(fill = guide_legend(title = "Phylum")) +
  
  coord_cartesian(clip = "off") +
  
  theme_minimal(base_family = "Arial", base_size = text_size) +
  
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = text_size + 2),
    axis.text.x  = element_text(angle = 20,vjust = 0.5, hjust = 0.5, size = text_size + 2),
    axis.text.y  = element_text(size = text_size + 2),
    
    axis.ticks = element_line(color = "black", size = 0.3),
    
    panel.grid = element_blank(),
    
    legend.title = element_text(size = text_size + 2),
    legend.text  = element_text(size = text_size + 2),
    legend.key.size = unit(3, "mm"),
    legend.position = "right",
    legend.box.margin = margin(t = 40)
  )

print(tax_plot)

# ==============================
# 17. Export
# ==============================
ggsave(
  "Plots/tax_phyla.tiff",
  tax_plot,
  width = 183,
  height = 68.5,
  units = "mm",
  dpi = 600
)

ggsave(
  "Plots/tax_phyla.svg",
  tax_plot,
  width = 183,
  height = 68.5,
  units = "mm"
)
