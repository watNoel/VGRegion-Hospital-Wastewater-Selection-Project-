# ==============================
# 0. Required libraries
# ==============================
library(vegan)          
library(ggplot2)        
library(dplyr)          
library(tibble)         
library(pairwiseAdonis) 
library(tidyverse)      
library(ggpubr)         
library(grid)           
library(svglite)        

setwd("path/to/file")

# ==============================
# 1. File paths & parameters
# ==============================
taxa_file <- "merged_taxonomy_abundance.tsv"

meta_df <- read.csv(
  "metadata.csv",
  stringsAsFactors = FALSE,
  fileEncoding = "ISO-8859-1"
)

tax_level <- "g__"  

# ==============================
# 2. Load taxonomy table
# ==============================
taxa_df <- read_tsv(taxa_file, col_names = TRUE)
colnames(taxa_df)[1] <- "clade_name"

sample_cols <- colnames(taxa_df)[-1]
clean_names <- basename(sample_cols)
colnames(taxa_df)[-1] <- clean_names

# ==============================
# 3. Separate taxonomy levels
# ==============================
taxa_df <- taxa_df %>%
  separate(
    clade_name,
    into = paste0("Level", 1:10),
    sep = "\\|",
    fill = "right",
    extra = "drop"
  )

level_col <- taxa_df %>%
  select(starts_with("Level")) %>%
  summarise(across(everything(), ~ any(str_starts(.x, tax_level)))) %>%
  pivot_longer(everything(), names_to = "level", values_to = "has_tax") %>%
  filter(has_tax) %>%
  slice(1) %>%
  pull(level)

level_num <- as.integer(str_extract(level_col, "\\d+"))

# ==============================
# 4. Filter to genus only
# ==============================
taxa_df <- taxa_df %>%
  filter(
    !is.na(.data[[level_col]]) &
      is.na(.data[[paste0("Level", level_num + 1)]])
  ) %>%
  mutate(taxon = str_replace(.data[[level_col]], paste0("^", tax_level), "")) %>%
  select(taxon, all_of(clean_names))

# ==============================
# 5. Metadata preparation
# ==============================
meta_df <- meta_df %>%
  mutate(filename = paste0(Sample, "_R1_001.fastq.gz")) %>%
  select(filename, Site, Type_Site, Type, Group)

# ==============================
# 6. Match samples
# ==============================
matched_samples <- intersect(meta_df$filename, colnames(taxa_df))

abund_matrix <- taxa_df %>%
  select(taxon, all_of(matched_samples)) %>%
  pivot_longer(-taxon, names_to = "filename", values_to = "abundance") %>%
  pivot_wider(names_from = taxon, values_from = abundance, values_fill = 0) %>%
  column_to_rownames("filename")

meta_df <- meta_df %>% filter(filename %in% rownames(abund_matrix))
abund_matrix <- abund_matrix[meta_df$filename, ]

# ==============================
# 7. Bray–Curtis distance
# ==============================
bray_dist <- vegdist(abund_matrix, method = "bray")

# ==============================
# 8. PCoA
# ==============================
pcoa <- cmdscale(bray_dist, k = 2, eig = TRUE)
pcoa_df <- as.data.frame(pcoa$points)
colnames(pcoa_df) <- c("PC1", "PC2")

pcoa_df <- pcoa_df %>%
  mutate(
    filename  = rownames(.),
    Site      = meta_df$Site,
    Type_Site = meta_df$Type_Site,
    Type      = meta_df$Type,
    Group     = meta_df$Group
  )

# ==============================
# 9. Convex hulls
# ==============================
hulls <- pcoa_df %>%
  group_by(Group) %>%
  slice(chull(PC1, PC2))

# ==============================
# 10. Factor ordering
# ==============================
desired_order <- c(
  "Hospital Wastewater",
  "Hospital Biofilm",
  "Municipal Influent Wastewater",
  "Municipal Influent Biofilm",
  "Municipal Effluent Biofilm"
)

pcoa_df$Group <- factor(pcoa_df$Group, levels = desired_order)
hulls$Group   <- factor(hulls$Group,   levels = desired_order)

# ==============================
# 11. Variance explained
# ==============================
eig_vals <- pcoa$eig
var_explained <- eig_vals / sum(eig_vals[eig_vals > 0])
pc1_var <- round(var_explained[1] * 100, 2)
pc2_var <- round(var_explained[2] * 100, 2)

# ==============================
# 12. PERMANOVA
# ==============================
permanova_res <- adonis2(bray_dist ~ Group, data = meta_df, permutations = 999)
permanova_p <- permanova_res$`Pr(>F)`[1]
sig_label <- ifelse(permanova_p < 0.001, "***",
                    ifelse(permanova_p < 0.01, "**",
                           ifelse(permanova_p < 0.05, "*", "ns")))

# ==============================
# 13. Beta dispersion
# ==============================
bd_group <- betadisper(bray_dist, meta_df$Group)
beta_disp_anova <- anova(bd_group)

# ==============================
# 14. Pairwise PERMANOVA
# ==============================
pairwise_res <- NULL
if(permanova_p < 0.05){
  pairwise_res <- pairwise.adonis2(bray_dist ~ Group, data = meta_df)
}

# ==============================
# 15. PCoA Plot
# ==============================
text_size <- 6

p_beta <- ggplot(pcoa_df, aes(PC1, PC2)) +
  geom_point(aes(color = Site, shape = Type_Site), size = 2, alpha = 0.8) +
  geom_polygon(data = hulls, aes(fill = Group, group = Group), alpha = 0.2, color = NA) +
  scale_shape_manual(values = c(16, 17)) +  # circles & triangles
  guides(
    shape = guide_legend(
      title = "Type",
      override.aes = list(size = 3, alpha = 1),
      keyheight = unit(3, "mm"),
      keywidth  = unit(3, "mm")
    ),
    color = guide_legend(
      title = "Site",
      override.aes = list(shape = 16, size = 3, alpha = 1),
      keyheight = unit(3, "mm"),
      keywidth  = unit(3, "mm")
    ),
    fill = guide_legend(
      title = "Group",
      override.aes = list(alpha = 1),
      keyheight = unit(3, "mm"),
      keywidth  = unit(3, "mm")
    )
  ) +
  labs(
    x = paste0("PC1 (", pc1_var, "%)"),
    y = paste0("PC2 (", pc2_var, "%)")
  ) +
  annotate(
    "text",
    x = min(pcoa_df$PC1),
    y = max(pcoa_df$PC2),
    label = paste0("PERMANOVA p = ", signif(permanova_p, 2), " (", sig_label, ")"),
    hjust = 0, vjust = 1,
    size = text_size / .pt
  ) +
  theme_minimal(base_family = "Arial", base_size = text_size) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    legend.spacing.y = unit(0, "mm"),
    legend.spacing.x = unit(0, "mm"),
    plot.margin = margin(10, 10, 5, 5)
  ) +
  coord_cartesian(clip = "off")

print(p_beta)

# ==============================
# 16. Export
# ==============================
ggsave(
  "Plots/beta_div.tiff",
  p_beta,
  width = 183, height = 100,
  units = "mm",
  dpi = 600
)

ggsave(
  "Plots/beta_div.svg",
  p_beta,
  width = 183, height = 100,
  units = "mm",
  dpi = 600
)
