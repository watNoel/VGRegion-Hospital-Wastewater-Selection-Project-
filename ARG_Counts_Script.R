# ==============================
# 0. Required libraries
# ==============================
library(ggplot2)   # Plotting
library(dplyr)     # Data wrangling
library(readr)     # File import
library(tidyverse) # Pipes + tidy helpers
library(stringi)   # String cleaning
library(RColorBrewer) # Color palettes
library(showtext)  # Proper Arial embedding
library(svglite)   # SVG export

# ==============================
# 1. Working directory
# ==============================
setwd("path/to/file")

# ==============================
# 2. Font handling (critical for journals)
# ==============================
font_add("Arial", regular = "Arial.ttf")
showtext_auto()

# ==============================
# 3. Load metadata
# ==============================
metadata <- read.csv(
  "metadata.csv",
  stringsAsFactors = FALSE,
  fileEncoding = "ISO-8859-1"
)

# ==============================
# 4. Import & process ResFinder data
# ==============================
import_resfinder <- function(filter_class = FALSE, filter_gene = FALSE) {
  
  # Load raw ARG counts
  args_df_RF <- read_tsv(
    "ResFinder_DB_counts_by_group.tsv",
    col_types = cols(),
    locale = locale(encoding = "UTF-8")
  )
  
  # Load phenotype mapping
  phenotypes_RF <- read_tsv("phenotypes.txt")
  
  # ------------------------------
  # Clean ARG naming artifacts
  # ------------------------------
  args_df_RF <- args_df_RF %>%
    mutate(
      Class = str_replace(Class, "_1$", "")  # Remove duplicate suffixes
    )
  
  # ------------------------------
  # Standardize gene names
  # ------------------------------
  args_df_RF <- args_df_RF %>%
    mutate(Class_clean = str_replace(Class, "_.*$", ""))
  
  phenotypes_RF <- phenotypes_RF %>%
    mutate(Gene_clean = str_replace(Gene, "_.*$", ""))
  
  phenotypes_RF_unique <- phenotypes_RF %>%
    group_by(Gene_clean) %>%
    slice(1) %>%
    ungroup()
  
  # ------------------------------
  # Join phenotype classification
  # ------------------------------
  args_df_RF <- left_join(
    args_df_RF,
    phenotypes_RF_unique,
    by = c("Class_clean" = "Gene_clean")
  )
  
  # ------------------------------
  # Fix known ResFinder quirks
  # ------------------------------
  args_df_RF$Class_clean <- gsub("^tet\\(O$", "tet(O)", args_df_RF$Class_clean)
  args_df_RF$Phenotype[args_df_RF$Class_clean == "tet(O)"] <- "Tetracycline"
  args_df_RF$Class.y[args_df_RF$Class_clean == "tet(O)"] <- "Tetracycline"
  
  args_df_RF$Class_clean <- gsub("^tet\\(S$", "tet(S)", args_df_RF$Class_clean)
  args_df_RF$Phenotype[args_df_RF$Class_clean == "tet(S)"] <- "Tetracycline"
  args_df_RF$Class.y[args_df_RF$Class_clean == "tet(S)"] <- "Tetracycline"
  
  # ------------------------------
  # Optional filtering
  # ------------------------------
  if (!identical(filter_class, FALSE)) {
    args_df_RF <- args_df_RF %>% filter(Class.y %in% filter_class)
  }
  
  if (!identical(filter_gene, FALSE)) {
    args_df_RF <- args_df_RF %>%
      filter(str_detect(Class, paste(filter_gene, collapse = "|")))
  }
  
  return(args_df_RF)
}

# Load processed ARG dataset
args_df_RF <- import_resfinder(FALSE, FALSE)

# ==============================
# 5. Normalize counts (CPM)
# ==============================
process_microbial_data <- function(
    count_df,
    meta_df,
    group_column = "Class.y",
    x_var = "Site",
    top_n = 10
) {
  
  # ------------------------------
  # Join counts with metadata
  # ------------------------------
  normalized_df <- meta_df %>%
    left_join(count_df, by = "Sample") %>%
    mutate(
      Count = ifelse(is.na(Count), 0, Count),
      
      # Counts Per Million (CPM)
      CPM = (Count / `Total.Reads_After`) * 1e6,
      
      # Defensive cleaning (very important)
      CPM = ifelse(is.na(CPM) | is.infinite(CPM) | CPM < 0, 0, CPM)
    ) %>%
    filter(!is.na(Group))
  
  # ------------------------------
  # Average CPM per grouping
  # ------------------------------
  avg_data <- normalized_df %>%
    group_by(!!sym(x_var), !!sym(group_column), Group) %>%
    summarise(avg_CPM = mean(CPM), .groups = "drop")
  
  # ------------------------------
  # Top-N ARG classes
  # ------------------------------
  top_groups <- avg_data %>%
    group_by(!!sym(group_column)) %>%
    summarise(total_CPM = sum(avg_CPM), .groups = "drop") %>%
    arrange(desc(total_CPM)) %>%
    slice(1:top_n) %>%
    pull(!!sym(group_column))
  
  avg_data <- avg_data %>%
    mutate(
      !!sym(group_column) :=
        ifelse(!!sym(group_column) %in% top_groups,
               as.character(!!sym(group_column)),
               "Other"),
      
      !!sym(group_column) :=
        factor(!!sym(group_column),
               levels = c(sort(top_groups), "Other"))
    ) %>%
    group_by(across(-avg_CPM)) %>%
    summarise(avg_CPM = sum(avg_CPM), .groups = "drop")
  
  # ------------------------------
  # Color palette
  # ------------------------------
  groups <- levels(avg_data[[group_column]])
  num_main_groups <- length(groups) - 1
  
  custom_palette <- c(
    colorRampPalette(brewer.pal(12, "Paired"))(num_main_groups),
    "gray70"
  )
  names(custom_palette) <- groups
  
  # ==============================
  # 6. Plot
  # ==============================
  text_size <- 7
  
  avg_data$Group <- factor(
    avg_data$Group,
    levels = c(
      "Hospital Wastewater",
      "Hospital Biofilm",
      "Municipal Influent Wastewater",
      "Municipal Influent Biofilm",
      "Municipal Effluent Biofilm"
    )
  )
  
  p <- ggplot(
    avg_data,
    aes(x = !!sym(x_var), y = avg_CPM, fill = !!sym(group_column))
  ) +
    geom_col(position = "stack") +
    
    scale_fill_manual(
      values = custom_palette,
      name = "ARG Class",
      labels = function(x) stringr::str_wrap(x, width = 18)
    ) +
    
    labs(
      y = "Average CPM",
      x = NULL
    ) +
    
    theme_minimal(base_family = "Arial") +
    theme(
      text = element_text(family = "Arial", size = text_size),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = text_size),
      axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 0, size = text_size),
      axis.text.y  = element_text(size = text_size),
      legend.title = element_text(size = text_size),
      legend.text  = element_text(size = text_size),
      legend.key.size = unit(3, "mm"),
      strip.text = element_text(size = text_size, margin = margin(t = 6, b = 6, l = 16, r = 16)),
      panel.grid = element_blank(),
      
      panel.background = element_rect(fill = NA, color = "grey70", linewidth = 0.8),
      panel.spacing.x = unit(2, "mm"),
      
      legend.position = "right",
      strip.background = element_rect(fill = "grey90")
    ) +
  
    
    facet_grid(
      ~ Group,
      scales = "free_x",
      space = "free_x",
      labeller = labeller(Group = function(x) stringr::str_wrap(x, width = 10))
    )
  
  list(
    normalized_data = normalized_df,
    summary_data    = avg_data,
    plot            = p
  )
}

# ==============================
# 7. Generate final object
# ==============================
arg_results_RF <- process_microbial_data(
  args_df_RF,
  metadata,
  group_column = "Class.y",
  x_var = "Site",
  top_n = 10
)

print(arg_results_RF$plot)

# ==============================
# 8. Export
# ==============================
ggsave(
  "Plots/arg_counts.tiff",
  arg_results_RF$plot,
  width = 183,
  height = 100,
  units = "mm",
  dpi = 600
)

ggsave(
  "Plots/arg_counts.svg",
  arg_results_RF$plot,
  width = 183,
  height = 100,
  units = "mm"
)
