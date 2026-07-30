# Figure 2: RVFV genomic surveillance coverage and completeness across Africa


#   A = Distribution of publicly available RVFV genomic surveillance relative 
       #to reported RVF occurrence
#   B = Percentage of isolates represented by complete L, M and S genomes
# Input file: sequences.csv

# 0. Setup
# --------
pkgs <- c(
  "dplyr", "readr", "sf", "ggplot2", "rnaturalearth",
  "rnaturalearthdata", "stringr", "tidyr", "ggrepel", "patchwork", "ggnewscale"
)



to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)

library(dplyr)
library(readr)
library(sf)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(stringr)
library(tidyr)
library(ggrepel)
library(patchwork)
library(ggnewscale)

sf_use_s2(FALSE)

getwd()
setwd("/Users/IsaacOmara/Desktop/Manuscripts/Review_Article/Africa_Map/")

# 1. Read metadata
# ----------------
meta <- read_csv("RVFV_clean_unique_isolate_segment_dataset.csv", show_col_types = FALSE) # downloaded from ncbi

print(names(meta))
print(dim(meta))


# 2. Clean metadata
# -----------------
meta <- meta %>%
  mutate(
    Country = str_trim(Country),
    Host = str_trim(Host),
    Species = str_trim(Species),
    Segment = str_trim(Segment),
    Collection_Date = str_trim(Collection_Date),
    Country = ifelse(is.na(Country) | Country == "", "unknown_country", Country),
    Host = ifelse(is.na(Host) | Host == "", "unknown_host", Host),
    Species = ifelse(is.na(Species) | Species == "", "unknown_species", Species),
    Collection_Date = ifelse(is.na(Collection_Date) | Collection_Date == "", "unknown_date", Collection_Date)
  ) %>%
  filter(Segment %in% c("L", "M", "S"))


# 3. Create isolate-level identifier
# 3. Use the cleaned isolate identifier
# ----------

meta <- meta %>%
  mutate(
    isolate_id = case_when(
      "isolate_key" %in% names(meta) ~ isolate_key,
      !is.na(Isolate) & Isolate != "" ~ paste(Country, Isolate, sep = " | "),
      TRUE ~ paste0("ACCESSION_", Accession)
    )
  )

# 4. Collapse to isolate level across all segments
# ------------------------------------------------
isolates <- meta %>%
  group_by(isolate_id, Country, Species, Host, Collection_Date) %>%
  summarise(
    n_accessions = n(),
    segments_present = paste(sort(unique(Segment)), collapse = ","),
    n_segments = n_distinct(Segment),
    has_L = any(Segment == "L"),
    has_M = any(Segment == "M"),
    has_S = any(Segment == "S"),
    complete_records = sum(Nuc_Completeness == "complete", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    segment_pattern = case_when(
      has_L & has_M & has_S ~ "L+M+S",
      has_L & has_M & !has_S ~ "L+M",
      has_L & !has_M & has_S ~ "L+S",
      !has_L & has_M & has_S ~ "M+S",
      has_L & !has_M & !has_S ~ "L only",
      !has_L & has_M & !has_S ~ "M only",
      !has_L & !has_M & has_S ~ "S only",
      TRUE ~ "Other"
    )
  )


# 5. Country-level isolate summary
# --------------------------------
country_summary <- isolates %>%
  filter(Country != "unknown_country") %>%
  group_by(Country) %>%
  summarise(
    n_isolates = n(),
    n_full_tripartite = sum(segment_pattern == "L+M+S"),
    n_partial = sum(segment_pattern != "L+M+S"),
    pct_full_tripartite = round(100 * n_full_tripartite / n_isolates, 1),
    .groups = "drop"
  ) %>%
  rename(admin = Country)

print(country_summary)


# 6. Countries with genomic data
# Automatically derived from the cleaned dataset
# ----------------------------------------------
genomic_countries <- country_summary %>%
  filter(n_isolates > 0) %>%
  pull(admin)

genomic_df <- country_summary %>%
  filter(n_isolates > 0) %>%
  transmute(
    admin,
    genomic_data = "Yes"
  )

print(genomic_countries)



# 7. Reported RVF occurrence but no eligible genomic data
# ----------------------------------------------------------

reported_rvf_countries <- c(
  "South Sudan",
  "Comoros",
  "Gambia"
)

reported_no_genomes <- setdiff(
  reported_rvf_countries,
  genomic_countries
)

reported_df <- tibble(
  admin = reported_no_genomes,
  reported_rvf = "Yes"
)

print(reported_no_genomes)


# 8. Load Africa map and harmonize names
# --------------------------------------
africa <- ne_countries(continent = "Africa", returnclass = "sf")

africa <- africa %>%
  mutate(
    admin = case_when(
      admin == "Central African Rep." ~ "Central African Republic",
      admin == "Dem. Rep. Congo" ~ "Democratic Republic of the Congo",
      admin == "Congo" ~ "Republic of the Congo",
      admin == "Ivory Coast" ~ "Côte d'Ivoire",
      admin == "United Republic of Tanzania" ~ "Tanzania",
      
      admin %in% c("The Gambia", "Gambia") ~ "Gambia",
      admin %in% c("Comoro Is.", "Comoros") ~ "Comoros",
      
      TRUE ~ admin
    )
  )


# 9. Join map, genomic metadata and outbreak-only countries
# ----------------------------------------------------------

reported_no_genomes <- c(
  "South Sudan",
  "Gambia",
  "Comoros"
)

africa_gap <- africa %>%
  left_join(country_summary, by = "admin") %>%
  mutate(
    
    n_isolates = ifelse(is.na(n_isolates), 0, n_isolates),
    
    n_full_tripartite = ifelse(is.na(n_full_tripartite), 0, n_full_tripartite),
    
    n_partial = ifelse(is.na(n_partial), 0, n_partial),
    
    pct_full_tripartite = ifelse(
      n_isolates == 0,
      NA_real_,
      pct_full_tripartite
    ),
    
    genomic_data = ifelse(
      n_isolates > 0,
      "Yes",
      "No"
    ),
    
    reported_rvf = ifelse(
      admin %in% reported_no_genomes,
      "Yes",
      "No"
    ),
    
    no_genome_outbreak = ifelse(
      n_isolates == 0 & reported_rvf == "No",
      "No public genomes /\nno reported outbreaks",
      NA_character_
    )
  )

# 10. Define final map categories
# -------------------------------
africa_gap <- africa_gap %>%
  mutate(
    gap_class = case_when(
      genomic_data == "Yes" & n_isolates >= 15 ~
        "Strong genomic coverage",
      
      genomic_data == "Yes" &
        n_isolates >= 5 &
        n_isolates < 15 ~
        "Moderate genomic coverage",
      
      genomic_data == "Yes" &
        n_isolates >= 1 &
        n_isolates < 5 ~
        "Limited genomic coverage",
      
      genomic_data == "No" &
        reported_rvf == "Yes" ~
        "Reported RVF\n(no public genomic data)",
      
      TRUE ~
        "No reported RVF occurrence identified / no public genomic data"
    )
  )


# 11. Labels for genomic-data countries
# -------------------------------------
genomic_label_df <- africa_gap %>%
  filter(
    admin %in% setdiff(genomic_countries, "Comoros")
  ) %>%
  st_point_on_surface()

genomic_label_coords <- genomic_label_df %>%
  mutate(
    X = st_coordinates(geometry)[,1],
    Y = st_coordinates(geometry)[,2],
    label = admin
  ) %>%
  st_drop_geometry()



# 12. Manual labels for reported-no-genomes countries
# ---------------------------------------------------

# ----------------------------------------------------------
# 12. Manual labels for reported-RVF countries
# Comoros is excluded because it will use a custom brown dot
# ----------------------------------------------------------

reported_label_coords <- tibble(
  admin = c(
    "South Sudan",
    "Gambia"
  ),
  
  label = c(
    "South Sudan",
    "Gambia"
  ),
  
  X = c(
    31.5,
    -15.4
  ),
  
  Y = c(
    7.3,
    13.4
  ),
  
  X_lab = c(
    38.0,
    -21.0
  ),
  
  Y_lab = c(
    10.5,
    16.5
  )
)


# 13. Plot settings
# -----------------
gap_colors <- c(

  "Strong genomic coverage" = "#238B45",

  "Moderate genomic coverage" = "#66BD63",

  "Limited genomic coverage" = "#D9F0D3", 

  "Genomic data country with no isolate count" = "#FDB863",

  "Reported RVF\n(no public genomic data)" = "#A6611A",

  "No reported RVF occurrence identified / no public genomic data" = "grey96"
)

map_theme <- theme_minimal(base_size = 11) +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 9),
    plot.margin = margin(8, 20, 8, 8)
  )


# 14A. Panel A: Surveillance gap map
# ---------------------------------
p_gap <- ggplot(africa_gap) +
  geom_sf(
    aes(fill = gap_class),
    color = "grey15",
    linewidth = 0.45
  ) +
  geom_text_repel(
    data = genomic_label_coords,
    aes(x = X, y = Y, label = label),
    size = 2.2,
    lineheight = 0.9,
    min.segment.length = 0,
    box.padding = 0.20,
    point.padding = 0.12,
    segment.color = "grey30",
    max.overlaps = Inf
  ) +
  geom_segment(
    data = reported_label_coords,
    aes(
      x = X,
      y = Y,
      xend = X_lab,
      yend = Y_lab
    ),
    color = "grey30",
    linewidth = 0.25
  ) +
  
  geom_text(
    data = reported_label_coords,
    aes(
      x = X_lab,
      y = Y_lab,
      label = label
    ),
    size = 2.2
  ) +
  
  
  
# Comoros annotation
# -----------------

# Comoros: brown dot with label directly above
geom_point(
  data = NULL,
  aes(x = 43.3, y = -11.7),
  shape = 21,
  fill = "#A6611A",
  colour = "black",
  size = 3.2,
  stroke = 0.3,
  inherit.aes = FALSE
) +
  
  geom_text(
    data = NULL,
    aes(x = 43.3, y = -10.5),
    label = "Comoros",
    size = 2.4,
    hjust = 0.5,
    vjust = 0,
    inherit.aes = FALSE
  ) +
  
  scale_fill_manual(
    values = gap_colors,
    breaks = c(
      "Strong genomic coverage",
      "Moderate genomic coverage",
      "Limited genomic coverage",
      "Genomic data country with no isolate count",
      "Reported RVF\n(no public genomic data)"
    ),
    name = "Genomic surveillance\n coverage"
  ) +
  labs(
    title = "(2A).",
    #subtitle = "Distinct isolate-level genomic records integrated across L, M, and S segments"
  ) +
  coord_sf(
    xlim = c(-20, 55),
    ylim = c(-38, 38),
    clip = "off"
  ) +
  map_theme +
  theme(
    legend.position = c(0.14, 0.20),
    legend.justification = c(0, 0),
    
    legend.background = element_rect(
      fill = "white",
      colour = NA
    ),
    
    legend.key.height = unit(0.45, "cm"),
    legend.key.width  = unit(0.45, "cm"),
    legend.title = element_text(
      size = 13,
      face = "bold"
    ),
    
    legend.text = element_text(size = 9)
  )

# 14B. Panel B: Completeness map
# ------------------------------
# 14B. Panel B: Completeness map
# Green palette matching Panel 2a
# ------------------------------

p_tripartite <- ggplot(africa_gap) +
  
  geom_sf(
    aes(fill = pct_full_tripartite),
    color = "grey15",
    linewidth = 0.45
  ) +
  
  geom_text_repel(
    data = genomic_label_coords,
    aes(x = X, y = Y, label = label),
    size = 2.2,
    lineheight = 0.9,
    min.segment.length = 0,
    box.padding = 0.20,
    point.padding = 0.12,
    segment.color = "grey30",
    max.overlaps = Inf
  ) +
  
  # Comoros marker
  geom_point(
    data = NULL,
    aes(x = 43.3, y = -11.7),
    shape = 21,
    fill = "white",
    colour = "black",
    size = 3.2,
    stroke = 0.5,
    inherit.aes = FALSE
  ) +
  
  geom_text(
    data = NULL,
    aes(x = 43.3, y = -10.5),
    label = "Comoros",
    size = 2.4,
    hjust = 0.5,
    vjust = 0,
    inherit.aes = FALSE
  ) +
  
  # Green completeness gradient
  scale_fill_gradientn(
    colours = c(
      "#D9F0D3",   # 0%
      "#A1D99B",
      "#74C476",
      "#31A354",
      "#006D2C"    # 100%
    ),
    values = scales::rescale(c(0, 20, 40, 70, 100)),
    limits = c(0, 100),
    breaks = c(0, 20, 40, 60, 80, 100),
    na.value = "white",
    name = "Tripartite genome\ncompleteness (%)"
  ) +
  
  labs(
    title = "(2B)."
  ) +
  
  coord_sf(
    xlim = c(-20, 55),
    ylim = c(-38, 38),
    clip = "off"
  ) +
  
  map_theme +
  
  theme(
    legend.title = element_text(
      face = "bold",
      size = 13
    ),
    legend.text = element_text(
      size = 10
    )
  )


# 15. Combine side by side
# ------------------------
combined_plot <- p_gap + p_tripartite +
  plot_layout(ncol = 2, widths = c(1, 1)) +
  plot_annotation(
    #title = "Figure 2. Geographic distribution and completeness of publicly available RVFV genomic surveillance across Africa",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
    )
  )

print(combined_plot)


# 16. Save outputs
# ---------------
ggsave(
  filename = "Figure2_combined_RVFV_maps.png",
  plot = combined_plot,
  width = 18,
  height = 8,
  dpi = 600
)


write_csv(isolates, "RVFV_collapsed_isolate_table_integrated.csv")
write_csv(country_summary, "RVFV_country_summary_integrated.csv")
