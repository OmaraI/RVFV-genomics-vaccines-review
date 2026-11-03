# FIGURE 1: RVFV GENOMIC DATA AVAILABILITY IN AFRICA
#
# Panel 1A:
# Publicly available RVFV sequences by country
# Country names and sequence counts in attached boxes
#
# Panel 1B:
# Percentage of isolates represented by L, M and S

# 0. LOAD REQUIRED PACKAGES
# -------------------------

pkgs <- c(
  "dplyr",
  "readr",
  "sf",
  "ggplot2",
  "rnaturalearth",
  "rnaturalearthdata",
  "stringr",
  "tidyr",
  "ggrepel",
  "patchwork",
  "tibble",
  "scales"
)

to_install <- pkgs[
  !vapply(
    pkgs,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(to_install) > 0) {
  install.packages(to_install)
}

invisible(
  lapply(pkgs, library, character.only = TRUE)
)

sf_use_s2(FALSE)

setwd(
  "/Users/IsaacOmara/Desktop/Manuscripts/Review_Article/Africa_Map/"
)



# 1. READ RVFV METADATA
# ---------------------

meta <- read_csv(
  "RVFV_clean_unique_isolate_segment_dataset.csv",
  show_col_types = FALSE
)

print(names(meta))
print(dim(meta))

required_cols <- c(
  "Country",
  "Segment",
  "Accession",
  "Nuc_Completeness"
)

missing_cols <- setdiff(
  required_cols,
  names(meta)
)

if (length(missing_cols) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

for (col in c(
  "Host",
  "Species",
  "Collection_Date"
)) {
  if (!col %in% names(meta)) {
    meta[[col]] <- NA_character_
  }
}


# 2. HARMONIZE COUNTRY NAMES
# --------------------------

harmonize_country <- function(x) {
  
  x <- stringr::str_squish(
    as.character(x)
  )
  
  dplyr::case_when(
    
    x %in% c(
      "DRC",
      "DR Congo",
      "Dem. Rep. Congo",
      "Democratic Republic of Congo"
    ) ~ "Democratic Republic of the Congo",
    
    x %in% c(
      "Congo",
      "Republic of Congo"
    ) ~ "Republic of the Congo",
    
    x == "Central African Rep." ~
      "Central African Republic",
    
    x %in% c(
      "Ivory Coast",
      "Cote d'Ivoire"
    ) ~ "Côte d'Ivoire",
    
    x == "United Republic of Tanzania" ~
      "Tanzania",
    
    x == "The Gambia" ~ "Gambia",
    
    x == "Swaziland" ~ "Eswatini",
    
    x == "Comoro Is." ~ "Comoros",
    
    TRUE ~ x
  )
}



# 3. CLEAN METADATA
# -----------------

meta <- meta %>%
  mutate(
    
    Country = harmonize_country(Country),
    
    Host = str_trim(
      as.character(Host)
    ),
    
    Species = str_trim(
      as.character(Species)
    ),
    
    Segment = str_to_upper(
      str_trim(as.character(Segment))
    ),
    
    Collection_Date = str_trim(
      as.character(Collection_Date)
    ),
    
    Accession = str_trim(
      as.character(Accession)
    ),
    
    Nuc_Completeness = str_to_lower(
      str_trim(as.character(Nuc_Completeness))
    ),
    
    Country = if_else(
      is.na(Country) | Country == "",
      "unknown_country",
      Country
    ),
    
    Host = if_else(
      is.na(Host) | Host == "",
      "unknown_host",
      Host
    ),
    
    Species = if_else(
      is.na(Species) | Species == "",
      "unknown_species",
      Species
    ),
    
    Collection_Date = if_else(
      is.na(Collection_Date) |
        Collection_Date == "",
      "unknown_date",
      Collection_Date
    )
  ) %>%
  filter(
    Segment %in% c("L", "M", "S")
  )



# 4. CREATE ISOLATE IDENTIFIER
# ----------------------------

if ("isolate_key" %in% names(meta)) {
  
  meta <- meta %>%
    mutate(
      isolate_id = as.character(isolate_key)
    )
  
} else if ("Isolate" %in% names(meta)) {
  
  meta <- meta %>%
    mutate(
      isolate_id = if_else(
        
        !is.na(Isolate) &
          str_trim(as.character(Isolate)) != "",
        
        paste(
          Country,
          str_trim(as.character(Isolate)),
          sep = " | "
        ),
        
        paste0(
          "ACCESSION_",
          Accession
        )
      )
    )
  
} else {
  
  meta <- meta %>%
    mutate(
      isolate_id = paste0(
        "ACCESSION_",
        Accession
      )
    )
}

meta <- meta %>%
  mutate(
    isolate_id = if_else(
      
      is.na(isolate_id) |
        str_trim(isolate_id) == "",
      
      paste0(
        "ACCESSION_",
        Accession
      ),
      
      isolate_id
    )
  )



# 5. COUNT PUBLICLY AVAILABLE RVFV SEQUENCES
# ------------------------------------------

# Counts represent distinct accession records.
# L, M and S sequences are counted separately.
# They do not directly measure surveillance capacity.

sequence_counts <- meta %>%
  filter(
    Country != "unknown_country",
    !is.na(Accession),
    Accession != ""
  ) %>%
  group_by(Country) %>%
  summarise(
    
    n_sequences = n_distinct(Accession),
    
    n_L = n_distinct(
      Accession[Segment == "L"]
    ),
    
    n_M = n_distinct(
      Accession[Segment == "M"]
    ),
    
    n_S = n_distinct(
      Accession[Segment == "S"]
    ),
    
    .groups = "drop"
  ) %>%
  rename(
    admin = Country
  )

print(sequence_counts)



# 6. COLLAPSE RECORDS TO ISOLATE LEVEL
# -----------------------------------

isolates <- meta %>%
  group_by(
    isolate_id,
    Country,
    Species,
    Host,
    Collection_Date
  ) %>%
  summarise(
    
    n_accessions = n(),
    
    segments_present = paste(
      sort(unique(Segment)),
      collapse = ","
    ),
    
    n_segments = n_distinct(Segment),
    
    has_L = any(Segment == "L"),
    
    has_M = any(Segment == "M"),
    
    has_S = any(Segment == "S"),
    
    complete_records = sum(
      Nuc_Completeness == "complete",
      na.rm = TRUE
    ),
    
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



# 7. COUNTRY-LEVEL ISOLATE SUMMARY
# --------------------------------

country_summary <- isolates %>%
  filter(
    Country != "unknown_country"
  ) %>%
  group_by(Country) %>%
  summarise(
    
    n_isolates = n(),
    
    n_full_tripartite = sum(
      segment_pattern == "L+M+S"
    ),
    
    n_partial = sum(
      segment_pattern != "L+M+S"
    ),
    
    pct_full_tripartite = round(
      100 * n_full_tripartite / n_isolates,
      1
    ),
    
    .groups = "drop"
  ) %>%
  rename(
    admin = Country
  )

print(country_summary)


# 8. REPORTED RVF COUNTRIES WITHOUT SEQUENCE RECORDS
# --------------------------------------------------

reported_rvf_countries <- c(
  "South Sudan",
  "Comoros",
  "Gambia"
)

genomic_countries <- sequence_counts %>%
  filter(n_sequences > 0) %>%
  pull(admin)

reported_no_genomes <- setdiff(
  reported_rvf_countries,
  genomic_countries
)



# 9. LOAD AFRICA MAP
# ------------------

africa <- ne_countries(
  continent = "Africa",
  returnclass = "sf"
) %>%
  mutate(
    admin = harmonize_country(admin)
  )

unmatched_countries <- sequence_counts %>%
  anti_join(
    africa %>%
      st_drop_geometry() %>%
      select(admin),
    by = "admin"
  )

if (nrow(unmatched_countries) > 0) {
  
  print(unmatched_countries)
  
  stop(
    "Unmatched country names detected. ",
    "Correct the names before plotting."
  )
}



# 10. JOIN SEQUENCE DATA TO MAP
# -----------------------------

africa_gap <- africa %>%
  left_join(
    country_summary,
    by = "admin"
  ) %>%
  left_join(
    sequence_counts,
    by = "admin"
  ) %>%
  mutate(
    
    n_isolates = replace_na(
      n_isolates, 0L
    ),
    
    n_sequences = replace_na(
      n_sequences, 0L
    ),
    
    n_full_tripartite = replace_na(
      n_full_tripartite, 0L
    ),
    
    n_partial = replace_na(
      n_partial, 0L
    ),
    
    pct_full_tripartite = if_else(
      n_isolates == 0,
      NA_real_,
      pct_full_tripartite
    ),
    
    reported_rvf = if_else(
      admin %in% reported_no_genomes,
      "Yes",
      "No"
    )
  )


# 11. FIGURE 1A SEQUENCE-COUNT CATEGORIES
# ---------------------------------------

africa_gap <- africa_gap %>%
  mutate(
    
    gap_class = case_when(
      
      n_sequences >= 15 ~
        "15 or more sequences",
      
      n_sequences >= 5 ~
        "5-14 sequences",
      
      n_sequences >= 1 ~
        "1-4 sequences",
      
      reported_rvf == "Yes" ~
        "Reported RVF\n(no public sequences)",
      
      TRUE ~
        "No public sequences identified"
    ),
    
    gap_class = factor(
      gap_class,
      
      levels = c(
        "15 or more sequences",
        "5-14 sequences",
        "1-4 sequences",
        "Reported RVF\n(no public sequences)",
        "No public sequences identified"
      )
    )
  )



# 12. CREATE COUNTRY LABEL COORDINATES
# ------------------------------------

genomic_label_df <- africa_gap %>%
  filter(
    n_sequences > 0,
    admin != "Comoros"
  ) %>%
  st_point_on_surface()

label_xy <- st_coordinates(
  genomic_label_df
)

genomic_label_coords <- genomic_label_df %>%
  st_drop_geometry() %>%
  mutate(
    
    X = label_xy[, 1],
    
    Y = label_xy[, 2],
    
    label_country = admin,
    
    label_count = as.character(
      n_sequences
    ),
    
    label_B = admin
  )

stopifnot(
  exists("genomic_label_coords"),
  nrow(genomic_label_coords) > 0
)


# 13. MANUAL COUNTRY LABEL POSITIONS
# ----------------------------------

# Adjust these coordinates if larger boxes overlap.
# X_lab controls horizontal position.
# Y_lab controls vertical position.


manual_positions <- tibble(
  
  admin = c(
    "Uganda",
    "Kenya",
    "Rwanda",
    "Burundi",
    "Senegal",
    "Guinea",
    "Burkina Faso"
  ),
  
  X_lab = c(
    31.0,    # Uganda: move west (left)
    41.5,    # Kenya: move east (right)
    27.5,
    27.5,
    -19.5,
    -14.0,
    -1.5
  ),
  
  Y_lab = c(
    4.5,     # Uganda: move north (up)
    -1.0,    # Kenya
    -1.5,
    -5.0,
    10.0,
    8.0,
    14.0
  )
)

# 14. PREPARE ATTACHED BOX LABELS FOR FIGURE 1A
# ---------------------------------------------

boxed_labels <- genomic_label_coords %>%
  left_join(
    manual_positions,
    by = "admin"
  ) %>%
  mutate(
    
    X_lab = coalesce(
      X_lab, X
    ),
    
    Y_lab = coalesce(
      Y_lab, Y
    )
  )

print(
  boxed_labels %>%
    select(
      admin,
      n_sequences,
      X_lab,
      Y_lab
    )
)



# 15. PREPARE FIGURE 1B LABELS
# ----------------------------

manual_label_coords <- genomic_label_coords %>%
  inner_join(
    manual_positions,
    by = "admin"
  )

auto_label_coords <- genomic_label_coords %>%
  filter(
    !admin %in% manual_positions$admin
  )



# 16. REPORTED-RVF COUNTRY ANNOTATIONS
# ------------------------------------

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
) %>%
  filter(
    admin %in% reported_no_genomes
  )


# 17. ORIGINAL FIGURE 1A COLOURS
# ------------------------------

gap_colors <- c(
  
  "15 or more sequences" = "#238B45",
  
  "5-14 sequences" = "#66BD63",
  
  "1-4 sequences" = "#D9F0D3",
  
  "Reported RVF\n(no public sequences)" = "#A6611A",
  
  "No public sequences identified" = "grey96"
)


# ----------------------------------------------------------
# 18. MAP THEME
# ----------------------------------------------------------

map_theme <- theme_minimal(
  base_size = 11
) +
  theme(
    
    axis.title = element_blank(),
    
    axis.text = element_blank(),
    
    panel.grid = element_blank(),
    
    plot.title = element_text(
      face = "bold",
      size = 11
    ),
    
    plot.margin = margin(
      8, 20, 8, 8
    )
  )


# 19. CUSTOM ATTACHED LABEL GEOMETRY
# =================================

# This section controls the size of the country
# rectangles and sequence-number boxes.
#
# Country name: 10 pt bold
# Country rectangle: 6 mm high
# Sequence count: 10 pt bold
# Count rectangle: 5 mm high
#
# Both boxes are attached at the same horizontal edge.


# 19A. ADJUSTABLE LABEL SIZES
# ---------------------------

country_font_size <- 10

count_font_size <- 10

country_box_height <- 6.0

count_box_height <- 5.0

country_horizontal_padding <- 4.0

count_horizontal_padding <- 3.0



# 19B. CREATE CUSTOM GGPLOT GEOMETRY
# ----------------------------------

GeomAttachedLabel <- ggplot2::ggproto(
  
  "GeomAttachedLabel",
  
  ggplot2::Geom,
  
  required_aes = c(
    "x",
    "y",
    "country",
    "count"
  ),
  
  default_aes = ggplot2::aes(
    colour = "black",
    fill = "white",
    linewidth = 0.3,
    alpha = 1
  ),
  
  draw_key = ggplot2::draw_key_blank,
  
  draw_panel = function(
    data,
    panel_params,
    coord
  ) {
    
    coords <- coord$transform(
      data,
      panel_params
    )
    
    grobs <- lapply(
      seq_len(nrow(coords)),
      
      function(i) {
        
        country <- as.character(
          coords$country[i]
        )
        
        count <- as.character(
          coords$count[i]
        )
        
        x <- coords$x[i]
        y <- coords$y[i]
        
      
        # COUNTRY NAME TEXT AND BOX DIMENSIONS
        # ----------------------------------------------
        
        country_text <- grid::textGrob(
          country,
          
          gp = grid::gpar(
            fontsize = country_font_size,
            fontface = "bold",
            col = "black"
          )
        )
        
        country_width <-
          grid::grobWidth(country_text) +
          grid::unit(
            country_horizontal_padding,
            "mm"
          )
        
        country_height <- grid::unit(
          country_box_height,
          "mm"
        )
        
        
        # SEQUENCE COUNT TEXT AND BOX DIMENSIONS
        # --------------------------------------
        
        count_text <- grid::textGrob(
          count,
          
          gp = grid::gpar(
            fontsize = count_font_size,
            fontface = "bold",
            col = "black"
          )
        )
        
        count_width <-
          grid::grobWidth(count_text) +
          grid::unit(
            count_horizontal_padding,
            "mm"
          )
        
        count_height <- grid::unit(
          count_box_height,
          "mm"
        )
        
        
       
        # ATTACH BOTH RECTANGLES
        # ----------------------
        
        # The lower edge of the country box
        # and upper edge of the count box
        # meet exactly at the same anchor point.
        
        country_y <- grid::unit(
          y, "npc"
        ) + country_height / 2
        
        count_y <- grid::unit(
          y, "npc"
        ) - count_height / 2
        
        
        
        # DRAW BOTH BOXES AND TEXT
        # ------------------------
        
        grid::grobTree(
          
          # Country-name rectangle
          grid::rectGrob(
            
            x = grid::unit(x, "npc"),
            
            y = country_y,
            
            width = country_width,
            
            height = country_height,
            
            gp = grid::gpar(
              fill = "white",
              col = "black",
              lwd = 0.7
            )
          ),
          
          # Country name
          grid::textGrob(
            
            country,
            
            x = grid::unit(x, "npc"),
            
            y = country_y,
            
            gp = grid::gpar(
              fontsize = country_font_size,
              fontface = "bold",
              col = "black"
            )
          ),
          
          # Attached sequence-count rectangle
          grid::rectGrob(
            
            x = grid::unit(x, "npc"),
            
            y = count_y,
            
            width = count_width,
            
            height = count_height,
            
            gp = grid::gpar(
              fill = "white",
              col = "black",
              lwd = 0.7
            )
          ),
          
          # Sequence count
          grid::textGrob(
            
            count,
            
            x = grid::unit(x, "npc"),
            
            y = count_y,
            
            gp = grid::gpar(
              fontsize = count_font_size,
              fontface = "bold",
              col = "black"
            )
          )
        )
      }
    )
    
    do.call(
      grid::grobTree,
      grobs
    )
  }
)



# 19C. CREATE ATTACHED-LABEL FUNCTION
# -----------------------------------

geom_attached_label <- function(
    mapping = NULL,
    data = NULL,
    ...
) {
  
  ggplot2::layer(
    
    geom = GeomAttachedLabel,
    
    stat = "identity",
    
    position = "identity",
    
    mapping = mapping,
    
    data = data,
    
    inherit.aes = FALSE,
    
    params = list(...)
  )
}



# 20. FIGURE 1A: PUBLIC SEQUENCE AVAILABILITY
# ===========================================

p_gap <- ggplot(africa_gap) +
  
  geom_sf(
    aes(fill = gap_class),
    
    color = "grey15",
    
    linewidth = 0.45
  )



# 20A. CONNECTOR LINES FOR DISPLACED LABELS
# -----------------------------------------

p_gap <- p_gap +
  
  geom_segment(
    
    data = boxed_labels %>%
      filter(
        abs(X_lab - X) > 0.01 |
          abs(Y_lab - Y) > 0.01
      ),
    
    aes(
      x = X,
      y = Y,
      xend = X_lab,
      yend = Y_lab
    ),
    
    inherit.aes = FALSE,
    
    color = "grey35",
    
    linewidth = 0.3
  )

# 20B. ADD ENLARGED ATTACHED COUNTRY LABELS
# -----------------------------------------

p_gap <- p_gap +
  
  geom_attached_label(
    
    data = boxed_labels,
    
    mapping = aes(
      x = X_lab,
      y = Y_lab,
      country = label_country,
      count = label_count
    )
  )



# 20C. REPORTED-RVF COUNTRY ANNOTATIONS
# -------------------------------------

p_gap <- p_gap +
  
  geom_segment(
    
    data = reported_label_coords,
    
    aes(
      x = X,
      y = Y,
      xend = X_lab,
      yend = Y_lab
    ),
    
    inherit.aes = FALSE,
    
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
    
    inherit.aes = FALSE,
    
    size = 2.2,
    
    fontface = "bold"
  )


# 20D. COMOROS ANNOTATION
# -----------------------

if ("Comoros" %in% reported_no_genomes) {
  
  p_gap <- p_gap +
    
    annotate(
      "point",
      
      x = 43.3,
      y = -11.7,
      
      shape = 21,
      
      fill = "#A6611A",
      
      colour = "black",
      
      size = 3.2,
      
      stroke = 0.3
    )
  
}

p_gap <- p_gap +
  
  annotate(
    "text",
    
    x = 43.3,
    y = -10.5,
    
    label = "Comoros",
    
    size = 2.4,
    
    fontface = "bold",
    
    hjust = 0.5,
    
    vjust = 0
  )



# 20E. FIGURE 1A COLOURS AND LEGEND
# ---------------------------------

p_gap <- p_gap +
  
  scale_fill_manual(
    
    values = gap_colors,
    
    breaks = c(
      "15 or more sequences",
      "5-14 sequences",
      "1-4 sequences",
      "Reported RVF\n(no public sequences)",
      "No public sequences identified"
    ),
    
    name = "Publicly available\nRVFV sequences",
    
    drop = FALSE
  ) +
  
  labs(
    title = "(1A)."
  ) +
  
  coord_sf(
    
    xlim = c(-20, 55),
    
    ylim = c(-38, 38),
    
    clip = "off"
  ) +
  
  map_theme +
  
  theme(
    
    # Original Figure 1A legend position
    legend.position = c(0.14, 0.20),
    
    legend.justification = c(0, 0),
    
    legend.background = element_rect(
      fill = "white",
      colour = NA
    ),
    
    legend.key.height = grid::unit(
      0.45, "cm"
    ),
    
    legend.key.width = grid::unit(
      0.45, "cm"
    ),
    
    legend.title = element_text(
      size = 13,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 9
    )
  )


# 21. FIGURE 1B: TRIPARTITE GENOME REPRESENTATION
#
# Country names in rectangular boxes
# No connector lines
# No sequence numbers
# Original green gradient retained



# 21A. PREPARE COUNTRY LABEL POSITIONS
# ------------------------------------

# Reuse manually adjusted positions from Figure 1A,
# including Uganda and Kenya.

boxed_labels_B <- genomic_label_coords %>%
  left_join(
    manual_positions,
    by = "admin"
  ) %>%
  mutate(
    X_lab = coalesce(X_lab, X),
    Y_lab = coalesce(Y_lab, Y)
  )


# 21B. CREATE FIGURE 1B MAP
# -------------------------

p_tripartite <- ggplot(africa_gap) +
  
  geom_sf(
    aes(fill = pct_full_tripartite),
    color = "grey15",
    linewidth = 0.45
  )


# 21C. ADD COUNTRY-NAME RECTANGLES ONLY
# -------------------------------------

# Uses the same 10-point font as Figure 1A.
# Rectangular white boxes with black borders.
# No connector lines or sequence counts.

p_tripartite <- p_tripartite +
  
  geom_label(
    data = boxed_labels_B,
    
    aes(
      x = X_lab,
      y = Y_lab,
      label = admin
    ),
    
    inherit.aes = FALSE,
    
    size = 10 / ggplot2::.pt,
    
    fontface = "bold",
    
    fill = "white",
    
    color = "black",
    
    label.size = 0.25,
    
    label.padding = grid::unit(
      0.20, "lines"
    ),
    
    label.r = grid::unit(
      0, "lines"
    )
  )

# 21D. COMOROS ANNOTATION
# -----------------------

# Keep Comoros visible as a labelled island.
# No connector line.

p_tripartite <- p_tripartite +
  
  annotate(
    "point",
    x = 43.3,
    y = -11.7,
    shape = 21,
    fill = "white",
    colour = "black",
    size = 3.2,
    stroke = 0.5
  ) +
  
  geom_label(
    data = data.frame(
      X = 43.3,
      Y = -10.5,
      label = "Comoros"
    ),
    
    aes(
      x = X,
      y = Y,
      label = label
    ),
    
    inherit.aes = FALSE,
    
    size = 10 / ggplot2::.pt,
    
    fontface = "bold",
    
    fill = "white",
    
    color = "black",
    
    label.size = 0.25,
    
    label.padding = grid::unit(
      0.20, "lines"
    ),
    
    label.r = grid::unit(
      0, "lines"
    )
  )


# 21E. IMPROVED FIGURE 1B LEGEND
# ------------------------------

# Green gradient: 0–100%
# White square: No isolate data available
# White square directly below the colour bar
# ==========================================

library(ggplot2)
library(cowplot)
library(grid)

# 1. APPLY THE ORIGINAL GREEN GRADIENT
# ------------------------------------

p_tripartite <- p_tripartite +
  
  scale_fill_gradientn(
    colours = c(
      "#D9F0D3",
      "#A1D99B",
      "#74C476",
      "#31A354",
      "#006D2C"
    ),
    
    values = seq(0, 1, length.out = 5),
    
    limits = c(0, 100),
    
    breaks = seq(0, 100, by = 20),
    
    labels = function(x) paste0(x, "%"),
    
    na.value = "white",
    
    name = "Isolates with all three\nRVFV segments (%)",
    
    guide = guide_colorbar(
      order = 1,
      barheight = unit(5, "cm"),
      barwidth = unit(0.6, "cm"),
      ticks = TRUE,
      frame.colour = "grey30"
    )
  ) +
  
  labs(title = "(1B).") +
  
  coord_sf(
    xlim = c(-20, 55),
    ylim = c(-38, 38),
    clip = "off"
  ) +
  
  map_theme +
  
  theme(
    legend.position = "right",
    
    legend.title = element_text(
      face = "bold",
      size = 11
    ),
    
    legend.text = element_text(
      size = 10
    ),
    
    legend.margin = margin(0, 0, 0, 0)
  )

# 2. EXTRACT THE GREEN GRADIENT LEGEND
# ------------------------------------

green_legend <- cowplot::get_legend(
  p_tripartite
)


# 3. CREATE THE WHITE-SQUARE LEGEND
# ---------------------------------

missing_data <- data.frame(
  status = "No isolate data available",
  x = 1,
  y = 1
)

p_missing <- ggplot(
  missing_data,
  aes(x = x, y = y, shape = status)
) +
  
  geom_point(
    size = 4,
    fill = "white",
    colour = "black",
    stroke = 0.5
  ) +
  
  scale_shape_manual(
    values = c(
      "No isolate data available" = 22
    ),
    name = NULL
  ) +
  
  guides(
    shape = guide_legend(
      override.aes = list(
        shape = 22,
        fill = "white",
        colour = "black",
        size = 4
      )
    )
  ) +
  
  theme_void() +
  
  theme(
    legend.position = "right",
    legend.justification = "left",
    legend.key.width = unit(0.6, "cm"),
    legend.key.height = unit(0.6, "cm"),
    legend.text = element_text(size = 10),
    legend.margin = margin(0, 0, 0, 0)
  )

white_legend <- cowplot::get_legend(p_missing)


# 4. STACK THE TWO LEGENDS
# ------------------------

combined_legend <- cowplot::plot_grid(
  green_legend,
  white_legend,
  
  ncol = 1,
  
  align = "v",
  
  axis = "l",
  
  rel_heights = c(5, 1)
)


# 5. REMOVE ORIGINAL LEGEND FROM MAP
# ----------------------------------

p_tripartite_map <- p_tripartite +
  theme(legend.position = "none")


# 6. PLACE THE COMBINED LEGEND BESIDE THE MAP
# -------------------------------------------
p_tripartite_final <- cowplot::plot_grid(
  p_tripartite_map,
  combined_legend,
  
  ncol = 2,
  
  rel_widths = c(4.5, 1),
  
  align = "h"
)


# 7. DISPLAY FIGURE 1B
# --------------------

print(p_tripartite_final)


# 21F. DISPLAY FIGURE 1B
# ----------------------

print(p_tripartite)


# 22. COMBINE FIGURE 1A AND FIGURE 1B
# ===================================

combined_plot <- p_gap + p_tripartite +
  
  plot_layout(
    ncol = 2,
    widths = c(1, 1)
  )

print(combined_plot)



# 23. SAVE PUBLICATION-QUALITY FIGURES
# ====================================

ggsave(
  
  filename = "Figure1_combined_RVFV_maps.png",
  
  plot = combined_plot,
  
  width = 18,
  
  height = 8,
  
  dpi = 600,
  
  bg = "white"
)

ggsave(
  
  filename = "Figure1_combined_RVFV_maps.pdf",
  
  plot = combined_plot,
  
  width = 18,
  
  height = 8,
  
  bg = "white"
)



# 24. EXPORT SUPPORTING DATA
# ==========================

write_csv(
  isolates,
  "RVFV_collapsed_isolate_table_integrated.csv"
)

write_csv(
  country_summary,
  "RVFV_country_summary_integrated.csv"
)

write_csv(
  sequence_counts,
  "RVFV_country_public_sequence_counts.csv"
)

write_csv(
  
  africa_gap %>%
    st_drop_geometry() %>%
    select(
      admin,
      n_sequences,
      n_isolates,
      n_full_tripartite,
      pct_full_tripartite,
      gap_class
    ),
  
  "RVFV_Figure1_country_mapping_data.csv"
)

write_csv(
  unmatched_countries,
  "RVFV_unmatched_country_names.csv"
)


# =============
# END OF SCRIPT
# =============
