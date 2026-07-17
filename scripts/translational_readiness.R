# Translational readiness of RVFV vaccine platforms
# Figure 3: From empirical vaccines to precision vaccinology 

# 0. Setup
# --------
pkgs <- c(
  "dplyr",
  "tidyr",
  "ggplot2",
  "readr",
  "stringr"
)

to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)

library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(stringr)

# -----------
# 1. Vaccine readiness matrix
# Scores:
# 1 = Low
# 2 = Moderate
# 3 = High

vaccine_data <- tribble(
  ~Platform, ~Category, ~Safety, ~Immunogenicity,
  ~Field_deployability, ~Cold_chain_feasibility,
  ~Human_readiness, ~Livestock_readiness,
  ~Africa_manufacturing_readiness, ~Evidence_level,
  
  "Smithburn", "Classical empirical vaccines",
  1, 3, 2, 2, 1, 3, 2, 3,
  
  "Clone 13", "Classical empirical vaccines",
  2, 3, 2, 2, 1, 3, 2, 3,
  
  "MP-12", "Classical empirical vaccines",
  2, 3, 2, 2, 2, 2, 1, 3,
  
  "RVFV-4s", "Rationally engineered vaccines",
  3, 2, 2, 2, 2, 2, 1, 2,
  
  "DDvax", "Rationally engineered vaccines",
  3, 2, 2, 2, 2, 2, 1, 2,
  
  "Viral-vectored vaccines", "Advanced platform technologies",
  3, 2, 2, 2, 2, 2, 1, 2,
  
  "mRNA vaccines", "Advanced platform technologies",
  3, 2, 1, 1, 2, 1, 1, 2,
  
  "Subunit / epitope-based vaccines",
  "Precision vaccinology approaches",
  3, 2, 2, 2, 1, 2, 1, 2
)

# -----------------------
# 2. Convert to long format

vaccine_long <- vaccine_data %>%
  pivot_longer(
    cols = c(
      Safety,
      Immunogenicity,
      Field_deployability,
      Cold_chain_feasibility,
      Human_readiness,
      Livestock_readiness,
      Africa_manufacturing_readiness,
      Evidence_level
    ),
    names_to = "Criterion",
    values_to = "Score"
  ) %>%
  mutate(
    
    Criterion = recode(
      Criterion,
      "Safety" = "Safety",
      "Immunogenicity" = "Immunogenicity",
      "Field_deployability" = "Field deployability",
      "Cold_chain_feasibility" = "Cold-chain feasibility",
      "Human_readiness" = "Human-use readiness",
      "Livestock_readiness" = "Livestock-use readiness",
      "Africa_manufacturing_readiness" = "African manufacturing readiness",
      "Evidence_level" = "Evidence level"
    ),
    
    Readiness = case_when(
      Score == 1 ~ "Low",
      Score == 2 ~ "Moderate",
      Score == 3 ~ "High"
    ),
    
    Platform = factor(
      Platform,
      levels = c(
        "Smithburn",
        "Clone 13",
        "MP-12",
        "RVFV-4s",
        "DDvax",
        "Viral-vectored vaccines",
        "mRNA vaccines",
        "Subunit / epitope-based vaccines"
      )
    ),
    
    Criterion = factor(
      Criterion,
      levels = c(
        "Safety",
        "Immunogenicity",
        "Field deployability",
        "Cold-chain feasibility",
        "Human-use readiness",
        "Livestock-use readiness",
        "African manufacturing readiness",
        "Evidence level"
      )
    ),
    
    Category = factor(
      Category,
      levels = c(
        "Classical empirical vaccines",
        "Rationally engineered vaccines",
        "Advanced platform technologies",
        "Precision vaccinology approaches"
      )
    )
  )

# ----------------------
# 3. Define colors

readiness_colors <- c(
  "Low" = "#D73027",
  "Moderate" = "#FDAE61",
  "High" = "#1A9850"
)

category_colors <- c(
  "Classical empirical vaccines" = "#FEE8C8",
  "Rationally engineered vaccines" = "#FDBB84",
  "Advanced platform technologies" = "#B3CDE3",
  "Precision vaccinology approaches" = "#CCEBC5"
)

# ---------------------
# 4. Plot

p_vaccine <- ggplot(
  vaccine_long,
  aes(x = Criterion, y = Platform)
) +
  geom_tile(
    aes(fill = Category),
    alpha = 0.12,
    color = NA
  ) +
  scale_fill_manual(
    values = category_colors,
    guide = "none"
  ) +
  geom_point(
    aes(size = Score, color = Readiness),
    shape = 16
  ) +
  scale_color_manual(
    values = readiness_colors,
    breaks = c("High", "Moderate", "Low"),
    name = "Readiness"
  ) +
  scale_size_continuous(
    range = c(3, 8),
    breaks = c(1, 2, 3),
    labels = c("Low", "Moderate", "High"),
    name = "Score"
  ) +
  scale_x_discrete(
    labels = function(x) stringr::str_wrap(x, width = 18)
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  guides(
    color = guide_legend(
      override.aes = list(size = 7)
    ),
    size = guide_legend(
      override.aes = list(size = c(4, 7, 10))
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(
  angle = 35,
  hjust = 1,
  vjust = 1,
  size = 14,
  face = "bold"
),

axis.text.y = element_text(
  size = 14,
  face = "bold"
),
    panel.grid.major = element_line(
      color = "grey88",
      linewidth = 0.3
    ),
    panel.grid.minor = element_blank(),
    
    legend.position = "right",
    
    legend.title = element_text(
      size = 14,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 12
    ),
    
    legend.key.size = unit(
      1.2,
      "cm"
    ),
    
    legend.spacing.y = unit(
      0.3,
      "cm"
    )
  )
# ---------------
# 5. Print figure


print(p_vaccine)

# ---------------
# 6. Save outputs

ggsave(
  filename = "Figure3_precision_vaccinology_landscape.png",
  plot = p_vaccine,
  width = 14,
  height = 8,
  dpi = 600
)

# Optional export table
write_csv(
  vaccine_data,
  "RVFV_vaccine_readiness_precision_vaccinology.csv"
)
