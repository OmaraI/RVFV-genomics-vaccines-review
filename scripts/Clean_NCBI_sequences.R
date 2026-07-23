# ==========================================================
# Clean RVFV NCBI Virus metadata
#
# Purpose:
# 1. Remove exact duplicate accession records
# 2. Exclude vaccine, laboratory-passaged, recombinant,
#    synthetic and experimentally manipulated sequences
# 3. Create a reliable isolate identifier
# 4. Retain one best sequence per isolate and genome segment
# 5. Prefer complete sequences over partial sequences
# 6. Prefer the longest sequence when records have equal
#    completeness
# 7. Export retained and removed records with summaries
# ==========================================================


# ----------------------------------------------------------
# 0. Setup
# ----------------------------------------------------------

required_packages <- c(
  "dplyr",
  "readr",
  "stringr",
  "tidyr"
)

packages_to_install <- required_packages[
  !sapply(required_packages, requireNamespace, quietly = TRUE)
]

if (length(packages_to_install) > 0) {
  install.packages(packages_to_install)
}

library(dplyr)
library(readr)
library(stringr)
library(tidyr)


# ----------------------------------------------------------
# 1. Set working directory
# ----------------------------------------------------------

# Change this path to the folder containing sequences(3).csv

setwd(
  "/Users/IsaacOmara/Desktop/Manuscripts/Review_Article/Africa_Map/"
)


# ----------------------------------------------------------
# 2. Read the dataset
# ----------------------------------------------------------

rvfv_raw <- read_csv(
  "sequences.csv",
  show_col_types = FALSE
)

cat("\nColumns found in the input file:\n")
print(names(rvfv_raw))

cat("\nInitial number of records:", nrow(rvfv_raw), "\n")


# ----------------------------------------------------------
# 3. Confirm that essential columns are present
# ----------------------------------------------------------

required_columns <- c(
  "Accession",
  "Isolate",
  "Segment",
  "Length",
  "Nuc_Completeness",
  "Country"
)

missing_columns <- setdiff(
  required_columns,
  names(rvfv_raw)
)

if (length(missing_columns) > 0) {
  stop(
    paste(
      "The following required columns are missing:",
      paste(missing_columns, collapse = ", ")
    )
  )
}


# ----------------------------------------------------------
# 4. Standardise metadata fields
# ----------------------------------------------------------

rvfv_standardised <- rvfv_raw %>%
  mutate(
    across(
      where(is.character),
      ~ str_squish(.x)
    ),
    
    Accession = na_if(Accession, ""),
    Isolate = na_if(Isolate, ""),
    Segment = str_to_upper(Segment),
    Country = na_if(Country, ""),
    Nuc_Completeness = str_to_lower(Nuc_Completeness),
    
    Length = suppressWarnings(
      as.numeric(Length)
    )
  )


# ----------------------------------------------------------
# 5. Restrict the analysis to RVFV L, M and S segments
# ----------------------------------------------------------

rvfv_segments <- rvfv_standardised %>%
  filter(
    Segment %in% c("L", "M", "S")
  )


# ----------------------------------------------------------
# 6. Remove records without accession numbers
# ----------------------------------------------------------

records_missing_accession <- rvfv_segments %>%
  filter(
    is.na(Accession)
  )

rvfv_with_accession <- rvfv_segments %>%
  filter(
    !is.na(Accession)
  )


# ----------------------------------------------------------
# 7. Remove repeated accession numbers
#
# If the same accession occurs more than once, retain the
# record with the greatest amount of usable metadata.
# ----------------------------------------------------------

rvfv_accession_ranked <- rvfv_with_accession %>%
  mutate(
    metadata_score =
      if_else(!is.na(Isolate), 1L, 0L) +
      if_else(!is.na(Country), 1L, 0L) +
      if_else(!is.na(Segment), 1L, 0L) +
      if_else(!is.na(Nuc_Completeness), 1L, 0L) +
      if_else(!is.na(Length), 1L, 0L)
  ) %>%
  arrange(
    Accession,
    desc(metadata_score),
    desc(Length)
  )

rvfv_unique_accessions <- rvfv_accession_ranked %>%
  group_by(Accession) %>%
  slice(1) %>%
  ungroup()

duplicate_accessions_removed <- rvfv_accession_ranked %>%
  group_by(Accession) %>%
  filter(row_number() > 1) %>%
  ungroup()


# ----------------------------------------------------------
# 8. Create one searchable metadata string per row
#
# This permits vaccine and laboratory terms to be searched
# across all available character metadata fields.
# ----------------------------------------------------------

character_columns <- names(
  rvfv_unique_accessions
)[
  vapply(
    rvfv_unique_accessions,
    is.character,
    logical(1)
  )
]

rvfv_searchable <- rvfv_unique_accessions %>%
  unite(
    col = "search_text",
    all_of(character_columns),
    sep = " | ",
    remove = FALSE,
    na.rm = TRUE
  )


# ----------------------------------------------------------
# 9. Define exclusion terms
# ----------------------------------------------------------

exclusion_pattern <- paste(
  c(
    "\\bsmithburn\\b",
    "\\bclone[ -]?13\\b",
    "\\bmp[ -]?12\\b",
    "\\brvfv[ -]?4s\\b",
    "\\bddvax\\b",
    "\\bvaccine\\b",
    "\\bvaccinal\\b",
    "\\battenuated\\b",
    "\\blaboratory[- ]adapted\\b",
    "\\blab[- ]adapted\\b",
    "\\blaboratory[- ]passaged\\b",
    "\\blab[- ]passaged\\b",
    "\\bpassaged\\b",
    "\\bpassage\\b",
    "\\bcell[- ]culture\\b",
    "\\btissue[- ]culture\\b",
    "\\brecombinant\\b",
    "\\bsynthetic\\b",
    "\\bengineered\\b",
    "\\breverse[- ]genetics\\b",
    "\\bplasmid\\b",
    "\\bconstruct\\b"
  ),
  collapse = "|"
)


# ----------------------------------------------------------
# 10. Flag and remove excluded records
# ----------------------------------------------------------

rvfv_flagged <- rvfv_searchable %>%
  mutate(
    exclusion_flag = str_detect(
      search_text,
      regex(
        exclusion_pattern,
        ignore_case = TRUE
      )
    )
  )

excluded_vaccine_lab_records <- rvfv_flagged %>%
  filter(exclusion_flag) %>%
  select(-search_text)

rvfv_field_sequences <- rvfv_flagged %>%
  filter(!exclusion_flag) %>%
  select(
    -search_text,
    -exclusion_flag
  )


# ----------------------------------------------------------
# 11. Create a robust isolate identifier
#
# If an isolate ID is present, use it.
# If it is missing, use the accession number so that unrelated
# unnamed sequences are not grouped together.
#
# Country is included to avoid merging isolates with the same
# name reported from different countries.
# ----------------------------------------------------------

rvfv_field_sequences <- rvfv_field_sequences %>%
  mutate(
    isolate_clean = case_when(
      !is.na(Isolate) &
        !str_to_lower(Isolate) %in% c(
          "unknown",
          "not available",
          "not applicable",
          "na",
          "n/a",
          "none"
        ) ~ Isolate,
      
      TRUE ~ paste0(
        "ACCESSION_",
        Accession
      )
    ),
    
    country_clean = if_else(
      is.na(Country),
      "Unknown country",
      Country
    ),
    
    isolate_key = paste(
      country_clean,
      isolate_clean,
      sep = " | "
    )
  )


# ----------------------------------------------------------
# 12. Rank sequence completeness
#
# Exact matching is used so that "partial" is never mistaken
# for "complete".
# ----------------------------------------------------------

rvfv_ranked <- rvfv_field_sequences %>%
  mutate(
    completeness_rank = case_when(
      str_detect(
        Nuc_Completeness,
        regex(
          "^complete$",
          ignore_case = TRUE
        )
      ) ~ 3L,
      
      str_detect(
        Nuc_Completeness,
        regex(
          "complete",
          ignore_case = TRUE
        )
      ) &
        !str_detect(
          Nuc_Completeness,
          regex(
            "partial|incomplete",
            ignore_case = TRUE
          )
        ) ~ 3L,
      
      str_detect(
        Nuc_Completeness,
        regex(
          "partial|incomplete",
          ignore_case = TRUE
        )
      ) ~ 2L,
      
      TRUE ~ 1L
    ),
    
    length_rank = if_else(
      is.na(Length),
      -1,
      Length
    )
  )


# ----------------------------------------------------------
# 13. Identify repeated isolate-segment combinations
# ----------------------------------------------------------

isolate_segment_duplicates <- rvfv_ranked %>%
  add_count(
    isolate_key,
    Segment,
    name = "records_per_isolate_segment"
  ) %>%
  filter(
    records_per_isolate_segment > 1
  ) %>%
  arrange(
    isolate_key,
    Segment,
    desc(completeness_rank),
    desc(length_rank)
  )


# ----------------------------------------------------------
# 14. Retain the best sequence per isolate × segment
#
# Priority:
# 1. Complete sequence
# 2. Longest sequence
# 3. Accession number as a reproducible tie-breaker
# ----------------------------------------------------------

rvfv_clean <- rvfv_ranked %>%
  arrange(
    isolate_key,
    Segment,
    desc(completeness_rank),
    desc(length_rank),
    Accession
  ) %>%
  group_by(
    isolate_key,
    Segment
  ) %>%
  slice(1) %>%
  ungroup()


# ----------------------------------------------------------
# 15. Identify records removed during isolate-segment
# deduplication
# ----------------------------------------------------------

duplicate_isolate_segment_records_removed <- rvfv_ranked %>%
  anti_join(
    rvfv_clean %>%
      select(Accession),
    by = "Accession"
  ) %>%
  mutate(
    removal_reason =
      "Additional record for the same isolate and segment"
  )


# ----------------------------------------------------------
# 16. Create isolate-level completeness table
# ----------------------------------------------------------

isolate_summary <- rvfv_clean %>%
  group_by(
    isolate_key,
    country_clean,
    isolate_clean
  ) %>%
  summarise(
    n_segment_records = n(),
    segments_present = paste(
      sort(unique(Segment)),
      collapse = "+"
    ),
    
    has_L = any(Segment == "L"),
    has_M = any(Segment == "M"),
    has_S = any(Segment == "S"),
    
    complete_L = any(
      Segment == "L" &
        completeness_rank == 3
    ),
    
    complete_M = any(
      Segment == "M" &
        completeness_rank == 3
    ),
    
    complete_S = any(
      Segment == "S" &
        completeness_rank == 3
    ),
    
    .groups = "drop"
  ) %>%
  mutate(
    has_all_three_segments =
      has_L & has_M & has_S,
    
    has_complete_LMS =
      complete_L &
      complete_M &
      complete_S,
    
    segment_pattern = case_when(
      has_L & has_M & has_S ~ "L+M+S",
      has_L & has_M ~ "L+M",
      has_L & has_S ~ "L+S",
      has_M & has_S ~ "M+S",
      has_L ~ "L only",
      has_M ~ "M only",
      has_S ~ "S only",
      TRUE ~ "No recognised segment"
    )
  )


# ----------------------------------------------------------
# 17. Country-level sequence summary
# ----------------------------------------------------------

country_sequence_summary <- rvfv_clean %>%
  count(
    country_clean,
    Segment,
    name = "n_sequences"
  ) %>%
  pivot_wider(
    names_from = Segment,
    values_from = n_sequences,
    values_fill = 0
  ) %>%
  mutate(
    total_segment_sequences =
      coalesce(L, 0L) +
      coalesce(M, 0L) +
      coalesce(S, 0L)
  ) %>%
  arrange(
    desc(total_segment_sequences)
  )


# ----------------------------------------------------------
# 18. Country-level isolate and completeness summary
# ----------------------------------------------------------

country_isolate_summary <- isolate_summary %>%
  group_by(country_clean) %>%
  summarise(
    n_unique_isolates = n(),
    n_isolates_with_LMS = sum(
      has_all_three_segments,
      na.rm = TRUE
    ),
    n_isolates_with_complete_LMS = sum(
      has_complete_LMS,
      na.rm = TRUE
    ),
    
    pct_isolates_with_LMS = round(
      100 *
        n_isolates_with_LMS /
        n_unique_isolates,
      1
    ),
    
    pct_isolates_with_complete_LMS = round(
      100 *
        n_isolates_with_complete_LMS /
        n_unique_isolates,
      1
    ),
    
    .groups = "drop"
  ) %>%
  arrange(
    desc(n_unique_isolates)
  )


# ----------------------------------------------------------
# 19. Segment summary
# ----------------------------------------------------------

segment_summary <- rvfv_clean %>%
  count(
    Segment,
    Nuc_Completeness,
    name = "n_sequences"
  ) %>%
  arrange(
    Segment,
    desc(n_sequences)
  )


# ----------------------------------------------------------
# 20. Create filtering summary
# ----------------------------------------------------------

filtering_summary <- tibble(
  filtering_stage = c(
    "Input records",
    "L, M and S segment records",
    "Records without accession numbers",
    "Repeated accession records removed",
    "Vaccine/laboratory/recombinant records removed",
    "Additional isolate-segment records removed",
    "Final retained segment sequences",
    "Final unique isolates",
    "Countries represented"
  ),
  
  n_records = c(
    nrow(rvfv_raw),
    nrow(rvfv_segments),
    nrow(records_missing_accession),
    nrow(duplicate_accessions_removed),
    nrow(excluded_vaccine_lab_records),
    nrow(duplicate_isolate_segment_records_removed),
    nrow(rvfv_clean),
    n_distinct(rvfv_clean$isolate_key),
    n_distinct(rvfv_clean$country_clean)
  )
)


# ----------------------------------------------------------
# 21. Remove temporary ranking variables from final dataset
# ----------------------------------------------------------

rvfv_clean_export <- rvfv_clean %>%
  select(
    -metadata_score,
    -completeness_rank,
    -length_rank
  )


# ----------------------------------------------------------
# 22. Export output files
# ----------------------------------------------------------

write_csv(
  rvfv_clean_export,
  "RVFV_clean_unique_isolate_segment_dataset.csv"
)

write_csv(
  duplicate_accessions_removed,
  "RVFV_removed_duplicate_accessions.csv"
)

write_csv(
  duplicate_isolate_segment_records_removed,
  "RVFV_removed_duplicate_isolate_segment_records.csv"
)

write_csv(
  excluded_vaccine_lab_records,
  "RVFV_excluded_vaccine_lab_recombinant_records.csv"
)

write_csv(
  records_missing_accession,
  "RVFV_records_missing_accession.csv"
)

write_csv(
  isolate_segment_duplicates,
  "RVFV_all_repeated_isolate_segment_records.csv"
)

write_csv(
  isolate_summary,
  "RVFV_unique_isolate_summary.csv"
)

write_csv(
  country_sequence_summary,
  "RVFV_country_sequence_summary.csv"
)

write_csv(
  country_isolate_summary,
  "RVFV_country_isolate_completeness_summary.csv"
)

write_csv(
  segment_summary,
  "RVFV_segment_completeness_summary.csv"
)

write_csv(
  filtering_summary,
  "RVFV_filtering_summary.csv"
)


# ----------------------------------------------------------
# 23. Print final summary
# ----------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("RVFV DATASET CLEANING SUMMARY\n")
cat("====================================================\n")

print(filtering_summary)

cat("\nSequences retained by segment:\n")

print(
  rvfv_clean %>%
    count(
      Segment,
      name = "n_sequences"
    )
)

cat("\nTop countries by unique isolates:\n")

print(
  country_isolate_summary %>%
    slice_head(n = 30)
)

cat("\nOutput files have been written to:\n")
cat(getwd(), "\n")
