#!/usr/bin/env Rscript

# Load necessary libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

# Function to calculate averages with specific logic
calculate_avg <- function(df, columns) {
  df %>%
    summarise(across(all_of(columns), ~ {
      mean_val <- mean(.x, na.rm = TRUE)
      total_sum <- sum(.x, na.rm = TRUE)
      if (total_sum <= 0) {
        return(total_sum)
      } else {
        return(round(mean_val, 3))
      }
    }))
}

# Function to extract PatID and Type from SampleID
extract_patid_type <- function(sample_id) {
  # Split Sample ID by '~'
  parts <- strsplit(sample_id, "~")[[1]]
  
  # Extract PatID (first part)
  patid <- parts[1]
  
  # Determine Type based on the third part
  # Mapping:
  # - Contains "germline" -> "germline"
  # - Contains "originator" -> "originator"
  # - Contains "pdorg" -> "pdorg"
  # - Contains "pdc" -> "pdc"
  # - Contains "normal" -> "pdx"
  # - Else -> "pdx"
  
  if (length(parts) >= 3) {
    type_part <- tolower(parts[3])
    if (grepl("germline", type_part)) {
      type <- "germline"
    } else if (grepl("ORIGINATOR", type_part)) {
      type <- "originator"
    } else if (grepl("organoid", type_part)) {
      type <- "pdorg"
    } else if (grepl("PDC", type_part)) {
      type <- "pdc"
    } else if (grepl("Normal", type_part)) {
      type <- "germline"
    } else {
      type <- "pdx" # Default to 'pdx' if no specific type found
    }
  } else {
    type <- "pdx" # Default to 'pdx' if insufficient parts
  }
  
  return(list(patid = patid, type = type))
}

# Main processing function
process_ancestry <- function(input_files, output_file) {
  
  # Read all input files, skipping lines starting with '#'
  data_list <- lapply(input_files, function(f) {
    read_tsv(f, comment = "#", col_names = FALSE, col_types = cols())
  })
  
  # Combine all data into one data frame
  data <- bind_rows(data_list)
  
  # Assign column names based on the expected structure
  colnames(data) <- c("SampleID", "Populationlabel", "NumberofSNPs", 
                      "PredictedPC1", "PredictedPC2", "PredictedPC3", "PredictedPC4",
                      "EURancestry", "AMRancestry", "AFRancestry", "SASAncestry",
		      "EASAncestry")
  
  # Extract PatID and Type from SampleID
  extracted <- lapply(data$SampleID, extract_patid_type)
  patid <- sapply(extracted, function(x) x$patid)
  type <- sapply(extracted, function(x) x$type)
  
  # Add PatID and Type columns to the data frame
  data <- data %>%
    mutate(PatID = patid,
           Type = type)
  
  # Group data by PatID and PopulationLabel
  processed_data <- data %>%
    group_by(PatID, Populationlabel) %>%
    group_modify(~ {
      patient_data <- .x
      
      if ("germline" %in% patient_data$Type) {
        # Process as 'germline'
        germline_data <- patient_data %>% filter(Type == "germline")
        avg_values <- calculate_avg(germline_data, c("NumberofSNPs",
                      "PredictedPC1", "PredictedPC2", "PredictedPC3", "PredictedPC4",
                      "EURancestry", "AMRancestry", "AFRancestry", "SASAncestry",
                      "EASAncestry"))
        # Create output row
        output_row <- avg_values %>%
          mutate(PatID = unique(germline_data$PatID),
                 PopulationLabel = unique(germline_data$Populationlabel),
                 Source = "germline")
        return(output_row)
        
      } else if ("originator" %in% patient_data$Type) {
        # Process as 'originator'
        originator_data <- patient_data %>% filter(Type == "originator")
        avg_values <- calculate_avg(originator_data, c("NumberofSNPs",
                      "PredictedPC1", "PredictedPC2", "PredictedPC3", "PredictedPC4",
                      "EURancestry", "AMRancestry", "AFRancestry", "SASAncestry",
                      "EASAncestry"))
        # Create output row
        output_row <- avg_values %>%
          mutate(PatID = unique(originator_data$PatID),
                 PopulationLabel = unique(originator_data$Populationlabel),
                 Source = "originator")
        return(output_row)
        
      } else {
        # Process as 'PDX'
        pdx_data <- patient_data %>% filter(Type %in% c("pdx", "pdorg", "pdc"))
        avg_values <- calculate_avg(pdx_data, c("NumberofSNPs",
                      "PredictedPC1", "PredictedPC2", "PredictedPC3", "PredictedPC4",
                      "EURancestry", "AMRancestry", "AFRancestry", "SASAncestry",
                      "EASAncestry"))

        # Create output row
        output_row <- avg_values %>%
          mutate(PatID = unique(pdx_data$PatID),
                 PopulationLabel = unique(pdx_data$Populationlabel),
                 Source = "PDX")
        return(output_row)
      }
    }) %>%
    ungroup() %>%
    select(PatID, Populationlabel, NumberofSNPs, PredictedPC1, 
           PredictedPC2, PredictedPC3, PredictedPC4, EURancestry, 
	   AMRancestry, AFRancestry, SASAncestry, EASAncestry, Source)
  
  # Prepare header
  header <- "#Sample ID\tPopulation label\tNumber of SNPs\tPredicted PC1\tPredicted PC2\tPredicted PC3\tPredicted PC4\t% EUR ancestry\t% AMR ancestry\t% AFR ancestry\t% SAS Ancestry\t% EAS Ancestry\tSource"
  
  # Write header to the output file
  write_lines(header, output_file)
  
  # Write processed data to the output file, appending after the header
  write_tsv(processed_data, output_file, append = TRUE, col_names = FALSE)
}

# Entry point of the script
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: process_patient_ancestry.R <output_file> <input_file1> <input_file2> ...")
}

output_file <- args[1]
input_files <- args[-1]

# Execute the processing function
process_ancestry(input_files, output_file)

