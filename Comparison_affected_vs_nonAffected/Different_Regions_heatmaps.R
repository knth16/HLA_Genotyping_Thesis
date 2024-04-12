###############################################################################
# Script: Difference_magnitude.R
# Author: Kenneth Valerio Aguilar
# Date: 01/29/2023
# Version: 1.0
# Purpose: Get the Difference magnitude of alleles
# Input Requirements: None
# Usage: run the script 
# Output:
################################################################################

#Call the required packages
library(tidyverse)
library(plotrix)
library(readxl)
library(tidyr)

#import the data
whole<-read.csv('/Documents and Settings/Kenneth/Documents/HLA/HLA/summary_complete_table.csv',sep = '\t',header = T)
Pedigree<-read_xlsx('/Documents and Settings/Kenneth/Documents/HLA/HLA/full-phenotype-table.xlsx')

#transform - for actual values
transform_table <- function(data) {
  for (col_name in names(data)) {
    if (endsWith(col_name, "_Gtype2")) {
      # Extract the corresponding _Gtype1 column name
      col_name_gtype1 <- sub("_Gtype2$", "_Gtype1", col_name)
      
      # Replace '-' values in _Gtype2 column with corresponding values from _Gtype1 column
      data[[col_name]][data[[col_name]] == "-"] <- data[[col_name_gtype1]][data[[col_name]] == "-"]
    }
  }
  return(data)
}

#apply the function
data<-transform_table(whole)

#trim the columns
extract_after_star <- function(x) {
  sub(".*\\*", "", x)
}

# for the pedigree, first filter the useful rows
Pedigree2<-Pedigree%>%drop_na(`Sample name`)%>%rename('Basename'='Sample name')

#clean the data the same way as before
whole2<-sapply(data, extract_after_star)
whole2<-data.frame(whole2)
#combine with the pedigree data
whole3<-left_join(whole2,Pedigree2,"Basename")

#create a funciton for heatmap creation
cross_regions_heatmap <- function(DF, R1_GT1, R1_GT2, R2_GT1, R2_GT2) {
  # Select relevant columns and gather data
  firsth <- DF %>%
    select(R1_GT1, R1_GT2, R2_GT1, R2_GT2, Affected) %>%
    gather(key = GTR2s, value = GTR2allele, -R1_GT1, -R1_GT2, -Affected) %>%
    gather(key = DR1s, value = GTR1allele, -Affected, -GTR2s, -GTR2allele) %>%
    count(Affected, GTR2allele, GTR1allele) %>%
    filter(!GTR2allele == '-') %>%
    filter(!GTR1allele == '-')
  
  # Convert to tibble
  firsth <- as_tibble(firsth)
  
  # Get unique values and generate combinations
  unique_vals_GTR2 <- sort(unique(firsth$GTR2allele))
  unique_vals_GTR1 <- sort(unique(firsth$GTR1allele))
  combinations <- expand.grid(GTR2allele = unique_vals_GTR2, GTR1allele = unique_vals_GTR1)
  
  # Merge with original data to identify missing combinations
  df_merged <- merge(firsth, combinations, by = c("GTR2allele", "GTR1allele"), all = TRUE)
  
  # Replace NA values in "Affected" and "n" columns with "N" and 0 respectively
  df_merged$Affected[is.na(df_merged$Affected)] <- "N"
  df_merged$n[is.na(df_merged$n)] <- 0
  
  # Add missing combinations with 0 count
  df_merged <- df_merged[order(df_merged$Affected, df_merged$GTR2allele, df_merged$GTR1allele), ]
  missing_combinations <- df_merged[df_merged$Affected == "N" & df_merged$n == 0, ]
  
  # Set the "Affected" column to "Y" for missing combinations
  missing_combinations$Affected <- "Y"
  # Set count values to 0 for missing combinations
  missing_combinations$n <- 0
  
  # Add missing combinations to the original data
  df_complete <- rbind(df_merged, missing_combinations)
  
  # Check each combination and add missing rows with "Affected" values of 'Y' and 'N'
  for (i in 1:nrow(combinations)) {
    row <- combinations[i, ]
    matched_rows <- df_complete[df_complete$GTR2allele == row$GTR2allele & df_complete$GTR1allele == row$GTR1allele, ]
    
    # If no matches, continue
    if (nrow(matched_rows) == 0) {
      next
    }
    
    # If there is a match, check the value of the "Affected" column
    unique_affected <- unique(matched_rows$Affected)
    
    # If more than one match and values of the "Affected" column are different, continue
    if (length(unique_affected) > 1) {
      next
    }
    
    # If only one match and the value of the "Affected" column is the same as 'N', add a row with 'Y'
    if (unique_affected == 'N') {
      df_complete <- rbind(df_complete, data.frame(Affected = 'Y', GTR2allele = row$GTR2allele, GTR1allele = row$GTR1allele, n = 0))
    }
    
    # If only one match and the value of the "Affected" column is the same as 'Y', add a row with 'N'
    if (unique_affected == 'Y') {
      df_complete <- rbind(df_complete, data.frame(Affected = 'N', GTR2allele = row$GTR2allele, GTR1allele = row$GTR1allele, n = 0))
      
      
      
    }
  }
  df_complete<-df_complete%>%group_by(Affected)%>%mutate(Percentage = (n/sum(n))*100)
  view(df_complete)
  # Plot heatmap
  plot <- df_complete %>% 
    ggplot(aes(GTR2allele, GTR1allele, fill = Percentage)) +
    geom_tile() +
    facet_grid(~Affected) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      legend.text = element_text(face = "bold"),
      axis.ticks = element_line(size = 0.5),
      plot.background = element_blank(),
      # remove plot border
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      panel.grid.major = element_line(colour = "grey70", linewidth = 0),
      axis.title.y = element_text(size = 9),
      axis.title.x = element_text(size = 9)
      ) +
    scale_fill_viridis_c()+
    labs(x = substr(R2_GT1, 1, 4), y = substr(R1_GT1, 1, 4))
  
  return(plot)
}

#test different results
result1 <- cross_regions_heatmap(whole3, "DRB1_Gtype1", "DRB1_Gtype2", "DQA1_Gtype1", "DQA1_Gtype2")
result2 <- cross_regions_heatmap(whole3, "DQB1_Gtype1", "DQB1_Gtype2", "DQA1_Gtype1", "DQA1_Gtype2")
result3 <- cross_regions_heatmap(whole3, "DPB1_Gtype1", "DPB1_Gtype2", "DQA1_Gtype1", "DQA1_Gtype2")
result4 <- cross_regions_heatmap(whole3, "A_Gtype1", "A_Gtype2", "DQA1_Gtype1", "DQA1_Gtype2")
result5 <- cross_regions_heatmap(whole3, "B_Gtype1", "B_Gtype2", "DQA1_Gtype1", "DQA1_Gtype2")
result6 <- cross_regions_heatmap(whole3, "C_Gtype1", "C_Gtype2", "DQA1_Gtype1", "DQA1_Gtype2")
result7 <- cross_regions_heatmap(whole3, "DPA1_Gtype1", "DPA1_Gtype2", "DQA1_Gtype1", "DQA1_Gtype2")
result8 <- cross_regions_heatmap(whole3, "DMA_Gtype1", "DMA_Gtype2", "DQA1_Gtype1", "DQA1_Gtype2")

