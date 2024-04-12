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

#For each HLA Region In the affected individuals, extract the Alleles which differ more (by its count) from the rest.
#For this first task I First created a count dataframe for each HLA region using the table 'summary_complete_table.csv':

#import the data
whole1<-read.csv('/Documents and Settings/Kenneth/Documents/HLA/HLA/summary_complete_table.csv',sep = '\t',header = T)
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

#apply the transformation function
whole<-transform_table(whole1)

#trim the columns
extract_after_star <- function(x) {
  sub(".*\\*", "", x)
}

# for the pedigree, first filter the useful rows
Pedigree2<-Pedigree%>%drop_na(`Sample name`)%>%rename('Basename'='Sample name')

#clean the data 
whole2<-sapply(whole, extract_after_star)
whole2<-data.frame(whole2)

#combine with the pedigree data
whole3<-left_join(whole2,Pedigree2,"Basename")



#prueba:
whole4<-whole3%>%select(-c('M/F','DQ','DQ alleles','Father','Mother','Family','Sequenced'))
names(whole4)
View()

whole4%>%pivot_longer(cols = -c(Basename, Affected),names_to =c(".value", "Gtype"),names_pattern = "(.*?)_Gtype(\\d+)")

merge_columns <- function(data) {
  # Pivot the columns _G1 and _G2 into long format
  data_long <- pivot_longer(data, 
                            cols = -c(Basename, Affected), 
                            names_to = c(".value", "Gtype"), 
                            names_pattern = "(.*?)_Gtype(\\d+)")
  
  
  return(data_long)
}

# Example usage
# Assuming 'data' is your original data frame
transformed_data <- merge_columns(whole4)
transformed_data<-transformed_data%>%select(-'Gtype')
view(transformed_data)



#Get the region names
Regions<-names(transformed_data)[-1]

#Create a list with all data frames for each region and remove the non genotyped
create_count_dataframes <- function(df, column_names) {
  result_list <- list()
  for (col_name in column_names) {
    # Get unique values and their counts in the column
    count_df <- df %>%
      count(Affected, !!sym(col_name)) %>%
      mutate(Affected = as.factor(Affected), count = n) %>%
      select(Affected, count,!!sym(col_name))%>%filter(!(!!sym(col_name) %in% c('-','Not_typed')))
    
    result_list[[col_name]] <- count_df
  }
  
  return(result_list)
}

#apply
result<-create_count_dataframes(transformed_data,Regions)
head(result)
result
#For this task it is first necessary to calculate the the difference
#in counts inside each group.
#Then in order to compare the counts between group it is
#also necessary to normalize the data by the number of individuals 
#inside each group. After that it is possible to obtain
#the alleles that are expressed with most difference between both Y and N groups:

Affected_vs_not <- function(list) {
  comparison_list <- list()
  
  for (x in names(list)) {
    # Check if 'Affected' column exists with factors 'Y' and 'N'
    if (!('Affected' %in% names(list[[x]])) || 
        !(all(c('Y', 'N') %in% unique(list[[x]]$Affected)))) {
      # Skip current element of the list if the column or factors are missing
      next
    }
    third_col_name <- names(list[[x]])[3]
    cat("Third column name:", third_col_name, "\n")
    
    # Check if third_col_name is not NA
    if (is.na(third_col_name)) {
      cat("Skipping current element as third_col_name is NA.\n")
      next
    }
    # Obtain the number of individuals in each category
    df1 <- list[[x]] %>% count(Affected)
    df2 <- left_join(list[[x]], df1, by = 'Affected')
    # Calculate the most different alleles between affected vs non-affected individuals
    df3 <- df2 %>% group_by(Affected)%>%mutate(Mean = mean(count),
                          Diff_inside_group = ((count - Mean)^2 / n)) %>%
      arrange(third_col_name, desc(Diff_inside_group)) %>%
      filter(Affected %in% c('Y', 'N')) %>%
      pivot_wider(names_from = Affected, values_from = Diff_inside_group) %>%
      select(third_col_name, Y, N) %>%
      group_by(!!sym(third_col_name)) %>%
      mutate_all(~ifelse(is.na(.), "", as.character(.))) %>%
      filter(!(!!sym(third_col_name) %in% c('-', 'Not_typed'))) %>%
      summarise_all(funs(trimws(paste(., collapse = '')))) %>%
      mutate_at(c('Y', 'N'), as.numeric) %>%
      mutate_all(~replace(., is.na(.), 0)) %>%
      mutate(Diff_between_groups = (Y - N)^2) %>%
      arrange(desc(Diff_between_groups))
    
    comparison_list[[x]] <- df3
  }
  return(comparison_list)
}


#apply the function
F1<-Affected_vs_not(result)

#create another function to pivot the count data for later merging
pivot_tables <- function(list_of_dataframes) {
  pivoted_list <- lapply(list_of_dataframes, function(df) {
    pivoted_df <- pivot_wider(df, names_from = "Affected", values_from = "count")
    
    # Replace NA values with 0
    if ("Y" %in% colnames(pivoted_df)) {
      pivoted_df[is.na(pivoted_df$Y), "Y"] <- 0
    }
    if ("N" %in% colnames(pivoted_df)) {
      pivoted_df[is.na(pivoted_df$N), "N"] <- 0
    }
    
    # Calculate percentages and add new columns if 'Y' and 'N' exist
    if ("Y" %in% colnames(pivoted_df) & "N" %in% colnames(pivoted_df)) {
      pivoted_df$percentage_Y <- pivoted_df$Y / colSums(pivoted_df["Y"], na.rm = TRUE)
      pivoted_df$percentage_N <- pivoted_df$N / colSums(pivoted_df["N"], na.rm = TRUE)
    }
    
    # Remove count columns if 'Y' and 'N' exist
    if ("Y" %in% colnames(pivoted_df)) {
      pivoted_df <- pivoted_df[, !(names(pivoted_df) %in% c("Y"))]
    }
    if ("N" %in% colnames(pivoted_df)) {
      pivoted_df <- pivoted_df[, !(names(pivoted_df) %in% c("N"))]
    }
    
    return(pivoted_df)
  })
  return(pivoted_list)
}


#apply the funcion
pivotedtable<-pivot_tables(result)

#create another funciton to join both tables and arrange the final table based on the Difference magnitude column

left_join_by_column <- function(list1, list2) {
  merged_list <- list()
  
  for (i in seq_along(list1)) {
    name <- names(list1)[i]
    if (name %in% names(list2)) {
      merged_df <- merge(list1[[i]], list2[[name]], by.x = name, by.y = name, all.x = TRUE)
      merged_list[[name]] <- merged_df
    } else {
      merged_list[[name]] <- list1[[i]]
    }
  }
  
  # Add data frames from list2 that are not present in list1
  for (j in seq_along(list2)) {
    name <- names(list2)[j]
    if (!(name %in% names(merged_list))) {
      merged_list[[name]] <- list2[[name]]
    }
  }
  #arrange
  merged_list <- lapply(merged_list, function(df) {
    if ("Diff_between_groups" %in% colnames(df)) {
      df <- df[order(-df$Diff_between_groups), ]
    }
    return(df)
  })
     
    
  return(merged_list)
}

#Create the list of final tables
Final_table_list<-left_join_by_column(F1,pivotedtable)

#from this final tbale list i will take into account the most relevant results
Final_table_list