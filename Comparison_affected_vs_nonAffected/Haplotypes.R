###############################################################################
# Script: freq_analysis_haplotypes.R
# Author: Kenneth Valerio Aguilar
# Date: 01/29/2023
# Version: 1.0
# Purpose: Get the frequency of each haplotype inside each group
# Input Requirements: None
# Usage: run the script 
# Output:
################################################################################

#Call the required packages
library(tidyverse)
library(plotrix)
library(readxl)
library(tidyr)

#import the data.
HLADQ<-read.csv('/Documents and Settings/Kenneth/Documents/HLA/HLA/summary_table.csv',header = T)
Pedigree<-read_xlsx('/Documents and Settings/Kenneth/Documents/HLA/HLA/full-phenotype-table.xlsx')

#First transform the '-' values for the exact values
HLADQ$DQA1_Gtype2[HLADQ$DQA1_Gtype2 == "-"] <- HLADQ$DQA1_Gtype1[HLADQ$DQA1_Gtype2 == "-"]
HLADQ$DQB1_Gtype2[HLADQ$DQB1_Gtype2 == "-"] <- HLADQ$DQB1_Gtype1[HLADQ$DQB1_Gtype2 == "-"]

#In order to analyse both tables, first they need to be cleaned and joined together:
#trim the HLADQ columns
extract_after_star <- function(x) {
  sub(".*\\*", "", x)
}

# Use mutate_at to apply the function to specified columns
HLADQ2 <- HLADQ %>%mutate_at(vars(-Basename), funs(extract_after_star))

# for the pedigree, first filter the useful rows
Pedigree2<-Pedigree%>%drop_na(`Sample name`)%>%rename('Basename'='Sample name')

#join both tables by individual
data1<-left_join(HLADQ2,Pedigree2,"Basename")

#Now, in this particular case, we know that the alleles of the HLA-DQ region are expressed as haplotypes, so the analysis must be carried using the haplotypes instead of the single alleles. 
#The next code will combine the alleles and create 4 columns for the 4 possible combination of these alleles:
#create new column with the 4 possible halpotypes

data2<-data1%>%mutate(haplotype_1=(str_c(data1$DQA1_Gtype1,'-',data1$DQB1_Gtype1)))%>%
  mutate(haplotype_2=(str_c(data1$DQA1_Gtype1,'-',data1$DQB1_Gtype2)))%>%
  mutate(haplotype_3=(str_c(data1$DQA1_Gtype2,'-',data1$DQB1_Gtype1)))%>%
  mutate(haplotype_4=(str_c(data1$DQA1_Gtype2,'-',data1$DQB1_Gtype2)))

# Define a function to extract values after the '-' symbol and remove consecutive '-'
extract_after_dash_and_remove_consecutive <- function(x) {
  ifelse(x == '---', '', str_replace_all(x, '--', '-'))}
data3 <- data2 %>%mutate_all(~extract_after_dash_and_remove_consecutive(.))

#Now using the cleaned and transformed data, It is possible to create a new data which contains the count of appearance of each haplotype for Susceptible vs not-suceptible individuals. (The 'DQ' column of the pedigree table is used instead of the 'Affected' one)
#colapse the haplotypes in a single column
data4<-data3%>%select(Basename,Affected,haplotype_1,haplotype_2,haplotype_3,haplotype_4)%>%
  gather(key = 'Haplotype',value = 'values',-Basename,-Affected)%>%count(Affected,values)%>%
  mutate(Affected=as.factor(Affected))%>%mutate(values=as.character(values))%>%group_by(Affected)%>%mutate(Percentage=(n/sum(n))*100)%>%filter(nchar(values)>8)
view(data4)

#Using this table is posible to create the heatmap:
library(ggtext)



p1 <- data4 %>%
  ggplot(aes(x = Affected, y = values, fill = Percentage)) +
  geom_tile(colour = "white", size = 0.25) +
  labs(x = "Affected", y = "HLADQA-HLADQB Combinations") + 
  scale_y_discrete(expand = c(0, 0)) +
  theme_grey(base_size = 8) +
  theme(
    # bold font for legend text
    legend.text = element_text(face = "bold"),
    # set thickness of axis ticks
    axis.ticks = element_line(size = 0.5),
    # remove plot background
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
  
  
  # Add annotation
  annotate(
    geom = "text", 
    x = "Y",  # Adjust the x-value as needed
    y = "05:01:01-02:01:01",    # Adjust the y-value as needed
    label = "40.27%", 
    hjust = "rigth",
    colour = "black",
    size = 3
  )+
  annotate(
    geom = "text", 
    x = "N",  # Adjust the x-value as needed
    y = "05:01:01-02:01:01",    # Adjust the y-value as needed
    label = "18.58%", 
    hjust = "rigth",
    colour = "black",
    size = 3
  )+
  annotate(
    geom = "text", 
    x = "N",  # Adjust the x-value as needed
    y = "05:05:01-03:01:01",    # Adjust the y-value as needed
    label = "11.53%", 
    hjust = "rigth",
    colour = "white",
    size = 3
  )+
  annotate(
    geom = "text", 
    x = "Y",  # Adjust the x-value as needed
    y = "05:01:01-05:01:01",    # Adjust the y-value as needed
    label = "9.72%", 
    hjust = "rigth",
    colour = "white",
    size = 3
  )+
  annotate(
    geom = "text", 
    x = "Y",  # Adjust the x-value as needed
    y = "01:01:01-05:01:01",    # Adjust the y-value as needed
    label = "8.33%", 
    hjust = "rigth",
    colour = "white",
    size = 3
  )+
  annotate(
    geom = "text", 
    x = "N",  # Adjust the x-value as needed
    y = "02:01-02:02:01",    # Adjust the y-value as needed
    label = "7.05%", 
    hjust = "rigth",
    colour = "white",
    size = 3
  )
  

png("/Documents and Settings/Kenneth/Documents/HLA/plot12.png", units="in", width=5, height=5, res=300)
p1
dev.off()