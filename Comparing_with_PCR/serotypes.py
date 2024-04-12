###############################################################################
# Script: serotype_table.py
# Author: Kenneth Valerio Aguilar
# Date: 12/17/2023
# Version: 1.0
# Purpose: Create table with general HLA_DQ serotype variants information
# Input Requirements: None
# Usage: run the script like: python3 serotype_table.py
# Output: summary serotype_table.tsv 
################################################################################

#import required packages
import pandas as pd
import os

#create table with the isomform information from multiple sources
serotypes = pd.DataFrame(
    {"Serotype": ['DQ2','DQ2','DQ2','DQ4','DQ4',"DQ4","DQ5","DQ5","DQ5","DQ5","DQ5","DQ5","DQ5","DQ5","DQ6","DQ6","DQ6","DQ6","DQ6","DQ6","DQ6","DQ6","DQ7","DQ7","DQ7","DQ7","DQ7","DQ7","DQ7","DQ7","DQ8","DQ8","DQ9","DQ9"],
     "Subtype": [2.5,2.2,2.3,4.3,4.3,4.4,5.1,5.1,5.1,5.1,5.2,5.2,5.3,5.4,6.1,6.2,6.2,6.2,6.3,6.3,6.4,6.9,7.2,7.3,7.3,7.3,7.3,7.4,7.5,7.6,8.1,8.1,9.2,9.3],
     "DQA1": ["05:01","02:01","03:02","03:01","03:02","04:01","01:01","01:02","01:03","01:04","01:02","01:03","01:04","01:02","01:03","01:02","01:03","01:04","01:02","01:03","01:02","01:02","02:01","03:01","03:03","03:01","03:02","04:01","05:05","06:01","03:01","03:02","02:01","03:02"],
     "DQB1": ["02:01","02:02","02:02","04:02","04:02","04:02","05:01","05:01","05:01","05:01","05:02","05:02","05:03","05:04","06:01","06:02","06:02","06:02","06:03","06:03","06:04","06:09","03:01","03:01","03:01","03:04","03:04","03:01","03:01","03:01","03:02","03:02","03:03","03:03"]}
)

serotypes.to_csv('/molbio/projects/hla_genotyping/hlahd/CD_analysis/serotype_table.tsv',sep='\t',index=False)