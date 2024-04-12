###############################################################################
# Script: matching_table.py
# Author: Kenneth Valerio Aguilar
# Date: 12/20/2023
# Version: 1.0
# Purpose: Create table with matching information between the genotyping and pedigree
# Input Requirements: serotype_table.tsv and summary_table.csv
# Usage: run the script like: python3 serotype_table.py
# Output: summary matching_table.csv 
################################################################################

#Import required packages
import pandas as pd
import os
import re
import numpy as np

#read in data
serotypes = pd.read_csv('/molbio/projects/hla_genotyping/hlahd/CD_analysis/serotype_table.tsv', sep='\t')
data = pd.read_csv('/molbio/projects/hla_genotyping/hlahd/CD_analysis/summary_table.csv',sep=',')


# Replace '-' for actual values
data.loc[data['DQA1_Gtype2'] == '-', 'DQA1_Gtype2'] = data['DQA1_Gtype1']
data.loc[data['DQB1_Gtype2'] == '-', 'DQB1_Gtype2'] = data['DQB1_Gtype1']

#create empty lists to fill with for loop
serotypes_dict = {}

for index1, row1 in data.iterrows():
    for index2, row2 in serotypes.iterrows():
        # Check if DQA1_Gtype1 equals DQA1_Gtype2 and DQB1_Gtype1 equals DQB1_Gtype2
        if row1['DQA1_Gtype1'] == row1['DQA1_Gtype2'] and row1['DQB1_Gtype1'] == row1['DQB1_Gtype2']:
            if (row2['DQA1'] in row1['DQA1_Gtype1'][9:14]) and (row2['DQB1'] in row1['DQB1_Gtype1'][9:14]):
                if row1['Basename'] not in serotypes_dict.keys():
                    serotypes_dict[row1['Basename']] = [row2['Subtype'], row2['Subtype']]
                else:
                    serotypes_dict[row1['Basename']].extend([row2['Subtype'], row2['Subtype']])
        else:
            # Check for DQA1_Gtype1 and DQB1_Gtype1
            if (row2['DQA1'] in row1['DQA1_Gtype1'][9:14]) and (row2['DQB1'] in row1['DQB1_Gtype1'][9:14]):
                if row1['Basename'] not in serotypes_dict.keys():
                    serotypes_dict[row1['Basename']] = [row2['Subtype']]
                else:
                    serotypes_dict[row1['Basename']].append(row2['Subtype'])
            # Check for DQA1_Gtype1 and DQB1_Gtype2
            elif (row2['DQA1'] in row1['DQA1_Gtype1'][9:14]) and (row2['DQB1'] in row1['DQB1_Gtype2'][9:14]):
                if row1['Basename'] not in serotypes_dict.keys():
                    serotypes_dict[row1['Basename']] = [row2['Subtype']]
                else:
                    serotypes_dict[row1['Basename']].append(row2['Subtype'])
            # Check for DQA1_Gtype2 and DQB1_Gtype1
            elif (row2['DQA1'] in row1['DQA1_Gtype2'][9:14]) and (row2['DQB1'] in row1['DQB1_Gtype1'][9:14]):
                if row1['Basename'] not in serotypes_dict.keys():
                    serotypes_dict[row1['Basename']] = [row2['Subtype']]
                else:
                    serotypes_dict[row1['Basename']].append(row2['Subtype'])
            # Check for DQA1_Gtype2 and DQB1_Gtype2
            elif (row2['DQA1'] in row1['DQA1_Gtype2'][9:14]) and (row2['DQB1'] in row1['DQB1_Gtype2'][9:14]):
                if row1['Basename'] not in serotypes_dict.keys():
                    serotypes_dict[row1['Basename']] = [row2['Subtype']]
                else:
                    serotypes_dict[row1['Basename']].append(row2['Subtype'])
 
for basename, subtypes in serotypes_dict.items():
    while len(subtypes) < 4:
        subtypes.append("NA")



result_df = pd.DataFrame.from_dict(serotypes_dict, orient='index')
# Rename columns
result_df.columns = ['Serotype_1', 'Serotype_2', 'Serotype_3', 'Serotype_4']

# Reset index to move 'Basename' from index to a regular column
result_df.reset_index(inplace=True)
result_df.rename(columns={'index': 'Basename'}, inplace=True)

#lets add column with Korponay coding
Program_obtained_coding = [] 
Korponay_coding = {'1':2.5,'2':2.2,'3':7.5,'8':8.1}

#lets add column with Korponay coding
Korponay_coding = {'1':2.5,'2':2.2,'3':7.5,'8':8.1}
#create function to classify
def map_serotype_to_coding(row):
    program_coding = ''
    for serotype in ['Serotype_1', 'Serotype_2', 'Serotype_3', 'Serotype_4']:
        if row[serotype] != 'NA':
            if row[serotype] in Korponay_coding.values():
                for key, value in Korponay_coding.items():
                    if row[serotype] == value:
                        program_coding += key
                        break
            else:
                program_coding += 'X'
    
    # Reorder the characters with 'X' values at the right
    program_coding = ''.join(sorted(program_coding, key=lambda x: (x == 'X', x)))

        # Cut out the last character if the length is greater than 2
    if len(program_coding) > 2:
        program_coding = program_coding[:2]
    
    return program_coding
    
 
# Apply the function to create a new column
result_df['program_Coding'] = result_df.apply(map_serotype_to_coding, axis=1)

#add data of serotypes determined in the lab
Lab_serotypes = {'CS01':'3X','CS11':'1X', 'CS12':'3X','CS13':'1X',\
                'CS21':'XX', 'CS22':'1X','CS23':'1X','CS24':'1X',\
                'CS25':'1X','CS32':'1X','CS41':'38','CS51':'12',\
                'CS52':'13','CS53':'12','CS54':'12','CS55':'11',\
                'CS56':'23','CS57':'11','CS61':'3X','CS62':'8X',\
                    'CS63':'8X', 'BS051':'13','BS052':'2X','BS053':'12',\
                        'BS054':'1X','BS055':'23','BS056':'23','BS061':'3X',\
                            'BS062':'1X','BS063':'1X','BS064':'1X','BS065':'1X',\
                                'BS066':'1X','BS067':'1X','BS068':'3X','BS069':'3X',\
                                    'BS0691':'XX','BS070':'12','BS071':'1X','BS073':'13',\
                                        'BS074':'13','BS075':'13','BS081':'13','BS091':'12',\
                                            'BS092':'18','BS093':'12','BS094':'11','BS095':'11',\
                                                'BS096':'12','BS097':'11','BS098':'28','BS099':'12',\
                                                    'BS100':'28','BS101':'12','BS102':'12','BS103':'11',\
                                                        'BS104':'28','BS105':'18'}


#create new column to compare the Lab_serotypes with the program_serotypes
result_df['Lab_serotypes'] = result_df['Basename'].map(Lab_serotypes)

#create comparison column
result_df['New_Column'] = np.where(result_df['program_Coding'] == result_df['Lab_serotypes'], 'Yes', 'No')

#combine the table with the 
result_df = pd.merge(data, result_df, how = 'left', on= 'Basename')

#create a no only table
no_match = result_df[result_df['New_Column'].str.contains('No')]

# export results
result_df.to_csv('/molbio/projects/hla_genotyping/hlahd/CD_analysis/matching_table.csv', sep='\t', index=False)
no_match.to_csv('/molbio/projects/hla_genotyping/hlahd/CD_analysis/no_matching_table.csv', sep='\t', index=False)