#!/bin/bash

################################################################################
# Script: hlahd_analysis_script.sh
# Author: Kenneth Valerio Aguilar
# Date: 11/18/2023
# Version: 1.0
# Purpose: Perform HLA genotyping analysis using hlahd.sh on multiple individuals
# Input Requirements: Requires fastq files with specific naming convention
# Usage: run the script in teh terminal with any of two arguments (BS or CS)
# Output: Log files, result files, and summary table in the results directory
################################################################################

# Set the project root directory and other variables
root="/molbio/projects/hla_genotyping/hlahd"
freq_data="/usr/local/molbio/hlahd/freq_data"
HLA_genesplit="/usr/local/molbio/hlahd/HLA_gene.split.txt"
Dictionary="/usr/local/molbio/hlahd/dictionary"
Outputdir="${root}/CD_analysis"

#Iterate over the files
for i in $root/$1*1.fastq.gz; do
    #Extract the files basenames
    if [[ $i =~ ($1[[:digit:]]+)_(.+)\.fastq\.gz ]]; then
        bs_number="${BASH_REMATCH[1]}"
        #Uncompres the files into fastqtemp directory
        gunzip -c "${root}/${bs_number}_1.fastq.gz" > "${root}/fastqtmp/${bs_number}_1.fastq"
        gunzip -c "${root}/${bs_number}_2.fastq.gz" > "${root}/fastqtmp/${bs_number}_2.fastq"
        #run the hla code
        hlahd.sh -t 20 -m 100 -c 0.95 -f "${freq_data}" "${root}/fastqtmp/${bs_number}_1.fastq" "${root}/fastqtmp/${bs_number}_2.fastq" "${HLA_genesplit}" "${Dictionary}" "${bs_number}" "${Outputdir}" > "${root}/logs/${bs_number}.log" 2> "${root}/logs/${bs_number}.err"
        #check creation of results file and run code for table creation
        if [ -s "${Outputdir}/${bs_number}/result/${bs_number}_final.result.txt" ]; then
            #add time of completion to the log file
            echo "Analysis completed successfully at $(date +"%Y-%m-%d %H:%M:%S")" >> "${root}/logs/${bs_number}.log"
            #remove unziped files
            rm -f "${root}/fastqtmp/${bs_number}_1.fastq" "${root}/fastqtmp/${bs_number}_2.fastq"
        else #if the results file is empty or inexistent 
            #Add time
            echo "Analysis failed at $(date +"%Y-%m-%d %H:%M:%S")" >> "${root}/logs/${bs_number}.log"
            #add files size to log file
            echo "Uncompressed file sizes:" >> "${root}/logs/${bs_number}.log"
            du -h "${root}/fastqtmp/${bs_number}_1.fastq" "${root}/fastqtmp/${bs_number}_2.fastq" >> "${root}/logs/${bs_number}.log"
            #remove files
            rm -f "${root}/fastqtmp/${bs_number}_1.fastq" "${root}/fastqtmp/${bs_number}_2.fastq"
        fi
    fi
done 
