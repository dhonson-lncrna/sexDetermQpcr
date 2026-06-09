# Analysis and data visualization of sex-determination qPCR for the Hawaiian bobtail squid

Note: The code, primer sequences, and methodology are the intellectual property of Drew Honson. Please acknowledge Dr. Honson in any publication using any portion of this repository. 

## Introduction

This repository contains two files for analyzing qPCR data to determine the sex of juvenile bobtail squid. 

1. An excel spreadsheet for matching samples and primers to a 384-well plate.
2. A python script for analyzing the data output from a QuantStudio 5 qPCR thermocycler. The script can easily be adapted for other machines by changing the collect_cts() function to match your instrument's output.

The protocol uses four sets of primers for each animal, two to the Z chromosome and two to autosomes. Cephalopods are diploid and use a ZZ/ZO sex determination system. As such, assuming perfect efficiency, sex can be determined from the following formula:

$$Ct_{chrZ} - Ct_{chrA} = \Delta Ct$$

$$\text{if } \Delta Ct \approx 0 \rightarrow Z:A = 1:1 \rightarrow \text{female}$$

$$\text{if } \Delta Ct \approx 1 \rightarrow Z:A = 1:2 \rightarrow \text{male}$$

The primers used for the protocol are:

| Name | Alias | Sequence | Product Length | Target | Efficiency |
|------|-------|----------|---------------|--------|------------|
| chrZ\_9\_FWD | chrZ\_1\_FWD | AACCAGGTCTTTCCCATGCC | 106 | LOC145268276 | 0.97 |
| chrZ\_9\_RVS | chrZ\_1\_RVS | GCTGGTGGACTCTTGCTAGG | | | |
| chrZ\_10\_FWD | chrZ\_2\_FWD | TCTGTGTGAACCCTTACGCA | 106 | LOC145268176 | 1.01 |
| chrZ\_10\_RVS | chrZ\_2\_RVS | TGAATCGGAGAGTAAATCCCTCA | | | |
| chr20\_2\_FWD | chrA\_1\_FWD | CTTTTCGCTTGCTTCGTGCT | 143 | LOC145241333 | 0.96 |
| chr20\_2\_RVS | chrA\_1\_RVS | CAGACCTTCACTAGCGCCTT | | | |
| chr37\_1\_FWD | chrA\_2\_FWD | GACCTTGCACCGTATCGTGG | 143 | LOC145261195 | 1.01 |
| chr37\_1\_RVS | chrA\_2\_RVS | TCCATCTGGGTATGCTGGCT | | | |

## Requirements

- Python >= 3.12
- polars
- fastexcel
- numpy
- matplotlib
- seaborn

## Usage

Fill out the 384-well plate layout sheets for samples and primers. The maximum number of samples is 24 if doing three replicate wells in a triangle, 32 if filling the entire plate using an adjustable width multichannel. Perform the qPCR on genomic DNA (1-10 ng per well, optimally ~5 ng) using SYBR Green chemistry.

Once the data have been exported, the basic execution of the script is:

'''
python analyzeQpcr.py --layout /path/to/layout.xlsx --results /path/to/results.xlsx
'''

Optionally, the script accepts two other arguments:

'''
--thresh / -t Maximum acceptable standard deviation of Ct replicates.
              Wells exceeding this threshold will have their worst outlier 
              removed before re-evaluation. If the deviation remains above 
              the threshold, the sample will be discarded and recorded in 
              "<prefix>_failedWells.log". (default: 0.1)

--prefix / -p 'Prefix for output files (<prefix>_failedWells.log, 
              <prefix>_results.csv, <prefix>_qpcrResults.png). Defaults to 
              the current date (YYYYMMDD). 
''' 
