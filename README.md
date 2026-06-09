# Analysis and data visualization of sex-determination qPCR for the Hawaiian bobtail squid

This repository contains two files for analyzing qPCR data to determine the sex of juvenile bobtail squid. The protocol uses four sets of primers for each animal, two to the Z chromosome and two to autosomes. Cephalopods are diploid and use a ZZ/ZO sex determination system. As such, assuming perfect efficiency, sex can be determined from the following formula:

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

