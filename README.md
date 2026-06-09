# Analysis and data visualization of sex-determination qPCR for the Hawaiian bobtail squid

This repository contains two files for analyzing qPCR data to determine the sex of juvenile bobtail squid. The protocol uses four sets of primers for each animal, two to the Z chromosome and two to autosomes. Cephalopods are diploid and use a ZZ/ZO sex determination system. As such, assuming perfect efficiency, sex can be determined from the following formula:

$$ Ct_{chrZ} - Ct_{chrA} &= \Delta Ct $$

$$ \text{if } \Delta Ct \approx 0 &\rightarrow Z:A = 1:1 \rightarrow \text{female} $$

$$ \text{if } \Delta Ct \approx 1 &\rightarrow Z:A = 1:2 \rightarrow \text{male} $$
