# DMS degradation from _Methanolobus_
Code for the paper:
S. L. Tsola, Y. Zhu, Y. Chen, I. A. Sanders, C. K. Economou, V. Brüchert, Ö. Eyice (2023). _Methanolobus_ use unspecific methyltransferases to produce methane from dimethylsulfide. 
https://zenodo.org/badge/629950360.svg

This repo contains the code for analysing the amplicon data presented in the paper using QIIME2 as well as code for creating a custom _mcrA_ database. Furthermore, it has R scripts detailing the bioinformatic analysis of the amplicon data and how the various graphs were made using ggplot2. Lastly, program usage information is given for PAST (used for the Spearman correlation analysis).

Any help with improving the scripts or any other comments are greatly appreciated.  


## Content
1) mcrAseq_analysis-QIIME2.sh - Script for the qiime2 analysis of _mcrA_ sequences from the Baltic Sea samples
2) mcrA_sample-metadata.tsv - Metadata file of Baltic Sea _mcrA_ samples
3) mcrA_analysis-microeco.R - Script for microeco _mcrA_ analysis 
4) mcrAseq_Database - Script for creating the _mcrA_ database
5) Graphs_AverageGas_Depth_qPCR_heatmaps-ggplot2.R - Script containing information on the creation of all R graphs (avg gases, metag/T heatmap, qPCR)
7) Spearman_correlation-PAST.txt - Step-by-step guide to using PAST for the Spearman's correlation analysis 
8) Baltic_mcrA_PCoA_scores.csv - csv file imported into PAST  
9) LICENSE.md

## License
This project is licensed under the terms of the MIT License. See the LICENSE.md for further information.

