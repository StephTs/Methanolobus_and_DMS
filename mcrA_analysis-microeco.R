##### Microeco - Exploration of the Baltic _mcrA_ amplicon data from QIIME2 #####
# Author: Stephania L. Tsola
# Date: September 2022
# R version 4.2.1
# RStudio 2022.07.1

library(file2meco)
library(microeco)
library(qiime2R)
library(phyloseq)
library(ggplot2)

#If unable to install qiime2R:
#download.file("https://github.com/jbisanz/qiime2R/Methhive/master.zip", "source.zip")
#unzip("source.zip")
#install.packages("qiime2R-master", repos = NULL, type="source")

# Covert the qiime2 files to phyloseq objects
BaltmcrAphyseq<-qza_to_phyloseq(
  features="mcrA-table.qza",
  taxonomy="mcrA-taxonomy.qza",
  metadata = "mcrA-metadata.tsv"
)
BaltmcrAphyseq

# Insert phyloseq object to microeco
BaltMeth_mec <- phyloseq2meco(BaltmcrAphyseq)
BaltMeth_mec


# Analysis starts here
BaltMeth_mec$cal_abund() # Abundance of mcrA although for the paper I used excel for the relative abundance graphs due to inexperience
# save abundance to a directory
BaltMeth_mec$save_abund(dirpath = "BaltMeth_taxa_abund")
BaltMeth_mec$cal_alphadiv(measure = "Shannon") #Calculate alpha diversity using the Shannon index
# save alpha diversity to a directory
BaltMeth_mec$save_alphadiv(dirpath = "BaltMeth_alpha_diversity")
# unifrac = FALSE means do not calculate unifrac metric - requires GUniFrac package installed
BaltMeth_mec$cal_betadiv(unifrac = FALSE)
# save beta diversity to a directory
BaltMeth_mec$save_betadiv(dirpath = "BaltMeth_beta_diversity")

# create trans_abund object
Meth_t1 <- trans_alpha$new(dataset = BaltMeth_mec, group = "Location")
Meth_t1$data_stat
Meth_t1$cal_diff(method = "anova", anova_set = "Treatment+Location+Depth")
Meth_t1$res_diff

#$Shannon
#    Df        Sum   Sq     Mean     Sq   F value Pr(>F)    
#  Treatment     1 12.726  12.726  45.470 1.05e-09 ***
#  Location    2  4.060   2.030   7.254  0.00115 **    
#  Depth       6  1.992   0.332   1.186  0.31992    
#  Residuals   99 27.708   0.280                       
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


Meth_t1$cal_diff(method = "t.test")
Meth_t1$res_diff
#    Comparison Measure  Test_method  Group  P.unadj    P.adj        Significance
#  1    H2 - H3 Shannon      t.test    H3 0.170235279 0.170235279           ns
#  2    H2 - H5 Shannon      t.test    H5 0.003222301 0.009666904           **
#  3    H3 - H5 Shannon      t.test    H5 0.030424309 0.045636463            *


Meth_t2 <- trans_alpha$new(dataset = BaltMeth_mec, group = "Treatment")
Meth_t2$cal_diff(method = "t.test")
Meth_t2$res_diff
#Comparison Measure Test_method    Group      P.unadj        P.adj Significance
#1 Original - DMS Shannon      t.test Original 7.618871e-09 7.618871e-09          ***

Meth_t3 <- trans_alpha$new(dataset = BaltMeth_mec, group = "Depth") 
Meth_t3$cal_diff(method = "KW_dunn")
Meth_t3$res_diff
Meth_t3$plot_alpha(pair_compare = TRUE, measure = "Shannon", shape = "Treatment")



# Methanogens grouped by Treatment
t1_BaltMeth_treatment <- trans_beta$new(dataset = BaltMeth_mec, group = "Treatment", measure = "bray")
t1_BaltMeth_treatment$cal_manova(manova_all = TRUE)
t1_BaltMeth_treatment$res_manova
#Permutation test for adonis under reduced model
#Terms added sequentially (first to last)
#Permutation: free
#Number of permutations: 999

#adonis2(formula = use_formula, data = metadata)
#Df SumOfSqs      R2      F Pr(>F)    
#Treatment  1    7.629 0.17702 23.015  0.001 ***
#  Residual  107   35.469 0.82298                    
#Total     108   43.099 1.00000                
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

t1_BaltMeth_treatment$cal_manova(manova_all = FALSE) #manova_all = FALSE can be used to calculate significance for each paired group
t1_BaltMeth_treatment$res_manova
#         Groups           measure  F        R2      p.value  p.adjusted    Significance
#  1 Original vs DMS    bray 23.01502 0.1770182   0.001      0.001          ***


t1_BaltMeth_treatment$cal_betadisper()
t1_BaltMeth_treatment$res_betadisper
#Permutation test for homogeneity of multivariate dispersions
#Permutation: free
#Number of permutations: 999
#
#Response: Distances
#           Df  Sum       Sq    Mean      Sq   F N.Perm  Pr(>F)    
#Groups     1 0.04729 0.047293 2.9224    999  0.107
#Residuals 107 1.73156 0.016183            
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Pairwise comparisons:
#  (Observed p-value below diagonal, permuted p-value above diagonal)
#           DMS    Original
#DMS                0.111
#Original 0.090255          

# Methanogens grouped by Location
t1_BaltMeth_location <- trans_beta$new(dataset = BaltMeth_mec, group = "Location", measure = "bray")
t1_BaltMeth_location$cal_manova(manova_all = TRUE)
t1_BaltMeth_location$res_manova
#Permutation test for adonis under reduced model
#Terms added sequentially (first to last)
#Permutation: free
#Number of permutations: 999

#adonis2(formula = use_formula, data = metadata)
#Df SumOfSqs      R2      F Pr(>F)   
#Location   2    4.471 0.10373 6.1339  0.001 ***
#  Residual 106   38.628 0.89627                         
#Total    108   43.099 1.00000                    
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

t1_BaltMeth_location$cal_manova(manova_all = FALSE) #manova_all = FALSE can be used to calculate significance for each paired group
t1_BaltMeth_location$res_manova
#Groups measure        F         R2 p.value p.adjusted Significance
#  1 H2 vs H3    bray 9.574689 0.12185463   0.001      0.001          ***
#  2 H2 vs H5    bray 3.436835 0.04496308   0.001      0.001          ***
#  3 H3 vs H5    bray 5.916329 0.07793223   0.001      0.001          ***

# for the whole comparison and for each paired groups
t1_BaltMeth_location$cal_betadisper()
t1_BaltMeth_location$res_betadisper 
#Permutation test for homogeneity of multivariate dispersions
#Permutation: free
#Number of permutations: 999
#
#Response: Distances
#           Df    Sum Sq   Mean Sq      F       N.Perm    Pr(>F)
#Groups     2 0.27922 0.139609 7.2592    999  0.001 ***
#Residuals 106 2.03859 0.019232                   
#
#Pairwise comparisons:
#  (Observed p-value below diagonal, permuted p-value above diagonal)
#H2         H3    H5
#H2            0.0100000 0.744
#H3 0.0101104           0.005
#H5 0.7511369 0.0036268      

# Methanogens grouped by Depth
t1_BaltMeth_depth <- trans_beta$new(dataset = BaltMeth_mec, group = "Depth", measure = "bray")
t1_BaltMeth_depth$cal_manova(manova_all = TRUE) #Grouped by Location
t1_BaltMeth_depth$res_manova
#Permutation test for adonis under reduced model
#Terms added sequentially (first to last)
#Permutation: free
#Number of permutations: 999

#adonis2(formula = use_formula, data = metadata)
#Df SumOfSqs      R2      F Pr(>F)   
#Location   6    7.959 0.18467 3.8504  0.001 ***
#  Residual 102   35.140 0.81533                  
#Total    108   43.099 1.00000                           
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

t1_BaltMeth_depth$cal_manova(manova_all = FALSE) #manova_all = FALSE can be used to calculate significance for each paired group
t1_BaltMeth_depth$res_manova
#Groups measure        F         R2 p.value  p.adjusted Significance
#  1  D1 vs D2    bray 4.522379 0.18441844   0.004 0.007636364           **
#  2  D1 vs D3    bray 1.320003 0.05015207   0.237 0.237000000             
#  3  D1 vs D4    bray 2.627566 0.08053208   0.030 0.035000000            *
#  4  D1 vs D5    bray 4.714598 0.13983847   0.001 0.002333333           **
#  5  D1 vs D6    bray 4.519849 0.13484096   0.001 0.002333333           **
#  6  D1 vs D7    bray 3.412507 0.10528366   0.008 0.011200000            *
#  7  D2 vs D3    bray 3.883297 0.15606038   0.007 0.011200000            *
#  8  D2 vs D4    bray 5.004851 0.16142153   0.002 0.004200000           **
#  9  D2 vs D5    bray 7.084768 0.22081407   0.001 0.002333333           **
#  10 D2 vs D6    bray 6.657991 0.21030997   0.001 0.002333333           **
#  11 D2 vs D7    bray 6.811338 0.21411667   0.001 0.002333333           **
#  12 D3 vs D4    bray 2.224556 0.06695517   0.047 0.051947368             
#  13 D3 vs D5    bray 5.095638 0.14519293   0.001 0.002333333           **
#  14 D3 vs D6    bray 4.862903 0.13948646   0.001 0.002333333           **
#  15 D3 vs D7    bray 3.212245 0.09671870   0.008 0.011200000            *
#  16 D4 vs D5    bray 4.331662 0.11013167   0.001 0.002333333           **
#  17 D4 vs D6    bray 4.224549 0.10770167   0.001 0.002333333           **
#  18 D4 vs D7    bray 2.790359 0.07383785   0.010 0.013125000            *
#  19 D5 vs D6    bray 1.750818 0.04897281   0.082 0.086100000             
#  20 D5 vs D7    bray 2.371411 0.06519985   0.022 0.027176471            *
#  21 D6 vs D7    bray 2.823897 0.07668654   0.007 0.011200000            *

t1_BaltMeth_depth$cal_betadisper()
t1_BaltMeth_depth$res_betadisper 
#Permutation test for homogeneity of multivariate dispersions
#Permutation: free
#Number of permutations: 999
#
#Response: Distances
#           Df    Sum Sq   Mean Sq      F       N.Perm    Pr(>F)
#Groups      6 0.2419 0.040316 7.2322    999  0.001 ***
#Residuals 102 0.5686 0.005575                              
#
#Pairwise comparisons:
#  (Observed p-value below diagonal, permuted p-value above diagonal)
#D1         D2         D3         D4         D5         D6    D7
#D1            1.0000e-03 1.0900e-01 3.1200e-01 7.0000e-02 2.3100e-01 0.262
#D2 9.3779e-04            1.0000e-03 2.0000e-03 1.0000e-03 1.0000e-03 0.001
#D3 1.0279e-01 1.0807e-04            9.8600e-01 7.4600e-01 9.2000e-01 0.901
#D4 3.1545e-01 4.2368e-04 9.7601e-01            8.7800e-01 9.1500e-01 0.942
#D5 6.6647e-02 2.0067e-05 7.6160e-01 8.6135e-01            7.1300e-01 0.948
#D6 2.0062e-01 7.0304e-05 9.2194e-01 9.1769e-01 7.0669e-01            0.843
#D7 2.4478e-01 2.7854e-04 8.9716e-01 9.3386e-01 9.3588e-01 8.3606e-01      


# PCoA plot showing the treatment as colours and the location as the point shape

t1_BaltMeth_treatment$cal_ordination(ordination = "PCoA")

# The code below exports the PCo-axes scores into a txt file that I can then import into the program PAST to run the Spearman correlation
# PAST: https://palaeo-electronica.org/2001_1/past/issue1_01.htm

BaltMeth_PCoA_scores_treatment <- t1_BaltMeth_treatment$res_ordination$scores
BaltMeth_PCoA_scores_treatment

write.table(BaltMeth_PCoA_scores_treatment, file = "BaltMeth_PCoA_scores_treatment.txt", sep = "\t",
            row.names = TRUE, col.names = NA)

# Making the PCoA plot using ggplot
BaltMeth_PCoA_treat_colTreat <- t1_BaltMeth_treatment$plot_ordination(plot_color = "Treatment",
                                                                    plot_shape = "Location", 
                                                                    plot_type = c("point", "ellipse"))
BaltMeth_PCoA_treat_colTreat + theme_bw() +
  theme(text = element_text(size = 21, face="bold", colour = "black"),
        axis.text.x = element_text(face="bold", colour = "black"), 
        axis.text.y = element_text(face="bold", colour = "black"))+ 
  geom_point(size = 4)

