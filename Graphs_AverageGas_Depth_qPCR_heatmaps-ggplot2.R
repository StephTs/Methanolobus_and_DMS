##### Graphs - Making the various graphs in the Baltic Sea paper #####
# Author: Stephania L. Tsola
# Date: September 2022

library("dplyr")
library("ggplot2")

# Graph that shows the total Methane and CO2 produced and DMS degraded during the incubation experiments per depth in each sampling station
DMS_graph_facet <- ggplot(Methane_CO2_DMS, aes(x=Depth, y=Mean, fill=Gas)) +
  geom_bar(stat="identity", colour="black", position="dodge")+
  geom_errorbar(aes(ymin=Mean-se, ymax=Mean+se), size=.3, width=.2, position=position_dodge(.9))+ 
  labs (y= "Gas concentration (µmol/g)", x="") +
  guides(fill=guide_legend(title=""))+
  theme(text=element_text(size=22, face = "bold", colour = "black"),
        axis.text.x = element_text(face="bold", colour = "black"), 
        axis.text.y = element_text(face="bold", colour = "black")) + 
  facet_wrap(~ Site, scales = "free_x")

DMS_graph_facet


#Graph showing the concentration of methane and sulfate per depth in the three sampling stations

ggplot(Baltic_2019_t0)+
  #set up the aesthetics, the line must be broken because this is technically discrete data
  geom_line(data=Baltic_2019_t0[!is.na(Baltic_2019_t0$Sulfate_mM), ],aes(x=Depth, y=Sulfate_mM, color="blue"),size=1, linetype=1)+
  geom_point(aes(x=Depth,y=Sulfate_mM))+
  geom_line(data=Baltic_2019_t0[!is.na(Baltic_2019_t0$Methane_µM), ], aes(x=Depth, y=Methane_µM, color="red"),size=1, linetype=1)+
  geom_point(aes(x=Depth,y=Methane_µM))+
  #reverse depth so it starts at zero and goes up by 10
  scale_x_reverse(breaks=pretty(Baltic_2019_t0$Depth, n=10))+
  #put the y axis labs on the opposite side so when its flipped it will appear at top
  scale_y_continuous(position="right", limits = c(0, 8), breaks = seq(0, 8, by = 1))+
  #this is how you reverse the look and order of the coordinates for the graph
  coord_flip()+
  theme(axis.ticks=element_line(linewidth=0.6),
        panel.border=element_blank(),
        axis.text.x = element_text(face="bold", size=11, colour = "black"),
        axis.text.y = element_text(face="bold", size=11, colour = "black"),
        axis.title = element_text(size = 14),
        legend.position = "none")+
  labs (y= expression(bold("Methane (µM) and Sulfate (mM) concentrations")), x=expression(bold("Depth (cmbsf)")))+
  facet_wrap(~Station)


# Graph showing the qPCR results per station and depth

## First import csv which contains Gene_Copy, Replicate, Treatment, Depth, Site

##This gives you the average and standard error between the triplicates
mcrA_sum_SD <- mcrA_Baltic_qPCR %>%
  group_by(Sample) %>%
  summarise(Mean= mean(Gene_Copy, na.rm=T), se = sd(Gene_Copy, na.rm=T)/sqrt(3), Replicate, Treatment, Depth, Site)

write.csv(mcrA_sum_SD,"mcrA_sum_SD.csv", row.names = FALSE) #Export summary to csv so can change table if needed

mcrA_qPCR_graph <- ggplot(mcrA_sum_SD, aes(x=Depth, y=Mean, fill=Treatment)) +
  geom_bar(stat="identity", colour="black", position="dodge")+
  geom_errorbar(aes(ymin=Mean-se, ymax=Mean+se), size=.3, width=.2, position=position_dodge(.9)) +
  labs (y= expression(bold(~bolditalic("mcrA")~ "copies/g")), x="") +
  guides(fill=guide_legend(title=""))+
  theme(text=element_text(size=22, face = "bold", colour = "black"),
        axis.text.x = element_text(face="bold", colour = "black"),
        axis.text.y = element_text(face="bold", colour = "black"),
        panel.grid.minor = element_blank()) +
  facet_wrap(~ Site, scales = "free_x")
mcrA_qPCR_graph

# Heatmaps showing the normalised data from the metagenomics (metaG) and metatranscriptomics (metaT) sequencing as well as the _Methanolobus_ MAGS

heatmap_MetaT_genes_log_normalised <- ggplot(data= Baltic_MetaT_Genes, mapping = aes(x=Sample ,y= Genes, fill = Log_FPKM)) +
  geom_tile(color="white", size=0.25) +
  labs(x="", y="")+
  scale_y_discrete(expand=c(0,0),limits = rev)+
  theme_grey(base_size=8)+
  theme(legend.text=element_text(size=14),
        legend.title=element_text(size=14, face="bold"),
        axis.ticks=element_line(linewidth =0.2),
        panel.border=element_blank(),
        axis.text.x = element_text(face="bold", size=14, angle = 90, hjust = 0, vjust = 0.5),
        axis.text.y = element_text(face="bold", size=14))+
  scale_fill_gradient(low = "darkorange", high = "blue3")
heatmap_MetaT_genes_log_normalised

heatmap_MetaG_genes_log_normalised <- ggplot(data= Baltic_MetaG_Genes, mapping = aes(x=Sample ,y= Genes, fill = Log_CPM)) +
  geom_tile(color="white", size=0.25) +
  labs(x="", y="")+
  scale_y_discrete(expand=c(0,0),limits = rev)+
  theme_grey(base_size=8)+
  theme(legend.text=element_text(size=14),
        legend.title=element_text(size=14, face="bold"),
        axis.ticks=element_line(linewidth =0.2),
        panel.border=element_blank(),
        axis.text.x = element_text(face="bold", size=14, angle = 90, hjust = 0, vjust = 0.5),
        axis.text.y = element_text(face="bold", size=14))+
  scale_fill_gradient(low = "darkorange", high = "blue3")
heatmap_MetaG_genes_log_normalised

heatmap_Methanolobus_MAG <- ggplot(data= MAGs_Methanolobus, mapping = aes(x=Sample ,y= Genes, fill = Presence)) +
  geom_tile(color="white", size=0.25) +
  coord_equal(ratio=0.8)+
  labs(x="", y="")+
  scale_y_discrete(expand=c(0,0),limits = rev)+
  theme_grey(base_size=8)+
  theme(legend.text=element_text(size=10),
        legend.title=element_text(size=10, face="bold"),
        axis.ticks=element_line(size=0.2),
        panel.border=element_blank(),
        axis.text.x = element_text(face="bold", size=14, angle = 90, hjust = 0, vjust = 0.5),
        axis.text.y = element_text(face="bold", size=14))+
  scale_fill_gradient(low = "grey63", high = "blue3")



