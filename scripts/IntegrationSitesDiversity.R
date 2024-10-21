working_dir = "~/git/htlv_integration_sites/data/"
setwd(working_dir)
dir()

sites = read.table(file="ALL.markedsorted_with_hits_to_viral_filtered_merged.bed.tsv", header=T, sep="\t")

sites$integration_site <- paste(sites$chromosome, sites$start_pos, sep = "_")

library(ggplot2)
library(RColorBrewer)
library(viridis)
library(scico)
library(gridExtra)

# some color palette options
# scale_fill_brewer(palette = "Paired") ++
# scale_fill_brewer(palette = "Dark2") ++
# scale_fill_brewer(palette = "Spectral")
# scale_fill_scico_d(palette = "hawaii")
# scale_fill_viridis_d(option = "inferno")

samples = unique(sites$sample)
b_adjust = -15
tf_size = 10
lw = 0.05
lsize = 2
  
#CTCF-7 plot
sites_sample = sites[which(sites$sample=="CTCF-7_SIC_934_SIC2"),]
sample_name = "CTCF-7"
gini_index = 0.69
n = 5 
o = order(sites_sample$count, decreasing = TRUE)
sites_sample = sites_sample[o,]
sites_sample$integration_site <- factor(sites_sample$integration_site, levels = sites_sample$integration_site[order(sites_sample$count, decreasing = TRUE)])
unique_sites = length(sites_sample$count)
total_count = sum(as.numeric(sites_sample$count))
sites_sample$label = ""
sites_sample$label[1:n] = sites_sample$count[1:n]
title = paste(sample_name, " (sites = ", unique_sites, ", count = ", total_count, ", GINI = ", gini_index, ")", sep="")
p1 = ggplot(data=sites_sample, aes(x = "", y = count, fill = integration_site)) +
  geom_bar(stat = "identity", width = 1, color = "black", size = lw) +        
  coord_polar("y") +                                                            
  theme_void() +                                                                
  theme(legend.position = "none") +                                             
  geom_text(aes(label=label), position=position_stack(vjust = 0.5), size = lsize) + 
  scale_fill_brewer(palette = "Spectral") +                                     
  ggtitle(title) +                                                              
  theme(plot.title = element_text(hjust = 0.5, size = tf_size, face = "bold", margin = margin(b = b_adjust)))
print(p1)

#CTCF-8 plot
sites_sample = sites[which(sites$sample=="CTCF-8_SIC_935_SIC2"),]
sample_name = "CTCF-8"
gini_index = 0.65
n = 5 
o = order(sites_sample$count, decreasing = TRUE)
sites_sample = sites_sample[o,]
sites_sample$integration_site <- factor(sites_sample$integration_site, levels = sites_sample$integration_site[order(sites_sample$count, decreasing = TRUE)])
unique_sites = length(sites_sample$count)
total_count = sum(as.numeric(sites_sample$count))
sites_sample$label = ""
sites_sample$label[1:n] = sites_sample$count[1:n]
title = paste(sample_name, " (sites = ", unique_sites, ", count = ", total_count, ", GINI = ", gini_index, ")", sep="")
p2 = ggplot(data=sites_sample, aes(x = "", y = count, fill = integration_site)) +
  geom_bar(stat = "identity", width = 1, color = "black", size = lw) +        
  coord_polar("y") +                                                            
  theme_void() +                                                                
  theme(legend.position = "none") +                                             
  geom_text(aes(label=label), position=position_stack(vjust = 0.5), size = lsize) + 
  scale_fill_brewer(palette = "Spectral") +                                     
  ggtitle(title) +                                                              
  theme(plot.title = element_text(hjust = 0.5, size = tf_size, face = "bold", margin = margin(b = b_adjust)))
print(p2)

#P12-10B plot
sites_sample = sites[which(sites$sample=="P12-10B_SIC_936_SIC2"),]
sample_name = "P12-10B"
gini_index = 0.50
n = 1 
o = order(sites_sample$count, decreasing = TRUE)
sites_sample = sites_sample[o,]
sites_sample$integration_site <- factor(sites_sample$integration_site, levels = sites_sample$integration_site[order(sites_sample$count, decreasing = TRUE)])
unique_sites = length(sites_sample$count)
total_count = sum(as.numeric(sites_sample$count))
sites_sample$label = ""
sites_sample$label[1:n] = sites_sample$count[1:n]
title = paste(sample_name, " (sites = ", unique_sites, ", count = ", total_count, ", GINI = ", gini_index, ")", sep="")
p3 = ggplot(data=sites_sample, aes(x = "", y = count, fill = integration_site)) +
  geom_bar(stat = "identity", width = 1, color = "black", size = lw) +        
  coord_polar("y") +                                                            
  theme_void() +                                                                
  theme(legend.position = "none") +                                             
  geom_text(aes(label=label), position=position_stack(vjust = 0.5), size = lsize) + 
  scale_fill_brewer(palette = "Spectral") +                                     
  ggtitle(title) +                                                              
  theme(plot.title = element_text(hjust = 0.5, size = tf_size, face = "bold", margin = margin(b = b_adjust)))
print(p3)

#P12-14 plot
sites_sample = sites[which(sites$sample=="P12-14_SIC_937_SIC2"),]
sample_name = "P12-14"
gini_index = 0.90
n = 3 
o = order(sites_sample$count, decreasing = TRUE)
sites_sample = sites_sample[o,]
sites_sample$integration_site <- factor(sites_sample$integration_site, levels = sites_sample$integration_site[order(sites_sample$count, decreasing = TRUE)])
unique_sites = length(sites_sample$count)
total_count = sum(as.numeric(sites_sample$count))
sites_sample$label = ""
sites_sample$label[1:n] = sites_sample$count[1:n]
title = paste(sample_name, " (sites = ", unique_sites, ", count = ", total_count, ", GINI = ", gini_index, ")", sep="")
p4 = ggplot(data=sites_sample, aes(x = "", y = count, fill = integration_site)) +
  geom_bar(stat = "identity", width = 1, color = "black", size = lw) +
  coord_polar("y") +                                                            
  theme_void() +                                                                
  theme(legend.position = "none") +                                             
  geom_text(aes(label=label), position=position_stack(vjust = 0.5), size = lsize) + 
  scale_fill_brewer(palette = "Spectral") +                                     
  ggtitle(title) +                                                              
  theme(plot.title = element_text(hjust = 0.5, size = tf_size, face = "bold", margin = margin(b = b_adjust)))
print(p4)

grid.arrange(p1, p2, p3, p4, nrow = 2, ncol = 2)


