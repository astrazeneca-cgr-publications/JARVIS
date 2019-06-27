library(ggplot2)
library(patchwork)
library(tidyverse)
library(ggridges)
library(scales)
library(RColorBrewer)


args = commandArgs(trailingOnly = TRUE)
results_dir = args[1]
out_file = args[2]
annot = args[3]

print(paste('results_dir:', results_dir))	
print(paste('out_file:', out_file))	
print(paste('annot:', annot))	



df <- read.csv(paste(results_dir, out_file, sep='/'), stringsAsFactors=FALSE)
print(head(df))

class_medians = sort(apply(df, 2, median, na.rm=TRUE))
print(class_medians)


functional_classes = names(class_medians) #c('intergenic', 'lincrna', 'tolerant', 'intron', 'mature_mirna', 'ccds', 'utr', 'pathogenic', 'vista', 'omim.HI', 'intolerant', 'ucne')

df = df[ , colnames(df) %in% functional_classes ]
df = df[ , match(functional_classes, colnames(df)) ]
print(tail(df))


#genomic_colors = c(brewer.pal(5, 'Set1'), brewer.pal(7, 'Set2')) 
genomic_colors = c( '#FC8D62', '#8DA0CB', '#E5C494', '#4DAF4A', '#984EA3', '#FF7F00', '#737373', '#A6D854', '#FFD92F', '#377EB8', '#E41A1C', '#E78AC3')
names(genomic_colors) = c('intergenic', 'lincrna', 'tolerant', 'intron', 'mature_mirna', 'ccds', 'utr', 'pathogenic', 'vista', 'omim.HI', 'intolerant', 'ucne')

#Melt data into column format.
melt_df <- gather(df, "variable", "value", 1:ncol(df))
print(head(melt_df))



melt_df$variable <- factor(melt_df$variable, colnames(df))
#melt_df = melt_df[ order(functional_classes, melt_df$variable), ]

# remove NA values
melt_df = melt_df[ !is.na(melt_df$value), ]

print(unique(melt_df$variable))


#Modify Theme:
source("z_theme.R")

# =======  BETA  ========
#melt_df = melt_df[ melt_df$value > 0,  ]
#melt_df$value = log(abs(melt_df$value))
#print(head(melt_df))
# =======================


# Boxplot
p1 = ggplot(melt_df,aes(variable,value))+
  geom_boxplot(aes(fill=variable), alpha=.8, na.rm=TRUE, outlier.colour="#535353", outlier.shape=16, outlier.size=1, notch=FALSE)+
  #geom_jitter(aes(color=variable),size=3,alpha=.2)+
  #scale_y_continuous(breaks=seq(0,1,.1), labels=scales::percent)+
  scale_fill_manual(values = genomic_colors )+
  guides(fill=FALSE,color=FALSE)+
  labs(title="Intolerance scores across different genomic classes",
       x="gwRVIS",
       y="Density")+
  coord_flip()+
  z_theme()
p1 = p1 + geom_vline(xintercept = 0, colour='#636363')
#ggsave(paste(results_dir, "/gwRVIS_joy_boxplot.", annot, ".pdf", sep=''), plot=p1, width=4, height=3, units="in", scale=3)
  

bw = 0.4



# Joyplot for probly
p2 = ggplot(melt_df,aes(y=variable,x=value))+
  geom_density_ridges(scale=15, aes(fill=variable), alpha=.8, na.rm=TRUE, bandwidth=bw)+
  #scale_x_continuous(breaks=seq(0,1,.1), labels=scales::percent)+
  #scale_fill_manual(values = c(brewer.pal(5, 'Set1'), brewer.pal(7, 'Set2')) )+
  scale_fill_manual(values = genomic_colors )+
  guides(fill=TRUE,color=FALSE)+
  labs(title="Intolerance scores across different genomic classes",
       y="Densities",
       x="gwRVIS") +
  z_theme()

#for(cl in unique(melt_df$variable)){
for(cl in c('intergenic', 'ucne')){
	cur_median = median(melt_df[ melt_df$variable == cl, 'value'])
	p2 = p2 + geom_vline(xintercept = cur_median, colour=genomic_colors[cl], linetype='dashed')
}
p2 = p2 + geom_vline(xintercept = 0, colour='#636363')
#ggsave(paste(results_dir, "/gwRVIS_density_ggjoy_plots.", annot, ".pdf", sep=''), plot=p2, width=4, height=3, units="in", scale=3)

p = p1 + p2
ggsave(paste(results_dir, "/gwRVIS_ggridges_Boxplot_n_Density_Plots.", annot, ".pdf", sep=''), plot=p, width=7, height=3, units="in", scale=3)
