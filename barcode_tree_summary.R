library(cowplot)
library(ape)
library(reshape2)
library(ggplot2)

# Set directory
# Microsoft
setwd('C:/Users/newUser/Dropbox/Barcoding review')
# Mac
setwd('~/Dropbox/Barcoding review')

##################################
### Tree distance calculations ###
##################################

rpoB = read.nexus('barcode_rpoB_MAFFT_GTR+G4+I_ML.nwk')
d.rpoB = cophenetic(rpoB)
d.rpoB[upper.tri(d.rpoB)] = NA
d.rpoB[which(d.rpoB==0)] = NA
write.csv(d.rpoB, 'barcode_rpoB_MAFFT_GTR+G4+I_ML_branch.csv')

ftsZ = read.nexus('barcode_ftsZ_MAFFT_GTR+G4+I_ML.nwk')
d.ftsZ = cophenetic(ftsZ)
d.ftsZ[upper.tri(d.ftsZ)] = NA
d.ftsZ[which(d.ftsZ==0)] = NA
write.csv(d.ftsZ, 'barcode_ftsZ_MAFFT_GTR+G4+I_ML_branch.csv')

gltA = read.nexus('barcode_gltA_MAFFT_GTR+G4+I_ML.nwk')
d.gltA = cophenetic(gltA)
d.gltA[upper.tri(d.gltA)] = NA
d.gltA[which(d.gltA==0)] = NA
write.csv(d.gltA, 'barcode_gltA_MAFFT_GTR+G4+I_ML_branch.csv')

ribC = read.nexus('barcode_ribC_MAFFT_GTR+G4+I_ML.nwk')
d.ribC = cophenetic(ribC)
d.ribC[upper.tri(d.ribC)] = NA
d.ribC[which(d.ribC==0)] = NA
write.csv(d.ribC, 'barcode_ribC_MAFFT_GTR+G4+I_ML_branch.csv')

gyrB = read.nexus('barcode_gyrB_MAFFT_GTR+G4+I_ML.nwk')
d.gyrB = cophenetic(gyrB)
d.gyrB[upper.tri(d.gyrB)] = NA
d.gyrB[which(d.gyrB==0)] = NA
write.csv(d.gyrB, 'barcode_gyrB_MAFFT_GTR+G4+I_ML_branch.csv')

groEL = read.nexus('barcode_groEL_MAFFT_GTR+G4+I_ML.nwk')
d.groEL = cophenetic(groEL)
d.groEL[upper.tri(d.groEL)] = NA
d.groEL[which(d.groEL==0)] = NA
write.csv(d.groEL, 'barcode_groEL_MAFFT_GTR+G4+I_ML_branch.csv')

nuoG = read.nexus('barcode_nuoG_MAFFT_GTR+G4+I_ML.nwk')
d.nuoG = cophenetic(nuoG)
d.nuoG[upper.tri(d.nuoG)] = NA
d.nuoG[which(d.nuoG==0)] = NA
write.csv(d.nuoG, 'barcode_nuoG_MAFFT_GTR+G4+I_ML_branch.csv')

ITS = read.nexus('barcode_ITS_MAFFT_GTR+G4+I_ML.nwk')
d.ITS = cophenetic(ITS)
d.ITS[upper.tri(d.ITS)] = NA
d.ITS[which(d.ITS==0)] = NA
write.csv(d.ITS, 'barcode_ITS_MAFFT_GTR+G4+I_ML_branch.csv')

x16SrRNA = read.nexus('barcode_16SrRNA_MAFFT_GTR+G4+I_ML.nwk')
d.16SrRNA = cophenetic(x16SrRNA)
d.16SrRNA[upper.tri(d.16SrRNA)] = NA
d.16SrRNA[which(d.16SrRNA==0)] = NA
write.csv(d.16SrRNA, 'barcode_16SrRNA_MAFFT_GTR+G4+I_ML_branch.csv')

x1R_2078C_membrane_protein = read.nexus('barcode_1R_2078C_MAFFT_GTR+G4+I_ML.nwk')
d.1R_2078C_membrane_protein = cophenetic(x1R_2078C_membrane_protein)
d.1R_2078C_membrane_protein[upper.tri(d.1R_2078C_membrane_protein)] = NA
d.1R_2078C_membrane_protein[which(d.1R_2078C_membrane_protein==0)] = NA
write.csv(d.1R_2078C_membrane_protein, 'barcode_1R_2078C_MAFFT_GTR+G4+I_ML_branch.csv')

x2R_2290C_nuclease = read.nexus('barcode_2R_2290C_MAFFT_GTR+G4+I_ML.nwk')
d.2R_2290C_nuclease = cophenetic(x2R_2290C_nuclease)
d.2R_2290C_nuclease[upper.tri(d.2R_2290C_nuclease)] = NA
d.2R_2290C_nuclease[which(d.2R_2290C_nuclease==0)] = NA
write.csv(d.2R_2290C_nuclease, 'barcode_2R_2290C_MAFFT_GTR+G4+I_ML_branch.csv')

x3R_781C = read.nexus('barcode_3R_781C_MAFFT_GTR+G4+I_ML.nwk')
d.3R_781C = cophenetic(x3R_781C)
d.3R_781C[upper.tri(d.3R_781C)] = NA
d.3R_781C[which(d.3R_781C==0)] = NA
write.csv(d.3R_781C, 'barcode_3R_781C_MAFFT_GTR+G4+I_ML_branch.csv')

x4R_3019C = read.nexus('barcode_4R_3019C_MAFFT_GTR+G4+I_ML.nwk')
d.4R_3019C = cophenetic(x4R_3019C)
d.4R_3019C[upper.tri(d.4R_3019C)] = NA
d.4R_3019C[which(d.4R_3019C==0)] = NA
write.csv(d.4R_3019C, 'barcode_4R_3019C_MAFFT_GTR+G4+I_ML_branch.csv')

x5R_675C_nuclease = read.nexus('barcode_5R_675C_MAFFT_GTR+G4+I_ML.nwk')
d.5R_675C_nuclease = cophenetic(x5R_675C_nuclease)
d.5R_675C_nuclease[upper.tri(d.5R_675C_nuclease)] = NA
d.5R_675C_nuclease[which(d.5R_675C_nuclease==0)] = NA
write.csv(d.5R_675C_nuclease, 'barcode_5R_675C_MAFFT_GTR+G4+I_ML_branch.csv')

x6R_3036C = read.nexus('barcode_6R_3036C_GTR+G4+I_ML.nwk')
d.6R_3036C = cophenetic(x6R_3036C)
d.6R_3036C[upper.tri(d.6R_3036C)] = NA
d.6R_3036C[which(d.6R_3036C==0)] = NA
write.csv(d.6R_3036C, 'barcode_6R_3036C_GTR+G4+I_ML_branch.csv')

x7R_1152C = read.nexus('barcode_7R_1152C_GTR+G4+I_ML.nwk')
d.7R_1152C = cophenetic(x7R_1152C)
d.7R_1152C[upper.tri(d.7R_1152C)] = NA
d.7R_1152C[which(d.7R_1152C==0)] = NA
write.csv(d.7R_1152C, 'barcode_7R_1152C_GTR+G4+I_ML_branch.csv')

x8R_2995C = read.nexus('barcode_8R_2995C_GTR+G4+I_ML.nwk')
d.8R_2995C = cophenetic(x8R_2995C)
d.8R_2995C[upper.tri(d.8R_2995C)] = NA
d.8R_2995C[which(d.8R_2995C==0)] = NA
write.csv(d.8R_2995C, 'barcode_8R_2995C_GTR+G4+I_ML_branch.csv')

x9R_3109C = read.nexus('barcode_9R_3109C_GTR+G4+I_ML.nwk')
d.9R_3109C = cophenetic(x9R_3109C)
d.9R_3109C[upper.tri(d.9R_3109C)] = NA
d.9R_3109C[which(d.9R_3109C==0)] = NA
write.csv(d.9R_3109C, 'barcode_9R_3109C_GTR+G4+I_ML_branch.csv')

x10R_1071C = read.nexus('barcode_10R_1071C_GTR+G4+I_ML.nwk')
d.10R_1071C = cophenetic(x10R_1071C)
d.10R_1071C[upper.tri(d.10R_1071C)] = NA
d.10R_1071C[which(d.10R_1071C==0)] = NA
write.csv(d.10R_1071C, 'barcode_10R_1071C_GTR+G4+I_ML_branch.csv')

x654R_2105C = read.nexus('barcode_654R_2105C_GTR+G4+I_ML.nwk')
d.654R_2105C = cophenetic(x654R_2105C)
d.654R_2105C[upper.tri(d.654R_2105C)] = NA
d.654R_2105C[which(d.654R_2105C==0)] = NA
write.csv(d.654R_2105C, 'barcode_654R_2105C_GTR+G4+I_ML_branch.csv')

x655R_4017C = read.nexus('barcode_655R_4017C_GTR+G4+I_ML.nwk')
d.655R_4017C = cophenetic(x655R_4017C)
d.655R_4017C[upper.tri(d.655R_4017C)] = NA
d.655R_4017C[which(d.655R_4017C==0)] = NA
write.csv(d.655R_4017C, 'barcode_655R_4017C_GTR+G4+I_ML_branch.csv')

x656R_2585C = read.nexus('barcode_656R_2585C_GTR+G4+I_ML.nwk')
d.656R_2585C = cophenetic(x656R_2585C)
d.656R_2585C[upper.tri(d.656R_2585C)] = NA
d.656R_2585C[which(d.656R_2585C==0)] = NA
write.csv(d.656R_2585C, 'barcode_656R_2585C_GTR+G4+I_ML_branch.csv')

x658R_1101C = read.nexus('barcode_658R_1101C_MAFFT_GTR+G4+I_ML.nwk')
d.658R_1101C = cophenetic(x658R_1101C)
d.658R_1101C[upper.tri(d.658R_1101C)] = NA
d.658R_1101C[which(d.658R_1101C==0)] = NA
write.csv(d.658R_1101C, 'barcode_658R_1101C_MAFFT_GTR+G4+I_ML_branch.csv')

x659R_1992C = read.nexus('barcode_659R_1992C_MAFFT_GTR+G4+I_ML.nwk')
d.659R_1992C = cophenetic(x659R_1992C)
d.659R_1992C[upper.tri(d.659R_1992C)] = NA
d.659R_1992C[which(d.659R_1992C==0)] = NA
write.csv(d.659R_1992C, 'barcode_659R_1992C_MAFFT_GTR+G4+I_ML_branch.csv')

x660R_3144C = read.nexus('barcode_660R_3144C_MAFFT_GTR+G4+I_ML.nwk')
d.660R_3144C = cophenetic(x660R_3144C)
d.660R_3144C[upper.tri(d.660R_3144C)] = NA
d.660R_3144C[which(d.660R_3144C==0)] = NA
write.csv(d.660R_3144C, 'barcode_660R_3144C_MAFFT_GTR+G4+I_ML_branch.csv')

x661R_3468C = read.nexus('barcode_661R_3468C_MAFFT_GTR+G4+I_ML.nwk')
d.661R_3468C = cophenetic(x661R_3468C)
d.661R_3468C[upper.tri(d.661R_3468C)] = NA
d.661R_3468C[which(d.661R_3468C==0)] = NA
write.csv(d.661R_3468C, 'barcode_661R_3468C_MAFFT_GTR+G4+I_ML_branch.csv')

x662R_2356C = read.nexus('barcode_662R_2356C_MAFFT_GTR+G4+I_ML.nwk')
d.662R_2356C = cophenetic(x662R_2356C)
d.662R_2356C[upper.tri(d.662R_2356C)] = NA
d.662R_2356C[which(d.662R_2356C==0)] = NA
write.csv(d.662R_2356C, 'barcode_662R_2356C_MAFFT_GTR+G4+I_ML_branch.csv')

x663R_4130C = read.nexus('barcode_663R_4130C_MAFFT_GTR+G4+I_ML.nwk')
d.663R_4130C = cophenetic(x663R_4130C)
d.663R_4130C[upper.tri(d.663R_4130C)] = NA
d.663R_4130C[which(d.663R_4130C==0)] = NA
write.csv(d.663R_4130C, 'barcode_663R_4130C_MAFFT_GTR+G4+I_ML_branch.csv')

###############################
### Figures and comparisons ###
###############################

# Import data for boxplots
sorted = read.csv('single_copy_clusters_sorted.csv', header=T)
group = c(rep('top10', 10), rep('remaining', 18))
sorted$group = group

# Create boxplots
a = ggplot(data=sorted, aes(x=factor(group), y=proportion.unique.32mers)) +
  geom_boxplot() +
  xlab('') +
  ylab('Proportion of unique 32-mers') +
  ylim(0, 1) +
  theme_classic(base_size=12) +
  theme(axis.line.x=element_blank(),
        axis.line.y=element_blank(),
        axis.ticks.x=element_line(colour='black'),
        axis.ticks.y=element_line(colour='black'),
        axis.title.x=element_text(colour='black'),
        axis.text.x=element_text(colour='black'),
        axis.title.y=element_text(colour='black'),
        axis.text.y=element_text(colour='black'),
        panel.border=element_rect(colour='black', fill=NA))

b = ggplot(data=sorted, aes(x=factor(group), y=proportion.segregating.sites)) +
  geom_boxplot() +
  xlab('') +
  ylab('Proportion of segregating sites') +
  ylim(0, 0.8) +
  theme_classic(base_size=12) +
  theme(axis.line.x=element_blank(),
        axis.line.y=element_blank(),
        axis.ticks.x=element_line(colour='black'),
        axis.ticks.y=element_line(colour='black'),
        axis.title.x=element_text(colour='black'),
        axis.text.x=element_text(colour='black'),
        axis.title.y=element_text(colour='black'),
        axis.text.y=element_text(colour='black'),
        panel.border=element_rect(colour='black', fill=NA))

c = ggplot(data=sorted, aes(x=factor(group), y=Watterson.estimator)) +
  geom_boxplot() +
  xlab('') +
  ylab('Watterson estimator') +
  ylim(0, 0.2) +
  theme_classic(base_size=12) +
  theme(axis.line.x=element_blank(),
        axis.line.y=element_blank(),
        axis.ticks.x=element_line(colour='black'),
        axis.ticks.y=element_line(colour='black'),
        axis.title.x=element_text(colour='black'),
        axis.text.x=element_text(colour='black'),
        axis.title.y=element_text(colour='black'),
        axis.text.y=element_text(colour='black'),
        panel.border=element_rect(colour='black', fill=NA))

d = ggplot(data=sorted, aes(x=factor(group), y=nucleotide.diversity)) +
  geom_boxplot() +
  xlab('') +
  ylab('Nucleotide diversity') +
  ylim(0, 0.3) +
  theme_classic(base_size=12) +
  theme(axis.line.x=element_blank(),
        axis.line.y=element_blank(),
        axis.ticks.x=element_line(colour='black'),
        axis.ticks.y=element_line(colour='black'),
        axis.title.x=element_text(colour='black'),
        axis.text.x=element_text(colour='black'),
        axis.title.y=element_text(colour='black'),
        axis.text.y=element_text(colour='black'),
        panel.border=element_rect(colour='black', fill=NA))

e = ggplot(data=sorted, aes(x=factor(group), y=Tamura.Nei.min)) +
  geom_boxplot() +
  xlab('') +
  ylab('Minimum Tamura-Nei distance') +
  ylim(0, 0.1) +
  theme_classic(base_size=12) +
  theme(axis.line.x=element_blank(),
        axis.line.y=element_blank(),
        axis.ticks.x=element_line(colour='black'),
        axis.ticks.y=element_line(colour='black'),
        axis.title.x=element_text(colour='black'),
        axis.text.x=element_text(colour='black'),
        axis.title.y=element_text(colour='black'),
        axis.text.y=element_text(colour='black'),
        panel.border=element_rect(colour='black', fill=NA))

f = ggplot(data=sorted, aes(x=factor(group), y=Tamura.Nei.median)) +
  geom_boxplot() +
  xlab('Group') +
  ylab('Median Tamura-Nei distance') +
  ylim(0, 0.5) +
  theme_classic(base_size=12) +
  theme(axis.line.x=element_blank(),
        axis.line.y=element_blank(),
        axis.ticks.x=element_line(colour='black'),
        axis.ticks.y=element_line(colour='black'),
        axis.title.x=element_text(colour='black'),
        axis.text.x=element_text(colour='black'),
        axis.title.y=element_text(colour='black'),
        axis.text.y=element_text(colour='black'),
        panel.border=element_rect(colour='black', fill=NA))

g = ggplot(data=sorted, aes(x=factor(group), y=Tamura.Nei.max)) +
  geom_boxplot() +
  xlab('') +
  ylab('Maximum Tamura-Nei distance') +
  ylim(0, 0.8) +
  theme_classic(base_size=12) +
  theme(axis.line.x=element_blank(),
        axis.line.y=element_blank(),
        axis.ticks.x=element_line(colour='black'),
        axis.ticks.y=element_line(colour='black'),
        axis.title.x=element_text(colour='black'),
        axis.text.x=element_text(colour='black'),
        axis.title.y=element_text(colour='black'),
        axis.text.y=element_text(colour='black'),
        panel.border=element_rect(colour='black', fill=NA))

# Combined plot and output to PDF
bottom = plot_grid(b, c, d, e, f, g, labels=c('b', 'c', 'd', 'e', 'f', 'g'), nrow=2, ncol=3)
pdf('barcode_boxplots.pdf', height=10, width=10)
plot_grid(a, bottom, labels=c('a', NA), nrow=2, ncol=1, rel_heights=c(0.33, 0.67))
dev.off()

# Import data for lineplots
sorted = read.csv('single_copy_clusters_sorted.csv', header=T)
msorted = melt(sorted, id.var='tag')

genes = c('1',
          '2',
          '3',
          '4',
          '5',
          '6',
          '7',
          '8',
          '9',
          '10',
          'ribC - 171',
          'gyrB - 386',
          'nuoG - 543',
          'gltA - 565',
          'ftsZ - 611',
          'rpoB - 635',
          'ITS - 652',
          '655',
          '656',
          '657',
          'groEL - 658',
          '659',
          '660',
          '661',
          '662',
          '663',
          '664',
          '16S rRNA - 665')

# Create lineplots
h = ggplot(msorted[1:112,], aes(x=tag, y=value, group=variable, colour=variable)) +
  geom_point() +
  geom_line() +
  scale_x_discrete(label=genes) +
  scale_colour_manual(name='Variables',
                      values=c('#e41a1c', '#377eb8', '#4daf4a', '#984ea3'),
                      labels=c('Proportion of\nunique 32-mers', 'Proportion of\nsegregating sites',
                               'Watterson estimator', 'Nucleotide diversity')) +
  ylab('Sequence diversity') +
  xlab('') +
  ylim(0, 1) +
  theme_classic(base_size=12) +
  theme(axis.line.x=element_blank(),
        axis.line.y=element_blank(),
        axis.ticks.x=element_line(colour='black'),
        axis.ticks.y=element_line(colour='black'),
        legend.title=element_text(size=10),
        legend.text=element_text(size=8),
        axis.title.x=element_text(colour='black'),
        axis.text.x=element_text(angle=60, hjust=1, size=8, colour='black'),
        axis.title.y=element_text(colour='black'),
        axis.text.y=element_text(size=8, colour='black'),
        panel.border=element_rect(colour='black', fill=NA))

i = ggplot(msorted[113:196,], aes(x=tag, y=value, group=variable, colour=variable)) +
  geom_point() +
  geom_line() +
  scale_x_discrete(label=genes) +
  scale_colour_manual(name='Variables',
                      values=c('#a65628', '#f781bf', '#999999'),
                      labels=c('Minimum Tamura-Nei\ndistance', 'Median Tamura-Nei\ndistance', 'Maximum Tamura-Nei\ndistance')) +
  ylab('Phylogenetic distance') +
  xlab('Ranking/Gene name') +
  ylim(0, 1) +
  theme_classic(base_size=12) +
  theme(axis.line.x=element_blank(),
        axis.line.y=element_blank(),
        axis.ticks.x=element_line(colour='black'),
        axis.ticks.y=element_line(colour='black'),
        legend.title=element_text(size=10),
        legend.text=element_text(size=8),
        axis.title.x=element_text(colour='black'),
        axis.text.x=element_text(angle=60, hjust=1, size=8, colour='black'),
        axis.title.y=element_text(colour='black'),
        axis.text.y=element_text(size=8, colour='black'),
        panel.border=element_rect(colour='black', fill=NA))

# Combined plot and output to PDF
pdf('barcode_lineplots.pdf', height=10, width=10)
plot_grid(h, i, labels=c('a', 'b'), nrow=2, ncol=1)
dev.off()

# Import data for marker frequency
markers = read.csv('Bartonella_Gene_Markers_sorted.csv', header=T)

# Create marker frequency barplot
j = ggplot(markers, aes(x=Rank, y=Publications.using.marker)) +
  geom_bar(stat = 'identity', colour=NA, fill='black') +
  scale_x_discrete(labels=markers$Gene.marker) +
  ylab('Publications using marker') +
  xlab('Marker') +
  ylim(0, 200) +
  theme_classic(base_size=12) +
  theme(axis.line.x=element_blank(),
        axis.line.y=element_blank(),
        axis.ticks.x=element_line(colour='black'),
        axis.ticks.y=element_line(colour='black'),
        axis.title.x=element_text(colour='black'),
        axis.text.x=element_text(angle=60, hjust=1, size=8, colour='black'),
        axis.title.y=element_text(colour='black'),
        axis.text.y=element_text(size=8, colour='black'),
        panel.border=element_rect(colour='black', fill=NA))

# Import and rearrange data for marker number scatterplot
numbers = read.csv('Bartonella_Gene_Markers_number.csv', header=T)
tab.num = as.data.frame(table(numbers))
tab.num.n0 = tab.num[which(tab.num$Freq>0),]
tab.num.n0$Year = as.numeric(as.character(tab.num.n0$Year))
tab.num.n0$Number.of.markers = as.numeric(as.character(tab.num.n0$Number.of.markers))
tab.num.hist = as.data.frame(table(numbers$Number.of.markers))
tab.num.hist$Var1 = as.numeric(as.character(tab.num.hist$Var1))
tab.num.hist$Freq = as.numeric(as.character(tab.num.hist$Freq))

# Create marker number scatterplot with trendline
k = ggplot(numbers, aes(x=Year, y=Number.of.markers)) +
  geom_smooth(method='lm', formula=y~x, colour='grey30', linetype=1, alpha=0) +
  geom_count(fill=NA, colour='black', alpha=0.8) +
  scale_size_area(name='Frequency') +
  scale_x_continuous(breaks=seq(1994, 2017, 1)) +
  ylab('Number of markers') +
  xlab('Year') +
  ylim(0, 15) +
  theme_classic(base_size=12) +
  theme(axis.line.x=element_blank(),
        axis.line.y=element_blank(),
        axis.ticks.x=element_line(colour='black'),
        axis.ticks.y=element_line(colour='black'),
        legend.title=element_text(size=10),
        legend.text=element_text(size=8),
        axis.title.x=element_text(colour='black'),
        axis.text.x=element_text(angle=60, hjust=1, size=8, colour='black'),
        axis.title.y=element_text(colour='black'),
        axis.text.y=element_text(size=8, colour='black'),
        panel.border=element_rect(colour='black', fill=NA))

l = ggplot(tab.num.hist, aes(x=Var1, y=Freq)) +
  geom_bar(stat='identity', colour=NA, fill='black') +
  scale_x_continuous(breaks=seq(0, 15, 1)) +
  ylab('Frequency') +
  xlab('Number of markers') +
  ylim(0, 150) +
  theme_classic(base_size=12) +
  theme(axis.line.x=element_blank(),
        axis.line.y=element_blank(),
        axis.ticks.x=element_line(colour='black'),
        axis.ticks.y=element_line(colour='black'),
        axis.title.x=element_text(colour='black'),
        axis.text.x=element_text(angle=0, size=8, colour='black'),
        axis.title.y=element_text(colour='black'),
        axis.text.y=element_text(size=8, colour='black'),
        panel.border=element_rect(colour='black', fill=NA))

# Combined plot and output to PDF
bottom = plot_grid(k, l, labels=c('b', 'c'), nrow=1, ncol=2, rel_widths=c(0.6, 0.4))
pdf('barcode_barplots.pdf', height=10, width=10)
plot_grid(j, bottom, labels=c('a',NA), nrow=2, ncol=1)
dev.off()

# Calculate trendline statistics with and without frequency
summary(lm(numbers$Number.of.markers~numbers$Year))
summary(lm(tab.num.n0$Number.of.markers~tab.num.n0$Year))