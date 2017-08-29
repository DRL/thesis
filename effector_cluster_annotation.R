# Analysis of effector clusters
library(readr)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(RColorBrewer)
library(ggthemes)
library(tools)
library(plyr)
library(stringr)
library(gtools)
library(gridExtra)
library(grid)
library(ggExtra)
library(scales)
library(cowplot)
read_delim_f <- function(filename){
  ret <- read_delim(file=filename, delim = "\t", escape_double = FALSE, trim_ws = TRUE,col_names = c('clusterID', 'proteins', 'proteomes', 'RFA'))
  ret
}
# IPR RFA
# EFFECTORS 
effector_clusters_import_f <- mixedsort(list.files(path=".", pattern="PCN_effectors.cluster_proteins_proteomes_RFA.tsv", all.files=T, full.names=T ))
effector_clusters_import <- data.frame(read_delim_f(effector_clusters_import_f))
effector_clusters_import$RFA <- factor(effector_clusters_import$RFA, levels = c('RFA', 'No RFA'))
# EFFECTORS RFA
effector_clusters_import_RFA <- effector_clusters_import[effector_clusters_import$RFA == "RFA",]
data <- effector_clusters_import_RFA
cluster_count <- paste("n =", length(data$clusterID))
colour = "orange"
rightbinwidth = 5
min_proteins = -5
max_proteins = max(data$proteins) + 20
min_proteomes = -1
max_proteomes = max(data$proteomes) + 2
RFA_hist_top <- ggplot(data, aes(x = proteomes)) + geom_histogram(aes(y=(..count..)/sum(..count..)), fill=colour, binwidth = 0.5, show.legend = F) + theme_minimal() + xlab("")  + theme(axis.text.x=element_blank(), plot.margin = unit(c(0.5, 0, 0, 0), "cm")) + xlim(min_proteomes, max_proteomes) + ylab("Proportion")
RFA_hist_right <- ggplot(data, aes(x = proteins)) + geom_histogram(aes(y=(..count..)/sum(..count..)),fill=colour, binwidth = rightbinwidth, show.legend = F) + theme_minimal() + xlab("") + coord_flip() + theme(axis.text.y=element_blank(), plot.margin = unit(c(0, 0.5, 0, 0), "cm"))  + ylab("Proportion") + xlim(min_proteins, max_proteins)
grob <- grobTree(textGrob(cluster_count, x=0.7,  y=0.95, hjust=0, gp=gpar(col ='grey3', fontsize=13)))
RFA_scatter <- ggplot(data) + geom_point(aes(x = proteomes, y = proteins), colour=colour, show.legend = F, alpha=0.25) + theme_minimal() + xlab("Proteomes") + ylab("Proteins") + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) + xlim(min_proteomes, max_proteomes) + ylim(min_proteins, max_proteins) + annotation_custom(grob)
RFA_hist_right_empty <- ggplot() + aes(1) + theme(axis.line.x = element_blank(),axis.line.y = element_blank(), axis.ticks=element_blank(), panel.background=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank()) + coord_flip()
p1 <- ggplot_gtable(ggplot_build(RFA_hist_top))
p2 <- ggplot_gtable(ggplot_build(RFA_hist_right))
p3 <- ggplot_gtable(ggplot_build(RFA_scatter))
p4 <- ggplot_gtable(ggplot_build(RFA_hist_right_empty))
first_row = plot_grid(p1, p4, rel_widths = c(2, 1), rel_heights = c(1, 1))
second_row = plot_grid(p3, p2, nrow = 1, rel_widths = c(2, 1), rel_heights = c(2, 1))
plot_grid(first_row, second_row, ncol=1, rel_widths = c(2, 2), rel_heights = c(1, 2))
# EFFECTORS NO RFA
effector_clusters_import_NORFA <- effector_clusters_import[effector_clusters_import$RFA == "No RFA",]
data <- effector_clusters_import_NORFA
cluster_count <- paste("n =", length(data$clusterID))
colour = "grey"
rightbinwidth = 5
#min_proteins = min(data$proteins) - 1.1
#max_proteins = max(data$proteins) + 5
#max_proteomes = max(data$proteomes) + 2
#min_proteomes = min(data$proteomes) - 1.1
RFA_hist_top <- ggplot(data, aes(x = proteomes)) + geom_histogram(aes(y=(..count..)/sum(..count..)), fill=colour, binwidth = 0.5, show.legend = F) + theme_minimal() + xlab("")  + theme(axis.text.x=element_blank(), plot.margin = unit(c(0.5, 0, 0, 0), "cm")) + xlim(min_proteomes, max_proteomes) + ylab("Proportion")
RFA_hist_right <- ggplot(data, aes(x = proteins)) + geom_histogram(aes(y=(..count..)/sum(..count..)),fill=colour, binwidth = rightbinwidth, show.legend = F) + theme_minimal() + xlab("") + coord_flip() + theme(axis.text.y=element_blank(), plot.margin = unit(c(0, 0.5, 0, 0), "cm"))  + ylab("Proportion") + xlim(min_proteins, max_proteins)
grob <- grobTree(textGrob(cluster_count, x=0.7,  y=0.95, hjust=0, gp=gpar(col ='grey3', fontsize=13)))
RFA_scatter <- ggplot(data) + geom_point(aes(x = proteomes, y = proteins), colour=colour, show.legend = F, alpha=0.25) + theme_minimal() + xlab("Proteomes") + ylab("Proteins") + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) + xlim(min_proteomes, max_proteomes) + ylim(min_proteins, max_proteins) + annotation_custom(grob)
RFA_hist_right_empty <- ggplot() + aes(1) + theme(axis.line.x = element_blank(),axis.line.y = element_blank(), axis.ticks=element_blank(), panel.background=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank()) + coord_flip()
p1 <- ggplot_gtable(ggplot_build(RFA_hist_top))
p2 <- ggplot_gtable(ggplot_build(RFA_hist_right))
p3 <- ggplot_gtable(ggplot_build(RFA_scatter))
p4 <- ggplot_gtable(ggplot_build(RFA_hist_right_empty))
first_row = plot_grid(p1, p4, rel_widths = c(2, 1), rel_heights = c(1, 1))
second_row = plot_grid(p3, p2, nrow = 1, rel_widths = c(2, 1), rel_heights = c(2, 1))
plot_grid(first_row, second_row, ncol=1, rel_widths = c(2, 2), rel_heights = c(1, 2))

# ALL clusters
all_clusters_import_f <- mixedsort(list.files(path=".", pattern="all_clusters.cluster_proteins_proteomes_RFA.tsv", all.files=T, full.names=T ))
all_clusters_import <- data.frame(read_delim_f(all_clusters_import_f))
all_clusters_import$RFA <- factor(all_clusters_import$RFA, levels = c('RFA', 'No RFA'))
summary(all_clusters_import)
# ALL RFA
all_clusters_import_RFA <- all_clusters_import[all_clusters_import$RFA == "RFA",]
data <- all_clusters_import_RFA
cluster_count <- paste("n =", length(data$clusterID))
colour = "orange"
rightbinwidth = 5
#min_proteins = min(data$proteins) - 1.1
#max_proteins = max(data$proteins) + 5
#min_proteomes = min(data$proteomes) - 1.1
#max_proteomes = max(data$proteomes) + 2
min_proteins = -5
max_proteins = max(data$proteins) + 20
min_proteomes = -1
max_proteomes = max(data$proteomes) + 2
RFA_hist_top <- ggplot(data, aes(x = proteomes)) + geom_histogram(aes(y=(..count..)/sum(..count..)), fill=colour, binwidth = 0.5, show.legend = F) + theme_minimal() + xlab("")  + theme(axis.text.x=element_blank(), plot.margin = unit(c(0.5, 0, 0, 0), "cm")) + xlim(min_proteomes, max_proteomes) + ylab("Proportion")
RFA_hist_right <- ggplot(data, aes(x = proteins)) + geom_histogram(aes(y=(..count..)/sum(..count..)),fill=colour, binwidth = rightbinwidth, show.legend = F) + theme_minimal() + xlab("") + coord_flip() + theme(axis.text.y=element_blank(), plot.margin = unit(c(0, 0.5, 0, 0), "cm"))  + ylab("Proportion") + xlim(min_proteins, max_proteins)
grob <- grobTree(textGrob(cluster_count, x=0.7,  y=0.95, hjust=0, gp=gpar(col ='grey3', fontsize=13)))
RFA_scatter <- ggplot(data) + geom_point(aes(x = proteomes, y = proteins), colour=colour, show.legend = F, alpha=0.25) + theme_minimal() + xlab("Proteomes") + ylab("Proteins") + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) + xlim(min_proteomes, max_proteomes) + ylim(min_proteins, max_proteins) + annotation_custom(grob)
RFA_hist_right_empty <- ggplot() + aes(1) + theme(axis.line.x = element_blank(),axis.line.y = element_blank(), axis.ticks=element_blank(), panel.background=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank()) + coord_flip()
p1 <- ggplot_gtable(ggplot_build(RFA_hist_top))
p2 <- ggplot_gtable(ggplot_build(RFA_hist_right))
p3 <- ggplot_gtable(ggplot_build(RFA_scatter))
p4 <- ggplot_gtable(ggplot_build(RFA_hist_right_empty))
first_row = plot_grid(p1, p4, rel_widths = c(2, 1), rel_heights = c(1, 1))
second_row = plot_grid(p3, p2, nrow = 1, rel_widths = c(2, 1), rel_heights = c(2, 1))
plot_grid(first_row, second_row, ncol=1, rel_widths = c(2, 2), rel_heights = c(1, 2))
# ALL NO RFA
all_clusters_import_NORFA <- all_clusters_import[all_clusters_import$RFA == "No RFA",]
data <- all_clusters_import_NORFA
cluster_count <- paste("n =", length(data$clusterID))
colour = "grey"
colour = "grey"
rightbinwidth = 5
#min_proteins = min(data$proteins) - 1.1
#max_proteins = max(data$proteins) + 5
#max_proteomes = max(data$proteomes) + 2
#min_proteomes = min(data$proteomes) - 1.1
RFA_hist_top <- ggplot(data, aes(x = proteomes)) + geom_histogram(aes(y=(..count..)/sum(..count..)), fill=colour, binwidth = 0.5, show.legend = F) + theme_minimal() + xlab("")  + theme(axis.text.x=element_blank(), plot.margin = unit(c(0.5, 0, 0, 0), "cm")) + xlim(min_proteomes, max_proteomes) + ylab("Proportion")
RFA_hist_right <- ggplot(data, aes(x = proteins)) + geom_histogram(aes(y=(..count..)/sum(..count..)),fill=colour, binwidth = rightbinwidth, show.legend = F) + theme_minimal() + xlab("") + coord_flip() + theme(axis.text.y=element_blank(), plot.margin = unit(c(0, 0.5, 0, 0), "cm"))  + ylab("Proportion") + xlim(min_proteins, max_proteins)
grob <- grobTree(textGrob(cluster_count, x=0.7,  y=0.95, hjust=0, gp=gpar(col ='grey3', fontsize=13)))
RFA_scatter <- ggplot(data) + geom_point(aes(x = proteomes, y = proteins), colour=colour, show.legend = F, alpha=0.25) + theme_minimal() + xlab("Proteomes") + ylab("Proteins") + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) + xlim(min_proteomes, max_proteomes) + ylim(min_proteins, max_proteins) + annotation_custom(grob)
RFA_hist_right_empty <- ggplot() + aes(1) + theme(axis.line.x = element_blank(),axis.line.y = element_blank(), axis.ticks=element_blank(), panel.background=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank()) + coord_flip()
p1 <- ggplot_gtable(ggplot_build(RFA_hist_top))
p2 <- ggplot_gtable(ggplot_build(RFA_hist_right))
p3 <- ggplot_gtable(ggplot_build(RFA_scatter))
p4 <- ggplot_gtable(ggplot_build(RFA_hist_right_empty))
first_row = plot_grid(p1, p4, rel_widths = c(2, 1), rel_heights = c(1, 1))
second_row = plot_grid(p3, p2, nrow = 1, rel_widths = c(2, 1), rel_heights = c(2, 1))
plot_grid(first_row, second_row, ncol=1, rel_widths = c(2, 2), rel_heights = c(1, 2))


# signalP RFA
# EFFECTORS 
effector_clusters_import_f <- mixedsort(list.files(path=".", pattern="PCN_effectors.cluster_proteins_proteomes_signalp.tsv", all.files=T, full.names=T ))
effector_clusters_import <- data.frame(read_delim_f(effector_clusters_import_f))
summary(effector_clusters_import)
effector_clusters_import$RFA <- factor(effector_clusters_import$RFA, levels = c('SignalP-TM;SignalP-noTM', 'SignalP-noTM', 'SignalP-noTM;SignalP-TM', 'None', 'SignalP-TM'))
# EFFECTORS RFA
effector_clusters_import_RFA <- effector_clusters_import[effector_clusters_import$RFA == "SignalP-noTM",]
effector_clusters_import_RFA <- rbind(effector_clusters_import_RFA, effector_clusters_import[effector_clusters_import$RFA == "SignalP-TM;SignalP-noTM",])
effector_clusters_import_RFA <- rbind(effector_clusters_import_RFA, effector_clusters_import[effector_clusters_import$RFA == "SignalP-noTM;SignalP-TM",])
data <- effector_clusters_import_RFA
cluster_count <- paste("n =", length(data$clusterID))
colour = "green3"
rightbinwidth = 5
min_proteins = -5
max_proteins = max(data$proteins) + 20
min_proteomes = -1
max_proteomes = max(data$proteomes) + 2
RFA_hist_top <- ggplot(data, aes(x = proteomes)) + geom_histogram(aes(y=(..count..)/sum(..count..)), fill=colour, binwidth = 0.5, show.legend = F) + theme_minimal() + xlab("")  + theme(axis.text.x=element_blank(), plot.margin = unit(c(0.5, 0, 0, 0), "cm")) + xlim(min_proteomes, max_proteomes) + ylab("Proportion")
RFA_hist_right <- ggplot(data, aes(x = proteins)) + geom_histogram(aes(y=(..count..)/sum(..count..)),fill=colour, binwidth = rightbinwidth, show.legend = F) + theme_minimal() + xlab("") + coord_flip() + theme(axis.text.y=element_blank(), plot.margin = unit(c(0, 0.5, 0, 0), "cm"))  + ylab("Proportion") + xlim(min_proteins, max_proteins)
grob <- grobTree(textGrob(cluster_count, x=0.7,  y=0.95, hjust=0, gp=gpar(col =colour, fontsize=13)))
RFA_scatter <- ggplot(data) + geom_point(aes(x = proteomes, y = proteins), colour=colour, show.legend = F, alpha=0.25) + theme_minimal() + xlab("Proteomes") + ylab("Proteins") + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) + xlim(min_proteomes, max_proteomes) + ylim(min_proteins, max_proteins) + annotation_custom(grob)
RFA_hist_right_empty <- ggplot() + aes(1) + theme(axis.line.x = element_blank(),axis.line.y = element_blank(), axis.ticks=element_blank(), panel.background=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank()) + coord_flip()
p1 <- ggplot_gtable(ggplot_build(RFA_hist_top))
p2 <- ggplot_gtable(ggplot_build(RFA_hist_right))
p3 <- ggplot_gtable(ggplot_build(RFA_scatter))
p4 <- ggplot_gtable(ggplot_build(RFA_hist_right_empty))
first_row = plot_grid(p1, p4, rel_widths = c(2, 1), rel_heights = c(1, 1))
second_row = plot_grid(p3, p2, nrow = 1, rel_widths = c(2, 1), rel_heights = c(2, 1))
plot_grid(first_row, second_row, ncol=1, rel_widths = c(2, 2), rel_heights = c(1, 2))
# EFFECTORS NO RFA
effector_clusters_import_NORFA <- effector_clusters_import[effector_clusters_import$RFA == "None",]
effector_clusters_import_NORFA <- rbind(effector_clusters_import_NORFA, effector_clusters_import[effector_clusters_import$RFA == "SignalP-TM",])
data <- effector_clusters_import_NORFA
cluster_count <- paste("n =", length(data$clusterID))
colour = "violet"
rightbinwidth = 5
#min_proteins = min(data$proteins) - 1.1
#max_proteins = max(data$proteins) + 5
#max_proteomes = max(data$proteomes) + 2
#min_proteomes = min(data$proteomes) - 1.1
RFA_hist_top <- ggplot(data, aes(x = proteomes)) + geom_histogram(aes(y=(..count..)/sum(..count..)), fill=colour, binwidth = 0.5, show.legend = F) + theme_minimal() + xlab("")  + theme(axis.text.x=element_blank(), plot.margin = unit(c(0.5, 0, 0, 0), "cm")) + xlim(min_proteomes, max_proteomes) + ylab("Proportion")
RFA_hist_right <- ggplot(data, aes(x = proteins)) + geom_histogram(aes(y=(..count..)/sum(..count..)),fill=colour, binwidth = rightbinwidth, show.legend = F) + theme_minimal() + xlab("") + coord_flip() + theme(axis.text.y=element_blank(), plot.margin = unit(c(0, 0.5, 0, 0), "cm"))  + ylab("Proportion") + xlim(min_proteins, max_proteins)
grob <- grobTree(textGrob(cluster_count, x=0.7,  y=0.95, hjust=0, gp=gpar(col =colour, fontsize=13)))
RFA_scatter <- ggplot(data) + geom_point(aes(x = proteomes, y = proteins), colour=colour, show.legend = F, alpha=0.25) + theme_minimal() + xlab("Proteomes") + ylab("Proteins") + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) + xlim(min_proteomes, max_proteomes) + ylim(min_proteins, max_proteins) + annotation_custom(grob)
RFA_hist_right_empty <- ggplot() + aes(1) + theme(axis.line.x = element_blank(),axis.line.y = element_blank(), axis.ticks=element_blank(), panel.background=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank()) + coord_flip()
p1 <- ggplot_gtable(ggplot_build(RFA_hist_top))
p2 <- ggplot_gtable(ggplot_build(RFA_hist_right))
p3 <- ggplot_gtable(ggplot_build(RFA_scatter))
p4 <- ggplot_gtable(ggplot_build(RFA_hist_right_empty))
first_row = plot_grid(p1, p4, rel_widths = c(2, 1), rel_heights = c(1, 1))
second_row = plot_grid(p3, p2, nrow = 1, rel_widths = c(2, 1), rel_heights = c(2, 1))
plot_grid(first_row, second_row, ncol=1, rel_widths = c(2, 2), rel_heights = c(1, 2))

# ALL clusters
all_clusters_import_f <- mixedsort(list.files(path=".", pattern="all_clusters.cluster_proteins_proteomes_signalp.tsv", all.files=T, full.names=T ))
all_clusters_import <- data.frame(read_delim_f(all_clusters_import_f))
all_clusters_import$RFA <- factor(all_clusters_import$RFA, levels = c('SignalP-TM;SignalP-noTM', 'SignalP-TM', 'SignalP-noTM', 'SignalP-noTM;SignalP-TM', 'None'))
summary(all_clusters_import)
# ALL RFA
all_clusters_import_RFA <- all_clusters_import[all_clusters_import$RFA == "SignalP-noTM",]
all_clusters_import_RFA <- rbind(all_clusters_import_RFA, all_clusters_import[all_clusters_import$RFA == "SignalP-TM;SignalP-noTM",])
all_clusters_import_RFA <- rbind(all_clusters_import_RFA, all_clusters_import[all_clusters_import$RFA == "SignalP-noTM;SignalP-TM",])
data <- all_clusters_import_RFA
cluster_count <- paste("n =", length(data$clusterID))
colour = "green3"
rightbinwidth = 5
#min_proteins = min(data$proteins) - 1.1
#max_proteins = max(data$proteins) + 5
#min_proteomes = min(data$proteomes) - 1.1
#max_proteomes = max(data$proteomes) + 2
min_proteins = -5
max_proteins = max(data$proteins) + 20
min_proteomes = -1
max_proteomes = max(data$proteomes) + 2
RFA_hist_top <- ggplot(data, aes(x = proteomes)) + geom_histogram(aes(y=(..count..)/sum(..count..)), fill=colour, binwidth = 0.5, show.legend = F) + theme_minimal() + xlab("")  + theme(axis.text.x=element_blank(), plot.margin = unit(c(0.5, 0, 0, 0), "cm")) + xlim(min_proteomes, max_proteomes) + ylab("Proportion")
RFA_hist_right <- ggplot(data, aes(x = proteins)) + geom_histogram(aes(y=(..count..)/sum(..count..)),fill=colour, binwidth = rightbinwidth, show.legend = F) + theme_minimal() + xlab("") + coord_flip() + theme(axis.text.y=element_blank(), plot.margin = unit(c(0, 0.5, 0, 0), "cm"))  + ylab("Proportion") + xlim(min_proteins, max_proteins)
grob <- grobTree(textGrob(cluster_count, x=0.7,  y=0.95, hjust=0, gp=gpar(col =colour, fontsize=13)))
RFA_scatter <- ggplot(data) + geom_point(aes(x = proteomes, y = proteins), colour=colour, show.legend = F, alpha=0.25) + theme_minimal() + xlab("Proteomes") + ylab("Proteins") + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) + xlim(min_proteomes, max_proteomes) + ylim(min_proteins, max_proteins) + annotation_custom(grob)
RFA_hist_right_empty <- ggplot() + aes(1) + theme(axis.line.x = element_blank(),axis.line.y = element_blank(), axis.ticks=element_blank(), panel.background=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank()) + coord_flip()
p1 <- ggplot_gtable(ggplot_build(RFA_hist_top))
p2 <- ggplot_gtable(ggplot_build(RFA_hist_right))
p3 <- ggplot_gtable(ggplot_build(RFA_scatter))
p4 <- ggplot_gtable(ggplot_build(RFA_hist_right_empty))
first_row = plot_grid(p1, p4, rel_widths = c(2, 1), rel_heights = c(1, 1))
second_row = plot_grid(p3, p2, nrow = 1, rel_widths = c(2, 1), rel_heights = c(2, 1))
plot_grid(first_row, second_row, ncol=1, rel_widths = c(2, 2), rel_heights = c(1, 2))
# ALL NO RFA
all_clusters_import_NORFA <- all_clusters_import[all_clusters_import$RFA == "SignalP-TM",]
all_clusters_import_NORFA <- rbind(all_clusters_import_NORFA, all_clusters_import[all_clusters_import$RFA == "None",])
data <- all_clusters_import_NORFA
cluster_count <- paste("n =", length(data$clusterID))
colour = "violet"
rightbinwidth = 5
#min_proteins = min(data$proteins) - 1.1
#max_proteins = max(data$proteins) + 5
#max_proteomes = max(data$proteomes) + 2
#min_proteomes = min(data$proteomes) - 1.1
RFA_hist_top <- ggplot(data, aes(x = proteomes)) + geom_histogram(aes(y=(..count..)/sum(..count..)), fill=colour, binwidth = 0.5, show.legend = F) + theme_minimal() + xlab("")  + theme(axis.text.x=element_blank(), plot.margin = unit(c(0.5, 0, 0, 0), "cm")) + xlim(min_proteomes, max_proteomes) + ylab("Proportion")
RFA_hist_right <- ggplot(data, aes(x = proteins)) + geom_histogram(aes(y=(..count..)/sum(..count..)),fill=colour, binwidth = rightbinwidth, show.legend = F) + theme_minimal() + xlab("") + coord_flip() + theme(axis.text.y=element_blank(), plot.margin = unit(c(0, 0.5, 0, 0), "cm"))  + ylab("Proportion") + xlim(min_proteins, max_proteins)
grob <- grobTree(textGrob(cluster_count, x=0.7,  y=0.95, hjust=0, gp=gpar(col =colour, fontsize=13)))
RFA_scatter <- ggplot(data) + geom_point(aes(x = proteomes, y = proteins), colour=colour, show.legend = F, alpha=0.25) + theme_minimal() + xlab("Proteomes") + ylab("Proteins") + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) + xlim(min_proteomes, max_proteomes) + ylim(min_proteins, max_proteins) + annotation_custom(grob)
RFA_hist_right_empty <- ggplot() + aes(1) + theme(axis.line.x = element_blank(),axis.line.y = element_blank(), axis.ticks=element_blank(), panel.background=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank()) + coord_flip()
p1 <- ggplot_gtable(ggplot_build(RFA_hist_top))
p2 <- ggplot_gtable(ggplot_build(RFA_hist_right))
p3 <- ggplot_gtable(ggplot_build(RFA_scatter))
p4 <- ggplot_gtable(ggplot_build(RFA_hist_right_empty))
first_row = plot_grid(p1, p4, rel_widths = c(2, 1), rel_heights = c(1, 1))
second_row = plot_grid(p3, p2, nrow = 1, rel_widths = c(2, 1), rel_heights = c(2, 1))
plot_grid(first_row, second_row, ncol=1, rel_widths = c(2, 2), rel_heights = c(1, 2))