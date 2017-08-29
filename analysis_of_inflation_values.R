# Analysis of inflation values
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

read_delim_filename <- function(filename){
  ret <- read_delim(file=filename, delim = "\t", escape_double = FALSE, trim_ws = TRUE, na = 'None')
  print(ret)
  ret$iv <- sub("./Orthogroups_I(\\d+\\.\\d).*", "\\1", filename ,perl=TRUE)
  ret$params <- sub("./Orthogroups_I\\d+\\.\\d.only_cladeIV.longest_isoform.cluster_functional_annotation.all.p(\\d+)\\.x(\\d+).no_none.tsv", "p=0.\\1 x=0.\\2", filename ,perl=TRUE)
  ret$source <- filename
  ret
}

files <- mixedsort(list.files(path=".", pattern="*.no_none.tsv", all.files=T, full.names=T))
unsorted_import <- data.frame(ldply(files, read_delim_filename))
summary(unsorted_import)
import <- unsorted_import[order(unsorted_import$iv),] 
import$X.cluster_id <- NULL
import$protein_count <- NULL
import$domain_description <- NULL
import$source <- NULL
import$domain_ids <- NULL
xtabs(~iv+params+taxon_count, import)
colourCount = length(unique(import$iv))
getPalette = colorRampPalette(brewer.pal(9, "Paired"))
ggplot(as.data.frame(xtabs(~iv+params+taxon_count, import)), aes(x=taxon_count, y = Freq, colour = iv, group=iv)) + geom_line(size=1, show.legend = T, alpha=1.0) + facet_wrap(ncol = 3, ~params) + theme_minimal() + scale_color_manual(name = "Inflation value", limits = mixedsort(unique(import$iv)), values = getPalette(colourCount)) + ylab("Clusters with represenative functional annotation") + xlab('Number of taxa in cluster') + scale_y_log10() + annotation_logticks(side="l", size = 0.2, colour = "grey") + scale_x_discrete(breaks=c(1,5,10,15,20))#+ theme(axis.text.x = element_text(angle = 45, hjust = 1))

library(UpSetR)
read_delim_filename_taxon_count <- function(filename){
  ret <- read_delim(file=filename, delim = "\t", escape_double = FALSE, trim_ws = TRUE, na = 'None')
  ret
}
set_order = rev(c("SCARP", "SSCAP", "SGLAS", "SMONT", "SFELT", "RKR3021", "PTRIC", "SSTER", "SRATT", "SVENE", "SPAPI", "PREDI", "BXYLO", "GPALL", "GROST", "MHAPL", "MINCO", "MFLOR", "MAREN", "CELEG", "CBRIG"))
set_size <- 6

taxon_count_files <- mixedsort(list.files(path=".", pattern="Orthogroups_I1.5.only_cladeIV.longest_isoform.RFA.p0.x75.txt", all.files=T, full.names=T ))
taxon_count_import <- data.frame(ldply(taxon_count_files, read_delim_filename_taxon_count))
taxon_count_import[taxon_count_import > 1] <- 1
upset(taxon_count_import[rowSums(taxon_count_import)== set_size,], text.scale = 1, sets.x.label = "Occurence in clusters", mainbar.y.label = "Number of clusters", sets=set_order, mb.ratio = c(0.3, 0.7), nintersects = 10, group.by = "degree", keep.order = TRUE, order.by = c("freq"))

taxon_count_files <- mixedsort(list.files(path=".", pattern="Orthogroups_I4.0.only_cladeIV.longest_isoform.RFA.p0.x75.txt", all.files=T, full.names=T))
taxon_count_import <- data.frame(ldply(taxon_count_files, read_delim_filename_taxon_count))
taxon_count_import[taxon_count_import > 1] <- 1
upset(taxon_count_import[rowSums(taxon_count_import)== set_size,], text.scale = 1, sets.x.label = "Occurence in clusters", mainbar.y.label = "Number of clusters", sets=set_order, mb.ratio = c(0.3, 0.7), nintersects = 10, group.by = "degree", keep.order = TRUE, order.by = c("freq"))

taxon_count_files <- mixedsort(list.files(path=".", pattern="Orthogroups_I9.0.only_cladeIV.longest_isoform.RFA.p0.x75.txt", all.files=T, full.names=T))
taxon_count_import <- data.frame(ldply(taxon_count_files, read_delim_filename_taxon_count))
taxon_count_import[taxon_count_import > 1] <- 1
upset(taxon_count_import[rowSums(taxon_count_import)== set_size,], text.scale = 1, sets.x.label = "Occurence in clusters", mainbar.y.label = "Number of clusters", sets=set_order, mb.ratio = c(0.3, 0.7), nintersects = 10, group.by = "degree", keep.order = TRUE, order.by = c("freq"))

set_size <- 5

taxon_count_files <- mixedsort(list.files(path=".", pattern="Orthogroups_I1.5.only_cladeIV.longest_isoform.RFA.p0.x75.txt", all.files=T, full.names=T ))
taxon_count_import <- data.frame(ldply(taxon_count_files, read_delim_filename_taxon_count))
taxon_count_import[taxon_count_import > 1] <- 1
upset(taxon_count_import[rowSums(taxon_count_import)== set_size,], text.scale = 1, sets.x.label = "Occurence in clusters", mainbar.y.label = "Number of clusters", sets=set_order, mb.ratio = c(0.3, 0.7), nintersects = 10, group.by = "degree", keep.order = TRUE, order.by = c("freq"))

taxon_count_files <- mixedsort(list.files(path=".", pattern="Orthogroups_I4.0.only_cladeIV.longest_isoform.RFA.p0.x75.txt", all.files=T, full.names=T))
taxon_count_import <- data.frame(ldply(taxon_count_files, read_delim_filename_taxon_count))
taxon_count_import[taxon_count_import > 1] <- 1
upset(taxon_count_import[rowSums(taxon_count_import)== set_size,], text.scale = 1, sets.x.label = "Occurence in clusters", mainbar.y.label = "Number of clusters", sets=set_order, mb.ratio = c(0.3, 0.7), nintersects = 10, group.by = "degree", keep.order = TRUE, order.by = c("freq"))

taxon_count_files <- mixedsort(list.files(path=".", pattern="Orthogroups_I9.0.only_cladeIV.longest_isoform.RFA.p0.x75.txt", all.files=T, full.names=T))
taxon_count_import <- data.frame(ldply(taxon_count_files, read_delim_filename_taxon_count))
taxon_count_import[taxon_count_import > 1] <- 1
upset(taxon_count_import[rowSums(taxon_count_import)== set_size,], text.scale = 1, sets.x.label = "Occurence in clusters", mainbar.y.label = "Number of clusters", sets=set_order, mb.ratio = c(0.3, 0.7), nintersects = 10, group.by = "degree", keep.order = TRUE, order.by = c("freq"))

set_size <- 15

taxon_count_files <- mixedsort(list.files(path=".", pattern="Orthogroups_I1.5.only_cladeIV.longest_isoform.RFA.p0.x75.txt", all.files=T, full.names=T ))
taxon_count_import <- data.frame(ldply(taxon_count_files, read_delim_filename_taxon_count))
taxon_count_import[taxon_count_import > 1] <- 1
upset(taxon_count_import[rowSums(taxon_count_import)== set_size,], text.scale = 1, sets.x.label = "Occurence in clusters", mainbar.y.label = "Number of clusters", sets=set_order, mb.ratio = c(0.3, 0.7), nintersects = 10, group.by = "degree", keep.order = TRUE, order.by = c("freq"))

taxon_count_files <- mixedsort(list.files(path=".", pattern="Orthogroups_I4.0.only_cladeIV.longest_isoform.RFA.p0.x75.txt", all.files=T, full.names=T))
taxon_count_import <- data.frame(ldply(taxon_count_files, read_delim_filename_taxon_count))
taxon_count_import[taxon_count_import > 1] <- 1
upset(taxon_count_import[rowSums(taxon_count_import)== set_size,], text.scale = 1, sets.x.label = "Occurence in clusters", mainbar.y.label = "Number of clusters", sets=set_order, mb.ratio = c(0.3, 0.7), nintersects = 10, group.by = "degree", keep.order = TRUE, order.by = c("freq"))

taxon_count_files <- mixedsort(list.files(path=".", pattern="Orthogroups_I9.0.only_cladeIV.longest_isoform.RFA.p0.x75.txt", all.files=T, full.names=T))
taxon_count_import <- data.frame(ldply(taxon_count_files, read_delim_filename_taxon_count))
taxon_count_import[taxon_count_import > 1] <- 1
upset(taxon_count_import[rowSums(taxon_count_import)== set_size,], text.scale = 1, sets.x.label = "Occurence in clusters", mainbar.y.label = "Number of clusters", sets=set_order, mb.ratio = c(0.3, 0.7), nintersects = 10, group.by = "degree", keep.order = TRUE, order.by = c("freq"))

set_size <- 21

taxon_count_files <- mixedsort(list.files(path=".", pattern="Orthogroups_I1.5.only_cladeIV.longest_isoform.RFA.p0.x75.txt", all.files=T, full.names=T ))
taxon_count_import <- data.frame(ldply(taxon_count_files, read_delim_filename_taxon_count))
taxon_count_import[taxon_count_import > 1] <- 1
upset(taxon_count_import[rowSums(taxon_count_import)== set_size,], text.scale = 1, sets.x.label = "Occurence in clusters", mainbar.y.label = "Number of clusters", sets=set_order, mb.ratio = c(0.3, 0.7), nintersects = 10, group.by = "degree", keep.order = TRUE, order.by = c("freq"))

taxon_count_files <- mixedsort(list.files(path=".", pattern="Orthogroups_I4.0.only_cladeIV.longest_isoform.RFA.p0.x75.txt", all.files=T, full.names=T))
taxon_count_import <- data.frame(ldply(taxon_count_files, read_delim_filename_taxon_count))
taxon_count_import[taxon_count_import > 1] <- 1
upset(taxon_count_import[rowSums(taxon_count_import)== set_size,], text.scale = 1, sets.x.label = "Occurence in clusters", mainbar.y.label = "Number of clusters", sets=set_order, mb.ratio = c(0.3, 0.7), nintersects = 10, group.by = "degree", keep.order = TRUE, order.by = c("freq"))

taxon_count_files <- mixedsort(list.files(path=".", pattern="Orthogroups_I9.0.only_cladeIV.longest_isoform.RFA.p0.x75.txt", all.files=T, full.names=T))
taxon_count_import <- data.frame(ldply(taxon_count_files, read_delim_filename_taxon_count))
taxon_count_import[taxon_count_import > 1] <- 1
upset(taxon_count_import[rowSums(taxon_count_import)== set_size,], text.scale = 1, sets.x.label = "Occurence in clusters", mainbar.y.label = "Number of clusters", sets=set_order, mb.ratio = c(0.3, 0.7), nintersects = 10, group.by = "degree", keep.order = TRUE, order.by = c("freq"))
################################### all clusters
read_delim_cluster_counts <- function(filename){
  ret <- read_delim(file=filename, delim = "\t", escape_double = FALSE, trim_ws = TRUE, na = 'None')
  ret
}
set_order = rev(c("SCARP", "SSCAP", "SGLAS", "SMONT", "SFELT", "RKR3021", "PTRIC", "SSTER", "SRATT", "SVENE", "SPAPI", "PREDI", "BXYLO", "GPALL", "GROST", "MHAPL", "MINCO", "MFLOR", "MAREN", "CELEG", "CBRIG"))
set_size <- 19

taxon_count_import <- data.frame(ldply("../Orthogroups_I4.0.only_cladeIV.longest_isoform.kinfin_results/cluster_counts_by_taxon.txt", read_delim_cluster_counts))
taxon_count_import[taxon_count_import > 1] <- 1
taxon_count_import$X.ID <- NULL
upset(taxon_count_import[rowSums(taxon_count_import)== set_size,], text.scale = 1, sets.x.label = "Occurence in clusters", mainbar.y.label = "Number of clusters", sets=set_order, mb.ratio = c(0.3, 0.7), nintersects = 10, group.by = "degree", keep.order = TRUE, order.by = c("freq"))
