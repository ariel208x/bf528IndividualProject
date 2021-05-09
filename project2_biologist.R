#libraries
library(tidyverse)
library(ggplot2)
library(RColorBrewer)

#load the data
fpkm.df <- read.table("/Users/arielxue/Documents/bf528/project5/genes.fpkm_tracking")
fpkm.matrix <- read.csv("/Users/arielxue/Documents/bf528/project5/fpkm_matrix.csv", sep =  "\t")
gene_exp_diff = read.table("/Users/arielxue/Documents/bf528/project5/gene_exp.diff")
gene_exp_diff <- gene_exp_diff %>% row_to_names(row_number = 1)
DGE.significant <- subset(gene_exp_diff,significant == "yes")
#rename the columns with names on first row
fpkm.df <- fpkm.df %>% row_to_names(row_number = 1)

#combine the two data frames together
fpkm.combined <- merge(fpkm.df, fpkm.matrix, by = 'tracking_id')
str(fpkm.combined)

#select needed columns as a new df
fpkm.selected <- fpkm.combined %>% select(gene_short_name,FPKM,Ad_1_FPKM,Ad_2_FPKM,P0_2_FPKM,P4_1_FPKM,P4_2_FPKM,P7_1_FPKM,P7_2_FPKM)

#change to numeric
fpkm.selected$FPKM <- as.numeric(fpkm.selected$FPKM)

###plot graphs for 7.1
##for SARCOMERE
gene_list_sarcomere = c('Pdlim5', 'Pygm', 'Myoz2', 'Des', 'Csrp3', 'Tcap', 'Cryab')
fpkm.selected.sarcomere <- fpkm.selected[fpkm.selected$gene_short_name %in% gene_list_sarcomere, ]

#get average
fpkm.selected.sarcomere$P0 <- rowMeans(fpkm.selected.sarcomere[ , c("FPKM","P0_2_FPKM")], na.rm=TRUE)
fpkm.selected.sarcomere$P4 <- rowMeans(fpkm.selected.sarcomere[ , c("P4_1_FPKM","P4_2_FPKM")], na.rm=TRUE)
fpkm.selected.sarcomere$P7 <- rowMeans(fpkm.selected.sarcomere[ , c("P7_1_FPKM","P7_2_FPKM")], na.rm=TRUE)
fpkm.selected.sarcomere$Ad <- rowMeans(fpkm.selected.sarcomere[ , c("Ad_1_FPKM","Ad_2_FPKM")], na.rm=TRUE)
fpkm.selected.sarcomere <- fpkm.selected.sarcomere %>% select(gene_short_name,P0,P4,P7,Ad)

#plot
fpkm.selected.sarcomere.plot <- fpkm.selected.sarcomere %>%
  select(gene_short_name,P0,P4,P7,Ad) %>%
  gather(key = "Time", value = "FPKM", -gene_short_name)

level_order <- factor(fpkm.selected.sarcomere.plot$Time, level = c('P0', 'P4', 'P7','Ad'))

ggplot(fpkm.selected.sarcomere.plot, aes(x = level_order, y = FPKM, group=gene_short_name),sort= False) +
  geom_line(aes(linetype=gene_short_name,color=gene_short_name))+
  geom_point(aes(shape=gene_short_name,color=gene_short_name))+ 
  theme_classic()+ 
  ggtitle("Sarcomere")+
  theme(axis.title.x=element_blank())

##for MITOCHONDRIA
gene_list_mitochondria = c("Mpc1", "Prdx3", "Acat1", "Echs1", "Slc25a11", "Phyh")
fpkm.selected.mitochondria <- fpkm.selected[fpkm.selected$gene_short_name %in% gene_list_mitochondria, ]

#get average
fpkm.selected.mitochondria$P0 <- rowMeans(fpkm.selected.mitochondria[ , c("FPKM","P0_2_FPKM")], na.rm=TRUE)
fpkm.selected.mitochondria$P4 <- rowMeans(fpkm.selected.mitochondria[ , c("P4_1_FPKM","P4_2_FPKM")], na.rm=TRUE)
fpkm.selected.mitochondria$P7 <- rowMeans(fpkm.selected.mitochondria[ , c("P7_1_FPKM","P7_2_FPKM")], na.rm=TRUE)
fpkm.selected.mitochondria$Ad <- rowMeans(fpkm.selected.mitochondria[ , c("Ad_1_FPKM","Ad_2_FPKM")], na.rm=TRUE)
fpkm.selected.mitochondria <- fpkm.selected.mitochondria %>% select(gene_short_name,P0,P4,P7,Ad)

#plot
fpkm.selected.mitochondria.plot <- fpkm.selected.mitochondria %>%
  select(gene_short_name,P0,P4,P7,Ad) %>%
  gather(key = "Time", value = "FPKM", -gene_short_name)

level_order.mitochondria <- factor(fpkm.selected.mitochondria.plot$Time, level = c('P0', 'P4', 'P7','Ad'))

ggplot(fpkm.selected.mitochondria.plot, aes(x = level_order.mitochondria, y = FPKM, group=gene_short_name),sort= False) +
  geom_line(aes(linetype=gene_short_name,color=gene_short_name))+
  geom_point(aes(shape=gene_short_name,color=gene_short_name))+ 
  theme_classic()+ 
  ggtitle("Mitochondria")+
  theme(axis.title.x=element_blank())

##for Cell_cycle
gene_list_Cell_cycle = c("Cdc7", "E2f8", "Cdk7", "Cdc26", "Cdc6", "E2f1", "Cdc27", "Bora", "Cdc45", "Rad51", "Aurkb", "Cdc23")
fpkm.selected.Cell_cycle <- fpkm.selected[fpkm.selected$gene_short_name %in% gene_list_Cell_cycle, ]

#get average
fpkm.selected.Cell_cycle$P0 <- rowMeans(fpkm.selected.Cell_cycle[ , c("FPKM","P0_2_FPKM")], na.rm=TRUE)
fpkm.selected.Cell_cycle$P4 <- rowMeans(fpkm.selected.Cell_cycle[ , c("P4_1_FPKM","P4_2_FPKM")], na.rm=TRUE)
fpkm.selected.Cell_cycle$P7 <- rowMeans(fpkm.selected.Cell_cycle[ , c("P7_1_FPKM","P7_2_FPKM")], na.rm=TRUE)
fpkm.selected.Cell_cycle$Ad <- rowMeans(fpkm.selected.Cell_cycle[ , c("Ad_1_FPKM","Ad_2_FPKM")], na.rm=TRUE)
fpkm.selected.Cell_cycle <- fpkm.selected.Cell_cycle %>% select(gene_short_name,P0,P4,P7,Ad)

#plot
fpkm.selected.Cell_cycle.plot <- fpkm.selected.Cell_cycle %>%
  select(gene_short_name,P0,P4,P7,Ad) %>%
  gather(key = "Time", value = "FPKM", -gene_short_name)

level_order.Cell_cycle <- factor(fpkm.selected.Cell_cycle.plot$Time, level = c('P0', 'P4', 'P7','Ad'))

ggplot(fpkm.selected.Cell_cycle.plot, aes(x = level_order.Cell_cycle, y = FPKM, group=gene_short_name),sort= False) +
  geom_line(aes(linetype=gene_short_name,color=gene_short_name))+
  geom_point(aes(shape=gene_short_name,color=gene_short_name))+ 
  theme_classic()+ 
  ggtitle("Cell Cycle")+
  theme(axis.title.x=element_blank())

###7.3
#get the combined matrix
fpkm.selected.matrix <- fpkm.combined %>% select(gene_short_name,FPKM,P0_2_FPKM,P4_1_FPKM,P4_2_FPKM,P7_1_FPKM,P7_2_FPKM,Ad_1_FPKM,Ad_2_FPKM)
fpkm.selected.matrix$FPKM <- as.numeric(fpkm.selected.matrix$FPKM)

#select top 1000 DEG
top1000 <- DGE.significant %>% top_n(n = 1000, wt = q_value)
fpkm.selected.matrix.top1000 <- fpkm.selected.matrix[fpkm.selected.matrix$gene_short_name %in% top1000$gene, ]

#get rid of dupication
fpkm.selected.matrix.top1000 <- fpkm.selected.matrix.top1000[!duplicated(fpkm.selected.matrix.top1000$gene_short_name), ]
rownames(fpkm.selected.matrix.top1000)<-fpkm.selected.matrix.top1000$gene_short_name

#make clustered heatmap with top1000 genes
fpkm.selected.matrix.top1000.fpkm <- fpkm.selected.matrix.top1000 %>% select(FPKM:Ad_2_FPKM)

fpkm.selected.matrix.top1000.fpkm <- fpkm.selected.matrix.top1000.fpkm %>% rename(P0_1_FPKM = FPKM)

fpkm.selected.matrix.top1000.fpkm <- data.matrix(fpkm.selected.matrix.top1000.fpkm)

fpkm.selected.matrix.top1000.fpkm <- fpkm.selected.matrix.top1000.fpkm[!(apply(fpkm.selected.matrix.top1000.fpkm, 1, function(y) any(y == 0))),]


pheatmap(fpkm.selected.matrix.top1000.fpkm, scale = "row", fontsize_row = 4, color = colorRampPalette(c("navy","yellow2"))(50), 
         cluster_row = T,
         cluster_col = T,
         show_rownames = F,
         show_colnames = T,
         legend = T,
         border_color = NA,
         fontsize = 10,
         cutree_rows = 3,
         column_title = "Samples")

