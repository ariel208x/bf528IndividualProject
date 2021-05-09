#libraries
library(janitor)

#load the data
gene_exp <- read.table("/Users/arielxue/Documents/bf528/project5/gene_exp.diff")
gene_exp_diff = read.table("/Users/arielxue/Documents/bf528/project5/gene_exp.diff")
str(gene_exp_diff)
#rename the columns with names on first row
gene_exp_diff <- gene_exp_diff %>% row_to_names(row_number = 1)

#sort the data with lowest q_value first
gene_exp_diff <- gene_exp_diff[order(gene_exp_diff$q_value),]
gene_exp_diff$value_1 <- as.numeric(gene_exp_diff$value_1)
gene_exp_diff$value_2 <- as.numeric(gene_exp_diff$value_2)
str(gene_exp_diff)
#create pseudo counts
gene_exp_diff$value_1.pseudo <- with(gene_exp_diff, value_1 + 1)
gene_exp_diff$value_2.pseudo <- with(gene_exp_diff, value_2 + 1)
gene_exp_diff$log2FC.pseudo <- with(gene_exp_diff, log2(value_2.pseudo/value_1.pseudo))

#make histogram for log2FC
DGE.log2FC <- as.numeric(gene_exp_diff$log2FC.pseudo)

hist(DGE.log2FC, breaks=50, main="Log2FC with pseudo-counts for all genes", xlab = "log2FC")

#df for genes with significant = yes
DGE.significant <- subset(gene_exp_diff,significant == "yes")

#second histogram for significant genes only
DGE.significant.log2FC <- as.numeric(DGE.significant$log2FC.pseudo)
hist(DGE.significant.log2FC, breaks=50, main="Log2FC with pseudo-counts for significant genes", xlab = "log2FC")

#genes with p-value < 0.01
significant <- subset(gene_exp, p_value < 0.01)
significant.up <- subset(significant,`log2(fold_change)` > 0)
significant.down <- subset(significant,`log2(fold_change)` < 0)

#sub df for up and down regulated genes
DGE.significant.up <- subset(DGE.significant, `log2(fold_change)` > 0)
DGE.significant.down <- subset(DGE.significant, `log2(fold_change)` < 0)

#write up and down regulated genes into seperate files
write.table(DGE.significant.up, file.path("/Users/arielxue/Documents/bf528/project5", "up_regulated_genes.txt"), row.names=F, quote=F)
write.table(DGE.significant.up$gene, file.path("/Users/arielxue/Documents/bf528/project5", "up_regulated_gene_names.txt"), col.names=F, row.names=F, quote=F)

write.table(DGE.significant.down, file.path("/Users/arielxue/Documents/bf528/project5", "down_regulated_genes.txt"), row.names=F, quote=F)
write.table(DGE.significant.down$gene, file.path("/Users/arielxue/Documents/bf528/project5", "down_regulated_gene_names.txt"), col.names=F, row.names=F, quote=F)




