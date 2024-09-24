output_dir = "./"
table_count_dir = "../4.TableCount/"

library(stringr)
library(DESeq2)
library(EnhancedVolcano)
library(dplyr)

## DESeq analysis of NH4 starvation

# 1. Read data and process the data ####
table_count = read.table(paste0(table_count_dir, "Table_count.txt"), header=TRUE, row.names=1)
colnames(table_count)[6:11] = c("Control_1", "Control_2", "Control_3", "Starv_1", "Starv_2", "Starv_3")
row.names(table_count) = str_remove(row.names(table_count), "gene-")
table_count_short = table_count[, 6:ncol(table_count)]
table_count_short = as.matrix(table_count_short)

# 2. Condition assignment ####
condition = factor(c("Control", "Control", "Control", "Starv", "Starv", "Starv"))
coldata = data.frame(row.names=colnames(table_count_short), condition)
dds = DESeqDataSetFromMatrix(countData=table_count_short[rowSums(table_count_short) > 6, ], 
                             colData=coldata, 
                             design = ~condition)

# 3. Running DESeq ####
dds = DESeq(dds)
resultsNames(dds)  # "condition_Starv_vs_Control"

# PCA plot
vst = varianceStabilizingTransformation(dds)
png(paste0(output_dir, "PCA.png"), width=763, height=442)
plotPCA(vst, intgroup = c("condition"))
dev.off()

# 4. Result ####
res = results(dds, name="condition_Starv_vs_Control")
write.table(res, file=paste0(output_dir, "ResFull_Starv_vs_Control_stat.tsv"),
            sep="\t", col.names=NA, quote=FALSE)

res = as.data.frame(res)
res$gene = rownames(res)


#### Make a list with valuable genes (according to p.adj.value) ####
sorted_res = res[with(res, order(padj, -log2FoldChange)), ]
sorted_res_df = data.frame("id"=row.names(sorted_res), sorted_res)
p_adj_below_0.05 = sorted_res_df %>% filter(padj < 0.05)

### ______
# Starv vs Control => Starv / Control
# log2(100 / 10) # +  -> Starv > NH4
# log2(5 / 500) # -   -> Starv < NH4
### ______

Starv_activated = p_adj_below_0.05[p_adj_below_0.05$log2FoldChange > 0.5, ]  # Starv > Control
Starv_repressed = p_adj_below_0.05[p_adj_below_0.05$log2FoldChange < -0.5, ]  # Starv < Control

write.table(Starv_activated,
            file=paste0(output_dir, "Result_Starv_0.05_activated.tsv"),
            sep="\t", col.names=NA, quote=FALSE)

write.table(Starv_repressed,
            file=paste0(output_dir, "Result_Starv_0.05_repressed.tsv"),
            sep="\t", col.names=NA, quote=FALSE)

# Volcano plot ####
res_witout_na = na.omit(res)
png(paste0(output_dir, "Volcano.png"), width=1422, height=582)
EnhancedVolcano(res_witout_na,
                lab = rownames(res_witout_na),
                subtitle = "Nitrogen starvation",
                x = 'log2FoldChange',
                y = 'padj',
                FCcutoff = 0.5,
                pCutoff = 0.05,
                selectLab = NA,
)
dev.off()

# MA plot ####
png(paste0(output_dir, "MAPlot.png"), width=1532, height=730)
plotMA(na.omit(results(dds, name="condition_Starv_vs_Control")),
       ylim = c(-2.5, 6),
       alpha=0.05,
       colSig = "red",
       cex=0.6,
       colLine = "black",
       ylab="Log fold change",
       xlab="Mean of normalized counts"
)
dev.off()


# Check
table_count_short["PAS_chr1-1_0030", ]  # FoldChange "+" => Starv > Control
Starv_activated["PAS_chr1-1_0030",]
res["PAS_chr1-1_0030", ]

table_count_short["PAS_chr3_1071", ]  # FoldChange "-" => Starv < Control
Starv_repressed["PAS_chr3_1071",]
res["PAS_chr3_1071", ]
