## Author: Dr Susan Corley through bioinformatics consultaton
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("gplots")
BiocManager::install("ggplot2")
BiocManager::install("dplyr")
BiocManager::install("rvg")
BiocManager::install("officer")
BiocManager::install("ggrepel")
BiocManager::install("tibble")
BiocManager::install("openxlsx")

library("DESeq2")
library("gplots")
library("ggplot2")
library("dplyr")
library("rvg")
library("officer")
library("ggrepel")
library("tibble")
library("openxlsx")


### Import input data (Input data: featureCounts output)

cts <- as.data.frame(read.table("featureCounts", header = TRUE))
cts <- cts[, -c(2:6)]  # remove unnecessary columns
cts <- column_to_rownames(cts, var = "Geneid")



### Import sample information
coldata <- read.table("coldata.txt", header = T, sep = "\t")
coldata <- column_to_rownames(coldata, var = "sample")
colnames(cts) <- rownames(coldata)
all(rownames(coldata) == colnames(cts))  # Check whether the sample order of coldata and cts is correct.



### DEG analysis
dds <- DESeqDataSetFromMatrix(countData = round(cts),
                              colData = coldata,
                              design = ~ genotype)
dds
keep <- rowSums(counts(dds)) >= 10
dds_keep <- dds[keep, ]
dds_keep
dds_keep$genotype <- relevel(dds_keep$genotype, ref = "WT")
dds_keep <- DESeq(dds_keep)

result <- results(dds_keep,
                        contrast = c("genotype", "GW", "WT"),
                        alpha = 0.05)
summary(result)

result_merge <- merge(as.data.frame(result),
                      as.data.frame(counts(dds_keep, normalized = TRUE)),
                      by = "row.names",
                      sort = FALSE)
rownames(result_merge) <- result_merge[,1]
result_merge <- result_merge[,c(2:13)]

result_merge <- merge(result_merge,
                      cts,
                      by = "row.names",
                      sort = FALSE)

names(result_merge)[1] <- "Gene"

write.xlsx(result_merge, file = "DESeq2_result.xlsx",
           colNames = TRUE, rowNames = FALSE)



### Make PCA plot
vsd <- vst(dds_keep, blind = FALSE)
zscore <- assay(vsd)
zscore <- t(scale(t(zscore)))
zscore <- as.data.frame(zscore)

z = plotPCA(vsd, intgroup = c("genotype"))
z = z + geom_label(aes(label = name))

pcaData <- plotPCA(vsd,
                   intgroup = c("genotype"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color="genotype", shape = "cell")) +
  geom_point(size=2) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed(ratio = 1) +
  theme_bw()
plotPCA(vsd, intgroup=c("genotype"))-> pca

PCAplot <- read_pptx() %>%
  add_slide() %>% ph_with(dml(ggobj = pca), location = ph_location_type(type = "body"))
print(PCAplot, target = "DESeq2_PCA.pptx")



### Make volcano plot
data <- result[,c("Gene", "log2FoldChange", "padj")]
dim(data)

data <- data %>% 
  mutate(
    Expression = case_when(log2FoldChange >= log2(2) & padj < 0.05 ~ "up-regulated",
                           log2FoldChange <= -log2(2) & padj < 0.05 ~ "down-regulated",
                           TRUE ~ "Unchanged")
  )
head(data)

data %>% 
  count(Expression)

ggplot(data, aes(log2FoldChange, -log(padj, 10))) +
  geom_point(aes(color = Expression), size = 0.5) +
  xlab(expression("log"[2]*"FC")) +
  ylab(expression("-log"[10]*"adj.p")) +
  scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
  guides(colour = guide_legend(override.aes = list(size = 1.5))) +
  theme_bw() -> volcano

# add Gene symbol in the volcano plot
top <- 10
top_logFC <- bind_rows(
  data %>% 
    filter(Expression == 'up-regulated') %>% 
    arrange(desc(log2FoldChange), padj) %>% 
    head(top),
  data %>% 
    filter(Expression == 'down-regulated') %>% 
    arrange(log2FoldChange, padj) %>% 
    head(top)
)

top_adjp <- bind_rows(
  data %>% 
    filter(Expression == 'up-regulated') %>% 
    arrange(padj) %>% 
    head(top),
  data %>% 
    filter(Expression == 'down-regulated') %>% 
    arrange(padj) %>% 
    head(top)
)

top_genes <- rbind(top_logFC, top_adjp)
top_genes <- distinct(top_genes, Gene, .keep_all = T)

volcano +
  geom_text_repel(data = top_genes,
                   mapping = aes(log2FoldChange, -log(padj, 10), label = Gene),
                   size = 3) -> volcano_2

volcanoplot <- read_pptx() %>%
  add_slide() %>% ph_with(dml(ggobj = volcano), location = ph_location_type(type = "body")) %>%
  add_slide() %>% ph_with(dml(ggobj = volcano_2), location = ph_location_type(type = "body"))

print(volcanoplot, target = "DESeq2_volcano.pptx")



### Make Heatmap
select_genes <- rownames(subset(result, padj < 0.05 & abs(log2FoldChange) > 1))

rld <- rlog(dds_keep)
rld_df <- as.data.frame(rld@assays@data)

heatmap.2(assay(rld)[select_genes,], 
          scale = "row",
          trace = "none",
          Colv = FALSE,
          Rowv = TRUE,
          density.info = "none",
          dendrogram = "both",
          col = greenred(75),
          ColSideColors = c(WT = "gray", GW = "black")[
            colData(rld)$genotype],
          labRow = FALSE)
