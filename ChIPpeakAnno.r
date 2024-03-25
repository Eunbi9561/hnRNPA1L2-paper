## Author: Dr. Eunbi LEE
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ChIPpeakAnno")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
BiocManager::install("Guitar")
BiocManager::install("rvg")
BiocManager::install("officer")


library("ChIPpeakAnno")
library("org.Hs.eg.db")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
library("Guitar")
library("rvg")
library("officer")



### Preparing the annotation data
annoData <- toGRanges(TxDb.Hsapiens.UCSC.hg19.knownGene, feature = "gene")

### Importing the data
SampleA = read.delim("SampleA.narrowPeak", header = FALSE)
SampleB = read.delim("SampleB.narrowPeak", header = FALSE)
SampleC = read.delim("SampleC.narrowPeak", header = FALSE)

colnames(SampleA) <- c("seqnames", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak")
colnames(SampleB) <- c("seqnames", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak")
colnames(SampleC) <- c("seqnames", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak")

SampleA_GR <- toGRanges(SampleA, format = "narrowPeak")
SampleB_GR <- toGRanges(SampleB, format = "narrowPeak")
SampleC_GR <- toGRanges(SampleC, format = "narrowPeak")



### Finding overlaps
ol_all <- findOverlapsOfPeaks(SampleA_GR, SampleB_GR, SampleC_GR, 
                              connectedPeaks = "keepAll")
ol_all <- addMetadata(ol_all, colNames = "score", FUN = mean)
makeVennDiagram(ol_all, by = "feature",
                NameOfPeaks = c("SampleA", "SampleB", "SampleC"),
                main = "hg19",
                fill = c("#440154ff", '#21908dff', '#fde725ff'),
                col = c("#440154ff", '#21908dff', '#fde725ff'),
                cat.col = c("#000000", "#000000", "#000000"),
                cat.cex = 1.3,
                cex = 1.5)
ol_target_GR <- ol_all$peaklist[["SampleB_GR///SampleC_GR"]]
export.bed(ol_target_GR, "bed_overlap.bed") # Export bed file



### Genomic distribution
overlap_dist <- genomicElementDistribution(GRangesList(overlaps = ol_target_GR), 
                                           TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                                           promoterRegion = c(upstream = 2000, downstream = 500),
                                           geneDownstream = c(upstream = 0, downstream = 2000))
overlap_dist <- overlap_dist$plot$data

tiff("overlap_distribution.tiff",
     width = 1000, height = 700, units = "px", bg = "white")
genomicElementDistribution(GRangesList(overlaps = ol_target_GR), 
                           TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                           promoterRegion = c(upstream = 2000, downstream = 500),
                           geneDownstream = c(upstream = 0, downstream = 2000))
dev.off()



### Annotation
ol_target_GR_anno <- annotatePeakInBatch(ol_target_GR,
                                         AnnotationData = annoData,
                                         output = "nearestLocation")
ol_target_GR_anno <- addGeneIDs(ol_target_GR_anno,
                           orgAnn = "org.Hs.eg.db",
                           IDs2Add = c("symbol"),
                           feature_id_type = "entrez_id")
# Gene list
ol_target_GR_anno_df <- as.data.frame(unname(ol_target_GR_anno))
ol_target_GR_anno_df_freq <- data.frame(table(ol_target_GR_anno_df$symbol))
names(ol_target_GR_anno_df_freq)[1] <- "symbol"
ol_target_GR_anno_df_freq <- ol_target_GR_anno_df_freq[order(-ol_target_GR_anno_df_freq$Freq), ]
write.xlsx(ol_target_GR_anno_df_freq,
           "overlap_freq.xlsx")



### Guitar package
## Importing input bed file
input <- list("SampleA.narrowPeak",
              "SampleB.narrowPeak",
              "SampleC.narrowPeak",
              "bed_overlap.bed")

## Build Guitar Coordinate
txdb_file <- system.file("extdata", "hg19_toy.sqlite",
                         package="Guitar")
txdb <- loadDb(txdb_file)

## Plotting
GuitarPlot(txTxdb = txdb,
           stBedFiles = input,
           headOrtail = TRUE,
           enableCI = FALSE,
           mapFilterTranscript = TRUE,
           pltTxType = c("mrna"),
           stGroupName = c("SampleA",
                           "SampleB",
                           "SampleC",
                           "Overlap")) -> p1

plot <- read_pptx() %>%
  add_slide() %>% ph_with(dml(ggobj = p1), location = ph_location_type(type = "body")) 
print(plot, target = "density plot.pptx")
