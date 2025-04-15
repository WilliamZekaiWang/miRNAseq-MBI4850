library(DESeq2)
library(readr)
library(dplyr)

# load files
count_files <- list.files(pattern = "*_miRNA_counts.txt$", full.names = TRUE)

# Initialize an empty list to store data
count_list <- list()

# Read each file and store it in the list
for (file in count_files) {
  sample_name <- gsub("*_miRNA_counts.txt", "", basename(file))
  df <- read.delim(file, header=TRUE, row.names=1, comment.char = "#")
  count_list[[sample_name]] <- setNames(df[, ncol(df)], rownames(df))
}

# Combine into a single data frame
final_counts <- as.data.frame(do.call(cbind, count_list))
colnames(final_counts) <- names(count_list)  # Set column names to sample names

# remove low counts
min_count <- 10
min_samples <- ceiling(0.3 * ncol(final_counts))  # 30% of total samples

# keep only miRNAs with ≥10 counts in at least 30% of samples (exceRpt steps)
filtered_counts <- final_counts[rowSums(final_counts >= min_count) >= min_samples, ]

# read meta and change rownames to Runs
metadata <- read.csv("SraRunTable.csv", header=TRUE, row.names=1)
metadata$cell_type <- as.factor(metadata$cell_type)

# load into dds (already noranmlized internall)
dds <- DESeqDataSetFromMatrix(countData = filtered_counts, 
                              colData = metadata, # row names look good
                              design = ~ cell_type)  # mdd/not is specified in this 

# PCA
vsd <- varianceStabilizingTransformation(dds, blind=F)
pca <- prcomp(t(assay(vsd)))

pdf('scree.pdf')
plot(summary(pca)$importance[2,1:10]*100, ylab="% Variance Explained", xlab="PC", type="b")
dev.off()

# Extract PC1-PC5 which had vairance > ~5%
pca_data <- pca$x[, 1:5] 

# Add PCs to the colData of the DESeq2 object
colData(dds) <- cbind(colData(dds), pca_data)

# remove outlier samples based on PCA

# test for autocorrelation
pc_data <- as.data.frame(colData(dds)[, c("cell_type", paste0("PC", 1:5))])
pdf("pca.pdf")

twolines = function(x,y) {
  points(x,y,pch=20)
  abline(lm(y~x),col="red")
  #legend("bottomright", paste("R=",prettyNum(cor(x,y), digits=3)),bty="n" ,cex=1.5)
}

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y, use="complete.obs"))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}
mydiag.panel <- function( x,  labels, ...){
  ll <- par("usr")
  rect(ll[1], ll[3], ll[2], ll[4], col="darkolivegreen1")
}

diag.labels=c("MDD/Ctrl", "PC1", "PC2", "PC3", 'PC4','PC5')
plot.formula=as.formula(~cell_type + PC1 + PC2 + PC3 + PC4 + PC5)
pairs(plot.formula, data=pc_data, upper.panel=twolines, labels=diag.labels, diag.panel=mydiag.panel, lower.panel=panel.cor, label.pos=0.5, main="Correlation Between Variables")
dev.off()



# plot PCA
library(ggplot2)
pdf('pca.pdf', width = 5, height = 5)
ggplot(data = pc_data, aes(x = PC1, y = PC2, colour = cell_type)) +
  geom_point() +
  theme_minimal() +
  xlab("PC1 (15.33%)") +              
  ylab("PC2 (12.23%)") +             
  labs(colour = "Case")              # Changes the legend title to "Case"
dev.off()


# remove any with greater 2sd from norm (95% interval)
pdf('density.pdf')
densityplot(PC1, pch=19)
densityplot(PC2, pch=19)
densityplot(PC3, pch=19)
densityplot(PC4, pch=19)
densityplot(PC5, pch=19)
dev.off()

# plot sample respect to PC
PC1 <- as.numeric(pca$x[,1])
PC2 <- as.numeric(pca$x[,2])
PC3 <- as.numeric(pca$x[,3])
PC4 <- as.numeric(pca$x[,4])
PC5 <- as.numeric(pca$x[,5])

sample_outliers=c()
alloutliers=c()
for(i in 1:5){
  a<-subset(rownames(pca$x), pca$x[,i] > (mean(pca$x[,i])+3*sd(pca$x[,i])))
  b<-subset(rownames(pca$x), pca$x[,i] < (mean(pca$x[,i])-3*sd(pca$x[,i])))
  out<-c(a,b)
  sample_outliers <- c(sample_outliers,out)
  print(paste("outliers in PCA",i,":",sep=""))
  print(sample_outliers)
  alloutliers=c(alloutliers,sample_outliers)
  sample_outliers=c()
}

# drop these then re-run
dds_filtered <- dds[, !colnames(dds) %in% alloutliers]

# run with covariates 
design(dds_filtered) <- formula(~ cell_type + PC1 + PC2 + PC3 + PC4 + PC5)

ddsRes <- DESeq(dds_filtered)
res <- results(ddsRes, contrast = c("cell_type", "MDD", "control"))
res_sig <- res[res$pvalue == min(res$pvalue),] # Adjusted Pvalues were not significant
res_df <- as.data.frame(res)


# make volcano plot
res_df$col <- ifelse(res$pvalue < 0.05, "red", "black")
res2 <- subset(res, res$col =="red")

volcano <- function(res,filename, sig=NULL){
    pdf(paste("volcano_", filename, ".pdf",sep=""),width = 5,height = 5)
    plot(-log10(res$pvalue)~res$log2FoldChange,
        xlab="Log",
        ylab=expression(-log[10]*P),
        main="Volcano Plot", ylim=c(0,10), xlim=c(-0.5, 0.5), pch=20, col=res$col)
        par(new=T)
        plot(-log10(res2$pvalue)~res2$log2FoldChange, axes=F,
        xlab="",
        ylab='', 
     ylim=c(0,10), xlim=c(-0.5, 0.5), pch=20, col=res2$col)
        if (is.null(sig)==FALSE){
        abline(h=sig,lty="dotted",lwd=2,col="chartreuse4")
    }
    dev.off()   
}

volcano(res_df,"mdd")

#make venn
# Define sets
miRNA_passed_QC <- c(
  "hsa-miR-9-5p", "hsa-miR-9-3p", "hsa-miR-181b-5p", "hsa-miR-1307-3p",
  "hsa-miR-708-3p", "hsa-miR-125b-5p", "hsa-let-7a-5p", "hsa-miR-26a-5p",
  "hsa-miR-323a-3p", "hsa-miR-543", "hsa-miR-7-5p", "hsa-miR-132-3p",
  "hsa-miR-1180-3p", "hsa-miR-99b-5p", "hsa-miR-99a-5p", "hsa-miR-30a-5p",
  "hsa-miR-29a-3p", "hsa-miR-126-3p", "hsa-miR-181a-5p", "hsa-miR-1307-5p",
  "hsa-miR-100-5p", "hsa-miR-16-5p", "hsa-miR-342-3p", "hsa-miR-127-3p",
  "hsa-miR-411-5p", "hsa-miR-128-3p", "hsa-miR-124-3p", "hsa-miR-425-5p",
  "hsa-miR-191-5p", "hsa-miR-30d-5p", "hsa-let-7f-5p", "hsa-miR-222-3p"
)

author_significant_miRNAs <- c(
  "miR-132-5p", "let-7e-5p", "miR-126-5p", "miR-320b", "miR-421",
  "miR-92a-3p", "miR-885-5p", "miR-708-5p"
)



pdf('venn.pdf', width=7.2, height = 7.2)
# Draw the Venn diagram
venn.plot <- venn.diagram(
  x = list(
    "miRNAs that passed QC" = miRNA_passed_QC,
    "Paper's miRNAs" = author_significant_miRNAs
  ),
  filename = NULL,
  fill = c("skyblue", "tomato"),
  alpha = 0.6,
  cex = 1.5,
  cat.cex = 1.2,
  cat.pos = c(-20, 20),
  cat.dist = 0.05,
  main = "Overlap Between QC-Passed and Author-Significant miRNAs"
)
grid::grid.draw(venn.plot)

dev.off()

