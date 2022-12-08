# Last checked: 7th December 2022
# R version 4.0.5

# Load packages
library(DESeq2) # version 1.30.1
library(pheatmap) # version 1.0.12

# Read data object
path <- "/Users/seanwen/Documents/Eleni/Eleni_2020/Github/Louka_JExpMed_2021/data/"
file <- "data.rdata"
load(paste(path, file, sep=""))

# Read feature counts
    # Read file
    df <- object$bulk_rnaseq$counts.matrix
    
    # Retrieve gene IDs and length
    geneids <- df[,1]
    genelength <- df[,6]
    
    # Convert data frame to array
    df.array <- array(df[,-c(1:6)])

    # Create row names using gene IDs
    rownames(df.array) <- geneids
    
    # Provide sample names for array
    . <- gsub(pattern="input.", replacement="", x=colnames(df[,-c(1:6)]), fixed=TRUE)
    . <- gsub(pattern=".Aligned.out.bam", replacement="", x=.)
    colnames(df.array) <- .

# Read sample info file
    # Read file
    sample.info <- object$bulk_rnaseq$sample.metadata
    
    # Recode Sample levles to only include letters, numbers, '_' and '.'
    sample.info$Sample[which(sample.info$Sample=="JMML ++")] <- "JMML DP"
    sample.info$Sample <- gsub(" ", "_", sample.info$Sample)
    
    # Check sample alignment between both files
    table(colnames(df.array)==rownames(sample.info))
    
# Subset relevant cell type
sample.info.HSC <- sample.info[which(sample.info$Population=="HSC"), ]
df.array.HSC <- df.array[, which(names(df.array) %in% row.names(sample.info.HSC))]
table(rownames(sample.info.HSC)==colnames(df.array.HSC))
    
# Differential analysis
    # Deseq2
    deseq2 <- DESeqDataSetFromMatrix(countData=df.array.HSC,
                                     colData=sample.info.HSC,
                                     design=~Type)
    out <- DESeq(deseq2)
    out <- results(out, contrast=c("Type", "JMML", "CB"))
    out.sig <- out[which(out$padj < 0.01 ), ]
    nrow(out.sig)
    
# Convert counts to TPM
    # Normalize gene length for each gene
    df.array.scaled <- sweep(x=df.array, MARGIN=1, FUN='/', STATS=genelength)
    
    # Normalize by library size for each sample
    size.factor <- colSums(df.array.scaled)
    df.array.scaled <- sweep(x=df.array.scaled, MARGIN=2, FUN='/', STATS=size.factor)
    
    # Convert values to per million
    df.array.scaled <- df.array.scaled * 1e6
    
# Subset significant genes from expression array
out.sig <- as.data.frame(out.sig)
sig.genes <- row.names(out.sig)
df.array.scaled <- df.array.scaled[sig.genes,]
    
# Subset relevant samples from expression array
df.array.scaled <- df.array.scaled[, which(names(df.array.scaled) %in% names(df.array.HSC))]

# K-means clustering
    # Plot
    pal2 <- colorRampPalette(list('#313695','#313695','#313695','#ffffbf','#d73027','#a50026','#a50026'))(100)
    
    path <- "/Users/seanwen/Documents/Eleni/Eleni_2020/Github/Louka_JExpMed_2021/figures/"
    file <- "figure_5H.pdf"
    pdf(paste(path, file, sep=""), width=2.5, height=5)
    
    set.seed(123)
    
    pheatmap(log2(df.array.scaled + 1),
             kmeans_k=5,
             clustering_method="ward.D",
             scale="row",
             color=pal2,
             show_colnames=TRUE,
             show_rownames=FALSE,
             treeheight_row=0,
             legend=TRUE,
             fontsize=10
             )
             
    dev.off()

    # Retrieve no. of genes in each cluster
    set.seed(123)
    
    plot <- pheatmap(log2(df.array.scaled + 1),
                     kmeans_k=5,
                     clustering_method="ward.D",
                     scale="row",
                     color=pal2,
                     show_colnames=TRUE,
                     show_rownames=FALSE,
                     treeheight_row=0,
                     legend=TRUE,
                     fontsize=10
                     )

    clusters <- plot$kmeans$size
    
    clusters <- data.frame("Cluster"=c(5:1), "Size"=clusters, stringsAsFactors=FALSE)
    
    clusters <- clusters[order(clusters$Cluster), ]
    clusters
    
    #path <- "/Users/seanwen/Documents/Eleni_2020/Figures/Figure 5G/"
    #file <- "Figure 5G.txt"
    #write.table(clusters, paste(path, file, sep=""), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
    
