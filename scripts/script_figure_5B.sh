# Last checked: 7th December 2022
# R version 4.0.5

# Load packages
library(DESeq2) # version 1.30.1
library(plyr) # version 1.8.4
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
    
# Split by disease type and cell population
    # Disease type: JMML
    sample.info.JMML <- sample.info[which(sample.info$Type=="JMML"), ]
    df.array.JMML <- df.array[, which(names(df.array) %in% row.names(sample.info.JMML))]
    table(rownames(sample.info.JMML)==colnames(df.array.JMML))

# Differential analysis
    # HSC vs. GMP
    deseq2 <- DESeqDataSetFromMatrix(countData=df.array.JMML,
                                     colData=sample.info.JMML,
                                     design=~Population)
    out <- DESeq(deseq2)
    out <- results(out, contrast=c("Population", "HSC", "GMP"))
    out.sig <- out[which(out$padj < 0.01 ), ]
    out.sig.1 <- out.sig
    
    # HSC vs. DP
    deseq2 <- DESeqDataSetFromMatrix(countData=df.array.JMML,
                                     colData=sample.info.JMML,
                                     design=~Population)
    out <- DESeq(deseq2)
    out <- results(out, contrast=c("Population", "HSC", "DP"))
    out.sig <- out[which(out$padj < 0.01 ), ]
    nrow(out.sig)
    out.sig.2 <- out.sig
    
    # GMP vs. DP
    deseq2 <- DESeqDataSetFromMatrix(countData=df.array.JMML,
                                     colData=sample.info.JMML,
                                     design=~Population)
    out <- DESeq(deseq2)
    out <- results(out, contrast=c("Population", "GMP", "DP"))
    out.sig <- out[which(out$padj < 0.01 ), ]
    nrow(out.sig)
    out.sig.3 <- out.sig
    
    # Merge all significant genes
    sig.genes <- unique(c(row.names(out.sig.1), row.names(out.sig.2), row.names(out.sig.3)))
    
# Convert counts to TPM
    # Normalize gene length for each gene
    df.array.scaled <- sweep(x=df.array, MARGIN=1, FUN='/', STATS=genelength)
    
    # Normalize by library size for each sample
    size.factor <- colSums(df.array.scaled)
    df.array.scaled <- sweep(x=df.array.scaled, MARGIN=2, FUN='/', STATS=size.factor)
    
    # Convert values to per million
    df.array.scaled <- df.array.scaled * 1e6
    
# Subset significant genes
    # Subset
    df.array.scaled <- df.array.scaled[sig.genes, ]
    
    # Read gene name-id file
    gene.info <- object$bulk_rnaseq$gene.metadata
    
    # Replace gene IDs with gene names
    . <- data.frame("gene_id"=row.names(df.array.scaled), stringsAsFactors=FALSE)
    . <- join(., gene.info, by="gene_id", type="left")
    row.names(df.array.scaled) <- .$gene_name
    
# Collapse expression by group
    # Define groups
    groups <- c("JMML_HSC", "JMML_GMP", "JMML_DP")
    
    # Compute mean for each group for each gene
    results.list <- NULL
    
    for(i in 1:length(groups)) {
    
        # Subset relevant samples
        samples <- row.names(sample.info[which(sample.info$Sample==groups[i]), ])
        df.array.scaled.small <- df.array.scaled[, which(names(df.array.scaled) %in% samples)]
        
        # Compute mean for each gene
        means <- apply(df.array.scaled.small, 1, function(x) {mean(x)})
        
        # Save results as data frame
        results <- data.frame(means)
        names(results) <- groups[i]
        results.list[[i]] <- results
        
    }
    
    df.array.scaled.collapsed <- do.call(cbind.data.frame, results.list)

# Heatmap
    # Plot
    data <- df.array.scaled.collapsed
    
    pal2 <- colorRampPalette(list('#313695','#313695','#313695','#ffffbf','#d73027','#a50026','#a50026'))(100)

    path <- "/Users/seanwen/Documents/Eleni/Eleni_2020/Github/Louka_JExpMed_2021/figures/"
    file <- "figure_5B.pdf"
    pdf(paste(path, file, sep=""), width=2, height=4)
    
    pheatmap(log2(data + 1),
             clustering_method="complete",
             cluster_rows=TRUE,
             cluster_cols=TRUE,
             scale="row",
             color=pal2,
             show_colnames=TRUE,
             show_rownames=FALSE,
             legend=TRUE,
             fontsize=5,
             border_color=NA
             )
                        
    dev.off()
