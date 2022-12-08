# Last checked: 7th December 2022
# R version 4.0.5

# Load packages
library(DESeq2) # version 1.30.1

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
    
    # Disease type: CB
    sample.info.CB <- sample.info[which(sample.info$Type=="CB"), ]
    df.array.CB <- df.array[, which(names(df.array) %in% row.names(sample.info.CB))]
    table(rownames(sample.info.CB)==colnames(df.array.CB))

    # Cell type: HSC
    sample.info.HSC <- sample.info[which(sample.info$Population=="HSC"), ]
    df.array.HSC <- df.array[, which(names(df.array) %in% row.names(sample.info.HSC))]
    table(rownames(sample.info.HSC)==colnames(df.array.HSC))
    
# Differential analysis
    n.sig <- NULL

    # Subgroup: JMML
        # HSC vs. GMP
        deseq2 <- DESeqDataSetFromMatrix(countData=df.array.JMML,
                                         colData=sample.info.JMML,
                                         design=~Population)
        out <- DESeq(deseq2)
        out <- results(out, contrast=c("Population", "HSC", "GMP"))
        out.sig <- out[which(out$padj < 0.01 ), ]
        n.sig[1] <- nrow(out.sig)
        
        # HSC vs. DP
        deseq2 <- DESeqDataSetFromMatrix(countData=df.array.JMML,
                                         colData=sample.info.JMML,
                                         design=~Population)
        out <- DESeq(deseq2)
        out <- results(out, contrast=c("Population", "HSC", "DP"))
        out.sig <- out[which(out$padj < 0.01 ), ]
        nrow(out.sig)
        n.sig[2] <- nrow(out.sig)
        
        # GMP vs. DP
        deseq2 <- DESeqDataSetFromMatrix(countData=df.array.JMML,
                                         colData=sample.info.JMML,
                                         design=~Population)
        out <- DESeq(deseq2)
        out <- results(out, contrast=c("Population", "GMP", "DP"))
        out.sig <- out[which(out$padj < 0.01 ), ]
        nrow(out.sig)
        n.sig[3] <- nrow(out.sig)
        
    # Subgroup: CB
        # HSC vs. GMP
        deseq2 <- DESeqDataSetFromMatrix(countData=df.array.CB,
                                         colData=sample.info.CB,
                                         design=~Population)
        out <- DESeq(deseq2)
        out <- results(out, contrast=c("Population", "HSC", "GMP"))
        out.sig <- out[which(out$padj < 0.01 ), ]
        nrow(out.sig)
        n.sig[4] <- nrow(out.sig)
        
    # Subgroup: HSC
        # JMML HSC vs. CB HSC
        deseq2 <- DESeqDataSetFromMatrix(countData=df.array.HSC,
                                         colData=sample.info.HSC,
                                         design=~Type)
        out <- DESeq(deseq2)
        out <- results(out, contrast=c("Type", "JMML", "CB"))
        out.sig <- out[which(out$padj < 0.01 ), ]
        nrow(out.sig)
        n.sig[5] <- nrow(out.sig)
    
    # Subgroup: None
        # JMML HSC vs. CB GMP
        deseq2 <- DESeqDataSetFromMatrix(countData=df.array,
                                         colData=sample.info,
                                         design=~Sample)
        out <- DESeq(deseq2)
        out <- results(out, contrast=c("Sample", "JMML_HSC", "CB_GMP"))
        out.sig <- out[which(out$padj < 0.01 ), ]
        nrow(out.sig)
        n.sig[6] <- nrow(out.sig)
        
    # Tabulate results
    results <- data.frame("Comparison"=c("JMML HSC vs. JMML GMP",
                                         "JMML HSC vs. JMML DP",
                                         "JMML GMP vs. JMML DP",
                                         "CB HSC vs. CB GMP",
                                         "JMML HSC vs. CB HSC",
                                         "JMML HSC vs. CB GMP"),
                          "significant.genes"=n.sig,
                          stringsAsFactors=FALSE
                          )
                          
    # Save file
    path <- "/Users/seanwen/Documents/Eleni/Eleni_2020/Github/Louka_JExpMed_2021/figures/"
    file <- "figure_5A.txt"
    write.table(results, paste(path, file, sep=""), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
