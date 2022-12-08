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
    
    # Cell type: HSC
    sample.info.HSC <- sample.info[which(sample.info$Population=="HSC"), ]
    df.array.HSC <- df.array[, which(names(df.array) %in% row.names(sample.info.HSC))]
    table(rownames(sample.info.HSC)==colnames(df.array.HSC))
    
# Retrieve relevant genes: Differential analysis
    # JMML HSC vs. GMP
        # DESeq2
        deseq2 <- DESeqDataSetFromMatrix(countData=df.array.JMML,
                                         colData=sample.info.JMML,
                                         design=~Population)
        out <- DESeq(deseq2)
        out <- results(out, contrast=c("Population", "HSC", "GMP"))
        
        # Subset significant genes
        out.sig <- out[which(out$padj < 0.25), ]
        nrow(out.sig)
        
        # Save as new object
        out.sig.1 <- out.sig
        
    # JMML HSC vs. CB HSC
        # DESeq2
        deseq2 <- DESeqDataSetFromMatrix(countData=df.array.HSC,
                                         colData=sample.info.HSC,
                                         design=~Type)
        out <- DESeq(deseq2)
        out <- results(out, contrast=c("Type", "JMML", "CB"))
        
        # Subset significant genes
        out.sig <- out[which(out$padj < 0.25), ]
        nrow(out.sig)
        
        # Save as new object
        out.sig.2 <- out.sig
        
        # Check overlap with previous
        length(intersect(row.names(out.sig.1), row.names(out.sig.2)))
                
# Retrieve relevant genes: Highly expressed in JMML HSC
    # Convert counts to TPM
        # Normalize gene length for each gene
        df.array.scaled <- sweep(x=df.array, MARGIN=1, FUN='/', STATS=genelength)
        
        # Normalize by library size for each sample
        size.factor <- colSums(df.array.scaled)
        df.array.scaled <- sweep(x=df.array.scaled, MARGIN=2, FUN='/', STATS=size.factor)
        
        # Convert values to per million
        df.array.scaled <- df.array.scaled * 1e6
    
    # Subset relevant samples
        # Define samples
        samples <- row.names(sample.info[which(sample.info$Sample=="JMML_HSC"), ])
        
        # Subset
        df.array.scaled.small <- df.array.scaled[, which(names(df.array.scaled) %in% samples)]
        
    # Subset highly expressed genes
        # Tabulate n cells in which each gene is expressed
        n.cells <- apply(df.array.scaled.small, 1, function(x) {sum(x > 1)})
        
        # Subset genes expressed above threshold
        genes.high <- names(n.cells[which(n.cells > 3)])
        
# Retrieve JMML HSC gene signature
    # Find intersection
    overlap <- intersect(row.names(out.sig.1), intersect(row.names(out.sig.2), genes.high))
    
    # Read gene name-id file
    gene.info <- object$bulk_rnaseq$gene.metadata
        
    # Retrieve gene names
    overlap <- gene.info[which(gene.info$gene_id %in% overlap), ]
    
    # Recode selected gene names
    overlap$gene_name[which(overlap$gene_name=="C14orf169")] <- "RIOX1"
    overlap$gene_name[which(overlap$gene_name=="RP11-732A19.1")] <- "AC091564.1"
    overlap$gene_name[which(overlap$gene_name=="RP11-196G11.2")] <- "AC135050.3"
    
# Heatmap
    # Create sample ID column for metadata
    sample.info$Sample.ID <- row.names(sample.info)
    
    # Reorder columns by samples
    table(sample.info$Sample)
    index.CB_HSC <- which(colnames(df.array.scaled) %in% row.names(sample.info[which(sample.info$Sample=="CB_HSC"),]))
    index.CB_GMP <- which(colnames(df.array.scaled) %in% row.names(sample.info[which(sample.info$Sample=="CB_GMP"),]))
    index.JMML_HSC <- which(colnames(df.array.scaled) %in% row.names(sample.info[which(sample.info$Sample=="JMML_HSC"),]))
    index.JMML_DP <- which(colnames(df.array.scaled) %in% row.names(sample.info[which(sample.info$Sample=="JMML_DP"),]))
    index.JMML_GMP <- which(colnames(df.array.scaled) %in% row.names(sample.info[which(sample.info$Sample=="JMML_GMP"),]))
    df.array.scaled <- df.array.scaled[, c(index.CB_HSC, index.CB_GMP, index.JMML_HSC, index.JMML_DP, index.JMML_GMP)]
    
    # Subset signature genes
    df.array.scaled <- df.array.scaled[overlap$gene_id,]
    table(overlap$gene_id==row.names(df.array.scaled))
    row.names(df.array.scaled) <- overlap$gene_name
    
    # Save as new object
    data <- df.array.scaled
    
    # Define colors
        # Columns
            # Create group labels
            annotation.col.df <- data.frame("Sample.ID"=colnames(data), stringsAsFactors=FALSE)
            annotation.col.df <- join(annotation.col.df, sample.info[,c("Sample.ID", "Sample")], by="Sample.ID", type="left")
            row.names(annotation.col.df) <- annotation.col.df$Sample.ID
            annotation.col.df$Sample.ID <- NULL
            names(annotation.col.df) <- "Group"
            annotation.col.df$Group <- as.character(annotation.col.df$Group)
            annotation.col.df$Group <- factor(annotation.col.df$Group, levels=c("CB_HSC", "CB_GMP", "JMML_HSC", "JMML_DP", "JMML_GMP"))
            
            # Create group colors
            annotation.col.color <- c("black", "grey","blue", "orange", "red")
            names(annotation.col.color) <- levels(annotation.col.df$Group)
            
        # Create final list
        annotation.colors <- list("Group"=annotation.col.color)
    
    # Plot
    pal2 <- colorRampPalette(list('#313695','#313695','#313695','#ffffbf','#d73027','#a50026','#a50026'))(100)
    
    pheatmap(log2(df.array.scaled + 1),
             clustering_method="ward.D",
             cluster_rows=TRUE,
             cluster_cols=FALSE,
             scale="row",
             color=pal2,
             show_colnames=FALSE,
             show_rownames=TRUE,
             legend=TRUE,
             fontsize=10,
             annotation_col=annotation.col.df,
             annotation_colors=annotation.colors
             )
             
    # Note genes also highly expressed in >1 non-JMML HSC groups
    genes.high.JMML_GMP <- c("PRAM1", "HMGB2", "UHRF1", "DNMT1")
    genes.high.CB_HSC <- c("MYCT1", "OXT")
    
    # Subset JMML HSC-specific genes
    overlap <- overlap[-which(overlap$gene_name %in% c(genes.high.JMML_GMP, genes.high.CB_HSC)), ]

# Summarize numbers
results <- data.frame("Category"=c("JMML HSC vs. JMML GMP",
                                 "JMML HSC vs. CB HSC",
                                 ">1 TPM in >50% JMML HSC"
                                 ),
                      "Genes"=c(nrow(out.sig.1),
                                length(intersect(row.names(out.sig.1), row.names(out.sig.2))),
                                nrow(overlap)
                                ),
                      stringsAsFactors=FALSE
                      )
                      
# Save file
path <- "/Users/seanwen/Documents/Eleni/Eleni_2020/Github/Louka_JExpMed_2021/figures/"
file <- "figure_5J.txt"
write.table(results, paste(path, file, sep=""), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

# Heatmap
    # Subset signature genes
    df.array.scaled <- df.array.scaled[overlap$gene_name,]
    
    # Save as new object
    data <- df.array.scaled
    
    # Define colors
        # Columns
            # Create group labels
            annotation.col.df <- data.frame("Sample.ID"=colnames(data), stringsAsFactors=FALSE)
            annotation.col.df <- join(annotation.col.df, sample.info[,c("Sample.ID", "Sample")], by="Sample.ID", type="left")
            row.names(annotation.col.df) <- annotation.col.df$Sample.ID
            annotation.col.df$Sample.ID <- NULL
            names(annotation.col.df) <- "Group"
            annotation.col.df$Group <- as.character(annotation.col.df$Group)
            annotation.col.df$Group <- factor(annotation.col.df$Group, levels=c("CB_HSC", "CB_GMP", "JMML_HSC", "JMML_DP", "JMML_GMP"))
            
            # Create group colors
            annotation.col.color <- c("black", "grey","blue", "orange", "red")
            names(annotation.col.color) <- levels(annotation.col.df$Group)
            
        # Create final list
        annotation.colors <- list("Group"=annotation.col.color)
    
    # Plot
    pal2 <- colorRampPalette(list('#313695','#313695','#313695','#ffffbf','#d73027','#a50026','#a50026'))(100)
    
    path <- "/Users/seanwen/Documents/Eleni/Eleni_2020/Github/Louka_JExpMed_2021/figures/"
    file <- "figure_5L.pdf"
    pdf(paste(path, file, sep=""), width=5, height=3)
    
    pheatmap(log2(df.array.scaled + 1),
             clustering_method="complete",
             cluster_rows=TRUE,
             cluster_cols=FALSE,
             scale="row",
             color=pal2,
             show_colnames=FALSE,
             show_rownames=TRUE,
             legend=TRUE,
             fontsize=10,
             annotation_col=annotation.col.df,
             annotation_colors=annotation.colors
             )
             
    dev.off()
             
    

          
          
          
