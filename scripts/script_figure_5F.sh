# Last checked: 7th December 2022
# R version 4.0.5

# Load packages
library(ggplot2) # version 3.3.6

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
    
# Convert counts to TPM
    # Normalize gene length for each gene
    df.array.scaled <- sweep(x=df.array, MARGIN=1, FUN='/', STATS=genelength)
    
    # Normalize by library size for each sample
    size.factor <- colSums(df.array.scaled)
    df.array.scaled <- sweep(x=df.array.scaled, MARGIN=2, FUN='/', STATS=size.factor)
    
    # Convert values to per million
    df.array.scaled <- df.array.scaled * 1e6
    
    # Save as matrix
    df.mat <- as.matrix(df.array.scaled)
    
# Set factor levels
table(sample.info$Sample)
sample.info$Sample <- factor(sample.info$Sample, levels=c("CB_HSC", "CB_GMP", "JMML_HSC", "JMML_DP", "JMML_GMP"))

# Plot
    # Read gene name-id file
    gene.info <- object$bulk_rnaseq$gene.metadata
    
    # Define genes
    gene.names <- c("HMGA2", "CNN3", "VNN2")
    gene.info.small <- gene.info[which(gene.info$gene_name %in% gene.names), ]

    for(i in 1:nrow(gene.info.small)) {
    
        # Retrieve expression values
        exp <- df.mat[gene.info.small$gene_id[i], , drop=FALSE]
        exp <- as.data.frame(t(exp))
        names(exp) <- "TPM"
        
        # Annotate with sample type
        table(row.names(exp)==row.names(sample.info))
        exp$Sample <- sample.info$Sample
        
        # Definition
        data <- exp
        x <- data$Sample
        y <- log2(data$TPM + 1)
        z <- data$Sample
        maintitle <- gene.info.small$gene_name[i]
        ytitle <- "log2(TPM + 1)"
        xtitle <- ""
        
        # Plot
        plot <- ggplot() +
            geom_boxplot(data, mapping=aes(x=x, y=y, fill=z)) +
            stat_summary(data, mapping=aes(x=x, y=y), fun.y=mean, colour="black", geom="point", shape=23, size=5, stroke=2, fill="white") +
            scale_fill_manual(values=c("black", "grey","blue", "orange", "red")) +
            labs(title=maintitle, x=xtitle, y=ytitle) +
            theme_bw() +
            theme(panel.background = element_blank(),
                  panel.border = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.grid.major = element_blank(),
                  axis.ticks.x=element_blank(),
                  axis.line = element_line(colour="black"),
                  axis.text.x=element_text(size=5, colour="black"),
                  axis.text.y=element_text(colour="black"),
                  axis.text=element_text(size=16),
                  axis.title=element_text(size=16),
                  plot.title=element_text(hjust = 0.5, size=16),
                  legend.position="none"
                  )
                  
        # Save file
        path <- "/Users/seanwen/Documents/Eleni/Eleni_2020/Github/Louka_JExpMed_2021/figures/"
        file <- paste("figure_5F_", gene.info.small$gene_name[i], ".pdf", sep="")
        ggsave(paste(path, file, sep=""), plot, width=3, height=5)
        
    }
