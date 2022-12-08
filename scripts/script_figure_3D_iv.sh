# Last checked: 7th December 2022
# R version 4.0.5

# Load packages
library(ggplot2) # version 3.3.6

# Read data object
path <- "/Users/seanwen/Documents/Eleni/Eleni_2020/Github/Louka_JExpMed_2021/data/"
file <- "data.rdata"
load(paste(path, file, sep=""))

# Retrieve genotype proportion
    # Retrieve
    df.1 <- object$single_cell_genotyping$cd38neg
    df.2 <- object$single_cell_genotyping$cd38pos
    
    # Merge
    df <- rbind.data.frame(df.1, df.2)
    
    # Subset relevant columns
    df <- df[, c("celltype", "genotype")]
    
# Set factor levels
    # Cell type
        # Recode levels
        table(df$celltype)
        df$celltype[which(df$celltype=="DP")] <- "+/+"
        
        # Set factor levels
        df$celltype <- factor(df$celltype, levels=c("HSC", "MPP", "LMPP",
                                                "+/+", "CMP", "MEP",
                                                "GMP"))
        table(df$celltype)
                                                
    # Genotypes
        # Recode levels
        table(df$genotype)
        df$genotype[which(df$genotype=="Mon7_WT_NRAS_MUT_SETBP1_WT")] <- "NRAS"
        df$genotype[which(df$genotype=="Mon7_WT_NRAS_MUT_SETBP1_MUT")] <- "NRAS_SETBP1"
        df$genotype[which(df$genotype=="Mon7_MUT_NRAS_MUT_SETBP1_WT")] <- "NRAS_mono7"
        df$genotype[which(df$genotype=="Mon7_MUT_NRAS_WT_SETBP1_WT")] <- "Others"
        df$genotype[which(df$genotype=="Mon7_MUT_NRAS_MUT_SETBP1_MUT")] <- "Others"
    
        # Set factor levels
        df$genotype <- factor(df$genotype, levels=c("Others", "NRAS_mono7", "NRAS_SETBP1", "NRAS"))
        table(df$genotype)
    
    
# Cross tabulate
table(df$genotype, df$celltype)

# Convert from n to % by cell type
    # Reformat data frame
    df <- as.data.frame(table(df$celltype, df$genotype))
    names(df) <- c("celltype", "genotype", "Freq")
    
    # Compute %
    cell.types <- levels(df$celltype)
    
    df.list <- list()

    for(i in 1:length(cell.types)) {

        df.small <- df[which(df$celltype==cell.types[i]), ]
        df.small$Pct <- df.small$Freq/sum(df.small$Freq)*100
        df.list[[i]] <- df.small
        
    }
    
    df <- do.call(rbind.data.frame, df.list)

# Stacked barplots
    # Definition
    data <- df
    x <- data$celltype
    y <- data$Pct
    z <- data$genotype
    maintitle <- ""
    ytitle <- "% of Population"
    xtitle <- ""
    legendtitle <- "Genotype"
    ymin <- 0; ymax <- 100; yinterval <- 20
    
    # Plot
    plot <- ggplot() +
            geom_bar(data, mapping=aes(x=x, y=y, fill=z), color="black", position="stack", stat="identity") +
            scale_y_continuous(limits=c(ymin, ymax), breaks=seq(ymin, ymax, by=yinterval)) +
            scale_fill_manual(values=c("gray", "blue", "purple", "darkorange")) +
            labs(title=maintitle, x=xtitle, y=ytitle, fill=legendtitle) +
            theme(panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank(),
                  panel.background=element_blank(),
                  axis.line=element_line(colour="black"),
                  axis.text=element_text(size=13),
                  axis.title=element_text(size=18),
                  axis.text.x=element_text(colour="black"),
                  plot.title=element_text(size=18)
                  )

              
    # Save file
    path <- "/Users/seanwen/Documents/Eleni/Eleni_2020/Github/Louka_JExpMed_2021/figures/"
    file <- "figure_3D_iv.png"
    ggsave(paste(path, file, sep=""), plot, width=8, height=8)

