#set working directory
#setwd()

# Load necessary libraries
library(GEOquery)
library(limma)
library(EnhancedVolcano) # for enhanced volcano plot

# Load series and platform data from GEO
gset <- getGEO("GSE52509", GSEMatrix = TRUE, AnnotGPL = TRUE)
if (length(gset) > 1) idx <- grep("GPL6885", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# Make proper column names to match topTable
fvarLabels(gset) <- make.names(fvarLabels(gset))

# Group membership for all samples
gsms <- "000111222333"
sml <- strsplit(gsms, split = "")[[1]]

# Log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE))
LogC <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)
if (LogC) {
  ex[which(ex <= 0)] <- NaN
  exprs(gset) <- log2(ex)
}

# Assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("4 month CS", "4 month control", "6 month CS", "6 month control"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

gset <- gset[complete.cases(exprs(gset)), ] # Skip missing values

fit <- lmFit(gset, design)  # Fit linear model

# Set up contrasts of interest and recalculate model coefficients
cts <- c(paste(groups[1], "-", groups[2], sep = ""), paste(groups[3], "-", groups[4], sep = ""))
cont.matrix <- makeContrasts(contrasts = cts, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)

# Compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)

# Extract data for volcano plot
tT2 <- topTable(fit2, adjust = "fdr", sort.by = "B", number = Inf)

#create volcano plots

#4 months
library(EnhancedVolcano)

# Create keyvals for coloring
keyvals_4 <- ifelse(
  tT2$adj.P.Val < 0.05 & tT2$X4.month.CS.X4.month.control > 0.58, "red", 
  ifelse(tT2$adj.P.Val < 0.05 & tT2$X4.month.CS.X4.month.control < -0.58, "blue", "black")
)

# Assign names to keyvals (required by EnhancedVolcano)
names(keyvals_4) <- rep("NS", length(keyvals_4))
names(keyvals_4)[tT2$adj.P.Val < 0.05 & tT2$X4.month.CS.X4.month.control > 0.58] <- "Upregulated"
names(keyvals_4)[tT2$adj.P.Val < 0.05 & tT2$X4.month.CS.X4.month.control < -0.58] <- "Downregulated"


# Calculate upregulated and downregulated genes
upregulated_4 <- tT2[tT2$adj.P.Val < 0.05 & tT2$X4.month.CS.X4.month.control > 0.58, ]
downregulated_4 <- tT2[tT2$adj.P.Val < 0.05 & tT2$X4.month.CS.X4.month.control < -0.58, ]

num_up_4 <- nrow(upregulated_4)   # Number of upregulated genes
num_down_4 <- nrow(downregulated_4) # Number of downregulated genes


# Create Enhanced Volcano Plot
volcano_plot_4 <- EnhancedVolcano(tT2,
                                lab = tT2$Gene.ID, 
                                x = 'X4.month.CS.X4.month.control', # Log fold change column
                                y = 'adj.P.Val', # Adjusted p-value column
                                title = '4 month CS-exposed vs. Control',
                                xlab = expression(bold(log[2] ~ fold ~ change)),
                                ylab = expression(bold('-log'[10] ~ adj.P.Val)),
                                pCutoff = 0.05, # P-value cutoff
                                FCcutoff = 0.58, # Fold change cutoff
                                pointSize = 2.0, # Point size
                                colAlpha = 0.8, # Transparency
                                colCustom = keyvals_4, # Use keyvals for custom colors
                                legendPosition = 'top', # Legend position
                                legendLabSize = 12, # Legend text size
                                drawConnectors = FALSE)
volcano_plot_4

# Add annotation for upregulated and downregulated counts
volcano_plot_4 +
  annotate("text", 
           x = max(tT2$X4.month.CS.X4.month.control),
           y = max(-log10(tT2$adj.P.Val)) + 5, # Offset above max y-axis value
           label = paste("Upregulated: ", num_up_4), 
           color = "red", 
           size = 5, 
           hjust = 1) +
  annotate("text", 
           x = min(tT2$X4.month.CS.X4.month.control) , 
           y = max(-log10(tT2$adj.P.Val)) + 5, # Offset above max y-axis value
           label = paste("Downregulated: ", num_down_4), 
           color = "blue", 
           size = 5, 
           hjust = 0)

#TOsave images without changing up and down colors
ggsave("4monthCS_volcanoplot.jpg", units="in", width=10, height=10, dpi=300)
jpeg("4monthCS_volcanoplot.jpg", units="in", width=10, height=10, res=300)
dev.off()


#########################################################################################################

#6 months

# Create keyvals for coloring
keyvals_6 <- ifelse(
  tT2$adj.P.Val < 0.05 & tT2$X6.month.CS.X6.month.control > 0.58, "red", 
  ifelse(tT2$adj.P.Val < 0.05 & tT2$X6.month.CS.X6.month.control < -0.58, "blue", "black")
)

# Assign names to keyvals (required by EnhancedVolcano)
names(keyvals_6) <- rep("NS", length(keyvals_6))
names(keyvals_6)[tT2$adj.P.Val < 0.05 & tT2$X6.month.CS.X6.month.control > 0.58] <- "Upregulated"
names(keyvals_6)[tT2$adj.P.Val < 0.05 & tT2$X6.month.CS.X6.month.control < -0.58] <- "Downregulated"


# Calculate upregulated and downregulated genes
upregulated_6 <- tT2[tT2$adj.P.Val < 0.05 & tT2$X6.month.CS.X6.month.control > 0.58, ]
downregulated_6 <- tT2[tT2$adj.P.Val < 0.05 & tT2$X6.month.CS.X6.month.control < -0.58, ]

num_up_6 <- nrow(upregulated_6)   # Number of upregulated genes
num_down_6 <- nrow(downregulated_6) # Number of downregulated genes


# Create Enhanced Volcano Plot
volcano_plot_6 <- EnhancedVolcano(tT2,
                                  lab = tT2$Gene.ID, 
                                  x = 'X6.month.CS.X6.month.control', # Log fold change column
                                  y = 'adj.P.Val', # Adjusted p-value column
                                  title = '6 month CS-exposed vs. Control',
                                  xlab = expression(bold(log[2] ~ fold ~ change)),
                                  ylab = expression(bold('-log'[10] ~ adj.P.Val)),
                                  pCutoff = 0.05, # P-value cutoff
                                  FCcutoff = 0.58, # Fold change cutoff
                                  pointSize = 2.0, # Point size
                                  colAlpha = 0.8, # Transparency
                                  colCustom = keyvals_6, # Use keyvals for custom colors
                                  legendPosition = 'top', # Legend position
                                  legendLabSize = 12, # Legend text size
                                  drawConnectors = FALSE)
volcano_plot_6

# Add annotation for upregulated and downregulated counts
volcano_plot_6 +
  annotate("text", 
           x = max(tT2$X6.month.CS.X6.month.control),
           y = max(-log10(tT2$adj.P.Val)) + 5, # Offset above max y-axis value
           label = paste("Upregulated: ", num_up_6), 
           color = "red", 
           size = 5, 
           hjust = 1) +
  annotate("text", 
           x = min(tT2$X6.month.CS.X6.month.control)-1 , 
           y = max(-log10(tT2$adj.P.Val)) + 5, # Offset above max y-axis value
           label = paste("Downregulated: ", num_down_6), 
           color = "blue", 
           size = 5, 
           hjust = 0)

#TOsave images without changing up and down colors
ggsave("6monthCS_volcanoplot.jpg", units="in", width=10, height=10, dpi=300)
jpeg("6monthCS_volcanoplot.jpg", units="in", width=10, height=10, res=300)
dev.off()


# Filter significant genes for 4 months
significant_genes_4 <- tT2[tT2$adj.P.Val < 0.05 & abs(tT2$X4.month.CS.X4.month.control) > 0.58, ]

# Filter significant genes for 6 months
significant_genes_6 <- tT2[tT2$adj.P.Val < 0.05 & abs(tT2$X6.month.CS.X6.month.control) > 0.58, ]

# Write significant genes to CSV files
write.csv(significant_genes_4, "significant_genes_4_months.csv", row.names = TRUE)
write.csv(significant_genes_6, "significant_genes_6_months.csv", row.names = TRUE)

