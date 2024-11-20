
#### Make data frame#######

#install.packages("/home2/BI311F/scripts/plyr_1.8.6.tar.gz", repos = NULL, type="source")
#install.packages("/home2/BI311F/scripts/dplyr_0.8.5.tar.gz", repos = NULL, type="source")
#install.packages("/home2/BI311F/scripts/rlang_0.4.10.tar.gz", repos = NULL, type="source"

install.packages("dplyr")
library(dplyr)
library(plyr)
library(rlang)

#file_list = list.files(pattern="*[0-9].txt$") 
file_list = list.files(pattern="count.*.txt$")
 
###Les inn data fra filnavnene
datalist = lapply(file_list, function(x){
  dat = read.table(file=x, header=T, sep = "\t")
  names(dat)[2] = x
  return(dat)
})
#library(plyr); library(dplyr)

library(plyr, lib.loc = "/home2/BI311F//R/x86_64-pc-linux-gnu-library/3.6/")
library(dplyr, lib.loc = "/home2/BI311F//R/x86_64-pc-linux-gnu-library/3.6/")

joined <- join_all(dfs = datalist,by = "Geneid",type ="full" )  

tail(joined,5)
head(joined,5)
dim(joined)
colnames(joined)
joined2 <- dplyr::select(joined,contains("sam"))
head(joined2)
dim(joined2)
rownames(joined2) <- joined$Geneid 
rownames(joined2) <- gsub(" ","_",rownames(joined2))
write.table(joined2, file="Count_DataFrame.csv",quote=F,sep="\t")

######LIMMA######################

# Load necessary libraries
library(BiocManager)
library(edgeR)
library(limma)
library(reshape2)
library(ggplot2)

# Install edgeR if not already installed (only needs to be done once)
if (!requireNamespace("edgeR", quietly = TRUE)) {
  BiocManager::install("edgeR")
}

# Load and prepare data
df <- read.table(file="Count_DataFrame.csv", header=TRUE, row.names=1)
df[is.na(df)] <- 0
head(df)

# Create DGEList object
df.mat.exp1.dge <- DGEList(df)
plotMDS(df.mat.exp1.dge)

# Filter low counts
keep <- rowSums(df.mat.exp1.dge$counts > 10) >= ncol(df.mat.exp1.dge)
df.mat.exp1.dge <- df.mat.exp1.dge[keep, ]
dim(df.mat.exp1.dge)

# Define group factors
group <- factor(ifelse(grepl("LC_37|LC_38|LC_39", colnames(df.mat.exp1.dge$counts)), 
                       "Lethal", "control"))

# Recalculate library sizes
df.mat.exp1.dge$samples$lib.size <- colSums(df.mat.exp1.dge$counts)

# Normalize counts
method <- "TMM"
df.mat.exp1.dge <- calcNormFactors(df.mat.exp1.dge, method=method)

# Plot MDS
cc <- ifelse(group == "Lethal", "pink", "green3")
pdf(file="MDS.pdf", width=6, height=6)
plotMDS(df.mat.exp1.dge, col=cc)
dev.off()

# Design matrix and model fitting
des <- model.matrix(~0 + group)
colnames(des) <- levels(group)
v <- voom(df.mat.exp1.dge, design=des, plot=FALSE)
fit <- lmFit(v, design=des)

# Contrast and fit
contrasts <- makeContrasts(Cond = Lethal - control, levels=des)
fit2 <- contrasts.fit(fit, contrasts=contrasts)
fit2 <- eBayes(fit2)

# Summarize results
topTable_Lethal_vs_control <- topTable(fit2, coef="Cond", sort.by="P", adjust.method="BH", n=Inf)
write.table(topTable_Lethal_vs_control, file="topTable_Lethal_vs_control.csv", quote=FALSE, col.names=NA)
head(topTable_Lethal_vs_control, 10)




############################### New code for volcano plot###############


# Load data
data <- read.table(file="topTable_Lethal_vs_control.csv", header=TRUE, sep=",", check.names=FALSE)

# Check data structure
head(data)

# Rename the first column to "Gene"
colnames(data)[1] <- "Gene"

# Load required libraries
library(ggplot2)
library(ggrepel)
library(calibrate)

# Mark significant genes
data$Significant <- as.factor(data$adj.P.Val < 0.05)

# Filter for significant genes only
data <- data[data$adj.P.Val < 0.05,]

# Create PDF for the plot
pdf(file="Lethal_vs_control.pdf", width=8, height=6)

# Generate the plot
ggplot(data, aes(logFC, -log10(adj.P.Val), colour=Significant)) + 
  geom_point() +
  theme_bw(base_size=8) +
  labs(x="logFC (Lethal vs Control)", y="Adjusted p-value (-log10)") +
  geom_text_repel(data=subset(data, -log10(adj.P.Val) > 3.5), segment.size=0.2,
                  aes(logFC, label=Gene), size=2, colour="black") +
  scale_colour_manual(values=c("TRUE"="#FF3333", "FALSE"="black")) +
  geom_vline(xintercept=1, colour="grey") + 
  geom_vline(xintercept=-1, colour="grey") 

# Close the PDF
dev.off()


