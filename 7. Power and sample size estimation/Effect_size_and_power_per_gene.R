###### Environmental Setup Steps #####

library(DESeq2)
library(gplots)
library(amap)
#install.packages('gplots')
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
#if (!requireNamespace('BiocManager', quietly = TRUE))
#install.packages('BiocManager')

#BiocManager::install('EnhancedVolcano')
#devtools::install_github('kevinblighe/EnhancedVolcano')
library(EnhancedVolcano)
library(ggrepel)
library(data.table)
library(writexl)
library(readr)
library(dplyr)
library(ggplot2)
library(ggplot.multistats)

#   On/Off Switches   #
#       - Analysis 1
#       - Analysis 2


##Set your working directory:
setwd("C:/Users/singh.1815/OneDrive - The Ohio State University/Gurkirat Singh/MS Degree outline/Aim1/DESeq/20240405 Power Statistics")


##Feed in your data

#GREM4
#data <- read.table("countmatrix_GREM.txt", header=T, row.names=1, com='', quote='', check.names=F, sep="\t")

#CabSav
#data <- read.table("countmatrix_CabSav.txt", header=T, row.names=1, com='', quote='', check.names=F, sep="\t")


#Inter-species Comparisons
data <- read.table("countmatrix_orthology_GREM4_CabSav_sorted.txt", header=T, row.names=1, com='', quote='', check.names=F, sep="\t")



## Print to screen the number of lines in the dataset and the number of columns
dim(data)
#37443 24 for GREM
#36688 24 for CabSav



##Set Your Parameters of the Run
lfcvalue <- 1  #this creates a variable which will set the log2foldchange value where inserting a '0'=lfc of 1, '1'=lfc of 2 (what we want)
padjvalue <- 0.05
pvaluevalue <- 0.05



###### Analysis 2: Pair-wise Analysis with DEseq2 ######

#### Preparatory Work ####

if(TRUE) {
  
  ## First you have to open up the count matrix and condition matrix files and make them yourself by copying and pasting the columns you need
  ## You cannot simply write a code to do this because we have different numbers of reps for each tissue type so simply pulling out, for example, every 6 columns won't work
  
  # Insert the Names of the Pairs of Data You Wish to Analyze Here \/ make a nice list for a for loop
  #allname <- list("leaf_root", "preveraisonberry_veraisonberry", "veraisonberry_matureberry")
  #     ^ example from the transcriptome by tissue analysis
  #                 \/ must be written exactly as done in the '...Cond.Input.File.txt' files
  
  #Do just one
  #allname <- list("GREM_Chill_v_GREM_Control")
  #allname <- list("GREM_Freeze_v_GREM_Control")
  #allname <- list("CabSav_Chill_v_CabSav_Control")
  #allname <- list("CabSav_Freeze_v_CabSav_Control")
  
  
  
  #Do a subset
  #allname <- list("GREM_Chill_v_GREM_Control", "GREM_Freeze_v_GREM_Control")
  #allname <- list("CabSav_Chill_v_CabSav_Control", "CabSav_Freeze_v_CabSav_Control")
  
  #INter-species comparison
  allname <- list("GREM_Control_v_CabSav_Control")
  
  
  for (p in allname){
    print (p)
    
    ## Read in the Count Matrix Input File (the data which you will do the comparisons on):
    data <- read.table(paste (p,".countmatrix.txt",sep=""), header=T, row.names=1, com='', quote='', check.names=F, sep="\t")
    
    
    ## Round the data using the round function to get it into integers and no decimal places (since the data may or many not be in integer format if you use CoCo for example)
    data <- round(data)
    
    
    ## Read in the Condition Input File (the comparison you want to make):
    sample <- read.table(paste (p,".Cond.Input.File.txt",sep=""), header=T, row.names=1, com='',
                         quote='', check.names=F, sep="\t", colClasses="factor")
    sample <- sample[match(colnames(data), rownames(sample)),, drop=F]
    sample_rowname <- rownames(sample)
    
    
    #### Variability Check with PCA Plot ####
    
    
    ## Read Files In to DESeq2-usable Variables with Specific Format Parameters Dictated
    # (DESeqDataSetFromMatrix is a DESeq2 parameter we are currently setting)
    ddsFullCountTable <- DESeqDataSetFromMatrix(countData = data, 
                                                colData = sample,  
                                                design = ~ condition)
    head(assay(ddsFullCountTable), 3)
    
    
    
    
    
    ## Normalize the 'ddsFullCountTable' input count data so that is can be more visually amiable to plotting
    #       Using the 'rlog' transformation
    rlogddsFullCountTable <- rlog(ddsFullCountTable, blind = FALSE)
    head(assay(rlogddsFullCountTable), 3)
    
    
    ## PCA plot with rlog normalized data - the program will ONLY allow this to be run on normalized data
    
    #Create the PCA Plot (dataset that is)
    #PCAplotoutputdata <- plotPCA(rlogddsFullCountTable,
    #                             ntop = 30000,  #how many of the top most expressed genes should be included in the analysis? (default appears to be 500 of them)
    #                             intgroup = "condition",  #how should you sort the samples for the PCA plot? (by condition in my case)
    #                             returnData = TRUE  #return the output data that would originally be graphed as data and not a graphical output (FALSE produces a graphic)
    #)
    
    #Plot the PCA Data with Sample Labels
    #pdf(paste(p,".rlogPCAPlotwlabels.pdf",sep=""),wi = 8, he = 6)
    #print(ggplot(PCAplotoutputdata,
    #             aes(x = PC1, y = PC2, color = condition)) +
    #       geom_point(size = 5) +
    #        geom_label_repel(aes(label = name)) +
    #        ggtitle(paste(p," PCA Plot of Sample Variance After rlog Transformation",sep="")) +
    #        theme(text = element_text(size = 10)
    #        ))
    #dev.off()
    
    #Plot the PCA Data without Sample Labels
    #pdf(paste(p,".rlogPCAPlot.pdf",sep=""),wi = 8, he = 6)
    #print(ggplot(PCAplotoutputdata,
    #             aes(x = PC1, y = PC2, color = condition)) +
    #        geom_point(size = 5) +
            #geom_label_repel(aes(label = name)) +
    #        ggtitle(paste(p," PCA Plot of Sample Variance After rlog Transformation",sep="")) +
    #        theme(text = element_text(size = 10)
    #        ))
    #dev.off()
    
    
    #### Run DESeq2 ####
    
    
    ## Run DESeq2
    dds <- DESeq(ddsFullCountTable)
    colnames( colData( dds ) )
    
    
    ## Create the naming conventions based on the file name needed for use throughout the rest of this script
    #                               This \/ is where you dictate how to read in your input data (so that is to say, what is your divider in your file name)
    tmpname <- matrix(unlist(strsplit(p,"_v_",2)),ncol = 2)
    sampleA = tmpname[1,1]
    sampleB = tmpname[1,2]  #these two extract sample names and give them to the next function which creates the variable contrastv
    
    # Specify the contrast for differential expression analysis
    contrastV <- c("condition", sampleA, sampleB)
    
    
    
    #Results
    results_dds <- results(dds, contrast = contrastV, lfcThreshold = lfcvalue, alpha = padjvalue, saveCols = "dispersion")
    
    
    # Calculate the Average Expression for Each Gene in Sample A
    baseA <- counts(dds, normalized=TRUE)[, colData(dds)$condition == sampleA]
    if (is.vector(baseA)){
      baseMeanA <- as.data.frame(baseA)
    } else {
      baseMeanA <- as.data.frame(rowMeans(baseA))
    }
    colnames(baseMeanA) <- sampleA
    
    # Calculate the Average Expression for Each Gene in Sample B
    baseB <- counts(dds, normalized=TRUE)[, colData(dds)$condition == sampleB]
    if (is.vector(baseB)){
      baseMeanB <- as.data.frame(baseB)
    } else {
      baseMeanB <- as.data.frame(rowMeans(baseB))
    }
    colnames(baseMeanB) <- sampleB
    
    
    
    # Calculate the standard deviation of expression for each gene in Sample A
    baseSdA <- apply(baseA, 1, sd)
    
    # Calculate the standard deviation of expression for each gene in Sample B
    baseSdB <- apply(baseB, 1, sd)
    
    # Add the mean and standard deviation of expression for Sample A to results_dds
    results_dds$Mean_Expression_SampleA <- baseMeanA[rownames(results_dds), sampleA]
    results_dds$SD_Expression_SampleA <- baseSdA[rownames(results_dds)]
    
    # Add the mean and standard deviation of expression for Sample B to results_dds
    results_dds$Mean_Expression_SampleB <- baseMeanB[rownames(results_dds), sampleB]
    results_dds$SD_Expression_SampleB <- baseSdB[rownames(results_dds)]
    
  
    
    # Add the Gene IDs in Again but with a Column Header of ID this time
    results_dds <- cbind(ID=rownames(results_dds), as.data.frame(results_dds))
    
    
    ## change padj as NA to 1
    results_dds$padj[is.na(results_dds$padj)] <- 1
    
    
    ## Order the listing by p-value and padj value
    results_dds <- results_dds[order(results_dds$pvalue),]
    results_dds <- results_dds[order(results_dds$padj),]
    
    # Calculate pooled standard deviation
    results_dds <- results_dds %>%
      mutate(SD_pooled = sqrt((SD_Expression_SampleA^2 + SD_Expression_SampleB^2) / 2))
    
    # Calculate Cohen's d or efect size
    results_dds <- results_dds %>%
      mutate(cohens_d = (Mean_Expression_SampleA - Mean_Expression_SampleB) / SD_pooled)
    results_dds$cohens_d_absolute <- abs(results_dds$cohens_d)
    
    # Calculate CV for SampleA
    results_dds <- results_dds %>%
      mutate(CV_SampleA = (SD_Expression_SampleA / Mean_Expression_SampleA))
    
    # Calculate CV for SampleB
    results_dds <- results_dds %>%
      mutate(CV_SampleB = (SD_Expression_SampleB / Mean_Expression_SampleB))
    
    effect <- results_dds$cohens_d_absolute
    cv1 <- results_dds$CV_SampleA
    cv2 <- results_dds$CV_SampleB
    n = 4
    
    # Define a function to calculate power for each row
    power_calculation <- function(cv1, cv2, effect, alpha = 0.05, n = 4) {
      pwr.t.test(n = n, d = effect, sig.level = 0.05, alternative = "two.sided")$power
    }
    
    # Add power calculation as a new column
    results_dds <- results_dds %>%
      mutate(power = power_calculation(CV_SampleA, CV_SampleB, cohens_d_absolute))
    
    
    
    summary(results_dds)
    View(results_dds)
    
   
    
  
    # Create an output excel file for DESeq2 Output
    write.xlsx(as.data.frame(results_dds), file = paste0(p, "_results.xlsx"))
    
   
    
  } #Closes the name loop
  
} #closes the on/off switch


##### END #####




##### Interaction Analysis #####


#Inter-species Comparisons
data <- read.table("countmatrix_orthology_GREM4_CabSav_sorted.txt", header=T, row.names=1, com='', quote='', check.names=F, sep="\t")

## Print to screen the number of lines in the dataset and the number of columns
dim(data)

##Set Your Parameters of the Run
lfcvalue <- 1  #this creates a variable which will set the log2foldchange value where inserting a '0'=lfc of 1, '1'=lfc of 2 (what we want)
padjvalue <- 0.05
pvaluevalue <- 0.05

#Do multi-level comparisons
allname <- list("GREM_Chill_v_GREM_Control_x_CabSav_Chill_v_CabSav_Control")
#allname <- list("GREM_Freeze_v_GREM_Control_x_CabSav_Freeze_v_CabSav_Control")


for (p in allname){
  print (p)
  
  ## Read in the Count Matrix Input File (the data which you will do the comparisons on):
  data <- read.table(paste (p,".countmatrix.txt",sep=""), header=T, row.names=1, com='', quote='', check.names=F, sep="\t")
  
  
  ## Round the data using the round function to get it into integers and no decimal places (since the data may or many not be in integer format if you use CoCo for example)
  data <- round(data)
  
  
  ## Read in the Condition Input File (the comparison you want to make):
  sample <- read.table(paste (p,".condition.txt",sep=""), header=T, row.names=1, com='',
                       quote='', check.names=F, sep="\t", colClasses="factor")
  sample <- sample[match(colnames(data), rownames(sample)),, drop=F]
  sample_rowname <- rownames(sample)
  

  ## Read Files In to DESeq2-usable Variables with Specific Format Parameters Dictated
  # (DESeqDataSetFromMatrix is a DESeq2 parameter we are currently setting)
  ddsFullCountTable <- DESeqDataSetFromMatrix(countData = data, 
                                              colData = sample,  
                                              design = ~ species + treatment + species:treatment)
  
  
  ## Ensure that the data was read in correctly by visual inspection
  head(assay(ddsFullCountTable), 3)
  
  ## Ensure that the conditions have been read in my the program correctly
  ddsFullCountTable$species 
  #       The above should output the reference level first, so in this case that should be CabSav
  #       because I am asking the question "Does GREM4 differ from CabSav?", so CabSav should be the reference
  
  ddsFullCountTable$treatment
  #       The above should output the reference level first, so in this case that should be control

  # If there is a reference not set correctly, you can set it using the below
  ddsFullCountTable$species <- relevel(ddsFullCountTable$species, ref="CabSav")
  ddsFullCountTable$treatment <- relevel(ddsFullCountTable$treatment, ref="Control")

  
  #### Run DESeq2 ####
  
  ## Run DESeq2
  dds <- DESeq(ddsFullCountTable)
  
  
  ## Create an Output File of the DESeq2 Results
  
  # Normalize the DESeq2 output data (dds) and convert to a file-printable format
  normalized_counts <- counts(dds, normalized=TRUE)

  #### Prepare DESeq2 Output for Downstream Analysis ####
  
  
  ## Create the naming conventions based on the file name needed for use throughout the rest of this script
  #                               This \/ is where you dictate how to read in your input data (so that is to say, what is your divider in your file name)
  tmpname <- matrix(unlist(strsplit(p,"_x_",2)),ncol = 2)
  sampleA = tmpname[1,1]
  sampleB = tmpname[1,2]  #these two extract sample names and give them to the next function which creates the variable contrastmulti
  
  
  ## Determine what comparisons were conducted in the analysis; determine which comparison you would like to carry forward and use to create the object of 'res'
  #       (That is to say, determine what data (comparisons' data) you wish to use moving forward)
  resultsNames(dds)
  #       This will print to screen the various comparisons conducted - pick one to use moving forward.
  #       In all likelihood, you will want to pick the one using the ':', and there should only be one of these
  #       as this is the multi-level comparison comparison data (interaction comparison data)
  #"Intercept", "species_GREM_vs_CabSav", "treatment_Chill_vs_Control", "speciesGREM.treatmentChill"
  
  ## Create the 'res' object; Dictate what comparison data you wish to use moving forward
  res <- results(dds, name="speciesGREM.treatmentChill", lfcThreshold = lfcvalue, alpha = padjvalue, saveCols = "dispersion")
  #res <- results(dds, name="speciesGREM.treatmentFreeze", lfcThreshold = lfcvalue, alpha = padjvalue, saveCols = "dispersion")
  
  ## Check the 'res' object to ensure that it indeed contains the comparison data you dictated it to contain
  head(res, 10) #display the first 10 lines of the object, check the very first one to ensure it matches the 'name' given above
  
  
  ## Print DEG Summary Statistics to Screen for the Selected Comparison (contained w/i the 'res' object)
  summary(res)
  

  # Calculate the Average Expression for Each Gene in Sample A
  baseA <- counts(dds, normalized=TRUE)[, dds@colData$species == "CabSav"] 
  if (is.vector(baseA)){
    baseMeanA <- as.data.frame(baseA)
  } else {
    baseMeanA <- as.data.frame(rowMeans(baseA))
  }
  colnames(baseMeanA) <- sampleA
  
  # Calculate the Average Expression for Each Gene in Sample B
  baseB <- counts(dds, normalized=TRUE)[, dds@colData$species == "GREM"]
  if (is.vector(baseB)){
    baseMeanB <- as.data.frame(baseB)
  } else {
    baseMeanB <- as.data.frame(rowMeans(baseB))
  }
  colnames(baseMeanB) <- sampleB
  
  # Add the Average Expression For Each Gene for Samples A and B Calculated Just Above to 'res'
  res <- cbind(baseMeanA, baseMeanB, as.data.frame(res))
  
  # Add the Gene IDs in Again but with a Column Header of ID this time
  res <- cbind(ID=rownames(res), as.data.frame(res))
  
  
  
  # Calculate the standard deviation of expression for each gene in Sample A
  baseSdA <- apply(baseA, 1, sd)
  
  # Calculate the standard deviation of expression for each gene in Sample B
  baseSdB <- apply(baseB, 1, sd)
  
  # Add the standard deviation of expression for Sample A to results_dds
  res$SD_Expression_SampleA <- baseSdA[rownames(res)]
  
  # Add the standard deviation of expression for Sample B to results_dds
  res$SD_Expression_SampleB <- baseSdB[rownames(res)]


  
  
  ## change padj as NA to 1
  res$padj[is.na(res$padj)] <- 1
  
  
  ## Order the listing by p-value and padj value
  res <- res[order(res$pvalue),]
  res <- res[order(res$padj),]
  
  # Calculate pooled standard deviation
  res <- res %>%
    mutate(SD_pooled = sqrt((SD_Expression_SampleA^2 + SD_Expression_SampleB^2) / 2))
  
  # Calculate Cohen's d or effect size
  res <- res %>%
    mutate(cohens_d = (GREM_Chill_v_GREM_Control - CabSav_Chill_v_CabSav_Control) / SD_pooled)
  res$cohens_d_absolute <- abs(res$cohens_d)
  
  # Calculate CV for SampleA
  res <- res %>%
    mutate(CV_SampleA = (SD_Expression_SampleA / GREM_Chill_v_GREM_Control))
  
  # Calculate CV for SampleB
  res <- res %>%
    mutate(CV_SampleB = (SD_Expression_SampleB / CabSav_Chill_v_CabSav_Control))
  
  effect <- res$cohens_d_absolute
  cv1 <- res$CV_SampleA
  cv2 <- res$CV_SampleB
  n = 4
  
  # Define a function to calculate power for each row
  power_calculation <- function(cv1, cv2, effect, alpha = 0.05, n = 4) {
    pwr.anova.test(n = n, d = effect, sig.level = 0.05, alternative = "two.sided")$power
  }
  
  # Add power calculation as a new column
  res <- res %>%
    mutate(power = power_calculation(CV_SampleA, CV_SampleB, cohens_d_absolute))
  
  
  
  summary(res)
  View(res)
  
  
  
  
  # Create an output excel file for DESeq2 Output
  write.xlsx(as.data.frame(res), file = paste0(p, "_results.xlsx"))





