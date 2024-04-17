#####   Environmental Setup   #####

library("clusterProfiler")
library("dplyr")
library("tidyr")
library("GOplot")
library("enrichplot")
library("fgsea")
library("data.table")
library("GO.db")
library("AnnotationForge")
library("AnnotationHub")
library("writexl")
library("ggnewscale")
library("readxl")
library("dplyr")
library("ggplot2")
library("DOSE")

#   Set Your Working Directory   #
#setwd("C:/Users/singh.1815/OneDrive - The Ohio State University/Gurkirat Singh/MS Degree outline/Aim1/DESeq/Freezing Tolerance Gene Study")
#setwd("C:/Users/singh.1815/OneDrive - The Ohio State University/Gurkirat Singh/MS Degree outline/Aim1/DESeq/Interspecies comparisons/nostress")
#setwd("C:/Users/singh.1815/OneDrive - The Ohio State University/Gurkirat Singh/MS Degree outline/Aim1/DESeq/Interspecies comparisons/chill")
setwd("C:/Users/singh.1815/OneDrive - The Ohio State University/Gurkirat Singh/MS Degree outline/Aim1/DESeq/Interspecies comparisons/freeze")


#   On/Off Switches   #
#       - ORA
#       - GSEA
#       - treeplot in GSEA

#   Important Selections to Turn On or Off
#       - TERM2GENE file selection based on genome running analysis upon [w/i Input File Prep]




##### Setting Up the Loop #####



### ORA ###




#ORA Enrichments in Genome Specific Genes
#allname <- list("GREM_uniq_gene_list")

#On-Deck
#"PN_uniq_genes_list"



#do some specific ones 
#allname <- list("GREM_Chill_v_GREM_Control")
#allname <- list("GREM_Freeze_v_GREM_Control")
#allname <- list("CabSav_Chill_v_CabSav_Control")
#allname <- list("CabSav_Freeze_v_CabSav_Control")

#allname <- list("GREM_Control_v_CabSav_Control")
#allname <- list("GREM_Chill_v_GREM_Control_x_CabSav_Chill_v_CabSav_Control")
allname <- list("GREM_Freeze_v_GREM_Control_x_CabSav_Freeze_v_CabSav_Control")


for (p in allname){
  print (p)
  
  tmpname <- matrix(unlist(strsplit(p,"_x_",2)),ncol = 2)
  sampleA = tmpname[1,1]
  sampleB = tmpname[1,2]  #these two extract sample names and give them to the next function which creates the variable contrastv
  file_name_higher <- paste(sampleA,"_higherThan_",sampleB, sep="")
  file_name_lower <- paste(sampleA,"_lowerThan_",sampleB, sep="")
  
  #}
  # ^ End the loop with the curly closed bracket
  #         (you can see me place this down below the GSEA analysis)
  
  
  
  
  
  
  
  
  ##### Input File Preparation #####
  
  
  
  ###   Prepare and Provide the Program your GO terms containing file   ###
  
  ## GREM4 TERM2GENE file -
  #First, you have to take the 'C:\Users\dixon.778\OneDrive - The Ohio State University\Gschwend Lab Shared Folder\AAA_Datasets\Genome_Library\V_labrusca\Construction 2 - Anno CD\Functional Anno Output Files\annotations.iprscan.txt'
  #     file and open it in excel to format it to work with the following analyses.
  
  #Add a header row containing the headers 'geneID' and 'GO Terms' over the 1st and 3rd column
  
  #Next, delete out the 2nd column containing information on which database the information in the third column originated from
  
  #Then filter the second column with condition 'Show only those that do not start with the text 'InterPro''
  
  #Then copy those results out (not by using the whole column click to select, because that will
  #       take everything, just use the ctrl+shift+down to grab them all) and paste them into a new tab
  
  #Make a new column in the third column and start to type in 2 or 3 of the GO term numbers from the 
  #       second column.  This gets flash fill working and then flash fill in the rest of the thousands of rows
  
  #Now, make a fourth column to add in the 'GO:' in front of the number.  Use flash fill to do this.
  #       Ensure that the header for this column is 'GO Terms'.
  
  #Delete out columns 2 and 3 as they are no longer needed.  Also, delete out the original tab leaving only this tab
  
  #Add a new column between columns 1 and 2; use flash fill to write in the gene names there without
  #       the transcript tag on the ends of them (remove the '...-T2', etc.)
  
  #Delete out column 1 leaving only the gene names column without the transcript tags now and the GO Term column
  
  #Swap columns 1 and 2 as the gene name needs to be in the 2nd column and GO terms in the 1st
  
  #Save as a .csv file formatted file
  
  # The GREM4 TERM2GENE file should now be ready to go now...
  
  
  ## PN TERM2GENE file - 
  #Start with the file of 'PN40024.v4.1.REF.b2g.annot' found at location 'C:\Users\dixon.778\OneDrive - The Ohio State University\Gschwend Lab Shared Folder\AAA_Datasets\Genome_Library\Vv_Pinot_Noir\Grapedia 2-28-23 12x.2-v4.1\PN40024.v4_11_05_21\Blast2GO_results_PN40024_september_2021\'
  
  #Read it into excel
  
  #Delete out the third column as it contains protein naming information we don't need
  
  #Add a header row containing the headers 'geneID' and 'GO Terms' over the 1st and 2nd column
  
  #Then filter the second column with condition 'Show only those that start with the text 'GO'
  
  #Then copy those results (including the gene names next to them) out (not by using the whole column 
  #       click to select, because that will take everything, just use the ctrl+shift+down to grab 
  #       them all) and paste them into a new tab
  
  #Add a new column to the left of columns 1; use flash fill to write in the gene names there without
  #       the transcript demarcator on the ends of them (remove the '...-T2', etc.)
  
  #Delete out column 2 leaving only the gene names (column without the transcript demarcators now) and 
  #       the GO Term column
  
  #Swap columns 1 and 2 as the gene name needs to be in the 2nd column and GO terms in the 1st
  
  #Save as a .csv file formatted file
  
  #The PN TERM2GENE file should now be ready to go now...
  
  
  
  
  
  ## Importing the file into the program is csv with two columns: 'geneID' and corresponding 'GO Terms'
  
  #GREM4 TERM2GENE - Use this for interspecies comaprisons as well bevcause we have GREM gene IDs
  TERM2GENE <- read.csv("annotations.iprscan.GSEAinputformat.csv", header = TRUE)
  
  #CAbSav TERM2GENE - 
  #TERM2GENE <- read.csv("CabSav_GOterms.csv", header = TRUE)
  
  #Notes:
  #Note- if there is more than one GO Term available for an individual gene, make a new row for it
  #Example code for making new rows for each GO Term
  #TERM2GENE <- tidyr::separate_rows(data = TERM2GENE,ID,sep = ",")
  #Example code for removing duplicate rows (if IPR gave you repetitive results)
  #TERM2GENE <- TERM2GENE %>% distinct()
  #export as csv for future reference
  #write.csv(TERM2GENE,"~/R_dir/Pennycress RNASeq/ClusterProfiler\\GO_annotation.csv", row.names = TRUE)
  #IF YOU FEED IN THE WRONG TERM2GENE FILE FOR THE GENOME YOU ARE CURRENTLY WORKING ON
  #     you will get the following error "No gene can be mapped"
  
  
  ###   Present the Program the Reference Guide for GO Number to Functional Name   ###
  
  #Read in the .tsv GO Term Reference file from QuickGO
  #           The 'QuickGOannoAllArabGO2_13_23.tsv' can be opened in excel and just saved as a .csv - this is the easiest way to do this.  Then read in the .csv below...
  TERM2NAME <- read.csv("QuickGOannoAllArabGO2_13_23.csv", header = TRUE)
  
  #View the original file was imported correctly
  head(TERM2NAME)
  
  #Retain only the GO Term ('GO.TERM')(col 3) and name ('GO.NAME')(col 4) column
  #       (although we will also keep around the 'GO.ASPECT' column for use with the GObubble plot for later)
  TERM2NAME <- TERM2NAME[c("GO.TERM", "GO.NAME", "GO.ASPECT")]
  
  #View the new object looks correct
  head(TERM2NAME)
  
  
  
  
  
  
  
  
  ##### Over-Representation Analysis (ORA) with your own GO annotation File #####
  
  if(TRUE) {
    
    
    
    ###   Provide the Program your list of DEGS   ###
    
    ## START WITH UP REGULATED   /\ /\ /\   ##
    
    
    #Provide the DEG gene list, provide up regulated
    #read in the .csv DEG up or down file
    gene <- read_xlsx(paste(paste(file_name_higher, ".xlsx", sep=""), sep=""))  
    
    #carry forward only the gene ID column
    gene <- gene[["ID"]]  #retain only the 'ID' column from the above file
    
    #Make sure gene list is in character vector format!
    #Verify that visually using the below command
    head(gene)  #view to make sure it looks right
    
    
    ###   Run Over-Representation Analysis   ###
    
    #If using your own GO annotation file as input, use the enricher function
    GO_overrep <- clusterProfiler::enricher(
      gene = gene,
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      #universe = GO_annotation,
      minGSSize = 10, #this sets 'M' in the equation of the program.  It spiders out all over to include GO terms that are only associated with the term to to other genes, thusly, that is why this will return enriched GO terms sets that have counts lower than 10, becuase it found a minimum of 10 'associated' genes that had GO terms associated with the GO term in question for the enrichmed term query
      maxGSSize = 500,
      qvalueCutoff = 0.2,
      TERM2GENE = TERM2GENE,
      TERM2NAME = TERM2NAME
    )
    
    #Filter the results so that the number of counts meets a certain threshold
    GO_overrep <- gsfilter(GO_overrep, #this function sets 'k' in the equation of the program, that is to say, it filters out any enrichments that are not over or under your threshold based on count alone, not using the associated method employed above in the min and maxGSsize option
                           by='Count', #tells it to filter by count, not by 'M' (default) in the equation
                           min=5, 
                           max=500)
    
    #print the output to a file
    write.csv(GO_overrep@result, paste(file_name_higher, "_ORA.csv", sep=""))
    
    
    
    ###   Graph Your Results   ###
    
    
    ## ORA Horizontal Bar Plot
    pdf(paste(file_name_higher, "_ORA_Barplot.pdf", sep=""), 
        width = 12,  #currently dictating that the width must be 12 inches but...
        #height = 7,   #the height is at the discretion of the barplot command
    )
    print(barplot(GO_overrep, #need to do 'print' to get it to work inside a for loop, otherwise, just using the barplot command would work on its own
                  #cex.axis = 3,  #supposed to increase font size of name labels but does not for me
                  #cex.names = 20, #supposed to increase font size of name labels but does not for me
                  drop = TRUE, #not sure what this does
                  showCategory = 10,  #dictate how many to show maximum
                  title = paste(file_name_higher, "Over-Represented GO Terms"), font.title = 50,
                  #font.size = 20  #increases font size of x and y axis labels but not title or legend
    ))
    dev.off()
    
    
    ## ORA Dot Plot
    pdf(paste(file_name_higher, "_ORA_Dotplot.pdf", sep=""), 
        width = 9, 
        height = 8
    ) 
    print(dotplot(GO_overrep, #data
                  title = paste(file_name_higher, "Over-Represented GO Terms") #title
    ))
    dev.off()
    
    
    ## ORA Linkage Plot    XXX DOES NOT WORK FOR SOME COMPARISONS - TURNED OFF
    
    #run the linkage-creation command
    #linkage <- pairwise_termsim(GO_overrep)
    
    #create the linkage plot
    #pdf(paste(file_name_higher, "_ORA_Linkage_Plot.pdf", sep=""), 
    #    #width = 12, 
    #    #height = 7
    #)
    #print(emapplot(linkage #data
    #               )) + ggtitle(paste(file_name_higher, "Over-Represented GO Terms"))
    #dev.off()
    
    #Repeat for your list of down-regulated DEGs
    
    
    
    
    
    ## NOW DO DOWN REGULATED   \/ \/ \/   ##
    
    #Provide the DEG gene list, provide down regulated
    #read in the .csv DEG up or down file
    gene <- read_xlsx(paste(paste(file_name_lower, ".xlsx", sep=""), sep=""))  
    
    #carry forward only the gene ID column
    gene <- gene[["ID"]]  #retain only the 'ID' column from the above file
    
    #Make sure gene list is in character vector format!
    #Verify that visually using the below command
    head(gene)  #view to make sure it looks right
    
    
    
    ###   Run Over-Representation Analysis   ###
    
    #If using your own GO annotation file as input, use the enricher function
    GO_overrep <- enricher(
      gene = gene,
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      #universe = GO_annotation,
      minGSSize = 10, #this sets 'M' in the equation of the program.  It spiders out all over to include GO terms that are only associated with the term to to other genes, thusly, that is why this will return enriched GO terms sets that have counts lower than 10, becuase it found a minimum of 10 'associated' genes that had GO terms associated with the GO term in question for the enrichmed term query
      maxGSSize = 500,
      qvalueCutoff = 0.2,
      TERM2GENE = TERM2GENE,
      TERM2NAME = TERM2NAME
    )
    
    #Filter the results so that the number of counts meets a certain threshold
    GO_overrep <- gsfilter(GO_overrep, #this function sets 'k' in the equation of the program, that is to say, it filters out any enrichments that are not over or under your threshold based on count alone, not using the associated method employed above in the min and maxGSsize option
                           by='Count',
                           min=5, 
                           max=500)
    
    #print the output to a file
    write.csv(GO_overrep@result, paste(file_name_lower, "_ORA.csv", sep=""))
    
    
    
    ###   Graph Your Results   ###
    
    
    ## ORA Horizontal Bar Plot
    pdf(paste(file_name_lower, "_ORA_Barplot.pdf", sep=""), 
        width = 12,  #currently dictating that the width must be 12 inches but...
        #height = 7,   #the height is at the discretion of the barplot command
    )
    print(barplot(GO_overrep, #need to do 'print' to get it to work inside a for loop, otherwise, just using the barplot command would work on its own
                  #cex.axis = 3,  #supposed to increase font size of name labels but does not for me
                  #cex.names = 20, #supposed to increase font size of name labels but does not for me
                  drop = TRUE, #not sure what this does
                  showCategory = 10,  #dictate how many to show maximum
                  title = paste(file_name_lower, "Over-Represented GO Terms"), font.title = 50,
                  #font.size = 20  #increases font size of x and y axis labels but not title or legend
    ))
    dev.off()
    
    
    ## ORA Dot Plot
    pdf(paste(file_name_lower, "_ORA_Dotplot.pdf", sep=""), 
        width = 9, 
        height = 8
    )
    print(dotplot(GO_overrep, #data
                  title = paste(file_name_lower, "Over-Represented GO Terms") #title
    ))
    dev.off()
    
    
    ## ORA Linkage Plot    XXX DOES NOT WORK FOR SOME COMPARISONS - TURNED OFF
    
    #run the linkage-creation command
    #linkage <- pairwise_termsim(GO_overrep)
    
    #create the linkage plot
    #pdf(paste(file_name_lower, "_ORA_Linkage_Plot.pdf", sep=""), 
    #    #width = 12, 
    #    #height = 7
    #)
    #print(emapplot(linkage,
    #               showCategory = 15
    #               ) #data
    #      + ggtitle(paste(file_name_lower, "Over-Represented GO Terms"))
    #)
    #dev.off()
    
  }
  
  
  
  
  
  
  #####  Conducting a GSEA with Your Own GO Anno. File  #####
  
  if(FALSE) {
    
    
    
    ###   Prepare the DEG List Again, but in a slightly different format   ###
    
    #   Re-set Your Working Directory Back to the GSEA Folder (if you every left in the first place)   #
    setwd("C:/Users/dixon.778/OneDrive - The Ohio State University/School/College - Graduate School/19-20/Gschwend Lab/AAA   Science/Coding/Aim 3 IH Timecourse V1/GSEA")
    
    
    ##Feed the DEG list file into R > Prepare geneList variable in decreasing sorted data 
    ##        frame format
    
    #First, feed in the list of all genes and their expression information from DESeq2
    #                             ... the "... DESeq2.results.csv" file
    #                       \/ this should be the data for all genes from from DESeq2 run
    geneList <- read.csv(paste(p, ".coco.count.matrix.DESeq2.results.csv", sep=""), header = TRUE)
    
    
    #Pull out only the two columns we need
    geneList <- geneList %>%
      pull(log2FoldChange, ID)
    
    #Sort the geneList (necessary)(descending order - high on top, low on bottom)
    #     It would appear this function also incidentally removes all genes which possess a log2FoldChange of 'NA' as well, which is convenient
    geneList <- sort(geneList, decreasing = TRUE)
    
    #Remove any rows with an N/A listed (that was what Rachel's file looked like, in my files
    #       I will actually have nothing filling that cell instead, but fortunately, this command
    #       still works for the purpose of removing those rows
    #       (That is, to ensure that all of the 'NA's were removed by the previous sort as well)
    geneList <- na.omit(geneList)
    
    #View the geneList to ensure accuracy
    #       High to low, 'snap' genes are still there (NA didn't accidentally remove them),
    #       only two columns exist of log2FoldChange and ID in that order
    head(geneList)
    
    
    
    ###   Run GSEA Using GSEA Function   ###
    
    ##Run the GSEA analysis itself
    #check manual for more options to specify below
    
    #set the seed to ensure consistent output between runs 
    #       The program requests a random starting seed within the run and this random seed can affect the marginal calls to such an extent that they may or may not be present between runs (their statistics reported are also affected but those enriched terms very high on the list are not affected much at all)
    set.seed(1234)
    
    #Run the GSEA Analysis
    GSEAoutput <- GSEA(geneList = geneList,
                       minGSSize = 10,         #dictate the minimum number of genes necessary within a bin to then run analysis on to determine if it is enriched or not
                       maxGSSize = 500,        #dictate the maximum number of genes necessary within a bin to then run analysis on to determine if it is enriched or not
                       eps = 0,                #dictate the eps, a measure of how stringent you want the pvalue calculations to be run in the fgsea function within the GSEA function 
                       nPerm = 1000000,        #dictate how many permutations you want the program to run to identify the enriched terms
                       pvalueCutoff = 0.05,    #dictate the padj-value cutoff (I know it says 'pvalue' but in looking at the output files, this actually dictates the padj value cutoff)
                       TERM2GENE = TERM2GENE, 
                       TERM2NAME = TERM2NAME,
                       seed = TRUE,            #dictate if you want to dictate in a specific seed that you are providing (TRUE) or if you would like to have it pick one at random (FALSE).  You want to dictate in a consistent seed value in order to attain reproducible results, otherwise, a randomly selected seed value will produce different enrichment results each time (those very clearly enriched will not change (aside from some of their statistics marginally), however, those that are marginally enriched will be heavily affected by this aspect)
                       by = "fgsea"            #dictate the program of GSEA you wish to run (fgsea or DOSE) - see elec lab notebook for more information on which of these you should select - you should go with fgsea though FYI
    )
    
    #Export results to data frame and then to excel file
    GSEA_output_writable <- as.data.frame(GSEAoutput)
    write_xlsx(GSEA_output_writable, paste(p, "_GSEA.xlsx", sep=""))
    
    
    
    ###   Graph GSEA Results   ###
    
    
    ## GSEA Horizontal Bar Plot
    ## I don't think it is possible to run this type of plot for this analysis but that is fine
    ## My attempts at it are below... 
    
    #We have to reformat the data a bit to allow it to be used in this graphic
    
    #GSEAoutputformattedforBarplot <- GSEAoutput$setSize
    
    #names(GSEAoutputformattedforBarplot) <- GSEAoutput$Description
    
    #GSEAoutputformattedforBarplot
    
    #class(GSEAoutput$setSize)
    
    
    #GSEAoutputformattedforBarplot = subset(as.data.frame(GSEAoutput), select = c(ID, Description, setSize, enrichmentScore, pvalue, p.adjust))
    
    #names(GSEAoutputformattedforBarplot)[names(GSEAoutputformattedforBarplot) == 'setSize'] <- 'Count'
    
    #GSEAoutputformattedforBarplot <- GSEAoutput[c("ID", "Description", "setSize", "enrichmentScore", "pvalue", "p.adjust")]
    
    
    
    ## GSEA Dot Plot
    pdf(paste(p, "_GSEA_Dotplot.pdf", sep=""), 
        width = 8, 
        height = 9
    )
    print(dotplot(GSEAoutput, 
                  title = paste(p, "Gene Set Enrichment Analysis") #title
    ))
    dev.off()
    
    
    ## GSEA Heatmap-like Functional Classification
    pdf(paste(p, "_GSEA_Heatmap_like.pdf", sep=""), 
        width = 10, 
        #height = 6
    )
    print(heatplot(GSEAoutput, showCategory = 5,
                   foldChange = geneList, #specify where to pull the FC data from for the heatmap function
                   pvalue = NULL, #place pvalue here (but I could never get it to work FYI)
                   symbol = "rect" #rectangle ["rect"] markers or dots ["dot"]
                   #label_format = "3" #how much space to give the y axis labels
    ) +
      ggplot2::ylab("Enriched GO Terms") + 
      ggplot2::xlab("Genes Implicated") +
      ggplot2::ggtitle(paste("Genes Implicated in Top 5 Enriched GO Terms for", p))
    )
    dev.off()
    
    
    ## GSEA Gene Concept Network
    pdf(paste(p, "_GSEA_Gene_Concept_Network.pdf", sep=""),
        width = 8, 
        height = 8
    )
    print(cnetplot(GSEAoutput, #data
                   layout = "dh", #in order of how much I like them = dh, kk, fr, star, graphopt, grid, mds, gem, circle, drl, lgl, randomly
                   showCategory = 3, #how many terms to display (works from most enriched down)
                   color.params = list(foldChange = NULL, 
                                       edge = TRUE
                                       #category = _, #use to change the color to a specific color for a specific category
                                       #gene = _ #use to change the color to a specific color for a specific gene
                   ),
                   cex.params = list(category_node = 1, #scale by which the node sizes are created
                                     category_label = 0.5, #scale by which the category labels are created
                                     gene_label = 0.3 #scale by which the gene labels are created
                   )
    ) + ggplot2::labs(title = paste("Gene Interaction Network of Top 3 Enriched GO Terms in", p)
                      #subtitle = paste("place your subtitle here")
    )
    )
    dev.off()
    
    
    ## GSEA Linkage Plot
    
    #run the linkage-creation command
    linkage <- pairwise_termsim(GSEAoutput)
    
    #In order to avoid 'subscript out of bounds error' from the linkage plot graphics below
    #     you need to run the following commands to resolve this issue.  To see the explanation of
    #     why this is necessary, please see my electronic lab notebook day 2/24/23 and 3/14/23.
    
    #Determine which rows in the human-readable data frame within 'linkage' contain blanks (*NA*)
    #     These row numbers will be used to correct for the blank column and row headers in the 
    #     non-human-readable 'termsim' object within the 'linkage' object in the next step.
    #     These numbers do correspond 1:1 so this is an appropriate way to do this.
    NAnumbers <- which(is.na(linkage@result$Description))
    
    #Replace any blanks (*NA*) with 'NA' in the 'results' object within 'linkage', but only in column 'Description'
    linkage@result$Description[is.na(linkage@result$Description)] <- "NA"
    
    #Replace any blank (*NA*) row or column names with 'NA' in object 'termsim' within 'linkage'
    row.names(linkage@termsim)[c(NAnumbers)] <- c("NA")
    colnames(linkage@termsim)[c(NAnumbers)] <- c("NA")
    
    #create the linkage plot
    pdf(paste(p, "_GSEA_Linkage_Plot.pdf", sep=""), 
        #width = 12, 
        #height = 7
    )
    print(emapplot(linkage, #data
                   showCategory = 50,
                   cex.params = list(category_label=0.5)
    )
    + ggplot2::labs(title = paste(p, "Gene Set Enrichment Analysis Linkage"),
                    subtitle = paste("(Capped at Top 50 Terms)")
    )
    )
    dev.off()
    
    
    ## GSEA Tree Plot
    
    if(FALSE) {
      
      
      pdf(paste(p, "_GSEA_Tree_Plot.pdf", sep=""), 
          width = 12, 
          height = 6
      )
      print(treeplot(linkage, #your data
                     #showCategory = _, #how many enriched terms to display
                     #label_format = 2, #how many words to allow the outmost label to have before it wraps (although I am not sure if it is really doing this persay idk)
                     hilight.params = list(hilight = TRUE, 
                                           align = "left"), #sets the highlighting of the highlighted region 
                     offset.params = list(bar_tree = rel(1), #bar_tree = rel(1) is a standard node length multiplier with  greater numbers resulting in shorter nodes
                                          tiplab = rel(1.2), #tiplab = rel(1) is a standard setting to set the distance between the leaf label and the dot
                                          extend = 0.3, #sets the distance between branches and leaves on the y-axis
                                          hexpand = 0.3 #sets how squished or expanded the overall graphic is on the x-axis including the tree but not including the outermost label level
                     ),  
                     cluster.params = list(method = "average", #method by which to create the phylogenetic tree
                                           #label_words = "1", #supposed to set how many words are allowed to be listed in the outermost label level but it does not appear to work
                                           label_format = 50 #how many characters are allowed to be printed in the outermost label level before it wraps
                     )
      ) + ggplot2::labs(title = paste("Hierarchial Clustering of Enriched GO Terms in", p),
                        subtitle = paste("(Calculated via Shared Genes Between Terms)")
      )
      )
      dev.off()
      
    } #close the on/off switch for the treeplot
    
  } #close the on/off switch for the GSEA
  
} #close the loop


##### END #####

