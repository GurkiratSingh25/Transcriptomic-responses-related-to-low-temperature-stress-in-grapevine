#### Environmental Setup Steps ####

library(dplyr)
library(forcats)
library(mjcbase)
library(stringr)
library(tidyverse)



#### Set your working directory: ####
setwd("C:/Users/singh.1815/OneDrive - The Ohio State University/Gurkirat Singh/MS Degree outline/Aim1/DESeq/Diamond/orthodiamond")


#### Read in Count Matrices ####

### CabSav
#Read in data
dataCabSav <- read.table("countmatrix_CabSav.txt", header=T,  com='', quote='', check.names=F, sep="\t")

#Name Column 1
names(dataCabSav)[1] <- 'ID_CabSav'

#Only retain columns seen below (those which were PN transcripts aligned to the PN genome)
dataCabSav <- dataCabSav %>% 
  select(ID_CabSav, c("CS_Chl_2":"CS_Fz_8")) 



### GREM4 ---
#Read in data
dataGREM4 <- read.table("countmatrix_GREM.txt", header=T, com='', quote='', check.names=F, sep="\t")

#Name Column 1
names(dataGREM4)[1] <- "ID_GREM"

#Only retain columns seen below (those which were GREM transcripts aligned to the GREM genome)
dataGREM4 <- dataGREM4 %>% 
  select(-c("CS_Chl_2":"CS_Fz_8")) 




#### Read in the Orthology Table ####

### Orthology (diamond BLAST results) ---

#Read in the table (simplified table FYI)
orthology <- read.table("diamondout4.0.blast.orthology_list.txt", header=F, com='', quote='', check.names=F, sep="\t")

#Name Column 1
names(orthology)[1] <- "ID_GREM"

#Name Column 2
names(orthology)[2] <- "ID_CabSav"


#### Create the Inter-Species Comparison Count Table ####

#Take each GREM4 gene name in the orthology table as input individually and search the count matrix GREM4 ID column for that name, if it is found, print to a new object - tmp
#Take each PN gene name in the tmp object as input individually and search the count matrix PN ID column for that name, if it is found, print to the right in a new object - tmp2
#Remove the 'ID_PN' column from the tmp2 object
#Remove the name of 'ID_GREM" from the first column from the tmp2 object
#Remove any rows containing 'NA's in the count matrix portion (does not search the gene names FYI)
tmp <- left_join(orthology, dataGREM4, by = "ID_GREM")
tmp2 <- left_join(tmp, dataCabSav, by = "ID_CabSav") %>% 
  select(-ID_CabSav)
names(tmp2)[1] <- ""
tmp3 <- na.omit(tmp2)

#Export the finalized inter-species count matrix for future use
write.table(tmp3, 
            file = "countmatrix_orthology_GREM4_CabSav_sorted.txt", 
            row.names = FALSE, 
            col.names = TRUE, #necessary to extract the column names from the object, w/o which, it does not
            quote = FALSE, #says to not retain the quotation marks all over
            sep = '\t' #says to make it tab delimited
) 
