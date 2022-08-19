install.packages()
library(SKAT)
#https://cran.r-project.org/web/packages/SKAT/vignettes/SKAT.pdf
#testing SKAT
data(SKAT.example)
names(SKAT.example)
attach(SKAT.example)

memory.limit(size=7000)
library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(readxl)
library(splitstackshape)
library(PheWAS)
library(broom)
library(lubridate)
library(cobalt)
library(MatchIt)
library(DMwR)
library(MASS)
library(ISLR)
library(ggrepel)
library(hrbrthemes)
library(survival)
library(survminer)
library(gridExtra)
library(forestplot)

library(GMMAT)

library('SKAT')
#load the selected patient ID for analysis
#this file was created from 09162021_comorbidity_matrix_06192020
load("selectedsample_EXOME_ratio5_final.RData")
write.table(selectedsample_EXOME_ratio5_final, file = "selectedsample_EXOME_ratio5_final.txt", col.names = F, row.names = F, quote = F)

#four patients removed from the final they are GHS_PT840201_195992214, GHS_PT474208_182822353, GHS_PT228258_180711474, GHS_PT1493628_212823593

library(splitstackshape)


cat("WARNING: this script doesn't perform a sanity check on your data and therefore it's your own responsibility to provide correct input data.\n")
args = unlist(strsplit(commandArgs(trailingOnly = TRUE)," "))
basename = strsplit(basename(args), "\\.")[[1]][1]
dat <- read.table(args[1], header=TRUE, stringsAsFactors=FALSE) 

FAM_90K <- read.table("GHS_Freeze_90.concov.GL.splitMulti.fam", header = FALSE, stringsAsFactors=FALSE)
head(FAM_90K)
FAM_90K <- cSplit(FAM_90K, "V1", "_")
ID_90K <- as.integer(gsub("PT", "", FAM_90K$V1_2))

FAM_175K <- read.table("GHS_Freeze_175.GL.splitmulti_sevengene_12032021.fam", header = FALSE, stringsAsFactors=FALSE)
head(FAM_175K)
FAM_175K <- cSplit(FAM_175K, "V1", "_")
FAM_175K$V1_2 <- gsub("PT", "", FAM_175K$V1_2)
FAM_90K_discovery <- FAM_175K[FAM_175K$V1_2 %in% phenotypes_ratio2$id, ]
FAM_90K_discovery_select <- data.frame(FAM_90K_discovery$V2, FAM_90K_discovery$V2)
colnames(FAM_90K_discovery_select) <- c("FID", "IID")
write.table(FAM_90K_discovery_select, file = "FAM_90K_discovery_select.txt", col.names=T, quote=F, row.names=F)


FAM_175K_discovery <- data.frame(FAM_175K[FAM_175K$V1_2 %in% ID_90K_p123, ]$V2, FAM_175K[FAM_175K$V1_2 %in% ID_90K_p123, ]$V2)
colnames(FAM_175K_discovery) <- c("FID", "IID")

FAM_175K_replication <- data.frame(FAM_175K[FAM_175K$V1_2 %in% ID_85K_p123, ]$V2, FAM_175K[FAM_175K$V1_2 %in% ID_85K_p123, ]$V2)
colnames(FAM_175K_replication) <- c("FID", "IID")

write.table(FAM_175K_discovery, file = "FAM_175K_discovery_select.txt", col.names=T, quote=F, row.names=F)
write.table(FAM_175K_replication, file = "FAM_175K_replication_select.txt", col.names=T, quote=F, row.names=F)

head(dat_RECUR_casecontrol_select_combined.matched) #13350*14
head(dat_RECUR_casecontrol_select_combined) #26593*14

#for pipeline 123 discovery and replication
load("dat_RECUR_casecontrol_select_combined_discoveryreplication_12032021_1to2.RData")
dat_RECUR_casecontrol_select_combined_85K
dat_RECUR_casecontrol_select_combined_90K
ID_90K_p123 <- as.integer(gsub("PT", "", dat_RECUR_casecontrol_select_combined_90K$PT_ID))
ID_85K_p123 <- as.integer(gsub("PT", "", dat_RECUR_casecontrol_select_combined_85K$PT_ID))

#load genotyping file
dat <- read.table("GHS_Freeze_175.GL.splitmulti_sevengene_12032021.raw", header =T, stringsAsFactors = F)

str(dat)
dat$PTID <- dat$FID
dat <- cSplit(dat, "PTID", "_")
dat$PTID_2 <- gsub("PT", "", dat$PTID_2)
dat$PTID_2 <- as.integer(dat$PTID_2)


ID_90K_p123 <- as.integer(gsub("PT", "", ID_90K_p123)) #8235
dat <- dat[dat$PTID_2 %in% ID_90K_p123, ] #8235*78 for discovery p123
dat <- dat[order(dat$PTID_2), ]
colnames(dat)
length(unique(dat$FID))

ID_85K_p123 <- as.integer(gsub("PT", "", ID_85K_p123)) #5115
dat <- dat[dat$PTID_2 %in% ID_85K_p123, ] #5115*78 for replication p123
dat <- dat[order(dat$PTID_2), ]
colnames(dat) 


dat <- dat[, 1:74] #12032021 ratio 2

#make sure the order of PTID is the same as EMMAX kinship
tFAM <- read.table("GHS_Freeze_175.GL.splitmulti_FAM_175K_discovery_select_EXOME2_mind01_hwe2_emmax_2.tfam", header = F, stringsAsFactors = F) #p123 discovery
tFAM <- read.table("GHS_Freeze_175.GL.splitmulti_FAM_175K_replication_select_EXOME2_mind01_hwe2_emmax_2.tfam", header = F, stringsAsFactors = F) #p123 discovery

colnames(tFAM) <- c("FID", "IID", "F", "M", "SEX", "PHENOTYPE")
tFAM$INDEX <- as.integer(rownames(tFAM)) #make sure change chr to integer
tFAM <- cSplit(tFAM, "FID", "_")
tFAM$FID_2 <- as.integer(gsub("PT", "", tFAM$FID_2))
tFAM <- tFAM[order(tFAM$FID_2), ]

dat <- cbind(dat, tFAM[, c("INDEX")])
dat <- dat[order(dat$INDEX), ]

dat <- dat[, 1:74] #12032021 ratio 2
dat <- dat[, 1:39] #p13
colnames(dat)

#real phenotype file
load("phenotypes_matched_ratio2_12032021.RData") #phenotypes_ratio2
str(phenotypes_ratio2$id)

phenotypes_ratio2 <- phenotypes109_age_ratio2

phenotypes_ratio2 <- phenotypes[phenotypes$id %in% ID_90K_p123, ] #8235*1858 for p123 discovery
phenotypes_ratio2 <- phenotypes[phenotypes$id %in% ID_85K_p123, ] #5115*1858 for p123 replication

table(phenotypes_ratio2$group)

str(phenotypes_ratio2)
phenotypes_ratio2 <- phenotypes_ratio2[order(phenotypes_ratio2$id), ]
head(phenotypes_ratio2$id, 10)

phenotypes_ratio2 <- cbind(phenotypes_ratio2, tFAM[, c("INDEX")])
phenotypes_ratio2 <- phenotypes_ratio2[order(phenotypes_ratio2$INDEX), ]

colnames(phenotypes_ratio2)[1859] #group
phenotypes_ratio2 <- phenotypes_ratio2[, 1:1858]
colnames(phenotypes_ratio2)
phenotypes_ratio2$id

#working (as.list)
phen = phenotypes_ratio2$'X496.21'
phen = phenotypes_ratio2$'442.1'

phen = phenotypes_ratio2$X442.1

table(phen)
colnames(dat)
geno = as.matrix(dat[,-c(1:6)]) #Z
dim(geno) 
sum(is.na(geno)) #10

geno
geno[is.na(geno)] <- 2 #change from 0 to 2
obj <- SKAT_Null_Model(phen~1, out_type="D")
obj <- SKAT_NULL_emmaX(phen~dat_RECUR_casecontrol_select_combined_90K$PT_SEX + dat_RECUR_casecontrol_select_combined_90K$PT_RACE, K=K)
obj <- SKAT_NULL_emmaX(phen~dat_RECUR_casecontrol_select_combined_90K$PT_SEX, K=K)

obj <- SKAT_NULL_emmaX(phen~1, K=K)
obj <- SKAT_Null_Model(phen~dat_RECUR_casecontrol_select_combined_90K$PT_SEX + dat_RECUR_casecontrol_select_combined_90K$PT_RACE, out_type="D")



# K: kinship matrix

K <- read.table("GHS_Freeze_175.GL.splitmulti_FAM_175K_discovery_select_EXOME2_mind01_hwe2_emmax_2.hIBS.kinf", header = F, stringsAsFactors = F) #discovery for p123
#K <- read.table("GHS_Freeze_175.GL.splitmulti_FAM_175K_replication_select_EXOME2_mind01_hwe2_emmax_2.hIBS.kinf", header = F, stringsAsFactors = F) #replication for p123
K <- as.matrix(K)
dim(K)
obj <- SKAT_NULL_emmaX(phen~1, K=K)
#SKAT(geno, obj, method="optimal.adj")$p.value
SKAT(geno, obj, method="SKATO")$p.value

##########################################################################################
#create loop for ratio 2 matched
#make sure using data.frame()
phenotypes_ratio2 <- data.frame(phenotypes_ratio2)
phecode <- colnames(phenotypes_ratio2)[2:1858]

#p123
dat_RECUR_casecontrol_select_combined$id <- as.integer(gsub("PT", "", dat_RECUR_casecontrol_select_combined$PT_ID)) 
str(dat_RECUR_casecontrol_select_combined_90K)
dat_RECUR_casecontrol_select_combined_90K <- dat_RECUR_casecontrol_select_combined[dat_RECUR_casecontrol_select_combined$id %in% ID_90K_p123, ]
dat_RECUR_casecontrol_select_combined_90K <- data.frame(dat_RECUR_casecontrol_select_combined_90K)
dat_RECUR_casecontrol_select_combined_90K <- dat_RECUR_casecontrol_select_combined_90K[order(dat_RECUR_casecontrol_select_combined_90K$id), ]
dat_RECUR_casecontrol_select_combined_90K <- cbind(dat_RECUR_casecontrol_select_combined_90K, tFAM[, c("INDEX")])
dat_RECUR_casecontrol_select_combined_90K <- dat_RECUR_casecontrol_select_combined_90K %>%
  arrange(INDEX)
head(dat_RECUR_casecontrol_select_combined_90K)

dat_RECUR_casecontrol_select_combined$id <- as.integer(gsub("PT", "", dat_RECUR_casecontrol_select_combined$PT_ID)) 
str(dat_RECUR_casecontrol_select_combined_85K)
dat_RECUR_casecontrol_select_combined_85K <- dat_RECUR_casecontrol_select_combined[dat_RECUR_casecontrol_select_combined$id %in% ID_85K_p123, ]
dat_RECUR_casecontrol_select_combined_85K <- data.frame(dat_RECUR_casecontrol_select_combined_85K)
dat_RECUR_casecontrol_select_combined_85K <- dat_RECUR_casecontrol_select_combined_85K[order(dat_RECUR_casecontrol_select_combined_85K$id), ]
dat_RECUR_casecontrol_select_combined_85K <- cbind(dat_RECUR_casecontrol_select_combined_85K, tFAM[, c("INDEX")])
dat_RECUR_casecontrol_select_combined_85K <- dat_RECUR_casecontrol_select_combined_85K %>%
  arrange(INDEX)
head(dat_RECUR_casecontrol_select_combined_85K)

save(dat, geno, phecode, phenotypes_ratio2, dat_RECUR_casecontrol_select_combined_90K, file = "12092021_SKAT_inputfile_p123_90K_updated.RData", version = 2)
save(dat, geno, phecode, phenotypes_ratio2, dat_RECUR_casecontrol_select_combined_85K, file = "12092021_SKAT_inputfile_p123_85K_updated.RData", version = 2)

save(dat, geno, phecode, phenotypes_ratio2, dat_RECUR_casecontrol_select_combined_85K, file = "12092021_SKAT_inputfile_p123_85K_updated_phenotypes109.RData", version = 2)
save(dat, geno, phecode, phenotypes_ratio2, dat_RECUR_casecontrol_select_combined_90K, file = "12092021_SKAT_inputfile_p123_90K_updated_phenotypes109.RData", version = 2)

save(dat, geno, phecode, phenotypes_ratio2, dat_RECUR_casecontrol_select_combined_90K, file = "12092021_SKAT_inputfile_p123_90K_updated_phenotypes109_age.RData", version = 2)
save(dat, geno, phecode, phenotypes_ratio2, dat_RECUR_casecontrol_select_combined_85K, file = "12092021_SKAT_inputfile_p123_85K_updated_phenotypes109_age.RData", version = 2)
save(dat, geno, phecode, phenotypes_ratio2, dat_RECUR_casecontrol_select_combined_90K, file = "12092021_SKAT_inputfile_p123_90K_updated_phenotypes109_age60.RData", version = 2)
save(dat, geno, phecode, phenotypes_ratio2, dat_RECUR_casecontrol_select_combined_85K, file = "12092021_SKAT_inputfile_p123_85K_updated_phenotypes109_age60.RData", version = 2)
save(dat, geno, phecode, phenotypes_ratio2, dat_RECUR_casecontrol_select_combined_90K, file = "12092021_SKAT_inputfile_p123_90K_updated_phenotypes109_age55.RData", version = 2)
save(dat, geno, phecode, phenotypes_ratio2, dat_RECUR_casecontrol_select_combined_85K, file = "12092021_SKAT_inputfile_p123_85K_updated_phenotypes109_age55.RData", version = 2)

load("12092021_SKAT_inputfile_p123_90K_updated.RData")
load("12092021_SKAT_inputfile_p123_90K_updated_phenotypes109.RData")
load("12092021_SKAT_inputfile_p123_85K_updated_phenotypes109.RData")
load("12092021_SKAT_inputfile_p123_90K_updated_phenotypes109_age.RData")
load("12092021_SKAT_inputfile_p123_85K_updated_phenotypes109_age.RData")

#do not include any PheCodes with no presence in carriers
SKATO <- list()
SKAT <- list()
BURDEN <- list()
summary_final_SKATO <- NULL
plot_list = list()
phen <- NULL
dim(geno)

for (i in phecode) { 
  print(i)
  phen = phenotypes_ratio2[, i] #change from phenotypes_ratio5
  print(table(phen))
  if(table(phen) == 8235) next # change from 9249 to 4713 for discovery p13; to 4536 for replication; 8235 for discovery p123; 5115 for replicaiton p123
  #obj <- SKAT_NULL_emmaX(phen~dat_RECUR_casecontrol_select_combined_90K$PT_SEX + dat_RECUR_casecontrol_select_combined_90K$PT_RACE, K=K)
  #obj <- SKAT_NULL_emmaX(phen~dat_RECUR_casecontrol_select_combined_85K$PT_SEX + dat_RECUR_casecontrol_select_combined_85K$PT_RACE, K=K)
  #obj <- SKAT_NULL_emmaX(phen~dat_RECUR_casecontrol_select_combined_90K$PT_SEX, K=K)
  obj <- SKAT_Null_Model(phen~dat_RECUR_casecontrol_select_combined_90K$PT_SEX, out_type="D")
  
  tryCatch({
    SKATO[i] <- SKAT(geno, obj, method="SKATO")$p.value
    SKAT[i] <- SKAT(geno, obj, method="SKAT")$p.value
    BURDEN[i] <- SKAT(geno, obj, method="Burden")$p.value
  }, error=function(e){})
}

summary_final_SKATO <- data.frame(do.call(rbind,SKATO))
summary_final_SKAT <- data.frame(do.call(rbind,SKAT))
summary_final_BURDEN <- data.frame(do.call(rbind,BURDEN))

X <- colnames(geno) #68
X

#remove X19.15174391.G.A_G
X <- X[X != "X19.15174391.G.A_G"]  #67

#individal gene
#geno_GENE <- geno[, colnames(geno) %in% X[grepl("^X23.*", X)]] #COL4A1/2
#dim(geno_GENE)

#all gene
geno_GENE <- geno[, colnames(geno) %in% X] 

phecode_select <- c("X425.11", "X425.1", "X425", "X433", "X433.1", "X433.2", "X433.3", "X286.8", "X286.81")

#do not include any PheCodes with no presence in carriers
SKATO <- list()
SKAT <- list()
BURDEN <- list()
summary_final_SKATO <- NULL
plot_list = list()
phen <- NULL
dim(geno_GENE)

for (i in phecode) { 
  print(i)
  phen = phenotypes_ratio2[, i] #change from phenotypes_ratio5
  print(table(phen))
  if(table(phen) == 8235) next # change from 9249 to 4713 for discovery p13; to 4536 for replication; 8235 for discovery p123; 5115 for replicaiton p123
  #obj <- SKAT_NULL_emmaX(phen~dat_RECUR_casecontrol_select_combined_90K$PT_SEX + dat_RECUR_casecontrol_select_combined_90K$PT_RACE, K=K)
  #obj <- SKAT_NULL_emmaX(phen~dat_RECUR_casecontrol_select_combined_85K$PT_SEX + dat_RECUR_casecontrol_select_combined_85K$PT_RACE, K=K)
  obj <- SKAT_NULL_emmaX(phen~dat_RECUR_casecontrol_select_combined_90K$PT_SEX, K=K)
  #obj <- SKAT_Null_Model_ChrX(phen~dat_RECUR_casecontrol_select_combined_85K$PT_SEX, SexVar = "PT_SEX", out_type="D")
  #obj <- SKAT_Null_Model(phen~dat_RECUR_casecontrol_select_combined_85K$PT_SEX, out_type="D")
  
  tryCatch({
    SKATO[i] <- SKAT(geno_GENE, obj, kernel = "linear.weighted", method="SKATO")$p.value
    SKAT[i] <- SKAT(geno_GENE, obj, kernel = "linear.weighted", method="SKAT")$p.value
    BURDEN[i] <- SKAT(geno_GENE, obj, kernel = "linear.weighted", method="Burden")$p.value
  }, error=function(e){})
}

summary_final_SKATO <- data.frame(do.call(rbind,SKATO))
summary_final_SKAT <- data.frame(do.call(rbind,SKAT))
summary_final_BURDEN <- data.frame(do.call(rbind,BURDEN))

save(summary_final_BURDEN, summary_final_SKAT, summary_final_SKATO, file = "summary_final_12032021_matched_ratio2_85K_ALLSIXGENES_linearweighted_ALL.RData", version = 2)
#save(summary_final_BURDEN, summary_final_SKAT, summary_final_SKATO, file = "summary_final_12032021_matched_ratio2_90K_age_chrX_linearweighted.RData", version = 2)

#do not include any PheCodes with no presence in carriers
phenotypes_ratio2 <- data.frame(phenotypes_ratio2)

SKATO <- list()
SKAT <- list()
BURDEN <- list()
summary_final_SKATO <- NULL
plot_list = list()
phen <- NULL
dim(geno_GENE)

for (i in phecode) { 
  print(i)
  phen = phenotypes_ratio2[, i] #change from phenotypes_ratio5
  print(table(phen))
  if(table(phen) == 5115) next # change from 9249 to 4713 for discovery p13; to 4536 for replication; 8235 for discovery p123; 5115 for replicaiton p123
  #obj <- SKAT_NULL_emmaX(phen~dat_RECUR_casecontrol_select_combined_90K$PT_SEX + dat_RECUR_casecontrol_select_combined_90K$PT_RACE, K=K)
  #obj <- SKAT_NULL_emmaX(phen~dat_RECUR_casecontrol_select_combined_85K$PT_SEX + dat_RECUR_casecontrol_select_combined_85K$PT_RACE, K=K)
  obj <- SKAT_NULL_emmaX(phen~dat_RECUR_casecontrol_select_combined_85K$PT_SEX, K=K)
  #obj <- SKAT_Null_Model_ChrX(phen~dat_RECUR_casecontrol_select_combined_85K$PT_SEX, SexVar = "PT_SEX", out_type="D")
  #obj <- SKAT_Null_Model(phen~dat_RECUR_casecontrol_select_combined_85K$PT_SEX, out_type="D")
  
  tryCatch({
    #SKATO[i] <- SKAT(geno_GENE, obj, kernel = "linear.weighted", weights.beta = c(1,25), weights=NULL, method ="SKATO")$p.value
    #SKATO[i] <- SKAT(geno_GENE, obj, kernel = "linear.weighted", weights.beta = c(1,50), weights=NULL, method ="SKATO")$p.value
    #SKATO[i] <- SKAT(geno_GENE, obj, kernel = "linear.weighted", weights.beta = c(1,1), weights=NULL, method ="SKATO")$p.value
    #SKATO[i] <- SKAT(geno_GENE, obj, kernel = "linear.weighted", weights.beta = c(0.5,0.5), weights=NULL, method ="SKATO")$p.value
    SKATO[i] <- SKAT(geno_GENE, obj, kernel = "linear.weighted", weights=weights_55, method ="SKATO")$p.value
    #SKATO[i] <- SKAT(geno_GENE, obj, kernel = "linear.weighted", weights.beta = c(1,200), weights=NULL, method ="SKATO")$p.value
    #SKATO[i] <- SKAT(geno_GENE, obj, kernel = "linear.weighted", weights=weights_CADD, method ="SKATO")$p.value
    #SKATO[i] <- SKAT(geno_GENE, obj, kernel = "linear.weighted", method="SKATO")$p.value
  }, error=function(e){})
}

summary_final_SKATO <- data.frame(do.call(rbind,SKATO))
#summary_final_SKATO_beta11 <- summary_final_SKATO
#summary_final_SKATO_beta125 <- summary_final_SKATO
#summary_final_SKATO_beta150 <- summary_final_SKATO
summary_final_SKATO_beta55 <- summary_final_SKATO
#summary_final_SKATO_beta1200 <- summary_final_SKATO
#summary_final_SKATO_CADD <- summary_final_SKATO


#save(summary_final_SKATO_beta55, file = "summary_final_12032021_matched_ratio2_85K_age_ALLSIXGENES_linearweighted_beta55_ALL_fixedweights.RData", version = 2)
save(summary_final_SKATO_beta55, file = "summary_final_12032021_matched_ratio2_85K_age_ALLSIXGENES_ALLPHENCODE_linearweighted_beta55_ALL_fixedweights.RData", version = 2)

#save(summary_final_SKATO_beta55, file = "summary_final_12032021_matched_ratio2_90K_age_ALLSIXGENES_ALLPHENCODE_linearweighted_beta55_ALL.RData", version = 2)
#save(summary_final_SKATO_beta55, file = "summary_final_12032021_matched_ratio2_90K_age_ALLSIXGENES_linearweighted_beta55_ALL.RData", version = 2)

#save(summary_final_SKATO_beta55, file = "summary_final_12032021_matched_ratio2_90K_chr23_linearweighted_beta55_ALL_fixedweights.RData", version = 2)

#save(summary_final_SKATO_beta125, file = "summary_final_12032021_matched_ratio2_85K_ALLSIXGENES_linearweighted_beta125_ALL_fixedweights.RData", version = 2)


#save(summary_final_SKATO_CADD, file = "summary_final_12032021_matched_ratio2_85K_ALLSIXGENES_linearweighted_CADD_ALL.RData", version = 2)
#save(summary_final_SKATO_CADD, file = "summary_final_12032021_matched_ratio2_90K_ALLSIXGENES_linearweighted_CADD_ALL.RData", version = 2)

#save(summary_final_SKATO_beta1200, file = "summary_final_12032021_matched_ratio2_85K_ALLSIXGENES_linearweighted_beta1200_ALL.RData", version = 2)
#save(summary_final_SKATO_beta1200, file = "summary_final_12032021_matched_ratio2_90K_ALLSIXGENES_linearweighted_beta1200_ALL.RData", version = 2)



#save(summary_final_SKATO_beta55, file = "summary_final_12032021_matched_ratio2_85K_ALLSIXGENES_ALLPHENCODE_linearweighted_beta55_ALL_fixedweights.RData", version = 2)
#save(summary_final_SKATO_beta11, file = "summary_final_12032021_matched_ratio2_90K_chr23_linearweighted_beta11_ALL.RData", version = 2)

#save(summary_final_SKATO_beta150, file = "summary_final_12032021_matched_ratio2_90K_chr23_linearweighted_beta150_ALL.RData", version = 2)
#save(summary_final_SKATO_beta150, file = "summary_final_12032021_matched_ratio2_85K_ALLGENES_ALLPHECODES_linearweighted_beta150_ALL.RData", version = 2)

save(summary_final_SKATO_beta55, file = "summary_final_12032021_matched_ratio2_85K_ALLSIXGENES_linearweighted_beta55_ALL_fixedweights.RData", version = 2)
#save(summary_final_SKATO_beta55, file = "summary_final_12032021_matched_ratio2_90K_ALLSIXGENES_linearweighted_beta55_ALL.RData", version = 2)

save(summary_final_SKATO, summary_final_SKATO_beta125, summary_final_SKATO_beta150, file = "summary_final_12032021_matched_ratio2_90K_ALLSIXGENES_linearweighted_beta_ALL.RData", version = 2)
save(summary_final_SKATO, summary_final_SKATO_beta125, summary_final_SKATO_beta150, file = "summary_final_12032021_matched_ratio2_85K_ALLSIXGENES_linearweighted_beta_ALL.RData", version = 2)

save(summary_final_SKATO, summary_final_SKATO_beta125, summary_final_SKATO_beta150, file = "summary_final_12032021_matched_ratio2_90K_age_ALLSIXGENES_linearweighted_beta_ALL.RData", version = 2)
save(summary_final_SKATO, summary_final_SKATO_beta125, summary_final_SKATO_beta150, file = "summary_final_12032021_matched_ratio2_85K_age_ALLSIXGENES_linearweighted_beta_ALL.RData", version = 2)


load("summary_final_12032021_matched_ratio2_90K_ALLSIXGENES_linearweighted_beta_ALL.RData")
load("summary_final_12032021_matched_ratio2_85K_ALLSIXGENES_linearweighted_beta_ALL.RData")

load("summary_final_12032021_matched_ratio2_85K_chr3_linearweighted_beta55_ALL_fixedweights.RData")
summary_final_SKATO_beta55_TREX1 <- summary_final_SKATO_beta55
colnames(summary_final_SKATO_beta55_TREX1) <- "pvalue"
summary_final_SKATO_beta55_TREX1$phecode <- rownames(summary_final_SKATO_beta55_TREX1)
summary_final_SKATO_beta55_TREX1$gene <- "TREX1"
summary_final_SKATO_beta55_TREX1$distribution <- "alpha=0.5, beta=0.5"
summary_final_SKATO_beta55_TREX1


load("summary_final_12032021_matched_ratio2_85K_chr10_linearweighted_beta55_ALL_fixedweights.RData")
summary_final_SKATO_beta55_HTRA1 <- summary_final_SKATO_beta55
colnames(summary_final_SKATO_beta55_HTRA1) <- "pvalue"
summary_final_SKATO_beta55_HTRA1$phecode <- rownames(summary_final_SKATO_beta55_HTRA1)
summary_final_SKATO_beta55_HTRA1$gene <- "HTRA1"
summary_final_SKATO_beta55_HTRA1$distribution <- "alpha=0.5, beta=0.5"
summary_final_SKATO_beta55_HTRA1

load("summary_final_12032021_matched_ratio2_85K_chr13_linearweighted_beta55_ALL_fixedweights.RData")
summary_final_SKATO_beta55_COL4A1 <- summary_final_SKATO_beta55
colnames(summary_final_SKATO_beta55_COL4A1) <- "pvalue"
summary_final_SKATO_beta55_COL4A1$phecode <- rownames(summary_final_SKATO_beta55_COL4A1)
summary_final_SKATO_beta55_COL4A1$gene <- "COL4A1/2"
summary_final_SKATO_beta55_COL4A1$distribution <- "alpha=0.5, beta=0.5"
summary_final_SKATO_beta55_COL4A1

load("summary_final_12032021_matched_ratio2_85K_chr17_linearweighted_beta55_ALL_fixedweights.RData")
summary_final_SKATO_beta55_CTC1 <- summary_final_SKATO_beta55
colnames(summary_final_SKATO_beta55_CTC1) <- "pvalue"
summary_final_SKATO_beta55_CTC1$phecode <- rownames(summary_final_SKATO_beta55_CTC1)
summary_final_SKATO_beta55_CTC1$gene <- "CTC1"
summary_final_SKATO_beta55_CTC1$distribution <- "alpha=0.5, beta=0.5"
summary_final_SKATO_beta55_CTC1

load("summary_final_12032021_matched_ratio2_85K_chr19_linearweighted_beta55_ALL_fixedweights.RData")
summary_final_SKATO_beta55_NOTCH3 <- summary_final_SKATO_beta55
colnames(summary_final_SKATO_beta55_NOTCH3) <- "pvalue"
summary_final_SKATO_beta55_NOTCH3$phecode <- rownames(summary_final_SKATO_beta55_NOTCH3)
summary_final_SKATO_beta55_NOTCH3$gene <- "NOTCH3"
summary_final_SKATO_beta55_NOTCH3$distribution <- "alpha=0.5, beta=0.5"
summary_final_SKATO_beta55_NOTCH3

load("summary_final_12032021_matched_ratio2_85K_chr23_linearweighted_beta55_ALL_fixedweights.RData")
summary_final_SKATO_beta55_GLA <- summary_final_SKATO_beta55
colnames(summary_final_SKATO_beta55_GLA) <- "pvalue"
summary_final_SKATO_beta55_GLA$phecode <- rownames(summary_final_SKATO_beta55_GLA)
summary_final_SKATO_beta55_GLA$gene <- "GLA"
summary_final_SKATO_beta55_GLA$distribution <- "alpha=0.5, beta=0.5"
summary_final_SKATO_beta55_GLA

#load("summary_final_12032021_matched_ratio2_90K_ALLSIXGENES_linearweighted_beta55_ALL.RData")
load("summary_final_12032021_matched_ratio2_85K_ALLSIXGENES_linearweighted_beta55_ALL_fixedweights.RData")
summary_final_SKATO_beta55_ALL_fixedweights <- summary_final_SKATO_beta55
colnames(summary_final_SKATO_beta55_ALL_fixedweights) <- "pvalue"
summary_final_SKATO_beta55_ALL_fixedweights$phecode <- rownames(summary_final_SKATO_beta55_ALL_fixedweights)
summary_final_SKATO_beta55_ALL_fixedweights$gene <- "ALLGENES"
summary_final_SKATO_beta55_ALL_fixedweights$distribution <- "alpha=0.5, beta=0.5"
summary_final_SKATO_beta55_ALL_fixedweights

foo <- list(summary_final_SKATO_beta55_ALL_fixedweights, summary_final_SKATO_beta55_COL4A1, summary_final_SKATO_beta55_NOTCH3, summary_final_SKATO_beta55_HTRA1, summary_final_SKATO_beta55_TREX1, summary_final_SKATO_beta55_CTC1, summary_final_SKATO_beta55_GLA)
summary_final_SKATO_beta55_ALL_summary <- do.call(rbind, foo)

#summary_final_SKATO_beta55_ALL_summary$cohort <- "discovery"
#summary_final_SKATO_beta55_ALL_summary_discovery <- summary_final_SKATO_beta55_ALL_summary

summary_final_SKATO_beta55_ALL_summary$cohort <- "replication"
summary_final_SKATO_beta55_ALL_summary_replication <- summary_final_SKATO_beta55_ALL_summary

summary_final_SKATO_beta55_ALL_summary<- rbind(summary_final_SKATO_beta55_ALL_summary_discovery, summary_final_SKATO_beta55_ALL_summary_replication)
summary_final_SKATO_beta55_ALL_summary

#for v1.2
library(dplyr)
phecode_mapping <- read.csv("Phecode_map_v1_2_icd10cm_beta.csv", header = T, stringsAsFactors = F)
phecode_mapping_unique <- phecode_mapping %>%
  dplyr::select(phecode, phecode_str, exclude_range, exclude_name) %>%
  distinct(phecode, .keep_all = T)  #1755*4
str(phecode_mapping_unique)
phecode_mapping_unique$phecode <- as.character(paste0("X", phecode_mapping_unique$phecode))

unique(phecode_mapping_unique$exclude_name)
# [1] "injuries & poisonings"   "neoplasms"               "circulatory system"      "digestive"              
# [5] "genitourinary"           "symptoms"                "sense organs"            "mental disorders"       
# [9] "neurological"            "endocrine/metabolic"     "hematopoietic"           "musculoskeletal"        
# [13] "respiratory"             ""                        "infectious diseases"     "pregnancy complications"
# [17] "NULL"                    "dermatologic"            "congenital anomalies"  

phecode_mapping_unique$phecode
summary_final_SKATO_beta55_ALL_summary <- merge(summary_final_SKATO_beta55_ALL_summary, phecode_mapping_unique, by = "phecode")
xlabels <- unique(summary_final_SKATO_beta55_ALL_summary$phecode_str)

library(ggplot2)
# The palette with grey: , "#CC79A7"
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")
tiff("Summary_phenotypes_merge_anno_SKATO_rare_topbottom_union_ratio2_candidatephecodes_beta55.tiff", units="in", width=10, height=8, res=600)
ggplot(summary_final_SKATO_beta55_ALL_summary, aes(x = as.factor(phecode), y = -log10(pvalue))) + 
  #geom_pointrange(aes(col = PHENOTYPE, shape = as.factor(Direction)), position=position_dodge(width=0.8)) + 
  #geom_point(aes(shape = as.factor(Direction), color = PHENOTYPE, y = estimate, size = -log10(p.value)), stat="identity", position = position_dodge(0.8))  +
  geom_point(aes(shape = as.factor(cohort), color = gene, y = -log10(pvalue), size = -log10(pvalue)), stat="identity", position = position_dodge(0.8)) +
  #scale_shape_manual(values=c("\u25BA","\u25C4","\u25BC","\u25B2")) +
  scale_color_manual(values=cbPalette) +
  scale_shape_manual(values=c("\u25BA","\u25C4")) +
  scale_x_discrete(labels= xlabels) +
  geom_hline(aes(yintercept = -log10(0.05)), linetype='dashed', col = 'black', width = 2) +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme_bw() +
  geom_text(data = subset(summary_final_SKATO_beta55_ALL_summary, pvalue < 0.018), 
            aes(phecode, -log10(pvalue), label = gene), vjust = -1, size = 4, angle = 0) +
  #facet_wrap(. ~ Header.x, scales="free") + 
  #theme(strip.text.x = element_text(size = 6, colour = "black")) +
  theme(text = element_text(size = 12),
        axis.text.x=element_text(angle = 45, vjust = 1, hjust=1)) + 
  guides(color = guide_legend(title="Gene", override.aes = list(size = 5))) +
  guides(shape = guide_legend(title="Cohort", override.aes = list(size = 5))) +
  ylab("-log10(p.value)") +
  xlab("phecode")
  #ggtitle("Summary_phenotypes_merge_anno_SKATO_rare_topbottom_union_ratio2")
dev.off()


load("summary_final_12032021_matched_ratio2_90K_ALLSIXGENES_linearweighted_beta_ALL.RData")
load("summary_final_12032021_matched_ratio2_85K_ALLSIXGENES_linearweighted_beta150_ALL_fixedweights.RData")
summary_final_SKATO_beta150_ALL <- summary_final_SKATO_beta150
colnames(summary_final_SKATO_beta150_ALL) <- "pvalue"
summary_final_SKATO_beta150_ALL$phecode <- rownames(summary_final_SKATO_beta150_ALL)
summary_final_SKATO_beta150_ALL$gene <- "ALLGENES"
summary_final_SKATO_beta150_ALL$distribution <- "alpha=1, beta=50"
summary_final_SKATO_beta150_ALL

load("summary_final_12032021_matched_ratio2_85K_ALLSIXGENES_linearweighted_beta125_ALL_fixedweights.RData")
summary_final_SKATO_beta125_ALL <- summary_final_SKATO_beta125
colnames(summary_final_SKATO_beta125_ALL) <- "pvalue"
summary_final_SKATO_beta125_ALL$phecode <- rownames(summary_final_SKATO_beta125_ALL)
summary_final_SKATO_beta125_ALL$gene <- "ALLGENES"
summary_final_SKATO_beta125_ALL$distribution <- "alpha=1, beta=25"
summary_final_SKATO_beta125_ALL

load("summary_final_12032021_matched_ratio2_90K_ALLSIXGENES_linearweighted_beta11_ALL.RData")
load("summary_final_12032021_matched_ratio2_85K_ALLSIXGENES_linearweighted_beta11_ALL.RData")
summary_final_SKATO_beta11_ALL <- summary_final_SKATO_beta11
colnames(summary_final_SKATO_beta11_ALL) <- "pvalue"
summary_final_SKATO_beta11_ALL$phecode <- rownames(summary_final_SKATO_beta11_ALL)
summary_final_SKATO_beta11_ALL$gene <- "ALLGENES"
summary_final_SKATO_beta11_ALL$distribution <- "alpha=1, beta=1"
summary_final_SKATO_beta11_ALL

load("summary_final_12032021_matched_ratio2_90K_ALLSIXGENES_linearweighted_beta55_ALL.RData")
load("summary_final_12032021_matched_ratio2_85K_ALLSIXGENES_linearweighted_beta55_ALL_fixedweights.RData")
summary_final_SKATO_beta55_ALL <- summary_final_SKATO_beta55
colnames(summary_final_SKATO_beta55_ALL) <- "pvalue"
summary_final_SKATO_beta55_ALL$phecode <- rownames(summary_final_SKATO_beta55_ALL)
summary_final_SKATO_beta55_ALL$gene <- "ALLGENES"
summary_final_SKATO_beta55_ALL$distribution <- "alpha=0.5, beta=0.5"
summary_final_SKATO_beta55_ALL

load("summary_final_12032021_matched_ratio2_90K_ALLSIXGENES_linearweighted_beta1200_ALL.RData")
load("summary_final_12032021_matched_ratio2_85K_ALLSIXGENES_linearweighted_beta1200_ALL_fixedweights.RData")
summary_final_SKATO_beta1200_ALL <- summary_final_SKATO_beta1200
colnames(summary_final_SKATO_beta1200_ALL) <- "pvalue"
summary_final_SKATO_beta1200_ALL$phecode <- rownames(summary_final_SKATO_beta1200_ALL)
summary_final_SKATO_beta1200_ALL$gene <- "ALLGENES"
summary_final_SKATO_beta1200_ALL$distribution <- "alpha=1, beta=200"
summary_final_SKATO_beta1200_ALL

load("summary_final_12032021_matched_ratio2_90K_ALLSIXGENES_linearweighted_CADD_ALL.RData")
load("summary_final_12032021_matched_ratio2_85K_ALLSIXGENES_linearweighted_CADD_ALL.RData")
summary_final_SKATO_CADD_ALL <- summary_final_SKATO_CADD
colnames(summary_final_SKATO_CADD_ALL) <- "pvalue"
summary_final_SKATO_CADD_ALL$phecode <- rownames(summary_final_SKATO_CADD_ALL)
summary_final_SKATO_CADD_ALL$gene <- "ALLGENES"
summary_final_SKATO_CADD_ALL$distribution <- "CADD PHRED"
summary_final_SKATO_CADD_ALL

foo <- list(summary_final_SKATO_beta11_ALL, summary_final_SKATO_beta125_ALL, summary_final_SKATO_beta150_ALL, summary_final_SKATO_beta55_ALL, summary_final_SKATO_beta1200_ALL, summary_final_SKATO_CADD_ALL)
summary_final_SKATO_beta_ALL_summary <- do.call(rbind, foo)

#summary_final_SKATO_beta_ALL_summary$cohort <- "discovery"
#summary_final_SKATO_beta_ALL_summary_discovery <- summary_final_SKATO_beta_ALL_summary

summary_final_SKATO_beta_ALL_summary$cohort <- "replication"
summary_final_SKATO_beta_ALL_summary_replication <- summary_final_SKATO_beta_ALL_summary

summary_final_SKATO_beta_ALL_summary<- rbind(summary_final_SKATO_beta_ALL_summary_discovery, summary_final_SKATO_beta_ALL_summary_replication)
summary_final_SKATO_beta_ALL_summary

summary_final_SKATO_beta_ALL_summary <- merge(summary_final_SKATO_beta_ALL_summary, phecode_mapping_unique, by = "phecode")
summary_final_SKATO_beta_ALL_summary

save(summary_final_SKATO_beta150_ALL_summary, summary_final_SKATO_beta_ALL_summary, file = "summary_final_SKATO_beta_ALL_summary.RData", version = 2)
summary_final_SKATO_beta_ALL_summary

save(summary_final_SKATO_beta55_ALL_summary, summary_final_SKATO_beta_ALL_summary, file = "summary_final_SKATO_beta55_ALL_summary.RData", version = 2)
summary_final_SKATO_beta_ALL_summary

save(summary_final_SKATO_beta_ALL_summary, file = "summary_final_SKATO_beta_ALL_summary_addCADD_fixedweights.RData", version = 2)
xlabels <- unique(summary_final_SKATO_beta_ALL_summary$phecode_str)

tiff("Summary_phenotypes_merge_anno_SKATO_rare_topbottom_union_ratio2_candidatephecodes_betadistribution_addCADD_fixweights.tiff", units="in", width=10, height=8, res=600)
ggplot(summary_final_SKATO_beta_ALL_summary, aes(x = as.factor(phecode), y = -log10(pvalue))) + 
  geom_point(aes(shape = as.factor(cohort), color = distribution, y = -log10(pvalue), size = -log10(pvalue)), stat="identity", position = position_dodge(0.8)) +
  #scale_shape_manual(values=c("\u25BA","\u25C4","\u25BC","\u25B2")) +
  scale_shape_manual(values=c("\u25BA","\u25C4")) +
  scale_x_discrete(labels= xlabels) +
  geom_hline(aes(yintercept = -log10(0.05)), linetype='dashed', col = 'black', width = 1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme_bw() +
  geom_text(data = subset(summary_final_SKATO_beta_ALL_summary, pvalue < 0.0006), 
            aes(phecode, -log10(pvalue), label = distribution), position=position_jitter(width=0,height=0.4), vjust = 0, size = 3, angle = 0) +
  #facet_wrap(. ~ Header.x, scales="free") + 
  #theme(strip.text.x = element_text(size = 6, colour = "black")) +
  theme(text = element_text(size = 12),
        axis.text.x=element_text(angle = 45, vjust = 1, hjust=1)) + 
  guides(color = guide_legend(title="Weights", override.aes = list(size = 5))) +
  guides(shape = guide_legend(title="Cohort", override.aes = list(size = 5))) +
  ylab("-log10(p.value)") +
  xlab("phecode")
#ggtitle("Summary_phenotypes_merge_anno_SKATO_rare_topbottom_union_ratio2")
dev.off()

#this common disease is for all based on frequency in both discovery and replication
commondisease <- read.table(file = "Common_disease.txt", header = TRUE, stringsAsFactors=FALSE)

#for v1.2
phecode_mapping <- read.csv("Phecode_map_v1_2_icd10cm_beta.csv", header = T, stringsAsFactors = F)
phecode_mapping_unique <- phecode_mapping %>%
  dplyr::select(phecode, phecode_str, exclude_range, exclude_name) %>%
  distinct(phecode, .keep_all = T)  #1755*4
str(phecode_mapping_unique)
phecode_mapping_unique$phecodeX <- as.character(paste0("X", phecode_mapping_unique$phecode))

unique(phecode_mapping_unique$exclude_name)

#identify common disease from control group (n=711) remove phecodes not belonging to any category
phecode_mapping_unique_commondisease <- phecode_mapping_unique[!phecode_mapping_unique$exclude_name %in% c("", "NULL"), ]
phecode_mapping_unique_commondisease <- phecode_mapping_unique_commondisease[phecode_mapping_unique_commondisease$phecodeX %in% commondisease$phecodeX, ] #711
head(phecode_mapping_unique_commondisease)

#read in 'Summary_phenotypes_merge_anno_common' and 'Summary_phenotypes_merge_anno_rare', two objects
load("Summary_phenotypes109_merge_anno_commonrare_delta_p123_85K_phenotypes109_age.RData")
load("Summary_phenotypes109_merge_anno_commonrare_delta_p123_90K_phenotypes109_age.RData")
load("Summary_phenotypes109_merge_anno_commonrare_delta_p123_85K_phenotypes109.RData")
load("Summary_phenotypes109_merge_anno_commonrare_delta_p123_90K_phenotypes109.RData")
unique(Summary_phenotypes_merge_anno_common$exclude_name)
Summary_phenotypes_merge_anno_common <- Summary_phenotypes_merge_anno_common[!Summary_phenotypes_merge_anno_common$exclude_name %in% c("", "NULL", NA), ]  #695
Summary_phenotypes_merge_anno_common$phecode <- paste0("X", Summary_phenotypes_merge_anno_common$phecode)

#n=696 PheCodes having results
Summary_phenotypes_merge_anno_common <- Summary_phenotypes_merge_anno[Summary_phenotypes_merge_anno$phecode %in% phecode_mapping_unique_commondisease$phecode, ] #696
head(Summary_phenotypes_merge_anno_common)
Summary_phenotypes_merge_anno_common$phecode <- paste0("X", Summary_phenotypes_merge_anno_common$phecode)

#n=224/290 PheCodes having PRE > 0 for 90K/85K; 208/695 for discovery and 248/580 for replication
Summary_phenotypes_merge_anno_common_top <- Summary_phenotypes_merge_anno_common[Summary_phenotypes_merge_anno_common$delta_change >0, ]


#load("summary_final_12032021_matched_ratio2_90K_ALLGENES_ALLPHECODES_linearweighted_beta150_ALL.RData")
#load("summary_final_12032021_matched_ratio2_85K_ALLGENES_ALLPHECODES_linearweighted_beta150_ALL.RData")

load("summary_final_12032021_matched_ratio2_90K_ALLSIXGENES_ALLPHENCODE_linearweighted_beta55_ALL.RData")
load("summary_final_12032021_matched_ratio2_85K_ALLSIXGENES_ALLPHENCODE_linearweighted_beta55_ALL_fixedweights.RData")

colnames(summary_final_SKATO_beta55) <- "pvalue"
summary_final_SKATO_beta55$phecode <- rownames(summary_final_SKATO_beta55)
summary_final_SKATO_beta55 <- merge(Summary_phenotypes_merge_anno_common_top, summary_final_SKATO_beta55, by = "phecode")
summary_final_SKATO_beta55_85K <- summary_final_SKATO_beta55
summary_final_SKATO_beta55_90K <- summary_final_SKATO_beta55

#put results from both cohort together 
summary_final_SKATO_beta55_90K
summary_final_SKATO_beta55_90K$cohort <- "discovery"
summary_final_SKATO_beta55_85K
summary_final_SKATO_beta55_85K$cohort <- "replication"
Summary_phenotypes_merge_anno_SKATO_beta55_common_90K85K <- rbind(summary_final_SKATO_beta55_90K, summary_final_SKATO_beta55_85K) %>%
  arrange(exclude_name)
dim(Summary_phenotypes_merge_anno_SKATO_beta55_common_90K85K) #1275*15

unique(Summary_phenotypes_merge_anno_SKATO_beta55_common_90K85K$cohort)
summary_final_SKATO_beta55_final <- merge(Summary_phenotypes_merge_anno_SKATO_beta55_common_90K85K_top[Summary_phenotypes_merge_anno_SKATO_beta55_common_90K85K_top$cohort == "discovery", ], Summary_phenotypes_merge_anno_SKATO_beta55_common_90K85K_top[Summary_phenotypes_merge_anno_SKATO_beta55_common_90K85K_top$cohort == "replication", ], by = "phecode") #208/248

summary_final_SKATO_beta55_final_replicated <- summary_final_SKATO_beta55_final %>%
  filter(pvalue.x < 0.05) %>%
  filter(pvalue.y < 0.05)

summary_final_SKATO_beta55_final_replicated
0

Summary_phenotypes_merge_anno_SKATO_beta55_common_90K85K_top <- Summary_phenotypes_merge_anno_SKATO_beta55_common_90K85K[Summary_phenotypes_merge_anno_SKATO_beta55_common_90K85K$delta_change > 0, ]  #456*15
dim(Summary_phenotypes_merge_anno_SKATO_beta55_common_90K85K_top)

library(ggplot2)
tiff("Summary_phenotypes_merge_anno_SKATO_common_top_phencode_109_ratio2_90K85K_ALLSIXGENES_linearweighted_beta55_commonfromeach2.tiff", units="in", width=24, height=16, res=600)
ggplot(Summary_phenotypes_merge_anno_SKATO_beta55_common_90K85K_top, aes(x = as.factor(phecode), y = -log10(pvalue))) + 
  #geom_pointrange(aes(col = PHENOTYPE, shape = as.factor(cohort)), position=position_dodge(width=0.8)) + 
  #geom_point(aes(shape = as.factor(cohort), color = PHENOTYPE, y = estimate, size = -log10(pvalue)), stat="identity", position = position_dodge(0.8))  +
  geom_point(aes(shape = as.factor(cohort), color = exclude_name, y = -log10(pvalue), size = abs(delta_change)), stat="identity", position = position_dodge(0.8))  +
  #ylim(0, 0) +
  #scale_y_log10(limits = c(1,30)) + 
  #scale_shape_manual(values=c("\u25BA","\u25C4","\u25BC","\u25B2"))
  scale_shape_manual(values=c("\u25BA","\u25C4")) +
  #geom_hline(aes(yintercept = -log10(0.05)), linetype='dashed', col = 'black', width = 1) +
  geom_hline(aes(yintercept = -log10(0.006152585)), linetype='dashed', col = 'red', width = 1) +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme_bw() +
  geom_text(data = subset(Summary_phenotypes_merge_anno_SKATO_beta55_common_90K85K, pvalue < 0.006152585), 
            aes(phecode, -log10(pvalue), label = phecode_str), position=position_jitter(width=0,height=0.05), hjust = 0.5, vjust = 0, size = 4.5, angle = 0) +
  #facet_wrap(. ~ Header.x, scales="free") + 
  #theme(strip.text.x = element_text(size = 6, colour = "black")) +
  theme(text = element_text(size = 12),
        axis.text.x=element_text(angle = 45, vjust = 1, hjust = 1)) + 
  guides(color = guide_legend(title="Disease Category", override.aes = list(size = 4))) +
  guides(shape = guide_legend(title="Cohort", override.aes = list(size = 4))) +
  guides(size = guide_legend(title="Percent Relative Effect")) +
  ylab("-log10(pvalue) in log scale") +
  xlab("phecode")
  #ggtitle("Summary_phenotypes_merge_anno_SKATO_common_90K85K")
dev.off()

dev.off()

#BH correction
#https://sphweb.bumc.bu.edu/otlt/MPH-Modules/BS/R/R-Manual/R-Manual22.html

load("summary_final_12032021_matched_ratio2_90K_ALLSIXGENES_ALLPHENCODE_linearweighted_beta55_ALL.RData")
load("summary_final_12032021_matched_ratio2_85K_ALLSIXGENES_ALLPHENCODE_linearweighted_beta55_ALL_fixedweights.RData")

colnames(summary_final_SKATO_beta55) <- "pvalue"
summary_final_SKATO_beta55$phecode <- rownames(summary_final_SKATO_beta55)

#read in 'Summary_phenotypes_merge_anno_common' and 'Summary_phenotypes_merge_anno_rare', two objects
load("Summary_phenotypes109_merge_anno_commonrare_delta_p123_85K_phenotypes109_age.RData")
load("Summary_phenotypes109_merge_anno_commonrare_delta_p123_90K_phenotypes109_age.RData")
load("Summary_phenotypes109_merge_anno_commonrare_delta_p123_85K_phenotypes109.RData")
load("Summary_phenotypes109_merge_anno_commonrare_delta_p123_90K_phenotypes109.RData")
unique(Summary_phenotypes_merge_anno_common$exclude_name)
Summary_phenotypes_merge_anno_common <- Summary_phenotypes_merge_anno_common[!Summary_phenotypes_merge_anno_common$exclude_name %in% c("", "NULL", NA), ]  #695
Summary_phenotypes_merge_anno_common$phecode <- paste0("X", Summary_phenotypes_merge_anno_common$phecode)
# #Discovery cohort
# load("Summary_phenotypes109_merge_anno_v1.2b_matched_ratio2_p123_90K_phenotypes109.RData")
# #n=696 PheCodes having results
# Summary_phenotypes_merge_anno_common <- Summary_phenotypes_merge_anno[Summary_phenotypes_merge_anno$phecode %in% phecode_mapping_unique_commondisease$phecode, ] #696
# head(Summary_phenotypes_merge_anno_common)
# 

summary_final_SKATO_beta55 <- merge(Summary_phenotypes_merge_anno_common, summary_final_SKATO_beta55, by = "phecode")
summary_final_SKATO_beta55_90K <- summary_final_SKATO_beta55
dim(summary_final_SKATO_beta55_90K) #696*12
summary_final_SKATO_beta55_90K$pvalue
summary_final_SKATO_beta55_90K$qvalue <- p.adjust(summary_final_SKATO_beta55_90K$pvalue, method="BH")
summary_final_SKATO_beta55_90K$direction <- ifelse(summary_final_SKATO_beta55_90K$delta_change > 0, "increase", "decrease")
save(summary_final_SKATO_beta55_90K, file = "summmary_final_SKATO_beta55_90K_qvalue.RData", version = 2)


summary_final_SKATO_beta55 <- merge(Summary_phenotypes_merge_anno_common, summary_final_SKATO_beta55, by = "phecode")
summary_final_SKATO_beta55_85K <- summary_final_SKATO_beta55
dim(summary_final_SKATO_beta55_85K) #580*12
summary_final_SKATO_beta55_85K$pvalue
summary_final_SKATO_beta55_85K$qvalue <- p.adjust(summary_final_SKATO_beta55_85K$pvalue, method="BH")
summary_final_SKATO_beta55_85K$direction <- ifelse(summary_final_SKATO_beta55_85K$delta_change > 0, "increase", "decrease")
save(summary_final_SKATO_beta55_85K, file = "summmary_final_SKATO_beta55_85K_qvalue.RData", version = 2)
str(summary_final_SKATO_beta55_85K)

tiff("Summary_phenotypes_merge_anno_SKATO_common_top_phencode_109_ratio2_90K_ALLSIXGENES_linearweighted_beta55_commonfromeach1.tiff", units="in", width=24, height=16, res=600)
ggplot(summary_final_SKATO_beta55_90K, aes(x = as.factor(phecode), y = -log10(pvalue))) + 
  #geom_pointrange(aes(col = PHENOTYPE, shape = as.factor(cohort)), position=position_dodge(width=0.8)) + 
  #geom_point(aes(shape = as.factor(cohort), color = PHENOTYPE, y = estimate, size = -log10(pvalue)), stat="identity", position = position_dodge(0.8))  +
  geom_point(aes(shape = as.factor(direction), color = exclude_name, y = -log10(pvalue), size = abs(delta_change)), stat="identity", position = position_dodge(0.8))  +
  #ylim(0, 0) +
  #scale_y_log10(limits = c(1,30)) + 
  #scale_shape_manual(values=c("\u25BA","\u25C4","\u25BC","\u25B2"))
  scale_shape_manual(values=c("\u25BC","\u25B2")) +
  #geom_hline(aes(yintercept = -log10(0.05)), linetype='dashed', col = 'black', width = 1) +
  geom_hline(aes(yintercept = -log10(0.0065619)), linetype='dashed', col = 'red', width = 1) +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme_bw() +
  geom_text(data = subset(summary_final_SKATO_beta55_90K, pvalue < 0.0065619), 
            aes(phecode, -log10(pvalue), label = phecode_str), position=position_jitter(width=0,height=0), hjust = 0.5, vjust = 0, size = 4, angle = 0) +
  #facet_wrap(. ~ Header.x, scales="free") + 
  #theme(strip.text.x = element_text(size = 6, colour = "black")) +
  theme(text = element_text(size = 14),
        axis.text.x=element_text(angle = 45, vjust = 1, hjust = 1)) + 
  guides(color = guide_legend(title="Disease Category", override.aes = list(size = 5))) +
  guides(shape = guide_legend(title="Cohort", override.aes = list(size = 5))) +
  guides(size = guide_legend(title="Percent Relative Effect")) +
  ylab("-log10(pvalue)") +
  xlab("phecode")
#ggtitle("Summary_phenotypes_merge_anno_SKATO_common_90K85K")
dev.off()

table(summary_final_SKATO_beta55_85K$direction)
# decrease increase 
# 332      248
table(summary_final_SKATO_beta55_90K$direction)
# decrease increase 
# 487      208

tiff("Summary_phenotypes_merge_anno_SKATO_common_top_phencode_109_ratio2_85K_ALLSIXGENES_linearweighted_beta55_commonfromeach1.tiff", units="in", width=24, height=16, res=300)
ggplot(summary_final_SKATO_beta55_85K, aes(x = as.factor(phecode), y = -log10(pvalue))) + 
  #geom_pointrange(aes(col = PHENOTYPE, shape = as.factor(cohort)), position=position_dodge(width=0.8)) + 
  #geom_point(aes(shape = as.factor(cohort), color = PHENOTYPE, y = estimate, size = -log10(pvalue)), stat="identity", position = position_dodge(0.8))  +
  geom_point(aes(shape = as.factor(direction), color = exclude_name, y = -log10(pvalue), size = abs(delta_change)), stat="identity", position = position_dodge(0.8))  +
  #ylim(0, 0) +
  #scale_y_log10(limits = c(1,5)) + 
  #scale_shape_manual(values=c("\u25BA","\u25C4","\u25BC","\u25B2"))
  scale_shape_manual(values=c("\u25BC","\u25B2")) +
  geom_hline(aes(yintercept = -log10(0.1)), linetype='dashed', col = 'black', width = 1) +
  #geom_hline(aes(yintercept = -log10(0.1)), linetype='dashed', col = 'red', width = 1) +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme_bw() +
  geom_text(data = subset(summary_final_SKATO_beta55_85K, pvalue < 0.1), 
            aes(phecode, -log10(pvalue), label = phecode_str), position=position_jitter(width=0,height=0), hjust = 0.5, vjust = 0, size = 4, angle = 0) +
  #facet_wrap(. ~ Header.x, scales="free") + 
  #theme(strip.text.x = element_text(size = 6, colour = "black")) +
  theme(text = element_text(size = 14),
        axis.text.x=element_text(angle = 45, vjust = 1, hjust = 1)) + 
  guides(color = guide_legend(title="Disease Category", override.aes = list(size = 5))) +
  guides(shape = guide_legend(title="Cohort", override.aes = list(size = 5))) +
  guides(size = guide_legend(title="Percent Relative Effect")) +
  ylab("-log10(pvalue)") +
  xlab("phecode")
#ggtitle("Summary_phenotypes_merge_anno_SKATO_common_90K85K")
dev.off()


write.table(summary_final_SKATO_beta55_85K, file = "summary_final_SKATO_beta55_85K.txt", sep = "|", quote = F, row.name = F, col.name = T)
write.table(summary_final_SKATO_beta55_90K, file = "summary_final_SKATO_beta55_90K.txt", sep = "|", quote = F, row.name = F, col.name = T)
save.image("05202022_SKATO_final.RData", version = 2)

setwd("D:/Geisinger/X19981_backup_10072020/Desktop/MS_project/")
load("summary_final_12032021_matched_ratio2_90K_age_ALLSIXGENES_ALLPHENCODE_linearweighted_beta55_ALL.RData")
load("summary_final_12032021_matched_ratio2_85K_age_ALLSIXGENES_ALLPHENCODE_linearweighted_beta55_ALL_fixedweights.RData")
#####################################################################################################################
save(summary_final_SKATO, file = "summary_final_SKATO_12032021_matched_ratio2.RData", version = 2)
save(summary_final_SKAT, file = "summary_final_SKAT_12032021_matched_ratio2.RData", version = 2)
save(summary_final_BURDEN, file = "summary_final_BURDEN_12032021_matched_ratio2.RData", version = 2)
#transfer the above files to local laptop

save(summary_final_SKATO, file = "summary_final_SKATO_12032021_matched_ratio2_p13.RData", version = 2)
save(summary_final_SKAT, file = "summary_final_SKAT_12032021_matched_ratio2_p13.RData", version = 2)
save(summary_final_BURDEN, file = "summary_final_BURDEN_12032021_matched_ratio2_p13.RData", version = 2)

save(summary_final_SKATO, file = "summary_final_SKATO_12032021_matched_ratio2_p13_90K.RData", version = 2)
save(summary_final_SKAT, file = "summary_final_SKAT_12032021_matched_ratio2_p13_90K.RData", version = 2)
save(summary_final_BURDEN, file = "summary_final_BURDEN_12032021_matched_ratio2_p13_90K.RData", version = 2)

save(summary_final_SKATO, file = "summary_final_SKATO_12032021_matched_ratio2_p13_90K_Kin.RData", version = 2)
save(summary_final_SKAT, file = "summary_final_SKAT_12032021_matched_ratio2_p13_90K_Kin.RData", version = 2)
save(summary_final_BURDEN, file = "summary_final_BURDEN_12032021_matched_ratio2_p13_90K_Kin.RData", version = 2)

save(summary_final_SKATO, file = "summary_final_SKATO_12032021_matched_ratio2_p13_85K_Kin.RData", version = 2)
save(summary_final_SKAT, file = "summary_final_SKAT_12032021_matched_ratio2_p13_85K_Kin.RData", version = 2)
save(summary_final_BURDEN, file = "summary_final_BURDEN_12032021_matched_ratio2_p13_85K_Kin.RData", version = 2)

save(summary_final_SKATO, file = "summary_final_SKATO_12032021_matched_ratio2_p123_90K_Kin.RData", version = 2)
save(summary_final_SKAT, file = "summary_final_SKAT_12032021_matched_ratio2_p123_90K_Kin.RData", version = 2)
save(summary_final_BURDEN, file = "summary_final_BURDEN_12032021_matched_ratio2_p123_90K_Kin.RData", version = 2)

save(summary_final_SKATO, file = "summary_final_SKATO_12032021_matched_ratio2_p123_85K_Kin.RData", version = 2)
save(summary_final_SKAT, file = "summary_final_SKAT_12032021_matched_ratio2_p123_85K_Kin.RData", version = 2)
save(summary_final_BURDEN, file = "summary_final_BURDEN_12032021_matched_ratio2_p123_85K_Kin.RData", version = 2)

save(summary_final_SKATO, file = "summary_final_SKATO_12032021_matched_ratio2_p123_85K_Kin_ICD109.RData", version = 2)
save(summary_final_SKAT, file = "summary_final_SKAT_12032021_matched_ratio2_p123_85K_Kin_ICD109.RData", version = 2)
save(summary_final_BURDEN, file = "summary_final_BURDEN_12032021_matched_ratio2_p123_85K_Kin_ICD109.RData", version = 2)

save(summary_final_SKATO, file = "summary_final_SKATO_12032021_matched_ratio2_p123_90K_Kin_ICD109.RData", version = 2)
save(summary_final_SKAT, file = "summary_final_SKAT_12032021_matched_ratio2_p123_90K_Kin_ICD109.RData", version = 2)
save(summary_final_BURDEN, file = "summary_final_BURDEN_12032021_matched_ratio2_p123_90K_Kin_ICD109.RData", version = 2)

save(summary_final_SKATO, file = "summary_final_SKATO_12032021_matched_ratio2_p123_85K_Kin_ICD109_AGE.RData", version = 2)
save(summary_final_SKAT, file = "summary_final_SKAT_12032021_matched_ratio2_p123_85K_Kin_ICD109_AGE.RData", version = 2)
save(summary_final_BURDEN, file = "summary_final_BURDEN_12032021_matched_ratio2_p123_85K_Kin_ICD109_AGE.RData", version = 2)

save(summary_final_SKATO, file = "summary_final_SKATO_12032021_matched_ratio2_p123_85K_Kin_ICD109_AGE60.RData", version = 2)
save(summary_final_SKAT, file = "summary_final_SKAT_12032021_matched_ratio2_p123_85K_Kin_ICD109_AGE60.RData", version = 2)
save(summary_final_BURDEN, file = "summary_final_BURDEN_12032021_matched_ratio2_p123_85K_Kin_ICD109_AGE60.RData", version = 2)
save(summary_final_SKATO, file = "summary_final_SKATO_12032021_matched_ratio2_p123_90K_Kin_ICD109_AGE60.RData", version = 2)
save(summary_final_SKAT, file = "summary_final_SKAT_12032021_matched_ratio2_p123_90K_Kin_ICD109_AGE60.RData", version = 2)
save(summary_final_BURDEN, file = "summary_final_BURDEN_12032021_matched_ratio2_p123_90K_Kin_ICD109_AGE60.RData", version = 2)

save(summary_final_SKATO, file = "summary_final_SKATO_12032021_matched_ratio2_p123_90K_Kin_ICD109_AGE.RData", version = 2)
save(summary_final_SKAT, file = "summary_final_SKAT_12032021_matched_ratio2_p123_90K_Kin_ICD109_AGE.RData", version = 2)
save(summary_final_BURDEN, file = "summary_final_BURDEN_12032021_matched_ratio2_p123_90K_Kin_ICD109_AGE.RData", version = 2)

load("12092021_SKAT_inputfile_p123_85K_updated_phenotypes109.RData")

#RACE specific inputfile generation
dat_RECUR_casecontrol_select_combined_85K_WHITE <- dat_RECUR_casecontrol_select_combined_85K[dat_RECUR_casecontrol_select_combined_85K$PT_RACE == "White", ]
dat <- cSplit(dat, "FID", "_")
dat$FID_2

dat_WHITE <- dat[dat$FID_2 %in% dat_RECUR_casecontrol_select_combined_85K_WHITE$PT_ID, ]
colnames(dat_WHITE)
geno_WHITE <- as.matrix(dat_WHITE[,c(6:73)]) #Z
dim(geno_WHITE) 

phenotypes_ratio2_WHITE <- phenotypes_ratio2[phenotypes_ratio2$id %in% dat_RECUR_casecontrol_select_combined_85K_WHITE$id, ]
phenotypes_ratio2_WHITE <- data.frame(phenotypes_ratio2_WHITE)
save(dat_WHITE, geno_WHITE, phecode, phenotypes_ratio2_WHITE, dat_RECUR_casecontrol_select_combined_85K_WHITE, file = "12092021_SKAT_inputfile_p123_85K_updated_phenotypes109_WHITE.RData", version = 2)

load("12092021_SKAT_inputfile_p123_90K_updated_phenotypes109_WHITE.RData")
load("12092021_SKAT_inputfile_p123_85K_updated_phenotypes109_WHITE.RData")

SKATO <- list()
SKAT <- list()
BURDEN <- list()
summary_final_SKATO <- NULL
plot_list = list()
phen <- NULL
dim(geno_WHITE)

for (i in phecode) { 
  print(i)
  phen = phenotypes_ratio2_WHITE[, i] #change from phenotypes_ratio5
  print(table(phen))
  if(table(phen) == 4520) next # change from 9249 to 4713 for discovery p13; to 4536 for replication; 7834/342 for discovery p123 WHITE/WHITE; 4520/500 for replicaiton p123 WHITE/BLACK
  #obj <- SKAT_NULL_emmaX(phen~dat_RECUR_casecontrol_select_combined_90K_WHITE$PT_SEX, K=K)
  obj <- SKAT_Null_Model(phen~dat_RECUR_casecontrol_select_combined_85K_WHITE$PT_SEX, out_type="D")
  
  tryCatch({
    SKATO[i] <- SKAT(geno_WHITE, obj, method="SKATO")$p.value
    SKAT[i] <- SKAT(geno_WHITE, obj, method="SKAT")$p.value
    BURDEN[i] <- SKAT(geno_WHITE, obj, method="Burden")$p.value
  }, error=function(e){})
}


summary_final_SKATO <- data.frame(do.call(rbind,SKATO))
summary_final_SKAT <- data.frame(do.call(rbind,SKAT))
summary_final_BURDEN <- data.frame(do.call(rbind,BURDEN))

save(summary_final_SKATO, file = "summary_final_SKATO_12032021_matched_ratio2_p123_85K_WHITE.RData", version = 2)
save(summary_final_SKAT, file = "summary_final_SKAT_12032021_matched_ratio2_p123_85K_WHITE.RData", version = 2)
save(summary_final_BURDEN, file = "summary_final_BURDEN_12032021_matched_ratio2_p123_85K_WHITE.RData", version = 2)


X <- colnames(geno_WHITE)
geno_WHITE_GENE <- geno_WHITE[, colnames(geno_WHITE) %in% X[grepl("^X23.*", X)]] #COL4A1/2
dim(geno_WHITE_GENE)

phecode_select <- c("X425.11", "X425.1", "X425", "X433", "X433.1", "X433.2", "X433.3", "X286.8", "X286.81")

SKATO <- list()
SKAT <- list()
BURDEN <- list()
summary_final_SKATO <- NULL
plot_list = list()
phen <- NULL
dim(geno_WHITE_GENE)

for (i in phecode_select) { 
  print(i)
  phen = phenotypes_ratio2_WHITE[, i] #change from phenotypes_ratio5
  print(table(phen))
  if(table(phen) == 4520) next # change from 9249 to 4713 for discovery p13; to 4536 for replication; 7834/342 for discovery p123 WHITE/WHITE; 4520/500 for replicaiton p123 WHITE/BLACK
  #obj <- SKAT_NULL_emmaX(phen~dat_RECUR_casecontrol_select_combined_90K_WHITE$PT_SEX, K=K)
  obj <- SKAT_Null_Model(phen~dat_RECUR_casecontrol_select_combined_90K_WHITE$PT_SEX, out_type="D")
  
  tryCatch({
    SKATO[i] <- SKAT(geno_WHITE_GENE, obj, method="SKATO")$p.value
    SKAT[i] <- SKAT(geno_WHITE_GENE, obj, method="SKAT")$p.value
    BURDEN[i] <- SKAT(geno_WHITE_GENE, obj, method="Burden")$p.value
  }, error=function(e){})
}

summary_final_SKATO <- data.frame(do.call(rbind,SKATO))
summary_final_SKAT <- data.frame(do.call(rbind,SKAT))
summary_final_BURDEN <- data.frame(do.call(rbind,BURDEN))

save(summary_final_SKATO, file = "summary_final_SKATO_12032021_matched_ratio2_p123_85K_WHITE.RData", version = 2)
save(summary_final_SKAT, file = "summary_final_SKAT_12032021_matched_ratio2_p123_85K_WHITE.RData", version = 2)
save(summary_final_BURDEN, file = "summary_final_BURDEN_12032021_matched_ratio2_p123_85K_WHITE.RData", version = 2)

load("summary_final_SKATO_12032021_matched_ratio2_p123_85K_WHITE.RData")
load("summary_final_SKATO_12032021_matched_ratio2_p123_90K_WHITE.RData")

commondisease <- read.table(file = "Common_disease.txt", header = TRUE, stringsAsFactors=FALSE)
str(summary_final_SKATO)

#for v1.2
phecode_mapping <- read.csv("Phecode_map_v1_2_icd10cm_beta.csv", header = T, stringsAsFactors = F)
phecode_mapping_unique <- phecode_mapping %>%
  dplyr::select(phecode, phecode_str, exclude_range, exclude_name) %>%
  distinct(phecode, .keep_all = T)  #1755*4
str(phecode_mapping_unique)
phecode_mapping_unique$phecodeX <- as.character(paste0("X", phecode_mapping_unique$phecode))

unique(phecode_mapping_unique$exclude_name)
# [1] "injuries & poisonings"   "neoplasms"               "circulatory system"      "digestive"              
# [5] "genitourinary"           "symptoms"                "sense organs"            "mental disorders"       
# [9] "neurological"            "endocrine/metabolic"     "hematopoietic"           "musculoskeletal"        
# [13] "respiratory"             ""                        "infectious diseases"     "pregnancy complications"
# [17] "NULL"                    "dermatologic"            "congenital anomalies"  



################################################################################
load("summary_final_SKATO_12032021_matched_ratio2.RData")

load("summary_final_SKATO_12032021_matched_ratio2_p13.RData")
load("summary_final_SKAT_12032021_matched_ratio2_p13.RData")
load("summary_final_BURDEN_12032021_matched_ratio2_p13.RData")

load("summary_final_SKATO_12032021_matched_ratio2_p13_90K.RData")
load("summary_final_SKAT_12032021_matched_ratio2_p13_90K.RData")
load("summary_final_BURDEN_12032021_matched_ratio2_p13_90K.RData")

load("summary_final_SKATO_12032021_matched_ratio2_p13_90K_Kin.RData")
load("summary_final_SKAT_12032021_matched_ratio2_p13_90K_Kin.RData")
load("summary_final_BURDEN_12032021_matched_ratio2_p13_90K_Kin.RData")

load("summary_final_SKATO_12032021_matched_ratio2_p123_85K_Kin.RData")
load("summary_final_SKAT_12032021_matched_ratio2_p123_85K_Kin.RData")
load("summary_final_BURDEN_12032021_matched_ratio2_p123_85K_Kin.RData")

load("summary_final_SKATO_12032021_matched_ratio2_p123_90K_Kin.RData")
load("summary_final_SKAT_12032021_matched_ratio2_p123_90K_Kin.RData")
load("summary_final_BURDEN_12032021_matched_ratio2_p123_90K_Kin.RData")

load("summary_final_SKATO_12032021_matched_ratio2_p123_85K_Kin_ICD109.RData")
load("summary_final_SKAT_12032021_matched_ratio2_p123_85K_Kin_ICD109.RData")
load("summary_final_BURDEN_12032021_matched_ratio2_p123_85K_Kin_ICD109.RData")

load("summary_final_SKATO_12032021_matched_ratio2_p123_90K_Kin_ICD109.RData")
load("summary_final_SKAT_12032021_matched_ratio2_p123_90K_Kin_ICD109.RData")
load("summary_final_BURDEN_12032021_matched_ratio2_p123_90K_Kin_ICD109.RData")

load("summary_final_SKATO_12032021_matched_ratio2_p123_85K_Kin_ICD109_AGE.RData")
load("summary_final_SKAT_12032021_matched_ratio2_p123_85K_Kin_ICD109_AGE.RData")
load("summary_final_BURDEN_12032021_matched_ratio2_p123_85K_Kin_ICD109_AGE.RData")
load("summary_final_SKATO_12032021_matched_ratio2_p123_85K_Kin_ICD109_AGE60.RData")
load("summary_final_SKAT_12032021_matched_ratio2_p123_85K_Kin_ICD109_AGE60.RData")
load("summary_final_BURDEN_12032021_matched_ratio2_p123_85K_Kin_ICD109_AGE60.RData")


load("summary_final_SKATO_12032021_matched_ratio2_p123_90K_Kin_ICD109_AGE.RData")
load("summary_final_SKAT_12032021_matched_ratio2_p123_90K_Kin_ICD109_AGE.RData")
load("summary_final_BURDEN_12032021_matched_ratio2_p123_90K_Kin_ICD109_AGE.RData")
load("summary_final_SKATO_12032021_matched_ratio2_p123_90K_Kin_ICD109_AGE60.RData")
load("summary_final_SKAT_12032021_matched_ratio2_p123_90K_Kin_ICD109_AGE60.RData")
load("summary_final_BURDEN_12032021_matched_ratio2_p123_90K_Kin_ICD109_AGE60.RData")

summary_final_SKATO$phecode <- rownames(summary_final_SKATO)
#summary_final_SKATO$phecode <- as.numeric(gsub("X", "", summary_final_SKATO$phecode))
colnames(summary_final_SKATO)[1] <- "SKATO_pvalue"
str(summary_final_SKATO)

load("Summary_phenotypes_merge_anno_commonrare_delta_p123_85K.RData")
load("Summary_phenotypes_merge_anno_commonrare_delta_p123_90K.RData")

load("Summary_phenotypes_merge_anno_commonrare_delta_p13.RData") #load the final common and rare diseases for pipeline 13
load("Summary_phenotypes109_merge_anno_commonrare_delta_p123_85K_phenotypes109_age.RData")
load("Summary_phenotypes109_merge_anno_commonrare_delta_p123_90K_phenotypes109_age.RData")
load("Summary_phenotypes109_merge_anno_commonrare_delta_p123_85K_phenotypes109.RData")
load("Summary_phenotypes109_merge_anno_commonrare_delta_p123_90K_phenotypes109.RData")
load("Summary_phenotypes109_merge_anno_commonrare_delta_p123_90K_phenotypes109_age60.RData")
load("Summary_phenotypes109_merge_anno_commonrare_delta_p123_85K_phenotypes109_age60.RData")

head(Summary_phenotypes_merge_anno_common) 
head(Summary_phenotypes_merge_anno_rare) 
str(Summary_phenotypes_merge_anno_common)
str(Summary_phenotypes_merge_anno_rare)

#optional
Summary_phenotypes_merge_anno_common$phecode <- paste0("X", Summary_phenotypes_merge_anno_common$phecode)
#Summary_phenotypes_merge_anno_common$phecode <- gsub("X", "", Summary_phenotypes_merge_anno_common$phecode)
Summary_phenotypes_merge_anno_rare$phecode <- paste0("X", Summary_phenotypes_merge_anno_rare$phecode)
#Summary_phenotypes_merge_anno_rare$phecode <- gsub("X", "", Summary_phenotypes_merge_anno_rare$phecode)


summary_final_SKATO_common <- merge(Summary_phenotypes_merge_anno_common, summary_final_SKATO, by = "phecode")
summary_final_SKATO_increase_common <- summary_final_SKATO_common[summary_final_SKATO_common$delta_change > 0, ]
summary_final_SKATO_decrease_common <- summary_final_SKATO_common[summary_final_SKATO_common$delta_change <= 0, ]

summary_final_SKATO_rare <- merge(Summary_phenotypes_merge_anno_rare, summary_final_SKATO, by = "phecode")
summary_final_SKATO_increase_rare <- summary_final_SKATO_rare[summary_final_SKATO_rare$delta_change > 0, ]
summary_final_SKATO_decrease_rare <- summary_final_SKATO_rare[summary_final_SKATO_rare$delta_change <= 0, ]

summary_final_SKATO_increase <- rbind(summary_final_SKATO_increase_common, summary_final_SKATO_increase_rare)
summary_final_SKATO_decrease <- rbind(summary_final_SKATO_decrease_common, summary_final_SKATO_decrease_rare)

save(summary_final_SKATO_increase, summary_final_SKATO_decrease, file = "summary_final_SKATO_12032021_matched_ratio2_increasedecrease_p123_90K_phenotypes109.RData", version = 2)
save(summary_final_SKAT_increase, summary_final_SKAT_decrease, file = "summary_final_SKAT_12032021_matched_ratio2_increasedecrease_p123_90K_phenotypes109.RData", version = 2)
save(summary_final_BURDEN_increase, summary_final_BURDEN_decrease, file = "summary_final_BURDEN_12032021_matched_ratio2_increasedecrease_p123_90K_phenotypes109.RData", version = 2)

save(summary_final_SKATO_increase, summary_final_SKATO_decrease, file = "summary_final_SKATO_12032021_matched_ratio2_increasedecrease_p123_85K_phenotypes109.RData", version = 2)
save(summary_final_SKAT_increase, summary_final_SKAT_decrease, file = "summary_final_SKAT_12032021_matched_ratio2_increasedecrease_p123_85K_phenotypes109.RData", version = 2)
save(summary_final_BURDEN_increase, summary_final_BURDEN_decrease, file = "summary_final_BURDEN_12032021_matched_ratio2_increasedecrease_p123_85K_phenotypes109.RData", version = 2)

save(summary_final_SKATO_increase, summary_final_SKATO_decrease, file = "summary_final_SKATO_12032021_matched_ratio2_increasedecrease_p123_85K_phenotypes109_age.RData", version = 2)
save(summary_final_SKAT_increase, summary_final_SKAT_decrease, file = "summary_final_SKAT_12032021_matched_ratio2_increasedecrease_p123_85K_phenotypes109_age.RData", version = 2)
save(summary_final_BURDEN_increase, summary_final_BURDEN_decrease, file = "summary_final_BURDEN_12032021_matched_ratio2_increasedecrease_p123_85K_phenotypes109_age.RData", version = 2)
save(summary_final_SKATO_increase, summary_final_SKATO_decrease, file = "summary_final_SKATO_12032021_matched_ratio2_increasedecrease_p123_85K_phenotypes109_age60.RData", version = 2)
save(summary_final_SKAT_increase, summary_final_SKAT_decrease, file = "summary_final_SKAT_12032021_matched_ratio2_increasedecrease_p123_85K_phenotypes109_age60.RData", version = 2)
save(summary_final_BURDEN_increase, summary_final_BURDEN_decrease, file = "summary_final_BURDEN_12032021_matched_ratio2_increasedecrease_p123_85K_phenotypes109_age60.RData", version = 2)

save(summary_final_SKATO_increase, summary_final_SKATO_decrease, file = "summary_final_SKATO_12032021_matched_ratio2_increasedecrease_p123_90K_phenotypes109_age.RData", version = 2)
save(summary_final_SKAT_increase, summary_final_SKAT_decrease, file = "summary_final_SKAT_12032021_matched_ratio2_increasedecrease_p123_90K_phenotypes109_age.RData", version = 2)
save(summary_final_BURDEN_increase, summary_final_BURDEN_decrease, file = "summary_final_BURDEN_12032021_matched_ratio2_increasedecrease_p123_90K_phenotypes109_age.RData", version = 2)
save(summary_final_SKATO_increase, summary_final_SKATO_decrease, file = "summary_final_SKATO_12032021_matched_ratio2_increasedecrease_p123_90K_phenotypes109_age60.RData", version = 2)
save(summary_final_SKAT_increase, summary_final_SKAT_decrease, file = "summary_final_SKAT_12032021_matched_ratio2_increasedecrease_p123_90K_phenotypes109_age60.RData", version = 2)
save(summary_final_BURDEN_increase, summary_final_BURDEN_decrease, file = "summary_final_BURDEN_12032021_matched_ratio2_increasedecrease_p123_90K_phenotypes109_age60.RData", version = 2)

summary_final_SKATO_increasedecrease <- rbind(summary_final_SKATO_increase, summary_final_SKATO_decrease)
summary_final_SKAT_increasedecrease <- rbind(summary_final_SKAT_increase, summary_final_SKAT_decrease)
summary_final_BURDEN_increasedecrease <- rbind(summary_final_BURDEN_increase, summary_final_BURDEN_decrease)

write.table(summary_final_SKATO_increasedecrease, file = "summary_final_SKATO_increasedecrease_90K_phenotypes109.txt", sep = "|", quote =F, row.name = F, col.name = T)
write.table(summary_final_SKAT_increasedecrease, file = "summary_final_SKAT_increasedecrease_90K_phenotypes109.txt", sep = "|", quote =F, row.name = F, col.name = T)
write.table(summary_final_BURDEN_increasedecrease, file = "summary_final_BURDEN_increasedecrease_90K_phenotypes109.txt", sep = "|", quote =F, row.name = F, col.name = T)

write.table(summary_final_SKATO_increasedecrease, file = "summary_final_SKATO_increasedecrease_85K_phenotypes109.txt", sep = "|", quote =F, row.name = F, col.name = T)
write.table(summary_final_SKAT_increasedecrease, file = "summary_final_SKAT_increasedecrease_85K_phenotypes109.txt", sep = "|", quote =F, row.name = F, col.name = T)
write.table(summary_final_BURDEN_increasedecrease, file = "summary_final_BURDEN_increasedecrease_85K_phenotypes109.txt", sep = "|", quote =F, row.name = F, col.name = T)

write.table(summary_final_SKATO_increasedecrease, file = "summary_final_SKATO_increasedecrease_85K_phenotypes109_age.txt", sep = "|", quote =F, row.name = F, col.name = T)
write.table(summary_final_SKAT_increasedecrease, file = "summary_final_SKAT_increasedecrease_85K_phenotypes109_age.txt", sep = "|", quote =F, row.name = F, col.name = T)
write.table(summary_final_BURDEN_increasedecrease, file = "summary_final_BURDEN_increasedecrease_85K_phenotypes109_age.txt", sep = "|", quote =F, row.name = F, col.name = T)
write.table(summary_final_SKATO_increasedecrease, file = "summary_final_SKATO_increasedecrease_85K_phenotypes109_age60.txt", sep = "|", quote =F, row.name = F, col.name = T)
write.table(summary_final_SKAT_increasedecrease, file = "summary_final_SKAT_increasedecrease_85K_phenotypes109_age60.txt", sep = "|", quote =F, row.name = F, col.name = T)
write.table(summary_final_BURDEN_increasedecrease, file = "summary_final_BURDEN_increasedecrease_85K_phenotypes109_age60.txt", sep = "|", quote =F, row.name = F, col.name = T)

write.table(summary_final_SKATO_increasedecrease, file = "summary_final_SKATO_increasedecrease_90K_phenotypes109_age.txt", sep = "|", quote =F, row.name = F, col.name = T)
write.table(summary_final_SKAT_increasedecrease, file = "summary_final_SKAT_increasedecrease_90K_phenotypes109_age.txt", sep = "|", quote =F, row.name = F, col.name = T)
write.table(summary_final_BURDEN_increasedecrease, file = "summary_final_BURDEN_increasedecrease_90K_phenotypes109_age.txt", sep = "|", quote =F, row.name = F, col.name = T)
write.table(summary_final_SKATO_increasedecrease, file = "summary_final_SKATO_increasedecrease_90K_phenotypes109_age60.txt", sep = "|", quote =F, row.name = F, col.name = T)
write.table(summary_final_SKAT_increasedecrease, file = "summary_final_SKAT_increasedecrease_90K_phenotypes109_age60.txt", sep = "|", quote =F, row.name = F, col.name = T)
write.table(summary_final_BURDEN_increasedecrease, file = "summary_final_BURDEN_increasedecrease_90K_phenotypes109_age60.txt", sep = "|", quote =F, row.name = F, col.name = T)

library(dplyr)
Summary_phenotypes_merge_anno_BURDEN_common_top <- Summary_phenotypes_merge_anno_BURDEN %>%
  dplyr::filter(delta_change > 0) %>%
  dplyr::filter(frequency == "common")

Summary_phenotypes_merge_anno_BURDEN_rare_top <- Summary_phenotypes_merge_anno_BURDEN %>%
  dplyr::filter(delta_change > 0) %>%
  dplyr::filter(frequency == "rare")

Summary_phenotypes_merge_anno_BURDEN_common_bottom <- Summary_phenotypes_merge_anno_BURDEN %>%
  dplyr::filter(delta_change <= 0) %>%
  dplyr::filter(frequency == "common")

Summary_phenotypes_merge_anno_BURDEN_rare_bottom <- Summary_phenotypes_merge_anno_BURDEN %>%
  dplyr::filter(delta_change <= 0) %>%
  dplyr::filter(frequency == "rare")

save(Summary_phenotypes_merge_anno_BURDEN_common_top, Summary_phenotypes_merge_anno_BURDEN_rare_top, Summary_phenotypes_merge_anno_BURDEN_common_bottom, Summary_phenotypes_merge_anno_BURDEN_rare_bottom, file = "Summary_phenotypes_merge_anno_BURDEN_commonrare_topbottom_concur_ratio5.RData", version = 2)
save(Summary_phenotypes_merge_anno_SKAT_common_top, Summary_phenotypes_merge_anno_SKAT_rare_top, Summary_phenotypes_merge_anno_SKAT_common_bottom, Summary_phenotypes_merge_anno_SKAT_rare_bottom, file = "Summary_phenotypes_merge_anno_SKAT_commonrare_topbottom_concur_ratio2.RData", version = 2)

save(Summary_phenotypes_merge_anno_BURDEN_common_top, Summary_phenotypes_merge_anno_BURDEN_rare_top, Summary_phenotypes_merge_anno_BURDEN_common_bottom, Summary_phenotypes_merge_anno_BURDEN_rare_bottom, file = "Summary_phenotypes_merge_anno_BURDEN_commonrare_topbottom_union_ratio5.RData", version = 2)
save(Summary_phenotypes_merge_anno_BURDEN_common_top, Summary_phenotypes_merge_anno_BURDEN_rare_top, Summary_phenotypes_merge_anno_BURDEN_common_bottom, Summary_phenotypes_merge_anno_BURDEN_rare_bottom, file = "Summary_phenotypes_merge_anno_BURDEN_commonrare_topbottom_union_ratio2.RData", version = 2)

#top
load("Summary_phenotypes_merge_anno_SKAT_commonrare_topbottom_concur_ratio5.RData")
colnames(Summary_phenotypes_merge_anno_SKAT_common_top)
# [1] "phecode"               "group.x"               "count.x"               "group.y"              
# [5] "count.y"               "delta_change"          "delta"                 "phecode_str"          
# [9] "exclude_range"         "exclude_name"          "frequency"             "do.call.rbind..SKAT."
str(Summary_phenotypes_merge_anno_SKAT_common_top)
colnames(Summary_phenotypes_merge_anno_SKAT_common_top)[12] <- "p.value"
class(Summary_phenotypes_merge_anno_SKAT_common_top)
Summary_phenotypes_merge_anno_SKAT_common_top <- Summary_phenotypes_merge_anno_SKAT_common_top[!is.na(Summary_phenotypes_merge_anno_SKAT_common_top$exclude_name), ] %>%
  filter(exclude_name != "NULL") %>%
  arrange(exclude_name)
#Summary_phenotypes_merge_anno_SKAT_common_top <- Summary_phenotypes_merge_anno_SKAT_common_top[sort(as.factor(Summary_phenotypes_merge_anno_SKAT_common_top$exclude_name)), ]
library(ggplot2)
tiff("Summary_phenotypes_merge_anno_SKAT_common_top_concur_ratio5.tiff", units="in", width=24, height=10, res=600)
ggplot(Summary_phenotypes_merge_anno_SKAT_common_top, aes(x = as.factor(phecode), y = -log10(p.value))) + 
  #geom_pointrange(aes(col = PHENOTYPE, shape = as.factor(Direction)), position=position_dodge(width=0.8)) + 
  #geom_point(aes(shape = as.factor(Direction), color = PHENOTYPE, y = estimate, size = -log10(p.value)), stat="identity", position = position_dodge(0.8))  +
  geom_point(aes(color = exclude_name, y = -log10(p.value), size = delta_change), stat="identity", position = position_dodge(0.8))  +
  
  #scale_shape_manual(values=c("\u25BA","\u25C4","\u25BC","\u25B2"))
  #scale_shape_manual(values=c("\u25BC","\u25B2")) +
  #geom_hline(aes(yintercept = -log(0.001)), linetype='dashed', col = 'black', width = 1) +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme_bw() +
  geom_text(data = subset(Summary_phenotypes_merge_anno_SKAT_common_top, p.value < 0.005), 
           aes(phecode, -log10(p.value), label = phecode_str), vjust = -1, size = 4, angle = 30) +
  #facet_wrap(. ~ Header.x, scales="free") + 
  #theme(strip.text.x = element_text(size = 6, colour = "black")) +
  theme(text = element_text(size = 10),
        axis.text.x=element_text(angle = 45, vjust = 0.5, hjust=0.5)) + 
  ylab("-log10(p.value)") +
  xlab("phecode") + 
  ggtitle("Summary_phenotypes_merge_anno_SKAT_common_top_concur_ratio5")
dev.off()


#bottom
load("Summary_phenotypes_merge_anno_SKAT_commonrare_topbottom_concur_ratio5.RData")
colnames(Summary_phenotypes_merge_anno_SKAT_common_bottom)
# [1] "phecode"               "group.x"               "count.x"               "group.y"              
# [5] "count.y"               "delta_change"          "delta"                 "phecode_str"          
# [9] "exclude_range"         "exclude_name"          "frequency"             "do.call.rbind..SKAT."
str(Summary_phenotypes_merge_anno_SKAT_common_bottom)
colnames(Summary_phenotypes_merge_anno_SKAT_common_bottom)[12] <- "p.value"
class(Summary_phenotypes_merge_anno_SKAT_common_bottom)
Summary_phenotypes_merge_anno_SKAT_common_bottom <- Summary_phenotypes_merge_anno_SKAT_common_bottom[!is.na(Summary_phenotypes_merge_anno_SKAT_common_bottom$exclude_name), ] %>%
  filter(exclude_name != "NULL") %>%
  arrange(exclude_name)
#Summary_phenotypes_merge_anno_SKAT_common_bottom <- Summary_phenotypes_merge_anno_SKAT_common_bottom[sort(as.factor(Summary_phenotypes_merge_anno_SKAT_common_bottom$exclude_name)), ]
library(ggplot2)
tiff("Summary_phenotypes_merge_anno_SKAT_common_bottom_concur_ratio5.tiff", units="in", width=24, height=10, res=600)
ggplot(Summary_phenotypes_merge_anno_SKAT_common_bottom, aes(x = as.factor(phecode), y = -log10(p.value))) + 
  #geom_pointrange(aes(col = PHENOTYPE, shape = as.factor(Direction)), position=position_dodge(width=0.8)) + 
  #geom_point(aes(shape = as.factor(Direction), color = PHENOTYPE, y = estimate, size = -log10(p.value)), stat="identity", position = position_dodge(0.8))  +
  geom_point(aes(color = exclude_name, y = -log10(p.value), size = -delta_change), stat="identity", position = position_dodge(0.8))  +
  
  #scale_shape_manual(values=c("\u25BA","\u25C4","\u25BC","\u25B2"))
  #scale_shape_manual(values=c("\u25BC","\u25B2")) +
  #geom_hline(aes(yintercept = -log(0.001)), linetype='dashed', col = 'black', width = 1) +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme_bw() +
  geom_text(data = subset(Summary_phenotypes_merge_anno_SKAT_common_bottom, p.value < 0.005), 
            aes(phecode, -log10(p.value), label = phecode_str), vjust = -1, size = 4, angle = 30) +
  #facet_wrap(. ~ Header.x, scales="free") + 
  #theme(strip.text.x = element_text(size = 6, colour = "black")) +
  theme(text = element_text(size = 10),
        axis.text.x=element_text(angle = 45, vjust = 0.5, hjust=0.5)) + 
  ylab("-log10(p.value)") +
  xlab("phecode") + 
  ggtitle("Summary_phenotypes_merge_anno_SKAT_common_bottom_concur_ratio5")
dev.off()

#combining top and bottom
load("Summary_phenotypes_merge_anno_BURDEN_commonrare_topbottom_union_ratio5.RData")
colnames(Summary_phenotypes_merge_anno_BURDEN_common_top)
# [1] "phecode"               "group.x"               "count.x"               "group.y"              
# [5] "count.y"               "delta_change"          "delta"                 "phecode_str"          
# [9] "exclude_range"         "exclude_name"          "frequency"             "do.call.rbind..BURDEN."
str(Summary_phenotypes_merge_anno_BURDEN_common_top)
colnames(Summary_phenotypes_merge_anno_BURDEN_common_top)[12] <- "p.value"
class(Summary_phenotypes_merge_anno_BURDEN_common_top)
Summary_phenotypes_merge_anno_BURDEN_common_top <- Summary_phenotypes_merge_anno_BURDEN_common_top[!is.na(Summary_phenotypes_merge_anno_BURDEN_common_top$exclude_name), ] %>%
  filter(exclude_name != "NULL" & exclude_name != "") %>%
  arrange(exclude_name) %>%
  mutate(direction = 1)

colnames(Summary_phenotypes_merge_anno_BURDEN_common_bottom)
# [1] "phecode"               "group.x"               "count.x"               "group.y"              
# [5] "count.y"               "delta_change"          "delta"                 "phecode_str"          
# [9] "exclude_range"         "exclude_name"          "frequency"             "do.call.rbind..BURDEN."
str(Summary_phenotypes_merge_anno_BURDEN_common_bottom)
colnames(Summary_phenotypes_merge_anno_BURDEN_common_bottom)[12] <- "p.value"
class(Summary_phenotypes_merge_anno_BURDEN_common_bottom)
Summary_phenotypes_merge_anno_BURDEN_common_bottom <- Summary_phenotypes_merge_anno_BURDEN_common_bottom[!is.na(Summary_phenotypes_merge_anno_BURDEN_common_bottom$exclude_name), ] %>%
  filter(exclude_name != "NULL" & exclude_name != "") %>%
  arrange(exclude_name) %>%
  mutate(direction = -1)

Summary_phenotypes_merge_anno_BURDEN_common_topbottom <- rbind(Summary_phenotypes_merge_anno_BURDEN_common_top, Summary_phenotypes_merge_anno_BURDEN_common_bottom) %>%
  arrange(exclude_name)
Summary_phenotypes_merge_anno_BURDEN_common_bottom
tiff("Summary_phenotypes_merge_anno_BURDEN_common_topbottom_union_ratio5.tiff", units="in", width=24, height=10, res=600)
ggplot(Summary_phenotypes_merge_anno_BURDEN_common_topbottom, aes(x = as.factor(phecode), y = -log10(p.value))) + 
  #geom_pointrange(aes(col = PHENOTYPE, shape = as.factor(Direction)), position=position_dodge(width=0.8)) + 
  #geom_point(aes(shape = as.factor(Direction), color = PHENOTYPE, y = estimate, size = -log10(p.value)), stat="identity", position = position_dodge(0.8))  +
  geom_point(aes(shape = as.factor(direction), color = exclude_name, y = -log10(p.value), size = abs(delta_change)), stat="identity", position = position_dodge(0.8))  +
  
  #scale_shape_manual(values=c("\u25BA","\u25C4","\u25BC","\u25B2"))
  scale_shape_manual(values=c("\u25BC","\u25B2")) +
  geom_hline(aes(yintercept = -log10(0.005)), linetype='dashed', col = 'black', width = 1) +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme_bw() +
  geom_text(data = subset(Summary_phenotypes_merge_anno_BURDEN_common_topbottom, p.value < 0.005), 
            aes(phecode, -log10(p.value), label = phecode_str), vjust = -1, size = 4, angle = 30) +
  #facet_wrap(. ~ Header.x, scales="free") + 
  #theme(strip.text.x = element_text(size = 6, colour = "black")) +
  theme(text = element_text(size = 10),
        axis.text.x=element_text(angle = 45, vjust = 0.5, hjust=0.5)) + 
  ylab("-log10(p.value)") +
  xlab("phecode") + 
  ggtitle("Summary_phenotypes_merge_anno_BURDEN_common_topbottom_union_ratio5")
dev.off()

#rare diseases
#combining top and bottom
load("Summary_phenotypes_merge_anno_SKATO_commonrare_topbottom_union_ratio2.RData")
colnames(Summary_phenotypes_merge_anno_SKATO_rare_top)
# [1] "phecode"               "group.x"               "count.x"               "group.y"              
# [5] "count.y"               "delta_change"          "delta"                 "phecode_str"          
# [9] "exclude_range"         "exclude_name"          "frequency"             "do.call.rbind..SKATO."
str(Summary_phenotypes_merge_anno_SKATO_rare_top)
colnames(Summary_phenotypes_merge_anno_SKATO_rare_top)[12] <- "p.value"
class(Summary_phenotypes_merge_anno_SKATO_rare_top)
Summary_phenotypes_merge_anno_SKATO_rare_top <- Summary_phenotypes_merge_anno_SKATO_rare_top[!is.na(Summary_phenotypes_merge_anno_SKATO_rare_top$exclude_name), ] %>%
  filter(exclude_name != "NULL" & exclude_name != "") %>%
  arrange(exclude_name) %>%
  mutate(direction = 1)

colnames(Summary_phenotypes_merge_anno_SKATO_rare_bottom)
# [1] "phecode"               "group.x"               "count.x"               "group.y"              
# [5] "count.y"               "delta_change"          "delta"                 "phecode_str"          
# [9] "exclude_range"         "exclude_name"          "frequency"             "do.call.rbind..SKATO."
str(Summary_phenotypes_merge_anno_SKATO_rare_bottom)
colnames(Summary_phenotypes_merge_anno_SKATO_rare_bottom)[12] <- "p.value"
class(Summary_phenotypes_merge_anno_SKATO_rare_bottom)
Summary_phenotypes_merge_anno_SKATO_rare_bottom <- Summary_phenotypes_merge_anno_SKATO_rare_bottom[!is.na(Summary_phenotypes_merge_anno_SKATO_rare_bottom$exclude_name), ] %>%
  filter(exclude_name != "NULL" & exclude_name != "") %>%
  arrange(exclude_name) %>%
  mutate(direction = -1)

Summary_phenotypes_merge_anno_SKATO_rare_topbottom <- rbind(Summary_phenotypes_merge_anno_SKATO_rare_top, Summary_phenotypes_merge_anno_SKATO_rare_bottom) %>%
  arrange(exclude_name)
Summary_phenotypes_merge_anno_SKATO_rare_bottom
tiff("Summary_phenotypes_merge_anno_SKATO_rare_topbottom_union_ratio2.tiff", units="in", width=24, height=10, res=600)
ggplot(Summary_phenotypes_merge_anno_SKATO_rare_topbottom, aes(x = as.factor(phecode), y = -log10(p.value))) + 
  #geom_pointrange(aes(col = PHENOTYPE, shape = as.factor(Direction)), position=position_dodge(width=0.8)) + 
  #geom_point(aes(shape = as.factor(Direction), color = PHENOTYPE, y = estimate, size = -log10(p.value)), stat="identity", position = position_dodge(0.8))  +
  geom_point(aes(shape = as.factor(direction), color = exclude_name, y = -log10(p.value), size = abs(delta_change)), stat="identity", position = position_doage(0.8)) +
  #scale_shape_manual(values=c("\u25BA","\u25C4","\u25BC","\u25B2")) +
  
  scale_shape_manual(values=c("\u25BC","\u25B2")) +
  geom_hline(aes(yintercept = -log10(0.005)), linetype='dashed', col = 'black', width = 1) +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme_bw() +
  geom_text(data = subset(Summary_phenotypes_merge_anno_SKATO_rare_topbottom, p.value < 0.005), 
            aes(phecode, -log10(p.value), label = phecode_str), vjust = -1, size = 4, angle = 30) +
  #facet_wrap(. ~ Header.x, scales="free") + 
  #theme(strip.text.x = element_text(size = 6, colour = "black")) +
  theme(text = element_text(size = 10),
        axis.text.x=element_text(angle = 45, vjust = 0.5, hjust=0.5)) + 
  ylab("-log10(p.value)") +
  xlab("phecode") + 
  ggtitle("Summary_phenotypes_merge_anno_SKATO_rare_topbottom_union_ratio2")
dev.off()

#############################################################################################################################################

load("summary_final_BURDEN_12032021_matched_ratio2_increasedecrease_p123_90K_phenotypes109.RData")
load("summary_final_BURDEN_12032021_matched_ratio2_increasedecrease_p123_85K_phenotypes109.RData")
summary_final_BURDEN_increasedecrease <- rbind(summary_final_BURDEN_increase, summary_final_BURDEN_decrease)
head(summary_final_BURDEN_increasedecrease)
common_phecode <- summary_final_BURDEN_increasedecrease[summary_final_BURDEN_increasedecrease$frequency == "common", ]$phecode #782
rare_phecode <- summary_final_BURDEN_increasedecrease[summary_final_BURDEN_increasedecrease$frequency == "rare", ]$phecode #866

load("summary_final_BURDEN_12032021_matched_ratio2_increasedecrease_p123_90K_phenotypes109_age.RData")

Summary_phenotypes_merge_anno_BURDEN_common_top <- summary_final_BURDEN_increase[summary_final_BURDEN_increase$phecode %in% common_phecode, ]

Summary_phenotypes_merge_anno_BURDEN_rare_top <- summary_final_BURDEN_increase[summary_final_BURDEN_increase$phecode %in% rare_phecode, ]

Summary_phenotypes_merge_anno_BURDEN_common_bottom <- summary_final_BURDEN_decrease[summary_final_BURDEN_decrease$phecode %in% common_phecode, ]

Summary_phenotypes_merge_anno_BURDEN_rare_bottom <- summary_final_BURDEN_decrease[summary_final_BURDEN_decrease$phecode %in% rare_phecode, ]

colnames(Summary_phenotypes_merge_anno_BURDEN_common_top)
# [1] "phecode"               "group.x"               "count.x"               "group.y"              
# [5] "count.y"               "delta_change"          "delta"                 "phecode_str"          
# [9] "exclude_range"         "exclude_name"          "frequency"             "do.call.rbind..SKAT."
str(Summary_phenotypes_merge_anno_BURDEN_common_top)
colnames(Summary_phenotypes_merge_anno_BURDEN_common_top)[12] <- "p.value"
class(Summary_phenotypes_merge_anno_BURDEN_common_top)
Summary_phenotypes_merge_anno_BURDEN_common_top <- Summary_phenotypes_merge_anno_BURDEN_common_top[!is.na(Summary_phenotypes_merge_anno_BURDEN_common_top$exclude_name), ] %>%
  filter(exclude_name != "NULL") %>%
  arrange(exclude_name)  #201
#Summary_phenotypes_merge_anno_BURDEN_common_top <- Summary_phenotypes_merge_anno_BURDEN_common_top[sort(as.factor(Summary_phenotypes_merge_anno_BURDEN_common_top$exclude_name)), ]

library(ggplot2)
tiff("Summary_phenotypes_merge_anno_BURDEN_common_top_phencode_109_ratio2_85K.tiff", units="in", width=24, height=10, res=600)
ggplot(Summary_phenotypes_merge_anno_BURDEN_common_top, aes(x = as.factor(phecode), y = -log10(p.value))) + 
  #geom_pointrange(aes(col = PHENOTYPE, shape = as.factor(Direction)), position=position_dodge(width=0.8)) + 
  #geom_point(aes(shape = as.factor(Direction), color = PHENOTYPE, y = estimate, size = -log10(p.value)), stat="identity", position = position_dodge(0.8))  +
  geom_point(aes(color = exclude_name, y = -log10(p.value), size = delta_change), stat="identity", position = position_dodge(0.8))  +
  
  #scale_shape_manual(values=c("\u25BA","\u25C4","\u25BC","\u25B2"))
  #scale_shape_manual(values=c("\u25BC","\u25B2")) +
  #geom_hline(aes(yintercept = -log(0.001)), linetype='dashed', col = 'black', width = 1) +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme_bw() +
  geom_text(data = subset(Summary_phenotypes_merge_anno_BURDEN_common_top, p.value < 0.05), 
            aes(phecode, -log10(p.value), label = phecode_str), vjust = -1, size = 4, angle = 30) +
  #facet_wrap(. ~ Header.x, scales="free") + 
  #theme(strip.text.x = element_text(size = 6, colour = "black")) +
  theme(text = element_text(size = 10),
        axis.text.x=element_text(angle = 45, vjust = 0.5, hjust=0.5)) + 
  ylab("-log10(p.value)") +
  xlab("phecode") + 
  ggtitle("Summary_phenotypes_merge_anno_BURDEN_common_top_85K_ratio2")
dev.off()


#bottom
colnames(Summary_phenotypes_merge_anno_BURDEN_common_bottom)
# [1] "phecode"               "group.x"               "count.x"               "group.y"              
# [5] "count.y"               "delta_change"          "delta"                 "phecode_str"          
# [9] "exclude_range"         "exclude_name"          "frequency"             "do.call.rbind..BURDEN."
str(Summary_phenotypes_merge_anno_BURDEN_common_bottom)
colnames(Summary_phenotypes_merge_anno_BURDEN_common_bottom)[12] <- "p.value"
class(Summary_phenotypes_merge_anno_BURDEN_common_bottom)
Summary_phenotypes_merge_anno_BURDEN_common_bottom <- Summary_phenotypes_merge_anno_BURDEN_common_bottom[!is.na(Summary_phenotypes_merge_anno_BURDEN_common_bottom$exclude_name), ] %>%
  filter(exclude_name != "NULL") %>%
  arrange(exclude_name)
#Summary_phenotypes_merge_anno_BURDEN_common_bottom <- Summary_phenotypes_merge_anno_BURDEN_common_bottom[sort(as.factor(Summary_phenotypes_merge_anno_BURDEN_common_bottom$exclude_name)), ]
library(ggplot2)
tiff("Summary_phenotypes_merge_anno_BURDEN_common_bottom_phencode_109_ratio2_90K.tiff", units="in", width=24, height=10, res=600)
ggplot(Summary_phenotypes_merge_anno_BURDEN_common_bottom, aes(x = as.factor(phecode), y = -log10(p.value))) + 
  #geom_pointrange(aes(col = PHENOTYPE, shape = as.factor(Direction)), position=position_dodge(width=0.8)) + 
  #geom_point(aes(shape = as.factor(Direction), color = PHENOTYPE, y = estimate, size = -log10(p.value)), stat="identity", position = position_dodge(0.8))  +
  geom_point(aes(color = exclude_name, y = -log10(p.value), size = -delta_change), stat="identity", position = position_dodge(0.8))  +
  
  #scale_shape_manual(values=c("\u25BA","\u25C4","\u25BC","\u25B2"))
  #scale_shape_manual(values=c("\u25BC","\u25B2")) +
  #geom_hline(aes(yintercept = -log(0.001)), linetype='dashed', col = 'black', width = 1) +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme_bw() +
  geom_text(data = subset(Summary_phenotypes_merge_anno_BURDEN_common_bottom, p.value < 0.05), 
            aes(phecode, -log10(p.value), label = phecode_str), vjust = -1, size = 4, angle = 30) +
  #facet_wrap(. ~ Header.x, scales="free") + 
  #theme(strip.text.x = element_text(size = 6, colour = "black")) +
  theme(text = element_text(size = 10),
        axis.text.x=element_text(angle = 45, vjust = 0.5, hjust=0.5)) + 
  ylab("-log10(p.value)") +
  xlab("phecode") + 
  ggtitle("Summary_phenotypes_merge_anno_BURDEN_common_bottom_90K_ratio2")
dev.off()

##########################################################################################
#to identify common disease

#read in 'Summary_phenotypes_merge_anno_common' and 'Summary_phenotypes_merge_anno_rare', two objects
load("Summary_phenotypes109_merge_anno_commonrare_delta_p123_85K_phenotypes109_age.RData")
load("Summary_phenotypes109_merge_anno_commonrare_delta_p123_90K_phenotypes109_age.RData")
load("Summary_phenotypes109_merge_anno_commonrare_delta_p123_85K_phenotypes109.RData")
load("Summary_phenotypes109_merge_anno_commonrare_delta_p123_90K_phenotypes109.RData")
unique(Summary_phenotypes_merge_anno_common$exclude_name)
Summary_phenotypes_merge_anno_common <- Summary_phenotypes_merge_anno_common[!Summary_phenotypes_merge_anno_common$exclude_name %in% c("", "NULL", NA), ]  #695


#this common disease is for all based on frequency in both discovery and replication
commondisease <- read.table(file = "Common_disease.txt", header = TRUE, stringsAsFactors=FALSE)

#for v1.2
phecode_mapping <- read.csv("Phecode_map_v1_2_icd10cm_beta.csv", header = T, stringsAsFactors = F)
phecode_mapping_unique <- phecode_mapping %>%
  dplyr::select(phecode, phecode_str, exclude_range, exclude_name) %>%
  distinct(phecode, .keep_all = T)  #1755*4
str(phecode_mapping_unique)
phecode_mapping_unique$phecodeX <- as.character(paste0("X", phecode_mapping_unique$phecode))

unique(phecode_mapping_unique$exclude_name)

#identify common disease from control group (n=711) remove phecodes not belonging to any category
phecode_mapping_unique_commondisease <- phecode_mapping_unique[!phecode_mapping_unique$exclude_name %in% c("", "NULL"), ]
phecode_mapping_unique_commondisease <- phecode_mapping_unique_commondisease[phecode_mapping_unique_commondisease$phecodeX %in% commondisease$phecodeX, ] #711
head(phecode_mapping_unique_commondisease)
# phecode                                    phecode_str exclude_range       exclude_name phecodeX
# 4   426.32                       Left bundle branch block    426-427.99 circulatory system  X426.32
# 5   433.10 Occlusion and stenosis of precerebral arteries    430-438.99 circulatory system   X433.1
# 10  773.00                                   Pain in limb    773-773.99           symptoms     X773
# 13  447.00     Other disorders of arteries and arterioles    440-449.99 circulatory system     X447
# 15  172.20            Other non-epithelial cancer of skin    172-173.99          neoplasms   X172.2
# 16  316.00              Substance addiction and disorders    316-319.99   mental disorders     X316


#load("Summary_phenotypes109_merge_anno_v1.2b_matched_ratio2_p123_90K_phenotypes109_WHITE.RData")
#load("Summary_phenotypes109_merge_anno_v1.2b_matched_ratio2_p123_85K_phenotypes109_WHITE.RData")

load("Summary_phenotypes109_merge_anno_v1.2b_matched_ratio2_p123_90K_phenotypes109.RData")
#load("Summary_phenotypes109_merge_anno_v1.2b_matched_ratio2_p123_85K_phenotypes109.RData")



#n=696 PheCodes having results
Summary_phenotypes_merge_anno_common <- Summary_phenotypes_merge_anno[Summary_phenotypes_merge_anno$phecode %in% phecode_mapping_unique_commondisease$phecode, ] #696
head(Summary_phenotypes_merge_anno_common)
Summary_phenotypes_merge_anno_common$phecode <- paste0("X", Summary_phenotypes_merge_anno_common$phecode)
Summary_phenotypes_merge_anno_common$direction <- ifelse(Summary_phenotypes_merge_anno_common$delta_change >0, "top", "bottom")
#n=224/290 PheCodes having PRE > 0 for 90K/85K; 208/695 for discovery and 248/580 for replication
Summary_phenotypes_merge_anno_common_top <- Summary_phenotypes_merge_anno_common[Summary_phenotypes_merge_anno_common$delta_change >0, ]
Summary_phenotypes_merge_anno_common_top$phecode <- paste0("X", Summary_phenotypes_merge_anno_common_top$phecode)

#load("summary_final_SKATO_12032021_matched_ratio2_p123_85K_WHITE.RData") #this is not linearweighted
#load("summary_final_SKATO_12032021_matched_ratio2_p123_90K_WHITE.RData") #this is not linearweighted

load("summary_final_12032021_matched_ratio2_90K_ALLSIXGENES_linearweighted_ALL.RData")
load("summary_final_12032021_matched_ratio2_85K_ALLSIXGENES_linearweighted_ALL.RData")

colnames(summary_final_SKATO) <- "pvalue"
summary_final_SKATO$phecode <- rownames(summary_final_SKATO)
summary_final_SKATO <- merge(Summary_phenotypes_merge_anno_common_top, summary_final_SKATO, by = "phecode")
summary_final_SKATO_85K <- summary_final_SKATO
summary_final_SKATO_90K <- summary_final_SKATO

summary_final_SKATO_85K_common <- merge(Summary_phenotypes_merge_anno_common, summary_final_SKATO, by = "phecode")
summary_final_SKATO_90K_common <- merge(Summary_phenotypes_merge_anno_common, summary_final_SKATO, by = "phecode")



summary_final_SKATO_final <- merge(summary_final_SKATO_90K, summary_final_SKATO_85K, by = "phecode") #208/248

summary_final_SKATO_final_replicated <- summary_final_SKATO_final %>%
  filter(pvalue.x < 0.05) %>%
  filter(pvalue.y < 0.05)

#BH correction
#https://sphweb.bumc.bu.edu/otlt/MPH-Modules/BS/R/R-Manual/R-Manual22.html
str(summary_final_SKATO_90K)
summary_final_SKATO_90K$qvalue <- p.adjust(summary_final_SKATO_90K$pvalue, method="BH")
summary_final_SKATO_90K_common$qvalue <- p.adjust(summary_final_SKATO_90K_common$pvalue, method="BH")
save(summary_final_SKATO_90K, file = "summmary_final_SKATO_90K_qvalue.RData", version = 2)
save(summary_final_SKATO_90K_common, file = "summary_final_SKATO_90K_common_qvalue.RData", version = 2)

str(summary_final_SKATO_85K)
summary_final_SKATO_85K$qvalue <- p.adjust(summary_final_SKATO_85K$pvalue, method="BH")
save(summary_final_SKATO_85K, file = "summmary_final_SKATO_85K_qvalue.RData", version = 2)

summary_final_SKATO_85K_common$qvalue <- p.adjust(summary_final_SKATO_85K_common$pvalue, method="BH")

save(summary_final_SKATO_85K_common, file = "summary_final_SKATO_85K_common_qvalue.RData", version = 2)
library(ggplot2)
tiff("Summary_phenotypes_merge_anno_SKATO_common_top_phencode_109_ratio2_90K_WHITE.tiff", units="in", width=24, height=10, res=600)
ggplot(summary_final_SKATO_90K, aes(x = as.factor(phecode), y = -log10(pvalue))) + 
  #geom_pointrange(aes(col = PHENOTYPE, shape = as.factor(Direction)), position=position_dodge(width=0.8)) + 
  #geom_point(aes(shape = as.factor(Direction), color = PHENOTYPE, y = estimate, size = -log10(p.value)), stat="identity", position = position_dodge(0.8))  +
  geom_point(aes(color = exclude_name, y = -log10(pvalue), size = delta_change), stat="identity", position = position_dodge(0.8))  +
  
  #scale_shape_manual(values=c("\u25BA","\u25C4","\u25BC","\u25B2"))
  #scale_shape_manual(values=c("\u25BC","\u25B2")) +
  #geom_hline(aes(yintercept = -log(0.001)), linetype='dashed', col = 'black', width = 1) +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme_bw() +
  geom_text(data = subset(summary_final_SKATO, pvalue < 0.05), 
            aes(phecode, -log10(pvalue), label = phecode_str), vjust = -1, size = 4, angle = 30) +
  #facet_wrap(. ~ Header.x, scales="free") + 
  #theme(strip.text.x = element_text(size = 6, colour = "black")) +
  theme(text = element_text(size = 10),
        axis.text.x=element_text(angle = 45, vjust = 0.5, hjust=0.5)) + 
  ylab("-log10(pvalue)") +
  xlab("phecode") + 
  ggtitle("Summary_phenotypes_merge_anno_SKATO_common_top_90K_ratio2_WHITE")
dev.off()

tiff("Summary_phenotypes_merge_anno_SKATO_common_top_phencode_109_ratio2_85K_WHITE.tiff", units="in", width=24, height=10, res=600)
ggplot(summary_final_SKATO_85K, aes(x = as.factor(phecode), y = -log10(pvalue))) + 
  #geom_pointrange(aes(col = PHENOTYPE, shape = as.factor(Direction)), position=position_dodge(width=0.8)) + 
  #geom_point(aes(shape = as.factor(Direction), color = PHENOTYPE, y = estimate, size = -log10(p.value)), stat="identity", position = position_dodge(0.8))  +
  geom_point(aes(color = exclude_name, y = -log10(pvalue), size = delta_change), stat="identity", position = position_dodge(0.8))  +
  
  #scale_shape_manual(values=c("\u25BA","\u25C4","\u25BC","\u25B2"))
  #scale_shape_manual(values=c("\u25BC","\u25B2")) +
  #geom_hline(aes(yintercept = -log(0.001)), linetype='dashed', col = 'black', width = 1) +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme_bw() +
  geom_text(data = subset(summary_final_SKATO, pvalue < 0.05), 
            aes(phecode, -log10(pvalue), label = phecode_str), vjust = -1, size = 4, angle = 30) +
  #facet_wrap(. ~ Header.x, scales="free") + 
  #theme(strip.text.x = element_text(size = 6, colour = "black")) +
  theme(text = element_text(size = 10),
        axis.text.x=element_text(angle = 45, vjust = 0.5, hjust=0.5)) + 
  ylab("-log10(pvalue)") +
  xlab("phecode") + 
  ggtitle("Summary_phenotypes_merge_anno_SKATO_common_top_85K_ratio2_WHITE")
dev.off()

tiff("Summary_phenotypes_merge_anno_BURDEN_common_top_phencode_109_ratio2_90K_ALLSIXGENES_linearweighted.tiff", units="in", width=24, height=10, res=600)
ggplot(summary_final_BURDEN_90K, aes(x = as.factor(phecode), y = -log10(pvalue))) + 
  #geom_pointrange(aes(col = PHENOTYPE, shape = as.factor(Direction)), position=position_dodge(width=0.8)) + 
  #geom_point(aes(shape = as.factor(Direction), color = PHENOTYPE, y = estimate, size = -log10(p.value)), stat="identity", position = position_dodge(0.8))  +
  geom_point(aes(color = exclude_name, y = -log10(pvalue), size = delta_change), stat="identity", position = position_dodge(0.8))  +
  
  #scale_shape_manual(values=c("\u25BA","\u25C4","\u25BC","\u25B2"))
  #scale_shape_manual(values=c("\u25BC","\u25B2")) +
  #geom_hline(aes(yintercept = -log(0.001)), linetype='dashed', col = 'black', width = 1) +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme_bw() +
  geom_text(data = subset(summary_final_BURDEN, pvalue < 0.05), 
            aes(phecode, -log10(pvalue), label = phecode_str), vjust = -1, size = 4, angle = 30) +
  #facet_wrap(. ~ Header.x, scales="free") + 
  #theme(strip.text.x = element_text(size = 6, colour = "black")) +
  theme(text = element_text(size = 10),
        axis.text.x=element_text(angle = 45, vjust = 0.5, hjust=0.5)) + 
  ylab("-log10(pvalue)") +
  xlab("phecode") + 
  ggtitle("Summary_phenotypes_merge_anno_BURDEN_common_top_90K_ratio2_linearweighted")
dev.off()

#put results from both cohort together 
summary_final_SKATO_90K$cohort <- "discovery"
summary_final_SKATO_85K$cohort <- "replication"
Summary_phenotypes_merge_anno_SKATO_common_90K85K <- rbind(summary_final_SKATO_90K, summary_final_SKATO_85K) %>%
  arrange(exclude_name)
Summary_phenotypes_merge_anno_SKATO_common_90K85K

tiff("Summary_phenotypes_merge_anno_SKATO_common_top_phencode_109_ratio2_90K85K_ALLSIXGENES_linearweighted_commonfromeach.tiff", units="in", width=24, height=10, res=600)
ggplot(Summary_phenotypes_merge_anno_SKATO_common_90K85K, aes(x = as.factor(phecode), y = -log10(pvalue))) + 
  #geom_pointrange(aes(col = PHENOTYPE, shape = as.factor(cohort)), position=position_dodge(width=0.8)) + 
  #geom_point(aes(shape = as.factor(cohort), color = PHENOTYPE, y = estimate, size = -log10(pvalue)), stat="identity", position = position_dodge(0.8))  +
  geom_point(aes(shape = as.factor(cohort), color = exclude_name, y = -log10(pvalue), size = abs(delta_change)), stat="identity", position = position_dodge(0.8))  +
  ylim(0, 5) +
  #scale_shape_manual(values=c("\u25BA","\u25C4","\u25BC","\u25B2"))
  scale_shape_manual(values=c("\u25BA","\u25C4")) +
  geom_hline(aes(yintercept = -log10(0.05)), linetype='dashed', col = 'black', width = 1) +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme_bw() +
  geom_text(data = subset(Summary_phenotypes_merge_anno_SKATO_common_90K85K, pvalue < 0.05), 
            aes(phecode, -log10(pvalue), label = phecode_str), vjust = -1, size = 4, angle = 30) +
  #facet_wrap(. ~ Header.x, scales="free") + 
  #theme(strip.text.x = element_text(size = 6, colour = "black")) +
  theme(text = element_text(size = 10),
        axis.text.x=element_text(angle = 45, vjust = 0.5, hjust=0.5)) + 
  ylab("-log10(pvalue)") +
  xlab("phecode") + 
  ggtitle("Summary_phenotypes_merge_anno_SKATO_common_90K85K")
dev.off()

#narrow down to CS
unique(Summary_phenotypes_merge_anno_SKATO_common_90K85K$exclude_name)
Summary_phenotypes_merge_anno_SKATO_common_90K85K_CS <- Summary_phenotypes_merge_anno_SKATO_common_90K85K[Summary_phenotypes_merge_anno_SKATO_common_90K85K$exclude_name == "circulatory system", ]
tiff("Summary_phenotypes_merge_anno_SKATO_common_top_phencode_109_ratio2_90K85K_ALLSIXGENES_linearweighted_commonfromeach_CSonly.tiff", units="in", width=10, height=6, res=600)
ggplot(Summary_phenotypes_merge_anno_SKATO_common_90K85K_CS, aes(x = as.factor(phecode), y = -log10(pvalue))) + 
  #geom_pointrange(aes(col = PHENOTYPE, shape = as.factor(cohort)), position=position_dodge(width=0.8)) + 
  #geom_point(aes(shape = as.factor(cohort), color = PHENOTYPE, y = estimate, size = -log10(pvalue)), stat="identity", position = position_dodge(0.8))  +
  geom_point(aes(shape = as.factor(cohort), color = exclude_name, y = -log10(pvalue), size = abs(delta_change)), stat="identity", position = position_dodge(0.8))  +
  ylim(1.2, 2.1) +
  #scale_shape_manual(values=c("\u25BA","\u25C4","\u25BC","\u25B2"))
  scale_shape_manual(values=c("\u25BA","\u25C4")) +
  geom_hline(aes(yintercept = -log10(0.05)), linetype='dashed', col = 'black', width = 1) +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme_bw() +
  geom_text(data = subset(Summary_phenotypes_merge_anno_SKATO_common_90K85K_CS, pvalue < 0.05), 
            aes(phecode, -log10(pvalue), label = phecode_str), vjust = -1, size = 4, angle = 30) +
  #facet_wrap(. ~ Header.x, scales="free") + 
  #theme(strip.text.x = element_text(size = 6, colour = "black")) +
  theme(text = element_text(size = 10),
        axis.text.x=element_text(angle = 45, vjust = 0.5, hjust=0.5)) + 
  ylab("-log10(pvalue)") +
  xlab("phecode") + 
  ggtitle("Summary_phenotypes_merge_anno_SKATO_common_90K85K_CSonly")
dev.off()

##############################################################################
load("summary_final_12032021_matched_ratio2_90K_chr19_linearweighted.RData")

colnames(summary_final_SKATO) <- "pvalue"
summary_final_SKATO$phecode <- rownames(summary_final_SKATO)
summary_final_SKATO_90K <- summary_final_SKATO

load("summary_final_12032021_matched_ratio2_85K_chr19_linearweighted.RData")
colnames(summary_final_SKATO) <- "pvalue"
summary_final_SKATO$phecode <- rownames(summary_final_SKATO)
summary_final_SKATO_85K <- summary_final_SKATO

summary_final_SKATO_NOTCH3 <- cbind(summary_final_SKATO_90K, summary_final_SKATO_85K)
summary_final_SKATO_NOTCH3$Gene <- "NOTCH3"

load("summary_final_12032021_matched_ratio2_90K_chr17_linearweighted.RData")

colnames(summary_final_SKATO) <- "pvalue"
summary_final_SKATO$phecode <- rownames(summary_final_SKATO)
summary_final_SKATO_90K <- summary_final_SKATO

load("summary_final_12032021_matched_ratio2_85K_chr17_linearweighted.RData")
colnames(summary_final_SKATO) <- "pvalue"
summary_final_SKATO$phecode <- rownames(summary_final_SKATO)
summary_final_SKATO_85K <- summary_final_SKATO

summary_final_SKATO_CTC1 <- cbind(summary_final_SKATO_90K, summary_final_SKATO_85K)
summary_final_SKATO_CTC1$Gene <- "CTC1"

load("summary_final_12032021_matched_ratio2_90K_chr13_linearweighted.RData")

colnames(summary_final_SKATO) <- "pvalue"
summary_final_SKATO$phecode <- rownames(summary_final_SKATO)
summary_final_SKATO_90K <- summary_final_SKATO

load("summary_final_12032021_matched_ratio2_85K_chr13_linearweighted.RData")
colnames(summary_final_SKATO) <- "pvalue"
summary_final_SKATO$phecode <- rownames(summary_final_SKATO)
summary_final_SKATO_85K <- summary_final_SKATO

summary_final_SKATO_COL4A12 <- cbind(summary_final_SKATO_90K, summary_final_SKATO_85K)
summary_final_SKATO_COL4A12$Gene <- "COL4A1/2"

load("summary_final_12032021_matched_ratio2_90K_chr10_linearweighted.RData")

colnames(summary_final_SKATO) <- "pvalue"
summary_final_SKATO$phecode <- rownames(summary_final_SKATO)
summary_final_SKATO_90K <- summary_final_SKATO

load("summary_final_12032021_matched_ratio2_85K_chr10_linearweighted.RData")
colnames(summary_final_SKATO) <- "pvalue"
summary_final_SKATO$phecode <- rownames(summary_final_SKATO)
summary_final_SKATO_85K <- summary_final_SKATO

summary_final_SKATO_HTRA1 <- cbind(summary_final_SKATO_90K, summary_final_SKATO_85K)
summary_final_SKATO_HTRA1$Gene <- "HTRA1"

load("summary_final_12032021_matched_ratio2_90K_chr3_linearweighted.RData")

colnames(summary_final_SKATO) <- "pvalue"
summary_final_SKATO$phecode <- rownames(summary_final_SKATO)
summary_final_SKATO_90K <- summary_final_SKATO

load("summary_final_12032021_matched_ratio2_85K_chr3_linearweighted.RData")
colnames(summary_final_SKATO) <- "pvalue"
summary_final_SKATO$phecode <- rownames(summary_final_SKATO)
summary_final_SKATO_85K <- summary_final_SKATO

summary_final_SKATO_TREX1 <- cbind(summary_final_SKATO_90K, summary_final_SKATO_85K)
summary_final_SKATO_TREX1$Gene <- "TREX1"

load("summary_final_12032021_matched_ratio2_90K_chrX_linearweighted.RData")

colnames(summary_final_SKATO) <- "pvalue"
summary_final_SKATO$phecode <- rownames(summary_final_SKATO)
summary_final_SKATO_90K <- summary_final_SKATO

load("summary_final_12032021_matched_ratio2_85K_chrX_linearweighted.RData")
colnames(summary_final_SKATO) <- "pvalue"
summary_final_SKATO$phecode <- rownames(summary_final_SKATO)
summary_final_SKATO_85K <- summary_final_SKATO

summary_final_SKATO_GLA <- cbind(summary_final_SKATO_90K, summary_final_SKATO_85K)
summary_final_SKATO_GLA$Gene <- "GLA"

load("summary_final_12032021_matched_ratio2_90K_ALLSIXGENES_linearweighted.RData")

colnames(summary_final_SKATO) <- "pvalue"
summary_final_SKATO$phecode <- rownames(summary_final_SKATO)
summary_final_SKATO_90K <- summary_final_SKATO

load("summary_final_12032021_matched_ratio2_85K_ALLSIXGENES_linearweighted.RData")
colnames(summary_final_SKATO) <- "pvalue"
summary_final_SKATO$phecode <- rownames(summary_final_SKATO)
summary_final_SKATO_85K <- summary_final_SKATO

summary_final_SKATO_ALL <- cbind(summary_final_SKATO_90K, summary_final_SKATO_85K)
summary_final_SKATO_ALL$Gene <- "All"
foo <- list(summary_final_SKATO_ALL, summary_final_SKATO_COL4A12, summary_final_SKATO_NOTCH3, summary_final_SKATO_HTRA1, summary_final_SKATO_TREX1, summary_final_SKATO_CTC1, summary_final_SKATO_GLA)
summary_final_SKATO_ALL_summary <- do.call(rbind, foo)

unique(summary_final_SKATO_ALL_summary$phecode)
#[1] "X425.11" "X425.1"  "X425"    "X433"    "X433.1"  "X433.2"  "X433.3"  "X286.8"  "X286.81"

write.table(summary_final_SKATO_ALL_summary, file = "summary_final_SKATO_ALL_summary.txt",  quote = F, row.name = F, col.name = T)
#######################################################################################
# Shape of the logistic weight
MAF <- 1:1000/1000
W <- Get_Logistic_Weights_MAF(MAF, par1=0.07, par2=150)
par(mfrow=c(1,2))
plot(MAF,W,xlab="MAF",ylab="Weights",type="l")
plot(MAF[1:100],W[1:100],xlab="MAF",ylab="Weights",type="l")
par(mfrow=c(1,2))
# Use logistic weight
weights<-Get_Logistic_Weights(Z, par1=0.07, par2=150)
SKAT(Z, obj, kernel = "linear.weighted", weights=weights)$p.value

colnames(dat)

weights_CADD <- c(31, 27.8, 28.5, 29.3, 
32, 
0.594, 
18.02, 
25.2, 
5.108, 
25.5, 
26.9, 
37, 
25.4, 
32, 
32, 
35, 
35, 
21.9, 
33, 
30, 
0.326, 
2.981, 
25.8, 
33, 
2.309, 
31, 
44, 
1.476, 
34, 
33, 
3.883, 
34, 
32, 
37, 
34, 
34, 
9.169, 
40, 
22.6, 
41, 
45, 
49, 
22.3, 
8.543, 
26.2, 
24.2, 
25.5, 
24.7, 
23.3, 
5.635, 
44, 
50, 
28.6, 
31, 
23.3, 
20.4, 
34, 
26.4, 
22.2, 
25.1, 
29.4, 
23.1, 
24.9, 
28.6, 
20.9, 
23.7, 
24.3) 

weights_125 <- c(0.999870215199605, 
                 0.999740446541552, 
                 0.999870215199605, 
                 0.999870215199605, 
                 0.999870215199605, 
                 0.999740446541552, 
                 0.999870215199605, 
                 0.950589769792554, 
                 0.94504152128383, 
                 0.985820032092687, 
                 0.981221314057172, 
                 0.996242787568163, 
                 0.999870215199605, 
                 0.991597480312237, 
                 0.990310805561371, 
                 0.988640520223112, 
                 0.99844364736759, 
                 0.996889616192861, 
                 0.993659492841176, 
                 0.999870215199605, 
                 0.996630836461704, 
                 0.994175637684247, 
                 0.994433806446441, 
                 0.998573254707789, 
                 0.999870215199605, 
                 0.999870215199605, 
                 0.998832517759547, 
                 0.999480957644803, 
                 0.999870215199605, 
                 0.999740446541552, 
                 0.999610694023925, 
                 0.998832517759547, 
                 0.999870215199605, 
                 0.999740446541552, 
                 0.999870215199605, 
                 0.999351237402262, 
                 0.999221533294387, 
                 0.997795852446721, 
                 0.999740446541552, 
                 0.998314056148626, 
                 0.999480957644803, 
                 0.999870215199605, 
                 0.995079509468143, 
                 0.999480957644803, 
                 0.999870215199605, 
                 0.997925379199938, 
                 0.999740446541552, 
                 0.999870215199605, 
                 0.999610694023925, 
                 0.997795852446721, 
                 0.999480957644803, 
                 0.999351237402262, 
                 0.999870215199605, 
                 0.999870215199605, 
                 0.999740446541552, 
                 0.999870215199605, 
                 0.999610694023925, 
                 0.999870215199605, 
                 0.999610694023925, 
                 0.999610694023925, 
                 0.999480957644803, 
                 0.999870215199605, 
                 0.999870215199605, 
                 0.871688100313487, 
                 0.999870215199605, 
                 0.998573254707789, 
                 0.999740446541552)

weights_150 <- c(0.999735040610632, 
                 0.999470149992387, 
                 0.999735040610632, 
                 0.999735040610632, 
                 0.999735040610632, 
                 0.999470149992387, 
                 0.999735040610632, 
                 0.901715052692282, 
                 0.891002454764131, 
                 0.971263002935384, 
                 0.962035069718055, 
                 0.992344035023909, 
                 0.999735040610632, 
                 0.982919923128185, 
                 0.980317710735127, 
                 0.976944921735305, 
                 0.996825022378877, 
                 0.993659920156559, 
                 0.987097543669265, 
                 0.999735040610632, 
                 0.993133361114436, 
                 0.988144663526561, 
                 0.988668631423231, 
                 0.997089226246439, 
                 0.999735040610632, 
                 0.999735040610632, 
                 0.997617839788611, 
                 0.998940574999376, 
                 0.999735040610632, 
                 0.999470149992387, 
                 0.999205328127794, 
                 0.997617839788611, 
                 0.999735040610632, 
                 0.999470149992387, 
                 0.999735040610632, 
                 0.998675890589658, 
                 0.998411274881184, 
                 0.995505031379154, 
                 0.999470149992387, 
                 0.996560887090428, 
                 0.998940574999376, 
                 0.999735040610632, 
                 0.989979742181389, 
                 0.998940574999376, 
                 0.999735040610632, 
                 0.995768892525435, 
                 0.999470149992387, 
                 0.999735040610632, 
                 0.999205328127794, 
                 0.995505031379154, 
                 0.998940574999376, 
                 0.998675890589658, 
                 0.999735040610632, 
                 0.999735040610632, 
                 0.999470149992387, 
                 0.999735040610632, 
                 0.999205328127794, 
                 0.999735040610632, 
                 0.999205328127794, 
                 0.999205328127794, 
                 0.998940574999376, 
                 0.999735040610632, 
                 0.999735040610632, 
                 0.755504892742761, 
                 0.999735040610632, 
                 0.997089226246439, 
                 0.999470149992387)
weights_1200 <- c(0.998924376755707, 
                  0.997849904669091, 
                  0.998924376755707, 
                  0.998924376755707, 
                  0.998924376755707, 
                  0.997849904669091, 
                  0.998924376755707, 
                  0.656941080366981, 
                  0.625816433804458, 
                  0.888325418140416, 
                  0.854544114690136, 
                  0.969269847616002, 
                  0.998924376755707, 
                  0.932426223038577, 
                  0.922441503486914, 
                  0.909620249240217, 
                  0.987168228445943, 
                  0.974500294420047, 
                  0.948625911243018, 
                  0.998924376755707, 
                  0.972404754192026, 
                  0.952719398116843, 
                  0.954772727481804, 
                  0.988231256159793, 
                  0.998924376755707, 
                  0.998924376755707, 
                  0.990360729676044, 
                  0.99570440906714, 
                  0.998924376755707, 
                  0.997849904669091, 
                  0.996776582514398, 
                  0.990360729676044, 
                  0.998924376755707, 
                  0.997849904669091, 
                  0.998924376755707, 
                  0.994633383104102, 
                  0.993563503403434, 
                  0.981870131844117, 
                  0.997849904669091, 
                  0.986106338477418, 
                  0.99570440906714, 
                  0.998924376755707, 
                  0.959925335060158, 
                  0.99570440906714, 
                  0.998924376755707, 
                  0.982927482934323, 
                  0.997849904669091, 
                  0.998924376755707, 
                  0.996776582514398, 
                  0.981870131844117, 
                  0.99570440906714, 
                  0.994633383104102, 
                  0.998924376755707, 
                  0.998924376755707, 
                  0.997849904669091, 
                  0.998924376755707, 
                  0.996776582514398, 
                  0.998924376755707, 
                  0.996776582514398, 
                  0.996776582514398, 
                  0.99570440906714, 
                  0.998924376755707, 
                  0.998924376755707, 
                  0.320253780071286, 
                  0.998924376755707, 
                  0.988231256159793, 
                  0.997849904669091)
weights_55 <- c(430.012790513733, 
                304.06578237417, 
                430.012790513733, 
                430.012790513733, 
                430.012790513733, 
                304.06578237417, 
                430.012790513733, 
                21.7974842291542, 
                20.6417736111316, 
                41.0122005656621, 
                35.6020735282693, 
                79.8574170218983, 
                430.012790513733, 
                53.3457550552544, 
                49.6635386620202, 
                45.8503050822325, 
                124.13769263027, 
                87.7814528425171, 
                61.4383734684747, 
                430.012790513733, 
                84.3381480562136, 
                64.1101501677541, 
                65.5837845106486, 
                129.657240675143, 
                430.012790513733, 
                430.012790513733, 
                143.340697655028, 
                215.008139431122, 
                430.012790513733, 
                304.06578237417, 
                248.269343018608, 
                143.340697655028, 
                430.012790513733, 
                304.06578237417, 
                430.012790513733, 
                192.309646214226, 
                175.554360144608, 
                104.297933419002, 
                304.06578237417, 
                119.267959755148, 
                215.008139431122, 
                430.012790513733, 
                69.764266190698, 
                215.008139431122, 
                430.012790513733, 
                107.507558276319, 
                304.06578237417, 
                430.012790513733, 
                248.269343018608, 
                104.297933419002, 
                215.008139431122, 
                192.309646214226, 
                430.012790513733, 
                430.012790513733, 
                304.06578237417, 
                430.012790513733, 
                248.269343018608, 
                430.012790513733, 
                248.269343018608, 
                248.269343018608, 
                215.008139431122, 
                430.012790513733, 
                430.012790513733, 
                13.2768927237404, 
                430.012790513733, 
                129.657240675143, 
                304.06578237417)

weights_51 <- c(430.009534782537, 
                304.061178067486, 
                430.009534782537, 
                430.009534782537, 
                430.009534782537, 
                304.061178067486, 
                430.009534782537, 
                21.7331480443692, 
                20.5738220521353, 
                40.9780481274239, 
                35.5627251182721, 
                79.8398835762722, 
                430.009534782537, 
                53.3195037849309, 
                49.6353398165916, 
                45.8197593038643, 
                124.126414245286, 
                87.7655024869396, 
                61.4155815749281, 
                430.009534782537, 
                84.3215463458923, 
                64.0883084995565, 
                65.5624337984408, 
                129.646442461148, 
                430.009534782537, 
                430.009534782537, 
                143.330930334662, 
                215.001627937035, 
                430.009534782537, 
                304.061178067486, 
                248.263703908461, 
                143.330930334662, 
                430.009534782537, 
                304.061178067486, 
                430.009534782537, 
                192.302366130708, 
                175.546385199742, 
                104.284509346902, 
                304.061178067486, 
                119.256220820826, 
                215.001627937035, 
                430.009534782537, 
                69.7441953107558, 
                215.001627937035, 
                430.009534782537, 
                107.494535034582, 
                304.061178067486, 
                430.009534782537, 
                248.263703908461, 
                104.284509346902, 
                215.001627937035, 
                192.302366130708, 
                430.009534782537, 
                430.009534782537, 
                304.061178067486, 
                430.009534782537, 
                248.263703908461, 
                430.009534782537, 
                248.263703908461, 
                248.263703908461, 
                215.001627937035, 
                430.009534782537, 
                430.009534782537, 
                13.1709624435612, 
                430.009534782537, 
                129.646442461148, 
                304.061178067486)

weights_55 <- weights_55/max(weights_55)
weights_5 <- c(430.011627749762, 
               304.064137970922, 
               430.011627749762, 
               430.011627749762, 
               430.011627749762, 
               304.064137970922, 
               430.011627749762, 
               21.7744851862956, 
               20.6174794718005, 
               41, 
               35.5880155293193, 
               79.8511546350566, 
               430.011627749762, 
               53.3363781182153, 
               49.6534658071989, 
               45.8393935387457, 
               124.133664517997, 
               87.7757559542117, 
               61.4302325356803, 
               430.011627749762, 
               84.3322184987628, 
               64.1023487175869, 
               65.5761584582349, 
               129.653384066904, 
               430.011627749762, 
               430.011627749762, 
               143.337209249921, 
               215.005813874881, 
               430.011627749762, 
               304.064137970922, 
               248.267329035994, 
               143.337209249921, 
               430.011627749762, 
               304.064137970922, 
               430.011627749762, 
               192.307046152761, 
               175.551511908423, 
               104.29313890918, 
               304.064137970922, 
               119.263767145962, 
               215.005813874881, 
               430.011627749762, 
               69.7570973563189, 
               215.005813874881, 
               430.011627749762, 
               107.502906937441, 
               304.064137970922, 
               430.011627749762, 
               248.267329035994, 
               104.29313890918, 
               215.005813874881, 
               192.307046152761, 
               430.011627749762, 
               430.011627749762, 
               304.064137970922, 
               430.011627749762, 
               248.267329035994, 
               430.011627749762, 
               248.267329035994, 
               248.267329035994, 
               215.005813874881, 
               430.011627749762, 
               430.011627749762, 
               13.2389630326584, 
               430.011627749762, 
               129.653384066904, 
               304.064137970922)
#to visualize beta distribution
#https://homepage.divms.uiowa.edu/~mbognar/applets/beta.html
