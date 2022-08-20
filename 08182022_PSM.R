###outline:
#1. identifying patient ID for RVs carriers and noncarriers and split based on the candidate genes;
#2. getting sequencing QC data from carriers using bash inquiry of vcf.gz file;
#3. random sampling from MyCode EHR to get controls (1:5 ratio);
#4. propensity score matching based on index age, sex, and ethnicity from controls (1:2 ratio);
#5. principal component analysis of selected RVs carriers and noncarriers;

######################################################################################################
##identifying patient ID for RVs carriers and noncarriers and split based on the candidate genes
Genotype<- read.table("discovery.raw", header =T, stringsAsFactors = F) 
#Genotype<- read.table("replication.raw", header =T, stringsAsFactors = F)

library(tidyverse)
Genotype <- as_tibble(Genotype)
Genotype <- cSplit(Genotype, "FID", "_") #FID_2 for PT####

Genotype_TREX1 <- Genotype %>% 
  dplyr::select(starts_with("X3.")) #X3. for TREX1(15); X10. for HTRA1(14); X13 for COL4A1/2; X19 for NOTCH3; X17 for CTC1; X23 for GLA

#optional
Genotype_TREX1 <- Genotype %>%
  dplyr::select(starts_with("X"))  #select all variants

Genotype_TREX1$PT_ID <- Genotype$FID_2
colnames(Genotype_TREX1)

gathercols <- colnames(Genotype_TREX1)[1:(ncol(Genotype_TREX1)-1)] 
keycol <- "genotype"
valuecol <- "dosage"
Genotype_long <- data.frame(gather_(Genotype_TREX1, keycol, valuecol, gathercols, factor_key=TRUE))
head(Genotype_long)

#90K or 85K
Genotype_long_mutation <- Genotype_long[!is.na(Genotype_long$dosage),] %>%
  filter(dosage != 2) %>% #only select RV carriers
  distinct(PT_ID, .keep_all = T) #some patients may have multiple mutations, here distinct to get the unique individual ID

TREX1 <- as.character(unique(Genotype_long_mutation$PT_ID)) 
# HTRA1 <- as.character(unique(Genotype_long_mutation$PT_ID)) 
# COL4A1 <- as.character(unique(Genotype_long_mutation$PT_ID)) 
# NOTCH3 <- as.character(unique(Genotype_long_mutation$PT_ID)) 
# CTC1 <- as.character(unique(Genotype_long_mutation$PT_ID)) 
# GLA <- as.character(unique(Genotype_long_mutation$PT_ID)) 

ALL <- as.character(unique(Genotype_long_mutation$PT_ID))

save(TREX1, HTRA1, COL4A1, NOTCH3, CTC1, GLA, ALL, file = "Individualgene_PTID_discovery.RData", version = 2)

##Getting sequencing QC data from carriers
load("Individualgene_PTID_discovery.RData")
ID <- as.integer(gsub("PT", "", CTC1))

##create bash file to extract QC sequencing data from vcf.gz batch calling file
STLIST <- NULL
for (i in colnames(Genotype_TREX1)[1:(ncol(Genotype_TREX1)-1)]) {
  STLIST[[i]] <- as.integer(unlist(strsplit(i, "\\."))[2]) #Escape special characters
  print(STLIST[i])
  
}
STLIST_final <- data.frame(do.call(rbind,STLIST))
colnames(STLIST_final) <- "coordinate"
STLIST_final$keep <- STLIST_final$coordinate

STLIST_final$coordinate <- paste0("<path_to_directory_of_server>/GHS_Freeze_175.967_X_87625896_101892854.GL.vcf.gz | sed '/^#/d' | awk '$2 == ", STLIST_final$coordinate)
STLIST_final$coordinate <- paste0(STLIST_final$coordinate, "' > /<local_directory>/GHS_Freeze_175.967_X_87625896_101892854.GL.vcf.test.")
STLIST_final$coordinate <- paste0(STLIST_final$coordinate,STLIST_final$keep)

rownames(STLIST_final) <- NULL
STLIST_final <- STLIST_final[, 1]
head(STLIST_final)
write.table(STLIST_final, file = "STLIST_final_GLA.txt", col.names = F, row.names = F, quote = F)

#create bash file
#sed -i -e 's/\r$//' STLIST_final_NOTCH3.sh
#chmod +x STLIST_final_NOTCH3.sh
#./STLIST_final_NOTCH3.sh

################################################################################################################
##Ramdom sampling from MyCode database to get controls

FID_90K <- cSplit(FID_90K, "FID", "_")
head(FID_90K)
FID_90K_excludecase <- FID_90K[!(FID_90K$FID_2 %in% df_ALL_90K$ALL), ]  
set.seed(42)
FID_90K_excludecase <- FID_90K_excludecase[sample(nrow(FID_90K_excludecase), N), ]  #N here referred to the number of 5 time sampling from the entire MyCOde discovery or replication cohort.
head(FID_90K_excludecase)

df_ALL_90K_control <- data.frame(FID_90K_excludecase$FID_2)
colnames(df_ALL_90K_control) <- "ALL"
df_ALL_90K_control$Cohort <- "discovery"

#same way to get replication sampling
df_ALL_control <- rbind(df_ALL_90K_control, df_ALL_85K_control)
write.table(df_ALL_control, file = "df_ALL_control.txt", col.names = T, row.names = F, quote = F)

############################################################################################
##propensity score matching
library(cobalt)
library(MatchIt)
library(DMwR)
library(broom)

dat_RECUR_ICD109_demographics$group[dat_RECUR_ICD109_demographics$group == 'case'] <- 1
dat_RECUR_ICD109_demographics$group[dat_RECUR_ICD109_demographics$group == 'control'] <- 0

table(dat_RECUR_ICD109_demographics$group)

dat_RECUR_ICD109_demographics$PT_SEX <- as.factor(dat_RECUR_ICD109_demographics$PT_SEX)
dat_RECUR_ICD109_demographics$group <- as.factor(dat_RECUR_ICD109_demographics$group)
dat_RECUR_ICD109_demographics$PT_RACE <- as.factor(dat_RECUR_ICD109_demographics$PT_RACE)
str(dat_RECUR_ICD109_demographics)


#optional for p123 matching and 90K sample only
load("FAM_90K_FAM_85K.RData")
ID_90K_p123 <- dat_RECUR_ICD109_demographics[dat_RECUR_ICD109_demographics$PT_ID %in% FAM_90K$V1_2, ]$PT_ID
ID_85K_p123 <- dat_RECUR_ICD109_demographics[!(dat_RECUR_ICD109_demographics$PT_ID %in% FAM_90K$V1_2), ]$PT_ID
length(ID_90K_p123)

dat_RECUR_ICD109_demographics_90K <- dat_RECUR_ICD109_demographics[dat_RECUR_ICD109_demographics$PT_ID %in% ID_90K_p123, ]
dat_RECUR_ICD109_demographics_85K <- dat_RECUR_ICD109_demographics[dat_RECUR_ICD109_demographics$PT_ID %in% ID_85K_p123, ]

set.seed(1234)
m.out <- matchit(group ~ Index_age + PT_SEX + PT_RACE, data = dat_RECUR_ICD109_demographics_90K, method="nearest", ratio=2) #ratio change from 1 to 5 or 1 to 10

m.out
a <- summary(m.out)
a

bal.tab(m.out)

tiff(paste('loveplot_matching', 'discovery1', '.tiff', sep=''), units="in", width=14, height=10, res=300)
love.plot(m.out, binary = "std", thresholds = c(m = .1))
dev.off()

summary (m.out,standardize = T)
dim(m.out$match.matrix) 
match.matrix<-m.out[["match.matrix"]]
index <- unique(as.integer(as.vector(match.matrix)))

#for discovery dataset(90K)
dat_RECUR_ICD109_demographics_90K.control <- dat_RECUR_ICD109_demographics_90K[rownames(dat_RECUR_ICD109_demographics_90K) %in% index, ]
dat_RECUR_ICD109_demographics_90K.case <- dat_RECUR_ICD109_demographics_90K[dat_RECUR_ICD109_demographics_90K$group == 1, ]
dat_RECUR_ICD109_demographics_90K <- rbind(dat_RECUR_ICD109_demographics_90K.case, dat_RECUR_ICD109_demographics_90K.control)
head(dat_RECUR_ICD109_demographics_90K)
save(dat_RECUR_ICD109_demographics_90K, file = "dat_RECUR_ICD109_demographics_90K_ratio2.RData", version = 2)

#same way to create replication dataset(85K)

#check the mean difference in index_age between carrier and noncarrier group
dat_RECUR_ICD109_demographics_90K %>%
  dplyr::group_by(group) %>%
  dplyr::summarise(mean = mean(Index_age), SD = sd(Index_age)) %>% 
  ungroup()

# # A tibble: 2 x 3
# group  mean    SD
# <fct> <dbl> <dbl>
#   1 0      56.8  18.2
# 2 1      56.7  18.4

est = dat_RECUR_ICD109_demographics_90K %>% do(tidy(t.test(Index_age ~ group, data = .)))
est
# # A tibble: 1 x 10
# estimate estimate1 estimate2 statistic p.value parameter conf.low conf.high method                  alternative
# <dbl>     <dbl>     <dbl>     <dbl>   <dbl>     <dbl>    <dbl>     <dbl> <chr>                   <chr>      
#   1    0.166      56.8      56.7     0.378   0.705     5146.   -0.696      1.03 Welch Two Sample t-test two.sided  

##########################################################################################################################
##PCA analysis

load("dat_RECUR_ICD109_demographics_90K_ratio2.RData")
head(dat_RECUR_ICD109_demographics_90K)
PCA <- read.table("pca_discovery.eigenvec", header = F, stringsAsFactors = F)
head(PCA) #7812*12
dim(PCA)
PCA$cohort <- "discovery"
PCA <- cSplit(PCA, "V1", "_")
colnames(PCA)[which(names(PCA) == "V1_2")] <- "PT_ID"
dat_RECUR_ICD109_demographics_90K_PCA <- merge(dat_RECUR_ICD109_demographics_90K, PCA, by = "PT_ID")
FAM_90K_discovery_ratio2 <- data.frame(PCA$V2, PCA$V2)
save(dat_RECUR_ICD109_demographics_90K_PCA, file = "dat_RECUR_ICD109_demographics_90K_PCA.RData", version =2)


PCA <- read.table("discoveryreplication_ratio2_EXOME2_pca.eigenvec", header = F, stringsAsFactors = F)
nrow(PCA) #12726
PCA <- cSplit(PCA, "V1", "_")
colnames(PCA)[which(names(PCA) == "V1_2")] <- "PT_ID"
PCA$PT_ID <- as.character(PCA$PT_ID)
PCA$cohort <- data.frame(ifelse(PCA$PT_ID %in% dat_RECUR_ICD109_demographics_90K$PT_ID, "discovery", "replication"))
PCA <- PCA[order(PCA$PT_ID),]

dat_RECUR_ICD109_demographics_90K85K <- rbind(dat_RECUR_ICD109_demographics_90K, dat_RECUR_ICD109_demographics_85K)
dat_RECUR_ICD109_demographics_90K85K <- dat_RECUR_ICD109_demographics_90K85K[order(dat_RECUR_ICD109_demographics_90K85K$PT_ID), ]
PCA_final <- cbind(PCA, dat_RECUR_ICD109_demographics_90K85K[, c("PT_SEX", "PT_RACE", "Index_age", "group")])
PCA_final$PT_RACE <- as.character(PCA_final$PT_RACE)
str(PCA_final)

tiff(paste('PCA_discovery_replication_06062022', '.tiff', sep=''), units="in", width=14, height=10, res=300)
ggplot(PCA_final, aes(x=V3, y=V4, color=cohort, size=Index_age, alpha = 0.8)) +
  geom_point() + 
  theme(text = element_text(size = 12)) + 
  ylab("PC2") +
  xlab("PC1")
dev.off()

########################################################################################################
