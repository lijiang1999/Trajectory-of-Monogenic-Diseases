setwd("D:/Geisinger/X19981_backup_10072020/Desktop/MS_project")
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

getwd()
ls(pattern= "O", all.names = TRUE)
ls()

#spike-in control
Genotype[Genotype$X19.15174391.G.A_G == 1,]$FID #this SNP should be removed
# [1] "GHS-60K-IDT_PT1117592_180757625" "GHS-60K-IDT_PT1202698_180739665" "GHS-60K-IDT_PT678791_180720974" 
# [4] "GHS-60K-IDT_PT726999_185878107"  "GHS-60K-IDT_PT865829_195995451"  "GHS-HSH_PT1360022_311967805"    
# [7] "GHS_PT1154976_301385422"         "GHS_PT121012_243902320"          "GHS_PT1346346_233626944"        
# [10] "GHS_PT242673_246252547"          "GHS_PT259508_246175955"          NA                               
# [13] "GHS_PT377514_231581910"          "GHS_PT527399_221629809"          "GHS_PT635221_214323031"         
# [16] "GHS_PT753905_301387964"          "GHS_PT865414_209999228"          "GHS_PT928617_311877320"

spikein_ID <- c(1117592, 1202698, 678791, 726999, 865829, 1360022, 1154976, 121012, 1346346, 242673, 259508,
                377514, 527399, 635221, 753905, 865414, 928617)
phenotypes_long_spikein_discovery <- phenotypes_long[phenotypes_long$id %in% spikein_ID, ]
phenotypes_long_spikein_discovery$cohort <- "discovery"

phenotypes_long_spikein_replication <- phenotypes_long[phenotypes_long$id %in% spikein_ID, ]
phenotypes_long_spikein_replication$cohort <- "replication"

phenotypes_long_spikein_discoveryreplication <- rbind(phenotypes_long_spikein_discovery, phenotypes_long_spikein_replication)

phenotypes_long_spikein_discoveryreplication_age <- rbind(phenotypes_long_spikein_discovery, phenotypes_long_spikein_replication)

save(phenotypes_long_spikein_discoveryreplication_age, file = "phenotypes_long_spikein_discoveryreplication_age.RData", version = 2)
load("phenotypes_long_spikein_discoveryreplication_age.RData")
save(phenotypes_long_spikein_discoveryreplication, file = "phenotypes_long_spikein_discoveryreplication.RData", version = 2)

phenotypes_long_spikein_discoveryreplication_top <- phenotypes_long_spikein_discoveryreplication[phenotypes_long_spikein_discoveryreplication$phecode == "X433", ]
phenotypes_long_spikein_discoveryreplication_top_summary <- phenotypes_long_spikein_discoveryreplication_top %>%
  group_by(cohort) %>%
  summarise(sum = sum(phenotype, na.rm = TRUE)) %>%
  ungroup
 
P <- fisher.test(phenotypes_long_spikein_discoveryreplication_top$phenotype, phenotypes_long_spikein_discoveryreplication_top$cohort)$p.value
S <- fisher.test(phenotypes_long_spikein_discoveryreplication_top$phenotype, phenotypes_long_spikein_discoveryreplication_top$cohort)$estimate #this is odds ratio
x <- c(P, S)
x <- t(data.frame(x))
colnames(x) <- c("p.value", "estimate")

# load("Summary_phenotypes109_merge_anno_commonrare_delta_p123_90K_phenotypes109_WHITE.RData")
# Summary_phenotypes_merge_anno_common$phecode
# Summary_phenotypes_merge_anno_common$phecode <- paste0("X", Summary_phenotypes_merge_anno_common$phecode)


candidate_phecode <- c("X286.7", "X286.8", "X286.81", "X395", "X396", "X401", "X411.2", "X411.4", "X415", "X425",
    "X425.1", "X426", "X427.2", "X430", "X433", "X433.1", "X433.11", "X433.12", "X433.2", "X433.21", "X433.3", "X433.31",
    "X433.5", "X433.6", "X433.8", "X440", "X443.9", "X452.2")

#only the following works
candidate_phecode <- c("X395", "X396", "X401", "X415", "X425",
                       "X425.1", "X427.2", "X433", "X433.1",
                       "X452.2")

phenotypes_long_spikein_discoveryreplication <- phenotypes_long_spikein_discoveryreplication_age 
#fisher exact test
##make a function to get all for top 50 using fisherexact test
fisherexact <- NULL
P <- NULL
S <- NULL

for (phecode in candidate_phecode){  # the top 50 file
  cat(phecode)
  phenotypes_long_spikein_discoveryreplication_top <- phenotypes_long_spikein_discoveryreplication[phenotypes_long_spikein_discoveryreplication$phecode == phecode, ]
  phenotypes_long_spikein_discoveryreplication_top_summary <- phenotypes_long_spikein_discoveryreplication_top %>%
    group_by(cohort) %>%
    summarise(sum = sum(phenotype, na.rm = TRUE)) %>%
    ungroup
  
  P <- fisher.test(phenotypes_long_spikein_discoveryreplication_top$phenotype, phenotypes_long_spikein_discoveryreplication_top$cohort)$p.value
  S <- fisher.test(phenotypes_long_spikein_discoveryreplication_top$phenotype, phenotypes_long_spikein_discoveryreplication_top$cohort)$estimate #this is odds ratio
  x <- c(P, S, phenotypes_long_spikein_discoveryreplication_top_summary$sum[1], phenotypes_long_spikein_discoveryreplication_top_summary$sum[2])
  x <- t(data.frame(x))
  colnames(x) <- c("p.value", "estimate", "Discovery", "Replication")
  rownames(x) <- phecode
  fisherexact <- rbind(fisherexact, x)
  write.csv(fisherexact, "Spikein_fisherexact_phenotypes109_age_ratio2_WHITE.csv")
  saveRDS(fisherexact, file= "Spikein_fisherexact_phenotypes109_age_ratio2_WHITE.rds")
}

fisherexact <- data.frame(readRDS("Spikein_fisherexact_phenotypes109_age_ratio2_WHITE.rds"))
dim(fisherexact)
fisherexact$phecode <- rownames(fisherexact)
head(fisherexact, 12)

#comparing discovery(n=7) and replication (n=10) for EUR w/o age cutoff
#          p.value  estimate Discovery Replication phecode
# X395   0.4117647 0.0000000         1           0    X395
# X396   1.0000000       Inf         0           1    X396
# X401   1.0000000 0.6823114         6           8    X401
# X415   1.0000000 0.6831461         1           1    X415
# X425   0.4117647 0.0000000         1           0    X425
# X425.1 0.4117647 0.0000000         1           0  X425.1
# X427.2 1.0000000 1.4656064         1           2  X427.2
# X433   1.0000000 1.4656064         1           2    X433
# X433.1 1.0000000 1.4656064         1           2  X433.1
# X452.2 0.4117647 0.0000000         1           0  X452.2

#comparing discovery(n=7) and replication (n=10) for EUR w/ age cutoff
#          p.value  estimate Discovery Replication phecode
# X395   0.4117647 0.0000000         1           0    X395
# X396   1.0000000       Inf         0           1    X396
# X401   0.6029412 0.4099416         6           7    X401
# X415   1.0000000 0.6831461         1           1    X415
# X425   0.4117647 0.0000000         1           0    X425
# X425.1 0.4117647 0.0000000         1           0  X425.1
# X427.2 1.0000000 1.4656064         1           2  X427.2
# X433   1.0000000 0.6831461         1           1    X433
# X433.1 1.0000000 0.6831461         1           1  X433.1
# X452.2 0.4117647 0.0000000         1           0  X452.2

discovery <- dat_RECUR_casecontrol_select_combined_90K[dat_RECUR_casecontrol_select_combined_90K$id %in% spikein_ID, ]
mean(discovery$Index_age)
sd(discovery$Index_age)

replication <- dat_RECUR_casecontrol_select_combined_85K[dat_RECUR_casecontrol_select_combined_85K$id %in% spikein_ID, ]
mean(replication$Index_age)
sd(replication$Index_age)

dat_RECUR_casecontrol_select_combined_90K$PT_ID
####################################################################################################
#04102022 the following is required to start.
load("dat_RECUR_ICD109.RData")
dat_RECUR_ICD109$phecode <- paste0("X", dat_RECUR_ICD109$phecode)
str(dat_RECUR_ICD109)

dat_RECUR_ICD109 <- dat_RECUR_ICD109[dat_RECUR_ICD109$PT_RACE == "White", ]
#dat_RECUR_ICD109 <- dat_RECUR_ICD109[dat_RECUR_ICD109$PT_RACE == "Black Or African American", ]

load("Summary_phenotypes109_merge_anno_commonrare_delta_p123_90K_phenotypes109_age.RData")
load("Summary_phenotypes109_merge_anno_commonrare_delta_p123_85K_phenotypes109_age.RData")

load("Summary_phenotypes109_merge_anno_commonrare_delta_p123_90K_phenotypes109.RData")
load("Summary_phenotypes109_merge_anno_commonrare_delta_p123_85K_phenotypes109.RData")


load("Summary_phenotypes109_merge_anno_commonrare_delta_p123_90K_phenotypes109_age_WHITE.RData")
load("Summary_phenotypes109_merge_anno_commonrare_delta_p123_85K_phenotypes109_age_WHITE.RData")
load("Summary_phenotypes109_merge_anno_commonrare_delta_p123_90K_phenotypes109_age_BLACK.RData")
load("Summary_phenotypes109_merge_anno_commonrare_delta_p123_85K_phenotypes109_age_BLACK.RData")

load("Summary_phenotypes109_merge_anno_commonrare_delta_p123_90K_phenotypes109_WHITE.RData")
load("Summary_phenotypes109_merge_anno_commonrare_delta_p123_85K_phenotypes109_WHITE.RData")
head(Summary_phenotypes_merge_anno_common)
table_summary_90K <- data.frame(table(dat_RECUR_casecontrol_select_combined_90K$PT_SEX, dat_RECUR_casecontrol_select_combined_90K$PT_RACE, dat_RECUR_casecontrol_select_combined_90K$group))
colnames(table_summary_90K) <- c("SEX", "RACE", "GROUP", "FREQ")
table_summary_90K$COHORT <- "Discovery"

table_summary_85K <- data.frame(table(dat_RECUR_casecontrol_select_combined_85K$PT_SEX, dat_RECUR_casecontrol_select_combined_85K$PT_RACE, dat_RECUR_casecontrol_select_combined_85K$group))
colnames(table_summary_85K) <- c("SEX", "RACE", "GROUP", "FREQ")
table_summary_85K$COHORT <- "Replication"
table_summary <- cbind(table_summary_90K, table_summary_85K)
save(table_summary, file = "table_summary_demo.RData", version = 2)
write.table(table_summary, file = "table_summary_demo.txt", sep = "|", col.names = F, row.names = F, quote = F)


load("Summary_phenotypes_merge_anno_common_enriched_parameters_65yrs_90K_WHITE.RData")
write.csv(Summary_phenotypes_merge_anno_common_enriched_parameters, file = "Summary_phenotypes_merge_anno_common_enriched_parameters_65yrs_90K_WHITE.csv", col.names = F, row.names = F, quote = F)
#load("phenotypes109_age.RData")
load("phenotypes109.RData")
load("12092021_SKAT_inputfile_p123_90K_updated_phenotypes109.RData") #to get phenotypes_ratio2 file for discovery
load("12092021_SKAT_inputfile_p123_85K_updated_phenotypes109.RData") #to get phenotypes_ratio2 file for replication
load("12092021_SKAT_inputfile_p123_90K_updated_phenotypes109_age.RData") #to get phenotypes_ratio2 file for discovery
load("12092021_SKAT_inputfile_p123_85K_updated_phenotypes109_age.RData") #to get phenotypes_ratio2 file for replication
#w/ age cutoff
phenotypes109_ratio2 <- phenotypes109_age[phenotypes109_age$id %in% phenotypes_ratio2$id, ]
#w/o age cutoff
phenotypes109_ratio2 <- phenotypes109[phenotypes109$id %in% phenotypes_ratio2$id, ] 

phenotypes109_ratio2$id
table(phenotypes109_ratio2$group)
# case control 
# 2745    5490
# case control 
# 1705    3410
phenotypes_ratio2 <- phenotypes109_ratio2 #from discovery or replication
colnames(phenotypes_ratio2)
gathercols <- colnames(phenotypes_ratio2[, 2:1858]) 
keycol <- "phecode"
valuecol <- "phenotype"
phenotypes_long <- data.frame(gather_(phenotypes_ratio2, keycol, valuecol, gathercols, factor_key=TRUE))
head(phenotypes_long)
# id   group phecode phenotype
# 1  64 control     008         0
# 2 185 control     008         0
# 3 221    case     008         0
# 4 319 control     008         0
# 5 467 control     008         0
# 6 517 control     008         0

phenotypes_long$phecode <- paste0("X", phenotypes_long$phecode)
phecode <- unique(phenotypes_long$phecode)
phecode
unique(Summary_phenotypes_merge_anno_common_enriched$exclude_name)
# [1] "infectious diseases"     "neoplasms"               "endocrine & metabolic"   "hematopoietic"          
# [5] "mental disorders"        "neurological"            "sense organs"            "circulatory system"     
# [9] "respiratory"             "digestive"               "genitourinary"           "pregnancy complications"
# [13] "dermatologic"            "musculoskeletal"         "congenital anomalies"    "symptoms"               
# [17] "injuries & poisonings" 
Summary_phenotypes_merge_anno_common_enriched$phecode

fooCox = list()
#KM analysis
for (i in unique(Summary_phenotypes_merge_anno_common_enriched$exclude_name)) {
  print(i)
  Summary_phenotypes_merge_anno_common_enriched_select <- Summary_phenotypes_merge_anno_common_enriched[Summary_phenotypes_merge_anno_common_enriched$exclude_name == i, ]
  phecode_select <- as.character(paste0("X", Summary_phenotypes_merge_anno_common_enriched_select$phecode))
  print(phecode_select)
  
  phenotypes_long_AA <- phenotypes_long %>% #change to relatedness or nonrelatedness files
    #dplyr::filter(phecode == "X433.3") %>% #Cerebral ischemia for X433; X433.3; #Cardiomyopathy for X425; X425.1
    dplyr::filter(phecode %in% phecode_select) %>% #select of enriched phecode set
    dplyr::filter(phenotype == 1)
  AA_ID <- unique(paste0("PT", "", phenotypes_long_AA$id))
  
  phenotypes_long_AA <- phenotypes_long %>% #change to relatedness or nonrelatedness files
    #dplyr::filter(phecode == "X433.3") %>% #Cerebral ischemia for X433; #Cardiomyopathy for X425
    dplyr::filter(phecode %in% phecode_select) %>% #select of enriched phecode set
    group_by(id) %>%
    dplyr::filter(sum(phenotype) == 0) %>%
    ungroup()
  AA_WO <- unique(paste0("PT", "", phenotypes_long_AA$id))
  
  dat_RECUR_ICD109$disease <- ifelse(dat_RECUR_ICD109$PT_ID %in% AA_ID, 1, 0)
  dat_RECUR_ICD109_select <- dat_RECUR_ICD109[dat_RECUR_ICD109$PT_ID %in% AA_ID, ] %>%
    #dplyr::filter(phecode == "X433.3") %>%
    dplyr::filter(phecode %in% phecode_select) %>% #select of enriched phecode set
    group_by(PT_ID) %>%
    arrange(desc(ENC_DT)) %>%
    slice(n()) %>%
    ungroup()
  dat_RECUR_ICD109_select <- data.frame(dat_RECUR_ICD109_select)
  print(table(dat_RECUR_ICD109_select$group))
  
  dat_RECUR_ICD109_WO <- dat_RECUR_ICD109[dat_RECUR_ICD109$PT_ID %in% AA_WO, ] %>%
    dplyr::group_by(PT_ID) %>%
    dplyr::arrange(ENC_DT) %>%
    slice(n()) %>%
    #mutate(disease = 0) %>%
    ungroup()
  dat_RECUR_ICD109_WO <- data.frame(dat_RECUR_ICD109_WO)
  GNSIS <- rbind(dat_RECUR_ICD109_select, dat_RECUR_ICD109_WO)
  
  GNSIS <- GNSIS %>% 
    mutate(timediff_recur_censor = 
             case_when(disease == 1 ~ (timediff_occur_censor - 0)*365.25,
                       disease == 0 ~ (timediff_occur_censor - 0)*365.25)) 
  GNSIS <- GNSIS %>% 
    mutate(timediff_recur_censor_15yr =
             case_when(timediff_recur_censor < 15*365.25 ~ timediff_recur_censor,
                       timediff_recur_censor >= 15*365.25 ~ 15*365.25),
           event_recur_censor_15yr =
             case_when(disease == 0 ~ 0,
                       disease == 1 & timediff_recur_censor <= 15*365.25 ~1,
                       disease == 1 & timediff_recur_censor > 15*365.25 ~0),
           timediff_recur_censor_25yr =
             case_when(timediff_recur_censor < 25*365.25 ~ timediff_recur_censor,
                       timediff_recur_censor >= 25*365.25 ~ 25*365.25),
           event_recur_censor_25yr =
             case_when(disease == 0 ~ 0,
                       disease == 1 & timediff_recur_censor <= 25*365.25 ~1,
                       disease == 1 & timediff_recur_censor > 25*365.25 ~0),
           timediff_recur_censor_35yr =
             case_when(timediff_recur_censor < 35*365.25 ~ timediff_recur_censor,
                       timediff_recur_censor >= 35*365.25 ~ 35*365.25),
           event_recur_censor_35yr =
             case_when(disease == 0 ~ 0,
                       disease == 1 & timediff_recur_censor <= 35*365.25 ~1,
                       disease == 1 & timediff_recur_censor > 35*365.25 ~0),
           timediff_recur_censor_65yr =
             case_when(timediff_recur_censor < 65*365.25 ~ timediff_recur_censor,
                       timediff_recur_censor >= 65*365.25 ~ 65*365.25),
           event_recur_censor_65yr =
             case_when(disease == 0 ~ 0,
                       disease == 1 & timediff_recur_censor <= 65*365.25 ~1,
                       disease == 1 & timediff_recur_censor > 65*365.25 ~0)
    )
  
  
  #optional
  GNSIS_select <- GNSIS %>%
    filter(timediff_recur_censor_65yr >= 0)
  #GNSIS_select <- GNSIS
  # KM stratified by cohort
  #ANN <- phecode_mapping_unique[phecode_mapping_unique$phecodeX %in% phecode_select, ]$phecode_str
  scale(GNSIS_select$timediff_occur_censor)
  #cox_65yrrecur_AA <- coxph(Surv(timediff_recur_censor_65yr/365.25, event_recur_censor_65yr)~as.factor(group) + as.factor(PT_SEX) + as.factor(PT_RACE), data = GNSIS_select)
  cox_65yrrecur_AA <- coxph(Surv(timediff_recur_censor_65yr/365.25, event_recur_censor_65yr)~as.factor(group) + as.factor(PT_SEX), data = GNSIS_select)
  coxtidy_65yrrecur_AA <- tidy(cox_65yrrecur_AA, exponentiate = T, conf.int = T, conf.level = 0.95)
  coxtidy_65yrrecur_AA$Category <- i
  coxtidy_65yrrecur_AA <- data.frame(coxtidy_65yrrecur_AA)
  fooCox[[i]] <- coxtidy_65yrrecur_AA
}
result_final_Cox <- do.call(rbind, fooCox)
result_final_Cox$estimate_95CI <- paste0(round(1/result_final_Cox$estimate, 3), "(", round(1/result_final_Cox$conf.high, 3), "-", round(1/result_final_Cox$conf.low, 3), ")")

#save(result_final_Cox, file = "result_final_Cox_90K_65yr_BLACK.RData", version =2)
#save(result_final_Cox, file = "result_final_Cox_85K_65yr_BLACK.RData", version =2)
# save(result_final_Cox, file = "result_final_Cox_85K_woRACE_65yr.RData", version =2)
# save(result_final_Cox, file = "result_final_Cox_85K_90Kphecodes_65yr.RData", version =2)
save(result_final_Cox, file = "result_final_Cox_85K_90Kphecodes_BLACK_65yr.RData", version =2)
#save(result_final_Cox, file = "result_final_Cox_90K_85Kphecodes_65yr_BLACK.RData", version =2)
# save(result_final_Cox, file = "result_final_Cox_90K_85Kphecodes_woRACE_65yr.RData", version =2)
# save(result_final_Cox, file = "result_final_Cox_90K_woRACE_65yr.RData", version =2)


result_final_Cox_select <- result_final_Cox %>%
  filter(term == 'as.factor(group)control')
result_final_Cox_select$p.value <- round(result_final_Cox_select$p.value, 3)

colnames(result_final_Cox_select)
# [1] "term"          "estimate"      "std.error"     "statistic"     "p.value"       "conf.low"      "conf.high"    
# [8] "Category"      "estimate_95CI"

#create plot without pregnancy complication.
result_final_Cox_select <- result_final_Cox_select[result_final_Cox_select$Category != "pregnancy complications", ]
#create a forest plot
tabletext <- result_final_Cox_select %>%
  dplyr::select(Category, p.value, estimate_95CI)

head(tabletext)

header <- c("Disease Category", "P value", "Hazard Ratio (95%CI)")
tabletext <- rbind(header, tabletext)
forestdata <- result_final_Cox_select %>%
  dplyr::select(estimate, conf.high, conf.low,)
forestdata$estimate <- 1/forestdata$estimate
forestdata$conf.low <- 1/forestdata$conf.low
forestdata$conf.high <- 1/forestdata$conf.high
head(forestdata)
header2 <- c(1, 1, 1)
forestdata <- rbind(header2, forestdata)
forestdata_new <- forestdata[c(1:3, 4:18),] 
library(forestplot)
tabletext_new <- tabletext[c(1:3, 4:18),] 

head(tabletext_new)

#tiff("Rplot_forestplot_multivariate_coxregression_disease_category_85K_90Kphecodes_woRACE_65yr_wopregnancy.tiff", units = "in", width = 10, height = 8, res = 300)
tiff("Rplot_forestplot_multivariate_coxregression_disease_category_85K_90Kphecodes_BLACK_65yr_wopregnancy.tiff", units = "in", width = 10, height = 8, res = 300)

p <- forestplot(tabletext_new, 
                graph.pos = 3,
                forestdata_new,new_page = TRUE,
                #is.summary=c(TRUE,TRUE,rep(FALSE,8),TRUE),
                clip=c(-2,6),
                boxsize=0.25,
                zero = 1,
                xlog=FALSE, 
                col=fpColors(box="royalblue",line="darkblue", summary="royalblue"),
                xlab = "estimate(95%CI)",
                txt_gp = fpTxtGp(label = list(gpar(fontfamily = "Arial"),
                                              gpar(fontfamily = "",
                                                   col = "#660000")),
                                 ticks = gpar(fontfamily = "", cex = 1.0),
                                 xlab  = gpar(fontfamily = "Arial", cex = 1.0)))

grid::grid.text("Disease category ~ having rare variants + sex + PCs(1-5)", .5, 0.99, gp=gpar(cex=1.5))
dev.off()

################################################################################

#select phecodes from a candidate list
phecode_select <- c("X286.7", "X286.8", "X286.81")
phecode_select <- c("X286.8", "X286.81")
phecode_select <- c("X425", "X425.1")
phecode_select <- c("X425.1")
phecode_select <- c("X433.1", "X433.11", "X433.12", "X433.2", "X433.21", "X433.3", "X433.31", "X433.32", "X433.6", "X433.8")
phecode_select <- c("X433.2", "X433.21", "X433.3", "X433.31", "X433.32", "X433.6", "X433.8")
phecode_select <- c("X433.2", "X433.21", "X433.3", "X433.31", "X433.32", "X433.8")

#get
phenotypes_long_AA <- phenotypes_long %>% 
  dplyr::filter(phecode %in% phecode_select) %>% #select of enriched phecode set
  dplyr::filter(phenotype == 1)
AA_ID <- unique(paste0("PT", "", phenotypes_long_AA$id))

phenotypes_long_AA <- phenotypes_long %>% #change to relatedness or nonrelatedness files
  #dplyr::filter(phecode == "X433.3") %>% #Cerebral ischemia for X433; #Cardiomyopathy for X425
  dplyr::filter(phecode %in% phecode_select) %>% #select of enriched phecode set
  group_by(id) %>%
  dplyr::filter(sum(phenotype) == 0) %>%
  ungroup()
AA_WO <- unique(paste0("PT", "", phenotypes_long_AA$id))

dat_RECUR_ICD109$disease <- ifelse(dat_RECUR_ICD109$PT_ID %in% AA_ID, 1, 0)
dat_RECUR_ICD109_select <- dat_RECUR_ICD109[dat_RECUR_ICD109$PT_ID %in% AA_ID, ] %>%
  #dplyr::filter(phecode == "X433.3") %>%
  dplyr::filter(phecode %in% phecode_select) %>% #select of enriched phecode set
  group_by(PT_ID) %>%
  arrange(desc(ENC_DT)) %>%
  slice(n()) %>%
  ungroup()
dat_RECUR_ICD109_select <- data.frame(dat_RECUR_ICD109_select)
print(table(dat_RECUR_ICD109_select$group))

dat_RECUR_ICD109_WO <- dat_RECUR_ICD109[dat_RECUR_ICD109$PT_ID %in% AA_WO, ] %>%
  dplyr::group_by(PT_ID) %>%
  dplyr::arrange(ENC_DT) %>%
  slice(n()) %>%
  #mutate(disease = 0) %>%
  ungroup()
dat_RECUR_ICD109_WO <- data.frame(dat_RECUR_ICD109_WO)
print(table(dat_RECUR_ICD109_WO$group))

GNSIS <- rbind(dat_RECUR_ICD109_select, dat_RECUR_ICD109_WO)

GNSIS <- GNSIS %>% 
  mutate(timediff_recur_censor = 
           case_when(disease == 1 ~ (timediff_occur_censor - 0)*365.25,
                     disease == 0 ~ (timediff_occur_censor - 0)*365.25)) 
GNSIS <- GNSIS %>% 
  mutate(timediff_recur_censor_55yr =
           case_when(timediff_recur_censor < 55*365.25 ~ timediff_recur_censor,
                     timediff_recur_censor >= 55*365.25 ~ 55*365.25),
         event_recur_censor_55yr =
           case_when(disease == 0 ~ 0,
                     disease == 1 & timediff_recur_censor <= 55*365.25 ~1,
                     disease == 1 & timediff_recur_censor > 55*365.25 ~0),
         timediff_recur_censor_60yr =
           case_when(timediff_recur_censor < 60*365.25 ~ timediff_recur_censor,
                     timediff_recur_censor >= 60*365.25 ~ 60*365.25),
         event_recur_censor_60yr =
           case_when(disease == 0 ~ 0,
                     disease == 1 & timediff_recur_censor <= 60*365.25 ~1,
                     disease == 1 & timediff_recur_censor > 60*365.25 ~0),
         timediff_recur_censor_85yr =
           case_when(timediff_recur_censor < 85*365.25 ~ timediff_recur_censor,
                     timediff_recur_censor >= 85*365.25 ~ 85*365.25),
         event_recur_censor_85yr =
           case_when(disease == 0 ~ 0,
                     disease == 1 & timediff_recur_censor <= 85*365.25 ~1,
                     disease == 1 & timediff_recur_censor > 85*365.25 ~0),
         timediff_recur_censor_65yr =
           case_when(timediff_recur_censor < 65*365.25 ~ timediff_recur_censor,
                     timediff_recur_censor >= 65*365.25 ~ 65*365.25),
         event_recur_censor_65yr =
           case_when(disease == 0 ~ 0,
                     disease == 1 & timediff_recur_censor <= 65*365.25 ~1,
                     disease == 1 & timediff_recur_censor > 65*365.25 ~0)
  )
GNSIS_select <- GNSIS

#create conditional covariates
GNSIS_condition <- GNSIS
GNSIS_condition$disease -> GNSIS_condition$hypercoagulability

GNSIS <- merge(GNSIS_select, GNSIS_condition[, c("PT_ID", "hypercoagulability")], by = "PT_ID")



#specify the patients with gene of mutation
load("Individualgene_PTID_12032021.RData")
listOfID <- list(COL4A1, NOTCH3, TREX1, GLA, HTRA1, CTC1)
x = 2
nrow(GNSIS_select[GNSIS_select$PT_ID %in% listOfID[[x]], ])

Intersect_ID <- intersect(intersect(COL4A1,NOTCH3),TREX1) #0
Intersect_ID_COL4A1TREX1 <- intersect(COL4A1, TREX1) #2
#[1] "PT149790"  "PT1657533"
Intersect_ID_COL4A1NOTCH3 <- intersect(COL4A1, NOTCH3) #12
# [1] "PT1437014" "PT104150"  "PT119364"  "PT1340859" "PT377762"  "PT626645"  "PT712165"  
# "PT865414"  "PT925831" 
# [10] "PT122220"  "PT1782039" "PT1154655"
Intersect_ID_COL4A1HTRA1 <- intersect(COL4A1, HTRA1) #0
Intersect_ID_COL4A1CTC1 <- intersect(COL4A1, CTC1) #0
Intersect_ID_COL4A1GLA <- intersect(COL4A1, GLA) #5
#[1] "PT472491"  "PT535386"  "PT1007297" "PT369310"  "PT479128"
Intersect_ID_NOTCH3TREX1 <- intersect(NOTCH3, TREX1) #2
#[1] "PT483269" "PT840406"
Intersect_ID_NOTCH3HTRA1 <- intersect(NOTCH3, HTRA1) #1
#[1] "PT1057017"
Intersect_ID_NOTCH3CTC1 <- intersect(NOTCH3, CTC1) #0
Intersect_ID_NOTCH3GLA <- intersect(NOTCH3, GLA) #0
#[1] "PT572316" "PT839522"
Intersect_ID_TREX1GLA <- intersect(TREX1, GLA) #0
Intersect_ID_TREX1HTRA1 <- intersect(TREX1, HTRA1) #0
Intersect_ID_TREX1CTC1 <- intersect(TREX1, CTC1) #0
Intersect_ID_GLAHTRA1 <- intersect(GLA, HTRA1) #0
Intersect_ID_GLACTC1 <- intersect(GLA, CTC1) #0
Intersect_ID_HTRA1CTC1 <- intersect(HTRA1, CTC1) #0


x <- list(Intersect_ID_COL4A1TREX1, Intersect_ID_COL4A1NOTCH3, Intersect_ID_COL4A1HTRA1, 
          Intersect_ID_COL4A1CTC1, Intersect_ID_COL4A1GLA, Intersect_ID_NOTCH3TREX1, 
          Intersect_ID_NOTCH3HTRA1, Intersect_ID_NOTCH3CTC1, Intersect_ID_NOTCH3GLA,
          Intersect_ID_TREX1GLA, Intersect_ID_TREX1HTRA1, Intersect_ID_TREX1CTC1, 
          Intersect_ID_GLAHTRA1, Intersect_ID_GLACTC1, Intersect_ID_HTRA1CTC1)

x <- list(Intersect_ID_COL4A1TREX1, Intersect_ID_COL4A1NOTCH3, Intersect_ID_COL4A1HTRA1, 
          Intersect_ID_COL4A1GLA, Intersect_ID_NOTCH3TREX1, 
          Intersect_ID_NOTCH3HTRA1, Intersect_ID_NOTCH3GLA,
          Intersect_ID_TREX1GLA, Intersect_ID_TREX1HTRA1, 
          Intersect_ID_GLAHTRA1)
x_list <- Reduce(c, x) #24 all independent
save(x_list, file = "x_list.RData", version =2)
x_list

library(venn)
venn(5, ilab=TRUE, zcolor = "style")
install.packages("ggpolypath")   # Install ggpolypath package
library("ggpolypath")            # Load ggpolypath package
#Next, we can draw a ggplot2 venn diagram by setting the ggplot2 argument within the venn function to be equal to TRUE:
venn(5, ggplot = TRUE)   
dev.off()

library("ggVennDiagram")
library(ggplot2)
p <- ggVennDiagram(x, label_alpha = 0)
p + scale_fill_distiller(palette = "Reds", direction = 1) +
  labs(title = "Pipeline 1, 2, and 3",
       subtitle = "Patients with rare variants from Seven Genes in 90K ")

#listOfID[[1]]
GNSIS_select <- GNSIS
#GNSIS_select <- GNSIS_select[GNSIS_select$PT_RACE == "Black Or African American", ]

GNSIS_select$subgroup <- "W/O Mut"
GNSIS_select[GNSIS_select$PT_ID %in% COL4A1, ]$subgroup <- "COL4A1"
GNSIS_select[GNSIS_select$PT_ID %in% TREX1, ]$subgroup <- "TREX1"
GNSIS_select[GNSIS_select$PT_ID %in% HTRA1, ]$subgroup <- "HTRA1"
GNSIS_select[GNSIS_select$PT_ID %in% GLA, ]$subgroup <- "GLA"
GNSIS_select[GNSIS_select$PT_ID %in% CTC1, ]$subgroup <- "CTC1"
GNSIS_select[GNSIS_select$PT_ID %in% NOTCH3, ]$subgroup <- "NOTCH3"

table(GNSIS_select$group)
table(GNSIS_select$subgroup)
# COL4A1    CTC1     GLA   HTRA1  NOTCH3   TREX1 W/O Mut 
# 1698      35     174      58     585     195    5490 
# COL4A1    CTC1     GLA   HTRA1  NOTCH3   TREx1 W/O Mut 
# 1696      35     174      58     583     195    5483 (hematopoietic)

# COL4A1    CTC1     GLA   HTRA1  NOTCH3   TREX1 W/O Mut 
# 761      33     142      62     518     188    3410
# COL4A1    CTC1     GLA   HTRA1  NOTCH3   TREX1 W/O Mut 
# 761      33     142      62     518     188    3408 (85K using 90K enriched phecode set)
# COL4A1    CTC1     GLA   HTRA1  NOTCH3   TREX1 W/O Mut 
# 761      33     142      61     518     188    3400 (hematopoietic)

#white
# COL4A1    CTC1     GLA   HTRA1  NOTCH3   TREX1 W/O Mut 
# 1560      34     174      56     577     192    5241 
# 1560      34     173      56     580     190    5241 #correct
# COL4A1    CTC1     GLA   HTRA1  NOTCH3   TREX1 W/O Mut 
# 568      33     142      61     511     169    3036

#Black
# COL4A1    CTC1   HTRA1   NOTCH3 W/O Mut 
# 125       1       2       3     211 
# COL4A1     GLA  TREX1   NOTCH3 W/O Mut 
# 178       1       5      17     299 

c("COL4A1", "NOTCH3", "TREX1", "GLA", "HTRA1", "CTC1", "ALL")

listOfDataFrames <- NULL
x = 0

for (i in c("COL4A1", "NOTCH3", "TREX1", "GLA", "HTRA1", "CTC1", "ALL")) {
  print(i)
  GNSIS_select <- GNSIS
  GNSIS_select$subgroup <- "W/O Mut"
  
  if (i == "ALL"){
    GNSIS_select_COL4A1 <- GNSIS_select
  } else {
    x = x + 1
    GNSIS_select[GNSIS_select$PT_ID %in% listOfID[[x]], ]$subgroup <- i
    GNSIS_select_COL4A1 <- rbind(GNSIS_select[GNSIS_select$subgroup == i, ], GNSIS_select[GNSIS_select$group == "control", ]) 
  }
  print(table(GNSIS_select_COL4A1$group))
  
  #cox_65yrrecur_COL4A1 <- coxph(Surv(timediff_recur_censor_65yr/365.25, event_recur_censor_65yr)~as.factor(group) + as.factor(PT_SEX), data = GNSIS_select_COL4A1)
  #add conditional covariate
  cox_65yrrecur_COL4A1 <- coxph(Surv(timediff_recur_censor_65yr/365.25, event_recur_censor_65yr)~as.factor(group) + as.factor(PT_SEX) + as.factor(hypercoagulability), data = GNSIS_select_COL4A1)
  
  coxtidy_65yrrecur_COL4A1 <- tidy(cox_65yrrecur_COL4A1, exponentiate = T, conf.int = T, conf.level = 0.95)
  coxtidy_65yrrecur_COL4A1$Onset <- "Early"
  coxtidy_65yrrecur_COL4A1$Category <- "Cerebral ischemia"
  coxtidy_65yrrecur_COL4A1$Gene <- i
  coxtidy_65yrrecur_COL4A1$Sex <- "Combined"
  coxtidy_65yrrecur_COL4A1 <- data.frame(coxtidy_65yrrecur_COL4A1)
  
  #cox_85yrrecur_COL4A1 <- coxph(Surv(timediff_recur_censor_85yr/365.25, event_recur_censor_85yr)~as.factor(group) + as.factor(PT_SEX), data = GNSIS_select_COL4A1)
  #add conditional covariate
  cox_85yrrecur_COL4A1 <- coxph(Surv(timediff_recur_censor_85yr/365.25, event_recur_censor_85yr)~as.factor(group) + as.factor(PT_SEX) + as.factor(hypercoagulability), data = GNSIS_select_COL4A1)
  
  coxtidy_85yrrecur_COL4A1 <- tidy(cox_85yrrecur_COL4A1, exponentiate = T, conf.int = T, conf.level = 0.95)
  coxtidy_85yrrecur_COL4A1$Onset <- "Lifetime"
  coxtidy_85yrrecur_COL4A1$Category <- "Cerebral ischemia"
  coxtidy_85yrrecur_COL4A1$Gene <- i
  coxtidy_85yrrecur_COL4A1$Sex <- "Combined"
  coxtidy_85yrrecur_COL4A1 <- data.frame(coxtidy_85yrrecur_COL4A1)
  
  fooCox_COL4A1_90K_WHITE <- rbind(coxtidy_65yrrecur_COL4A1, coxtidy_85yrrecur_COL4A1)
  
  GNSIS_select_COL4A1_FEMALE <- GNSIS_select_COL4A1[GNSIS_select_COL4A1$PT_SEX == "Female", ]
  
  #cox_65yrrecur_COL4A1 <- coxph(Surv(timediff_recur_censor_65yr/365.25, event_recur_censor_65yr)~as.factor(group), data = GNSIS_select_COL4A1_FEMALE)
  #add conditional covariate
  cox_65yrrecur_COL4A1 <- coxph(Surv(timediff_recur_censor_65yr/365.25, event_recur_censor_65yr)~as.factor(group) + as.factor(hypercoagulability), data = GNSIS_select_COL4A1_FEMALE)
  
  coxtidy_65yrrecur_COL4A1 <- tidy(cox_65yrrecur_COL4A1, exponentiate = T, conf.int = T, conf.level = 0.95)
  coxtidy_65yrrecur_COL4A1$Onset <- "Early"
  coxtidy_65yrrecur_COL4A1$Category <- "Cerebral ischemia"
  coxtidy_65yrrecur_COL4A1$Gene <- i
  coxtidy_65yrrecur_COL4A1$Sex <- "Female"
  coxtidy_65yrrecur_COL4A1 <- data.frame(coxtidy_65yrrecur_COL4A1)
  
  
  #cox_85yrrecur_COL4A1 <- coxph(Surv(timediff_recur_censor_85yr/365.25, event_recur_censor_85yr)~as.factor(group), data = GNSIS_select_COL4A1_FEMALE)
  #add conditional covariate
  cox_85yrrecur_COL4A1 <- coxph(Surv(timediff_recur_censor_85yr/365.25, event_recur_censor_85yr)~as.factor(group) + as.factor(hypercoagulability), data = GNSIS_select_COL4A1_FEMALE)
  coxtidy_85yrrecur_COL4A1 <- tidy(cox_85yrrecur_COL4A1, exponentiate = T, conf.int = T, conf.level = 0.95)
  coxtidy_85yrrecur_COL4A1$Onset <- "Lifetime"
  coxtidy_85yrrecur_COL4A1$Category <- "Cerebral ischemia"
  coxtidy_85yrrecur_COL4A1$Gene <- i
  coxtidy_85yrrecur_COL4A1$Sex <- "Female"
  coxtidy_85yrrecur_COL4A1 <- data.frame(coxtidy_85yrrecur_COL4A1)
  
  fooCox_COL4A1_90K_WHITE_FEMALE <- rbind(coxtidy_65yrrecur_COL4A1, coxtidy_85yrrecur_COL4A1)
  
  GNSIS_select_COL4A1_MALE <- GNSIS_select_COL4A1[GNSIS_select_COL4A1$PT_SEX == "Male", ]
  
  #cox_65yrrecur_COL4A1 <- coxph(Surv(timediff_recur_censor_65yr/365.25, event_recur_censor_65yr)~as.factor(group), data = GNSIS_select_COL4A1_MALE)
  #add conditional covariate
  cox_65yrrecur_COL4A1 <- coxph(Surv(timediff_recur_censor_65yr/365.25, event_recur_censor_65yr)~as.factor(group) + as.factor(hypercoagulability), data = GNSIS_select_COL4A1_MALE)
  
  coxtidy_65yrrecur_COL4A1 <- tidy(cox_65yrrecur_COL4A1, exponentiate = T, conf.int = T, conf.level = 0.95)
  coxtidy_65yrrecur_COL4A1$Onset <- "Early"
  coxtidy_65yrrecur_COL4A1$Category <- "Cerebral ischemia"
  coxtidy_65yrrecur_COL4A1$Gene <- i
  coxtidy_65yrrecur_COL4A1$Sex <- "Male"
  coxtidy_65yrrecur_COL4A1 <- data.frame(coxtidy_65yrrecur_COL4A1)
  
  #cox_85yrrecur_COL4A1 <- coxph(Surv(timediff_recur_censor_85yr/365.25, event_recur_censor_85yr)~as.factor(group), data = GNSIS_select_COL4A1_MALE)
  #add conditional covariate
  cox_85yrrecur_COL4A1 <- coxph(Surv(timediff_recur_censor_85yr/365.25, event_recur_censor_85yr)~as.factor(group) + as.factor(hypercoagulability), data = GNSIS_select_COL4A1_MALE)
  
  coxtidy_85yrrecur_COL4A1 <- tidy(cox_85yrrecur_COL4A1, exponentiate = T, conf.int = T, conf.level = 0.95)
  coxtidy_85yrrecur_COL4A1$Onset <- "Lifetime"
  coxtidy_85yrrecur_COL4A1$Category <- "Cerebral ischemia"
  coxtidy_85yrrecur_COL4A1$Gene <- i
  coxtidy_85yrrecur_COL4A1$Sex <- "Male"
  coxtidy_85yrrecur_COL4A1 <- data.frame(coxtidy_85yrrecur_COL4A1)
  
  fooCox_COL4A1_90K_WHITE_MALE <- rbind(coxtidy_65yrrecur_COL4A1, coxtidy_85yrrecur_COL4A1)
  
  fooCox_COL4A1_90K_WHITE_MALE_FEMALE <- rbind(fooCox_COL4A1_90K_WHITE_MALE, fooCox_COL4A1_90K_WHITE_FEMALE)
  
  fooCox_COL4A1_90K_WHITE <- rbind(fooCox_COL4A1_90K_WHITE, fooCox_COL4A1_90K_WHITE_MALE_FEMALE)
  listOfDataFrames[[i]] <- fooCox_COL4A1_90K_WHITE
}

fooCox_STROKE_90K_WHITE <- do.call("rbind", listOfDataFrames)
fooCox_STROKE_90K_WHITE$Category[fooCox_STROKE_90K_WHITE$Category == "Cerebral ischemia"] <- "Hypercoagulability"
fooCox_STROKE_90K_WHITE$Category[fooCox_STROKE_90K_WHITE$Category == "Cerebral ischemia"] <- "Primary Cardiomyopathy"

save(fooCox_STROKE_90K_WHITE, file = "fooCox_stroke_85K_WHITE.RData", version = 2) 
save(fooCox_STROKE_90K_WHITE, file = "fooCox_hypercoagulability_85K_WHITE.RData", version = 2) 
save(fooCox_STROKE_90K_WHITE, file = "fooCox_cardiomyopathy_85K_WHITE.RData", version = 2) 

save(fooCox_STROKE_90K_WHITE, file = "fooCox_stroke_10phecodes_90K_WHITE.RData", version = 2) 
save(fooCox_STROKE_90K_WHITE, file = "fooCox_stroke_10phecodes_85K_WHITE.RData", version = 2) 
load("fooCox_stroke_10phecodes_90K_WHITE.RData")
load("fooCox_stroke_90K_WHITE.RData")

fooCox_STROKE_90K_WHITE_NOTCH3 <- fooCox_STROKE_90K_WHITE[fooCox_STROKE_90K_WHITE$Gene == "NOTCH3", ]
fooCox_STROKE_90K_WHITE_NOTCH3$estimate <- 1/fooCox_STROKE_90K_WHITE_NOTCH3$estimate
fooCox_STROKE_90K_WHITE_NOTCH3$conf.low <- 1/fooCox_STROKE_90K_WHITE_NOTCH3$conf.low
fooCox_STROKE_90K_WHITE_NOTCH3$conf.high <- 1/fooCox_STROKE_90K_WHITE_NOTCH3$conf.high
fooCox_STROKE_90K_WHITE_NOTCH3

fooCox_STROKE_90K_WHITE_TREX1 <- fooCox_STROKE_90K_WHITE[fooCox_STROKE_90K_WHITE$Gene == "TREX1", ]
fooCox_STROKE_90K_WHITE_TREX1$estimate <- 1/fooCox_STROKE_90K_WHITE_TREX1$estimate
fooCox_STROKE_90K_WHITE_TREX1$conf.low <- 1/fooCox_STROKE_90K_WHITE_TREX1$conf.low
fooCox_STROKE_90K_WHITE_TREX1$conf.high <- 1/fooCox_STROKE_90K_WHITE_TREX1$conf.high
fooCox_STROKE_90K_WHITE_TREX1

fooCox_STROKE_90K_WHITE_COL4A1 <- fooCox_STROKE_90K_WHITE[fooCox_STROKE_90K_WHITE$Gene == "COL4A1", ]
fooCox_STROKE_90K_WHITE_COL4A1$estimate <- 1/fooCox_STROKE_90K_WHITE_COL4A1$estimate
fooCox_STROKE_90K_WHITE_COL4A1$conf.low <- 1/fooCox_STROKE_90K_WHITE_COL4A1$conf.low
fooCox_STROKE_90K_WHITE_COL4A1$conf.high <- 1/fooCox_STROKE_90K_WHITE_COL4A1$conf.high
fooCox_STROKE_90K_WHITE_COL4A1

save(fooCox_STROKE_90K_WHITE, file = "fooCox_stroke_90K_WHITE_condition_hypercoagulability.RData", version = 2)
################
#_________________________________________________________________________
#this is for COL4A1 in BLACK only
GNSIS_select_COL4A1 <- GNSIS_select[GNSIS_select$subgroup %in% c("COL4A1", "W/O Mut"), ]
#GNSIS_select_ALL <- GNSIS_select

cox_65yrrecur_COL4A1 <- coxph(Surv(timediff_recur_censor_65yr/365.25, event_recur_censor_65yr)~as.factor(group) + as.factor(PT_SEX), data = GNSIS_select_COL4A1)
coxtidy_65yrrecur_COL4A1 <- tidy(cox_65yrrecur_COL4A1, exponentiate = T, conf.int = T, conf.level = 0.95)
coxtidy_65yrrecur_COL4A1$Onset <- "Early"
coxtidy_65yrrecur_COL4A1$Category <- "Cerebral ischemia"
coxtidy_65yrrecur_COL4A1$Gene <- "COL4A1"
coxtidy_65yrrecur_COL4A1$Sex <- "Combined"
coxtidy_65yrrecur_COL4A1 <- data.frame(coxtidy_65yrrecur_COL4A1)

cox_85yrrecur_COL4A1 <- coxph(Surv(timediff_recur_censor_85yr/365.25, event_recur_censor_85yr)~as.factor(group) + as.factor(PT_SEX), data = GNSIS_select_COL4A1)
coxtidy_85yrrecur_COL4A1 <- tidy(cox_85yrrecur_COL4A1, exponentiate = T, conf.int = T, conf.level = 0.95)
coxtidy_85yrrecur_COL4A1$Onset <- "Lifetime"
coxtidy_85yrrecur_COL4A1$Category <- "Cerebral ischemia"
coxtidy_85yrrecur_COL4A1$Gene <- "COL4A1"
coxtidy_85yrrecur_COL4A1$Sex <- "Combined"
coxtidy_85yrrecur_COL4A1 <- data.frame(coxtidy_85yrrecur_COL4A1)

fooCox_COL4A1_90K_WHITE <- rbind(coxtidy_65yrrecur_COL4A1, coxtidy_85yrrecur_COL4A1)

GNSIS_select_COL4A1_FEMALE <- GNSIS_select_COL4A1[GNSIS_select_COL4A1$PT_SEX == "Female", ]

cox_65yrrecur_COL4A1 <- coxph(Surv(timediff_recur_censor_65yr/365.25, event_recur_censor_65yr)~as.factor(group), data = GNSIS_select_COL4A1_FEMALE)
coxtidy_65yrrecur_COL4A1 <- tidy(cox_65yrrecur_COL4A1, exponentiate = T, conf.int = T, conf.level = 0.95)
coxtidy_65yrrecur_COL4A1$Onset <- "Early"
coxtidy_65yrrecur_COL4A1$Category <- "Cerebral ischemia"
coxtidy_65yrrecur_COL4A1$Gene <- "COL4A1"
coxtidy_65yrrecur_COL4A1$Sex <- "Female"
coxtidy_65yrrecur_COL4A1 <- data.frame(coxtidy_65yrrecur_COL4A1)


cox_85yrrecur_COL4A1 <- coxph(Surv(timediff_recur_censor_85yr/365.25, event_recur_censor_85yr)~as.factor(group), data = GNSIS_select_COL4A1_FEMALE)
coxtidy_85yrrecur_COL4A1 <- tidy(cox_85yrrecur_COL4A1, exponentiate = T, conf.int = T, conf.level = 0.95)
coxtidy_85yrrecur_COL4A1$Onset <- "Lifetime"
coxtidy_85yrrecur_COL4A1$Category <- "Cerebral ischemia"
coxtidy_85yrrecur_COL4A1$Gene <- "COL4A1"
coxtidy_85yrrecur_COL4A1$Sex <- "Female"
coxtidy_85yrrecur_COL4A1 <- data.frame(coxtidy_85yrrecur_COL4A1)

fooCox_COL4A1_90K_WHITE_FEMALE <- rbind(coxtidy_65yrrecur_COL4A1, coxtidy_85yrrecur_COL4A1)

GNSIS_select_COL4A1_MALE <- GNSIS_select_COL4A1[GNSIS_select_COL4A1$PT_SEX == "Male", ]

cox_65yrrecur_COL4A1 <- coxph(Surv(timediff_recur_censor_65yr/365.25, event_recur_censor_65yr)~as.factor(group), data = GNSIS_select_COL4A1_MALE)
coxtidy_65yrrecur_COL4A1 <- tidy(cox_65yrrecur_COL4A1, exponentiate = T, conf.int = T, conf.level = 0.95)
coxtidy_65yrrecur_COL4A1$Onset <- "Early"
coxtidy_65yrrecur_COL4A1$Category <- "Cerebral ischemia"
coxtidy_65yrrecur_COL4A1$Gene <- "COL4A1"
coxtidy_65yrrecur_COL4A1$Sex <- "Male"
coxtidy_65yrrecur_COL4A1 <- data.frame(coxtidy_65yrrecur_COL4A1)

cox_85yrrecur_COL4A1 <- coxph(Surv(timediff_recur_censor_85yr/365.25, event_recur_censor_85yr)~as.factor(group), data = GNSIS_select_COL4A1_MALE)
coxtidy_85yrrecur_COL4A1 <- tidy(cox_85yrrecur_COL4A1, exponentiate = T, conf.int = T, conf.level = 0.95)
coxtidy_85yrrecur_COL4A1$Onset <- "Lifetime"
coxtidy_85yrrecur_COL4A1$Category <- "Cerebral ischemia"
coxtidy_85yrrecur_COL4A1$Gene <- "COL4A1"
coxtidy_85yrrecur_COL4A1$Sex <- "Male"
coxtidy_85yrrecur_COL4A1 <- data.frame(coxtidy_85yrrecur_COL4A1)

fooCox_COL4A1_90K_WHITE_MALE <- rbind(coxtidy_65yrrecur_COL4A1, coxtidy_85yrrecur_COL4A1)

fooCox_COL4A1_90K_WHITE_MALE_FEMALE <- rbind(fooCox_COL4A1_90K_WHITE_MALE, fooCox_COL4A1_90K_WHITE_FEMALE)

fooCox_COL4A1_90K_WHITE <- rbind(fooCox_COL4A1_90K_WHITE, fooCox_COL4A1_90K_WHITE_MALE_FEMALE)




listOfDataFrames <- list(fooCox_ALL_90K_WHITE, 
                         fooCox_COL4A1_90K_WHITE,
                         fooCox_CTC1_90K_WHITE,
                         fooCox_GLA_90K_WHITE,
                         fooCox_HTRA1_90K_WHITE,
                         fooCox_TREX1_90K_WHITE,
                         fooCox_NOTCH3_90K_WHITE)

fooCox_STROKE_90K_WHITE <- do.call("rbind", listOfDataFrames)
#_______________________________________________________________________________________________
#fooCox_STROKE_90K_WHITE <- fooCox_COL4A1_90K_WHITE #used for only combined analysis
result_final_Cox_select <- fooCox_STROKE_90K_WHITE %>%
  filter(term == 'as.factor(group)control')
result_final_Cox_select$p.value <- round(result_final_Cox_select$p.value, 3)
result_final_Cox_select$estimate_95CI <- paste0(round(1/result_final_Cox_select$estimate, 3), "(", round(1/result_final_Cox_select$conf.high, 3), "-", round(1/result_final_Cox_select$conf.low, 3), ")")

colnames(result_final_Cox_select)
# [1] "term"      "estimate"  "std.error" "statistic" "p.value"   "conf.low"  "conf.high" "Onset"     "Category" 
# [10] "Gene"      "Sex"   "estimate_95CI"

#create plot without pregnancy complication.
result_final_Cox_select <- result_final_Cox_select[!result_final_Cox_select$Gene %in% c("CTC1", "HTRA1"), ]
#create a forest plot
tabletext <- result_final_Cox_select %>%
  dplyr::select(Gene, Sex, Onset, p.value, estimate_95CI)

head(tabletext)

header <- c("Gene", "Sex", "Onset", "P value", "Hazard Ratio (95%CI)")
tabletext <- rbind(header, tabletext)
forestdata <- result_final_Cox_select %>%
  dplyr::select(estimate, conf.high, conf.low,)
forestdata$estimate <- 1/forestdata$estimate
forestdata$conf.low <- 1/forestdata$conf.low
forestdata$conf.high <- 1/forestdata$conf.high
head(forestdata)
header2 <- c(1, 1, 1)
forestdata <- rbind(header2, forestdata)
forestdata_new <- forestdata[c(1:3, 4:7),] #43 or 31 or 13
library(forestplot)
tabletext_new <- tabletext[c(1:3, 4:7),] 

head(tabletext_new)

tiff("Rplot_forestplot_multivariate_coxregression_disease_category_85K_stroke_BLACK_combined_65yr_limited.tiff", units = "in", width = 8, height = 4, res = 300) #height change to 8
#tiff("Rplot_forestplot_multivariate_coxregression_disease_category_90K_cardiomyopathy_white_65yr.tiff", units = "in", width = 8, height = 15, res = 300)

p <- forestplot(tabletext_new, 
                graph.pos = 3,
                forestdata_new,new_page = TRUE,
                #is.summary=c(TRUE,TRUE,rep(FALSE,8),TRUE),
                clip=c(-2,6),
                boxsize=0.25,
                zero = 1,
                xlog=FALSE, 
                col=fpColors(box="royalblue",line="darkblue", summary="royalblue"),
                xlab = "estimate(95%CI)",
                txt_gp = fpTxtGp(label = list(gpar(fontfamily = "Arial"),
                                              gpar(fontfamily = "",
                                                   col = "#660000")),
                                 ticks = gpar(fontfamily = "", cex = 1.0),
                                 xlab  = gpar(fontfamily = "Arial", cex = 1.0)))

grid::grid.text("STROKE ~ having rare variants + (Sex) + PCs(1-5)", .5, 0.99, gp=gpar(cex=1.5))
dev.off()




###############################################
load("phenotypes109.RData")
load("12092021_SKAT_inputfile_p123_90K_updated_phenotypes109.RData") #to get phenotype_ratio2 file for discovery
load("12092021_SKAT_inputfile_p123_85K_updated_phenotypes109.RData") #to get phenotype_ratio2 file for replication
phenotypes109_ratio2 <- phenotypes109[phenotypes109$id %in% phenotypes_ratio2$id, ] 
phenotypes109_ratio2$id
table(phenotypes109_ratio2$group)
# case control 
# 2745    5490
# case control 
# 1705    3410
phenotypes_ratio2 <- phenotypes109_ratio2 #from discovery or replication
colnames(phenotypes_ratio2)
gathercols <- colnames(phenotypes_ratio2[, 2:1858]) 
keycol <- "phecode"
valuecol <- "phenotype"
phenotypes_long <- data.frame(gather_(phenotypes_ratio2, keycol, valuecol, gathercols, factor_key=TRUE))

#phenotypes_matched_ratio2
head(phenotypes_long)
# id   group phecode phenotype
# 1  64 control     008         0
# 2 185 control     008         0
# 3 221    case     008         0
# 4 319 control     008         0
# 5 467 control     008         0
# 6 517 control     008         0
# id   group phecode phenotype
# 1  287 control     008         0
# 2 1136 control     008         0
# 3 1214 control     008         0
# 4 1510 control     008         0
# 5 1604 control     008         0
# 6 1692 control     008         0
str(phenotypes_long)
phenotypes_long$phecode <- paste0("X", phenotypes_long$phecode)
phecode <- unique(phenotypes_long$phecode)
phecode

load("dat_RECUR_ICD109.RData")
dat_RECUR_ICD109$phecode <- paste0("X", dat_RECUR_ICD109$phecode)
str(dat_RECUR_ICD109)
getwd()
load("D:/Geisinger/X19981_backup_10072020/Desktop/MS_project/Summary_phenotypes_merge_anno_common_enriched_85K.RData")
load("D:/Geisinger/X19981_backup_10072020/Desktop/MS_project/Summary_phenotypes_merge_anno_common_enriched_90K.RData") #in 85K using 90K enriched phecodes

listOfDataFrames <- list(Summary_phenotypes_merge_anno_common_enriched_CS,
                         Summary_phenotypes_merge_anno_common_enriched_GU,
                         Summary_phenotypes_merge_anno_common_enriched_EM,
                         Summary_phenotypes_merge_anno_common_enriched_DM,
                         Summary_phenotypes_merge_anno_common_enriched_DI,
                         Summary_phenotypes_merge_anno_common_enriched_RS,
                         Summary_phenotypes_merge_anno_common_enriched_NE,
                         Summary_phenotypes_merge_anno_common_enriched_HE,
                         Summary_phenotypes_merge_anno_common_enriched_IP,
                         Summary_phenotypes_merge_anno_common_enriched_SO,
                         Summary_phenotypes_merge_anno_common_enriched_MD,
                         Summary_phenotypes_merge_anno_common_enriched_MS,
                         Summary_phenotypes_merge_anno_common_enriched_SY,
                         Summary_phenotypes_merge_anno_common_enriched_PC,
                         Summary_phenotypes_merge_anno_common_enriched_ID,
                         Summary_phenotypes_merge_anno_common_enriched_NEO,
                         Summary_phenotypes_merge_anno_common_enriched_CA)

Summary_phenotypes_merge_anno_common_enriched <- do.call("rbind", listOfDataFrames)

unique(Summary_phenotypes_merge_anno_common_enriched$exclude_name)
#select enriched phecodes
plot_list = list()
plot_listq = list()
p <- NULL
q <- NULL
GNSIS_select <- NULL
for (i in unique(Summary_phenotypes_merge_anno_common_enriched$exclude_name)) {
  print(i)
  Summary_phenotypes_merge_anno_common_enriched_select <- Summary_phenotypes_merge_anno_common_enriched[Summary_phenotypes_merge_anno_common_enriched$exclude_name == i, ]
  phecode_select <- as.character(paste0("X", Summary_phenotypes_merge_anno_common_enriched_select$phecode))
  print(phecode_select)
  
  phenotypes_long_AA <- phenotypes_long %>% #change to relatedness or nonrelatedness files
    #dplyr::filter(phecode == "X433.3") %>% #Cerebral ischemia for X433; X433.3; #Cardiomyopathy for X425; X425.1
    dplyr::filter(phecode %in% phecode_select) %>% #select of enriched phecode set
    dplyr::filter(phenotype == 1)
  AA_ID <- unique(paste0("PT", "", phenotypes_long_AA$id))
  
  phenotypes_long_AA <- phenotypes_long %>% #change to relatedness or nonrelatedness files
    #dplyr::filter(phecode == "X433.3") %>% #Cerebral ischemia for X433; #Cardiomyopathy for X425
    dplyr::filter(phecode %in% phecode_select) %>% #select of enriched phecode set
    group_by(id) %>%
    dplyr::filter(sum(phenotype) == 0) %>%
    ungroup()
  AA_WO <- unique(paste0("PT", "", phenotypes_long_AA$id))
  
  dat_RECUR_ICD109$disease <- ifelse(dat_RECUR_ICD109$PT_ID %in% AA_ID, 1, 0)
  dat_RECUR_ICD109_select <- dat_RECUR_ICD109[dat_RECUR_ICD109$PT_ID %in% AA_ID, ] %>%
    #dplyr::filter(phecode == "X433.3") %>%
    dplyr::filter(phecode %in% phecode_select) %>% #select of enriched phecode set
    group_by(PT_ID) %>%
    arrange(desc(ENC_DT)) %>%
    slice(n()) %>%
    ungroup()
  dat_RECUR_ICD109_select <- data.frame(dat_RECUR_ICD109_select)
  print(table(dat_RECUR_ICD109_select$group))
  
  dat_RECUR_ICD109_WO <- dat_RECUR_ICD109[dat_RECUR_ICD109$PT_ID %in% AA_WO, ] %>%
    dplyr::group_by(PT_ID) %>%
    dplyr::arrange(ENC_DT) %>%
    slice(n()) %>%
    #mutate(disease = 0) %>%
    ungroup()
  dat_RECUR_ICD109_WO <- data.frame(dat_RECUR_ICD109_WO)
  GNSIS <- rbind(dat_RECUR_ICD109_select, dat_RECUR_ICD109_WO)
  
  GNSIS <- GNSIS %>% 
    mutate(timediff_recur_censor = 
             case_when(disease == 1 ~ (timediff_occur_censor - 0)*365.25,
                       disease == 0 ~ (timediff_occur_censor - 0)*365.25)) 
  GNSIS <- GNSIS %>% 
    mutate(timediff_recur_censor_15yr =
             case_when(timediff_recur_censor < 15*365.25 ~ timediff_recur_censor,
                       timediff_recur_censor >= 15*365.25 ~ 15*365.25),
           event_recur_censor_15yr =
             case_when(disease == 0 ~ 0,
                       disease == 1 & timediff_recur_censor <= 15*365.25 ~1,
                       disease == 1 & timediff_recur_censor > 15*365.25 ~0),
           timediff_recur_censor_25yr =
             case_when(timediff_recur_censor < 25*365.25 ~ timediff_recur_censor,
                       timediff_recur_censor >= 25*365.25 ~ 25*365.25),
           event_recur_censor_25yr =
             case_when(disease == 0 ~ 0,
                       disease == 1 & timediff_recur_censor <= 25*365.25 ~1,
                       disease == 1 & timediff_recur_censor > 25*365.25 ~0),
           timediff_recur_censor_35yr =
             case_when(timediff_recur_censor < 35*365.25 ~ timediff_recur_censor,
                       timediff_recur_censor >= 35*365.25 ~ 35*365.25),
           event_recur_censor_35yr =
             case_when(disease == 0 ~ 0,
                       disease == 1 & timediff_recur_censor <= 35*365.25 ~1,
                       disease == 1 & timediff_recur_censor > 35*365.25 ~0),
           timediff_recur_censor_65yr =
             case_when(timediff_recur_censor < 65*365.25 ~ timediff_recur_censor,
                       timediff_recur_censor >= 65*365.25 ~ 65*365.25),
           event_recur_censor_65yr =
             case_when(disease == 0 ~ 0,
                       disease == 1 & timediff_recur_censor <= 65*365.25 ~1,
                       disease == 1 & timediff_recur_censor > 65*365.25 ~0)
           )
  
  #optional
  # GNSIS_select <- GNSIS %>%
  #   filter(timediff_recur_censor_35yr >= 0)
  
  GNSIS_select <- GNSIS
  # KM stratified by cohort
  ANN <- phecode_mapping_unique[phecode_mapping_unique$phecodeX %in% phecode_select, ]$phecode_str
  if (i == "endocrine/metabolic") {
    q <- ggsurvplot(
      surv_fit(Surv(timediff_recur_censor_65yr/365.25, event_recur_censor_65yr)~as.factor(group), data = GNSIS_select), # survfit object with calculated statistics.
      fun = function(x) {1-x},  #plot cumulative probability F(t) = 1 - S(t)
      risk.table = TRUE,       # show risk table.
      pval = TRUE,             # show p-value of log-rank test.
      pval.coord = c(2, 1),
      pval.size = 8,
      conf.int = TRUE,       # show confidence intervals for 
      censor.shape="|", censor.size = 4,
      fontsize = 4, #font size in the table
      font.legend = c(14, "plain", "black"),
      # point estimates of survival curves.
      xlim = c(0, 65),         # present narrower X axis, but not affect
      ylim = c(0, 1),
      # survival estimates.
      xlab = "Time in years",   # customize X axis label.
      ylab = "Cumulative Incidence",
      break.time.by = 1,     # break X axis in time intervals by 500.
      ggtheme = theme_light(), # customize plot and risk table with a theme.
      n.risk = FALSE,
      legend.labs =
        c("w/ rare alleles", "w/o rare alleles"),    # change legend labels.
      font.title    = c(14, "bold.italic", "darkgreen"),
      font.subtitle = c(14, "bold", "darkgreen"),
      font.caption  = c(14, "plain", "darkgreen"),
      font.x        = c(14, "bold.italic", "black"),
      font.y        = c(14, "bold.italic", "black"),
      font.xtickslab = c(14, "bold", "black"),
      font.ytickslab = c(14, "bold", "black")
    ) +
      labs(
        title    = "Cumulative probability for 65yrs",
        subtitle = paste('Disease = ', i, length(phecode_select), 'enriched phecodes'),
        #subtitle = paste('Disease = ', "Cerebral ischemia (X433.3)"),
        #subtitle = paste('Disease = ', "Cardiomyopathy (425)"),
        caption  = "Plotted with R survminer")
    q$plot <- q$plot + ggplot2::annotate("text", x = 15, y = seq(0,0.95,length=length(phecode_select)),# x and y coordinates of the text
                                         label = ANN[1:length(ANN)], size = 5)
    plot_listq[[i]] = q
  } else {
    p <- ggsurvplot(
      surv_fit(Surv(timediff_recur_censor_65yr/365.25, event_recur_censor_65yr)~as.factor(group), data = GNSIS_select), # survfit object with calculated statistics.
      fun = function(x) {1-x},  #plot cumulative probability F(t) = 1 - S(t)
      risk.table = TRUE,       # show risk table.
      pval = TRUE,             # show p-value of log-rank test.
      pval.coord = c(2, 0.5),
      pval.size = 8,
      conf.int = TRUE,       # show confidence intervals for 
      censor.shape="|", censor.size = 4,
      fontsize = 4, #font size in the table
      font.legend = c(14, "plain", "black"),
      # point estimates of survival curves.
      xlim = c(0, 65),         # present narrower X axis, but not affect
      ylim = c(0, 0.6),
      # survival estimates.
      xlab = "Time in years",   # customize X axis label.
      ylab = "Cumulative Incidence",
      break.time.by = 1,     # break X axis in time intervals by 500.
      ggtheme = theme_light(), # customize plot and risk table with a theme.
      n.risk = FALSE,
      legend.labs =
        c("w/ rare alleles", "w/o rare alleles"),    # change legend labels.
      font.title    = c(14, "bold.italic", "darkgreen"),
      font.subtitle = c(14, "bold", "darkgreen"),
      font.caption  = c(14, "plain", "darkgreen"),
      font.x        = c(14, "bold.italic", "black"),
      font.y        = c(14, "bold.italic", "black"),
      font.xtickslab = c(14, "bold", "black"),
      font.ytickslab = c(14, "bold", "black")
    ) +
      labs(
        title    = "Cumulative probability for 65yrs",
        subtitle = paste('Disease = ', i, length(phecode_select), 'enriched phecodes'),
        #subtitle = paste('Disease = ', "Cerebral ischemia (X433.3)"),
        #subtitle = paste('Disease = ', "Cardiomyopathy (425)"),
        caption  = "Plotted with R survminer")
    p$plot <- p$plot + ggplot2::annotate("text", x = 15, y = seq(0,0.45,length=length(phecode_select)),# x and y coordinates of the text
                                         label = ANN[1:length(ANN)], size = 5)
    plot_list[[i]] = p
  }
  # 
  # p <- ggsurvplot(
  #   surv_fit(Surv(timediff_recur_censor_25yr/365.25, event_recur_censor_25yr)~as.factor(group), data = GNSIS_select), # survfit object with calculated statistics.
  #   fun = function(x) {1-x},  #plot cumulative probability F(t) = 1 - S(t)
  #   risk.table = TRUE,       # show risk table.
  #   pval = TRUE,             # show p-value of log-rank test.
  #   pval.coord = c(2, 0.5),
  #   pval.size = 8,
  #   conf.int = TRUE,       # show confidence intervals for 
  #   censor.shape="|", censor.size = 4,
  #   fontsize = 4, #font size in the table
  #   font.legend = c(14, "plain", "black"),
  #   # point estimates of survival curves.
  #   xlim = c(0, 25),         # present narrower X axis, but not affect
  #   ylim = c(0, 0.5),
  #   # survival estimates.
  #   xlab = "Time in years",   # customize X axis label.
  #   ylab = "Cumulative Incidence",
  #   break.time.by = 1,     # break X axis in time intervals by 500.
  #   ggtheme = theme_light(), # customize plot and risk table with a theme.
  #   n.risk = FALSE,
  #   legend.labs =
  #     c("w/ rare alleles", "w/o rare alleles"),    # change legend labels.
  #   font.title    = c(14, "bold.italic", "darkgreen"),
  #   font.subtitle = c(14, "bold", "darkgreen"),
  #   font.caption  = c(14, "plain", "darkgreen"),
  #   font.x        = c(14, "bold.italic", "black"),
  #   font.y        = c(14, "bold.italic", "black"),
  #   font.xtickslab = c(14, "bold", "black"),
  #   font.ytickslab = c(14, "bold", "black")
  # ) +
  #   labs(
  #     title    = "Cumulative probability for 25yrs",
  #     subtitle = paste('Disease = ', i, length(phecode_select), 'enriched phecodes'),
  #     #subtitle = paste('Disease = ', "Cerebral ischemia (X433.3)"),
  #     #subtitle = paste('Disease = ', "Cardiomyopathy (425)"),
  #     caption  = "Plotted with R survminer")
  # p$plot <- p$plot + ggplot2::annotate("text", x = 5, y = seq(0,0.45,length=length(phecode_select)),# x and y coordinates of the text
  #   label = ANN[1:length(ANN)], size = 5)
  # plot_list[[i]] = p
}

#dir.create("D:/Geisinger/X19981_backup_10072020/Desktop/MS_project/survivalcurve_85K_90Kphecodes/")
dir.create("D:/Geisinger/X19981_backup_10072020/Desktop/MS_project/survivalcurve_90K_85Kphecodes/")
for (i in unique(Summary_phenotypes_merge_anno_common_enriched$exclude_name)) {
  tiff(paste('D:/Geisinger/X19981_backup_10072020/Desktop/MS_project/survivalcurve_90K_85Kphecodes/', 'Cerebralvascular_85K_', i, '_enriched_65yr', '.tiff', sep=''), units="in", width=15, height=15, res=600)
  if (i == "endocrine/metabolic") {
    print(plot_listq[[i]])
    dev.off()
  } else {
    print(plot_list[[i]])
    dev.off()
  }
}

fooCox = list()
#KM analysis
for (i in unique(Summary_phenotypes_merge_anno_common_enriched$exclude_name)) {
  print(i)
  Summary_phenotypes_merge_anno_common_enriched_select <- Summary_phenotypes_merge_anno_common_enriched[Summary_phenotypes_merge_anno_common_enriched$exclude_name == i, ]
  phecode_select <- as.character(paste0("X", Summary_phenotypes_merge_anno_common_enriched_select$phecode))
  print(phecode_select)
  
  phenotypes_long_AA <- phenotypes_long %>% #change to relatedness or nonrelatedness files
    #dplyr::filter(phecode == "X433.3") %>% #Cerebral ischemia for X433; X433.3; #Cardiomyopathy for X425; X425.1
    dplyr::filter(phecode %in% phecode_select) %>% #select of enriched phecode set
    dplyr::filter(phenotype == 1)
  AA_ID <- unique(paste0("PT", "", phenotypes_long_AA$id))
  
  phenotypes_long_AA <- phenotypes_long %>% #change to relatedness or nonrelatedness files
    #dplyr::filter(phecode == "X433.3") %>% #Cerebral ischemia for X433; #Cardiomyopathy for X425
    dplyr::filter(phecode %in% phecode_select) %>% #select of enriched phecode set
    group_by(id) %>%
    dplyr::filter(sum(phenotype) == 0) %>%
    ungroup()
  AA_WO <- unique(paste0("PT", "", phenotypes_long_AA$id))
  
  dat_RECUR_ICD109$disease <- ifelse(dat_RECUR_ICD109$PT_ID %in% AA_ID, 1, 0)
  dat_RECUR_ICD109_select <- dat_RECUR_ICD109[dat_RECUR_ICD109$PT_ID %in% AA_ID, ] %>%
    #dplyr::filter(phecode == "X433.3") %>%
    dplyr::filter(phecode %in% phecode_select) %>% #select of enriched phecode set
    group_by(PT_ID) %>%
    arrange(desc(ENC_DT)) %>%
    slice(n()) %>%
    ungroup()
  dat_RECUR_ICD109_select <- data.frame(dat_RECUR_ICD109_select)
  print(table(dat_RECUR_ICD109_select$group))
  
  dat_RECUR_ICD109_WO <- dat_RECUR_ICD109[dat_RECUR_ICD109$PT_ID %in% AA_WO, ] %>%
    dplyr::group_by(PT_ID) %>%
    dplyr::arrange(ENC_DT) %>%
    slice(n()) %>%
    #mutate(disease = 0) %>%
    ungroup()
  dat_RECUR_ICD109_WO <- data.frame(dat_RECUR_ICD109_WO)
  GNSIS <- rbind(dat_RECUR_ICD109_select, dat_RECUR_ICD109_WO)
  
  GNSIS <- GNSIS %>% 
    mutate(timediff_recur_censor = 
             case_when(disease == 1 ~ (timediff_occur_censor - 0)*365.25,
                       disease == 0 ~ (timediff_occur_censor - 0)*365.25)) 
  GNSIS <- GNSIS %>% 
    mutate(timediff_recur_censor_15yr =
             case_when(timediff_recur_censor < 15*365.25 ~ timediff_recur_censor,
                       timediff_recur_censor >= 15*365.25 ~ 15*365.25),
           event_recur_censor_15yr =
             case_when(disease == 0 ~ 0,
                       disease == 1 & timediff_recur_censor <= 15*365.25 ~1,
                       disease == 1 & timediff_recur_censor > 15*365.25 ~0),
           timediff_recur_censor_25yr =
             case_when(timediff_recur_censor < 25*365.25 ~ timediff_recur_censor,
                       timediff_recur_censor >= 25*365.25 ~ 25*365.25),
           event_recur_censor_25yr =
             case_when(disease == 0 ~ 0,
                       disease == 1 & timediff_recur_censor <= 25*365.25 ~1,
                       disease == 1 & timediff_recur_censor > 25*365.25 ~0),
           timediff_recur_censor_35yr =
             case_when(timediff_recur_censor < 35*365.25 ~ timediff_recur_censor,
                       timediff_recur_censor >= 35*365.25 ~ 35*365.25),
           event_recur_censor_35yr =
             case_when(disease == 0 ~ 0,
                       disease == 1 & timediff_recur_censor <= 35*365.25 ~1,
                       disease == 1 & timediff_recur_censor > 35*365.25 ~0),
           timediff_recur_censor_65yr =
             case_when(timediff_recur_censor < 65*365.25 ~ timediff_recur_censor,
                       timediff_recur_censor >= 65*365.25 ~ 65*365.25),
           event_recur_censor_65yr =
             case_when(disease == 0 ~ 0,
                       disease == 1 & timediff_recur_censor <= 65*365.25 ~1,
                       disease == 1 & timediff_recur_censor > 65*365.25 ~0)
           )
  
  GNSIS_select <- GNSIS %>%
    filter(timediff_recur_censor_35yr >= 0)
  #optional
  #GNSIS_select <- GNSIS
  # KM stratified by cohort
  ANN <- phecode_mapping_unique[phecode_mapping_unique$phecodeX %in% phecode_select, ]$phecode_str
  scale(GNSIS_select$timediff_occur_censor)
  cox_65yrrecur_AA <- coxph(Surv(timediff_recur_censor_65yr/365.25, event_recur_censor_65yr)~as.factor(group) + as.factor(PT_SEX) + as.factor(PT_RACE), data = GNSIS_select)
  #cox_65yrrecur_AA <- coxph(Surv(timediff_recur_censor_65yr/365.25, event_recur_censor_65yr)~as.factor(group) + as.factor(PT_SEX), data = GNSIS_select)
  coxtidy_65yrrecur_AA <- tidy(cox_65yrrecur_AA, exponentiate = T, conf.int = T, conf.level = 0.95)
  coxtidy_65yrrecur_AA$Category <- i
  coxtidy_65yrrecur_AA <- data.frame(coxtidy_65yrrecur_AA)
  fooCox[[i]] <- coxtidy_65yrrecur_AA
}
result_final_Cox <- do.call(rbind, fooCox)
result_final_Cox$estimate_95CI <- paste0(round(1/result_final_Cox$estimate, 3), "(", round(1/result_final_Cox$conf.high, 3), "-", round(1/result_final_Cox$conf.low, 3), ")")

save(result_final_Cox, file = "result_final_Cox_90K_65yr.RData", version =2)
save(result_final_Cox, file = "result_final_Cox_85K_65yr.RData", version =2)
save(result_final_Cox, file = "result_final_Cox_85K_woRACE_65yr.RData", version =2)
save(result_final_Cox, file = "result_final_Cox_85K_90Kphecodes_65yr.RData", version =2)
save(result_final_Cox, file = "result_final_Cox_85K_90Kphecodes_woRACE_65yr.RData", version =2)
save(result_final_Cox, file = "result_final_Cox_90K_85Kphecodes_65yr.RData", version =2)
save(result_final_Cox, file = "result_final_Cox_90K_85Kphecodes_woRACE_65yr.RData", version =2)
save(result_final_Cox, file = "result_final_Cox_90K_woRACE_65yr.RData", version =2)

load("result_final_Cox_90K_65yr.RData")
load("result_final_Cox_90K_woRACE_65yr.RData")
load("result_final_Cox_85K_65yr.RData")
load("result_final_Cox_85K_woRACE_65yr.RData")
load("result_final_Cox_90K_85Kphecodes_65yr.RData")
load("result_final_Cox_90K_85Kphecodes_woRACE_65yr.RData")
load("result_final_Cox_85K_90Kphecodes_65yr.RData")
load("result_final_Cox_85K_90Kphecodes_woRACE_65yr.RData")


result_final_Cox_select <- result_final_Cox %>%
  filter(term == 'as.factor(group)control')
result_final_Cox_select$p.value <- round(result_final_Cox_select$p.value, 3)

colnames(result_final_Cox_select)
# [1] "term"          "estimate"      "std.error"     "statistic"     "p.value"       "conf.low"      "conf.high"    
# [8] "Category"      "estimate_95CI"

#create plot without pregnancy complication.
result_final_Cox_select <- result_final_Cox_select[result_final_Cox_select$Category != "pregnancy complications", ]
#create a forest plot
tabletext <- result_final_Cox_select %>%
  dplyr::select(Category, p.value, estimate_95CI)

head(tabletext)

header <- c("Disease Category", "P value", "Hazard Ratio (95%CI)")
tabletext <- rbind(header, tabletext)
forestdata <- result_final_Cox_select %>%
  dplyr::select(estimate, conf.high, conf.low,)
forestdata$estimate <- 1/forestdata$estimate
forestdata$conf.low <- 1/forestdata$conf.low
forestdata$conf.high <- 1/forestdata$conf.high
head(forestdata)
header2 <- c(1, 1, 1)
forestdata <- rbind(header2, forestdata)
forestdata_new <- forestdata[c(1:3, 4:18),] 
library(forestplot)
tabletext_new <- tabletext[c(1:3, 4:18),] 

head(tabletext_new)

tiff("Rplot_forestplot_multivariate_coxregression_disease_category_85K_90Kphecodes_woRACE_65yr_wopregnancy.tiff", units = "in", width = 10, height = 8, res = 300)
p <- forestplot(tabletext_new, 
                graph.pos = 3,
                forestdata_new,new_page = TRUE,
                #is.summary=c(TRUE,TRUE,rep(FALSE,8),TRUE),
                clip=c(-2,6),
                boxsize=0.25,
                zero = 1,
                xlog=FALSE, 
                col=fpColors(box="royalblue",line="darkblue", summary="royalblue"),
                xlab = "estimate(95%CI)",
                txt_gp = fpTxtGp(label = list(gpar(fontfamily = "Arial"),
                                              gpar(fontfamily = "",
                                                   col = "#660000")),
                                 ticks = gpar(fontfamily = "", cex = 1.0),
                                 xlab  = gpar(fontfamily = "Arial", cex = 1.0)))

grid::grid.text("Disease category ~ having rare variants + sex + race", .5, 0.99, gp=gpar(cex=1.5))
dev.off()

###########
# custom pathway
#metabolic system

Summary_phenotypes_merge_anno_common_enriched_select <- Summary_phenotypes_merge_anno_common_enriched[Summary_phenotypes_merge_anno_common_enriched$exclude_name == "endocrine/metabolic", ]
phecode_select <- as.character(paste0("X", Summary_phenotypes_merge_anno_common_enriched_select$phecode))
print(phecode_select)

phenotypes_long_AA <- phenotypes_long %>% #change to relatedness or nonrelatedness files
  #dplyr::filter(phecode == "X433.3") %>% #Cerebral ischemia for X433; X433.3; #Cardiomyopathy for X425; X425.1
  dplyr::filter(phecode %in% phecode_select) %>% #select of enriched phecode set
  dplyr::filter(phenotype == 1)
AA_ID <- unique(paste0("PT", "", phenotypes_long_AA$id))

phenotypes_long_AA <- phenotypes_long %>% #change to relatedness or nonrelatedness files
  #dplyr::filter(phecode == "X433.3") %>% #Cerebral ischemia for X433; #Cardiomyopathy for X425
  dplyr::filter(phecode %in% phecode_select) %>% #select of enriched phecode set
  group_by(id) %>%
  dplyr::filter(sum(phenotype) == 0) %>%
  ungroup()
AA_WO <- unique(paste0("PT", "", phenotypes_long_AA$id))

dat_RECUR_ICD109$disease <- ifelse(dat_RECUR_ICD109$PT_ID %in% AA_ID, 1, 0)
dat_RECUR_ICD109_select <- dat_RECUR_ICD109[dat_RECUR_ICD109$PT_ID %in% AA_ID, ] %>%
  #dplyr::filter(phecode == "X433.3") %>%
  dplyr::filter(phecode %in% phecode_select) %>% #select of enriched phecode set
  group_by(PT_ID) %>%
  arrange(desc(ENC_DT)) %>%
  slice(n()) %>%
  ungroup()
dat_RECUR_ICD109_select <- data.frame(dat_RECUR_ICD109_select)
print(table(dat_RECUR_ICD109_select$group))

dat_RECUR_ICD109_WO <- dat_RECUR_ICD109[dat_RECUR_ICD109$PT_ID %in% AA_WO, ] %>%
  dplyr::group_by(PT_ID) %>%
  dplyr::arrange(ENC_DT) %>%
  slice(n()) %>%
  #mutate(disease = 0) %>%
  ungroup()
dat_RECUR_ICD109_WO <- data.frame(dat_RECUR_ICD109_WO)
GNSIS <- rbind(dat_RECUR_ICD109_select, dat_RECUR_ICD109_WO)

GNSIS <- GNSIS %>% 
  mutate(timediff_recur_censor = 
           case_when(disease == 1 ~ (timediff_occur_censor - 0)*365.25,
                     disease == 0 ~ (timediff_occur_censor - 0)*365.25)) 
GNSIS <- GNSIS %>% 
  mutate(timediff_recur_censor_15yr =
           case_when(timediff_recur_censor < 15*365.25 ~ timediff_recur_censor,
                     timediff_recur_censor >= 15*365.25 ~ 15*365.25),
         event_recur_censor_15yr =
           case_when(disease == 0 ~ 0,
                     disease == 1 & timediff_recur_censor <= 15*365.25 ~1,
                     disease == 1 & timediff_recur_censor > 15*365.25 ~0),
         timediff_recur_censor_25yr =
           case_when(timediff_recur_censor < 25*365.25 ~ timediff_recur_censor,
                     timediff_recur_censor >= 25*365.25 ~ 25*365.25),
         event_recur_censor_25yr =
           case_when(disease == 0 ~ 0,
                     disease == 1 & timediff_recur_censor <= 25*365.25 ~1,
                     disease == 1 & timediff_recur_censor > 25*365.25 ~0),
         timediff_recur_censor_35yr =
           case_when(timediff_recur_censor < 35*365.25 ~ timediff_recur_censor,
                     timediff_recur_censor >= 35*365.25 ~ 35*365.25),
         event_recur_censor_35yr =
           case_when(disease == 0 ~ 0,
                     disease == 1 & timediff_recur_censor <= 35*365.25 ~1,
                     disease == 1 & timediff_recur_censor > 35*365.25 ~0),
         timediff_recur_censor_65yr =
           case_when(timediff_recur_censor < 65*365.25 ~ timediff_recur_censor,
                     timediff_recur_censor >= 65*365.25 ~ 65*365.25),
         event_recur_censor_65yr =
           case_when(disease == 0 ~ 0,
                     disease == 1 & timediff_recur_censor <= 65*365.25 ~1,
                     disease == 1 & timediff_recur_censor > 65*365.25 ~0)
  )

#optional
# GNSIS_select <- GNSIS %>%
#   filter(timediff_recur_censor_35yr >= 0)

GNSIS_select <- GNSIS
ANN <- phecode_mapping_unique[phecode_mapping_unique$phecodeX %in% phecode_select, ]$phecode_str
q <- ggsurvplot(
  surv_fit(Surv(timediff_recur_censor_65yr/365.25, event_recur_censor_65yr)~as.factor(group), data = GNSIS_select), # survfit object with calculated statistics.
  fun = function(x) {1-x},  #plot cumulative probability F(t) = 1 - S(t)
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  pval.coord = c(2, 0.5),
  pval.size = 8,
  conf.int = TRUE,       # show confidence intervals for 
  censor.shape="|", censor.size = 4,
  fontsize = 4, #font size in the table
  font.legend = c(14, "plain", "black"),
  # point estimates of survival curves.
  xlim = c(0, 65),         # present narrower X axis, but not affect
  ylim = c(0, 0.6),
  # survival estimates.
  xlab = "Time in years",   # customize X axis label.
  ylab = "Cumulative Incidence",
  break.time.by = 1,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  n.risk = FALSE,
  legend.labs =
    c("w/ rare alleles", "w/o rare alleles"),    # change legend labels.
  font.title    = c(14, "bold.italic", "darkgreen"),
  font.subtitle = c(14, "bold", "darkgreen"),
  font.caption  = c(14, "plain", "darkgreen"),
  font.x        = c(14, "bold.italic", "black"),
  font.y        = c(14, "bold.italic", "black"),
  font.xtickslab = c(14, "bold", "black"),
  font.ytickslab = c(14, "bold", "black")
) +
  labs(
    title    = "Cumulative probability for 65yrs",
    subtitle = paste('Disease = ', 'endocrine_metabolic', length(phecode_select), 'enriched phecodes'),
    #subtitle = paste('Disease = ', "Cerebral ischemia (X433.3)"),
    #subtitle = paste('Disease = ', "Cardiomyopathy (425)"),
    caption  = "Plotted with R survminer")
q$plot <- q$plot + ggplot2::annotate("text", x = 15, y = seq(0,0.45,length=length(phecode_select)),# x and y coordinates of the text
                                     label = ANN[1:length(ANN)], size = 5)
tiff(paste('D:/Geisinger/X19981_backup_10072020/Desktop/MS_project/survivalcurve_90K/', 'Cerebralvascular_85K_', 'endocrine_metabolic', '_enriched_65yr', '.tiff', sep=''), units="in", width=15, height=15, res=600)
print(q)
dev.off()

##custom phecodes
#select phecodes from top ranked list of a specific disease category
unique(Summary_phenotypes_merge_anno_common_enriched$exclude_name)
# [1] "circulatory system"      "genitourinary"           "endocrine/metabolic"     "dermatologic"           
# [5] "digestive"               "respiratory"             "neurological"            "hematopoietic"          
# [9] "injuries & poisonings"   "sense organs"            "mental disorders"        "musculoskeletal"        
# [13] "symptoms"                "pregnancy complications" "infectious diseases"     "neoplasms"              
# [17] "congenital anomalies"

#using hematopoietic as an example
Summary_phenotypes_merge_anno_common_enriched_select <- Summary_phenotypes_merge_anno_common_enriched[Summary_phenotypes_merge_anno_common_enriched$exclude_name == "hematopoietic", ]
phecode_select <- as.character(paste0("X", Summary_phenotypes_merge_anno_common_enriched_select$phecode))
print(phecode_select)

#select phecodes from a list of candidate phecodes
phecode_select <- c("X286.7", "X286.8", "X286.81")
phecode_select <- c("X286.8", "X286.81")
phecode_select <- c("X425", "X425.1")
phecode_select <- c("X425.1")
phecode_select <- c("X433.1", "X433.11", "X433.12", "X433.2", "X433.21", "X433.3", "X433.31", "X433.32", "X433.6", "X433.8")
phecode_select <- c("X433.2", "X433.21", "X433.3", "X433.31", "X433.32", "X433.6", "X433.8")
phecode_select <- c("X433.2", "X433.21", "X433.3", "X433.31", "X433.32", "X433.8")
phecode_select <- c("X433.2")

#get patient ID with RVs from a specific phecode_select
phenotypes_long_AA <- phenotypes_long %>%
  dplyr::filter(phecode %in% phecode_select) %>%
  dplyr::filter(phenotype == 1)
AA_ID <- unique(paste0("PT", "", phenotypes_long_AA$id))

#get patient ID without RVs from a specific phecode_select
phenotypes_long_AA <- phenotypes_long %>%
  dplyr::filter(phecode %in% phecode_select) %>%
  group_by(id) %>%
  dplyr::filter(sum(phenotype) == 0) %>%
  ungroup()
AA_WO <- unique(paste0("PT", "", phenotypes_long_AA$id))

##
dat_RECUR_ICD109$disease <- ifelse(dat_RECUR_ICD109$PT_ID %in% AA_ID, 1, 0)
dat_RECUR_ICD109_select <- dat_RECUR_ICD109[dat_RECUR_ICD109$PT_ID %in% AA_ID, ] %>%
  dplyr::filter(phecode %in% phecode_select) %>%
  group_by(PT_ID) %>%
  arrange(desc(ENC_DT)) %>%
  slice(n()) %>%
  ungroup()
dat_RECUR_ICD109_select <- data.frame(dat_RECUR_ICD109_select)
print(table(dat_RECUR_ICD109_select$group))
#85K X433.3
# case control 
# 62     124 
#90K X433.3
# case control 
# 188     358/359 
#85K X433
# case control 
# 140     243
#90K X433
# case control 
# 350     663

#85K X425
# case control 
# 51      88
# 90K X425
# case control 
# 81     156
#85K X425.1
# case control 
# 49      75 
#90K X425.1
# case control 
# 73     137 
dat_RECUR_ICD109_WO <- dat_RECUR_ICD109[dat_RECUR_ICD109$PT_ID %in% AA_WO, ] %>%
  dplyr::group_by(PT_ID) %>%
  dplyr::arrange(ENC_DT) %>%
  slice(n()) %>%
  ungroup()
dat_RECUR_ICD109_WO <- data.frame(dat_RECUR_ICD109_WO)
print(table(dat_RECUR_ICD109_WO$group))
#85K X433.3
# case control 
# 1643    3285
#90K X433.3
# case control 
# 2557    5131
#85K X433
# case control 
# 1565    3167
#90K X433
# case control 
# 2395    4827
#90K X425
# case control 
# 2664    5333 
#85K X425
# case control 
# 1654    3322 
#85K X425.1
# case control 
# 1656    3335
# 90K X425.1
# case control 
# 2672    5353
GNSIS <- rbind(dat_RECUR_ICD109_select, dat_RECUR_ICD109_WO)

################################################
##calculate IR (incidence rate)
#remove spikein cases
spikein_ID <- c(1117592, 1202698, 678791, 726999, 865829, 1360022, 1154976, 121012, 1346346, 242673, 259508,
                377514, 527399, 635221, 753905, 865414, 928617)
spikein_ID <- paste0("PT", spikein_ID)
GNSIS <- GNSIS[!GNSIS$PT_ID %in% spikein_ID, ] #change from 7834 to 7827

sum(GNSIS$timediff_occur_censor) #449512.6
table(GNSIS$group, GNSIS$disease)
#            0    1
# case    2518   68
# control 5110  131
table(GNSIS$group, GNSIS$disease)[1,2] #68
table(GNSIS$group, GNSIS$disease)[2,2] #131

#alternative expression of the count
nrow(GNSIS[GNSIS$group == "case" & GNSIS$disease == 1, ])
nrow(GNSIS[GNSIS$group == "control" & GNSIS$disease == 1, ])

sum(GNSIS[GNSIS$group == "case", ]$timediff_occur_censor) #149183.6

IR_withmutation <- table(GNSIS$group, GNSIS$disease)[1,2]/sum(GNSIS[GNSIS$group == "case", ]$timediff_occur_censor)*1000  
IR_withoutmutation <- table(GNSIS$group, GNSIS$disease)[2,2]/sum(GNSIS[GNSIS$group == "control", ]$timediff_occur_censor)*1000  
#Calculation of incidence rate ratios and 95% confidence intervals
#Reference for all calculations: Rothman KJ, Greenland S. Modern epidemiology. 3rd ed.
#Philadelphia: Lippincott Williams & Wilkins, 2008.
IR <- IR_withmutation/IR_withoutmutation

SD <- (1/table(GNSIS$group, GNSIS$disease)[1,2] + 1/table(GNSIS$group, GNSIS$disease)[2,2])^0.5

upper <- exp(log(IR) + 1.96*SD) 
lower <- exp(log(IR) - 1.96*SD) 

LIST <- c(IR_withmutation, IR_withoutmutation, IR, SD, lower, upper)

GNSIS <- NULL
LIST <- NULL
SD <- NULL
IR_withmutation <- NULL
IR_withoutmutation <- NULL
RD <- NULL
upper <- NULL
lower <- NULL
withmutation <- NULL
withoutmutation <- NULL

Summary_phenotypes_merge_anno_common$phecode <- as.character(paste0("X", Summary_phenotypes_merge_anno_common$phecode))
# str(Summary_phenotypes_merge_anno_common)
phecodes <- Summary_phenotypes_merge_anno_common$phecode
phecodes <- c("X433.1", "X433.11", "X433.12", "X433.2", "X433.21", "X433.3", "X433.31", "X433.32", "X433.6", "X433.8")
phecodes <- c("X286.7",
              "X286.8",
              "X286.81",
              "X395",
              "X396",
              "X401",
              "X411.2",
              "X411.4",
              "X415",
              "X425",
              "X425.1",
              "X426",
              "X427.2",
              "X430",
              "X433",
              "X433.1",
              "X433.11",
              "X433.12",
              "X433.2",
              "X433.21",
              "X433.3",
              "X433.31",
              "X433.5",
              "X433.6",
              "X433.8",
              "X440",
              "X443.9",
              "X452.2")
              

#phecodes <- Summary_phenotypes_merge_anno_common$phecode
#calculate culmulative incidence
for (i in unique(phecodes)) {
  print(i)
  phecode_select <- i
  phenotypes_long_AA <- phenotypes_long %>% #change to relatedness or nonrelatedness files
    #dplyr::filter(phecode == "X433.3") %>% #Cerebral ischemia for X433; X433.3; #Cardiomyopathy for X425; X425.1
    dplyr::filter(phecode %in% phecode_select) %>% #select of enriched phecode set
    dplyr::filter(phenotype == 1)
  AA_ID <- unique(paste0("PT", "", phenotypes_long_AA$id))
  print(length(AA_ID))
  
  phenotypes_long_AA <- phenotypes_long %>% #change to relatedness or nonrelatedness files
    #dplyr::filter(phecode == "X433.3") %>% #Cerebral ischemia for X433; #Cardiomyopathy for X425
    dplyr::filter(phecode %in% phecode_select) %>% #select of enriched phecode set
    group_by(id) %>%
    dplyr::filter(sum(phenotype) == 0) %>%
    ungroup()
  AA_WO <- unique(paste0("PT", "", phenotypes_long_AA$id))
  print(length(AA_WO))
  
  dat_RECUR_ICD109$disease <- ifelse(dat_RECUR_ICD109$PT_ID %in% AA_ID, 1, 0)
  dat_RECUR_ICD109_select <- dat_RECUR_ICD109[dat_RECUR_ICD109$PT_ID %in% AA_ID, ] %>%
    #dplyr::filter(phecode == "X433.3") %>%
    dplyr::filter(phecode %in% phecode_select) %>% #select of enriched phecode set
    group_by(PT_ID) %>%
    arrange(desc(ENC_DT)) %>%
    slice(n()) %>%
    ungroup()
  dat_RECUR_ICD109_select <- data.frame(dat_RECUR_ICD109_select)
  
  dat_RECUR_ICD109_WO <- dat_RECUR_ICD109[dat_RECUR_ICD109$PT_ID %in% AA_WO, ] %>%
    dplyr::group_by(PT_ID) %>%
    dplyr::arrange(ENC_DT) %>%
    slice(n()) %>%
    #mutate(disease = 0) %>%
    ungroup()
  dat_RECUR_ICD109_WO <- data.frame(dat_RECUR_ICD109_WO)
  GNSIS <- rbind(dat_RECUR_ICD109_select, dat_RECUR_ICD109_WO)
  print(dim(GNSIS))
  #remove spikein cases
  spikein_ID <- c(1117592, 1202698, 678791, 726999, 865829, 1360022, 1154976, 121012, 1346346, 242673, 259508,
                  377514, 527399, 635221, 753905, 865414, 928617)
  spikein_ID <- paste0("PT", spikein_ID)
  GNSIS <- GNSIS[!GNSIS$PT_ID %in% spikein_ID, ] #change from 7834 to 7827
  print(table(GNSIS$group, GNSIS$disease))
  
  withmutation <- nrow(GNSIS[GNSIS$group == "case" & GNSIS$disease == 1, ])
  withoutmutation <- nrow(GNSIS[GNSIS$group == "control" & GNSIS$disease == 1, ])
  IR_withmutation <- nrow(GNSIS[GNSIS$group == "case" & GNSIS$disease == 1, ])/sum(GNSIS[GNSIS$group == "case", ]$timediff_occur_censor)*1000  
  IR_withoutmutation <- nrow(GNSIS[GNSIS$group == "control" & GNSIS$disease == 1, ])/sum(GNSIS[GNSIS$group == "control", ]$timediff_occur_censor)*1000  
  print(dim(GNSIS))
  #Calculation of incidence rate ratios and 95% confidence intervals
  #Reference for all calculations: Rothman KJ, Greenland S. Modern epidemiology. 3rd ed.
  #Philadelphia: Lippincott Williams & Wilkins, 2008.
  IR <- IR_withmutation/IR_withoutmutation
  RD <- IR_withmutation - IR_withoutmutation
  
  SD <- (1/(nrow(GNSIS[GNSIS$group == "case" & GNSIS$disease == 1, ]) + 0.01) + 1/(nrow(GNSIS[GNSIS$group == "control" & GNSIS$disease == 1, ]) + 0.01))^0.5
  
  upper <- exp(log(IR) + 1.96*SD) 
  lower <- exp(log(IR) - 1.96*SD) 
  
  LIST[[i]] <- c(IR_withmutation, IR_withoutmutation, IR, RD, SD, lower, upper, withmutation, withoutmutation)
}

IR_summary <- data.frame(do.call(rbind,LIST))
head(IR_summary)
colnames(IR_summary) <- c("IR_withmutation", "IR_withoutmutation", "IR", "RD", "SD", "lower", "upper", "withmutation", "withoutmutation")
#colnames(IR_summary) <- c("IR_withmutation", "IR_withoutmutation", "IR", "RD", "SD", "lower", "upper")
IR_summary$phecode <- rownames(IR_summary)

save(IR_summary, file = "IR_summary_discovery_age.RData", version = 2)
save(IR_summary, file = "IR_summary_discovery_common.RData", version = 2)
write.table(IR_summary, file = "IR_summary_discovery_age.txt", col.names = T, row.names = F, quote = F)

save(IR_summary, file = "IR_summary_replication_age.RData", version = 2)
save(IR_summary, file = "IR_summary_replication_common.RData", version = 2)
save(IR_summary, file = "IR_summary_replication_rare.RData", version = 2)
write.table(IR_summary, file = "IR_summary_replication_age.txt", col.names = T, row.names = F, quote = F)


#calculate 95%CI of incidence rate
GNSIS <- NULL
LIST <- NULL
SD <- NULL
IR_withmutation <- NULL
IR_withoutmutation <- NULL
RD <- NULL
upper <- NULL
lower <- NULL
withmutation <- NULL
withoutmutation <- NULL

for (i in unique(phecodes)) {
  print(i)
  phecode_select <- i
  phenotypes_long_AA <- phenotypes_long %>% #change to relatedness or nonrelatedness files
    #dplyr::filter(phecode == "X433.3") %>% #Cerebral ischemia for X433; X433.3; #Cardiomyopathy for X425; X425.1
    dplyr::filter(phecode %in% phecode_select) %>% #select of enriched phecode set
    dplyr::filter(phenotype == 1)
  AA_ID <- unique(paste0("PT", "", phenotypes_long_AA$id))
  print(length(AA_ID))
  
  phenotypes_long_AA <- phenotypes_long %>% #change to relatedness or nonrelatedness files
    #dplyr::filter(phecode == "X433.3") %>% #Cerebral ischemia for X433; #Cardiomyopathy for X425
    dplyr::filter(phecode %in% phecode_select) %>% #select of enriched phecode set
    group_by(id) %>%
    dplyr::filter(sum(phenotype) == 0) %>%
    ungroup()
  AA_WO <- unique(paste0("PT", "", phenotypes_long_AA$id))
  print(length(AA_WO))
  
  dat_RECUR_ICD109$disease <- ifelse(dat_RECUR_ICD109$PT_ID %in% AA_ID, 1, 0)
  dat_RECUR_ICD109_select <- dat_RECUR_ICD109[dat_RECUR_ICD109$PT_ID %in% AA_ID, ] %>%
    #dplyr::filter(phecode == "X433.3") %>%
    dplyr::filter(phecode %in% phecode_select) %>% #select of enriched phecode set
    group_by(PT_ID) %>%
    arrange(desc(ENC_DT)) %>%
    slice(n()) %>%
    ungroup()
  dat_RECUR_ICD109_select <- data.frame(dat_RECUR_ICD109_select)
  
  dat_RECUR_ICD109_WO <- dat_RECUR_ICD109[dat_RECUR_ICD109$PT_ID %in% AA_WO, ] %>%
    dplyr::group_by(PT_ID) %>%
    dplyr::arrange(ENC_DT) %>%
    slice(n()) %>%
    #mutate(disease = 0) %>%
    ungroup()
  dat_RECUR_ICD109_WO <- data.frame(dat_RECUR_ICD109_WO)
  GNSIS <- rbind(dat_RECUR_ICD109_select, dat_RECUR_ICD109_WO)
  print(dim(GNSIS))
  #remove spikein cases
  spikein_ID <- c(1117592, 1202698, 678791, 726999, 865829, 1360022, 1154976, 121012, 1346346, 242673, 259508,
                  377514, 527399, 635221, 753905, 865414, 928617)
  spikein_ID <- paste0("PT", spikein_ID)
  GNSIS <- GNSIS[!GNSIS$PT_ID %in% spikein_ID, ] #change from 7834 to 7827
  print(table(GNSIS$group, GNSIS$disease))
  
  withmutation <- nrow(GNSIS[GNSIS$group == "case" & GNSIS$disease == 1, ])
  withoutmutation <- nrow(GNSIS[GNSIS$group == "control" & GNSIS$disease == 1, ])
  
  withmutation_personyear <- sum(GNSIS[GNSIS$group == "case", ]$timediff_occur_censor)
  withoutmutation_personyear <- sum(GNSIS[GNSIS$group == "control", ]$timediff_occur_censor)
  IR_withmutation <- withmutation/withmutation_personyear*1000  
  IR_withoutmutation <- withoutmutation/withoutmutation_personyear*1000  
  RD <- IR_withmutation - IR_withoutmutation
  RRR <- RD/IR_withoutmutation
  print(dim(GNSIS))
  #Calculation of incidence rate and 95% confidence intervals
  #https://stats.stackexchange.com/questions/269613/how-to-calculate-95-ci-for-cumulative-incidence-rate-using-r
  #https://stats.stackexchange.com/questions/301777/how-to-calculate-confidence-interval-of-incidence-rate-under-the-poisson-distrib
  #use a poisson GLM to do maximum likelihood under the assumption of independent interarrival times:
  fit_withmutation <- glm(withmutation ~ offset(log(withmutation_personyear)), family=poisson)
  lower_withmutation <- exp(confint(fit_withmutation))[1]*1000
  upper_withmutation <- exp(confint(fit_withmutation))[2]*1000
  
  fit_withoutmutation <- glm(withoutmutation ~ offset(log(withoutmutation_personyear)), family=poisson)
  lower_withoutmutation <- exp(confint(fit_withoutmutation))[1]*1000
  upper_withoutmutation <- exp(confint(fit_withoutmutation))[2]*1000
  
  LIST[[i]] <- c(withmutation, withoutmutation, withmutation_personyear, withoutmutation_personyear, IR_withmutation, 
                 lower_withmutation, upper_withmutation, IR_withoutmutation, lower_withoutmutation, upper_withoutmutation, 
                 RD, RRR)
}

IR_summary <- data.frame(do.call(rbind,LIST))
head(IR_summary)
colnames(IR_summary) <- c("withmutation", "withoutmutation", "withmutation_personyear", "withoutmutation_personyear", "IR_withmutation", 
                 "lower_withmutation", "upper_withmutation", "IR_withoutmutation", "lower_withoutmutation", "upper_withoutmutation", 
                 "RD", "RRR")
IR_summary$phecode <- rownames(IR_summary)
save(IR_summary, file = "IR_summary_discovery_95CI.RData", version = 2)
write.table(IR_summary, file = "IR_summary_discovery_95CI.txt", col.names = T, row.names = F, quote = F)

save(IR_summary, file = "IR_summary_replication_95CI.RData", version = 2)
write.table(IR_summary, file = "IR_summary_replication_95CI.txt", col.names = T, row.names = F, quote = F)

#https://www.statsdirect.com/help/rates/compare_crude_incidence_rates.html

#for v1.2
phecode_mapping <- read.csv("Phecode_map_v1_2_icd10cm_beta.csv", header = T, stringsAsFactors = F)
phecode_mapping_unique <- phecode_mapping %>%
  dplyr::select(phecode, phecode_str, exclude_range, exclude_name) %>%
  distinct(phecode, .keep_all = T)  #1755*4
str(phecode_mapping_unique)

phecode_mapping_unique$phecode <- paste0("X", phecode_mapping_unique$phecode)
load("IR_summary_discovery_common.RData")
load("IR_summary_replication_common.RData")
IR_summary <- merge(IR_summary, phecode_mapping_unique, by = "phecode")
IR_summary$RRR <- IR_summary$RD/IR_summary$IR_withoutmutation
save(IR_summary, file = "IR_summary_discovery_common_1.RData", version = 2)
save(IR_summary, file = "IR_summary_replication_common_1.RData", version = 2)


#####################
GNSIS <- GNSIS %>% 
  mutate(timediff_recur_censor = 
           case_when(disease == 1 ~ (timediff_occur_censor - 0)*365.25,
                     disease == 0 ~ (timediff_occur_censor - 0)*365.25)) 
GNSIS <- GNSIS %>% 
  mutate(timediff_recur_censor_55yr =
           case_when(timediff_recur_censor < 55*365.25 ~ timediff_recur_censor,
                     timediff_recur_censor >= 55*365.25 ~ 55*365.25),
         event_recur_censor_55yr =
           case_when(disease == 0 ~ 0,
                     disease == 1 & timediff_recur_censor <= 55*365.25 ~1,
                     disease == 1 & timediff_recur_censor > 55*365.25 ~0),
         timediff_recur_censor_60yr =
           case_when(timediff_recur_censor < 60*365.25 ~ timediff_recur_censor,
                     timediff_recur_censor >= 60*365.25 ~ 60*365.25),
         event_recur_censor_60yr =
           case_when(disease == 0 ~ 0,
                     disease == 1 & timediff_recur_censor <= 60*365.25 ~1,
                     disease == 1 & timediff_recur_censor > 60*365.25 ~0),
         timediff_recur_censor_85yr =
           case_when(timediff_recur_censor < 85*365.25 ~ timediff_recur_censor,
                     timediff_recur_censor >= 85*365.25 ~ 85*365.25),
         event_recur_censor_85yr =
           case_when(disease == 0 ~ 0,
                     disease == 1 & timediff_recur_censor <= 85*365.25 ~1,
                     disease == 1 & timediff_recur_censor > 85*365.25 ~0),
         timediff_recur_censor_65yr =
           case_when(timediff_recur_censor < 65*365.25 ~ timediff_recur_censor,
                     timediff_recur_censor >= 65*365.25 ~ 65*365.25),
         event_recur_censor_65yr =
           case_when(disease == 0 ~ 0,
                     disease == 1 & timediff_recur_censor <= 65*365.25 ~1,
                     disease == 1 & timediff_recur_censor > 65*365.25 ~0)
  )

#optional
# GNSIS_select <- GNSIS %>%
#   filter(timediff_recur_censor_35yr >= 0)

GNSIS_select <- GNSIS
unique(GNSIS_select$PT_RACE)
GNSIS_select <- GNSIS_select[!GNSIS_select$PT_RACE %in% c("Unknown", "American Indian Or Alaska Native", "Native Hawaiian Or Other Pacific Islander"), ]
head(GNSIS_select)
ANN <- phecode_mapping_unique[phecode_mapping_unique$phecode == "X433.2", ]$phecode_str
ANN <- "circulatory system"
ANN <- "Cerebral ischemia"
ANN <- "Hypercoagulability"
ANN <- "Hematopoietic"
ANN <- "Cardiomyopathy"
q <- NULL
q <- ggsurvplot(
  #surv_fit(Surv(timediff_recur_censor_85yr/365.25, event_recur_censor_85yr)~as.factor(group), data = GNSIS_select), # survfit object with calculated statistics.
  surv_fit(Surv(timediff_recur_censor_85yr/365.25, event_recur_censor_85yr)~as.factor(group) + as.factor(PT_RACE), data = GNSIS_select), # survfit object with calculated statistics.
  
  #surv_fit(Surv(timediff_recur_censor_85yr/365.25, event_recur_censor_85yr)~as.factor(group) + as.factor(PT_SEX), data = GNSIS_select), # survfit object with calculated statistics.
  fun = function(x) {1-x},  #plot cumulative probability F(t) = 1 - S(t)
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  pval.coord = c(2, 0.15),
  pval.size = 8,
  conf.int = FALSE,       # show confidence intervals for 
  censor.shape="|", censor.size = 4,
  fontsize = 4, #font size in the table
  font.legend = c(14, "plain", "black"),
  # point estimates of survival curves.
  xlim = c(0, 85),         # present narrower X axis, but not affect
  ylim = c(0, 0.15),
  # survival estimates.
  xlab = "Time in years",   # customize X axis label.
  ylab = "Cumulative Incidence",
  break.time.by = 1,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  n.risk = FALSE,
  #legend.labs = c("w/ rare allele - female", "w/ rare allele - male", 
  #               "w/o rare alleles - female", "w/o rare allele - male"),    # change legend labels.
  # legend.labs = c("w/ rare allele - Asian", "w/ rare allele - African American", "w/ rare allele - Caucasian", 
  #                "w/o rare alleles - Asian", "w/o rare allele - African American", "w/o rare allele - Caucasian"),    # change legend labels.
  # 
  font.title    = c(14, "bold.italic", "darkgreen"),
  font.subtitle = c(14, "bold", "darkgreen"),
  font.caption  = c(14, "plain", "darkgreen"),
  font.x        = c(14, "bold.italic", "black"),
  font.y        = c(14, "bold.italic", "black"),
  font.xtickslab = c(14, "bold", "black"),
  font.ytickslab = c(14, "bold", "black")
) +
  labs(
    title    = "Cumulative probability for 85yrs",
    #subtitle = paste('Disease = ', "Circulatory system 20 phecodes"),
    #subtitle = paste('Disease = ', "Hematopoietic 7 phecodes"),
    #subtitle = paste('Disease = ', "Hypercoagulability 3 phecodes"),
    #subtitle = paste('Disease = ', "Cerebral ischemia 10 phecodes"),
    #subtitle = paste('Disease = ', "Cerebral ischemia (X433.3)"),
    #subtitle = paste('Disease = ', "Cerebrovascular disease (X433)"),
    subtitle = paste('Disease = ', "Cardiomyopathy(425.1)"),
    #subtitle = paste('Disease = ', "Primary/intrinsic cardiomyopathies (425.1)"),
    caption  = "Plotted with R survminer")
q$plot <- q$plot + ggplot2::annotate("text", x = 15, y = 0.10,# x and y coordinates of the text
                                     label = ANN, size = 5)

q$plot <- q$plot + 
  scale_x_continuous(breaks = sort(c(seq(0, 85, 5))))

#tiff(paste('D:/Geisinger/X19981_backup_10072020/Desktop/MS_project/survivalcurve_90K/', 'Cerebralvascular_90K_', 'X433_7genes', '_enriched_85yr', '.tiff', sep=''), units="in", width=15, height=10, res=600)
#tiff(paste('D:/Geisinger/X19981_backup_10072020/Desktop/MS_project/survivalcurve_90K/', 'Cerebralvascular_90K_', 'Circulatory_system_7genes', '_enriched_85yr', '.tiff', sep=''), units="in", width=15, height=10, res=600)
#tiff(paste('D:/Geisinger/X19981_backup_10072020/Desktop/MS_project/survivalcurve_90K/', 'Cerebralvascular_90K_', 'Cerebralvascular_6phecodes_7genes', '_enriched_85yr', '.tiff', sep=''), units="in", width=15, height=10, res=600)
#tiff(paste('D:/Geisinger/X19981_backup_10072020/Desktop/MS_project/survivalcurve_90K/', 'Cerebralvascular_90K_', 'Cerebralvascular_10phecodes_RACE_7genes', '_enriched_85yr', '.tiff', sep=''), units="in", width=15, height=10, res=600)
#tiff(paste('D:/Geisinger/X19981_backup_10072020/Desktop/MS_project/survivalcurve_90K/', 'Hypercoagulability_85K_', 'Hypercoagulability_3phecodes_RACE_7genes', '_enriched_85yr', '.tiff', sep=''), units="in", width=15, height=10, res=600)
#tiff(paste('D:/Geisinger/X19981_backup_10072020/Desktop/MS_project/survivalcurve_90K/', 'Hematopoietic_85K_', 'Hematopoietic_7phecodes_RACE_7genes', '_enriched_85yr', '.tiff', sep=''), units="in", width=15, height=10, res=600)
#tiff(paste('D:/Geisinger/X19981_backup_10072020/Desktop/MS_project/survivalcurve_90K/', 'Circulatorysystem_85K_', 'Circulatorysystem_20phecodes_RACE_7genes', '_enriched_85yr', '.tiff', sep=''), units="in", width=15, height=10, res=600)
tiff(paste('D:/Geisinger/X19981_backup_10072020/Desktop/MS_project/survivalcurve_90K_06062022/', 'Cardiomyopathy_90K_', 'Cardiomyopathy_1phecodes_RACE_7genes', '_enriched_85yr', '.tiff', sep=''), units="in", width=15, height=10, res=600)
print(q)
dev.off()

# Pairwise survdiff
library(survival)
library(survminer)
GNSIS_select$subgroup <- paste0(GNSIS_select$group, GNSIS_select$PT_RACE) 
res <- pairwise_survdiff(Surv(timediff_recur_censor_85yr/365.25, event_recur_censor_85yr)~subgroup, p.adjust.method = "none", data = GNSIS_select)
#res <- pairwise_survdiff(Surv(timediff_recur_censor_85yr/365.25, event_recur_censor_85yr)~group, p.adjust.method = "none", data = GNSIS_select)

res
# caseAsian caseBlack Or African American caseWhite controlAsian
# caseBlack Or African American    0.20432   -                             -         -           
#   caseWhite                        0.00245   0.00016                       -         -           
#   controlAsian                     0.16354   0.99739                       0.36442   -           
#   controlBlack Or African American 0.01835   0.36449                       0.00751   0.66334     
# controlWhite                     0.00040   2.7e-06                       0.01514   0.22507     
# controlBlack Or African American
# caseBlack Or African American    -                               
#   caseWhite                        -                               
#   controlAsian                     -                               
#   controlBlack Or African American -                               
#   controlWhite                     0.00024

# caseAsian caseBlack Or African American caseWhite controlAsian
# caseBlack Or African American    0.799     -                             -         -           
#   caseWhite                        0.748     0.840                         -         -           
#   controlAsian                     1.000     0.854                         0.688     -           
#   controlBlack Or African American 0.858     0.296                         0.138     0.821       
# controlWhite                     0.818     0.621                         0.011     0.789       
# controlBlack Or African American
# caseBlack Or African American    -                               
#   caseWhite                        -                               
#   controlAsian                     -                               
#   controlBlack Or African American -                               
#   controlWhite                     0.358                           

# caseAsian caseBlack Or African American caseWhite controlAsian
# caseBlack Or African American    0.487     -                             -         -           
#   caseWhite                        0.746     1.9e-05                       -         -           
#   controlAsian                     1.000     0.568                         0.794     -           
#   controlBlack Or African American 0.650     0.078                         0.157     0.744       
# controlWhite                     0.755     2.3e-06                       0.524     0.794       
# controlBlack Or African American
# caseBlack Or African American    -                               
#   caseWhite                        -                               
#   controlAsian                     -                               
#   controlBlack Or African American -                               
#   controlWhite                     0.094      

#specify the patients with gene of mutation
load("Individualgene_PTID_12032021.RData")
load("Individualgene_PTID_12032021_HTRA1_remove.RData") #HTRA1_remove

head(GNSIS_select)
#GNSIS_select <- data.frame(GNSIS_select)
GNSIS_select <- GNSIS
#GNSIS_select <- GNSIS_select[GNSIS_select$PT_RACE == "Black Or African American", ]
# GNSIS_select$subgroup <- "W/O Mut"
# GNSIS_select[GNSIS_select$PT_ID %in% COL4A1, ]$subgroup <- "COL4A1"
# GNSIS_select[GNSIS_select$PT_ID %in% NOTCH3, ]$subgroup <- "NOTCH3"
# GNSIS_select[GNSIS_select$PT_ID %in% HTRA1, ]$subgroup <- "HTRA1"
# GNSIS_select[GNSIS_select$PT_ID %in% GLA, ]$subgroup <- "GLA"
# GNSIS_select[GNSIS_select$PT_ID %in% CTC1, ]$subgroup <- "CTC1"
# GNSIS_select[GNSIS_select$PT_ID %in% NOTCH3, ]$subgroup <- "TREX1"
# 
# table(GNSIS_select$subgroup)
# COL4A1    CTC1     GLA   HTRA1  NOTCH3   NOTCH3 W/O Mut 
# 1698      35     174      58     585     195    5490 
# COL4A1    CTC1     GLA   HTRA1  NOTCH3   NOTCH3 W/O Mut 
# 1696      35     174      58     583     195    5483 (hematopoietic)

# COL4A1    CTC1     GLA   HTRA1  NOTCH3   NOTCH3 W/O Mut 
# 761      33     142      62     518     188    3410
# COL4A1    CTC1     GLA   HTRA1  NOTCH3   NOTCH3 W/O Mut 
# 761      33     142      62     518     188    3408 (85K using 90K enriched phecode set)
# COL4A1    CTC1     GLA   HTRA1  NOTCH3   NOTCH3 W/O Mut 
# 761      33     142      61     518     188    3400 (hematopoietic)

#white
# COL4A1    CTC1     GLA   HTRA1  NOTCH3   NOTCH3 W/O Mut 
# 1560      34     174      56     577     192    5241 
# COL4A1    CTC1     GLA   HTRA1  NOTCH3   TREX1 W/O Mut 
# 568      33     142      61     511     169    3036

#Black
# COL4A1    CTC1   HTRA1   TREX1 W/O Mut 
# 125       1       2       3     211 
# COL4A1     GLA  NOTCH3   NOTCH3 W/O Mut 
# 178       1       5      17     299 

GNSIS_select_COL4A1 <- GNSIS_select[GNSIS_select$PT_ID %in% COL4A1, ] #1574; 572; 125; 179
GNSIS_select_COL4A1$subgroup <- "COL4A1"
GNSIS_select_NOTCH3 <- GNSIS_select[GNSIS_select$PT_ID %in% NOTCH3, ] #580; 513; 0; 5
GNSIS_select_NOTCH3$subgroup <- "NOTCH3"
GNSIS_select_TREX1 <- GNSIS_select[GNSIS_select$PT_ID %in% TREX1, ] #192; 169; 3; 17
GNSIS_select_TREX1$subgroup <- "TREX1"
GNSIS_select_GLA <- GNSIS_select[GNSIS_select$PT_ID %in% GLA, ] #174; 142; 0; 1
GNSIS_select_GLA$subgroup <- "GLA"
GNSIS_select_HTRA1 <- GNSIS_select[GNSIS_select$PT_ID %in% HTRA1, ] #56; 61; 2; 0
GNSIS_select_HTRA1$subgroup <- "HTRA1"
GNSIS_select_CTC1 <- GNSIS_select[GNSIS_select$PT_ID %in% CTC1, ] #34; 33; 1; 0
GNSIS_select_CTC1$subgroup <- "CTC1"
PT_ID <- unique(c(COL4A1, NOTCH3, TREX1, GLA, HTRA1, CTC1)) #4463
GNSIS_select_WOMUT <- GNSIS_select[!(GNSIS_select$PT_ID %in% PT_ID), ] 
GNSIS_select_WOMUT$subgroup <- "W/O Mut"

#90K or 85K WHITE
listOfDataFrames <- list(GNSIS_select_WOMUT, 
                         GNSIS_select_COL4A1,
                         GNSIS_select_NOTCH3,
                         GNSIS_select_TREX1,
                         GNSIS_select_GLA,
                         GNSIS_select_HTRA1,
                         GNSIS_select_CTC1)
#90K BLACK
listOfDataFrames <- list(GNSIS_select_WOMUT, 
                         GNSIS_select_COL4A1,
                         GNSIS_select_TREX1,
                         GNSIS_select_HTRA1,
                         GNSIS_select_CTC1)
#85K BLACK
listOfDataFrames <- list(GNSIS_select_WOMUT, 
                         GNSIS_select_COL4A1,
                         GNSIS_select_NOTCH3,
                         GNSIS_select_TREX1,
                         GNSIS_select_GLA)

GNSIS_select <- do.call("rbind", listOfDataFrames) #7834 to 7851


GNSIS_select$COHORT <- "Discovery"
save(GNSIS_select, file = "D:/Geisinger/X19981_backup_10072020/Desktop/MS_project/GNSIS_select_white_gene_90K_cardiomyopathy2phecodes_remove10.RData", version = 2)
GNSIS_select$COHORT <- "Replication"
save(GNSIS_select, file = "GNSIS_select_white_gene_85K_duplicated_cardiomyopathy.RData", version = 2)

# load("GNSIS_select_black_gene_90K.RData")
# GNSIS_select_discovery <- GNSIS_select
# load("GNSIS_select_black_gene_85K.RData")
# GNSIS_select_replication <- GNSIS_select

#no racial stratification
load("GNSIS_select_ischemicstroke_90K.RData")  #8235
GNSIS_select_discovery <- GNSIS_select
load("GNSIS_select_ischemicstroke_85K.RData")  #5115
GNSIS_select_replication <- GNSIS_select
GNSIS_select <- rbind(GNSIS_select_discovery, GNSIS_select_replication)

Noncarrier <- as.character(GNSIS_select[GNSIS_select$group == "control", ]$PT_ID)
GNSIS_select$subgroup <- "W/ Mut"
GNSIS_select$subgroup <- ifelse(GNSIS_select$PT_ID %in% HTRA1, "HTRA1", GNSIS_select$subgroup)
GNSIS_select$subgroup <- ifelse(GNSIS_select$PT_ID %in% Noncarrier, "W/O Mut", GNSIS_select$subgroup)
unique(GNSIS_select$subgroup)
table(GNSIS_select$subgroup)

# HTRA1  W/ Mut W/O Mut 
# 120    4330    8900

GNSIS_select <- GNSIS_select[GNSIS_select$subgroup != "W/ Mut", ]

GNSIS_select$subgroup <- ifelse(GNSIS_select$subgroup == "W/O Mut", GNSIS_select$subgroup, 
                                ifelse(GNSIS_select$PT_ID %in% HTRA1_remove, GNSIS_select$subgroup, "remove"))
# HTRA1  remove W/O Mut  all race
# 104      16    8900
# HTRA1  remove W/O Mut  white only
# 103      14    8277 

GNSIS_select <- GNSIS_select[GNSIS_select$subgroup != "remove", ]

#with racial stratification
load("GNSIS_select_white_gene_90K_duplicated_stroke.RData")
GNSIS_select_discovery <- GNSIS_select
load("GNSIS_select_white_gene_85K_duplicated_stroke.RData")
GNSIS_select_replication <- GNSIS_select

# load("GNSIS_select_black_gene_90K_duplicated_hypercoagulability.RData")
# GNSIS_select_discovery <- GNSIS_select
# load("GNSIS_select_black_gene_85K_duplicated_hypercoagulability.RData")
# GNSIS_select_replication <- GNSIS_select


load("D:/Geisinger/X19981_backup_10072020/Desktop/MS_project/GNSIS_select_white_gene_90K_duplicated_cardiomyopathy.RData")
GNSIS_select_discovery <- GNSIS_select
load("D:/Geisinger/X19981_backup_10072020/Desktop/MS_project/GNSIS_select_white_gene_85K_duplicated_cardiomyopathy.RData")
GNSIS_select_replication <- GNSIS_select

GNSIS_select <- rbind(GNSIS_select_discovery, GNSIS_select_replication)
GNSIS_select_WHITE <- GNSIS_select
GNSIS_select_BLACK <- GNSIS_select
GNSIS_select_BLACKWHITE <- rbind(GNSIS_select_WHITE, GNSIS_select_BLACK)

GNSIS_select <- GNSIS_select_BLACKWHITE[GNSIS_select_BLACKWHITE$COHORT == "Discovery", ]
GNSIS_select <- GNSIS_select_BLACKWHITE[GNSIS_select_BLACKWHITE$COHORT == "Replication", ]

table(GNSIS_select$subgroup, GNSIS_select$COHORT)
WHITE_ID <- unique(GNSIS_select[GNSIS_select$subgroup != "W/O Mut", ]$PT_ID) #4077
#save(WHITE_ID, file = "WHITE_ID.RData", version = 2)

BLACK_ID <- unique(GNSIS_select[GNSIS_select$subgroup != "W/O Mut", ]$PT_ID) #332
#save(BLACK_ID, file = "BLACK_ID.RData", version = 2)

head(GNSIS_select)
unique(GNSIS_select$subgroup)
unique(GNSIS_select$PT_SEX)
GNSIS_select_Female <- GNSIS_select[GNSIS_select$PT_SEX == "Female", ]
GNSIS_select_Male <- GNSIS_select[GNSIS_select$PT_SEX == "Male", ]
#table(GNSIS_select$subgroup, GNSIS_select$COHORT)
GNSIS_select <- GNSIS_select_Female
GNSIS_select <- GNSIS_select_Male
GNSIS_select$subgroup[GNSIS_select$subgroup != "W/O Mut"] <- "W/ Mut"
unique(GNSIS_select$subgroup)
GNSIS_select <- GNSIS_select[GNSIS_select$subgroup %in% c("COL4A1", "W/O Mut"), ]
GNSIS_select <- GNSIS_select[GNSIS_select$subgroup %in% c("NOTCH3", "W/O Mut"), ]
GNSIS_select <- GNSIS_select[GNSIS_select$subgroup %in% c("TREX1", "W/O Mut"), ]
GNSIS_select <- GNSIS_select[GNSIS_select$subgroup %in% c("GLA", "W/O Mut"), ]
GNSIS_select <- GNSIS_select[GNSIS_select$subgroup %in% c("CTC1", "W/O Mut"), ]
GNSIS_select <- GNSIS_select[GNSIS_select$subgroup %in% c("HTRA1", "W/O Mut"), ]

GNSIS_select <- GNSIS_select[GNSIS_select$COHORT == "Discovery", ]
GNSIS_select <- GNSIS_select[GNSIS_select$COHORT == "Replication", ]


#Remove HTRA1 p.Gln151Lys
GNSIS_select <- GNSIS_select[!(GNSIS_select$PT_ID %in% HTRA1_remove), ] #change from 8394 to 8291
#include HTRA1 p.Gln151Lys
GNSIS_select <- GNSIS_select[GNSIS_select$PT_ID %in% HTRA1_remove, ] #change from 8394 to 8291

#GNSIS_select <- GNSIS_select[GNSIS_select$COHORT == "Replication", ]
#ANN <- phecode_mapping_unique[phecode_mapping_unique$phecodeX == "X433.3", ]$phecode_str
ANN <- "circulatory system"
ANN <- "ischemic stroke (6 PheCodes)"
ANN <- "hypercoagulability (3 PheCodes)"
ANN <- "cardiomyopathy (2 PheCodes)"
library(survival)
library(survminer)

q <- NULL
q <- ggsurvplot(
 # surv_fit(Surv(timediff_recur_censor_85yr/365.25, event_recur_censor_85yr)~as.factor(subgroup) + as.factor(COHORT) + as.factor(PT_SEX), data = GNSIS_select), # survfit object with calculated statistics.
  #surv_fit(Surv(timediff_recur_censor_85yr/365.25, event_recur_censor_85yr)~as.factor(subgroup) + as.factor(COHORT), data = GNSIS_select), # survfit object with calculated statistics.
  #surv_fit(Surv(timediff_recur_censor_85yr/365.25, event_recur_censor_85yr)~as.factor(subgroup) + as.factor(PT_SEX), data = GNSIS_select), # survfit object with calculated statistics.
  surv_fit(Surv(timediff_recur_censor_85yr/365.25, event_recur_censor_85yr)~as.factor(subgroup), data = GNSIS_select), # survfit object with calculated statistics.
  #surv_fit(Surv(timediff_recur_censor_85yr/365.25, event_recur_censor_85yr)~as.factor(subgroup) + as.factor(PT_SEX), data = GNSIS_select), # survfit object with calculated statistics.
  fun = function(x) {1-x},  #plot cumulative probability F(t) = 1 - S(t)
  risk.table = TRUE,  # show risk table.
  pval = FALSE,             # show p-value of log-rank test.
  pval.coord = c(2, 0.4),
  pval.size = 8,
  conf.int = FALSE,       # show confidence intervals for 
  #censor.shape="|", censor.size = 4,
  fontsize = 4, #font size in the table
  font.legend = c(14, "plain", "black"),
  # point estimates of survival curves.
  xlim = c(0, 80),         # present narrower X axis, but not affect
  ylim = c(0, 0.35),
  # survival estimates.
  xlab = "Time in years",   # customize X axis label.
  ylab = "Cumulative Incidence",
  break.time.by = 5,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  n.risk = FALSE,
  #legend.labs = c("Carrier-Discovery", "Carrier-Replication", "NonCarrier-Discovery", "NonCarrier-Replication"),
  #legend.labs = c("Carrier-Discovery-Female", "Carrier-Discovery-Male", "Carrier-Replication-Female", "Carrier-Replication-Male", "NonCarrier-Discovery-Female", "NonCarrier-Discovery-Male", "NonCarrier-Replication-Female", "NonCarrier-Replication-Male"),
  #legend.labs = c("NOTCH3-Discovery", "NOTCH3-Replication", "W/O Mut-Discovery", "W/O Mut-Replication"),
  #legend.labs = c("TREX1-Discovery", "TREX1-Replication", "W/O Mut-Discovery", "W/O Mut-Replication"),
  #legend.labs = c("GLA-Discovery", "GLA-Replication", "W/O Mut-Discovery", "W/O Mut-Replication"),
  #legend.labs = c("CTC1-Discovery", "CTC1-Replication", "W/O Mut-Discovery", "W/O Mut-Replication"),
  #legend.labs = c("HTRA1-Discovery", "HTRA1-Replication", "W/O Mut-Discovery", "W/O Mut-Replication"),
  #legend.labs = c("COL4A1-Female", "COL4A1-Male", "NonCarrier-Female", "NonCarrier-Male"),
  #legend.labs = c("HTRA1-Female", "HTRA1-Male", "NonCarrier-Female", "NonCarrier-Male"),
  #legend.labs = c("COL4A1", "CTC1", "GLA", "HTRA1", "NOTCH3", "TREX1", "Noncarrier"), # change legend labels.
  #legend.labs = c("COL4A1", "CTC1", "HTRA1", "TREX1", "W/O Mut"),  #90K
  #legend.labs = c("COL4A1", "GLA", "NOTCH3", "TREX1", "W/O Mut"),  #85K
  #legend.labs = c("TREX1-Female", "TREX1-Male", "NonCarrier-Female", "NonCarrier-Male"),
  #legend.labs = c("NOTCH3-Female", "NOTCH3-Male", "NonCarrier-Female", "NonCarrier-Male"),
  legend.labs = c("Carrier_HTRA1", "NonCarrier"),
  font.title    = c(14, "bold.italic", "darkgreen"),
  font.subtitle = c(14, "bold", "darkgreen"),
  font.caption  = c(14, "plain", "darkgreen"),
  font.x        = c(14, "bold.italic", "black"),
  font.y        = c(14, "bold.italic", "black"),
  font.xtickslab = c(14, "bold", "black"),
  font.ytickslab = c(14, "bold", "black")
) +
  labs(
    title    = "Cumulative probability for 85yrs",
    #subtitle = paste('Disease = ', "Circulatory system", " Top enriched 20 phecodes"),
    #subtitle = paste('Disease = ', "hematopoietic", " Top enriched 7 phecodes"),
    #subtitle = paste('Disease = ', "hypercoagulability", " 3 phecodes", "for the replication cohort with European or African Ancestry"),
    #subtitle = paste('Disease = ', "ischemic stroke", " 6 PheCodes", " Carriers vs NonCarriers in either discovery or replication cohorts stratified by sex"),
    #subtitle = paste('Disease = ', "ischemic stroke", " 6 PheCodes", " Carriers vs NonCarriers in either discovery or replication Female cohorts "),
    subtitle = paste('Disease = ', "ischemic stroke", " 6 phecodes", " in both discovery and replication cohorts"),
    #subtitle = paste('Disease = ', "Cerebral ischemia (X433.3)"),
    #subtitle = paste('Disease = ', "Cerebrovascular disease (X433)"),
    #subtitle = paste('Disease = ', "cardiomyopathy", " 2 phecodes", "for the replication cohort with European or African Ancestry"),
    #subtitle = paste('Disease = ', "Primary/intrinsic cardiomyopathies (425.1)"),
    caption  = "Plotted with R survminer")
q$plot <- q$plot + ggplot2::annotate("text", x = 15, y = 0.50,# x and y coordinates of the text
                                     label = ANN, size = 5)

# q$plot <- q$plot + 
#   scale_x_continuous(breaks = sort(c(seq(0, 85, 5))))

tiff(paste('D:/Geisinger/X19981_backup_10072020/Desktop/MS_project/survivalcurve_90K/', 'ischemicstroke_6phecodes', '_system_7genes_sep', '_enriched_85yr_discoveryreplication_HTRA1_test_white', '.tiff', sep=''), units="in", width=15, height=10, res=600)
print(q)
dev.off()

GNSIS_select$subgroup <- paste0(GNSIS_select$subgroup, GNSIS_select$PT_SEX)
#GNSIS_select$subgroup <- paste0(GNSIS_select$subgroup, GNSIS_select$COHORT)
#GNSIS_select$subgroup <- paste0(GNSIS_select$subgroup, GNSIS_select$COHORT, GNSIS_select$PT_SEX)
table(GNSIS_select$subgroup)
res <- pairwise_survdiff(Surv(timediff_recur_censor_85yr/365.25, event_recur_censor_85yr)~subgroup, p.adjust.method = "none", data = GNSIS_select)
#res <- pairwise_survdiff(Surv(timediff_recur_censor_80yr/365.25, event_recur_censor_80yr)~group, p.adjust.method = "none", data = GNSIS_select)

res


cox_65yrrecur_COL4A1 <- coxph(Surv(timediff_recur_censor_65yr/365.25, event_recur_censor_65yr)~as.factor(group) + as.factor(PT_SEX), data = GNSIS_select)

coxtidy_65yrrecur_COL4A1 <- tidy(cox_65yrrecur_COL4A1, exponentiate = T, conf.int = T, conf.level = 0.95)
coxtidy_65yrrecur_COL4A1$Onset <- "Early"
coxtidy_65yrrecur_COL4A1$Category <- "Cerebral ischemia"
coxtidy_65yrrecur_COL4A1$Sex <- "Combined"
coxtidy_65yrrecur_COL4A1 <- data.frame(coxtidy_65yrrecur_COL4A1)

cox_85yrrecur_COL4A1 <- coxph(Surv(timediff_recur_censor_85yr/365.25, event_recur_censor_85yr)~as.factor(group) + as.factor(PT_SEX), data = GNSIS_select)

coxtidy_85yrrecur_COL4A1 <- tidy(cox_85yrrecur_COL4A1, exponentiate = T, conf.int = T, conf.level = 0.95)
coxtidy_85yrrecur_COL4A1$Onset <- "Lifetime"
coxtidy_85yrrecur_COL4A1$Category <- "Cerebral ischemia"
coxtidy_85yrrecur_COL4A1$Sex <- "Combined"
coxtidy_85yrrecur_COL4A1 <- data.frame(coxtidy_85yrrecur_COL4A1)

fooCox_COL4A1_90K_WHITE <- rbind(coxtidy_65yrrecur_COL4A1, coxtidy_85yrrecur_COL4A1)

GNSIS_select_COL4A1_FEMALE <- GNSIS_select[GNSIS_select$PT_SEX == "Female", ]

cox_65yrrecur_COL4A1 <- coxph(Surv(timediff_recur_censor_65yr/365.25, event_recur_censor_65yr)~as.factor(group), data = GNSIS_select_COL4A1_FEMALE)

coxtidy_65yrrecur_COL4A1 <- tidy(cox_65yrrecur_COL4A1, exponentiate = T, conf.int = T, conf.level = 0.95)
coxtidy_65yrrecur_COL4A1$Onset <- "Early"
coxtidy_65yrrecur_COL4A1$Category <- "Cerebral ischemia"
coxtidy_65yrrecur_COL4A1$Sex <- "Female"
coxtidy_65yrrecur_COL4A1 <- data.frame(coxtidy_65yrrecur_COL4A1)

cox_85yrrecur_COL4A1 <- coxph(Surv(timediff_recur_censor_85yr/365.25, event_recur_censor_85yr)~as.factor(group), data = GNSIS_select_COL4A1_FEMALE)
coxtidy_85yrrecur_COL4A1 <- tidy(cox_85yrrecur_COL4A1, exponentiate = T, conf.int = T, conf.level = 0.95)
coxtidy_85yrrecur_COL4A1$Onset <- "Lifetime"
coxtidy_85yrrecur_COL4A1$Category <- "Cerebral ischemia"
coxtidy_85yrrecur_COL4A1$Sex <- "Female"
coxtidy_85yrrecur_COL4A1 <- data.frame(coxtidy_85yrrecur_COL4A1)

fooCox_COL4A1_90K_WHITE_FEMALE <- rbind(coxtidy_65yrrecur_COL4A1, coxtidy_85yrrecur_COL4A1)

GNSIS_select_COL4A1_MALE <- GNSIS_select[GNSIS_select$PT_SEX == "Male", ]

cox_65yrrecur_COL4A1 <- coxph(Surv(timediff_recur_censor_65yr/365.25, event_recur_censor_65yr)~as.factor(group), data = GNSIS_select_COL4A1_MALE)

coxtidy_65yrrecur_COL4A1 <- tidy(cox_65yrrecur_COL4A1, exponentiate = T, conf.int = T, conf.level = 0.95)
coxtidy_65yrrecur_COL4A1$Onset <- "Early"
coxtidy_65yrrecur_COL4A1$Category <- "Cerebral ischemia"
coxtidy_65yrrecur_COL4A1$Sex <- "Male"
coxtidy_65yrrecur_COL4A1 <- data.frame(coxtidy_65yrrecur_COL4A1)

cox_85yrrecur_COL4A1 <- coxph(Surv(timediff_recur_censor_85yr/365.25, event_recur_censor_85yr)~as.factor(group), data = GNSIS_select_COL4A1_MALE)

coxtidy_85yrrecur_COL4A1 <- tidy(cox_85yrrecur_COL4A1, exponentiate = T, conf.int = T, conf.level = 0.95)
coxtidy_85yrrecur_COL4A1$Onset <- "Lifetime"
coxtidy_85yrrecur_COL4A1$Category <- "Cerebral ischemia"
coxtidy_85yrrecur_COL4A1$Sex <- "Male"
coxtidy_85yrrecur_COL4A1 <- data.frame(coxtidy_85yrrecur_COL4A1)

fooCox_COL4A1_90K_WHITE_MALE <- rbind(coxtidy_65yrrecur_COL4A1, coxtidy_85yrrecur_COL4A1)

fooCox_COL4A1_90K_WHITE_MALE_FEMALE <- rbind(fooCox_COL4A1_90K_WHITE_MALE, fooCox_COL4A1_90K_WHITE_FEMALE)

fooCox_COL4A1_90K_WHITE <- rbind(fooCox_COL4A1_90K_WHITE, fooCox_COL4A1_90K_WHITE_MALE_FEMALE)
colnames(fooCox_COL4A1_90K_WHITE)
fooCox_COL4A1_90K_WHITE$Gene <- "HTRA1"
save(fooCox_COL4A1_90K_WHITE, file = "fooCox_stroke_90K85K_WHITE_HTRA1.RData", version = 2)

fooCox_COL4A1_90K_WHITE$estimate_95CI <- paste0(round(1/fooCox_COL4A1_90K_WHITE$estimate, 3), "(", round(1/fooCox_COL4A1_90K_WHITE$conf.high, 3), "-", round(1/fooCox_COL4A1_90K_WHITE$conf.low, 3), ")")
#fooCox_COL4A1_90K_WHITE$term[fooCox_COL4A1_90K_WHITE$term == 'as.factor(group)control'] <- "genetics"
fooCox_COL4A1_90K_WHITE_select <- fooCox_COL4A1_90K_WHITE %>%
  filter(term == 'as.factor(group)control')
fooCox_COL4A1_90K_WHITE_select$p.value <- round(fooCox_COL4A1_90K_WHITE_select$p.value, 3)

head(fooCox_COL4A1_90K_WHITE_select)
#create a forest plot
tabletext <- fooCox_COL4A1_90K_WHITE_select %>%
  dplyr::select(Gene, Sex, Onset, p.value, estimate_95CI)

head(tabletext)

header <- c("Gene", "Sex", "Onset", "P value", "Hazard Ratio (95%CI)")
tabletext <- rbind(header, tabletext)
forestdata <- fooCox_COL4A1_90K_WHITE_select %>%
  dplyr::select(estimate, conf.high, conf.low,)
forestdata$estimate <- 1/forestdata$estimate
forestdata$conf.low <- 1/forestdata$conf.low
forestdata$conf.high <- 1/forestdata$conf.high
head(forestdata)
header2 <- c(1, 1, 1)
forestdata <- rbind(header2, forestdata)
forestdata_new <- forestdata[c(1:2, 3:7),] #43 or 31 or 13
colnames(forestdata_new) <- c("mean", "lower", "upper")
library(forestplot)
library(dplyr)
tabletext_new <- as.matrix(tabletext[c(1:2, 3:7),])
tabletext_new <- as.matrix(tabletext)
dim(tabletext_new)

forestdata_new %>% 
  forestplot(labeltext = tabletext_new, 
             graph.pos = 3,
             #is.summary = c(rep(TRUE, 2), rep(FALSE, 8), TRUE),
             clip = c(0.1, 6), 
             xlog = TRUE, 
             col = fpColors(box = "royalblue",
                            line = "darkblue",
                            summary = "royalblue"),
             xlab = "estimate(95%CI)",
             txt_gp = fpTxtGp(label = list(gpar(fontfamily = "Arial"),
                                           gpar(fontfamily = "",
                                                col = "black")),
                              ticks = gpar(fontfamily = "", cex = 1.0),
                              xlab  = gpar(fontfamily = "Arial", cex = 1.0))
  )

tiff("Rplot_forestplot_multivariate_coxregression_disease_category_9085K_stroke_WHITE_combined_65yr_limited_HTRA1.tiff", units = "in", width = 8, height = 4, res = 300) #height change to 8
P <- forestdata_new %>% 
  forestplot(labeltext = tabletext_new, 
             graph.pos = 3,
             #is.summary = c(rep(TRUE, 2), rep(FALSE, 8), TRUE),
             clip = c(0.1, 6), 
             xlog = TRUE, 
             col = fpColors(box = "royalblue",
                            line = "darkblue",
                            summary = "royalblue"),
             xlab = "estimate(95%CI)",
             xlab = "estimate(95%CI)",
             txt_gp = fpTxtGp(label = list(gpar(fontfamily = "Arial"),
                                           gpar(fontfamily = "",
                                                col = "660000")),
                              ticks = gpar(fontfamily = "Arial", cex = 1.0),
                              xlab  = gpar(fontfamily = "Arial", cex = 1.0))
  )
#grid::grid.text("STROKE ~ having rare variants + (Sex) + PCs(1-5)", .5, 0.99, gp=gpar(cex=1.5))
dev.off()


#example code
library(forestplot)
library(dplyr)
# Cochrane data from the 'rmeta'-package
cochrane_from_rmeta <- structure(list(mean  = c(NA, NA, 0.578, 0.165, 0.246, 0.700, 0.348, 0.139, 1.017, NA, 0.531), 
                                      lower = c(NA, NA, 0.372, 0.018, 0.072, 0.333, 0.083, 0.016, 0.365, NA, 0.386),
                                      upper = c(NA, NA, 0.898, 1.517, 0.833, 1.474, 1.455, 1.209, 2.831, NA, 0.731)),
                                 .Names = c("mean", "lower", "upper"), 
                                 row.names = c(NA, -11L), 
                                 class = "data.frame")
dim(cochrane_from_rmeta)
tabletext <- cbind(c("", "Study", "Auckland", "Block", "Doran", "Gamsu", "Morrison", "Papageorgiou", "Tauesch", NA, "Summary"),
                   c("Deaths", "(steroid)", "36", "1", "4", "14", "3", "1", "8", NA, NA),
                   c("Deaths", "(placebo)", "60", "5", "11", "20", "7", "7", "10", NA, NA),
                   c("", "OR", "0.58", "0.16", "0.25", "0.70", "0.35", "0.14", "1.02", NA, "0.53"))
dim(tabletext)
head(tabletext)
cochrane_from_rmeta %>% 
  forestplot(labeltext = tabletext, 
             is.summary = c(rep(TRUE, 2), rep(FALSE, 8), TRUE),
             clip = c(0.1, 2.5), 
             xlog = TRUE, 
             col = fpColors(box = "royalblue",
                            line = "darkblue",
                            summary = "royalblue"))






















q <- NULL
q <- ggsurvplot(
  surv_fit(Surv(timediff_recur_censor_85yr/365.25, event_recur_censor_85yr)~as.factor(subgroup), data = GNSIS_select), # survfit object with calculated statistics.
  fun = function(x) {1-x},  #plot cumulative probability F(t) = 1 - S(t)
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  pval.coord = c(2, 0.15),
  pval.size = 8,
  conf.int = FALSE,       # show confidence intervals for 
  #censor.shape="|", censor.size = 4,
  fontsize = 4, #font size in the table
  font.legend = c(14, "plain", "black"),
  # point estimates of survival curves.
  xlim = c(0, 85),         # present narrower X axis, but not affect
  ylim = c(0, 0.15),
  # survival estimates.
  xlab = "Time in years",   # customize X axis label.
  ylab = "Cumulative Incidence",
  break.time.by = 1,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  n.risk = FALSE,
  #legend.labs = c("COL4A1", "CTC1", "GLA", "HTRA1", "NOTCH3", "TREX1", "W/O Mut"), # change legend labels.
  #legend.labs = c("COL4A1", "CTC1", "HTRA1", "TREX1", "W/O Mut"),  #90K
  #legend.labs = c("COL4A1", "GLA", "NOTCH3", "TREX1", "W/O Mut"),  #85K
  font.title    = c(14, "bold.italic", "darkgreen"),
  font.subtitle = c(14, "bold", "darkgreen"),
  font.caption  = c(14, "plain", "darkgreen"),
  font.x        = c(14, "bold.italic", "black"),
  font.y        = c(14, "bold.italic", "black"),
  font.xtickslab = c(14, "bold", "black"),
  font.ytickslab = c(14, "bold", "black")
) +
  labs(
    title    = "Cumulative probability for 85yrs",
    #subtitle = paste('Disease = ', "Circulatory system", " Top enriched 20 phecodes"),
    #subtitle = paste('Disease = ', "hematopoietic", " Top enriched 7 phecodes"),
    #subtitle = paste('Disease = ', "hypercoagulability", " 3 phecodes"),
    subtitle = paste('Disease = ', "Cerebral ischemia", " 6 phecodes"),
    #subtitle = paste('Disease = ', "Cerebral ischemia (X433.3)"),
    #subtitle = paste('Disease = ', "Cerebrovascular disease (X433)"),
    #subtitle = paste('Disease = ', "Cardiomyopathy (425)"),
    #subtitle = paste('Disease = ', "Primary/intrinsic cardiomyopathies (425.1)"),
    caption  = "Plotted with R survminer")
q$plot <- q$plot + ggplot2::annotate("text", x = 15, y = 0.10,# x and y coordinates of the text
                                     label = ANN, size = 5)

q$plot <- q$plot + 
  scale_x_continuous(breaks = sort(c(seq(0, 85, 5))))

#tiff(paste('D:/Geisinger/X19981_backup_10072020/Desktop/MS_project/survivalcurve_90K/', 'Cerebralvascular_90K_', 'X433.3_7genes', '_enriched_85yr', '.tiff', sep=''), units="in", width=15, height=10, res=600)
#tiff(paste('D:/Geisinger/X19981_backup_10072020/Desktop/MS_project/survivalcurve_90K/', 'Cerebralvascular_90K_', 'Circulatory_system_7genes_sep_90K', '_enriched_85yr', '.tiff', sep=''), units="in", width=15, height=10, res=600)
#tiff(paste('D:/Geisinger/X19981_backup_10072020/Desktop/MS_project/survivalcurve_90K/', 'Hypercoagulability_85K_black', '_system_7genes_sep_90K', '_enriched_85yr', '.tiff', sep=''), units="in", width=15, height=10, res=600)
#tiff(paste('D:/Geisinger/X19981_backup_10072020/Desktop/MS_project/survivalcurve_90K/', 'Cerebralischemia7_85K_black_', '_system_7genes_sep_90K', '_enriched_85yr', '.tiff', sep=''), units="in", width=15, height=10, res=600)
#tiff(paste('D:/Geisinger/X19981_backup_10072020/Desktop/MS_project/survivalcurve_90K/', 'Cardiomyopathy_85K_black_', '_system_7genes_sep_90K', '_enriched_85yr', '.tiff', sep=''), units="in", width=15, height=10, res=600)
tiff(paste('D:/Geisinger/X19981_backup_10072020/Desktop/MS_project/survivalcurve_90K/', 'Cerebralischemia7_90K_85K_black_', '_system_7genes_sep', '_enriched_85yr', '.tiff', sep=''), units="in", width=15, height=10, res=600)

print(q)
dev.off()

#select only one gene
GNSIS_select_onegene <- GNSIS_select[GNSIS_select$subgroup %in% c("NOTCH3", "W/O Mut"), ]
q <- ggsurvplot(
  surv_fit(Surv(timediff_recur_censor_85yr/365.25, event_recur_censor_85yr)~as.factor(subgroup), data = GNSIS_select_onegene), # survfit object with calculated statistics.
  fun = function(x) {1-x},  #plot cumulative probability F(t) = 1 - S(t)
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  pval.coord = c(2, 0.25),
  pval.size = 8,
  conf.int = FALSE,       # show confidence intervals for 
  censor.shape="|", censor.size = 4,
  fontsize = 4, #font size in the table
  font.legend = c(14, "plain", "black"),
  # point estimates of survival curves.
  xlim = c(0, 80),         # present narrower X axis, but not affect
  ylim = c(0, 0.30),
  # survival estimates.
  xlab = "Time in years",   # customize X axis label.
  ylab = "Cumulative Incidence",
  break.time.by = 1,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  n.risk = FALSE,
  #legend.labs = c("COL4A1", "CTC1", "GLA", "HTRA1", "NOTCH3", "NOTCH3", "W/O Mut"),    # change legend labels.
  font.title    = c(14, "bold.italic", "darkgreen"),
  font.subtitle = c(14, "bold", "darkgreen"),
  font.caption  = c(14, "plain", "darkgreen"),
  font.x        = c(14, "bold.italic", "black"),
  font.y        = c(14, "bold.italic", "black"),
  font.xtickslab = c(14, "bold", "black"),
  font.ytickslab = c(14, "bold", "black")
) +
  labs(
    title    = "Cumulative probability for 85yrs",
    subtitle = paste('Disease = ', "Cerebral ischemia (X433.3) for NOTCH3"),
    #subtitle = paste('Disease = ', "Cerebrovascular disease (X433)"),
    #subtitle = paste('Disease = ', "Cardiomyopathy (425)"),
    #subtitle = paste('Disease = ', "Primary/intrinsic cardiomyopathies (425.1)"),
    caption  = "Plotted with R survminer")
q$plot <- q$plot + ggplot2::annotate("text", x = 15, y = 0.25,# x and y coordinates of the text
                                     label = ANN, size = 5)

q$plot <- q$plot + 
  scale_x_continuous(breaks = sort(c(seq(0, 85, 5))))

tiff(paste('D:/Geisinger/X19981_backup_10072020/Desktop/MS_project/survivalcurve_90K/', 'Cerebralvascular_90K_', 'X433.3_NOTCH3', '_enriched_85yr', '.tiff', sep=''), units="in", width=15, height=10, res=600)
print(q)
dev.off()

# Pairwise survdiff
library(survival)
library(survminer)
res <- pairwise_survdiff(Surv(timediff_recur_censor_85yr/365.25, event_recur_censor_85yr)~subgroup, p.adjust.method = "none", data = GNSIS_select)
#res <- pairwise_survdiff(Surv(timediff_recur_censor_85yr/365.25, event_recur_censor_85yr)~group, p.adjust.method = "none", data = GNSIS_select)

res
# COL4A1  CTC1    GLA     HTRA1   NOTCH3  NOTCH3  
# CTC1    0.65735 -       -       -       -       -      
#   GLA     0.51524 0.97848 -       -       -       -      
#   HTRA1   0.64982 0.95690 0.96503 -       -       -      
#   NOTCH3  0.66385 0.69998 0.43875 0.42786 -       -      
#   NOTCH3   0.00369 0.29995 0.08654 0.44239 0.00195 -      
#   W/O Mut 0.66534 0.68696 0.41235 0.47185 0.87298 0.00019

# COL4A1  CTC1    GLA     HTRA1   NOTCH3  NOTCH3  
# CTC1    0.74135 -       -       -       -       -      
#   GLA     0.26751 0.45944 -       -       -       -      
#   HTRA1   0.77672 0.95085 0.38652 -       -       -      
#   NOTCH3  0.55944 0.81517 0.15476 0.93681 -       -      
#   NOTCH3   0.01586 0.12393 0.25158 0.20054 0.00351 -      
#   W/O Mut 0.28753 0.82213 0.08152 0.97926 0.88741 0.00035

# COL4A1 CTC1  GLA   HTRA1 NOTCH3 NOTCH3
# CTC1    0.283  -     -     -     -      -    
#   GLA     0.461  0.474 -     -     -      -    
#   HTRA1   0.711  0.434 0.419 -     -      -    
#   NOTCH3  0.297  0.257 0.206 0.898 -      -    
#   NOTCH3   0.609  0.560 0.851 0.608 0.312  -    
#   W/O Mut 0.476  0.353 0.603 0.598 0.091  0.787

# COL4A1 CTC1 GLA  HTRA1 NOTCH3 NOTCH3
# CTC1    0.41   -    -    -     -      -    
#   GLA     0.20   0.85 -    -     -      -    
#   HTRA1   0.94   0.56 0.51 -     -      -    
#   NOTCH3  0.49   0.36 0.12 0.77  -      -    
#   NOTCH3   0.66   0.35 0.17 0.98  0.96   -    
#   W/O Mut 0.49   0.49 0.27 0.85  0.21   0.46 

# COL4A1 CTC1  GLA   HTRA1 NOTCH3 NOTCH3
# CTC1    0.344  -     -     -     -      -    
#   GLA     0.531  0.223 -     -     -      -    
#   HTRA1   0.145  0.965 0.145 -     -      -    
#   NOTCH3  0.832  0.401 0.501 0.176 -      -    
#   NOTCH3   0.703  0.549 0.385 0.368 0.761  -    
#   W/O Mut 0.964  0.315 0.555 0.093 0.740  0.554

# COL4A1 CTC1  GLA   HTRA1 NOTCH3 NOTCH3
# CTC1    0.641  -     -     -     -      -    
#   GLA     0.562  0.470 -     -     -      -    
#   HTRA1   0.088  0.659 0.092 -     -      -    
#   NOTCH3  0.954  0.622 0.646 0.079 -      -    
#   NOTCH3   0.733  0.549 0.772 0.097 0.826  -    
#   W/O Mut 0.312  0.420 0.915 0.012 0.384  0.782

# COL4A1 CTC1  GLA   HTRA1 NOTCH3 NOTCH3
# CTC1    0.913  -     -     -     -      -    
#   GLA     0.208  0.423 -     -     -      -    
#   HTRA1   0.603  0.951 0.338 -     -      -    
#   NOTCH3  0.013  0.314 0.695 0.118 -      -    
#   NOTCH3   0.221  0.428 0.934 0.201 0.543  -    
#   W/O Mut 0.006  0.442 0.888 0.151 0.373  0.956

# COL4A1  CTC1    GLA     HTRA1   NOTCH3  NOTCH3  
# CTC1    0.53415 -       -       -       -       -      
#   GLA     0.13030 0.87020 -       -       -       -      
#   HTRA1   0.88029 0.66569 0.59669 -       -       -      
#   NOTCH3  0.04293 0.97376 0.80262 0.55621 -       -      
#   NOTCH3   0.42149 0.85428 0.57880 0.80827 0.64640 -      
#   W/O Mut 0.00027 0.99162 0.78275 0.52447 0.98620 0.59315

# COL4A1 CTC1  GLA   HTRA1 NOTCH3 NOTCH3
# CTC1    0.496  -     -     -     -      -    
#   GLA     0.055  0.872 -     -     -      -    
#   HTRA1   0.339  0.854 0.997 -     -      -    
#   NOTCH3  0.252  0.729 0.277 0.523 -      -    
#   NOTCH3   0.218  0.891 0.516 0.781 0.661  -    
#   W/O Mut 0.064  0.645 0.161 0.539 0.895  0.551

#No adjustment
# COL4A1 CTC1   GLA    HTRA1  NOTCH3 NOTCH3 
# CTC1    0.7602 -      -      -      -      -     
#   GLA     0.7472 0.5394 -      -      -      -     
#   HTRA1   0.5205 0.9411 0.2503 -      -      -     
#   NOTCH3  0.1498 0.3738 0.4939 0.1520 -      -     
#   NOTCH3   0.3650 0.5973 0.3408 0.9723 0.0624 -     
#   W/O Mut 0.0029 0.3085 0.2092 0.0618 0.4260 0.0090 

# COL4A1 CTC1  GLA   HTRA1 NOTCH3 NOTCH3
# CTC1    0.954  -     -     -     -      -    
#   GLA     0.684  0.772 -     -     -      -    
#   HTRA1   0.096  0.354 0.090 -     -      -    
#   NOTCH3  0.724  0.840 0.733 0.045 -      -    
#   NOTCH3   0.145  0.273 0.149 0.382 0.074  -    
#   W/O Mut 0.241  0.940 0.872 0.025 0.520  0.019

#BH
# COL4A1 CTC1  GLA   HTRA1 NOTCH3 NOTCH3
# CTC1    0.840  -     -     -     -      -    
#   GLA     0.840  0.708 -     -     -      -    
#   HTRA1   0.708  0.972 0.654 -     -      -    
#   NOTCH3  0.532  0.654 0.708 0.532 -      -    
#   NOTCH3   0.654  0.738 0.654 0.972 0.328  -    
#   W/O Mut 0.060  0.654 0.628 0.328 0.688  0.095
#P value adjustment method: BH 

# COL4A1 CTC1 GLA  HTRA1 NOTCH3 NOTCH3
# CTC1    0.95   -    -    -     -      -    
#   GLA     0.95   0.95 -    -     -      -    
#   HTRA1   0.33   0.67 0.33 -     -      -    
#   NOTCH3  0.95   0.95 0.95 0.32  -      -    
#   NOTCH3   0.39   0.57 0.39 0.67  0.33   -    
#   W/O Mut 0.56   0.95 0.95 0.27  0.84   0.27 

# COL4A1 CTC1   GLA    HTRA1  NOTCH3 NOTCH3 
# CTC1    0.9916 -      -      -      -      -     
#   GLA     0.9121 0.9916 -      -      -      -     
#   HTRA1   0.9916 0.9916 0.9916 -      -      -     
#   NOTCH3  0.4507 0.9916 0.9916 0.9916 -      -     
#   NOTCH3   0.9916 0.9916 0.9916 0.9916 0.9916 -     
#   W/O Mut 0.0056 0.9916 0.9916 0.9916 0.9916 0.9916

# COL4A1 CTC1 GLA  HTRA1 NOTCH3 NOTCH3
# CTC1    0.96   -    -    -     -      -    
#   GLA     0.66   0.71 -    -     -      -    
#   HTRA1   0.84   0.96 0.71 -     -      -    
#   NOTCH3  0.13   0.71 0.91 0.66  -      -    
#   NOTCH3   0.66   0.71 0.96 0.66  0.81   -    
#   W/O Mut 0.13   0.71 0.96 0.66  0.71   0.96 

#for v1.2
phecode_mapping <- read.csv("Phecode_map_v1_2_icd10cm_beta.csv", header = T, stringsAsFactors = F)
phecode_mapping_unique <- phecode_mapping %>%
  dplyr::select(phecode, phecode_str, exclude_range, exclude_name) %>%
  distinct(phecode, .keep_all = T)  #1755*4
str(phecode_mapping_unique)
phecode_mapping_unique$phecodeX <- as.character(paste0("X", phecode_mapping_unique$phecode))
####################################################################################################
colnames(GNSIS_select)
sum(is.na(GNSIS_select$timediff_occur_censor))
scale(GNSIS_select$timediff_occur_censor)
cox_25yrrecur_AA <- coxph(Surv(timediff_recur_censor_25yr/365.25, event_recur_censor_25yr)~as.factor(group) + as.factor(PT_SEX) + as.factor(PT_RACE), data = GNSIS_select)
#cox_25yrrecur_AA <- coxph(Surv(timediff_recur_censor_25yr/365.25, event_recur_censor_25yr)~PRS_norm_cat + as.factor(PT_SEX.x) + AGE_AT_INDEX, data = PRS_GNSIS)
coxtidy_25yrrecur_AA <- tidy(cox_25yrrecur_AA, exponentiate = T, conf.int = T, conf.level = 0.95)
coxtidy_25yrrecur_AA$Category <- "Aneurysm"

save(summary_final, anova_final, leven_final, plot_list, file = "dat_RECUR_all_matched_ratio2_phenotypes109_testresult_violinplot_85K.RData", version =2)

unique(Summary_phenotypes_merge_anno_common_enriched$exclude_name)
# [1] "circulatory system"      "genitourinary"           "endocrine/metabolic"     "dermatologic"           
# [5] "digestive"               "respiratory"             "neurological"            "hematopoietic"          
# [9] "injuries & poisonings"   "sense organs"            "mental disorders"        "musculoskeletal"        
# [13] "symptoms"                "pregnancy complications" "infectious diseases"     "neoplasms"              
# [17] "congenital anomalies" 

Summary_phenotypes_merge_anno_common_enriched_select <- Summary_phenotypes_merge_anno_common_enriched[Summary_phenotypes_merge_anno_common_enriched$exclude_name == "circulatory system", ]
phecode_select <- as.character(paste0("X", Summary_phenotypes_merge_anno_common_enriched_select$phecode))
print(phecode_select)

phecode_select
phenotypes_long_AA <- phenotypes_long %>% #change to relatedness or nonrelatedness files
  #dplyr::filter(phecode == "X433.3") %>% #Cerebral ischemia for X433; X433.3; #Cardiomyopathy for X425; X425.1
  dplyr::filter(phecode %in% phecode_select) %>% #select of enriched phecode set
  dplyr::filter(phenotype == 1)
AA_ID <- unique(paste0("PT", "", phenotypes_long_AA$id))


phenotypes_long_AA <- phenotypes_long %>% #change to relatedness or nonrelatedness files
  #dplyr::filter(phecode == "X433.3") %>% #Cerebral ischemia for X433; #Cardiomyopathy for X425
  dplyr::filter(phecode %in% phecode_select) %>% #select of enriched phecode set
  group_by(id) %>%
  dplyr::filter(sum(phenotype) == 0) %>%
  ungroup()
AA_WO <- unique(paste0("PT", "", phenotypes_long_AA$id))


head(phenotypes_long)
dat_RECUR_ICD109$disease <- ifelse(dat_RECUR_ICD109$PT_ID %in% AA_ID, 1, 0)
dat_RECUR_ICD109_select <- dat_RECUR_ICD109[dat_RECUR_ICD109$PT_ID %in% AA_ID, ] %>%
  #dplyr::filter(phecode == "X433.3") %>%
  dplyr::filter(phecode %in% phecode_select) %>% #select of enriched phecode set
  group_by(PT_ID) %>%
  arrange(desc(ENC_DT)) %>%
  slice(n()) %>%
  #mutate(disease = 1) %>%
  ungroup()
dat_RECUR_ICD109_select <- data.frame(dat_RECUR_ICD109_select)
print(table(dat_RECUR_ICD109_select$group))
#CS
# case control 
# 528     981 
#DI
# case control 
# 497     868 
# CA
# case control 
# 68     152
#GU
# case control 
# 494     969
#EM
# case control 
# 1165    2236

#433
# case control 
# 177     325 
# case control 
# 140     243

#433.3
# case control 
# 62     124 

#425
# case control 
# 51      88 

#425.1
# case control 
# 49      75

head(dat_RECUR_ICD109_select)

#wo phecode(s)
dat_RECUR_ICD109_WO <- dat_RECUR_ICD109[dat_RECUR_ICD109$PT_ID %in% AA_WO, ] %>%
  dplyr::group_by(PT_ID) %>%
  dplyr::arrange(ENC_DT) %>%
  slice(n()) %>%
  #mutate(disease = 0) %>%
  ungroup()
dat_RECUR_ICD109_WO <- data.frame(dat_RECUR_ICD109_WO)
GNSIS <- rbind(dat_RECUR_ICD109_select, dat_RECUR_ICD109_WO)

# load("dat_RECUR_casecontrol_select_combined_AA.RData")
# load("dat_RECUR_casecontrol_select_combined_bronchitis.RData")
# load("dat_RECUR_casecontrol_select_combined_latecvs.RData")
# load("dat_RECUR_casecontrol_select_combined_CKD.RData")
# load("dat_RECUR_casecontrol_select_combined_hypercoagulable.RData")
# load("dat_RECUR_casecontrol_select_combined_ulcer.RData")
# load("dat_RECUR_casecontrol_select_combined_HF.RData")
# load("dat_RECUR_casecontrol_select_combined_Stroke.RData")
# load("dat_RECUR_casecontrol_select_combined_DVT.RData")
# 
# colnames(dat_RECUR_casecontrol_select_combined)
# head(dat_RECUR_casecontrol_select_combined)
# [1] "PT_ID"                 "ENC_DT"                "phecode"               "phecode_str"          
# [5] "exclude_name"          "PT_BIRTH_DT"           "PT_DEATH_DT"           "PT_SEX"               
# [9] "PT_RACE"               "LAST_ACTIVE_DT"        "timediff_occur_censor" "group"                
# [13] "disease"    

# matched <- paste0("PT", phenotypes_matched$id) #selected matched data
# dat_RECUR_casecontrol_select_combined <- dat_RECUR_casecontrol_select_combined[dat_RECUR_casecontrol_select_combined$PT_ID %in% matched, ]

# GNSIS <- dat_RECUR_casecontrol_select_combined 
# head(GNSIS)
## Event: RECUR_STROKE (already available in data)
## Time to event: timediff_recur_censor
GNSIS <- GNSIS %>% 
  mutate(timediff_recur_censor = 
           case_when(disease == 1 ~ (timediff_occur_censor - 0)*365.25,
                     disease == 0 ~ (timediff_occur_censor - 0)*365.25)) 
#optional
GNSIS <- GNSIS %>% 
  mutate(timediff_recur_censor = 
           case_when(disease == 1 ~ (timediff_occur_censor - 50)*365.25,
                     disease == 0 ~ (timediff_occur_censor - 50)*365.25)) 

##   ____________________________________________________________________________
##   Limiting to 1/3/5 years for recurrence                                 ####

## Event: Limiting to RECUR STROKE in 5 year (event_recur_censor_5yr), 3 year (event_recur_censor_3yr), 1 year (event_recur_censor_1yr)
## Time to event: timediff_recur_censor_5yr, timediff_recur_censor_3yr, timediff_recur_censor_1yr
GNSIS <- GNSIS %>% 
  mutate(timediff_recur_censor_15yr =
           case_when(timediff_recur_censor < 15*365.25 ~ timediff_recur_censor,
                     timediff_recur_censor >= 15*365.25 ~ 15*365.25),
         event_recur_censor_15yr =
           case_when(disease == 0 ~ 0,
                     disease == 1 & timediff_recur_censor <= 15*365.25 ~1,
                     disease == 1 & timediff_recur_censor > 15*365.25 ~0),
         timediff_recur_censor_25yr =
           case_when(timediff_recur_censor < 25*365.25 ~ timediff_recur_censor,
                     timediff_recur_censor >= 25*365.25 ~ 25*365.25),
         event_recur_censor_25yr =
           case_when(disease == 0 ~ 0,
                     disease == 1 & timediff_recur_censor <= 25*365.25 ~1,
                     disease == 1 & timediff_recur_censor > 25*365.25 ~0),
         timediff_recur_censor_35yr =
           case_when(timediff_recur_censor < 35*365.25 ~ timediff_recur_censor,
                     timediff_recur_censor >= 35*365.25 ~ 35*365.25),
         event_recur_censor_35yr =
           case_when(disease == 0 ~ 0,
                     disease == 1 & timediff_recur_censor <= 35*365.25 ~1,
                     disease == 1 & timediff_recur_censor > 35*365.25 ~0),
         timediff_recur_censor_65yr =
           case_when(timediff_recur_censor < 65*365.25 ~ timediff_recur_censor,
                     timediff_recur_censor >= 65*365.25 ~ 65*365.25),
         event_recur_censor_65yr =
           case_when(disease == 0 ~ 0,
                     disease == 1 & timediff_recur_censor <= 65*365.25 ~1,
                     disease == 1 & timediff_recur_censor > 65*365.25 ~0)
         )
         

# GNSIS$timediff_recur_censor_15yr[GNSIS$timediff_recur_censor_15yr == 0] <- 0.5
# GNSIS$timediff_recur_censor_25yr[GNSIS$timediff_recur_censor_25yr == 0] <- 0.5
# GNSIS$timediff_recur_censor_35yr[GNSIS$timediff_recur_censor_35yr == 0] <- 0.5
colnames(GNSIS)
# [1] "PT_ID"                      "ENC_DT"                     "phecode"                   
# [4] "PT_BIRTH_DT"                "PT_DEATH_DT"                "PT_SEX"                    
# [7] "PT_RACE"                    "LAST_ACTIVE_DT"             "timediff_occur_censor"     
# [10] "group"                      "Index_age"                  "disease"                   
# [13] "timediff_recur_censor"      "timediff_recur_censor_15yr" "event_recur_censor_15yr"   
# [16] "timediff_recur_censor_25yr" "event_recur_censor_25yr"    "timediff_recur_censor_35yr"
# [19] "event_recur_censor_35yr" 

######################################################################################
library(survival)
library(survminer)
GNSIS_select <- GNSIS %>%
  filter(timediff_recur_censor_35yr >= 0)
#optional
GNSIS_select <- GNSIS
# KM stratified by cohort
tiff(paste('D:/Geisinger/X19981_backup_10072020/Desktop/MS_project/', 'Cerebralvascular_90K_EM_enriched_65yr', '.tiff', sep=''), units="in", width=15, height=15, res=600)

#tiff(paste('D:/Geisinger/X19981_backup_10072020/Desktop/MS_project/', 'Cerebralvascular_85K_X433.3_35yr_add', '.tiff', sep=''), units="in", width=15, height=15, res=600)
#tiff(paste('D:/Geisinger/X19981_backup_10072020/Desktop/MS_project/', 'Cardiomyopathy_85K_X425_add', '.tiff', sep=''), units="in", width=15, height=15, res=600)
ggsurvplot(
  surv_fit(Surv(timediff_recur_censor_65yr/365.25, event_recur_censor_65yr)~as.factor(group), data = GNSIS_select), # survfit object with calculated statistics.
  fun = function(x) {1-x},  #plot cumulative probability F(t) = 1 - S(t)
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  pval.coord = c(2, 0.5),
  pval.size = 8,
  conf.int = TRUE,       # show confidence intervals for 
  censor.shape="|", censor.size = 4,
  fontsize = 4, #font size in the table
  font.legend = c(14, "plain", "black"),
  # point estimates of survival curves.
  xlim = c(0, 65),         # present narrower X axis, but not affect
  ylim = c(0, 0.5),
  # survival estimates.
  xlab = "Time in years",   # customize X axis label.
  ylab = "Cumulative Incidence",
  break.time.by = 1,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  n.risk = FALSE,
  legend.labs =
    c("w/ rare alleles", "w/o rare alleles"),    # change legend labels.
  font.title    = c(14, "bold.italic", "darkgreen"),
  font.subtitle = c(14, "bold", "darkgreen"),
  font.caption  = c(14, "plain", "darkgreen"),
  font.x        = c(14, "bold.italic", "black"),
  font.y        = c(14, "bold.italic", "black"),
  font.xtickslab = c(14, "bold", "black"),
  font.ytickslab = c(14, "bold", "black")
) +
  labs(
    title    = "Cumulative probability for 65yrs",
    subtitle = paste('Disease = ', "20 enriched phecodes"),
    #subtitle = paste('Disease = ', "Cerebral ischemia (X433.3)"),
    #subtitle = paste('Disease = ', "Cardiomyopathy (425)"),
    caption  = "Plotted with R survminer")
dev.off()


# KM stratified by cohort
tiff(paste('D:/Geisinger/X19981_backup_10072020/Desktop/MS_project/', 'Stroke_35yrs', '.tiff', sep=''), units="in", width=15, height=15, res=600)
ggsurvplot(
  surv_fit(Surv(timediff_recur_censor_35yr/365.25, event_recur_censor_35yr)~as.factor(group), data = GNSIS_select), # survfit object with calculated statistics.
  fun = function(x) {1-x},  #plot cumulative probability F(t) = 1 - S(t)
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  #pval.coord = c(2, 0.4),
  pval.size = 8,
  conf.int = TRUE,       # show confidence intervals for 
  censor.shape="|", censor.size = 4,
  fontsize = 4, #font size in the table
  font.legend = c(14, "plain", "black"),
  # point estimates of survival curves.
  xlim = c(0, 35),         # present narrower X axis, but not affect
  ylim = c(0, 0.55),
  # survival estimates.
  xlab = "Time in years",   # customize X axis label.
  ylab = "Cumulative Incidence",
  break.time.by = 1,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  n.risk = FALSE,
  legend.labs =
    c("w/ rare alleles", "w/o rare alleles"),    # change legend labels.
  font.title    = c(14, "bold.italic", "darkgreen"),
  font.subtitle = c(14, "bold", "darkgreen"),
  font.caption  = c(14, "plain", "darkgreen"),
  font.x        = c(14, "bold.italic", "black"),
  font.y        = c(14, "bold.italic", "black"),
  font.xtickslab = c(14, "bold", "black"),
  font.ytickslab = c(14, "bold", "black")
) +
  labs(
    title    = "Cumulative probability for 35yrs",
    subtitle = paste('Disease = ', "Ischemic Stroke"),
    caption  = "Plotted with R survminer")
dev.off()

# KM stratified by cohort
tiff(paste('D:/Geisinger/X19981_backup_10072020/Desktop/MS_project/', 'Stroke_15yrs', '.tiff', sep=''), units="in", width=15, height=15, res=600)
ggsurvplot(
  surv_fit(Surv(timediff_recur_censor_15yr/365.25, event_recur_censor_15yr)~as.factor(group), data = GNSIS_select), # survfit object with calculated statistics.
  fun = function(x) {1-x},  #plot cumulative probability F(t) = 1 - S(t)
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  #pval.coord = c(2, 0.4),
  pval.size = 8,
  conf.int = TRUE,       # show confidence intervals for 
  censor.shape="|", censor.size = 4,
  fontsize = 4, #font size in the table
  font.legend = c(14, "plain", "black"),
  # point estimates of survival curves.
  xlim = c(0, 15),         # present narrower X axis, but not affect
  ylim = c(0, 0.40),
  # survival estimates.
  xlab = "Time in years",   # customize X axis label.
  ylab = "Cumulative Incidence",
  break.time.by = 1,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  n.risk = FALSE,
  legend.labs =
    c("w/ rare alleles", "w/o rare alleles"),    # change legend labels.
  font.title    = c(14, "bold.italic", "darkgreen"),
  font.subtitle = c(14, "bold", "darkgreen"),
  font.caption  = c(14, "plain", "darkgreen"),
  font.x        = c(14, "bold.italic", "black"),
  font.y        = c(14, "bold.italic", "black"),
  font.xtickslab = c(14, "bold", "black"),
  font.ytickslab = c(14, "bold", "black")
) +
  labs(
    title    = "Cumulative probability for 15yrs",
    subtitle = paste('Disease = ', "Ischemic Stroke"),
    caption  = "Plotted with R survminer")
dev.off()

######################################################################################

colnames(GNSIS_select)
sum(is.na(GNSIS_select$timediff_occur_censor))
scale(GNSIS_select$timediff_occur_censor)
cox_25yrrecur_AA <- coxph(Surv(timediff_recur_censor_25yr/365.25, event_recur_censor_25yr)~as.factor(group) + as.factor(PT_SEX) + as.factor(PT_RACE), data = GNSIS_select)
#cox_25yrrecur_AA <- coxph(Surv(timediff_recur_censor_25yr/365.25, event_recur_censor_25yr)~PRS_norm_cat + as.factor(PT_SEX.x) + AGE_AT_INDEX, data = PRS_GNSIS)
coxtidy_25yrrecur_AA <- tidy(cox_25yrrecur_AA, exponentiate = T, conf.int = T, conf.level = 0.95)
coxtidy_25yrrecur_AA$Category <- "Aneurysm"

cox_25yrrecur_COPD <- coxph(Surv(timediff_recur_censor_25yr/365.25, event_recur_censor_25yr)~as.factor(group) + as.factor(PT_SEX) + as.factor(PT_RACE), data = GNSIS_select)
#cox_25yrrecur_COPD <- coxph(Surv(timediff_recur_censor_25yr/365.25, event_recur_censor_25yr)~PRS_norm_cat + as.factor(PT_SEX.x) + AGE_AT_INDEX, data = PRS_GNSIS)
coxtidy_25yrrecur_COPD <- tidy(cox_25yrrecur_COPD, exponentiate = T, conf.int = T, conf.level = 0.95)
coxtidy_25yrrecur_COPD$Category <- "COPD"

cox_25yrrecur_lateCVS <- coxph(Surv(timediff_recur_censor_25yr/365.25, event_recur_censor_25yr)~as.factor(group) + as.factor(PT_SEX) + as.factor(PT_RACE), data = GNSIS_select)
#cox_25yrrecur_lateCVS <- coxph(Surv(timediff_recur_censor_25yr/365.25, event_recur_censor_25yr)~PRS_norm_cat + as.factor(PT_SEX.x) + AGE_AT_INDEX, data = PRS_GNSIS)
coxtidy_25yrrecur_lateCVS <- tidy(cox_25yrrecur_lateCVS, exponentiate = T, conf.int = T, conf.level = 0.95)
coxtidy_25yrrecur_lateCVS$Category <- "Sequelae of Stroke"

cox_25yrrecur_CKD <- coxph(Surv(timediff_recur_censor_25yr/365.25, event_recur_censor_25yr)~as.factor(group) + as.factor(PT_SEX) + as.factor(PT_RACE), data = GNSIS_select)
#cox_25yrrecur_CKD <- coxph(Surv(timediff_recur_censor_25yr/365.25, event_recur_censor_25yr)~PRS_norm_cat + as.factor(PT_SEX.x) + AGE_AT_INDEX, data = PRS_GNSIS)
coxtidy_25yrrecur_CKD <- tidy(cox_25yrrecur_CKD, exponentiate = T, conf.int = T, conf.level = 0.95)
coxtidy_25yrrecur_CKD$Category <- "Chronic Kidney Disease"

cox_25yrrecur_Hypercoagulable <- coxph(Surv(timediff_recur_censor_25yr/365.25, event_recur_censor_25yr)~as.factor(group) + as.factor(PT_SEX) + as.factor(PT_RACE), data = GNSIS_select)
#cox_25yrrecur_Hypercoagulable <- coxph(Surv(timediff_recur_censor_25yr/365.25, event_recur_censor_25yr)~PRS_norm_cat + as.factor(PT_SEX.x) + AGE_AT_INDEX, data = PRS_GNSIS)
coxtidy_25yrrecur_Hypercoagulable <- tidy(cox_25yrrecur_Hypercoagulable, exponentiate = T, conf.int = T, conf.level = 0.95)
coxtidy_25yrrecur_Hypercoagulable$Category <- "Hypercoagulable"

cox_25yrrecur_Ulcer <- coxph(Surv(timediff_recur_censor_25yr/365.25, event_recur_censor_25yr)~as.factor(group) + as.factor(PT_SEX) + as.factor(PT_RACE), data = GNSIS_select)
#cox_25yrrecur_Ulcer <- coxph(Surv(timediff_recur_censor_25yr/365.25, event_recur_censor_25yr)~PRS_norm_cat + as.factor(PT_SEX.x) + AGE_AT_INDEX, data = PRS_GNSIS)
coxtidy_25yrrecur_Ulcer <- tidy(cox_25yrrecur_Ulcer, exponentiate = T, conf.int = T, conf.level = 0.95)
coxtidy_25yrrecur_Ulcer$Category <- "Chronic Skin Ulcer"

cox_25yrrecur_HF <- coxph(Surv(timediff_recur_censor_25yr/365.25, event_recur_censor_25yr)~as.factor(group) + as.factor(PT_SEX) + as.factor(PT_RACE), data = GNSIS_select)
#cox_25yrrecur_HF <- coxph(Surv(timediff_recur_censor_25yr/365.25, event_recur_censor_25yr)~PRS_norm_cat + as.factor(PT_SEX.x) + AGE_AT_INDEX, data = PRS_GNSIS)
coxtidy_25yrrecur_HF <- tidy(cox_25yrrecur_HF, exponentiate = T, conf.int = T, conf.level = 0.95)
coxtidy_25yrrecur_HF$Category <- "Heart Failure"

cox_25yrrecur_Stroke <- coxph(Surv(timediff_recur_censor_25yr/365.25, event_recur_censor_25yr)~as.factor(group) + as.factor(PT_SEX) + as.factor(PT_RACE), data = GNSIS_select)
#cox_25yrrecur_Stroke <- coxph(Surv(timediff_recur_censor_25yr/365.25, event_recur_censor_25yr)~PRS_norm_cat + as.factor(PT_SEX.x) + AGE_AT_INDEX, data = PRS_GNSIS)
coxtidy_25yrrecur_Stroke <- tidy(cox_25yrrecur_Stroke, exponentiate = T, conf.int = T, conf.level = 0.95)
coxtidy_25yrrecur_Stroke$Category <- "Ischemic Stroke"
#####################################################################################

fooCox <- list(data.frame(coxtidy_25yrrecur_AA),
               data.frame(coxtidy_25yrrecur_COPD),
               data.frame(coxtidy_25yrrecur_lateCVS),
               data.frame(coxtidy_25yrrecur_CKD), 
               data.frame(coxtidy_25yrrecur_Hypercoagulable), 
               data.frame(coxtidy_25yrrecur_Ulcer), 
               data.frame(coxtidy_25yrrecur_HF), 
               data.frame(coxtidy_25yrrecur_Stroke))

result_final_Cox <- do.call(rbind, fooCox)
result_final_Cox$estimate_95CI <- paste0(round(1/result_final_Cox$estimate, 3), "(", round(1/result_final_Cox$conf.high, 3), "-", round(1/result_final_Cox$conf.low, 3), ")")
#result_final_Cox$term[result_final_Cox$term == 'as.factor(group)control'] <- "genetics"
result_final_Cox_select <- result_final_Cox %>%
  filter(term == 'as.factor(group)control')
result_final_Cox_select$p.value <- round(result_final_Cox_select$p.value, 3)
save.image("10042021_MS_project.RData", version = 2)

colnames(result_final_Cox_select)
#create a forest plot
tabletext <- result_final_Cox_select %>%
  dplyr::select(Category, p.value, estimate_95CI)

head(tabletext)

header <- c("Phecode Name", "P value", "Hazard Ratio (95%CI)")
tabletext <- rbind(header, tabletext)
forestdata <- result_final_Cox_select %>%
  dplyr::select(estimate, conf.high, conf.low,)
forestdata$estimate <- 1/forestdata$estimate
forestdata$conf.low <- 1/forestdata$conf.low
forestdata$conf.high <- 1/forestdata$conf.high
head(forestdata)
header2 <- c(1, 1, 1)
forestdata <- rbind(header2, forestdata)
forestdata_new <- forestdata[c(1:3, 4:9),] 
library(forestplot)
tabletext_new <- tabletext[c(1:3, 4:9),] 

head(tabletext_new)

tiff("Rplot_forestplot_multivariate_coxregression.tiff", units = "in", width = 10, height = 8, res = 300)
p <- forestplot(tabletext_new, 
                graph.pos = 3,
                forestdata_new,new_page = TRUE,
                #is.summary=c(TRUE,TRUE,rep(FALSE,8),TRUE),
                clip=c(-2,6),
                boxsize=0.25,
                zero = 1,
                xlog=FALSE, 
                col=fpColors(box="royalblue",line="darkblue", summary="royalblue"),
                xlab = "estimate(95%CI)",
                txt_gp = fpTxtGp(label = list(gpar(fontfamily = "Arial"),
                                              gpar(fontfamily = "",
                                                   col = "#660000")),
                                 ticks = gpar(fontfamily = "", cex = 1.0),
                                 xlab  = gpar(fontfamily = "Arial", cex = 1.0)))

grid::grid.text("Disease phenotype ~ rare variants of monogenic stroke risk genes", .5, 0.99, gp=gpar(cex=1.5))
dev.off()


#########################################################################################################
#calculate incidence rate difference and its confidence intervals
#Rothman KJ (2012) Epidemiology: An Introduction. 2nd Ed., Oxford University Press, Oxford.
library(fmsb)
res <- ratedifference(22, 25,	87946, 178884, CRC = TRUE, conf.level = 0.95)
res_result <- data.frame(tidy(res))

IR <- IRCI(22, 87946, conf.level = 0.95)
IR_result <- data.frame(t(c(IR$IR, IR$IRL, IR$IRU)))


dataFrame <- read.csv("table_summary_discovery.csv", header = T, stringsAsFactors = F)
head(dataFrame)
str(dataFrame)
#need to change from integer to numeric
dataFrame$carrier_case <- as.numeric(dataFrame$carrier_case)
dataFrame$noncarrier_case <- as.numeric(dataFrame$noncarrier_case)
dataFrame$carrier_personyear <- as.numeric(dataFrame$carrier_personyear)
dataFrame$noncarrier_personyear <- as.numeric(dataFrame$noncarrier_personyear)


# PheCodes carrier_case carrier_personyear noncarrier_case noncarrier_personyear
# 1 Other and unspecified coagulation defects(286.7)           22              87946              25                178884
# 2                     Hypercoagulable state(286.8)           22              87949              39                178755
# 3            Primary hypercoagulable state(286.81)           22              87949              38                178764
# 4                       Heart valve disorders(395)           89              87500             198                177623
# 5                       Abnormal heart sounds(396)           32              87888              93                178438
# 6                                Hypertension(401)          720              81756            1465                165418
# carrier_total noncarrier_total    onset
# 1          1695             3410 lifetime
# 2          1695             3410 lifetime
# 3          1695             3410 lifetime
# 4          1695             3410 lifetime
# 5          1695             3410 lifetime
# 6          1695             3410 lifetime

#rate difference
res <- ratedifference(dataFrame$carrier_case[1], dataFrame$noncarrier_case[1],	dataFrame$carrier_personyear[1], dataFrame$noncarrier_personyear[1], CRC = TRUE, conf.level = 0.95)

result <- NULL
for (i in 1:nrow(dataFrame)) {
  print(i)
  row <- dataFrame[i, ]
  #print(row)
  
  #result of IR
  IR <- IRCI(row$carrier_case, row$carrier_personyear, conf.level = 0.95)
  IR_result_case <- data.frame(t(c(IR$IR, IR$IRL, IR$IRU)))
  IR <- IRCI(row$noncarrier_case, row$noncarrier_personyear, conf.level = 0.95)
  IR_result_control <- data.frame(t(c(IR$IR, IR$IRL, IR$IRU)))
  IR_result <- cbind(IR_result_case, IR_result_control)
  IR_result$Phecode <- dataFrame$PheCodes[i]
  #print(IR_result)
  
  #calculate rate difference with 95%CI
  res <- ratedifference(row$carrier_case, row$noncarrier_case, row$carrier_personyear, row$noncarrier_personyear, CRC = TRUE, conf.level = 0.95)
  print(res)
  
  #result of RD
  res_result <- data.frame(tidy(res))
  res_result <- cbind(res_result, IR_result)
  result <- rbind(result, res_result)
  
}
write.csv(result, "result_for_ratedifference_discovery.csv")
saveRDS(result, file= "result_for_ratedifference_discovery.rds")
ratediff <- data.frame(readRDS("result_for_ratedifference_discovery.rds"))
dim(ratediff)

#risk difference
res <- riskdifference(dataFrame$carrier_case[1], dataFrame$noncarrier_case[1],	dataFrame$carrier_total[1], dataFrame$noncarrier_total[1], CRC = TRUE, conf.level = 0.95)

result <- NULL
for (i in 1:nrow(dataFrame)) {
  print(i)
  row <- dataFrame[i, ]
  #print(row)
  
  #result of IR
  R <- RCI(row$carrier_case, row$carrier_total, conf.level = 0.95)
  R_result_case <- data.frame(t(c(R$R, R$RL, R$RU)))
  R <- RCI(row$noncarrier_case, row$noncarrier_total, conf.level = 0.95)
  R_result_control <- data.frame(t(c(R$R, R$RL, R$RU)))
  R_result <- cbind(R_result_case, R_result_control)
  R_result$Phecode <- dataFrame$PheCodes[i]
  #print(IR_result)
  
  #calculate risk difference with 95%CI
  res <- riskdifference(row$carrier_case, row$noncarrier_case, row$carrier_total, row$noncarrier_total, CRC = TRUE, conf.level = 0.95)
  print(res)
  
  #result of RD
  res_result <- data.frame(tidy(res))
  res_result <- cbind(res_result, R_result)
  result <- rbind(result, res_result)
  
}
write.csv(result, "result_for_riskdifference_discovery.csv")
saveRDS(result, file= "result_for_riskdifference_discovery.rds")
riskdiff <- data.frame(readRDS("result_for_riskdifference_discovery.rds"))
dim(riskdiff)

