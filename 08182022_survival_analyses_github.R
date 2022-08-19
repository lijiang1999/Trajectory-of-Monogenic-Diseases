###outline:
#1. the following files created by other code files are required to start;
#2. preparing a longitudinal file for a specific PheCode set;
#3. creating censored table for KM analyses;
#4. Kaplan-Meier analysis and posthoc logRank test;
#5. performing coxph modeling;
#6. creating a forest plot for the coxph results

####################################################################################################
##the following files created by other code files are required to start.

load("dat_RECUR_ICD109.RData") #this file was created by '08182022_mapping_to_phecodes'
dat_RECUR_ICD109$phecode <- paste0("X", dat_RECUR_ICD109$phecode)
str(dat_RECUR_ICD109)

head(Summary_phenotypes_merge_anno_common) #this file was created by '08182022_create_RR'

load("phenotypes109.RData") #this file was created by '08182022_mapping_to_phecodes'
#load("phenotypes109_age.RData")

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

################################################################################
##preparing a longitudinal file for a specific phecode set
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
# phecode_select <- c("X286.8", "X286.81")
# phecode_select <- c("X425", "X425.1")
# phecode_select <- c("X425.1")
# phecode_select <- c("X433.1", "X433.11", "X433.12", "X433.2", "X433.21", "X433.3", "X433.31", "X433.32", "X433.6", "X433.8")
# phecode_select <- c("X433.2", "X433.21", "X433.3", "X433.31", "X433.32", "X433.6", "X433.8")
# phecode_select <- c("X433.2", "X433.21", "X433.3", "X433.31", "X433.32", "X433.8")
# phecode_select <- c("X433.2")

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

dat_RECUR_ICD109$disease <- ifelse(dat_RECUR_ICD109$PT_ID %in% AA_ID, 1, 0)
dat_RECUR_ICD109_select <- dat_RECUR_ICD109[dat_RECUR_ICD109$PT_ID %in% AA_ID, ] %>%
  dplyr::filter(phecode %in% phecode_select) %>%
  dplyr::group_by(PT_ID) %>%
  arrange(desc(ENC_DT)) %>%
  slice(n()) %>%
  ungroup()
dat_RECUR_ICD109_select <- data.frame(dat_RECUR_ICD109_select)
print(table(dat_RECUR_ICD109_select$group))

dat_RECUR_ICD109_WO <- dat_RECUR_ICD109[dat_RECUR_ICD109$PT_ID %in% AA_WO, ] %>%
  dplyr::group_by(PT_ID) %>%
  dplyr::arrange(ENC_DT) %>%
  slice(n()) %>%
  ungroup()
dat_RECUR_ICD109_WO <- data.frame(dat_RECUR_ICD109_WO)
print(table(dat_RECUR_ICD109_WO$group))

GNSIS <- rbind(dat_RECUR_ICD109_select, dat_RECUR_ICD109_WO)
##################################################################################
##Create censored table

# NN <- 'age for the follow-up starting from'
# GNSIS <- GNSIS %>% 
#   mutate(timediff_recur_censor = (timediff_occur_censor - NN)*365.25)

NN <- 0 #suggesting lifetime follow-up
GNSIS <- GNSIS %>% 
  mutate(timediff_recur_censor = (timediff_occur_censor - NN)*365.25)

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

#specify the patients with gene of mutation
load("Individualgene_PTID_12032021.RData")

head(GNSIS_select)
#GNSIS_select <- data.frame(GNSIS_select)
GNSIS_select <- GNSIS


# load("GNSIS_select_black_gene_90K.RData")
# GNSIS_select_discovery <- GNSIS_select
# load("GNSIS_select_black_gene_85K.RData")
# GNSIS_select_replication <- GNSIS_select

#no racial stratification
load("GNSIS_select_ischemicstroke_90K.RData")  
GNSIS_select_discovery <- GNSIS_select
load("GNSIS_select_ischemicstroke_85K.RData") 
GNSIS_select_replication <- GNSIS_select
GNSIS_select <- rbind(GNSIS_select_discovery, GNSIS_select_replication)

Noncarrier <- as.character(GNSIS_select[GNSIS_select$group == "control", ]$PT_ID)
GNSIS_select$subgroup <- "W/ Mut"
GNSIS_select$subgroup <- ifelse(GNSIS_select$PT_ID %in% HTRA1, "HTRA1", GNSIS_select$subgroup)
GNSIS_select$subgroup <- ifelse(GNSIS_select$PT_ID %in% Noncarrier, "W/O Mut", GNSIS_select$subgroup)
unique(GNSIS_select$subgroup)
table(GNSIS_select$subgroup)

GNSIS_select <- GNSIS_select[GNSIS_select$subgroup != "W/ Mut", ]

GNSIS_select$subgroup <- ifelse(GNSIS_select$subgroup == "W/O Mut", GNSIS_select$subgroup, 
                                ifelse(GNSIS_select$PT_ID %in% HTRA1_remove, GNSIS_select$subgroup, "remove"))

GNSIS_select <- GNSIS_select[GNSIS_select$subgroup != "remove", ]

#with racial stratification
load("GNSIS_select_white_gene_90K_duplicated_stroke.RData")
GNSIS_select_discovery <- GNSIS_select
load("GNSIS_select_white_gene_85K_duplicated_stroke.RData")
GNSIS_select_replication <- GNSIS_select

GNSIS_select <- rbind(GNSIS_select_discovery, GNSIS_select_replication)
GNSIS_select_WHITE <- GNSIS_select
GNSIS_select_BLACK <- GNSIS_select
GNSIS_select_BLACKWHITE <- rbind(GNSIS_select_WHITE, GNSIS_select_BLACK)

GNSIS_select <- GNSIS_select_BLACKWHITE[GNSIS_select_BLACKWHITE$COHORT == "Discovery", ]
GNSIS_select <- GNSIS_select_BLACKWHITE[GNSIS_select_BLACKWHITE$COHORT == "Replication", ]

table(GNSIS_select$subgroup, GNSIS_select$COHORT)
WHITE_ID <- unique(GNSIS_select[GNSIS_select$subgroup != "W/O Mut", ]$PT_ID) 
#save(WHITE_ID, file = "WHITE_ID.RData", version = 2)

BLACK_ID <- unique(GNSIS_select[GNSIS_select$subgroup != "W/O Mut", ]$PT_ID) 
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
GNSIS_select <- GNSIS_select[!(GNSIS_select$PT_ID %in% HTRA1_remove), ] 
#include HTRA1 p.Gln151Lys
GNSIS_select <- GNSIS_select[GNSIS_select$PT_ID %in% HTRA1_remove, ] 

####################################################################################
##Kaplan-Meier analysis and posthoc logRank test
ANN <- "ischemic stroke (6 PheCodes)"
#ANN <- "hypercoagulability (3 PheCodes)"
#ANN <- "cardiomyopathy (2 PheCodes)"

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

tiff(paste('<path_to_directory>/', 'ischemicstroke_6phecodes', '_system_7genes_sep', '_enriched_85yr_discoveryreplication_HTRA1_test_white', '.tiff', sep=''), units="in", width=15, height=10, res=600)
print(q)
dev.off()

GNSIS_select$subgroup <- paste0(GNSIS_select$subgroup, GNSIS_select$PT_SEX)
#GNSIS_select$subgroup <- paste0(GNSIS_select$subgroup, GNSIS_select$COHORT)
#GNSIS_select$subgroup <- paste0(GNSIS_select$subgroup, GNSIS_select$COHORT, GNSIS_select$PT_SEX)
table(GNSIS_select$subgroup)
res <- pairwise_survdiff(Surv(timediff_recur_censor_85yr/365.25, event_recur_censor_85yr)~subgroup, p.adjust.method = "none", data = GNSIS_select)
#res <- pairwise_survdiff(Surv(timediff_recur_censor_80yr/365.25, event_recur_censor_80yr)~group, p.adjust.method = "none", data = GNSIS_select)
res

###################################################################################################################

#performing coxph modeling
cox_65yrrecur_COL4A1 <- coxph(Surv(timediff_recur_censor_65yr/365.25, event_recur_censor_65yr)~as.factor(group) + as.factor(PT_SEX) + index_age + PC1 + PC2 + PC3 + PC4 + PC5, data = GNSIS_select)

coxtidy_65yrrecur_COL4A1 <- tidy(cox_65yrrecur_COL4A1, exponentiate = T, conf.int = T, conf.level = 0.95)
coxtidy_65yrrecur_COL4A1$Onset <- "Early"
coxtidy_65yrrecur_COL4A1$Category <- "Cerebral ischemia"
coxtidy_65yrrecur_COL4A1$Sex <- "Combined"
coxtidy_65yrrecur_COL4A1 <- data.frame(coxtidy_65yrrecur_COL4A1)

cox_85yrrecur_COL4A1 <- coxph(Surv(timediff_recur_censor_85yr/365.25, event_recur_censor_85yr)~as.factor(group) + as.factor(PT_SEX) + index_age + PC1 + PC2 + PC3 + PC4 + PC5, data = GNSIS_select)

coxtidy_85yrrecur_COL4A1 <- tidy(cox_85yrrecur_COL4A1, exponentiate = T, conf.int = T, conf.level = 0.95)
coxtidy_85yrrecur_COL4A1$Onset <- "Lifetime"
coxtidy_85yrrecur_COL4A1$Category <- "Cerebral ischemia"
coxtidy_85yrrecur_COL4A1$Sex <- "Combined"
coxtidy_85yrrecur_COL4A1 <- data.frame(coxtidy_85yrrecur_COL4A1)

fooCox_COL4A1_90K_WHITE <- rbind(coxtidy_65yrrecur_COL4A1, coxtidy_85yrrecur_COL4A1)

GNSIS_select_COL4A1_FEMALE <- GNSIS_select[GNSIS_select$PT_SEX == "Female", ]

cox_65yrrecur_COL4A1 <- coxph(Surv(timediff_recur_censor_65yr/365.25, event_recur_censor_65yr)~as.factor(group) + index_age + PC1 + PC2 + PC3 + PC4 + PC5, data = GNSIS_select_COL4A1_FEMALE)

coxtidy_65yrrecur_COL4A1 <- tidy(cox_65yrrecur_COL4A1, exponentiate = T, conf.int = T, conf.level = 0.95)
coxtidy_65yrrecur_COL4A1$Onset <- "Early"
coxtidy_65yrrecur_COL4A1$Category <- "Cerebral ischemia"
coxtidy_65yrrecur_COL4A1$Sex <- "Female"
coxtidy_65yrrecur_COL4A1 <- data.frame(coxtidy_65yrrecur_COL4A1)

cox_85yrrecur_COL4A1 <- coxph(Surv(timediff_recur_censor_85yr/365.25, event_recur_censor_85yr)~as.factor(group) + index_age + PC1 + PC2 + PC3 + PC4 + PC5, data = GNSIS_select_COL4A1_FEMALE)
coxtidy_85yrrecur_COL4A1 <- tidy(cox_85yrrecur_COL4A1, exponentiate = T, conf.int = T, conf.level = 0.95)
coxtidy_85yrrecur_COL4A1$Onset <- "Lifetime"
coxtidy_85yrrecur_COL4A1$Category <- "Cerebral ischemia"
coxtidy_85yrrecur_COL4A1$Sex <- "Female"
coxtidy_85yrrecur_COL4A1 <- data.frame(coxtidy_85yrrecur_COL4A1)

fooCox_COL4A1_90K_WHITE_FEMALE <- rbind(coxtidy_65yrrecur_COL4A1, coxtidy_85yrrecur_COL4A1)

GNSIS_select_COL4A1_MALE <- GNSIS_select[GNSIS_select$PT_SEX == "Male", ]

cox_65yrrecur_COL4A1 <- coxph(Surv(timediff_recur_censor_65yr/365.25, event_recur_censor_65yr)~as.factor(group) + index_age + PC1 + PC2 + PC3 + PC4 + PC5, data = GNSIS_select_COL4A1_MALE)

coxtidy_65yrrecur_COL4A1 <- tidy(cox_65yrrecur_COL4A1, exponentiate = T, conf.int = T, conf.level = 0.95)
coxtidy_65yrrecur_COL4A1$Onset <- "Early"
coxtidy_65yrrecur_COL4A1$Category <- "Cerebral ischemia"
coxtidy_65yrrecur_COL4A1$Sex <- "Male"
coxtidy_65yrrecur_COL4A1 <- data.frame(coxtidy_65yrrecur_COL4A1)

cox_85yrrecur_COL4A1 <- coxph(Surv(timediff_recur_censor_85yr/365.25, event_recur_censor_85yr)~as.factor(group) + index_age + PC1 + PC2 + PC3 + PC4 + PC5, data = GNSIS_select_COL4A1_MALE)

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


##create a forest plot for the coxph results
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
#grid::grid.text("STROKE ~ having rare variants + index_age + (Sex) + PCs(1-5)", .5, 0.99, gp=gpar(cex=1.5))
dev.off()

