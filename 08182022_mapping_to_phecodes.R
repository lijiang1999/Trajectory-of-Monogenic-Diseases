###outline:
#1. processing EHR longitudinal data pullout from case:control(here referred to as carriers:noncarriers) at 1:5 ratio without matching
#2. creating file for all cases or controls mapped to ICD10 or ICD9 and timediff_occur_censor for a specific ICD under a specific encounter
#3. using LAST_ACTIVE_DT make more sense for index date and creating the final demographic data;
#4. making files for PheCode mapping (including making files for equal or less than 65yrs. This is optional for creating early onset files for HheCode mapping);
#5. PheCode mapping using PheWAS package;
#6. post mapping file processing

####################################################################################################################
##processing EHR pullout from 5x case:control
library(broom)
library(car)

#175K datapull (this is to pulldown demographic data)
datapull_09162021_cases <- read.csv("<path_to_directory>/sevengenes_cases.csv", header = F, stringsAsFactors = F)
head(datapull_09162021_cases) 
datapull_09162021_cases <- datapull_09162021_cases[, c(2:4, 6:7, 13)]
colnames(datapull_09162021_cases) <- c("PT_BIRTH_DT", "PT_DEATH_DT", "PT_STATUS", "PT_SEX", "PT_RACE_1", "PT_ID")

#175K datapull (this is to pulldown ICD9 data)
datapull_09162021_cases_ICD9 <- read.csv("<path_to_directory>/sevengenes_cases_ICD9.csv", header = F, stringsAsFactors = F)
head(datapull_09162021_cases_ICD9)
datapull_09162021_cases_ICD9 <- datapull_09162021_cases_ICD9[, c(2:4, 7, 8, 12, 14, 16, 25, 27)]
colnames(datapull_09162021_cases_ICD9) <- c("PT_ID", "ENC_TYPE", "DX_CD", "DX_PRIMARY_YN", "DX_ON_ORD_YN", "CODE_SYSTEM", "CS_CD", "CONCEPT_NM", "ENC_ID", "ENC_DT")


#175K datapull (this is to pulldown ICD10 data)
datapull_09162021_cases_ICD10 <- read.csv("<path_to_directory>/sevengenes_cases_ICD10.csv", header = F, stringsAsFactors = F)
head(datapull_09162021_cases_ICD10) 
datapull_09162021_cases_ICD10 <- datapull_09162021_cases_ICD10[, c(2:4, 7, 8, 12, 14, 16, 25, 27)]
colnames(datapull_09162021_cases_ICD10) <- c("PT_ID", "ENC_TYPE", "DX_CD", "DX_PRIMARY_YN", "DX_ON_ORD_YN", "CODE_SYSTEM", "CS_CD", "CONCEPT_NM", "ENC_ID", "ENC_DT")

save(datapull_09162021_cases, datapull_09162021_cases_ICD9, datapull_09162021_cases_ICD10, file = "datapull_09162021_cases_datapull.RData", version = 2)

#same way to process controls

#ICD10 code mapping
phecode_mapping <- read.csv("Phecode_map_v1_2_icd10cm_beta.csv", header = T, stringsAsFactors = F)
str(phecode_mapping)
#for v1.2b
phecode_mapping_unique <- phecode_mapping %>%
  dplyr::select(phecode, phecode_str, exclude_range, exclude_name) %>%
  distinct(phecode, .keep_all = T)  #1755*4
str(phecode_mapping_unique)

#ICD9 code mapping
phecode_mapping_ICD9 <- read.csv("phecode_icd9_map_unrolled.csv", header = T, stringsAsFactors = F)
str(phecode_mapping_ICD9)

##creating file for all cases or controls mapped to ICD10 or ICD9
SelectedPT <- data.frame()
unique_file <- data.frame()
RECUR <- NULL
dat_RECUR <- data.frame()

Result_longitudinal_cases = function(unique_file, datapull_file, PT) {   #datapull_file  with demo information;
  unique_file <- unique_file %>%
    dplyr::select(PT_ID, CS_CD, ENC_DT) %>%
    group_by(PT_ID, CS_CD) %>%
    arrange(desc(ENC_DT)) %>%
    slice(n()) %>%
    ungroup()
  unique_file$icd10cm <- unique_file$CS_CD
  #head(unique_file, 5) #sanity check
  
  #provide selected patient ID
  for (i in PT){
    print(i)
    SelectedPT <- unique_file[unique_file$PT_ID == i, ]
    SelectedPT <- merge(SelectedPT, phecode_mapping, by = "icd10cm")
    SelectedPT$PT_BIRTH_DT <- datapull_file[datapull_file$PT_ID == i, ]$PT_BIRTH_DT[1]
    SelectedPT$PT_DEATH_DT <- sort(datapull_file[datapull_file$PT_ID == i, ]$PT_DEATH_DT, decreasing = TRUE)[1]
    SelectedPT$PT_SEX <- datapull_file[datapull_file$PT_ID == i, ]$PT_SEX[1]
    SelectedPT$PT_RACE <- datapull_file[datapull_file$PT_ID == i, ]$PT_RACE_1[1]
    SelectedPT$LAST_ACTIVE_DT <- sort(as.POSIXlt(unique_file[unique_file$PT_ID == i, ]$ENC_DT), decreasing = TRUE)[1]
    SelectedPT <- SelectedPT %>% 
      mutate(timediff_occur_censor = 
               case_when(!is.na(ENC_DT) ~ as.numeric(difftime(ENC_DT, PT_BIRTH_DT, units = "days")),
                         is.na(ENC_DT) ~ as.numeric(difftime(LAST_ACTIVE_DT, PT_BIRTH_DT, units = "days")))) #change to PT_BIRTH_DT not ENC_DT (04082022)
    SelectedPT$timediff_occur_censor <- SelectedPT$timediff_occur_censor/365.65
    SelectedPT$phecode <- as.numeric(SelectedPT$phecode)
    SelectedPT <- SelectedPT %>%
      dplyr::select(PT_ID, ENC_DT, phecode, phecode_str, exclude_name, PT_BIRTH_DT, PT_DEATH_DT, PT_SEX, PT_RACE, LAST_ACTIVE_DT, timediff_occur_censor)
    
    RECUR[[i]] <- SelectedPT
    #print(head(SelectedPT, 5)) #sanity check
  }
  dat_RECUR <- data.frame(do.call(rbind,RECUR))
  #print(dat_RECUR) #sanity check
  dat_RECUR$group <- "case"
  save(dat_RECUR, file="cases_trajectory_all_ratio2_ICD10.RData", version = 2)
  return(dat_RECUR)
}

Result_longitudinal_cases_result <- Result_longitudinal_cases(datapull_09162021_cases_ICD10, datapull_09162021_cases, unique(datapull_09162021_cases_ICD10$PT_ID))
head(Result_longitudinal_cases_result)

#same way to create file for all controls mapping with ICD10
#same way to create file for all cases/controls mapping with ICD9

load("cases_trajectory_all_ratio2_ICD10.RData")
dat_RECUR_case_ICD10 <- dat_RECUR 
load("controls_trajectory_all_ratio2_ICD10.RData")
dat_RECUR_control_ICD10 <- dat_RECUR


dat_RECUR_ICD10 <- rbind(dat_RECUR_case_ICD10, dat_RECUR_control_ICD10)
head(dat_RECUR_ICD10)
dat_RECUR_ICD10 <- dat_RECUR_ICD10[, !colnames(dat_RECUR_ICD10) %in% c("phecode_str", "exclude_name")]


load("cases_trajectory_all_ratio2_ICD9.RData")
dat_RECUR_case_ICD9 <- dat_RECUR
load("controls_trajectory_all_ratio2_ICD9.RData")
dat_RECUR_control_ICD9 <- dat_RECUR

dat_RECUR_ICD9 <- rbind(dat_RECUR_case_ICD9, dat_RECUR_control_ICD9)
head(dat_RECUR_ICD9)
dat_RECUR_ICD9 <- dat_RECUR_ICD9[, !colnames(dat_RECUR_ICD9) %in% c("phecode_str", "exclude_name")]


dat_RECUR_ICD109 <- rbind(dat_RECUR_ICD10, dat_RECUR_ICD9)
save(dat_RECUR_ICD109, file = "dat_RECUR_ICD109.RData", version = 2)

##using LAST_ACTIVE_DT make more sense for index date
library(lubridate)
colnames(dat_RECUR_ICD109)
# [1] "PT_ID"                 "ENC_DT"                "phecode"
# [4] "PT_BIRTH_DT"           "PT_DEATH_DT"           "PT_SEX"
# [7] "PT_RACE"               "LAST_ACTIVE_DT"        "timediff_occur_censor"
# [10] "group"

dat_RECUR_ICD109 <- dat_RECUR_ICD109 %>% 
  mutate(Index_age = 
           case_when(is.na(LAST_ACTIVE_DT) ~ as.numeric(difftime(ENC_DT, PT_BIRTH_DT, units = "days")),
                     !is.na(LAST_ACTIVE_DT) ~ as.numeric(difftime(LAST_ACTIVE_DT, PT_BIRTH_DT, units = "days")))) 

#Create the index age
dat_RECUR_ICD109$Index_age <- dat_RECUR_ICD109$Index_age/365.65

#Create final demographics file
dat_RECUR_ICD109_demographics <- dat_RECUR_ICD109 %>%
  dplyr::select(PT_ID, PT_BIRTH_DT, PT_DEATH_DT, PT_SEX, PT_RACE, LAST_ACTIVE_DT, group, Index_age) %>%
  group_by(PT_ID) %>%
  arrange(Index_age) %>%
  slice(n()) %>%
  ungroup()
dat_RECUR_ICD109_demographics <- data.frame(dat_RECUR_ICD109_demographics)

######################################################################################################
##make files for PheCode mapping
datapull_09162021_casescontrols <- rbind(datapull_09162021_cases, datapull_09162021_controls)
datapull_09162021_casescontrols <- datapull_09162021_casescontrols %>%
  distinct(PT_ID, .keep_all = T) 
datapull_09162021_casescontrols$id = as.integer(gsub("PT", "", datapull_09162021_casescontrols$PT_ID))
head(datapull_09162021_casescontrols)

#identify the encounter date for the diagnosis of the first time for either cases or controls
#ICD10
#for controls
datapull_09162021_controls_ICD10_unique <- datapull_09162021_controls_ICD10 %>%
  dplyr::select(PT_ID, CS_CD, ENC_DT) %>%
  group_by(PT_ID, CS_CD) %>%
  arrange(desc(ENC_DT)) %>%
  slice(n()) %>%
  ungroup() 

datapull_09162021_controls_ICD10_count <- datapull_09162021_controls_ICD10 %>%
  dplyr::select(PT_ID, CS_CD, ENC_DT) %>%
  group_by(PT_ID, CS_CD) %>%
  summarise(count=n()) %>%
  ungroup() 

datapull_09162021_controls_ICD10_count_unique <- inner_join(datapull_09162021_controls_ICD10_count, datapull_09162021_controls_ICD10_unique, by = c( "PT_ID", "CS_CD"))

length(unique(datapull_09162021_controls_ICD10_count_unique$PT_ID))

#ICD9
datapull_09162021_controls_ICD9_unique <- datapull_09162021_controls_ICD9 %>%
  dplyr::select(PT_ID, CS_CD, ENC_DT) %>%
  group_by(PT_ID, CS_CD) %>%
  arrange(desc(ENC_DT)) %>%
  slice(n()) %>%
  ungroup() 

datapull_09162021_controls_ICD9_count <- datapull_09162021_controls_ICD9 %>%
  dplyr::select(PT_ID, CS_CD, ENC_DT) %>%
  group_by(PT_ID, CS_CD) %>%
  summarise(count=n()) %>%
  ungroup() 

datapull_09162021_controls_ICD9_count_unique <- inner_join(datapull_09162021_controls_ICD9_count, datapull_09162021_controls_ICD9_unique, by = c( "PT_ID", "CS_CD"))
head(datapull_09162021_controls_ICD10_count_unique)

#for cases
datapull_09162021_cases_ICD10_unique <- datapull_09162021_cases_ICD10 %>%
  dplyr::select(PT_ID, CS_CD, ENC_DT) %>%
  group_by(PT_ID, CS_CD) %>%
  arrange(desc(ENC_DT)) %>%
  slice(n()) %>%
  ungroup() 

datapull_09162021_cases_ICD10_count <- datapull_09162021_cases_ICD10 %>%
  dplyr::select(PT_ID, CS_CD, ENC_DT) %>%
  group_by(PT_ID, CS_CD) %>%
  summarise(count=n()) %>%
  ungroup() 

datapull_09162021_cases_ICD10_count_unique <- inner_join(datapull_09162021_cases_ICD10_count, datapull_09162021_cases_ICD10_unique, by = c( "PT_ID", "CS_CD"))

length(unique(datapull_09162021_cases_ICD10_count_unique$PT_ID)) 

#ICD9 more patients mapping to ICD9 codes
datapull_09162021_cases_ICD9_unique <- datapull_09162021_cases_ICD9 %>%
  dplyr::select(PT_ID, CS_CD, ENC_DT) %>%
  group_by(PT_ID, CS_CD) %>%
  arrange(desc(ENC_DT)) %>%
  slice(n()) %>%
  ungroup()

datapull_09162021_cases_ICD9_count <- datapull_09162021_cases_ICD9 %>%
  dplyr::select(PT_ID, CS_CD, ENC_DT) %>%
  group_by(PT_ID, CS_CD) %>%
  summarise(count=n()) %>%
  ungroup()

datapull_09162021_cases_ICD9_count_unique <- inner_join(datapull_09162021_cases_ICD9_count, datapull_09162021_cases_ICD9_unique, by = c( "PT_ID", "CS_CD"))

Summary_ICD10_control <- datapull_09162021_controls_ICD10_unique %>%
  group_by(CS_CD) %>%
  summarise(count=n()) %>%
  ungroup() %>%
  arrange(desc(count))

Summary_ICD10_case <- datapull_09162021_cases_ICD10_unique %>%
  group_by(CS_CD) %>%
  summarise(count=n()) %>%
  ungroup() %>%
  arrange(desc(count))

#merged cases and controls
datapull_09162021_casescontrols_ICD10 <- rbind(datapull_09162021_controls_ICD10, datapull_09162021_cases_ICD10)
head(datapull_09162021_casescontrols_ICD10)

datapull_09162021_casescontrols_ICD10_count <- rbind(datapull_09162021_controls_ICD10_count_unique,datapull_09162021_cases_ICD10_count_unique)
head(datapull_09162021_casescontrols_ICD10_count)

################################################################################################################
# install.packages("devtools")
# install.packages(c("dplyr","tidyr","ggplot2","MASS","meta","ggrepel","DT"))
# devtools::install_github("PheWAS/PheWAS")

library(PheWAS)
?phewas

#ICD10
head(datapull_09162021_casescontrols_ICD10_count)

colnames(datapull_09162021_casescontrols_ICD10_count) <- c("id", "code", "count", "ENC_DT")
datapull_09162021_casescontrols_ICD10_count$vocabulary_id <- "ICD10CM"
str(datapull_09162021_casescontrols_ICD10_count)

datapull_09162021_casescontrols_ICD10_count$id <- gsub("PT", "", datapull_09162021_casescontrols_ICD10_count$id)
datapull_09162021_casescontrols_ICD10_count$id <- as.integer(datapull_09162021_casescontrols_ICD10_count$id)

head(datapull_09162021_casescontrols_ICD10_count)
datapull_09162021_casescontrols_ICD10_count1 <- datapull_09162021_casescontrols_ICD10_count %>%
  dplyr::select(id, vocabulary_id, code, count)

#ICD9

colnames(datapull_09162021_casescontrols_ICD9_count) <- c("id", "code", "count", "ENC_DT")
datapull_09162021_casescontrols_ICD9_count$vocabulary_id <- "ICD9CM"
str(datapull_09162021_casescontrols_ICD9_count)

datapull_09162021_casescontrols_ICD9_count$id <- gsub("PT", "", datapull_09162021_casescontrols_ICD9_count$id)
datapull_09162021_casescontrols_ICD9_count$id <- as.integer(datapull_09162021_casescontrols_ICD9_count$id)

head(datapull_09162021_casescontrols_ICD9_count)
datapull_09162021_casescontrols_ICD9_count1 <- datapull_09162021_casescontrols_ICD9_count %>%
  dplyr::select(id, vocabulary_id, code, count)

datapull_09162021_casescontrols_ICD9_count_age <- merge(datapull_09162021_casescontrols_ICD9_count, 
                                                        datapull_09162021_casescontrols, by = "id")

datapull_09162021_casescontrols_ICD9_count_age <- datapull_09162021_casescontrols_ICD9_count_age %>% 
  mutate(timediff_occur_censor = as.numeric(difftime(ENC_DT, PT_BIRTH_DT, units = "days"))) %>%
  mutate(timediff_occur_censor = timediff_occur_censor/365.65)

#equal or less than 65yrs(this is optional for creating early onset files for phecode mapping)
datapull_09162021_casescontrols_ICD9_count_age$count_age <- ifelse(datapull_09162021_casescontrols_ICD9_count_age$timediff_occur_censor <= 65, datapull_09162021_casescontrols_ICD9_count_age$count, 0)
datapull_09162021_casescontrols_ICD9_count1_age <- datapull_09162021_casescontrols_ICD9_count_age %>%
  dplyr::select(id, vocabulary_id, code, count_age) %>%
  rename(count = count_age)

str(datapull_09162021_casescontrols_ICD9_count1_age)
#age filtered at 65
save(datapull_09162021_casescontrols_ICD9_count_age, datapull_09162021_casescontrols_ICD9_count1_age, file = "datapull_05312022_casescontrols_ICD9_count1_age.RData", version = 2)
save(datapull_09162021_casescontrols_ICD9_count_age, datapull_09162021_casescontrols_ICD9_count1_age, file = "datapull_06062022_casescontrols_ICD9_ratio2_discovery_count1_age.RData", version = 2)


#adding sex to control mislabeling.
head(datapull_09162021_casescontrols)

id.sex.phenotype109 <- data.frame(datapull_09162021_casescontrols$id, datapull_09162021_casescontrols$PT_SEX)
colnames(id.sex.phenotype109) <- c("id", "sex")
dim(id.sex.phenotype109) #29708*2
unique(id.sex.phenotype109$sex)

id.sex.phenotype109$sex[id.sex.phenotype109$sex == "Female"] <- "F"
id.sex.phenotype109$sex[id.sex.phenotype109$sex == "Male"] <- "M"

#ICD9 and ICD10
datapull_09162021_casescontrols_ICD109_count1 <- rbind(datapull_09162021_casescontrols_ICD10_count1, datapull_09162021_casescontrols_ICD9_count1)

phenotypes=createPhenotypes(datapull_09162021_casescontrols_ICD109_count1, 
                            aggregate.fun=sum, id.sex=id.sex.phenotype109)

#also create early onset phecode mapping
datapull_09162021_casescontrols_ICD109_count1_age <- rbind(datapull_09162021_casescontrols_ICD10_count1_age, datapull_09162021_casescontrols_ICD9_count1_age)

phenotypes_age=createPhenotypes(datapull_09162021_casescontrols_ICD109_count1_age, 
                                aggregate.fun=sum, id.sex=id.sex.phenotype109)


##post-mapping file processing
dim(phenotypes) 
unique(phenotypes$id)

cols <- sapply(phenotypes, is.logical)
phenotypes[,cols] <- lapply(phenotypes[,cols], as.numeric)
phenotypes[is.na(phenotypes)] <- 0

#group column representing with RVs(case) or without RVs(control) 
cases_ID_final <- c(ALL_90K, ALL_85K)
cases_ID_final <- as.integer(gsub("PT", "", cases_ID_final))
phenotypes$group <- ifelse(phenotypes$id %in% cases_ID_final, "case", "control") #this case/control represents with/without RVs
table(phenotypes$group)

save(phenotypes, file = "phenotypes109_sexverified.RData", version = 2) #this is phenotypes file mapped by both ICD10 and ICD9


cols <- sapply(phenotypes_age, is.logical)
phenotypes_age[,cols] <- lapply(phenotypes_age[,cols], as.numeric)
phenotypes_age[is.na(phenotypes_age)] <- 0

#group column representing with RVs(case) or without RVs(control) 
cases_ID_final <- c(ALL_90K, ALL_85K)
cases_ID_final <- as.integer(gsub("PT", "", cases_ID_final))
phenotypes_age$group <- ifelse(phenotypes_age$id %in% cases_ID_final, "case", "control")
table(phenotypes_age$group)

save(phenotypes_age, file = "phenotypes109_age_sexverified.RData", version = 2) #this is phenotypes file mapped by both ICD10 and ICD9
