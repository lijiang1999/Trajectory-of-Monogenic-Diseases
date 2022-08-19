###outline:
#1. calculate risk difference and percent relative effect (RR-1);
#2. identify common and rare disease phenotypes in the cohort;
#3. identify the top ranked PheCodes based on percent relative effect;
#4. Fisher exact test to determine the association between RVs and phenotypes




#Phenotypes_ratio2 file was created by PheCode mapping procedure
colnames(phenotypes_ratio2)
gathercols <- colnames(phenotypes_ratio2[, 2:1858]) 
keycol <- "phecode"
valuecol <- "phenotype"
phenotypes_long <- data.frame(gather_(phenotypes_ratio2, keycol, valuecol, gathercols, factor_key=TRUE))

head(phenotypes_long)
# id group phecode phenotype
# 1 1000149     0    X008         0
# 2 1000246     1    X008         0
# 3 1001322     0    X008         0
# 4 1001423     0    X008         0
# 5 1001689     0    X008         0
# 6  100170     0    X008         0

phecode <- unique(phenotypes_long$phecode)
phecode

#remove some cases(optional) for sensitivity analsysis
# Remove_cases <- gsub("PT", "", NOTCH3)
# phenotypes_long <- phenotypes_long[!(phenotypes_long$id %in% Remove_cases), ]

Summary_phenotypes <- phenotypes_long %>%
  dplyr::group_by(group, phecode) %>%
  dplyr::summarise(count=sum(phenotype)) %>%
  dplyr::arrange(desc(count)) %>%
  ungroup()

#sanity check using one phecode (433.2) as example.
Summary_phenotypes[Summary_phenotypes$phecode == "X433.2", ]

Summary_phenotypes_merge <- merge(Summary_phenotypes[Summary_phenotypes$group == 1, ], Summary_phenotypes[Summary_phenotypes$group == 0, ], by = "phecode", all = TRUE)

sum(is.na(Summary_phenotypes_merge)) # 0
Summary_phenotypes_merge[Summary_phenotypes_merge == 0] <- 0.1
length(unique(Summary_phenotypes_merge$phecode)) #1857 phecodes

Summary_phenotypes_merge$delta_change <- (Summary_phenotypes_merge$count.x-(Summary_phenotypes_merge$count.y/2))/(Summary_phenotypes_merge$count.y/2) #2 fold
Summary_phenotypes_merge$delta <- (Summary_phenotypes_merge$count.x-(Summary_phenotypes_merge$count.y/2))

#for v1.2
phecode_mapping <- read.csv("Phecode_map_v1_2_icd10cm_beta.csv", header = T, stringsAsFactors = F)
phecode_mapping_unique <- phecode_mapping %>%
  dplyr::select(phecode, phecode_str, exclude_range, exclude_name) %>%
  distinct(phecode, .keep_all = T)
str(phecode_mapping_unique)
#optional
phecode_mapping_unique$phecode <- paste0("X", phecode_mapping_unique$phecode)
# phecode_mapping_unique$phecode <- gsub("X", "", phecode_mapping_unique$phecode)

Summary_phenotypes_merge_anno <- merge(Summary_phenotypes_merge, phecode_mapping_unique, by = "phecode", all.x = T) #1853

#for a specific ratio, 2 will be changed into the exact ratio between noncarriers and carriers
Summary_phenotypes_merge_anno$count.y <- Summary_phenotypes_merge_anno$count.y/2

N = '<total number of noncarriers>'
Summary_phenotypes_merge_anno$frequency <- ifelse(Summary_phenotypes_merge_anno$count.y*2/N > 0.01, "common", "rare")

#increase in common disease (>1%)
Summary_phenotypes_merge_anno_common <- Summary_phenotypes_merge_anno[Summary_phenotypes_merge_anno$frequency == "common", ]
Summary_phenotypes_merge_anno_common <- Summary_phenotypes_merge_anno_common[!grepl(",", Summary_phenotypes_merge_anno_common$exclude_range), ]
#increase in rare disease (<=1%)
Summary_phenotypes_merge_anno_rare <- Summary_phenotypes_merge_anno[Summary_phenotypes_merge_anno$frequency == "rare", ]
Summary_phenotypes_merge_anno_rare <- Summary_phenotypes_merge_anno_rare[!grepl(",", Summary_phenotypes_merge_anno_rare$exclude_range), ]

save(Summary_phenotypes_merge_anno_common, Summary_phenotypes_merge_anno_rare, file = "<path_to_directory>/Summary_phenotypes109_merge_anno_common_rare_delta_p123_90K_phenotypes109_black.RData", version =2)

Summary_phenotypes_merge_anno_commonrare <- rbind(Summary_phenotypes_merge_anno_common, Summary_phenotypes_merge_anno_rare)
Summary_phenotypes_merge_anno_commonrare <- data.frame(Summary_phenotypes_merge_anno_commonrare)

save(Summary_phenotypes_merge_anno_commonrare, file = "<path_to_directory>/Summary_phenotypes109_merge_anno_commonrare_delta_p123_90K_phenotypes109_black.RData", version =2)

Summary_phenotypes_merge_anno_delta <- head(arrange(Summary_phenotypes_merge_anno_common[!is.na(Summary_phenotypes_merge_anno_common$phecode_str),], desc(delta_change)), n = 50)
Summary_phenotypes_merge_anno_delta_bottom <- head(arrange(Summary_phenotypes_merge_anno_common[!is.na(Summary_phenotypes_merge_anno_common$phecode_str),], delta_change), n = 50)

#fisher exact test
##make a function to get all for top 50 using fisherexact test
fisherexact <- NULL
P <- NULL
S <- NULL

for (phecode in Summary_phenotypes_merge_anno_delta$phecode){  # the top 50 file
  cat(phecode)
  phenotypes_long_top <- phenotypes_long[phenotypes_long$phecode == phecode, ]
  phenotypes_long_top_summary <- phenotypes_long_top %>%
    group_by(id, group) %>%
    summarise(sum = sum(phenotype, na.rm = TRUE)) %>%
    ungroup
  phenotypes_long_top_summary$bin <- ifelse(phenotypes_long_top_summary$sum >= 1, "disease", "nondisease")
  #table(phenotypes_long_top_summary$group, phenotypes_long_top_summary$bin)
  P <- fisher.test(phenotypes_long_top_summary$group, phenotypes_long_top_summary$bin)$p.value
  S <- fisher.test(phenotypes_long_top_summary$group, phenotypes_long_top_summary$bin)$estimate #this is odds ratio
  x <- c(P, S)
  x <- t(data.frame(x))
  colnames(x) <- c("p.value", "estimate")
  rownames(x) <- phecode
  fisherexact <- rbind(fisherexact, x)
  write.csv(fisherexact, "<path_to_directory>/top_fisherexact_phenotypes109_ratio2_90K.csv")
  saveRDS(fisherexact, file= "<path_to_directory>/top_fisherexact_phenotypes109_ratio2_90K.rds")
}

fisherexact <- data.frame(readRDS("<path_to_directory>/top_fisherexact_phenotypes109_ratio2_90K.rds"))
dim(fisherexact)
fisherexact$phecode <- rownames(fisherexact)
head(fisherexact)
Summary_phenotypes_merge_anno_delta_top_fisherexact <- merge(Summary_phenotypes_merge_anno_delta, fisherexact, by = "phecode")

#for matched
save(Summary_phenotypes_merge_anno_delta_top_fisherexact, file = "<path_to_directory>/Summary_phenotypes_merge_anno_delta_top_fisherexact_phenotypes109_ratio2_90K.RData", version = 2)

