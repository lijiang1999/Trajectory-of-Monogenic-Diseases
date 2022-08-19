
###########################################################################################################
#load("Summary_phenotypes109_merge_anno_commonrare_delta_p123_90K_phenotypes109.RData")
#load("Summary_phenotypes109_merge_anno_commonrare_delta_p123_90K_phenotypes109_age.RData")
#load("Summary_phenotypes109_merge_anno_commonrare_delta_p123_85K_phenotypes109_age.RData")
load("Summary_phenotypes109_merge_anno_commonrare_delta_p123_85K_phenotypes109.RData")
Summary_phenotypes_merge_anno_common #this object loaded

#fisher exact test
##make a function to get all for top 50 using fisherexact test
fisherexact <- NULL
P <- NULL
S <- NULL

for (phecode in Summary_phenotypes_merge_anno_common$phecode){  # the top 50 phecodes
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
  write.csv(fisherexact, "ALL_fisherexact_phenotypes109_ratio2_85K.csv")
  saveRDS(fisherexact, file= "ALL_fisherexact_phenotypes109_ratio2_85K.rds")
}

fisherexact <- data.frame(readRDS("ALL_fisherexact_phenotypes109_ratio2_85K.rds"))
dim(fisherexact)

fisherexact$phecode <- rownames(fisherexact)
Summary_phenotypes_merge_anno_delta_top_fisherexact <- merge(Summary_phenotypes_merge_anno_common, fisherexact, by = "phecode")
write.table(Summary_phenotypes_merge_anno_delta_top_fisherexact, "ALL_fisherexact_phenotypes109_ratio2_85K_common.txt", sep = "|", col.names = T, row.names = F, quote = F)


######################################################################################################
load("phenotypes109.RData")
#load("12092021_SKAT_inputfile_p123_90K_updated_phenotypes109.RData") #to get phenotype_ratio2 file
load("12092021_SKAT_inputfile_p123_85K_updated_phenotypes109.RData") #to get phenotype_ratio2 file

phenotypes109_ratio2 <- phenotypes109[phenotypes109$id %in% phenotypes_ratio2$id, ] 

load("phenotypes109_age.RData")
#load("12092021_SKAT_inputfile_p123_90K_updated_phenotypes109_age.RData") #to get phenotype_ratio2 file
load("12092021_SKAT_inputfile_p123_85K_updated_phenotypes109_age.RData") #to get phenotype_ratio2 file
#optional for age cutoff
phenotypes109_ratio2 <- phenotypes109_age[phenotypes109_age$id %in% phenotypes_ratio2$id, ]

phenotypes109_ratio2$id

#optional race specific analyses
unique(dat_RECUR_casecontrol_select_combined_85K$PT_RACE)
# [1] White                                     Unknown                                  
# [3] Black Or African American                 Asian                                    
# [5] Native Hawaiian Or Other Pacific Islander American Indian Or Alaska Native         
# 7 Levels: American Indian Or Alaska Native Asian ... White

WHITE <- dat_RECUR_casecontrol_select_combined_85K[dat_RECUR_casecontrol_select_combined_85K$PT_RACE == "White", ]$PT_ID  #7834; #4520
WHITE <- as.integer(gsub("PT", "", WHITE))

BLACK <- dat_RECUR_casecontrol_select_combined_85K[dat_RECUR_casecontrol_select_combined_85K$PT_RACE == "Black Or African American", ]$PT_ID  #342; #500
BLACK <- as.integer(gsub("PT", "", BLACK))

phenotypes_ratio2 <- phenotypes109_ratio2[phenotypes109_ratio2$id %in% WHITE, ]  #change WHITE to BLACK
table(phenotypes_ratio2$group)
# case control 
# 2593    5241
5241/2593 #2.021211

# case control 
# 1484    3036 
3036/1484 #2.045822

# case control 
# 131     211
211/131  #1.610687

# case control 
# 201     299 
299/201 #1.487562

#optional non race specific analyses
phenotypes_ratio2 <- phenotypes109_ratio2 #from discovery or replication


colnames(phenotypes_ratio2)
gathercols <- colnames(phenotypes_ratio2[, 2:1858]) 
keycol <- "phecode"
valuecol <- "phenotype"
phenotypes_long <- data.frame(gather_(phenotypes_ratio2, keycol, valuecol, gathercols, factor_key=TRUE))
head(phenotypes_long) 
str(phenotypes_long)
unique(phenotypes_long$phenotype)


#phenotypes_matched_ratio2
head(phenotypes_long)
# id   group phecode phenotype
# 1  64 control     008         0
# 2 185 control     008         0
# 3 221    case     008         0
# 4 319 control     008         0
# 5 467 control     008         0
# 6 517 control     008         0

Summary_phenotypes <- phenotypes_long %>%
  dplyr::group_by(group, phecode) %>%
  dplyr::summarise(count=sum(phenotype)) %>%
  dplyr::arrange(desc(count)) %>%
  ungroup()

Summary_phenotypes_merge <- merge(Summary_phenotypes[Summary_phenotypes$group == "case", ], 
                                  Summary_phenotypes[Summary_phenotypes$group == "control", ], 
                                  by = "phecode", all = TRUE)

sum(is.na(Summary_phenotypes_merge)) # 0
Summary_phenotypes_merge[Summary_phenotypes_merge == 0] <- 0.1
length(unique(Summary_phenotypes_merge$phecode)) #1857 phecodes

#for two ratio or five ratio 
Summary_phenotypes_merge$delta_change <- (Summary_phenotypes_merge$count.x-(Summary_phenotypes_merge$count.y/2.045822))/(Summary_phenotypes_merge$count.y/2.045822) #2.045822 fold
Summary_phenotypes_merge$delta <- (Summary_phenotypes_merge$count.x-(Summary_phenotypes_merge$count.y/2.045822))
#optional
phecode_mapping_unique$phecode <- paste0("X", phecode_mapping_unique$phecode)
phecode_mapping_unique$phecode <- gsub("X", "", phecode_mapping_unique$phecode)

#for v1.2
phecode_mapping <- read.csv("Phecode_map_v1_2_icd10cm_beta.csv", header = T, stringsAsFactors = F)
phecode_mapping_unique <- phecode_mapping %>%
  dplyr::select(phecode, phecode_str, exclude_range, exclude_name) %>%
  distinct(phecode, .keep_all = T)  #1755*4
str(phecode_mapping_unique)

Summary_phenotypes_merge_anno <- merge(Summary_phenotypes_merge, phecode_mapping_unique, by = "phecode", all.x = T) #1853

#for matched ratio = 2
Summary_phenotypes_merge_anno$count.y <- Summary_phenotypes_merge_anno$count.y/2.045822


#race-specific
Summary_phenotypes_merge_anno$frequency <- ifelse(Summary_phenotypes_merge_anno$count.y*2.045822/3036 > 0.01, "common", "rare")


#increase in common disease (>1%)
Summary_phenotypes_merge_anno_common <- Summary_phenotypes_merge_anno[Summary_phenotypes_merge_anno$frequency == "common", ]
Summary_phenotypes_merge_anno_common <- Summary_phenotypes_merge_anno_common[!grepl(",", Summary_phenotypes_merge_anno_common$exclude_range), ]
#increase in rare disease (<=1%)
Summary_phenotypes_merge_anno_rare <- Summary_phenotypes_merge_anno[Summary_phenotypes_merge_anno$frequency == "rare", ]
Summary_phenotypes_merge_anno_rare <- Summary_phenotypes_merge_anno_rare[!grepl(",", Summary_phenotypes_merge_anno_rare$exclude_range), ]

save(Summary_phenotypes_merge_anno_common, Summary_phenotypes_merge_anno_rare, file = "Summary_phenotypes109_merge_anno_commonrare_delta_p123_85K_phenotypes109_age_WHITE.RData", version =2)
#save(Summary_phenotypes_merge_anno_common, Summary_phenotypes_merge_anno_rare, file = "Summary_phenotypes109_merge_anno_commonrare_delta_p123_85K_phenotypes109_age_BLACK.RData", version =2)

Summary_phenotypes_merge_anno_commonrare <- rbind(Summary_phenotypes_merge_anno_common, Summary_phenotypes_merge_anno_rare)
Summary_phenotypes_merge_anno_commonrare <- data.frame(Summary_phenotypes_merge_anno_commonrare)
Summary_phenotypes_merge_anno_commonrare$phecode <- paste0("X", Summary_phenotypes_merge_anno_commonrare$phecode)

write.table(Summary_phenotypes_merge_anno_commonrare, file = "Summary_phenotypes109_merge_anno_commonrare_delta_p123_85K_phenotypes109_age_WHITE.txt", sep = "|", quote = F, row.name = F, col.name = T)
#write.table(Summary_phenotypes_merge_anno_commonrare, file = "Summary_phenotypes109_merge_anno_commonrare_delta_p123_85K_phenotypes109_age_BLACK.txt", sep = "|", quote = F, row.name = F, col.name = T)

Summary_phenotypes_merge_anno_delta <- head(arrange(Summary_phenotypes_merge_anno_common[!is.na(Summary_phenotypes_merge_anno_common$phecode_str),], desc(delta_change)), n = 50)
Summary_phenotypes_merge_anno_delta_bottom <- head(arrange(Summary_phenotypes_merge_anno_common[!is.na(Summary_phenotypes_merge_anno_common$phecode_str),], delta_change), n = 50)

Summary_phenotypes_merge_anno_rare_delta <- head(arrange(Summary_phenotypes_merge_anno_rare[!is.na(Summary_phenotypes_merge_anno_rare$phecode_str),], desc(delta)), n = 50)
Summary_phenotypes_merge_anno_rare_delta_bottom <- head(arrange(Summary_phenotypes_merge_anno_rare[!is.na(Summary_phenotypes_merge_anno_rare$phecode_str),], delta), n = 50)

load("Summary_phenotypes109_merge_anno_commonrare_delta_p123_90K_phenotypes109.RData")
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
  write.csv(fisherexact, "top_fisherexact_phenotypes109_age_ratio2_WHITE_85K.csv")
  saveRDS(fisherexact, file= "top_fisherexact_phenotypes109_age_ratio2_WHITE_85K.rds")
}

fisherexact <- data.frame(readRDS("top_fisherexact_phenotypes109_age_ratio2_WHITE_85K.rds"))
dim(fisherexact)

fisherexact$phecode <- rownames(fisherexact)
Summary_phenotypes_merge_anno_delta_top_fisherexact <- merge(Summary_phenotypes_merge_anno_delta, fisherexact, by = "phecode")

#for matched
#save(Summary_phenotypes_merge_anno_delta_top_fisherexact, file = "Summary_phenotypes_merge_anno_delta_top_fisherexact_phenotypes109_ratio2.RData", version = 2)
save(Summary_phenotypes_merge_anno_delta_top_fisherexact, file = "Summary_phenotypes_merge_anno_delta_top_fisherexact_phenotypes109_age_ratio2_BLACK_85K.RData", version = 2)
# save(Summary_phenotypes_merge_anno_delta_top_fisherexact, file = "Summary_phenotypes_merge_anno_delta_top_fisherexact_phenotypes109_ratio2_WHITE_85K.RData", version = 2)
save(Summary_phenotypes_merge_anno_delta_top_fisherexact, file = "Summary_phenotypes_merge_anno_delta_top_fisherexact_phenotypes109_age_ratio2_WHITE_85K.RData", version = 2)

#create circoplot for using ggplot2

#load("Summary_phenotypes_merge_anno_delta_top_chisquare_phenotypes109_ratio2.RData")
data1 <- Summary_phenotypes_merge_anno_delta_top_fisherexact

colnames(data1)
# [1] "phecode"       "group.x"       "count.x"       "group.y"       "count.y"       "delta_change"  "delta"        
# [8] "phecode_str"   "exclude_range" "exclude_name"  "frequency"     "p.value"       "statistic"

data1 <- data1 %>% arrange(exclude_name, phecode)

library(tidyverse)

to_add <- matrix(NA, empty_bar, ncol(data1))
colnames(to_add) <- colnames(data1)
data1 <- rbind(data1, to_add)
data1$id <- seq(1, nrow(data1))

#skip this if no space optional step
empty_bar <- 4 #4
to_add <- data.frame(matrix(0, empty_bar*nlevels(data1$exclude_name), ncol(data1)))
colnames(to_add) <- colnames(data1)
to_add$exclude_name <- rep(levels(data1$exclude_name), each=empty_bar)
data1 <- rbind(data1, to_add)

data1 <- data1 %>% arrange(exclude_name)
data1$id <- seq(1, nrow(data1))

label_data <- data1
number_of_bar <- nrow(label_data)
# I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)
colnames(label_data)

# Make the plot
#matched ratio2 the height of y axis (projection)
data11 <- data1 #create a temp file
data11$delta_change[data11$delta_change >= 3] <- 3  #delta_change is percent relative effect (RR-1), which was controled less and equal to 3 for plotting.
tiff(paste('TRAJECTORY_Circoplot_top50_matched_phenotypes109_age_WHITE_ratio2_fisherexact_85K', '.tiff', sep=''), units="in", width=14, height=10, res=300)

p <- ggplot(data11, aes(x=as.factor(id), y=delta_change*100, fill = as.factor(exclude_name))) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  geom_bar(stat="identity") +
  geom_point(aes(x=as.factor(id), y=-10, size=-log10(p.value)*10), stat="identity", alpha = 0.2) +
  ylim(-100,300) +  #100 change to 150 depends to 200
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar(start = 0) + 
  geom_text(data=label_data, aes(x=id, y=-delta_change + 10, label=phecode_str, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=3.5, angle= label_data$angle, inherit.aes = FALSE) 

p 
dev.off()

#Create circo plot using circlize
#library(devtools)
#install_github("jokergoo/circlize")
library(circlize)

################################################################################################################
#create mat file
#select this function to get all cells filled to avoid error message
fishermatrix <- function(x) {
  names = colnames(x);  num = length(names)
  m = matrix(nrow=num,ncol=num,dimnames=list(names,names))
  for (i in 1:num) {
    for (j in 1:num) {
      #browser()
      m[j,i] = fisher.test(x[, i, drop = TRUE],x[, j, drop = TRUE])$p.value
    }
  }
  return (m)
}

#phenotypes109_age

#optional
phenotypes109_ratio2 <- phenotypes109_ratio2[phenotypes109_ratio2$id %in% WHITE, ]  #change from BLACK to WHITE

#phenotypes109_ratio2_top <- phenotypes109_ratio2[phenotypes109_ratio2$group == "case", colnames(phenotypes109_ratio2) %in% Summary_phenotypes_matched_merge_anno_top_chisquare$phecode]
phenotypes109_ratio2_top <- phenotypes109_ratio2[phenotypes109_ratio2$group == "case", colnames(phenotypes109_ratio2) %in% Summary_phenotypes_merge_anno_delta_top_fisherexact$phecode]

colnames(phenotypes109_ratio2_top)
phenotypes109_ratio2_top <- phenotypes109_ratio2_top %>%
  mutate_if(is.logical, as.factor)
mat <- fishermatrix(phenotypes109_ratio2_top)
diag(mat) <- 1
mat_transformed <- round(-log10(mat), 0)

#phenotypes109_ratio2_top_control <- phenotypes109_ratio2[phenotypes109_ratio2$group == "control", colnames(phenotypes109_ratio2) %in% Summary_phenotypes_matched_merge_anno_top_chisquare$phecode]
phenotypes109_ratio2_top_control <- phenotypes109_ratio2[phenotypes109_ratio2$group == "control", colnames(phenotypes109_ratio2) %in% Summary_phenotypes_merge_anno_delta_top_fisherexact$phecode]

phenotypes109_ratio2_top_control <- phenotypes109_ratio2_top_control %>%
  mutate_if(is.logical, as.factor)

#random select the same size of matched noncarriers to create correlation matrix
set.seed(1234)
N = 2738 #N should be the same number of carriers
picked = sample(seq_len(nrow(phenotypes109_ratio2_top_control)),size = N) 
development = phenotypes109_ratio2_top_control[picked,]
holdout = phenotypes109_ratio2_top_control[-picked,]
mat_control <- fishermatrix(development)
diag(mat_control) <- 1
mat_control_transformed <- round(-log10(mat_control), 0)

#create correlation links between pheocodes in carriers
colnames(data1_remove)
CORR <- data.frame(row=rownames(mat_transformed)[row(mat_transformed)], col=colnames(mat_transformed)[col(mat_transformed)], corr=c(mat_transformed))
#for circo.plot
CORR <- CORR %>%
  filter(corr >3) #p.value < 0.001 which is -log transformed to 3
#for cytoscape may not use
CORR <- CORR %>%
  filter(corr >1) #remove any p.value > 0.1 1;

colnames(CORR)[1] <- "phecode"
CORR_merge <- merge(CORR, data1_remove[, c("phecode", "exclude_name", "x")], by = "phecode")
colnames(CORR_merge)[1] <- "row"
colnames(CORR_merge)[2] <- "phecode"
CORR_merge <- merge(CORR_merge, data1_remove[, c("phecode", "exclude_name", "x")], by = "phecode")
colnames(CORR_merge)[1] <- "col"
head(CORR_merge)
# col    row corr        exclude_name.x x.x exclude_name.y x.y
# 1   195    198   45             neoplasms   4      neoplasms   2
# 2   195  288.1    6         hematopoietic   5      neoplasms   2
# 3   195  195.1   94             neoplasms   3      neoplasms   2
# 4   195 288.11    8         hematopoietic   6      neoplasms   2
# 5   195  963.1   15 injuries & poisonings   4      neoplasms   2
# 6 195.1    195   94             neoplasms   2      neoplasms   3

CORR_control <- data.frame(row=rownames(mat_control_transformed)[row(mat_control_transformed)], col=colnames(mat_control_transformed)[col(mat_control_transformed)], corr=c(mat_control_transformed))
#for circo.plot
CORR_control <- CORR_control %>%
  filter(corr >3) #p.value > 0.001 3
#for cytoscape may not use
CORR_control <- CORR_control %>%
  filter(corr >1) #remove any p.value > 0.1 1;

colnames(CORR_control)[1] <- "phecode"
CORR_control_merge <- merge(CORR_control, data1_remove[, c("phecode", "exclude_name", "x")], by = "phecode")
colnames(CORR_control_merge)[1] <- "row"
colnames(CORR_control_merge)[2] <- "phecode"
CORR_control_merge <- merge(CORR_control_merge, data1_remove[, c("phecode", "exclude_name", "x")], by = "phecode")
colnames(CORR_control_merge)[1] <- "col"

#exclude the connection in control samples (optional)
CORR_control_merge$connection <- paste0(CORR_control_merge$col,"-", CORR_control_merge$row)
CORR_merge$connection <- paste0(CORR_merge$col,"-", CORR_merge$row)
CORR_merge_update <- CORR_merge[!CORR_merge$connection %in% CORR_control_merge$connection, ]
########################################################################################################
#preparing files for creating circo plot
str(data1)
data1 <- data1 %>% group_by(exclude_name) %>% mutate(index = row_number())
normal_sector_index = unique(data1$exclude_name)
xrange = tapply(data1$phecode, data1$exclude_name, function(x) length(x))
xrange
# circulatory system    congenital anomalies            dermatologic               digestive 
# 6                       2                       1                       2 
# endocrine/metabolic           genitourinary           hematopoietic     infectious diseases 
# 3                       9                       7                       1 
# injuries & poisonings        mental disorders         musculoskeletal               neoplasms 
# 4                       3                       4                       6 
# pregnancy complications 
# 2 
sector.width = xrange[normal_sector_index] / sum(xrange[normal_sector_index])
sector.width
# circulatory system    congenital anomalies            dermatologic               digestive 
# 0.12                    0.04                    0.02                    0.04 
# endocrine/metabolic           genitourinary           hematopoietic     infectious diseases 
# 0.06                    0.18                    0.14                    0.02 
# injuries & poisonings        mental disorders         musculoskeletal               neoplasms 
# 0.08                    0.06                    0.08                    0.12 
# pregnancy complications 
# 0.04 
#par(mfrow = c(1, 2))

data1 <- data1 %>%
  group_by(exclude_name)%>%
  mutate(x=1,
         x=cumsum(x))

#remove any sections with only one item; 
#phenotypes109_age
data1_remove <- data1[!(data1$exclude_name %in% c("dermatologic", "digestive", "infectious diseases", "injuries & poisonings", "respiratory", "sense organs")), ]
normal_sector_index = unique(data1_remove$exclude_name)
xrange = tapply(data1_remove$phecode, data1_remove$exclude_name, function(x) length(x))
xrange
sector.width = xrange[normal_sector_index] / sum(xrange[normal_sector_index])
sector.width

labels=data1_remove$phecode_str

########################################################################################################
#now create circo plot
circos.clear()
tiff("Rplot_circoplot_phenotypes109_age_90K_WHITE_fisherexact.tiff", units="in", width=12, height=12, res=300)

circos.par(start.degree = 90, 
           points.overflow.warning = FALSE,
           "track.height" = 0.25, gap.degree = 8, 
           cell.padding = c(0, 0, 0, 0),
           circle.margin = c(0.2, 0.1, 0.2, 0.2)
)
#circos.initialize
color4 = c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF", 
           "#E6F598","#ABDDA4", "#66C2A5", "#3288BD" ,"#5E4FA2")

circos.initialize(data1_remove$exclude_name, x = data1_remove$x)
circos.track(
  sectors = data1_remove$exclude_name,
  x = data1_remove$x,
  y = data1_remove$delta_change - 0.14,   #need to play with this number (0.14) to get the better fit.
  panel.fun = function(x,y) {
    circos.barplot(
      value = y + 0.14, 
      pos = x+0.1,
      col = color4,
      bar_width=0.3
    )
    sector.index = get.cell.meta.data("sector.index")
    print(sector.index)
    xcenter = get.cell.meta.data("xcenter")
    #print(xcenter)
    ycenter = get.cell.meta.data("ycenter")
    circos.text(xcenter, ycenter, sector.index,
                cex = 0.9,
                niceFacing =T,
                facing="bending.inside",
                adj = c(0.5,-11))
    circos.text(x, y + 0.2,
                data1_remove[data1_remove$exclude_name == CELL_META$sector.index, ]$phecode_str,
                cex = 1.0,
                facing = "clockwise")
    circos.yaxis(side = "left", sector.index = "circulatory system")
  }
)

data1_remove <- data1_remove %>%
  arrange(exclude_name, x)

circos.track(sectors = data1_remove$exclude_name, 
             x = data1_remove$x, 
             y = -log10(data1_remove$p.value), #for Fisherexact test or Chisquare test 
             #y = -log10(data1_remove$BURDEN_pvalue), #for BURDEN test
             "track.height" = 0.15,
             panel.fun = function(x, y) {
               circos.points(x, 1.2, cex = y, pch = 16, col = color4)
             })

# Build the regions of track #1
#bg.col = rgb(0.1,0.1,seq(0,1,0.1),0.4)
circos.trackPlotRegion(data1_remove$exclude_name, y = data1_remove$delta_change, bg.border = NA)

#circos.link add-on
for(i in 1:nrow(CORR_merge)) {
  print(i)
  circos.link(CORR_merge$exclude_name.x[i], c(CORR_merge$x.x[i], CORR_merge$x.x[i]+0.1), CORR_merge$exclude_name.y[i], c(CORR_merge$x.x[i], CORR_merge$x.x[i]+0.1), lwd = log10(CORR_merge$corr[i]), h.ratio = 0.6, rou = 0.5, w2 = 0.6, col = "#FDAE61")
}
circos.clear()
dev.off()

