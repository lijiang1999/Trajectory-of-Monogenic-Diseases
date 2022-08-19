###outline:
#1. conducting fisher exact test to determine the association between RVs and the top 50 ranked PheCodes;
#2. creating circoplot using ggplot2;
#3. creating correlation/association matrix file between top ranked PheCodes in either RV carrrier or noncarriers;
#4. preparing files for creating circo plot using circlize;
#5. creating circo plot with multiple layers if nencessary


###########################################################################################################
#Conducting fisher exact test to determine the association between RVs and the top 50 ranked PheCodes;

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


####################################################################################################
##create circoplot using ggplot2

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


################################################################################################################
##create correlation/association matrix file between top ranked PheCodes
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
N = '<number of carriers>' #N should be the same number of carriers
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
#preparing files for creating circo plot using circlize

#Create circo plot using circlize
#library(devtools)
#install_github("jokergoo/circlize")
library(circlize)

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
##now creating circo plot
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
circos.trackPlotRegion(data1_remove$exclude_name, y = data1_remove$delta_change, bg.border = NA)

#circos.link add-on
for(i in 1:nrow(CORR_merge)) {
  print(i)
  circos.link(CORR_merge$exclude_name.x[i], c(CORR_merge$x.x[i], CORR_merge$x.x[i]+0.1), CORR_merge$exclude_name.y[i], c(CORR_merge$x.x[i], CORR_merge$x.x[i]+0.1), lwd = log10(CORR_merge$corr[i]), h.ratio = 0.6, rou = 0.5, w2 = 0.6, col = "#FDAE61")
}
circos.clear()
dev.off()

