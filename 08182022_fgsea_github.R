
#installation of fgsea
library(devtools)
install_github("ctlab/fgsea")

library(fgsea)
###############################################################################
#load the summary statistics of percent relative effect (RR-1) with phecode annotation from the previous analysis 
Summary_phenotypes_merge_anno_common <- NULL
load("Summary_phenotypes109_merge_anno_commonrare_delta_p123_90K_phenotypes109_age.RData")
load("Summary_phenotypes109_merge_anno_commonrare_delta_p123_85K_phenotypes109_age.RData")

load("Summary_phenotypes109_merge_anno_commonrare_delta_p123_90K_phenotypes109_age_WHITE.RData")
load("Summary_phenotypes109_merge_anno_commonrare_delta_p123_85K_phenotypes109_age_WHITE.RData")
load("Summary_phenotypes109_merge_anno_commonrare_delta_p123_90K_phenotypes109_age_BLACK.RData")
load("Summary_phenotypes109_merge_anno_commonrare_delta_p123_85K_phenotypes109_age_BLACK.RData")

load("Summary_phenotypes109_merge_anno_commonrare_delta_p123_90K_phenotypes109_WHITE.RData")
load("Summary_phenotypes109_merge_anno_commonrare_delta_p123_85K_phenotypes109_WHITE.RData")
head(Summary_phenotypes_merge_anno_common)

library(dplyr)
#need to make a change to avoid error
Summary_phenotypes_merge_anno_common$exclude_name[Summary_phenotypes_merge_anno_common$exclude_name == "endocrine/metabolic"] <- "endocrine & metabolic"
Summary_phenotypes_merge_anno_common <- Summary_phenotypes_merge_anno_common[!grepl(",", Summary_phenotypes_merge_anno_common$exclude_range), ]

#remove NA in phecode_string first, then select the top 50 in delta column
Summary_phenotypes_merge_anno_common <- Summary_phenotypes_merge_anno_common[!is.na(Summary_phenotypes_merge_anno_common$phecode_str),] %>%
  filter(!exclude_name %in% c("NULL", ""))
print(nrow(Summary_phenotypes_merge_anno_common))
#572 #for 90K; 478 #for 85K; 697/582 for 90K WHITE/BLACK; 591/522  for 85K WHITE/BLACK

Summary <- Summary_phenotypes_merge_anno_common %>%
  dplyr::count(exclude_name) %>%
  ungroup()
Summary
#90K phenotypes109_age
# exclude_name  n
# 1       circulatory system 69
# 2     congenital anomalies  6
# 3             dermatologic 43
# 4                digestive 56
# 5      endocrine/metabolic 65
# 6            genitourinary 66
# 7            hematopoietic 19
# 8      infectious diseases 11
# 9    injuries & poisonings 21
# 10        mental disorders 30
# 11         musculoskeletal 40
# 12               neoplasms 10
# 13            neurological 21
# 14 pregnancy complications 14
# 15             respiratory 40
# 16            sense organs 35
# 17                symptoms 26


exampleRanks_sevengenes <- Summary_phenotypes_merge_anno_common %>%
  dplyr::select(phecode, delta_change) %>%
  arrange(delta_change) %>%
  mutate(rank = 1:nrow(Summary_phenotypes_merge_anno_common))

exampleRanks_sevengenes$phecode <- as.character(exampleRanks_sevengenes$phecode)
ranklist <- with(exampleRanks_sevengenes,setNames(delta_change,phecode))
ranklist
df_ranklist <- data.frame(ranklist)
df_ranklist$phecode <- rownames(df_ranklist)
df_ranklist <- df_ranklist[order(-ranklist),]
df_ranklist$index <- c(1:nrow(df_ranklist)) 

#get the cut_point, enrichment score max 
#create a table for extracted parameters
Pathwaylist <- NULL
listOfDataFrames <- NULL
list_df <- NULL
df <- NULL
df_max <- NULL
X <- NULL
plot_list <- NULL

#to avoid error
#Summary_phenotypes_merge_anno_common$exclude_name[Summary_phenotypes_merge_anno_common$exclude_name == "endocrine/metabolic"] <- "endocrine & metabolic"

for (j in unique(Summary_phenotypes_merge_anno_common$exclude_name)){
  cat(j)
  Pathwaylist <- as.character(Summary_phenotypes_merge_anno_common[Summary_phenotypes_merge_anno_common$exclude_name == j, ]$phecode)
  
  p <- plotEnrichment(Pathwaylist,
                      ranklist)
  
  df <- p$data
  df_max <- df[df$y == max(df$y), ]
  X <- df_max$x #cutpoint of phecode ranks to the top
  colnames(df_max) <- c("Cut_point", "Enrichment_score_max")
  df_max$Category <- j
  df_max$TotalNum <- nrow(Summary_phenotypes_merge_anno_common[Summary_phenotypes_merge_anno_common$exclude_name == j, ])
   
  enriched_top <- df_ranklist[df_ranklist$index <= X, ]$phecode
  Summary_phenotypes_merge_anno_common_enriched_CS <- Summary_phenotypes_merge_anno_common[Summary_phenotypes_merge_anno_common$phecode %in% enriched_top, ]
  Summary_phenotypes_merge_anno_common_enriched_CS <- Summary_phenotypes_merge_anno_common_enriched_CS[Summary_phenotypes_merge_anno_common_enriched_CS$exclude_name == j, ]
  listOfDataFrames[[j]] <- Summary_phenotypes_merge_anno_common_enriched_CS
  
  df_max$Num <- nrow(Summary_phenotypes_merge_anno_common_enriched_CS[Summary_phenotypes_merge_anno_common_enriched_CS$exclude_name == j, ]) #total number of phecodes
  list_df[[j]] <- df_max
  
  #make plots
  plotEnrichment(Pathwaylist,
                 ranklist)
  plot_list[[j]] = p

}
Summary_phenotypes_merge_anno_common_enriched <- do.call("rbind", listOfDataFrames)
Summary_phenotypes_merge_anno_common_enriched_parameters <- do.call("rbind", list_df)
save(Summary_phenotypes_merge_anno_common_enriched, Summary_phenotypes_merge_anno_common_enriched_parameters, file = "Summary_phenotypes_merge_anno_common_enriched_parameters_65yrs_85K_WHITE.RData", version =2)


#create ranking plot per disease category for enrichment analyses
dir.create("path_to_directory/fgsea_85K_WHITE/")
for (i in unique(Summary_phenotypes_merge_anno_common$exclude_name)) {
  tiff(paste('path_to_directory/fgsea_85K_WHITE/', i, '.tiff', sep=''), units="in", width=10, height=8, res=300)
  print(plot_list[[i]])
  dev.off()
}

####################################################################################################
#Create candidate phecode sets for testing
Pathwaylist <- NULL
for (j in c('circulatory system', 'congenital anomalies', 'genitourinary', 'pregnancy complications', 'sense organs', 'symptoms', 'musculoskeletal', 'neoplasms', 'neurological', 'infectious diseases', 'injuries & poisonings', 'mental disorders', 'respiratory', 'dermatologic', 'endocrine/metabolic', 'hematopoietic', 'digestive')){
  cat(j)
  score <- as.character(Summary_phenotypes_merge_anno_common_enriched[Summary_phenotypes_merge_anno_common_enriched$exclude_name == j, ]$phecode)
  Pathwaylist[[j]] <- score
}
Pathwaylist


#using fgseaRes for enrichment analysis
fgseaRes <- fgsea(pathways = Pathwaylist, 
                  stats = ranklist,
                  minSize=1,
                  maxSize=60,
                  nperm=100000)
df_fgseaRes_90K <- data.frame(fgseaRes)
df_fgseaRes_85K <- data.frame(fgseaRes)
save(df_fgseaRes_90K, df_fgseaRes_85K, file = "df_fgseaRes_final_WHITE.RData", version = 2)


#using fgseaMultilevel for enrichment analysis (https://rdrr.io/bioc/fgsea/man/fgseaMultilevel.html)
fgseaMultilevelRes <- fgseaMultilevel(Pathwaylist, ranklist, minSize=1, maxSize=500)
df_fgseaMultilevelRes_90K <- fgseaMultilevelRes
#df_fgseaMultilevelRes_85K <- fgseaMultilevelRes
save(df_fgseaMultilevelRes_90K, df_fgseaMultilevelRes_85K, file = "df_fgseaMultilevelRes_final_WHITE.RData", version = 2)

topPathwaysUp <- fgseaMultilevelRes[ES > 0][head(order(pval), n=18), pathway]
topPathwaysDown <- fgseaMultilevelRes[ES < 0][head(order(pval), n=18), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

#plot the result from fgseaRes
tiff("fgseaplot_85K_phenotypes109_age_WHITE.tiff", units="in", width=10, height=10, res=600)
plotGseaTable(Pathwaylist[topPathways], ranklist, fgseaRes,
              gseaParam=0.5)
dev.off()

#plot the result from fgseaMultilevelRes
tiff("<path_to_directory/fgsea_90K_WHITE/fgseaMultilevelResplot_85K_90KsetWHITE_phenotypes109_age_WHITE.tiff", units="in", width=10, height=10, res=600)
plotGseaTable(Pathwaylist[topPathways], ranklist, fgseaMultilevelRes,
              gseaParam=0.5)
dev.off()


#To save the results in a text format data:table::fwrite function can be used:
fwrite(fgseaRes, file="df_fgseaRes_85K_WHITE.txt", sep="\t", sep2=c("", " ", ""))
fwrite(fgseaRes, file="df_fgseaRes_90K_WHITE.txt", sep="\t", sep2=c("", " ", ""))
fwrite(fgseaMultilevelRes, file="df_fgseaMultilevelRes_85K_WHITE.txt", sep="\t", sep2=c("", " ", ""))
fwrite(fgseaMultilevelRes, file="df_fgseaMultilevelRes_90K_WHITE.txt", sep="\t", sep2=c("", " ", ""))

#######################################################################################################
#create heatmap/bubble plot for PSEA

summary_statistics <- read.csv("06122022_PSEA_table.csv", header = T, stringsAsFactors = F)
head(summary_statistics)
unique(summary_statistics$PheCode_set_List)
bb <- c(1, 2, 4, 6, 8) # define breaks.
ll <- c("0.1", "0.01","0.0001","0.000001","0.00000001") # labels.

q <- summary_statistics %>%
  #arrange(desc(pop)) %>%
  #mutate(country = factor(country, country)) %>%
  ggplot(aes(x=PheCode_set_List, y=Rank_List, size=-log10(padj), fill=NES)) +
  geom_point(alpha=0.9, shape=21) +
  #scale_size(range = c(.1, 24), name="Population (M)") +
  scale_fill_viridis() +
  #scale_fill_viridis(discrete=FALSE, guide=TRUE, option="A") +
  scale_size_continuous(name = "p.adj",
                        breaks = bb,
                        limits = c(0, 10),
                        labels = ll,
                        range = c(0, 20) ) +
  theme_ipsum() +
  theme(legend.position="bottom") +
  xlab("PheCode Set List") +
  ylab("Rank List") +
  theme(legend.position = "right", text = element_text(size=10),) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.text.y = element_text(angle = 45, vjust = 1, hjust=1),
        axis.title.x = element_text(color = "grey20", size = 14, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 14, angle = 90, hjust = .5, vjust = .5, face = "plain")
        )
tiff("<path_to_directory>/fgsea_90K_WHITE/bubble_plot_06212022.tiff", units="in", width=10, height=10, res=600)
q
dev.off()



