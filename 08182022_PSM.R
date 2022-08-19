setwd("D:/Geisinger/X19981_backup_10072020/Desktop/MS_project")


save.image("05312022_sevengenes_files.RData", version = 2)
save(FID_90K, FID_85K, file = "FID_90K_FID_85K.RData", version = 2)
save(FAM_90K, FAM_85K, file = "FAM_90K_FAM_85K.RData", version = 2)
#05312022

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

FAM_85K <- FAM_175K[!(FAM_175K$V1_2 %in% FAM_90K$V1_2), ]
FID_90K <- merge(FAM_90K, FAM_175K, by = "V1_2") #90051
head(FID_90K)

FID_90K <- data.frame(FID_90K$V2.y, FID_90K$V2.y) #90051
colnames(FID_90K) <- c("FID", "IID")

FID_85K <- FAM_175K[!(FAM_175K$V2 %in% FID_90K$FID), ] #83534
FID_85K <- data.frame(FID_85K$V2, FID_85K$V2) #90051
colnames(FID_85K) <- c("FID", "IID")
head(FID_85K)
save(FID_90K, FID_85K, file = "FID_for_90K_85K.RData", version = 2)
write.table(FID_90K, file = "FID_90K.txt", col.names = T, row.names = F, quote = F)
write.table(FID_85K, file = "FID_85K.txt", col.names = T, row.names = F, quote = F)

FAM_175K_discovery <- data.frame(FAM_175K[FAM_175K$V1_2 %in% ID_90K_p123, ]$V2, FAM_175K[FAM_175K$V1_2 %in% ID_90K_p123, ]$V2)
colnames(FAM_175K_discovery) <- c("FID", "IID")

FAM_175K_replication <- data.frame(FAM_175K[FAM_175K$V1_2 %in% ID_85K_p123, ]$V2, FAM_175K[FAM_175K$V1_2 %in% ID_85K_p123, ]$V2)
colnames(FAM_175K_replication) <- c("FID", "IID")

write.table(FAM_175K_discovery, file = "FAM_175K_discovery_select.txt", col.names=T, quote=F, row.names=F)
write.table(FAM_175K_replication, file = "FAM_175K_replication_select.txt", col.names=T, quote=F, row.names=F)

######################################################################################################
#adding some patients with rare variants from COL4A1/2 in BLACK and NOTCH in WHITE
Genotype_12032021 <- read.table("GHS_Freeze_175.GL.splitmulti_sevengene_12032021.raw", header =T, stringsAsFactors = F)
colnames(Genotype_12032021)
# [1] "FID"                 "IID"                 "PAT"                 "MAT"                 "SEX"                
# [6] "PHENOTYPE"           "X3.48466629.G.A_G"   "X3.48466945.G.A_G"   "X3.48466995.C.T_C"   "X3.48466996.G.A_G"  
# [11] "X3.48467064.G.C_G"   "X3.48467073.A.G_A"   "X3.48467163.G.A_G"   "X3.48467322.G.A_G"   "X3.48467375.G.C_G"  
# [16] "X3.48467569.A.G_A"   "X10.122462103.C.A_C" "X10.122489438.C.T_C" "X10.122489616.T.C_T" "X10.122506796.G.A_G"
# [21] "X10.122506874.G.A_G" "X10.122506886.G.A_G" "X10.122512066.G.A_G" "X13.110167161.G.C_G" "X13.110192285.C.A_C"
# [26] "X13.110205493.C.T_C" "X13.110209941.T.A_T" "X13.110210116.A.G_A" "X13.110210163.C.T_C" "X13.110211923.C.A_C"
# [31] "X13.110242665.G.A_G" "X13.110307859.G.C_G" "X13.110308119.T.A_T" "X13.110357595.C.T_C" "X13.110424998.G.A_G"
# [36] "X13.110428584.G.A_G" "X13.110428632.T.C_T" "X13.110429953.T.A_T" "X13.110438001.G.T_G" "X13.110439828.C.T_C"
# [41] "X13.110462113.G.A_G" "X13.110462385.G.A_G" "X13.110465576.C.T_C" "X13.110469231.C.T_C" "X13.110491334.C.A_C"
# [46] "X13.110495443.C.T_C" "X13.110503936.C.T_C" "X13.110506597.C.T_C" "X13.110512120.G.A_G" "X17.8231334.C.T_C"  
# [51] "X17.8231981.C.A_C"   "X17.8237392.C.T_C"   "X17.8237487.G.A_G"   "X17.8248035.A.T_A"   "X17.8248036.T.A_T"  
# [56] "X19.15165393.A.G_A"  "X19.15165506.G.A_G"  "X19.15167362.C.T_C"  "X19.15174391.G.A_G"  "X19.15179079.A.C_A" 
# [61] "X19.15179115.G.A_G"  "X19.15180807.G.A_G"  "X19.15181023.T.G_T"  "X19.15181164.T.C_T"  "X19.15185371.G.A_G" 
# [66] "X19.15185632.C.A_C"  "X19.15187126.G.A_G"  "X19.15187273.G.A_G"  "X19.15187315.G.A_G"  "X19.15192130.T.C_T" 
# [71] "X19.15192449.G.A_G"  "X19.15192493.C.G_C"  "X23.101401752.C.T_C" "X23.101403828.G.A_G"

"X13.110357595.C.T_C"
"X13.110491334.C.A_C"
"X13.110512120.G.A_G"
"X19.15181023.T.G_T"
"X19.15192130.T.C_T"
"X3.48466996.G.A_G"
Genotype <- Genotype_12032021[, c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE", "X13.110357595.C.T_C", "X13.110491334.C.A_C", "X13.110512120.G.A_G", "X19.15181023.T.G_T", "X19.15192130.T.C_T", "X3.48466996.G.A_G")
                              ]
save(ALL, file = "Individualgene_PTID_05312022_additional.RData", version = 2)

#testing on X19.15181023.T.G_T, "X19.15192130.T.C_T"
Genotype_NOTCH3 <- Genotype_12032021[, c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE", "X19.15192130.T.C_T")]
                     
Genotype_NOTCH3 <- as_tibble(Genotype_NOTCH3)
Genotype_NOTCH3 <- cSplit(Genotype_NOTCH3, "FID", "_") #FID_2 for PT####
Genotype_NOTCH3 <- Genotype_NOTCH3 %>%
  dplyr::select(starts_with("X19"))  #select all variants

Genotype_NOTCH3$PT_ID <- Genotype$FID_2
head(Genotype_NOTCH3)
gathercols <- colnames(Genotype_NOTCH3)[1:(ncol(Genotype_NOTCH3)-1)] 
keycol <- "genotype"
valuecol <- "dosage"
Genotype_long <- data.frame(gather_(Genotype_NOTCH3, keycol, valuecol, gathercols, factor_key=TRUE))
head(Genotype_long)

#90K or 85K
Genotype_long_mutation <- Genotype_long[!is.na(Genotype_long$dosage),] %>%
  filter(dosage != 2) %>% 
  distinct(PT_ID, .keep_all = T) #
#some patients may have multiple mutations, here distinct to get the unique individual ID

#the number in bracket is same variant strategy
NOTCH3 <- as.character(unique(Genotype_long_mutation$PT_ID)) #245 vs 238(205); 746 vs 680

######################################################################################################
#05/31/2022
#top variants (#297)
Genotype<- read.table("GHS_Freeze_175.GL.splitmulti_sevengene_05312022_90K.raw", header =T, stringsAsFactors = F) #186 variants
Genotype<- read.table("GHS_Freeze_175.GL.splitmulti_sevengene_05312022_85K.raw", header =T, stringsAsFactors = F) #220 variants

colnames(Genotype)
ncol(Genotype)
#90K remove 111
# [1] "FID"                  "IID"                  "PAT"                  "MAT"                 
# [5] "SEX"                  "PHENOTYPE"            "X3.48466629.G.A_G"    "X3.48466750.C.T_C"   
# [9] "X3.48466792.I.1_G"    "X3.48466852.A.G_A"    "X3.48466881.I.8_G"    "X3.48466945.G.A_G"   
# [13] "X3.48466995.C.T_C"    "X3.48466996.G.A_G"    "X3.48467038.G.A_G"    "X3.48467043.G.A_G"   
# [17] "X3.48467128.C.T_C"    "X3.48467247.G.A_G"    "X3.48467322.G.A_G"    "X3.48467524.C.T_C"   
# [21] "X3.48467562.A.C_A"    "X3.48467569.A.G_A"    "X10.122462042.C.A_C"  "X10.122462103.C.A_C" 
# [25] "X10.122488925.C.T_C"  "X10.122489438.C.T_C"  "X10.122489616.T.C_T"  "X10.122506760.G.A_G" 
# [29] "X10.122506796.G.A_G"  "X10.122506874.G.A_G"  "X10.122506886.G.A_G"  "X10.122508766.I.2_C" 
# [33] "X10.122512065.C.T_C"  "X10.122512066.G.A_G"  "X10.122514204.I.13_G" "X10.122514264.G.C_G" 
# [37] "X13.110150412.I.1_T"  "X13.110162434.C.A_C"  "X13.110170705.C.G_C"  "X13.110174509.C.A_C" 
# [41] "X13.110176935.I.8_T"  "X13.110181326.C.T_C"  "X13.110187164.C.A_C"  "X13.110192285.C.G_C" 
# [45] "X13.110192285.C.A_C"  "X13.110201467.G.A_G"  "X13.110205493.C.T_C"  "X13.110210163.C.T_C" 
# [49] "X13.110211919.C.A_C"  "X13.110211923.C.A_C"  "X13.110211924.T.G_T"  "X13.110307859.G.C_G" 
# [53] "X13.110308119.T.A_T"  "X13.110357595.C.T_C"  "X13.110424998.G.A_G"  "X13.110428550.I.1_G" 
# [57] "X13.110428584.G.A_G"  "X13.110428585.T.G_T"  "X13.110429925.I.1_C"  "X13.110429953.T.A_T" 
# [61] "X13.110438001.G.T_G"  "X13.110439828.C.T_C"  "X13.110446797.G.T_G"  "X13.110450406.C.T_C" 
# [65] "X13.110457430.G.A_G"  "X13.110457524.C.A_C"  "X13.110462113.G.A_G"  "X13.110462385.G.A_G" 
# [69] "X13.110465492.G.T_G"  "X13.110465552.G.T_G"  "X13.110469216.G.A_G"  "X13.110469231.C.T_C" 
# [73] "X13.110478060.I.1_C"  "X13.110480375.I.1_A"  "X13.110482614.G.T_G"  "X13.110484910.C.T_C" 
# [77] "X13.110485837.G.A_G"  "X13.110491334.C.A_C"  "X13.110493246.G.A_G"  "X13.110495443.C.T_C" 
# [81] "X13.110503123.I.10_T" "X13.110503273.G.T_G"  "X13.110503873.G.A_G"  "X13.110503927.C.T_C" 
# [85] "X13.110503936.C.T_C"  "X13.110503972.C.T_C"  "X13.110506597.C.T_C"  "X13.110511932.A.G_A" 
# [89] "X13.110512120.G.A_G"  "X17.8228290.C.A_C"    "X17.8228302.G.A_G"    "X17.8228735.G.A_G"   
# [93] "X17.8228807.C.A_C"    "X17.8229306.I.1_C"    "X17.8229905.I.1_A"    "X17.8229943.G.A_G"   
# [97] "X17.8230304.T.C_T"    "X17.8230395.I.1_A"    "X17.8231274.A.C_A"    "X17.8231334.C.T_C"   
# [101] "X17.8231358.C.A_C"    "X17.8231364.C.A_C"    "X17.8231427.G.A_G"    "X17.8231749.G.A_G"   
# [105] "X17.8232162.G.C_G"    "X17.8232222.I.1_T"    "X17.8232360.C.A_C"    "X17.8232904.A.C_A"   
# [109] "X17.8234454.C.T_C"    "X17.8234520.G.A_G"    "X17.8235846.I.1_A"    "X17.8236228.G.A_G"   
# [113] "X17.8236238.C.T_C"    "X17.8236294.A.G_A"    "X17.8236302.C.A_C"    "X17.8236323.C.T_C"   
# [117] "X17.8237374.C.A_C"    "X17.8237392.C.T_C"    "X17.8237427.A.G_A"    "X17.8237487.G.A_G"   
# [121] "X17.8237497.G.A_G"    "X17.8238505.G.A_G"    "X17.8238575.I.4_G"    "X17.8243032.C.T_C"   
# [125] "X17.8248002.A.C_A"    "X17.8248018.G.A_G"    "X17.8248035.A.T_A"    "X17.8248035.A.G_A"   
# [129] "X17.8248036.T.A_T"    "X17.8248036.T.G_T"    "X19.15161059.G.A_G"   "X19.15161525.I.1_C"  
# [133] "X19.15161715.C.A_C"   "X19.15165506.G.A_G"   "X19.15165517.T.C_T"   "X19.15167362.C.T_C"  
# [137] "X19.15170669.A.C_A"   "X19.15170690.G.C_G"   "X19.15177977.I.1_T"   "X19.15178039.G.A_G"  
# [141] "X19.15179079.A.C_A"   "X19.15179115.G.A_G"   "X19.15179397.G.A_G"   "X19.15179498.T.C_T"  
# [145] "X19.15180703.G.T_G"   "X19.15180807.G.A_G"   "X19.15180960.C.T_C"   "X19.15181023.T.G_T"  
# [149] "X19.15181164.T.C_T"   "X19.15181604.C.A_C"   "X19.15181621.T.C_T"   "X19.15181802.C.A_C"  
# [153] "X19.15185256.C.G_C"   "X19.15185371.G.A_G"   "X19.15185393.A.T_A"   "X19.15185404.G.A_G"  
# [157] "X19.15185546.I.1_C"   "X19.15185632.C.A_C"   "X19.15187126.G.A_G"   "X19.15187171.G.A_G"  
# [161] "X19.15187213.G.A_G"   "X19.15187216.T.C_T"   "X19.15187273.G.A_G"   "X19.15187315.G.A_G"  
# [165] "X19.15187956.A.G_A"   "X19.15191496.C.T_C"   "X19.15191793.C.T_C"   "X19.15191838.C.T_C"  
# [169] "X19.15191868.C.G_C"   "X19.15192037.C.G_C"   "X19.15192130.T.C_T"   "X19.15192188.G.C_G"  
# [173] "X19.15192256.I.2_C"   "X19.15192449.G.A_G"   "X19.15192493.C.T_C"   "X19.15192493.C.G_C"  
# [177] "X23.101397874.G.T_G"  "X23.101397915.C.G_C"  "X23.101397942.T.G_T"  "X23.101398032.C.T_C" 
# [181] "X23.101398052.C.T_C"  "X23.101398504.T.C_T"  "X23.101398862.T.C_T"  "X23.101398906.C.T_C" 
# [185] "X23.101398942.T.C_T"  "X23.101400709.A.G_A"  "X23.101401724.T.C_T"  "X23.101401748.C.T_C" 
# [189] "X23.101403941.C.T_C"  "X23.101403983.T.C_T"  "X23.101403984.C.G_C"  "X23.101407795.C.T_C" 

#85K remove 77
# [1] "FID"                  "IID"                  "PAT"                  "MAT"                 
# [5] "SEX"                  "PHENOTYPE"            "X3.48466629.G.A_G"    "X3.48466656.A.G_A"   
# [9] "X3.48466687.T.G_T"    "X3.48466750.C.T_C"    "X3.48466792.I.1_G"    "X3.48466852.A.G_A"   
# [13] "X3.48466881.I.8_G"    "X3.48466945.G.A_G"    "X3.48466947.I.1_C"    "X3.48466995.C.T_C"   
# [17] "X3.48466996.G.A_G"    "X3.48467043.G.A_G"    "X3.48467128.C.T_C"    "X3.48467145.C.T_C"   
# [21] "X3.48467251.I.3_G"    "X3.48467277.I.4_T"    "X3.48467322.G.A_G"    "X3.48467524.C.T_C"   
# [25] "X3.48467562.A.C_A"    "X3.48467569.A.G_A"    "X3.48467572.G.C_G"    "X10.122462103.C.A_C" 
# [29] "X10.122462126.T.C_T"  "X10.122488925.C.T_C"  "X10.122489603.G.A_G"  "X10.122506734.G.A_G" 
# [33] "X10.122506796.G.A_G"  "X10.122506817.C.T_C"  "X10.122506818.G.A_G"  "X10.122506874.G.A_G" 
# [37] "X10.122506886.G.A_G"  "X10.122508758.C.T_C"  "X10.122512065.C.T_C"  "X10.122514204.I.13_G"
# [41] "X10.122514234.C.T_C"  "X10.122514264.G.C_G"  "X13.110155398.C.A_C"  "X13.110162263.G.A_G" 
# [45] "X13.110164886.I.2_G"  "X13.110178921.A.G_A"  "X13.110179352.C.T_C"  "X13.110187164.C.A_C" 
# [49] "X13.110187218.C.G_C"  "X13.110192285.C.G_C"  "X13.110192915.T.C_T"  "X13.110201437.C.A_C" 
# [53] "X13.110201467.G.A_G"  "X13.110205406.C.T_C"  "X13.110205493.C.T_C"  "X13.110209988.C.A_C" 
# [57] "X13.110211868.C.T_C"  "X13.110211923.C.A_C"  "X13.110211924.T.G_T"  "X13.110307858.A.G_A" 
# [61] "X13.110307859.G.C_G"  "X13.110307908.I.4_G"  "X13.110308119.T.A_T"  "X13.110357595.C.T_C" 
# [65] "X13.110424772.I.1_G"  "X13.110424957.I.1_C"  "X13.110424998.G.A_G"  "X13.110424999.T.C_T" 
# [69] "X13.110428584.G.A_G"  "X13.110429925.I.1_C"  "X13.110429939.G.T_G"  "X13.110429953.T.A_T" 
# [73] "X13.110429958.T.A_T"  "X13.110436359.C.T_C"  "X13.110438001.G.T_G"  "X13.110445835.C.T_C" 
# [77] "X13.110446797.G.T_G"  "X13.110450400.G.T_G"  "X13.110450406.C.T_C"  "X13.110450420.I.1_A" 
# [81] "X13.110457430.G.A_G"  "X13.110457524.C.T_C"  "X13.110462113.G.A_G"  "X13.110462385.G.A_G" 
# [85] "X13.110465492.G.T_G"  "X13.110469231.C.T_C"  "X13.110472943.C.T_C"  "X13.110484910.C.T_C" 
# [89] "X13.110485801.G.T_G"  "X13.110485837.G.A_G"  "X13.110491334.C.A_C"  "X13.110493219.G.A_G" 
# [93] "X13.110493246.G.A_G"  "X13.110495443.C.T_C"  "X13.110501673.C.T_C"  "X13.110503123.I.10_T"
# [97] "X13.110503283.G.A_G"  "X13.110503845.A.C_A"  "X13.110503855.G.A_G"  "X13.110503936.C.T_C" 
# [101] "X13.110506597.C.T_C"  "X13.110512120.G.A_G"  "X17.8228302.G.A_G"    "X17.8228726.C.T_C"   
# [105] "X17.8228807.C.A_C"    "X17.8229190.I.1_A"    "X17.8229905.I.1_A"    "X17.8229925.G.A_G"   
# [109] "X17.8229943.G.A_G"    "X17.8230292.A.T_A"    "X17.8230304.T.C_T"    "X17.8230395.I.1_A"   
# [113] "X17.8230469.C.T_C"    "X17.8231291.I.1_T"    "X17.8231364.C.A_C"    "X17.8231427.G.A_G"   
# [117] "X17.8231749.G.A_G"    "X17.8232038.I.1_T"    "X17.8232162.G.C_G"    "X17.8232175.I.1_G"   
# [121] "X17.8232360.C.A_C"    "X17.8232401.C.A_C"    "X17.8234570.I.1_T"    "X17.8234657.T.C_T"   
# [125] "X17.8234859.C.A_C"    "X17.8234889.G.A_G"    "X17.8235132.C.A_C"    "X17.8236189.G.A_G"   
# [129] "X17.8236228.G.A_G"    "X17.8236238.C.T_C"    "X17.8236255.G.A_G"    "X17.8236276.G.A_G"   
# [133] "X17.8236294.A.G_A"    "X17.8236323.C.T_C"    "X17.8237374.C.A_C"    "X17.8237374.C.G_C"   
# [137] "X17.8237392.C.T_C"    "X17.8237473.G.A_G"    "X17.8238505.G.A_G"    "X17.8238550.G.A_G"   
# [141] "X17.8238575.I.4_G"    "X17.8238616.G.A_G"    "X17.8243032.C.T_C"    "X17.8243073.G.A_G"   
# [145] "X17.8248018.G.A_G"    "X17.8248035.A.T_A"    "X17.8248036.T.A_T"    "X17.8248036.T.G_T"   
# [149] "X19.15160724.G.A_G"   "X19.15161059.G.A_G"   "X19.15165366.A.C_A"   "X19.15166049.I.1_G"  
# [153] "X19.15167362.C.T_C"   "X19.15170329.A.C_A"   "X19.15170420.C.T_C"   "X19.15170690.G.C_G"  
# [157] "X19.15174066.A.G_A"   "X19.15174391.G.C_G"   "X19.15177984.C.T_C"   "X19.15178003.G.A_G"  
# [161] "X19.15179079.A.C_A"   "X19.15179115.G.A_G"   "X19.15179397.G.A_G"   "X19.15180103.C.T_C"  
# [165] "X19.15180217.C.T_C"   "X19.15180717.G.A_G"   "X19.15180761.T.C_T"   "X19.15180807.G.A_G"  
# [169] "X19.15180964.G.T_G"   "X19.15181023.T.G_T"   "X19.15181604.C.A_C"   "X19.15181621.T.C_T"  
# [173] "X19.15181802.C.A_C"   "X19.15184922.G.T_G"   "X19.15185256.C.G_C"   "X19.15185350.G.A_G"  
# [177] "X19.15185371.G.A_G"   "X19.15185404.G.A_G"   "X19.15185633.G.T_G"   "X19.15186911.G.A_G"  
# [181] "X19.15187126.G.A_G"   "X19.15187186.G.A_G"   "X19.15187213.G.A_G"   "X19.15187216.T.C_T"  
# [185] "X19.15187273.G.A_G"   "X19.15189064.A.G_A"   "X19.15189114.I.1_T"   "X19.15191496.C.T_C"  
# [189] "X19.15191793.C.T_C"   "X19.15191817.G.A_G"   "X19.15191838.C.T_C"   "X19.15191968.C.T_C"  
# [193] "X19.15192022.C.T_C"   "X19.15192095.G.A_G"   "X19.15192130.T.C_T"   "X19.15192188.G.C_G"  
# [197] "X19.15192493.C.T_C"   "X19.15197503.C.G_C"   "X19.15200859.G.T_G"   "X19.15200905.T.G_T"  
# [201] "X23.101397810.T.A_T"  "X23.101397915.C.G_C"  "X23.101398011.C.T_C"  "X23.101398032.C.T_C" 
# [205] "X23.101398033.G.A_G"  "X23.101398044.G.C_G"  "X23.101398051.C.T_C"  "X23.101398368.A.T_A" 
# [209] "X23.101398431.T.C_T"  "X23.101398468.G.C_G"  "X23.101398495.C.T_C"  "X23.101398501.T.G_T" 
# [213] "X23.101398504.T.C_T"  "X23.101398786.A.G_A"  "X23.101398862.T.C_T"  "X23.101398942.T.C_T" 
# [217] "X23.101400692.G.T_G"  "X23.101400709.A.G_A"  "X23.101400712.A.G_A"  "X23.101403809.A.C_A" 
# [221] "X23.101403845.C.T_C"  "X23.101403846.G.A_G"  "X23.101403864.G.A_G"  "X23.101403941.C.T_C" 
# [225] "X23.101403983.T.C_T"  "X23.101403984.C.G_C" 

table(Genotype$X13.110462385.G.A_G)
# 1     2 
# 2 90049
# 1     2 
# 1 83533
table(Genotype$X19.15179115.G.A_G)
# 1     2 
# 17 90033
# 1     2 
# 7 83526
table(Genotype$X19.15179079.A.C_A)
# 1     2 
# 43 90007
# 1     2 
# 43 83491 

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
  filter(dosage != 2) %>% 
  distinct(PT_ID, .keep_all = T) #
#some patients may have multiple mutations, here distinct to get the unique individual ID

#the number in bracket is same variant strategy
TREX1 <- as.character(unique(Genotype_long_mutation$PT_ID)) #128 vs 118(110); 216 vs 207
HTRA1 <- as.character(unique(Genotype_long_mutation$PT_ID)) #81 vs 101(90); 81 vs 101
COL4A1 <- as.character(unique(Genotype_long_mutation$PT_ID)) #1206 vs 155(122); 1360 vs 447
NOTCH3 <- as.character(unique(Genotype_long_mutation$PT_ID)) #245 vs 238(205); 746 vs 680
CTC1 <- as.character(unique(Genotype_long_mutation$PT_ID)) #178 vs 182(140); 178 vs 182
GLA <- as.character(unique(Genotype_long_mutation$PT_ID)) #44 vs 42(24); 44 vs 42

ALL <- as.character(unique(Genotype_long_mutation$PT_ID)) #2606 vs 1647(689) total = 4253; 
2606+1647

save(TREX1, HTRA1, COL4A1, NOTCH3, CTC1, GLA, ALL, file = "Individualgene_PTID_05312022_90K.RData", version = 2)

df_ALL_90K <- data.frame(ALL)
df_ALL_90K$Cohort <- "discovery"
head(df_ALL_90K)

#save(TREX1, HTRA1, COL4A1, NOTCH3, CTC1, GLA, ALL, file = "Individualgene_PTID_05312022_85K.RData", version = 2) #same variants
save(TREX1, HTRA1, COL4A1, NOTCH3, CTC1, GLA, ALL, file = "Individualgene_PTID_05312022_85K_samegene.RData", version = 2)
df_ALL_85K <- data.frame(ALL)
df_ALL_85K$Cohort <- "replication"
head(df_ALL_85K)

df_ALL <- rbind(df_ALL_90K, df_ALL_85K) #4253
write.table(df_ALL, file = "df_ALL.txt", col.names = T, row.names = F, quote = F)



load("Individualgene_PTID_05312022_90K.RData")
ID <- as.integer(gsub("PT", "", CTC1))

#create bash file to extract QC sequencing data
STLIST <- NULL
for (i in colnames(Genotype_TREX1)[1:(ncol(Genotype_TREX1)-1)]) {
  STLIST[[i]] <- as.integer(unlist(strsplit(i, "\\."))[2]) #Escape special characters
  print(STLIST[i])
  
}
STLIST_final <- data.frame(do.call(rbind,STLIST))
colnames(STLIST_final) <- "coordinate"
STLIST_final$keep <- STLIST_final$coordinate

STLIST_final$coordinate <- paste0("zcat /media/12TB/MyCode_t1/data/pVCF/GL_by_segments/GHS_Freeze_175.967_X_87625896_101892854.GL.vcf.gz | sed '/^#/d' | awk '$2 == ", STLIST_final$coordinate)
STLIST_final$coordinate <- paste0(STLIST_final$coordinate, "' > /media/10TB/MS/test/GHS_Freeze_175.967_X_87625896_101892854.GL.vcf.test.")
STLIST_final$coordinate <- paste0(STLIST_final$coordinate,STLIST_final$keep)

rownames(STLIST_final) <- NULL
STLIST_final <- STLIST_final[, 1]
head(STLIST_final)
write.table(STLIST_final, file = "STLIST_final_GLA.txt", col.names = F, row.names = F, quote = F)
#create bash file
#sed -i -e 's/\r$//' STLIST_final_NOTCH3.sh
#chmod +x STLIST_final_NOTCH3.sh
#./STLIST_final_NOTCH3.sh

STLIST_final$coordinate <- paste0("zcat /media/12TB/MyCode_t1/data/pVCF/GL_by_segments/GHS_Freeze_175.852_19_15185138_15984381.GL.vcf.gz | sed '/^#/d' | awk '$2 == ", STLIST_final$coordinate)
STLIST_final$coordinate <- paste0(STLIST_final$coordinate, "' > /media/10TB/MS/test/GHS_Freeze_175.852_19_15185138_15984381.GL.vcf.test.")
STLIST_final$coordinate <- paste0(STLIST_final$coordinate,STLIST_final$keep)

rownames(STLIST_final) <- NULL
head(STLIST_final)
STLIST_final <- STLIST_final[, 1]
write.table(STLIST_final, file = "STLIST_final_NOTCH3.txt", col.names = F, row.names = F, quote = F)

STLIST_final$coordinate <- paste0("zcat /media/12TB/MyCode_t1/data/pVCF/GL_by_segments/GHS_Freeze_175.767_17_8146255_9245272.GL.vcf.gz | sed '/^#/d' | awk '$2 == ", STLIST_final$coordinate)
STLIST_final$coordinate <- paste0(STLIST_final$coordinate, "' > /media/10TB/MS/test/GHS_Freeze_175.767_17_8146255_9245272.GL.vcf.test.")
STLIST_final$coordinate <- paste0(STLIST_final$coordinate,STLIST_final$keep)

rownames(STLIST_final) <- NULL
head(STLIST_final)
STLIST_final <- STLIST_final[, 1]
write.table(STLIST_final, file = "STLIST_final_CTC1.txt", col.names = F, row.names = F, quote = F)


GHS_Freeze_175.767_17_8146255_9245272.GL.vcf.test.8231981
################################################################################################################
#Ramdom sampling to get controls

FID_90K <- cSplit(FID_90K, "FID", "_")
head(FID_90K)
FID_90K_excludecase <- FID_90K[!(FID_90K$FID_2 %in% df_ALL_90K$ALL), ]  #88181
set.seed(42)
FID_90K_excludecase <- FID_90K_excludecase[sample(nrow(FID_90K_excludecase), 18700), ]  #10times sampling
head(FID_90K_excludecase)

df_ALL_90K_control <- data.frame(FID_90K_excludecase$FID_2)
colnames(df_ALL_90K_control) <- "ALL"
df_ALL_90K_control$Cohort <- "discovery"

FID_85K <- cSplit(FID_85K, "FID", "_")
head(FID_85K)
FID_85K_excludecase <- FID_85K[!(FID_85K$FID_2 %in% df_ALL_85K$ALL), ]  #82702
set.seed(42)
FID_85K_excludecase <- FID_85K_excludecase[sample(nrow(FID_85K_excludecase), 8320), ]  #10times sampling
head(FID_85K_excludecase)

df_ALL_85K_control <- data.frame(FID_85K_excludecase$FID_2)
colnames(df_ALL_85K_control) <- "ALL"
df_ALL_85K_control$Cohort <- "replication"

df_ALL_control <- rbind(df_ALL_90K_control, df_ALL_85K_control) #27020
write.table(df_ALL_control, file = "df_ALL_control.txt", col.names = T, row.names = F, quote = F)

############################################################################################
#propensity score matching
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
ID_90K_p123 <- dat_RECUR_ICD109_demographics[dat_RECUR_ICD109_demographics$PT_ID %in% FAM_90K$V1_2, ]$PT_ID #21127
ID_85K_p123 <- dat_RECUR_ICD109_demographics[!(dat_RECUR_ICD109_demographics$PT_ID %in% FAM_90K$V1_2), ]$PT_ID #9787
length(ID_90K_p123)

dat_RECUR_ICD109_demographics_90K <- dat_RECUR_ICD109_demographics[dat_RECUR_ICD109_demographics$PT_ID %in% ID_90K_p123, ]
dat_RECUR_ICD109_demographics_85K <- dat_RECUR_ICD109_demographics[dat_RECUR_ICD109_demographics$PT_ID %in% ID_85K_p123, ]
save(dat_RECUR_ICD109_demographics_90K, dat_RECUR_ICD109_demographics_85K, file = "dat_RECUR_ICD109_demographics_90K_85K.RData", version = 2)

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
dat_RECUR_ICD109_demographics_90K <- rbind(dat_RECUR_ICD109_demographics_90K.case, dat_RECUR_ICD109_demographics_90K.control) #8235
head(dat_RECUR_ICD109_demographics_90K)
save(dat_RECUR_ICD109_demographics_90K, file = "dat_RECUR_ICD109_demographics_90K_ratio2.RData", version = 2)

#same way to create replication dataset(85K)
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
#PCA analysis

load("dat_RECUR_ICD109_demographics_90K_ratio2.RData") #7812
head(dat_RECUR_ICD109_demographics_90K)
PCA <- read.table("GHS_Freeze_175.GL.splitmulti_FAM_175K_discovery_ratio2_EXOME2_pca.eigenvec", header = F, stringsAsFactors = F)
head(PCA) #7812*12
dim(PCA)
PCA$cohort <- "discovery"
PCA <- cSplit(PCA, "V1", "_")
colnames(PCA)[which(names(PCA) == "V1_2")] <- "PT_ID"
dat_RECUR_ICD109_demographics_90K_PCA <- merge(dat_RECUR_ICD109_demographics_90K, PCA, by = "PT_ID")
FAM_90K_discovery_ratio2 <- data.frame(PCA$V2, PCA$V2)
save(dat_RECUR_ICD109_demographics_90K_PCA, file = "dat_RECUR_ICD109_demographics_90K_PCA.RData", version =2)


load("dat_RECUR_ICD109_demographics_85K_ratio2.RData")
head(dat_RECUR_ICD109_demographics_85K) #4914
PCA <- read.table("GHS_Freeze_175.GL.splitmulti_FAM_175K_replication_ratio2_EXOME2_pca.eigenvec", header = F, stringsAsFactors = F)
dim(PCA) #4914*12
head(PCA)
PCA$cohort <- "replication"
PCA <- cSplit(PCA, "V1", "_")
colnames(PCA)[which(names(PCA) == "V1_2")] <- "PT_ID"
dat_RECUR_ICD109_demographics_85K_PCA <- merge(dat_RECUR_ICD109_demographics_85K, PCA, by = "PT_ID")
FAM_85K_replication_ratio2 <- data.frame(PCA$V2, PCA$V2) 
FAM_90K_85K_ratio2 <- rbind(FAM_90K_discovery_ratio2, FAM_85K_replication_ratio2)
colnames(FAM_90K_85K_ratio2) <- c("FID", "IID")
save(dat_RECUR_ICD109_demographics_90K_PCA, dat_RECUR_ICD109_demographics_85K_PCA, file = "dat_RECUR_ICD109_demographics_90K_85K_PCA.RData", version =2)
write.table(FAM_90K_85K_ratio2, file = "FAM_90K_85K_ratio2.txt", col.names = T, row.names = F, quote = F)
nrow(FAM_90K_85K_ratio2)  #12726

PCA <- read.table("GHS_Freeze_175.GL.splitmulti_FAM_175K_discoveryreplication_ratio2_EXOME2_pca.eigenvec", header = F, stringsAsFactors = F)
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
library(ggplot2)
# Basic scatter plot
ggplot(PCA, aes(x=V3, y=V5)) + geom_point()
# Change the point size, and shape
ggplot(mtcars, aes(x=wt, y=mpg)) +
  geom_point(size=2, shape=23)

# Change point shapes, colors and sizes
tiff(paste('PCA_discovery_replication_06062022', '.tiff', sep=''), units="in", width=14, height=10, res=300)
ggplot(PCA_final, aes(x=V3, y=V4, color=cohort, size=Index_age, alpha = 0.8)) +
  geom_point() + 
  theme(text = element_text(size = 12)) + 
  ylab("PC2") +
  xlab("PC1")
dev.off()

PCA_final$PT_RACE[PCA_final$PT_RACE == "Unknown"] <- "RACE_Unknown"
tiff(paste('PCA_discovery_replication_byRACE_06062022', '.tiff', sep=''), units="in", width=14, height=10, res=300)
ggplot(PCA_final, aes(x=V3, y=V4, shape=as.factor(PT_RACE), color=cohort, size=Index_age, alpha = 0.8)) +
  geom_point() + 
  theme(text = element_text(size = 12)) + 
  ylab("PC2") +
  xlab("PC1")
dev.off()
save(PCA_final, file = "060622022_PCA_final.RData", version =2)
load("060622022_PCA_final.RData")

WHITE <- PCA_final[(PCA_final$V3 < 0.005 & PCA_final$V4 > -0.025), ]$PT_ID
length(WHITE) # 11431
BLACK <- PCA_final[(PCA_final$V3 > 0.01 & PCA_final$V4 > -0.025), ]$PT_ID
length(BLACK) #972
#this is genetics determined race
save(WHITE, BLACK, file = "Genetics_determined_race.RData", version = 2)
load("Genetics_determined_race.RData")
########################################################################################################

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

#remove some cases(optional)
# Remove_cases <- gsub("PT", "", NOTCH3)
# phenotypes_long <- phenotypes_long[!(phenotypes_long$id %in% Remove_cases), ]

Summary_phenotypes <- phenotypes_long %>%
  dplyr::group_by(group, phecode) %>%
  dplyr::summarise(count=sum(phenotype)) %>%
  dplyr::arrange(desc(count)) %>%#3703
  ungroup()

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
  distinct(phecode, .keep_all = T)  #1755*4
str(phecode_mapping_unique)
#optional
phecode_mapping_unique$phecode <- paste0("X", phecode_mapping_unique$phecode)
# phecode_mapping_unique$phecode <- gsub("X", "", phecode_mapping_unique$phecode)

Summary_phenotypes_merge_anno <- merge(Summary_phenotypes_merge, phecode_mapping_unique, by = "phecode", all.x = T) #1853

#for a specific ratio
Summary_phenotypes_merge_anno$count.y <- Summary_phenotypes_merge_anno$count.y/2
#4899 change to 4889for discovery white; 2975 for replication white: 258 for discovery black; 222 for replication black; 
N = '<total number of >'
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


#load("Summary_phenotypes109_merge_anno_common_rare_delta_p123_90K_phenotypes109_black_06062022.RData")
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
  write.csv(fisherexact, "<path_to_directory>/top_fisherexact_phenotypes109_ratio2_90K_black.csv")
  saveRDS(fisherexact, file= "<path_to_directory>/top_fisherexact_phenotypes109_ratio2_90K_black.rds")
}

fisherexact <- data.frame(readRDS("<path_to_directory>/top_fisherexact_phenotypes109_ratio2_90K_black_06062022.rds"))
dim(fisherexact)
fisherexact$phecode <- rownames(fisherexact)
head(fisherexact)
Summary_phenotypes_merge_anno_delta_top_fisherexact <- merge(Summary_phenotypes_merge_anno_delta, fisherexact, by = "phecode")

#for matched
save(Summary_phenotypes_merge_anno_delta_top_fisherexact, file = "<path_to_directory>/Summary_phenotypes_merge_anno_delta_top_fisherexact_phenotypes109_ratio2_90K_black.RData", version = 2)

#create mat file
#select this function to get all cells filled 
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

#phenotypes109
table(phenotypes_ratio2$group)
phenotypes_ratio2_top <- phenotypes_ratio2[phenotypes_ratio2$group == 1, colnames(phenotypes_ratio2) %in% Summary_phenotypes_merge_anno_delta_top_fisherexact$phecode]
dim(phenotypes_ratio2_top)
colnames(phenotypes_ratio2_top)
phenotypes_ratio2_top <- phenotypes_ratio2_top %>%
  mutate_if(is.logical, as.factor)
mat <- fishermatrix(phenotypes_ratio2_top)
diag(mat) <- 1
mat_transformed <- round(-log10(mat), 0)

#create correlation links between pheocodes
CORR <- data.frame(row=rownames(mat_transformed)[row(mat_transformed)], col=colnames(mat_transformed)[col(mat_transformed)], corr=c(mat_transformed))
#for circo.plot
CORR <- CORR %>%
  filter(corr >3) #p.value < 0.001 -log transform = 3

#for cytoscape may not use
# CORR <- CORR %>%
#   filter(corr >1) #remove any p.value > 0.1 1;
save(mat_transformed, CORR, file = "<path_to_directory>/Summary_phenotypes_merge_anno_delta_top_fisherexact_phenotypes109_black_ratio2_90K_mat_corr.RData", version = 2)

1/(50*49) #0.0004081633

colnames(data1_remove)
colnames(CORR)[1] <- "phecode"
CORR_merge <- merge(CORR, data1_remove[, c("phecode", "exclude_name", "x")], by = "phecode")
colnames(CORR_merge)[1] <- "row"
colnames(CORR_merge)[2] <- "phecode"
CORR_merge <- merge(CORR_merge, data1_remove[, c("phecode", "exclude_name", "x")], by = "phecode")
colnames(CORR_merge)[1] <- "col"
head(CORR_merge)

############################################################################################################
