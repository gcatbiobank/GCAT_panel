library(data.table)
library(dplyr)
library(xtable)
library(stargazer)

### GIAB biallelic SNVs #####

# read golden

golden = fread("./SNV/GIAB/Golden_Giab_ref_alt_SNVs")
golden$V1[golden$V1=="X"] = 23
golden = golden %>% filter(V1 %in% 1:23)

golden$ID = do.call(paste0,list(golden$V1,"_",golden$V2,"_",golden$V3,"_",golden$V4))
golden$ID2 = do.call(paste0,list(golden$V1,"_",golden$V2))

colnames(golden)[5] = "GIAB"

### Deepvariant #####

# read all output from the caller

caller_all = fread("./SNV/GIAB/deepvariant_snvs_bialelic_final",sep="\t")
caller_all$ID = do.call(paste0,list(caller_all$V1,"_",caller_all$V2,"_",caller_all$V3,"_",caller_all$V4))
caller_all$V1[caller_all$V1=="X"] = 23
caller_all = caller_all %>% filter(V1 %in% 1:23)

# read conservative zones

caller_cons = fread("./SNV/GIAB/Snvs_deepvariant_conservative")
caller_cons$ID = do.call(paste0,list(caller_cons$V1,"_",caller_cons$V2,"_",caller_cons$V4,"_",caller_cons$V5))
caller_cons$V1[caller_cons$V1=="X"] = 23

# intersect conservatives zones

caller = caller_all[caller_all$ID %in% caller_cons$ID,]

# build pass variable

caller$PASS = "NO"

caller$PASS[caller$ID %in% golden$ID] = "YES"

# name all variables

# PL

aux = tstrsplit(caller$V9,split=",")

caller$PL_AA_deepvariant = aux[[1]]
caller$PL_AB_deepvariant = aux[[2]]
caller$PL_BB_deepvariant = aux[[3]]

# AD

aux = tstrsplit(caller$V8,split=",")

caller$AD_REF_deepvariant = aux[[1]]
caller$AD_ALT_deepvariant = aux[[2]]

# chr pos REF ALT GT GQ

colnames(caller)[3:7] = c("REF_deepvariant","ALT_deepvariant","GT_deepvariant","GQ_deepvariant","DP_deepvariant")

deepvariant = caller %>% dplyr::select(c(11,10,3:7,15:16))

saveRDS(deepvariant,file = "./SNV/GIAB/outputs/deepvariant.rds")


### Gatk #####

# read all output from the caller

caller_all = fread("./SNV/GIAB/gatk_bialelic_SNVs_pass_final",sep="\t")
caller_all$ID = do.call(paste0,list(caller_all$V1,"_",caller_all$V2,"_",caller_all$V3,"_",caller_all$V4))
caller_all$V1[caller_all$V1=="X"] = 23
caller_all = caller_all %>% filter(V1 %in% 1:23)


# read conservative zones

caller_cons = fread("./SNV/GIAB/Snvs_gatk_conservative")
caller_cons$ID = do.call(paste0,list(caller_cons$V1,"_",caller_cons$V2,"_",caller_cons$V4,"_",caller_cons$V5))
caller_cons$V1[caller_cons$V1=="X"] = 23

# intersect conservatives zones

caller = caller_all[caller_all$ID %in% caller_cons$ID,]

# build pass variable

caller$PASS = "NO"

caller$PASS[caller$ID %in% golden$ID] = "YES"


# name all variables

# PL

aux = tstrsplit(caller$V9,split=",")

caller$PL_AA_gatk = aux[[1]]
caller$PL_AB_gatk = aux[[2]]
caller$PL_BB_gatk = aux[[3]]

# AD

aux = tstrsplit(caller$V6,split=",")

caller$AD_REF_gatk = aux[[1]]
caller$AD_ALT_gatk = aux[[2]]

# chr pos REF ALT GT GQ

colnames(caller)[c(3:5,7:8)] = c("REF_gatk","ALT_gatk",
                                 "GT_gatk","DP_gatk","GQ_gatk")

gatk = caller %>% dplyr::select(c(11,10,3:5,7:8,15:16))

saveRDS(gatk,file = "./SNV/GIAB/outputs/gatk.rds")


### Strelka #####

# read all output from the caller

caller_all = fread("./SNV/GIAB/strelka_snvs_bialellic_final",sep="\t")
caller_all$ID = do.call(paste0,list(caller_all$V1,"_",caller_all$V2,"_",caller_all$V3,"_",caller_all$V4))
caller_all$V1[caller_all$V1=="X"] = 23
caller_all = caller_all %>% filter(V1 %in% 1:23)


# read conservative zones

caller_cons = fread("./SNV/GIAB/Snvs_strelka_conservative")
caller_cons$ID = do.call(paste0,list(caller_cons$V1,"_",caller_cons$V2,"_",caller_cons$V4,"_",caller_cons$V5))
caller_cons$V1[caller_cons$V1=="X"] = 23

# intersect conservatives zones

caller = caller_all[caller_all$ID %in% caller_cons$ID,]


# build pass variable

caller$PASS = "NO"

caller$PASS[caller$ID %in% golden$ID] = "YES"

# name all variables

# AD

aux = tstrsplit(caller$V10,split=",")

caller$AD_REF_strelka = aux[[1]]
caller$AD_ALT_strelka = aux[[2]]

# ADF

aux = tstrsplit(caller$V11,split=",")

caller$ADF_REF_strelka = aux[[1]]
caller$ADF_ALT_strelka = aux[[2]]

# ADR

aux = tstrsplit(caller$V12,split=",")

caller$ADR_REF_strelka = aux[[1]]
caller$ADR_ALT_strelka = aux[[2]]

#PL 

aux = tstrsplit(caller$V15,split=",")

caller$PL_AA_strelka = as.numeric(aux[[1]])
caller$PL_AB_strelka = as.numeric(aux[[2]])
caller$PL_BB_strelka = as.numeric(aux[[3]])

# chr pos REF ALT GT GQ

colnames(caller)[c(3:9,13)] = c("REF_strelka","ALT_strelka",
                                "GT_strelka","GQ_strelka",
                                "GQx_strelka","DP_strelka",
                                "DPF_strelka","SB_strelka")

strelka = caller %>% dplyr::select(c(17:16,3:6,8,18:19))

saveRDS(strelka,file = "./SNV/GIAB/outputs/strelka.rds")


### Merge callers #####

deepvariant = readRDS("./SNV/GIAB/outputs/deepvariant.rds")
gatk = readRDS("./SNV/GIAB/outputs/gatk.rds")
strelka = readRDS("./SNV/GIAB/outputs/strelka.rds")
my_data = full_join(deepvariant,gatk)
my_data = full_join(my_data,strelka)


my_data$CHR = tstrsplit(my_data$ID,split="_")[[1]]
my_data$BP = tstrsplit(my_data$ID,split="_")[[2]]


saveRDS(my_data,file = "./SNV/GIAB/outputs/callers_data.rds")


### prepare model for all_genome ####

my_data = readRDS("./SNV/GIAB/outputs/callers_data.rds")

my_data$GT_deepvariant[my_data$GT_deepvariant=="0/0"] = "9/9"
my_data$GT_deepvariant[is.na(my_data$GT_deepvariant)] = "9/9"

my_data$GT_gatk[my_data$GT_gatk=="0/0"] = "9/9"
my_data$GT_gatk[is.na(my_data$GT_gatk)] = "9/9"


my_data$GT_strelka[is.na(my_data$GT_strelka)] = "9/9"
my_data$GT_strelka[my_data$GT_strelka==1] = "9/9"
my_data$GT_strelka[my_data$GT_strelka=="1/0"] = "0/1"
my_data$GT_strelka[my_data$GT_strelka=="0/0"] = "9/9"

my_data1=my_data
my_data$GT_deepvariant[my_data$GT_deepvariant=="0/1"] = "0/1-1/1"
my_data$GT_deepvariant[my_data$GT_deepvariant=="1/1"] = "0/1-1/1"

my_data$GT_gatk[my_data$GT_gatk=="0/1"] = "0/1-1/1"
my_data$GT_gatk[my_data$GT_gatk=="1/1"] = "0/1-1/1"

my_data$GT_strelka[my_data$GT_strelka=="0/1"] = "0/1-1/1"
my_data$GT_strelka[my_data$GT_strelka=="1/1"] = "0/1-1/1"

my_data$PASS = ifelse(my_data$PASS=="YES",0,1)

my_model = glm(PASS ~ GT_deepvariant + 
                 GT_gatk + GT_strelka,
               data=my_data,family = "binomial")

my_data$PASS = ifelse(my_data$PASS==0,"YES","NO")

summary(my_model)

my_pred = predict(my_model,my_data,type="response")

my_pred = ifelse(my_pred>0.5,1,0)

my_pred = ifelse(my_pred==1,"NO PASS","PASS")

my_data$prediction = my_pred
saveRDS(my_model,"./SNV/glm_model.rds")

### add genotype golden ####
my_data2 = my_data %>% dplyr::select("PASS","ID","REF_deepvariant","ALT_deepvariant","GQ_deepvariant","DP_deepvariant",
                                     "AD_REF_deepvariant","AD_ALT_deepvariant","REF_gatk","ALT_gatk",
                                     "DP_gatk","GQ_gatk","AD_REF_gatk","AD_ALT_gatk","REF_strelka","ALT_strelka",
                                     "GQ_strelka","DP_strelka","AD_REF_strelka","AD_ALT_strelka","CHR","BP","prediction")
my_data3 = my_data1 %>% dplyr::select("GT_deepvariant","GT_gatk","GT_strelka")

all_my_data <- cbind(my_data2,my_data3)
all_my_data1<- all_my_data[, c(1, 2, 3, 4, 24, 5, 6,7,8,9,10,25,11,12,13,14,15,16,26,17,18,19,20,21,22,23)]
giab  = left_join(all_my_data1,golden %>% dplyr::select(ID,GIAB))

giab = giab %>% dplyr::select(PASS,ID,GT_deepvariant,GT_gatk,GT_strelka,
                              GIAB,prediction)

saveRDS(giab,"./SNV/GIAB/outputs/GIAB_SNVs_biallelic.rds")

### genotype concordance ####

giab = readRDS("./SNV/GIAB/outputs/GIAB_SNVs_biallelic.rds")

giab$consensus = NA

giab$GT_deepvariant[giab$GT_deepvariant=="9/9"] = NA
giab$GT_gatk[giab$GT_gatk=="9/9"] = NA
giab$GT_strelka[giab$GT_strelka=="9/9"] = NA


giab_prediction_pass = giab %>% filter(prediction=="PASS")
giab_prediction_pass = as.matrix(giab_prediction_pass)
a <- Sys.time()
for(i in 1:nrow(giab_prediction_pass)){
  
  my_tab = table(giab_prediction_pass[i,3:5])
  
  if(length(which(my_tab==max(my_tab)))==1){
    giab_prediction_pass[i,8] = names(which.max(my_tab))
  }
  
  if(length(which(my_tab==max(my_tab)))>1){
    
    giab_prediction_pass[i,8] = NA
  }
}

b <- Sys.time()

b-a

giab_prediction_pass = as.data.frame(giab_prediction_pass)

giab_prediction_pass_golden = giab_prediction_pass %>% filter(PASS=="YES") %>% 
  filter(!is.na(consensus)) 

giab_prediction_pass_golden$GIAB[giab_prediction_pass_golden$GIAB=="0|1"] = "0/1"
giab_prediction_pass_golden$GIAB[giab_prediction_pass_golden$GIAB=="1|0"] = "0/1"
giab_prediction_pass_golden$GIAB[giab_prediction_pass_golden$GIAB=="1|1"] = "1/1"
giab_prediction_pass_golden$GIAB[giab_prediction_pass_golden$GIAB=="1/0"] = "0/1"
saveRDS(giab_prediction_pass_golden,
      file = "./SNV/GIAB/outputs/GIAB_SNVs_biallelic_consensus.rds")


# concordance genotype of the model with golden ######

giab_prediction_pass_golden$GIAB = as.character(giab_prediction_pass_golden$GIAB)
giab_prediction_pass_golden$consensus = as.character(giab_prediction_pass_golden$consensus)
model_g_error= 100 -(sum(giab_prediction_pass_golden$consensus==giab_prediction_pass_golden$GIAB)/
                       nrow(giab_prediction_pass_golden))*100

# concordance genotype with the rest of callers #####

deepvariant = giab_prediction_pass_golden %>% filter(!is.na(GT_deepvariant)) %>% filter(PASS=="YES")

deep_g_error = 100 - (sum(deepvariant$GIAB==deepvariant$GT_deepvariant)/
                        nrow(deepvariant))*100

gatk = giab_prediction_pass_golden %>% filter(!is.na(GT_gatk)) %>% filter(PASS=="YES")

gatk_g_error = 100 -(sum(gatk$GIAB==gatk$GT_gatk,na.rm = T)/
                       nrow(gatk))*100

strelka = giab_prediction_pass_golden %>% filter(!is.na(GT_strelka)) %>% filter(PASS=="YES")

strelka_g_error = 100 - (sum(strelka$GIAB==strelka$GT_strelka)/
                           nrow(strelka))*100


# non concordance ########

non_concordance = which(giab_prediction_pass_golden$consensus!=giab_prediction_pass_golden$GIAB)

non_concordance = giab_prediction_pass_golden[non_concordance,]

### sensitivity specificity for each caller #####

deepvariant = giab %>% filter(!is.na(GT_deepvariant))

deep_sens = (sum(deepvariant$ID %in% unique(golden$ID))/length(unique(golden$ID)))*100

deep_spec = (sum(deepvariant$ID %in% unique(golden$ID))/length(unique(deepvariant$ID)))*100


gatk = giab %>% filter(!is.na(GT_gatk))

gatk_sens = (sum(gatk$ID %in% unique(golden$ID))/length(unique(golden$ID)))*100

gatk_spec = (sum(gatk$ID %in% unique(golden$ID))/length(unique(gatk$ID)))*100

strelka = giab %>% filter(!is.na(GT_strelka))

strelka_sens = (sum(strelka$ID %in% unique(golden$ID))/length(unique(golden$ID)))*100

strelka_spec = (sum(strelka$ID %in% unique(golden$ID))/length(unique(strelka$ID)))*100

## table sens,spec, geno error #####

my_tab_metrics = matrix(c(deep_sens,gatk_sens,strelka_sens,model_sens,
                          deep_spec,gatk_spec,strelka_spec,model_spec,
                          deep_g_error,gatk_g_error,strelka_g_error,model_g_error),ncol = 3)

final_table = c("Deepvariant", "Haplotype caller", "Strelka2","Model")
final_table= as.data.frame(final_table)
colnames(final_table)[1] <- "Variant callers"
final_table$Recall= c(deep_sens,gatk_sens,strelka_sens,model_sens)
final_table$Precision = c(deep_spec,gatk_spec,strelka_spec,model_spec)
final_table$`Genotype error`= c(deep_g_error,gatk_g_error,strelka_g_error,model_g_error)

fwrite(final_table, "./SNV/table_snvs_model_giab_Lrm.csv", sep="\t")

xtable(my_tab_metrics,digits = 3)

stargazer(my_model)

### INSILICO 3 biallelic SNVs prediction #####

setwd("./SNV/insilico3/")

my_files = list.files(".")

callers = c("GT_deepvariant","GT_gatk","GT_strelka")

caller = fread(my_files[1])
colnames(caller)[5] = callers[1]
caller$ID = do.call(paste0,list(caller$V1,"_",caller$V2,"_",caller$V3,"_",caller$V4))
insilico  = caller %>% dplyr::select(ID,GT_deepvariant)

caller = fread(my_files[2])
colnames(caller)[5] = callers[2]
caller=subset(caller, caller$GT_gatk !="2/1")
caller$ID = do.call(paste0,list(caller$V1,"_",caller$V2,"_",caller$V3,"_",caller$V4))
aux=caller %>% dplyr::select(ID,GT_gatk)
insilico = full_join(insilico,aux)

caller = fread(my_files[5])
colnames(caller)[5] = callers[3]
caller$ID = do.call(paste0,list(caller$V1,"_",caller$V2,"_",caller$V3,"_",caller$V4))
aux=caller %>% dplyr::select(ID,GT_strelka)
insilico = full_join(insilico,aux)

insilico$GT_deepvariant[is.na(insilico$GT_deepvariant)] = "9/9"
insilico$GT_gatk[is.na(insilico$GT_gatk)] = "9/9"
insilico$GT_strelka[is.na(insilico$GT_strelka)] = "9/9"
insilico$GT_strelka[!(insilico$GT_strelka %in% c("0/0","0/1","1/1"))] = "9/9"
insilico$GT_deepvariant[insilico$GT_deepvariant=="1/0"] = "0/1"

golden = fread("SNVs_Insilico_3")
golden$V1 = gsub("chr","",golden$V1)
golden$ID = do.call(paste0,list(golden$V1,"_",golden$V3,"_",golden$V4,"_",golden$V5))
insilico$PASS = "NO"
insilico[insilico$ID %in% golden$ID,]$PASS = "YES"

### add genotype golden #####

golden$V1[golden$V1=="X"] = 23

golden$V3 = as.numeric(golden$V3)

golden_geno = NULL

for(i in 1:23){
  
  aux = golden[golden$V1==i,]
  
  which_duplicated = aux$V3[which(duplicated(aux$V3))]
  
  golden_dup = unique(aux[aux$V3 %in% which_duplicated,])
  golden_no_dup = aux[!aux$V3 %in% which_duplicated,]
  
  golden_dup$genotype = "1/1"
  golden_no_dup$genotype = "0/1"
  
  aux_geno = rbind(golden_dup,golden_no_dup)
  
  golden_geno = rbind(golden_geno,aux_geno)
}

insilico = left_join(insilico,golden_geno %>% dplyr::select(ID,genotype))

## prediction ######

insilico$CHR = tstrsplit(insilico$ID,split="_")[[1]]
insilico$CHR[insilico$CHR=="X"] = 23

insilico = insilico %>% filter(CHR %in% 1:23)
insilico1 = insilico
insilico$GT_deepvariant[insilico$GT_deepvariant=="0/1"] = "0/1-1/1"
insilico$GT_deepvariant[insilico$GT_deepvariant=="1/1"] = "0/1-1/1"

insilico$GT_gatk[insilico$GT_gatk=="0/1"] = "0/1-1/1"
insilico$GT_gatk[insilico$GT_gatk=="1/1"] = "0/1-1/1"

insilico$GT_strelka[insilico$GT_strelka=="0/1"] = "0/1-1/1"
insilico$GT_strelka[insilico$GT_strelka=="1/1"] = "0/1-1/1"

summary(my_model)

my_pred = predict(my_model,insilico,type="response")

my_pred = ifelse(my_pred>0.5,1,0)

my_pred = ifelse(my_pred==0,"PASS","NO PASS")

insilico$prediction = my_pred

insilico = insilico[,c(5,1,2:4,6,8)]
insilico =insilico %>% rename(chr_pos=ID)
colnames(insilico)[6] = "insilico3"
insilico2 = insilico %>% dplyr::select("PASS","chr_pos","insilico3","prediction")
insilico3 = insilico1 %>% dplyr::select("GT_deepvariant","GT_gatk","GT_strelka")
all_insilico <- cbind(insilico2,insilico3)
all_insilico2<- all_insilico[, c(1, 2, 5, 6, 7, 3, 4)]
saveRDS(all_insilico2,file = "./outputs/insilico3_SNVs_biallelic.rds")

### genotype concordance ####

insilico = readRDS("./outputs/insilico3_SNVs_biallelic.rds")

insilico$consensus = NA

insilico$GT_deepvariant[insilico$GT_deepvariant=="9/9"] = NA
insilico$GT_gatk[insilico$GT_gatk=="9/9"] = NA
insilico$GT_strelka[insilico$GT_strelka=="9/9"] = NA

insilico_prediction_pass = insilico %>% filter(prediction=="PASS")
insilico_prediction_pass = as.matrix(insilico_prediction_pass)

a <- Sys.time()
for(i in 1:nrow(insilico_prediction_pass)){
  
  my_tab = table(insilico_prediction_pass[i,3:5])
  
  if(length(which(my_tab==max(my_tab)))==1){
    insilico_prediction_pass[i,8] = names(which.max(my_tab))
  }
  
  if(length(which(my_tab==max(my_tab)))>1){
    
    insilico_prediction_pass[i,8] = NA
  }
}

b <- Sys.time()

b-a

insilico_prediction_pass = as.data.frame(insilico_prediction_pass)

insilico_prediction_pass_golden = insilico_prediction_pass %>% filter(PASS=="YES") %>%
  filter(prediction=="PASS") %>% filter(!is.na(consensus)) 

saveRDS(insilico_prediction_pass_golden,
        file = "./outputs/insilico3_SNVs_biallelic_consensus.rds")


# concordance genotype of the model with golden ######
model_g_error= 100 -(sum(insilico_prediction_pass_golden$consensus==insilico_prediction_pass_golden$insilico3)/
                       nrow(insilico_prediction_pass_golden))*100

# concordance genotype with the rest of callers #####

deepvariant = insilico %>% filter(!is.na(GT_deepvariant)) %>% filter(PASS=="YES")

deep_g_error = 100 - (sum(deepvariant$insilico3==deepvariant$GT_deepvariant)/
                        nrow(deepvariant))*100

gatk = insilico %>% filter(!is.na(GT_gatk)) %>% filter(PASS=="YES")

gatk_g_error = 100 -(sum(gatk$insilico3==gatk$GT_gatk,na.rm = T)/
                       nrow(gatk))*100

strelka = insilico %>% filter(!is.na(GT_strelka)) %>% filter(PASS=="YES")

strelka_g_error = 100 - (sum(strelka$insilico3==strelka$GT_strelka)/nrow(strelka))*100

# non concordance ########

non_concordance = which(insilico_prediction_pass_golden$consensus!=insilico_prediction_pass_golden$insilico3)

non_concordance = insilico_prediction_pass_golden[non_concordance,]

### sensitivity specificity for each caller in ref, alt, chrom, and position #####
deepvariant = insilico %>% filter(!is.na(GT_deepvariant))

deep_sens = (sum(deepvariant$chr_pos %in% unique(golden$ID))/length(unique(golden$ID)))*100

deep_spec = (sum(deepvariant$chr_pos %in% unique(golden$ID))/length(unique(deepvariant$chr_pos)))*100
gatk = insilico %>% filter(!is.na(GT_gatk))

gatk_sens = (sum(gatk$chr_pos %in% unique(golden$ID))/length(unique(golden$ID)))*100

gatk_spec = (sum(gatk$chr_pos %in% unique(golden$ID))/length(unique(gatk$chr_pos)))*100

strelka = insilico %>% filter(!is.na(GT_strelka))

strelka_sens = (sum(strelka$chr_pos %in% unique(golden$ID))/length(unique(golden$ID)))*100

strelka_spec = (sum(strelka$chr_pos %in% unique(golden$ID))/length(unique(strelka$chr_pos)))*100

## table sens,spec, geno error #####

my_tab_metrics = matrix(c(deep_sens,gatk_sens,strelka_sens,model_sens,
                          deep_spec,gatk_spec,strelka_spec,model_spec,
                          deep_g_error,gatk_g_error,strelka_g_error,model_g_error),ncol = 3)

final_table1 = c("Deepvariant", "Haplotype caller", "Strelka2","Model")
final_table1= as.data.frame(final_table1)
colnames(final_table1)[1] <- "Variant callers"
final_table1$Recall= c(deep_sens,gatk_sens,strelka_sens,model_sens)
final_table1$Precision = c(deep_spec,gatk_spec,strelka_spec,model_spec)
final_table1$`Genotype error`= c(deep_g_error,gatk_g_error,strelka_g_error,model_g_error)
setwd("../")
fwrite(final_table1, "./table_snvs_model_insilico3.csv", sep="\t")
