library(data.table)
library(dplyr)
library(xtable)
library(stargazer)

### GIAB biallelic SNVs #####

# read golden

golden = fread("./GIAB/Golden_Giab_ref_alt_SNVs.gz")
golden$V1[golden$V1=="X"] = 23
golden = golden %>% filter(V1 %in% 1:23)

golden$ID = do.call(paste0,list(golden$V1,"_",golden$V2,"_",golden$V3,"_",golden$V4))
golden$ID2 = do.call(paste0,list(golden$V1,"_",golden$V2))

colnames(golden)[5] = "GIAB"

### Deepvariant #####

# read all output from the caller

caller_all = fread("./GIAB/deepvariant_snvs_bialelic_final.gz",sep="\t")
caller_all$ID = do.call(paste0,list(caller_all$V1,"_",caller_all$V2,"_",caller_all$V3,"_",caller_all$V4))
caller_all$V1[caller_all$V1=="X"] = 23
caller_all = caller_all %>% filter(V1 %in% 1:23)

# read conservative zones

caller_cons = fread("./GIAB/Snvs_deepvariant_conservative.gz")
caller_cons$ID = do.call(paste0,list(caller_cons$V1,"_",caller_cons$V2,"_",caller_cons$V4,"_",caller_cons$V5))
caller_cons$V1[caller_cons$V1=="X"] = 23

# intersect conservatives zones

caller = caller_all[caller_all$ID %in% caller_cons$ID,]

# build pass variable

caller$PASS = "NO"

caller$PASS[caller$ID %in% golden$ID] = "YES"

# chr pos REF ALT GT GQ

colnames(caller)[3:5] = c("REF_deepvariant","ALT_deepvariant","GT_deepvariant")

deepvariant = caller %>% dplyr::select(c(7,6,1:5))

saveRDS(deepvariant,file = "./GIAB/outputs/deepvariant.rds")


### Gatk #####

# read all output from the caller

caller_all = fread("./GIAB/gatk_bialelic_SNVs_pass_final.gz",sep="\t")
caller_all$ID = do.call(paste0,list(caller_all$V1,"_",caller_all$V2,"_",caller_all$V3,"_",caller_all$V4))
caller_all$V1[caller_all$V1=="X"] = 23
caller_all = caller_all %>% filter(V1 %in% 1:23)


# read conservative zones

caller_cons = fread("./GIAB/Snvs_gatk_conservative.gz")
caller_cons$ID = do.call(paste0,list(caller_cons$V1,"_",caller_cons$V2,"_",caller_cons$V4,"_",caller_cons$V5))
caller_cons$V1[caller_cons$V1=="X"] = 23

# intersect conservatives zones

caller = caller_all[caller_all$ID %in% caller_cons$ID,]

# build pass variable

caller$PASS = "NO"

caller$PASS[caller$ID %in% golden$ID] = "YES"

colnames(caller)[c(3:5)] = c("REF_gatk","ALT_gatk",
                                 "GT_gatk")

gatk = caller %>% dplyr::select(c(7,6,1:5))

saveRDS(gatk,file = "./GIAB/outputs/gatk.rds")


### Strelka #####

# read all output from the caller

caller_all = fread("./GIAB/strelka_snvs_bialellic_final.gz",sep="\t")
caller_all$ID = do.call(paste0,list(caller_all$V1,"_",caller_all$V2,"_",caller_all$V3,"_",caller_all$V4))
caller_all$V1[caller_all$V1=="X"] = 23
caller_all = caller_all %>% filter(V1 %in% 1:23)


# read conservative zones

caller_cons = fread("./GIAB/Snvs_strelka_conservative.gz")
caller_cons$ID = do.call(paste0,list(caller_cons$V1,"_",caller_cons$V2,"_",caller_cons$V4,"_",caller_cons$V5))
caller_cons$V1[caller_cons$V1=="X"] = 23

# intersect conservatives zones

caller = caller_all[caller_all$ID %in% caller_cons$ID,]


# build pass variable

caller$PASS = "NO"

caller$PASS[caller$ID %in% golden$ID] = "YES"

# chr pos REF ALT GT GQ

colnames(caller)[c(3:5)] = c("REF_strelka","ALT_strelka",
                                "GT_strelka")

strelka = caller %>% dplyr::select(c(7,6,1:5))

saveRDS(strelka,file = "./GIAB/outputs/strelka.rds")


### Merge callers #####

deepvariant = readRDS("./GIAB/outputs/deepvariant.rds")
gatk = readRDS("./GIAB/outputs/gatk.rds")
strelka = readRDS("./GIAB/outputs/strelka.rds")
my_data = full_join(deepvariant,gatk)
my_data = full_join(my_data,strelka)


my_data$CHR = tstrsplit(my_data$ID,split="_")[[1]]
my_data$BP = tstrsplit(my_data$ID,split="_")[[2]]


saveRDS(my_data,file = "./GIAB/outputs/callers_data.rds")


### prepare model for all_genome ####

my_data = readRDS("./GIAB/outputs/callers_data.rds")

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
saveRDS(my_model,"./glm_model.rds")

### add genotype golden ####
my_data2 = my_data %>% dplyr::select("PASS","ID","REF_deepvariant","ALT_deepvariant",
                                     "REF_gatk","ALT_gatk",
                                     "REF_strelka","ALT_strelka",
                                     "CHR","BP","prediction")
my_data3 = my_data1 %>% dplyr::select("GT_deepvariant","GT_gatk","GT_strelka")

all_my_data <- cbind(my_data2,my_data3)
all_my_data1<- all_my_data[, c(1, 2, 3, 4, 12, 5, 6,13,7,8,14,9,10,11)]

giab  = left_join(all_my_data1,golden %>% dplyr::select(ID,GIAB))

giab = giab %>% dplyr::select(PASS,ID,GT_deepvariant,GT_gatk,GT_strelka,
                              GIAB,prediction)

saveRDS(giab,"./GIAB/outputs/GIAB_SNVs_biallelic.rds")

### genotype concordance ####

giab = readRDS("./GIAB/outputs/GIAB_SNVs_biallelic.rds")

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
      file = "./GIAB/outputs/GIAB_SNVs_biallelic_consensus.rds")


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

fwrite(final_table, "./table_snvs_model_giab_Lrm.csv", sep="\t")

xtable(my_tab_metrics,digits = 3)

stargazer(my_model)

### INSILICO 3 biallelic SNVs prediction #####

setwd("./insilico3/")

callers = c("GT_deepvariant","GT_gatk","GT_strelka")

caller = fread("Deepvariant_snvs_normalizados_Insilico3.gz")
colnames(caller)[5] = callers[1]
caller$ID = do.call(paste0,list(caller$V1,"_",caller$V2,"_",caller$V3,"_",caller$V4))
insilico  = caller %>% dplyr::select(ID,GT_deepvariant)

caller = fread("Gatk_snvs_normalizados_Insilico3.gz")
colnames(caller)[5] = callers[2]
caller=subset(caller, caller$GT_gatk !="2/1")
caller$ID = do.call(paste0,list(caller$V1,"_",caller$V2,"_",caller$V3,"_",caller$V4))
aux=caller %>% dplyr::select(ID,GT_gatk)
insilico = full_join(insilico,aux)

caller = fread("strelka_snvs_normalizados_Insilico3.gz")
colnames(caller)[5] = callers[3]
caller$ID = do.call(paste0,list(caller$V1,"_",caller$V2,"_",caller$V3,"_",caller$V4))
aux=caller %>% dplyr::select(ID,GT_strelka)
insilico = full_join(insilico,aux)

insilico$GT_deepvariant[is.na(insilico$GT_deepvariant)] = "9/9"
insilico$GT_gatk[is.na(insilico$GT_gatk)] = "9/9"
insilico$GT_strelka[is.na(insilico$GT_strelka)] = "9/9"
insilico$GT_strelka[!(insilico$GT_strelka %in% c("0/0","0/1","1/1"))] = "9/9"
insilico$GT_deepvariant[insilico$GT_deepvariant=="1/0"] = "0/1"

golden = fread("SNVs_Insilico_3.gz")
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
