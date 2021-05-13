library(data.table)
library(dplyr)

setwd("C:/Users/ivangalvan/Desktop/SV_gcat_panel_github/")

source("ext/functions_JV_IGF_IC.R")

# read golden ######

golden_del = fread("mid_DEL/data/Insilico_mid_DEL")

golden_del = golden_del %>% filter(V3>30 & V3<151)

golden_del$V1 = gsub("chr","",golden_del$V1)

golden_geno_del = NULL

for(i in names(table(golden_del$V1))){
  
  aux = golden_del[golden_del$V1==i,]
  
  which_duplicated = aux$V2[which(duplicated(aux$V2))]
  
  golden_dup = unique(aux[aux$V2 %in% which_duplicated,])
  golden_no_dup = aux[!aux$V2 %in% which_duplicated,]
  
  if(nrow(golden_dup)>0 & nrow(golden_no_dup)>0 ){
    golden_dup$genotype = "1/1"
    golden_no_dup$genotype = "0/1"
    
    aux_geno = rbind(golden_dup,golden_no_dup)
  }
  
  if(nrow(golden_dup)==0){
    golden_no_dup$genotype = "0/1"
    
    aux_geno = golden_no_dup
  }
  
  if(nrow(golden_no_dup)==0){
    golden_dup$genotype = "0/1"
    
    aux_geno = golden_dup
  }
  
  golden_geno_del = rbind(golden_geno_del,aux_geno)
}

golden = golden_geno_del
colnames(golden) = c("chr","start_golden","length_golden","GT_golden")
golden$lower_golden = golden$start_golden-0
golden$upper_golden = golden$start_golden+0
golden$chr_pos_golden = do.call(paste0,list(golden$chr,"_",golden$start_golden))


# read vcfs files ########

delly = fread("mid_DEL/data/Delly_mid_DEL")
delly$V4 = as.numeric(abs(delly$V4))
delly = delly %>% filter(V4>30 & V4<151)
colnames(delly) = c("chr","start_delly","end","length_delly","GT_delly")
delly$lower_delly = delly$start_delly-10
delly$upper_delly = delly$start_delly+10
delly$chr_pos_delly = do.call(paste0,list(delly$chr,"_",delly$start_delly))


lumpy = fread("mid_DEL/data/Lumpy_mid_DEL")
lumpy$V4 = as.numeric(abs(lumpy$V4))
lumpy = lumpy %>% filter(V4>30 & V4<151)
colnames(lumpy) = c("chr","start_lumpy","end","length_lumpy","GT_lumpy")
lumpy$lower_lumpy = lumpy$start_lumpy-10
lumpy$upper_lumpy = lumpy$start_lumpy+10
lumpy$chr_pos_lumpy = do.call(paste0,list(lumpy$chr,"_",lumpy$start_lumpy))


pindel = fread("mid_DEL/data/Pindel_mid_DEL")
colnames(pindel) = c("chr","start_pindel","end","length_pindel","GT_pindel")
pindel = pindel %>% filter(length_pindel>30 & length_pindel<151)
pindel$lower_pindel = pindel$start_pindel-10
pindel$upper_pindel = pindel$start_pindel+10
pindel$chr_pos_pindel = do.call(paste0,list(pindel$chr,"_",pindel$start_pindel))

whamg = fread("mid_DEL/data/Whamg_mid_DEL",fill = T)
whamg$V4 = as.numeric(abs(whamg$V4))
whamg = whamg %>% filter(V4>30 & V4<151)
colnames(whamg) = c("chr","start_whamg","end","length_whamg","GT_whamg")
whamg$lower_whamg = whamg$start_whamg-10
whamg$upper_whamg = whamg$start_whamg+10
whamg$chr_pos_whamg = do.call(paste0,list(whamg$chr,"_",whamg$start_whamg))

svaba = fread("mid_DEL/data/Svaba_mid_DEL")
svaba$length_svaba = abs(nchar(svaba$V3)-nchar(svaba$V4))
svaba = svaba %>% filter(length_svaba>30 & length_svaba<151)
colnames(svaba)[c(1,2,5)] = c("chr","start_svaba","GT_svaba")
svaba = svaba[,c(1,2,5,6)]
svaba$lower_svaba = svaba$start_svaba-10
svaba$upper_svaba = svaba$start_svaba+10
svaba$chr_pos_svaba = do.call(paste0,list(svaba$chr,"_",svaba$start_svaba))

gatk = fread("mid_DEL/data/Gatk_mid_DEL")
gatk$length_gatk = abs(nchar(gatk$V3)-nchar(gatk$V4))
gatk = gatk %>% filter(length_gatk>30 & length_gatk<151)
colnames(gatk)[c(1,2,5)] = c("chr","start_gatk","GT_gatk")
gatk = gatk[,c(1,2,5,6)]
gatk$end = gatk$length_gatk+gatk$start_gatk
gatk$lower_gatk = gatk$start_gatk-10
gatk$upper_gatk = gatk$start_gatk+10
gatk$chr_pos_gatk = do.call(paste0,list(gatk$chr,"_",gatk$start_gatk))

deepvariant = fread("mid_DEL/data/Deepvariant_mid_DEL")
deepvariant$length_deepvariant = abs(nchar(deepvariant$V3)-nchar(deepvariant$V4))
deepvariant = deepvariant %>% filter(length_deepvariant>30 & length_deepvariant<151)
colnames(deepvariant)[c(1,2,5)] = c("chr","start_deepvariant","GT_deepvariant")
deepvariant = deepvariant[,c(1,2,5,6)]
deepvariant$end = deepvariant$length_deepvariant+deepvariant$start_deepvariant
deepvariant$lower_deepvariant = deepvariant$start_deepvariant-10
deepvariant$upper_deepvariant = deepvariant$start_deepvariant+10
deepvariant$chr_pos_deepvariant = do.call(paste0,list(deepvariant$chr,"_",deepvariant$start_deepvariant))

strelka = fread("mid_DEL/data/Strelka_mid_DEL")
strelka$length_strelka = abs(nchar(strelka$V3)-nchar(strelka$V4))
strelka = strelka %>% filter(length_strelka>30 & length_strelka<151)
colnames(strelka)[c(1,2,5)] = c("chr","start_strelka","GT_strelka")
strelka = strelka[,c(1,2,5,6)]
strelka$end = strelka$length_strelka+strelka$start_strelka
strelka$lower_strelka = strelka$start_strelka-10
strelka$upper_strelka = strelka$start_strelka+10
strelka$chr_pos_strelka = do.call(paste0,list(strelka$chr,"_",strelka$start_strelka))

manta = fread("mid_DEL/data/Manta_mid_DEL_50_150",fill = T)
manta$V4 = as.numeric(abs(manta$V4))
manta = manta %>% filter(V4>30 & V4<151)
colnames(manta) = c("chr","start_manta","end","length_manta","GT_manta")
manta$lower_manta = manta$start_manta-10
manta$upper_manta = manta$start_manta+10
manta$chr_pos_manta = do.call(paste0,list(manta$chr,"_",manta$start_manta))

# we add "0/1" gewnotype because Manta for deletions from 30 to 50 doesn't give genotype
manta2 = fread("mid_DEL/data/Manta_mid_DEL_30_49",fill = T)
manta2$V4 = as.numeric(abs(manta2$V4))
manta2 = manta2 %>% filter(V4>30 & V4<151)
manta2$GT_manta = "0/1"
colnames(manta2) = c("chr","start_manta","end","length_manta","GT_manta")
manta2$lower_manta = manta2$start_manta-10
manta2$upper_manta = manta2$start_manta+10
manta2$chr_pos_manta = do.call(paste0,list(manta2$chr,"_",manta2$start_manta))

manta = rbind(manta,manta2)

### sensitivity-specifitiy-geno error for each caller ######

# gatk 

sens_gatk = sensitivity_precision_geno_error_deletions_duplications(gatk,golden,"deletion","gatk")

# deepvariant

sens_deep = sensitivity_precision_geno_error_deletions_duplications(deepvariant,golden,"deletion","deepvariant")

# strelka

sens_strelka = sensitivity_precision_geno_error_deletions_duplications(strelka,golden,"deletion","strelka")

# lumpy

sens_lumpy = sensitivity_precision_geno_error_deletions_duplications(lumpy,golden,"deletion","lumpy")

# pindel

sens_pindel = sensitivity_precision_geno_error_deletions_duplications(pindel,golden,"deletion","pindel")

# manta

manta[manta$length_manta<51,]$GT_manta = "9/9" # dels 30-50bp no genotype

sens_manta = sensitivity_precision_geno_error_deletions_duplications(manta,golden,"deletion","manta")

sens_manta$n_het = sum(sens_manta$all$GT_manta=="0/1" & !is.na(sens_manta$all$GT_golden))

sens_manta$n_homo = sum(sens_manta$all$GT_manta=="1/1" & !is.na(sens_manta$all$GT_golden))

sens_manta$g_error_het = 100 - 100*sum(sens_manta$all[GT_manta=="0/1",]$GT_manta==sens_manta$all[GT_manta=="0/1",]$GT_golden,na.rm = T)/sens_manta$n_het 

sens_manta$g_error_hom = 100 - 100*sum(sens_manta$all[GT_manta=="1/1",]$GT_manta==sens_manta$all[GT_manta=="1/1",]$GT_golden,na.rm = T)/sens_manta$n_homo

# svaba

sens_svaba = sensitivity_precision_geno_error_deletions_duplications(svaba,golden,"deletion","svaba")

# whamg

sens_whamg = sensitivity_precision_geno_error_deletions_duplications(whamg,golden,"deletion","whamg")

# delly

sens_delly = sensitivity_precision_geno_error_deletions_duplications(delly,golden,"deletion","delly")


metrics_callers = matrix(c(t(unlist(sens_deep))[c(3,4,1,12,2,13,5,6,7,9,8,10)],
                           t(unlist(sens_delly))[c(3,4,1,12,2,13,5,6,7,9,8,10)],
                           t(unlist(sens_gatk))[c(3,4,1,12,2,13,5,6,7,9,8,10)],
                           t(unlist(sens_lumpy))[c(3,4,1,12,2,13,5,6,7,9,8,10)],
                           t(unlist(sens_manta))[c(3,4,1,12,2,13,5,6,7,9,8,10)],
                           t(unlist(sens_pindel))[c(3,4,1,12,2,13,5,6,7,9,8,10)],
                           t(unlist(sens_strelka))[c(3,4,1,12,2,13,5,6,7,9,8,10)],
                           t(unlist(sens_svaba))[c(3,4,1,12,2,13,5,6,7,9,8,10)],
                           t(unlist(sens_whamg))[c(3,4,1,12,2,13,5,6,7,9,8,10)]),byrow=T,ncol = 12)

metrics_callers1 = as.data.frame(metrics_callers)

c =apply(metrics_callers, 2, as.numeric)

c = cbind(c("Deepvariant","Delly","Gatk","Lumpy","Manta","Pindel","Strelka","Svaba","Whamg"),
          round(c,2))


colnames(c) = c("Caller","TP","FP","Recall","CI Recall","Precision","CI Precision","F1-Score","Genotype error",
                "N 0/1","Genotype error 0/1","N 1/1","Genotype error 1/1")

c = as.data.frame(c)

c$`CI Recall` = metrics_callers1$V4
c$`CI Precision` = metrics_callers1$V6

c = c %>% as.data.frame() %>% arrange(desc(`F1-Score`))

c = cbind(N=c(nrow(golden),rep("",nrow(c)-1)),
          c)

c = cbind(SV_type=c("Deletions 30-150bp",rep("",nrow(c)-1)),
          c)


c[,1:3] = apply(c[,1:3], 2, function(x) as.character(x));

metrics_callers = c

write.csv(metrics_callers,"mid_DEL/outputs/metrics_mid_DEL.csv",row.names = F)



## prepare BBDD ordered by F-score ####

delly$ID = paste0("deletion_",1:nrow(delly))
manta$ID = "none"
svaba$ID = "none"
lumpy$ID = "none"
whamg$ID = "none"
gatk$ID = "none"
pindel$ID = "none"
deepvariant$ID = "none"
strelka$ID = "none"
golden$ID = "none"


my_data = merge_callers_deletions_duplications(manta,delly,callers_ref = c("delly"),
                                               caller_to_merge = c("manta"),svtype="deletion",repro = 0.8) 

my_data = merge_callers_deletions_duplications(svaba,my_data,callers_ref = c("delly","manta"),
                                               caller_to_merge = c("svaba"),svtype="deletion",repro = 0.8) 

my_data = merge_callers_deletions_duplications(lumpy,my_data,callers_ref = c("delly","manta","svaba"),
                                               caller_to_merge = c("lumpy"),svtype="deletion",repro = 0.8) 

my_data = merge_callers_deletions_duplications(whamg,my_data,callers_ref = c("delly","manta","svaba","lumpy"),
                                               caller_to_merge = c("whamg"),svtype="deletion",repro = 0.8) 

my_data = merge_callers_deletions_duplications(gatk,my_data,callers_ref = c("delly","manta","svaba","lumpy",
                                                                             "whamg"),
                                               caller_to_merge = c("gatk"),svtype="deletion",repro = 0.8) 

my_data = merge_callers_deletions_duplications(pindel,my_data,callers_ref = c("delly","manta","svaba","lumpy",
                                                                            "whamg","gatk"),
                                               caller_to_merge = c("pindel"),svtype="deletion",repro = 0.8) 

my_data = merge_callers_deletions_duplications(deepvariant,my_data,callers_ref = c("delly","manta","svaba","lumpy",
                                                                              "whamg","gatk","pindel"),
                                               caller_to_merge = c("deepvariant"),svtype="deletion",repro = 0.8) 

my_data = merge_callers_deletions_duplications(strelka,my_data,callers_ref = c("delly","manta","svaba","lumpy",
                                                                                "whamg","gatk","pindel","deepvariant"),
                                               caller_to_merge = c("strelka"),svtype="deletion",repro = 0.8) 

my_data2 = merge_callers_deletions_duplications(golden,my_data,callers_ref = c("delly","manta","svaba","lumpy",
                                                                               "whamg","gatk","pindel","deepvariant","strelka"),
                                                caller_to_merge = c("golden"),svtype="deletion",repro = 0.8)

saveRDS(my_data2,"mid_DEL/outputs/merge_callers_mid_DEL.rds")


### consensus length #######

my_data2$length = consensus_length(my_data2,
                                   callers = c("gatk","whamg","svaba","lumpy",
                                               "pindel","delly","strelka","deepvariant",
                                               "manta"))

summary(my_data2$length)

aux = cut(my_data2$length,
          breaks = c(30,50,75,100,125,151),right = F)

my_data2$length_stretch = factor(aux)


## add number of callers detected #####

my_data2$callers_detected = n_callers_detected(my_data2,
                                   callers = c("gatk","whamg","svaba","lumpy",
                                               "pindel","delly","strelka","deepvariant",
                                               "manta"))

table(my_data2$callers_detected)

my_data2$callers_detected[my_data2$callers_detected=="6"] = "6-7-8-9"
my_data2$callers_detected[my_data2$callers_detected=="7"] = "6-7-8-9"
my_data2$callers_detected[my_data2$callers_detected=="8"] = "6-7-8-9"
my_data2$callers_detected[my_data2$callers_detected=="9"] = "6-7-8-9"

## add reciprocity ###

my_data2$reciprocity = reprocicity(my_data2,
                                          callers = c("gatk","whamg","svaba","lumpy",
                                                      "pindel","delly","strelka","deepvariant",
                                                      "manta"))
summary(my_data2$reciprocity)


## add strategy ####

strategies = fread("ext/strategies.csv")

my_data2$strategy = strategy(my_data2,
                                callers = c("gatk","whamg","svaba","lumpy",
                                            "pindel","delly","strelka","deepvariant",
                                            "manta"),
                             strategies)

table(my_data2$strategy)


### colineality #####

table(my_data2$callers_detected,my_data2$strategy)
chisq.test(table(my_data2$callers_detected,my_data2$strategy))


## prepare model choose the best model####

my_data2$PASS = "YES"

my_data2$PASS[is.na(my_data2$length_golden)] = "NO"
my_data2 = my_data2 %>% filter(callers_detected!=0)

table(my_data2$PASS)

dim(golden)

my_data2$GT_deepvariant[is.na(my_data2$GT_deepvariant)] = "9/9"
my_data2$GT_gatk[is.na(my_data2$GT_gatk)] = "9/9"
my_data2$GT_strelka[is.na(my_data2$GT_strelka)] = "9/9"
my_data2$GT_manta[is.na(my_data2$GT_manta)] = "9/9"
my_data2$GT_whamg[is.na(my_data2$GT_whamg)] = "9/9"
my_data2$GT_delly[is.na(my_data2$GT_delly)] = "9/9"
my_data2$GT_lumpy[is.na(my_data2$GT_lumpy)] = "9/9"
my_data2$GT_pindel[is.na(my_data2$GT_pindel)] = "9/9"
my_data2$GT_svaba[is.na(my_data2$GT_svaba)] = "9/9"

my_data2$GT_deepvariant[my_data2$GT_deepvariant=="0/0"] = "9/9"
my_data2$GT_whamg[my_data2$GT_whamg==""] = "9/9"

my_data2$GT_deepvariant2 = my_data2$GT_deepvariant
my_data2$GT_gatk2 = my_data2$GT_gatk
my_data2$GT_strelka2 = my_data2$GT_strelka
my_data2$GT_manta2 = my_data2$GT_manta
my_data2$GT_whamg2 = my_data2$GT_whamg
my_data2$GT_delly2 = my_data2$GT_delly
my_data2$GT_lumpy2 = my_data2$GT_lumpy
my_data2$GT_pindel2 = my_data2$GT_pindel
my_data2$GT_svaba2 = my_data2$GT_svaba

my_data2$GT_deepvariant[my_data2$GT_deepvariant=="0/1"] = "0/1-1/1"
my_data2$GT_deepvariant[my_data2$GT_deepvariant=="1/1"] = "0/1-1/1"
my_data2$GT_gatk[my_data2$GT_gatk=="0/1"] = "0/1-1/1"
my_data2$GT_gatk[my_data2$GT_gatk=="1/1"] = "0/1-1/1"
my_data2$GT_strelka[my_data2$GT_strelka=="0/1"] = "0/1-1/1"
my_data2$GT_strelka[my_data2$GT_strelka=="1/1"] = "0/1-1/1"
my_data2$GT_manta[my_data2$GT_manta=="0/1"] = "0/1-1/1"
my_data2$GT_manta[my_data2$GT_manta=="1/1"] = "0/1-1/1"
my_data2$GT_whamg[my_data2$GT_whamg=="0/1"] = "0/1-1/1"
my_data2$GT_whamg[my_data2$GT_whamg=="1/1"] = "0/1-1/1"
my_data2$GT_delly[my_data2$GT_delly=="0/1"] = "0/1-1/1"
my_data2$GT_delly[my_data2$GT_delly=="1/1"] = "0/1-1/1"
my_data2$GT_lumpy[my_data2$GT_lumpy=="0/1"] = "0/1-1/1"
my_data2$GT_lumpy[my_data2$GT_lumpy=="1/1"] = "0/1-1/1"
my_data2$GT_pindel[my_data2$GT_pindel=="0/1"] = "0/1-1/1"
my_data2$GT_pindel[my_data2$GT_pindel=="1/1"] = "0/1-1/1"
my_data2$GT_svaba[my_data2$GT_svaba=="0/1"] = "0/1-1/1"
my_data2$GT_svaba[my_data2$GT_svaba=="1/1"] = "0/1-1/1"

table(my_data2$GT_deepvariant,useNA = "always")
table(my_data2$GT_gatk,useNA = "always")
table(my_data2$GT_strelka,useNA = "always")
table(my_data2$GT_manta,useNA = "always")
table(my_data2$GT_whamg,useNA = "always")
table(my_data2$GT_delly,useNA = "always")
table(my_data2$GT_lumpy,useNA = "always")
table(my_data2$GT_svaba,useNA = "always")
table(my_data2$GT_pindel,useNA = "always")

my_data2$length = as.factor(my_data2$length)

my_data2$PASS = ifelse(my_data2$PASS=="YES",0,1)

my_model = glm(PASS ~ factor(GT_deepvariant) + 
                 factor(GT_gatk) + 
                 factor(GT_strelka)+
                 factor(GT_manta)+
                 factor(GT_whamg)+
                 factor(GT_delly)+
                 factor(GT_lumpy)+
                 factor(GT_pindel)+
                 factor(GT_svaba) + 
                 factor(callers_detected),
               data=my_data2,family = "binomial")

my_data2$PASS = ifelse(my_data2$PASS==0,"YES","NO")

summary(my_model)

my_pred = predict(my_model,my_data2,type="response")

my_pred = ifelse(my_pred>0.5,1,0)

my_pred = ifelse(my_pred==1,"NO PASS","PASS")

#sensitivity
model_sens = (table(my_pred,my_data2$PASS)[2,2]/nrow(golden))*100
#specificity
model_spec = (table(my_pred,my_data2$PASS)[2,2])/(table(my_pred,my_data2$PASS)[2,2]+
                                                    table(my_pred,my_data2$PASS)[2,1])*100

model_tp = as.numeric(table(my_pred,my_data2$PASS)[2,2])
model_fp = as.numeric(table(my_pred)[2]-model_tp) 

N = nrow(golden)

IC_lower_sens = round(as.numeric(model_sens) - 1.96*sqrt(as.numeric(model_sens)*(100-as.numeric(model_sens))/N),digits=2)
IC_upper_sens = round(as.numeric(model_sens) + 1.96*sqrt(as.numeric(model_sens)*(100-as.numeric(model_sens))/N),digits=2)
IC_sens = paste0("[",IC_lower_sens,",",IC_upper_sens,"]")

IC_lower_spec = round(as.numeric(model_spec) - 1.96*sqrt(as.numeric(model_spec)*(100-as.numeric(model_spec))/(model_fp+model_tp)),digits=2)
IC_upper_spec = round(as.numeric(model_spec) + 1.96*sqrt(as.numeric(model_spec)*(100-as.numeric(model_spec))/(model_fp+model_tp)),digits=2)
IC_spec = paste0("[",IC_lower_spec,",",IC_upper_spec,"]")

f_score = ((2*model_sens*model_spec)/(model_sens + model_spec))/100

my_data2$prediction = my_pred

saveRDS(my_data2,"mid_DEL/outputs/merge_callers_mid_DEL_prediction.rds")


### genotype error in the model ######

my_data2$GT = NA

data_call = my_data2 %>% filter(prediction=="PASS") %>% filter(PASS=="YES") %>% 
  dplyr::select(GT_gatk2,
                GT_whamg2,
                GT_svaba2,
                GT_lumpy2,
                GT_pindel2,
                GT_delly2,
                GT_strelka2,
                GT_deepvariant2,
                GT_manta2,length_manta,
                GT,GT_golden,PASS,prediction)



data_call$GT_manta2[which(data_call$length_manta<51)] = NA # indels 30 a 50 no genotype

data_call$GT_deepvariant2[data_call$GT_deepvariant2=="9/9"] = NA
data_call$GT_gatk2[data_call$GT_gatk2=="9/9"] = NA
data_call$GT_strelka2[data_call$GT_strelka2=="9/9"] = NA
data_call$GT_delly2[data_call$GT_delly2=="9/9"] = NA
data_call$GT_pindel2[data_call$GT_pindel2=="9/9"] = NA
data_call$GT_svaba2[data_call$GT_svaba2=="9/9"] = NA
data_call$GT_whamg2[data_call$GT_whamg2=="9/9"] = NA
data_call$GT_lumpy2[data_call$GT_lumpy2=="9/9"] = NA
data_call$GT_manta2[data_call$GT_manta2=="9/9"] = NA


table(data_call$GT_gatk2)
table(data_call$GT_whamg2)
table(data_call$GT_svaba2)
table(data_call$GT_lumpy2)
table(data_call$GT_pindel2)
table(data_call$GT_delly2)
table(data_call$GT_strelka2)
table(data_call$GT_deepvariant2)
table(data_call$GT_manta2)


data_call = as.matrix(data_call)

for(i in 1:nrow(data_call)){
  
  my_tab = table(as.character(data_call[i,c(1:9)]))
  
  if(length(which(my_tab==max(my_tab)))==1){
    data_call[i,11] = names(which.max(my_tab))
  }
  
  if(length(which(my_tab==max(my_tab)))>1){
    
    data_call[i,11] = NA
  }
}


head(data_call)
dim(data_call)

data_call = as.data.frame(data_call)

table(data_call$GT,useNA = "always")

# concordance genotype of the model with golden ######

data_call = data_call[!is.na(data_call$GT),]

g.error = 100 -(sum(as.character(data_call$GT)==as.character(data_call$GT_golden))/
                  nrow(data_call))*100

g.error


homo_golden = nrow(data_call[data_call$GT=="1/1",])
het_golden = nrow(data_call[data_call$GT=="0/1",])

g_error_het = round(100-nrow(data_call[data_call$GT=="0/1" & data_call$GT_golden=="0/1",])/het_golden*100,digits=2)
g_error_hom = round(100-nrow(data_call[data_call$GT=="1/1" & data_call$GT_golden=="1/1",])/homo_golden*100,digits=2)

line = as.data.frame(t(c("","","LR",model_tp,model_fp,round(model_sens,2),IC_sens,round(model_spec,2),IC_spec,round(f_score,2),round(g.error,2),
                         het_golden,g_error_het,homo_golden,g_error_hom)))

colnames(line) <- c("SV_type","N","Caller","TP","FP","Recall","CI Recall","Precision","CI Precision","F1-Score","Genotype error",
                    "N 0/1","Genotype error 0/1","N 1/1","Genotype error 1/1")

metrics_callers = rbind(metrics_callers,line)


### name callers detected and prediction #####

my_data3 = my_data2 
my_data3$callers_detected <- gsub("6-7-8-9", 6,my_data3$callers_detected )

table(my_data3$callers_detected)
table(my_data3$strategy)

dim(my_data3)

# detected by at least 1 caller

det_call =  my_data3 %>% filter(callers_detected %in% 1:9 & reciprocity >=80)

callers1 = logical_sens_spec(det_call)

line_call1=as.data.frame(t(c("","",">=1 caller",callers1,NA,NA,NA,NA,NA)))
colnames(line_call1) <- c("SV_type","N","Caller","TP","FP","Recall","CI Recall","Precision","CI Precision","F1-Score","Genotype error",
                          "N 0/1","Genotype error 0/1","N 1/1","Genotype error 1/1")

metrics_callers = rbind(metrics_callers,line_call1)


# detected by at least 2 callers and 2 strategies

det_call = my_data3 %>% filter(callers_detected %in% 2:9 & strategy %in% 2:6 & reciprocity >=80)

callers2 = logical_sens_spec(det_call)

line_call2=as.data.frame(t(c("","",">=2 caller **",callers2,NA,NA,NA,NA,NA)))
colnames(line_call2) <- c("SV_type","N","Caller","TP","FP","Recall","CI Recall","Precision","CI Precision","F1-Score","Genotype error",
                          "N 0/1","Genotype error 0/1","N 1/1","Genotype error 1/1")

metrics_callers = rbind(metrics_callers,line_call2)


# detected by at least 3 callers and 2 strategies

det_call = my_data3 %>% filter(callers_detected %in% 3:5 & strategy %in% 2:6 & reciprocity >=80)

callers3 = logical_sens_spec(det_call)

line_call3=as.data.frame(t(c("","",">=3 caller **",callers3,NA,NA,NA,NA,NA)))
colnames(line_call3) <- c("SV_type","N","Caller","TP","FP","Recall","CI Recall","Precision","CI Precision","F1-Score","Genotype error",
                          "N 0/1","Genotype error 0/1","N 1/1","Genotype error 1/1")

metrics_callers = rbind(metrics_callers,line_call3)


# detected by at least 2 callers and 2 strategies

det_call = my_data3 %>% filter(callers_detected %in% 4:5 & strategy %in% 2:6 & reciprocity >80)

callers4 = logical_sens_spec(det_call)

line_call4=as.data.frame(t(c("","",">=4 caller **",callers4,NA,NA,NA,NA,NA)))
colnames(line_call4) <- c("SV_type","N","Caller","TP","FP","Recall","CI Recall","Precision","CI Precision","F1-Score","Genotype error",
                          "N 0/1","Genotype error 0/1","N 1/1","Genotype error 1/1")

metrics_callers = rbind(metrics_callers,line_call4)

write.csv(metrics_callers,"mid_DEL/outputs/metrics_mid_DEL.csv",row.names = F)


### CV ########

library(caret)
library(e1071)

ctrl <- trainControl(method = "repeatedcv", number = 10, 
                     savePredictions = TRUE)

## 70-30

0.7*nrow(my_data2)

set.seed(2589)

training = sample(1:nrow(my_data2),0.7*nrow(my_data2))

table(my_data2$PASS,useNA = "always")

my_data2$PASS = as.factor(my_data2$PASS)

my_data2_70 = my_data2[training,]

table(my_data2_70$PASS,useNA = "always")

mod_fit <- train(factor(PASS) ~ factor(GT_deepvariant) + factor(GT_gatk) + 
                   factor(GT_strelka) + factor(GT_manta) + factor(GT_whamg) + 
                   factor(GT_delly) + factor(GT_lumpy) + factor(GT_pindel) + 
                   factor(GT_svaba) + factor(callers_detected), 
                 data=my_data2_70, method="glm", family="binomial",
                 trControl = ctrl)

my_pred = predict(mod_fit, newdata=my_data2_70)

my_pred = ifelse(my_pred=="YES",0,1)

my_pred = ifelse(my_pred==0,"PASS","NO PASS")

table(my_data2_70$PASS)
table(my_pred,my_data2_70$PASS)
table(my_pred)

N = as.numeric(table(my_data2_70$PASS)[2])

tp_cv = as.numeric(table(my_pred,my_data2_70$PASS)[2,2])
fp_cv = as.numeric(table(my_pred)[2])-tp_cv

sens_call_cv = round(tp_cv/N*100,3)
spec_call_cv = round(tp_cv/(fp_cv+tp_cv)*100,3)

f_score_cv = round((2*sens_call_cv*spec_call_cv)/(sens_call_cv+spec_call_cv)/100,3)


# que lo detecte al menos 1 caller

det_call = my_data2_70 %>% filter(callers_detected %in% 1:9)

callers1 = logical_sens_spec(det_call)

# que lo detecte al menos 2 caller y 2 estrategias a la vez

det_call = my_data2_70 %>% filter(callers_detected %in% 2:9 & strategy %in% 2:6 & reciprocity >80)

callers2 = logical_sens_spec(det_call)

# que lo detecte al menos 3 caller

det_call = my_data2_70 %>% filter(callers_detected %in% 3:9 & strategy %in% 2:6 & reciprocity >80)

callers3 = logical_sens_spec(det_call)

# que lo detecte al menos 4 caller

det_call = my_data2_70 %>% filter(callers_detected %in% 4:9 & strategy %in% 2:6 & reciprocity >80)

callers4 = logical_sens_spec(det_call)


metrics_callers = rbind(metrics_callers,
                        c("",N,"LR 10-fold CV (training 70%)",
                          c(tp_cv,fp_cv,sens_call_cv,spec_call_cv,f_score_cv),""))

metrics_callers = rbind(metrics_callers,
                        c("","",">=1 caller",callers1,""))

metrics_callers = rbind(metrics_callers,
                        c("","",">=2 caller **",callers2,""))

metrics_callers = rbind(metrics_callers,
                        c("","",">=3 caller **",callers3,""))

metrics_callers = rbind(metrics_callers,
                        c("","",">=4 caller **",callers4,""))


my_data2_30 = my_data2[-training,]

my_pred = predict(mod_fit,my_data2_30)

my_pred = ifelse(my_pred=="YES",0,1)

my_pred = ifelse(my_pred==0,"PASS","NO PASS")

table(my_data2_30$PASS)
table(my_pred,my_data2_30$PASS)
table(my_pred)

N = as.numeric(table(my_data2_30$PASS)[2])

tp_cv = as.numeric(table(my_pred,my_data2_30$PASS)[2,2])
fp_cv = as.numeric(table(my_pred)[2])-tp_cv

sens_call_cv = round(tp_cv/N*100,3)
spec_call_cv = round(tp_cv/(fp_cv+tp_cv)*100,3)

f_score_cv = round((2*sens_call_cv*spec_call_cv)/(sens_call_cv+spec_call_cv)/100,3)


# que lo detecte al menos 1 caller

det_call = my_data2_30 %>% filter(callers_detected %in% 1:9)

callers1 = logical_sens_spec(det_call)

# que lo detecte al menos 2 caller y 2 estrategias a la vez

det_call = my_data2_30 %>% filter(callers_detected %in% 2:9 & strategy %in% 2:6 & reciprocity >80)

callers2 = logical_sens_spec(det_call)

# que lo detecte al menos 3 caller

det_call = my_data2_30 %>% filter(callers_detected %in% 3:9 & strategy %in% 2:6 & reciprocity >80)

callers3 = logical_sens_spec(det_call)

# que lo detecte al menos 4 caller

det_call = my_data2_30 %>% filter(callers_detected %in% 4:9 & strategy %in% 2:6 & reciprocity >80)

callers4 = logical_sens_spec(det_call)


metrics_callers = rbind(metrics_callers,
                        c("",N,"LR (test 30%)",
                          c(tp_cv,fp_cv,sens_call_cv,spec_call_cv,f_score_cv),""))

metrics_callers = rbind(metrics_callers,
                        c("","",">=1 caller",callers1,""))

metrics_callers = rbind(metrics_callers,
                        c("","",">=2 caller **",callers2,""))

metrics_callers = rbind(metrics_callers,
                        c("","",">=3 caller **",callers3,""))

metrics_callers = rbind(metrics_callers,
                        c("","",">=4 caller **",callers4,""))

