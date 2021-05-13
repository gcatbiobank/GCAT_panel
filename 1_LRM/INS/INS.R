library(data.table)
library(dplyr)

source("ext/functions.R")

setwd("1_LRM/INS/")

# read golden ######

golden_ins = fread("data/Insilico_INS",fill=T)

golden_ins = golden_ins[,c(1,3,4),with=F]

golden_ins$V4 = nchar(golden_ins$V4)

summary(golden_ins$V4)

golden_ins$V1 = gsub("chr","",golden_ins$V1)

golden_ins2 = fread("data/Insilico_INS",fill=T,skip = 180)

golden_ins2 = golden_ins2[,c(1,3,4),with=F]

golden_ins2$V4 = nchar(golden_ins2$V4)

summary(golden_ins2$V4)

golden_ins2$V1 = gsub("chr","",golden_ins2$V1)

golden_ins = rbind(golden_ins,golden_ins2)

golden_geno_ins = NULL

for(i in names(table(golden_ins$V1))){
  
  aux = golden_ins[golden_ins$V1==i,]
  
  which_duplicated = aux$V3[which(duplicated(aux$V3))]
  
  golden_dup = unique(aux[aux$V3 %in% which_duplicated,])
  golden_no_dup = aux[!aux$V3 %in% which_duplicated,]
  
  golden_dup$genotype = "1/1"
  golden_no_dup$genotype = "0/1"
  
  aux_geno = rbind(golden_dup,golden_no_dup)
  
  golden_geno_ins = rbind(golden_geno_ins,aux_geno)
}

golden = golden_geno_ins
colnames(golden) = c("chr","start_golden","length_golden","GT_golden")
golden$lower_golden = golden$start_golden-0
golden$upper_golden = golden$start_golden+0
golden$chr_pos_golden = do.call(paste0,list(golden$chr,"_",golden$start_golden))



# read vcfs files ########

delly = fread("data/Delly_INS_30_50")
delly$V4 = as.numeric(abs(delly$V4))
delly = delly %>% filter(V4>30)
colnames(delly) = c("chr","start_delly","end","length_delly","GT_delly")
delly$lower_delly = delly$start_delly-10
delly$upper_delly = delly$start_delly+10
delly$chr_pos_delly = do.call(paste0,list(delly$chr,"_",delly$start_delly))
summary(delly$length_delly)


pindel = fread("data/Pindel_INS",fill = T)
pindel = pindel[1:60,c(1,2,5,4)]
pindel$V4 = as.numeric(nchar(pindel$V4))
colnames(pindel) = c("chr","start_pindel","GT_pindel","length_pindel")

pindel$lower_pindel = pindel$start_pindel-10
pindel$upper_pindel = pindel$start_pindel+10
pindel$chr_pos_pindel = do.call(paste0,list(pindel$chr,"_",pindel$start_pindel))
summary(pindel$length_pindel)
pindel$chr = as.character(pindel$chr)


whamg = fread("data/Whamg_INS")
colnames(whamg) = c("chr","start_whamg","length_whamg","GT_whamg")
whamg$lower_whamg = whamg$start_whamg-10
whamg$upper_whamg = whamg$start_whamg+10
whamg$chr_pos_whamg = do.call(paste0,list(whamg$chr,"_",whamg$start_whamg))
summary(whamg$length_whamg)

whamg = whamg %>% arrange(chr,chr_pos_whamg)
head(whamg,40)
whamg$start_whamg = as.numeric(whamg$start_whamg)
whamg$diff = c(0,diff(whamg$start_whamg))
whamg = whamg %>% filter(diff<0 | diff>10)


svaba = fread("data/SVaBA_INS",fill=T,skip = 5289)
svaba$length_svaba = abs(nchar(svaba$V3)-nchar(svaba$V4))
colnames(svaba)[c(1,2,5)] = c("chr","start_svaba","GT_svaba")
svaba = svaba[,c(1,2,5,6)]
svaba$lower_svaba = svaba$start_svaba-10
svaba$upper_svaba = svaba$start_svaba+10
svaba$chr_pos_svaba = do.call(paste0,list(svaba$chr,"_",svaba$start_svaba))
summary(svaba$length_svaba)

svaba = svaba %>% arrange(chr,start_svaba)
head(svaba,40)
svaba$start_svaba = as.numeric(svaba$start_svaba)
svaba$diff = c(0,diff(svaba$start_svaba))
svaba = svaba %>% filter(diff<0 | diff>10)

gatk = fread("data/Gatk_INS_30_150")
gatk$length_gatk = abs(nchar(gatk$V3)-nchar(gatk$V4))
colnames(gatk)[c(1,2,5)] = c("chr","start_gatk","GT_gatk")
gatk = gatk[,c(1,2,5,6)]
gatk$end = gatk$length_gatk+gatk$start_gatk
gatk$lower_gatk = gatk$start_gatk-10
gatk$upper_gatk = gatk$start_gatk+10
gatk$chr_pos_gatk = do.call(paste0,list(gatk$chr,"_",gatk$start_gatk))
summary(gatk$length_gatk)

strelka = fread("data/Strelka_INS_30_50")
strelka$length_strelka = abs(nchar(strelka$V3)-nchar(strelka$V4))
strelka = strelka %>% filter(length_strelka>30 & length_strelka<151)
colnames(strelka)[c(1,2,5)] = c("chr","start_strelka","GT_strelka")
strelka = strelka[,c(1,2,5,6)]
strelka$end = strelka$length_strelka+strelka$start_strelka
strelka$lower_strelka = strelka$start_strelka-10
strelka$upper_strelka = strelka$start_strelka+10
strelka$chr_pos_strelka = do.call(paste0,list(strelka$chr,"_",strelka$start_strelka))
summary(strelka$length_strelka)


manta = fread("data/Manta_INS")
manta$V4 = as.numeric(abs(manta$V4))
manta$V3 = as.character(manta$V3)
manta$V3 = "0/1"
colnames(manta) = c("chr","start_manta","GT_manta","length_manta")
manta2 = fread("data/Manta_INS",skip=147)
head(manta2)
colnames(manta2) = c("chr","start_manta","GT_manta")
manta2$length_manta = 0
manta = rbind(manta,manta2)
manta$lower_manta = manta$start_manta-10
manta$upper_manta = manta$start_manta+10
manta$chr_pos_manta = do.call(paste0,list(manta$chr,"_",manta$start_manta))

manta = manta %>% arrange(chr,start_manta)
head(manta,40)
manta$start_manta = as.numeric(manta$start_manta)
manta$diff = c(0,diff(manta$start_manta))
manta = manta %>% filter(diff<0 | diff>10)

popins = fread("data/Popins_INS")
colnames(popins) = c("chr","start_popins","GT_popins")
popins$lower_popins = popins$start_popins-200
popins$upper_popins = popins$start_popins+200
popins$chr_pos_popins = do.call(paste0,list(popins$chr,"_",popins$start_popins))
popins$length_popins = 0
summary(popins$length_popins)

popins = popins %>% arrange(chr,start_popins)
head(popins,40)
popins$start_popins = as.numeric(popins$start_popins)
popins$diff = c(0,diff(popins$start_popins))
popins = popins %>% filter(diff<0 | diff>200)


### sensitivity-specifitiy-geno error for each caller ######

# gatk 

sens_gatk = sensitivity_precision_geno_error(gatk,golden,"insertion","gatk")

# strelka

strelka = strelka %>% as.data.table()

sens_strelka = sensitivity_precision_geno_error(strelka,golden,"insertion","strelka")

# popins

popins = popins %>% as.data.table()

sens_popins = sensitivity_precision_geno_error(popins,golden,"insertion","popins")

# pindel

pindel = pindel %>% as.data.table()

sens_pindel = sensitivity_precision_geno_error(pindel,golden,"insertion","pindel")

# manta

manta = manta %>% as.data.table()

sens_manta = sensitivity_precision_geno_error(manta,golden,"insertion","manta")

# whamg

whamg = whamg %>% as.data.table()

sens_whamg = sensitivity_precision_geno_error(whamg,golden,"insertion","whamg")

# delly

delly = delly %>% as.data.table()

sens_delly = sensitivity_precision_geno_error(delly,golden,"insertion","delly")

# svaba

svaba = svaba %>% as.data.table()

sens_svaba = sensitivity_precision_geno_error(svaba,golden,"insertion","svaba")


metrics_callers = matrix(c(t(unlist(sens_delly))[c(3,4,1,2,5,6)],
                           t(unlist(sens_gatk))[c(3,4,1,2,5,6)],
                           t(unlist(sens_popins))[c(3,4,1,2,5,6)],
                           t(unlist(sens_manta))[c(3,4,1,2,5,6)],
                           t(unlist(sens_pindel))[c(3,4,1,2,5,6)],
                           t(unlist(sens_strelka))[c(3,4,1,2,5,6)],
                           t(unlist(sens_svaba))[c(3,4,1,2,5,6)],
                           t(unlist(sens_whamg))[c(3,4,1,2,5,6)]),byrow=T,ncol = 6)

metrics_callers = cbind(c("Delly","Gatk","Popins","Manta","Pindel","Strelka","Svaba","Whamg"),
                        round(metrics_callers,2))

colnames(metrics_callers) = c("Caller","TP","FP","Sensitivity","Precision","F1-Score","Genotype error")

metrics_callers = metrics_callers %>% as.data.frame() %>% arrange(desc(`F1-Score`))

metrics_callers = cbind(N=c(nrow(golden),rep("",nrow(metrics_callers)-1)),
                        metrics_callers)

metrics_callers = cbind(SV_type=c("Insertions",rep("",nrow(metrics_callers)-1)),
                        metrics_callers)

metrics_callers[,1:3] = apply(metrics_callers[,1:3], 2, function(x) as.character(x));
metrics_callers[,4:9] = apply(metrics_callers[,4:9], 2, function(x) as.numeric(x));


## prepare BBDD ####

popins$ID = paste0("insertion_",1:nrow(popins))
pindel$ID = "none"
manta$ID = "none"
whamg$ID = "none"
gatk$ID = "none"
svaba$ID = "none"
strelka$ID = "none"
delly$ID = "none"
golden$ID = "none"


my_data = merge_callers(pindel,popins,callers_ref = c("popins"),
                        caller_to_merge = c("pindel")) 

my_data = merge_callers(manta,my_data,callers_ref = c("popins","pindel"),
                        caller_to_merge = c("manta")) 

my_data = merge_callers(whamg,my_data,callers_ref = c("popins","pindel","manta"),
                        caller_to_merge = c("whamg")) 

my_data = merge_callers(gatk,my_data,callers_ref = c("popins","pindel","manta","whamg"),
                        caller_to_merge = c("gatk")) 

my_data = merge_callers(svaba,my_data,callers_ref = c("popins","pindel","manta","whamg",
                                                      "gatk"),
                        caller_to_merge = c("svaba")) 

my_data = merge_callers(strelka,my_data,callers_ref = c("popins","pindel","manta","whamg",
                                                      "gatk","svaba"),
                        caller_to_merge = c("strelka")) 

my_data = merge_callers(delly,my_data,callers_ref = c("popins","pindel","manta","whamg",
                                                      "gatk","svaba","strelka"),
                        caller_to_merge = c("delly")) 


my_data2 = merge_callers(golden,my_data,callers_ref = c("popins","pindel","manta","whamg",
                                                        "gatk","svaba","strelka","delly"),
                         caller_to_merge = c("golden"))


dim(my_data)
dim(my_data2)

my_data2 = unique(my_data2)


## add number of callers detected #####

my_data2$length_svaba = 0

my_data2$callers_detected = n_callers_detected(my_data2,
                                               c("popins","pindel","manta","whamg",
                                                 "gatk","svaba","strelka","delly"))

table(my_data2$callers_detected)

## add strategy ####

strategies = fread("strategies.csv")

my_data2$strategy = strategy(my_data2,
                             c("popins","pindel","manta","whamg",
                               "gatk","svaba","strelka","delly"),strategies)

table(my_data2$strategy)

### colineality #####

table(my_data2$callers_detected,my_data2$strategy)
chisq.test(table(my_data2$callers_detected,my_data2$strategy))

## prepare model ####

my_data2$PASS = "YES"

my_data2$PASS[is.na(my_data2$length_golden)] = "NO"

dim(golden)

my_data2 = my_data2[-which(duplicated(my_data2$ID)),]

my_data2 = my_data2 %>% filter(callers_detected!=0)


table(my_data2$PASS)

my_data2$GT_gatk[is.na(my_data2$GT_gatk)] = "9/9"
my_data2$GT_strelka[is.na(my_data2$GT_strelka)] = "9/9"
my_data2$GT_manta[is.na(my_data2$GT_manta)] = "9/9"
my_data2$GT_whamg[is.na(my_data2$GT_whamg)] = "9/9"
my_data2$GT_delly[is.na(my_data2$GT_delly)] = "9/9"
my_data2$GT_popins[is.na(my_data2$GT_popins)] = "9/9"
my_data2$GT_pindel[is.na(my_data2$GT_pindel)] = "9/9"
my_data2$GT_svaba[is.na(my_data2$GT_svaba)] = "9/9"

my_data2$GT_whamg[my_data2$GT_whamg==""] = "9/9"
my_data2$GT_popins[my_data2$GT_popins=="0/0"] = "9/9"

table(my_data2$GT_gatk,useNA = "always")
table(my_data2$GT_strelka,useNA = "always")
table(my_data2$GT_manta,useNA = "always")
table(my_data2$GT_whamg,useNA = "always")
table(my_data2$GT_delly,useNA = "always")
table(my_data2$GT_popins,useNA = "always")
table(my_data2$GT_svaba,useNA = "always")

table(my_data2$PASS)

my_data2$PASS = ifelse(my_data2$PASS=="YES",0,1)

my_model = glm(PASS ~ factor(GT_gatk) + 
                 factor(GT_strelka)+
                 factor(GT_manta)+
                 factor(GT_whamg)+
                 factor(GT_delly)+
                 factor(GT_popins)+
                 factor(GT_pindel)+
                 factor(GT_svaba) + 
                 factor(callers_detected), 
                 #factor(strategy)+
               data=my_data2,family = "binomial")

my_data2$PASS = ifelse(my_data2$PASS==0,"YES","NO")

summary(my_model)

my_pred = predict(my_model,my_data2,type="response")

my_pred = ifelse(my_pred>0.5,1,0)

my_pred = ifelse(my_pred==1,"NO PASS","PASS")

table(my_data2$PASS)
table(my_pred,my_data2$PASS)
table(my_pred)

length(unique(golden$chr_pos))

# without CHR and BP
#sensitivity
model_sens = (table(my_pred,my_data2$PASS)[2,2]/nrow(golden))*100
#specificity
model_spec = (table(my_pred,my_data2$PASS)[2,2])/(table(my_pred,my_data2$PASS)[2,2]+
                                                    table(my_pred,my_data2$PASS)[2,1])*100

f_score = ((2*model_sens*model_spec)/(model_sens + model_spec))/100


table(my_pred,my_data2$callers_detected)


### choose best model ####

my_data2$PASS = ifelse(my_data2$PASS=="YES",0,1)

full.model =  glm(PASS ~ factor(GT_gatk) + 
                    factor(GT_strelka)+
                    factor(GT_manta)+
                    factor(GT_whamg)+
                    factor(GT_delly)+
                    factor(GT_popins)+
                    factor(GT_pindel)+
                    factor(GT_svaba) + 
                    factor(callers_detected)+ factor(strategy),
                  data=my_data2,family = "binomial")

base.model = glm(PASS ~ factor(GT_gatk) + 
                   factor(GT_strelka)+
                   factor(GT_manta)+
                   factor(GT_whamg)+
                   factor(GT_delly)+
                   factor(GT_popins)+
                   factor(GT_pindel)+
                   factor(GT_svaba),
                 data=my_data2,family = "binomial")

step(full.model,
     scope = list(lower = formula(base.model),
                  upper = formula(full.model)),
     direction = "backward")

best.model = glm(formula = PASS ~ factor(GT_gatk) + factor(GT_strelka) + factor(GT_manta) + 
                   factor(GT_whamg) + factor(GT_delly) + factor(GT_popins) + 
                   factor(GT_pindel) + factor(GT_svaba) + factor(callers_detected), 
                 family = "binomial", data = my_data2)

my_data2$PASS = ifelse(my_data2$PASS==0,"YES","NO")

my_pred = predict(best.model,my_data2,type="response")

my_pred = ifelse(my_pred>0.5,1,0)

my_pred = ifelse(my_pred==1,"NO PASS","PASS")

table(my_data2$PASS)
table(my_pred,my_data2$PASS)
table(my_pred)

N = length(unique(golden$chr_pos))

model_tp = as.numeric(table(my_pred,my_data2$PASS)[2,2])
model_fp = as.numeric(table(my_pred)[2]-model_tp)  
#sensitivity
model_sens = round(model_tp/N*100,3)
#specificity
model_spec = round(model_tp/(model_fp+model_tp)*100,3)

model_tp;model_fp;model_sens;model_spec

f_score = round(((2*model_sens*model_spec)/(model_sens + model_spec))/100,3)

my_data2$prediction = my_pred


# detected by at least 1 caller

det_call = my_data2 %>% filter(callers_detected %in% 1:9)

callers1 = logical_sens_spec(det_call)

# detected by at least 2 callers and 2 strategies

det_call = my_data2 %>% filter(callers_detected %in% 2:9 & strategy %in% 2:6)

callers2 = logical_sens_spec(det_call)

# detected by at least 3 callers and 2 strategies

det_call = my_data2 %>% filter(callers_detected %in% 3:9 & strategy %in% 2:6)

callers3 = logical_sens_spec(det_call)

# detected by at least 4 callers and 2 strategies

det_call = my_data2 %>% filter(callers_detected %in% 4:9 & strategy %in% 2:6)

callers4 = logical_sens_spec(det_call)

metrics_callers = rbind(metrics_callers,
                        c("","","LR",model_tp,model_fp,model_sens,model_spec,f_score,
                          ""))

metrics_callers = rbind(metrics_callers,
                        c("","",">=1 caller",callers1,""))

metrics_callers = rbind(metrics_callers,
                        c("","",">=2 caller **",callers2,""))

metrics_callers = rbind(metrics_callers,
                        c("","",">=3 caller **",callers3,""))

metrics_callers = rbind(metrics_callers,
                        c("","",">=4 caller **",callers4,""))


## genotype error in the model ######

data_call = my_data2 %>% filter(prediction=="PASS") %>% filter(PASS=="YES") %>% 
  dplyr::select(GT_whamg,
                GT_strelka,
                GT_pindel,
                GT_delly,
                GT_manta,
                GT_popins,
                GT_gatk,
                GT_svaba,
                GT_golden,PASS,prediction)

data_call$GT_delly[data_call$GT_delly=="9/9"] = NA
data_call$GT_pindel[data_call$GT_pindel=="9/9"] = NA
data_call$GT_whamg[data_call$GT_whamg=="9/9"] = NA
data_call$GT_strelka[data_call$GT_strelka=="9/9"] = NA
data_call$GT_manta[data_call$GT_manta=="9/9"] = NA
data_call$GT_svaba[data_call$GT_svaba=="9/9"] = NA
data_call$GT_gatk[data_call$GT_gatk=="9/9"] = NA
data_call$GT_popins[data_call$GT_popins=="9/9"] = NA

table(data_call$GT_whamg)
table(data_call$GT_strelka)
table(data_call$GT_pindel)
table(data_call$GT_delly)
table(data_call$GT_manta)
table(data_call$GT_svaba)
table(data_call$GT_gatk)
table(data_call$GT_popins)


data_call$GT = NA

data_call = as.matrix(data_call)

for(i in 1:nrow(data_call)){
  
  my_tab = table(as.character(data_call[i,c(1:8)]))
  
  if(length(which(my_tab==max(my_tab)))==1){
    data_call[i,12] = names(which.max(my_tab))
  }
  
  if(length(which(my_tab==max(my_tab)))>1){
    
    data_call[i,12] = NA
  }
}


head(data_call)
dim(data_call)

data_call = as.data.frame(data_call)

table(data_call$GT,useNA = "always")

data_call = data_call[!is.na(data_call$GT),]

model_g_error = 100 -(sum(as.character(data_call$GT)==as.character(data_call$GT_golden))/
                        nrow(data_call))*100

model_g_error



# caret package cv k-fold ####

library(caret)
library(e1071)

ctrl <- trainControl(method = "repeatedcv", number = 10, 
                     savePredictions = TRUE)

## 70-30

0.7*nrow(my_data2)

set.seed(2589)

training = sample(1:nrow(my_data2),1223)

table(my_data2$PASS,useNA = "always")

my_data2$PASS = as.factor(my_data2$PASS)

my_data2_70 = my_data2[training,]

table(my_data2_70$PASS,useNA = "always")

mod_fit <- train(factor(PASS) ~ factor(GT_gatk) + factor(GT_strelka) + factor(GT_manta) + 
                   factor(GT_whamg) + factor(GT_delly) + factor(GT_popins) + 
                   factor(GT_pindel) + factor(GT_svaba) + factor(callers_detected), 
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


# detected by at least 1 caller

det_call = my_data2_70 %>% filter(callers_detected %in% 1:9)

callers1 = logical_sens_spec(det_call)

# detected by at least 2 callers and 2 strategies

det_call = my_data2_70 %>% filter(callers_detected %in% 2:9 & strategy %in% 2:6)

callers2 = logical_sens_spec(det_call)

# detected by at least 3 callers and 2 strategies

det_call = my_data2_70 %>% filter(callers_detected %in% 3:9 & strategy %in% 2:6)

callers3 = logical_sens_spec(det_call)

# detected by at least 4 callers and 2 strategies

det_call = my_data2_70 %>% filter(callers_detected %in% 4:9 & strategy %in% 2:6)

callers4 = logical_sens_spec(det_call)


metrics_callers = rbind(metrics_callers,
                        c("",N,"LR (training 70%)",
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


# detected by at least 1 caller

det_call = my_data2_30 %>% filter(callers_detected %in% 1:9)

callers1 = logical_sens_spec(det_call)

# detected by at least 2 callers and 2 strategies

det_call = my_data2_30 %>% filter(callers_detected %in% 2:9 & strategy %in% 2:6)

callers2 = logical_sens_spec(det_call)

# detected by at least 3 callers and 2 strategies

det_call = my_data2_30 %>% filter(callers_detected %in% 3:9 & strategy %in% 2:6)

callers3 = logical_sens_spec(det_call)

# detected by at least 4 callers and 2 strategies

det_call = my_data2_30 %>% filter(callers_detected %in% 4:9 & strategy %in% 2:6)

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

metrics_callers$Sensitivity = as.numeric(metrics_callers$Sensitivity)
metrics_callers$Precision = as.numeric(metrics_callers$Precision)

metrics_callers$F1.Score = ((2*metrics_callers$Sensitivity*metrics_callers$Precision)/
                              (metrics_callers$Sensitivity + metrics_callers$Precision))
