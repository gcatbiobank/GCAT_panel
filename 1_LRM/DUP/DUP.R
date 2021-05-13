library(data.table)
library(dplyr)

source("ext/functions.R")

setwd("1_LRM/DUP/")

# read golden ######

golden_duplications = fread("data/Insilico_DUP")
head(golden_duplications)
golden_duplications$V1 = gsub("chr","",golden_duplications$V1)

golden_geno_dup = NULL

for(i in names(table(golden_duplications$V1))){
  
  aux = golden_duplications[golden_duplications$V1==i,]
  
  which_duplicated = aux$V2[which(duplicated(aux$V2))]
  which_duplicated2 = aux$V3[which(duplicated(aux$V3))]
  
  which_duplicated = unique(c(which_duplicated,which_duplicated2))
  
  golden_dup = unique(aux[aux$V2 %in% which_duplicated,])
  golden_no_dup = aux[!aux$V2 %in% which_duplicated,]
  
  if(nrow(golden_dup)>0){golden_dup$genotype = "1/1"
  golden_no_dup$genotype = "0/1"
  
  aux_geno = rbind(golden_dup,golden_no_dup)}
  
  
  if(nrow(golden_dup)==0){
  golden_no_dup$genotype = "0/1"
  
  aux_geno = golden_no_dup}
  
  golden_geno_dup = rbind(golden_geno_dup,aux_geno)
}

golden = golden_geno_dup

golden$diff = golden$V3-golden$V2

for(i in 1:nrow(golden)){
  
  if(golden$diff[i]<0){
    aux = golden$V2[i]
    aux2 = golden$V3[i]
    
    golden$V2[i] = aux2
    golden$V3[i] = aux
    
  }
  
}
golden$diff = golden$V3-golden$V2

summary(golden$diff)

golden$V3 = golden$diff
dim(golden)
golden = golden[,1:4,with=F]
head(golden)
colnames(golden) = c("chr","start_golden","length_golden","GT_golden")

golden$lower_golden = golden$start_golden-0
golden$upper_golden = golden$start_golden+0
golden$chr_pos_golden = do.call(paste0,list(golden$chr,"_",golden$start_golden))
golden= distinct(golden, chr_pos_golden,.keep_all = TRUE)

golden$end_golden = golden$start_golden + golden$length_golden
golden$diff = c(0,diff(golden$start_golden))


golden$diff = ifelse(golden$diff<0,0,golden$diff)

golden$window_bp = 10

golden = golden %>% filter(diff!=1)

golden = golden %>% as.data.table()

dim(golden)

# read vcfs files ########

delly = fread("data/Delly_DUP")
colnames(delly) = c("chr","start_delly","end","length_delly","GT_delly")
delly$lower_delly = delly$start_delly-10
delly$upper_delly = delly$start_delly+10
delly$chr_pos_delly = do.call(paste0,list(delly$chr,"_",delly$start_delly))
summary(delly$length_delly)

pindel = fread("data/Pindel_DUP",fill = T)
colnames(pindel) = c("chr","start_pindel","end","length_pindel","GT_pindel")
pindel$lower_pindel = pindel$start_pindel-10
pindel$upper_pindel = pindel$start_pindel+10
pindel$chr_pos_pindel = do.call(paste0,list(pindel$chr,"_",pindel$start_pindel))
summary(pindel$length_pindel)

whamg = fread("data/Whamg_DUP")
colnames(whamg) = c("chr","start_whamg","end","length_whamg","GT_whamg")
whamg$lower_whamg = whamg$start_whamg-10
whamg$upper_whamg = whamg$start_whamg+10
whamg$chr_pos_whamg = do.call(paste0,list(whamg$chr,"_",whamg$start_whamg))
summary(whamg$length_whamg)


svaba = fread("data/SVaBA_DUP")
svaba$length = svaba$V2-svaba$V3

# change start and end

svaba$length_aux = svaba$V3-svaba$V2

for(i in 1:nrow(svaba)){
  
  if(svaba$length_aux[i]<0){
    
    aux = svaba$V2[i]
    aux2 = svaba$V3[i]
    
    svaba$V2[i] = aux2
    svaba$V3[i] = aux
    
  }
  
}

svaba$length_aux = svaba$V3-svaba$V2
length(which(golden$length<=0))
svaba = svaba[,c(1:4,6)]
colnames(svaba) = c("chr","start_svaba","end_svaba","GT_svaba","length_svaba")
svaba$lower_svaba = svaba$start_svaba-10
svaba$upper_svaba = svaba$start_svaba+10
dim(svaba)
svaba = unique(svaba)
svaba$chr_pos_svaba = do.call(paste0,list(svaba$chr,"_",svaba$start_svaba))

manta = fread("data/Manta_DUP")
colnames(manta) = c("chr","start_manta","end","length_manta","GT_manta")
manta$lower_manta = manta$start_manta-20
manta$upper_manta = manta$start_manta+20
manta$chr_pos_manta = do.call(paste0,list(manta$chr,"_",manta$start_manta))

lumpy = fread("data/Lumpy_DUP")
colnames(lumpy) = c("chr","start_lumpy","end","length_lumpy","GT_lumpy")
lumpy$lower_lumpy = lumpy$start_lumpy-50
lumpy$upper_lumpy = lumpy$start_lumpy+50
lumpy$chr_pos_lumpy = do.call(paste0,list(lumpy$chr,"_",lumpy$start_lumpy))

cnvnator = fread("data/CNVnator_DUP")
colnames(cnvnator) = c("chr","start_cnvnator","end","length_cnvnator","GT_cnvnator")
cnvnator$lower_cnvnator = cnvnator$start_cnvnator-100
cnvnator$upper_cnvnator = cnvnator$start_cnvnator+100
cnvnator$chr_pos_cnvnator = do.call(paste0,list(cnvnator$chr,"_",cnvnator$start_cnvnator))

head(golden)
### sensitivity-specifitiy-geno error for each caller ######

## without translocations intra

# pindel

sens_pindel = sensitivity_precision_geno_error_deletions_duplications(pindel,golden,"duplication","pindel")

# manta

sens_manta = sensitivity_precision_geno_error_deletions_duplications(manta,
                                                                     golden,"duplication","manta")

# whamg

sens_whamg = sensitivity_precision_geno_error_deletions_duplications(whamg,golden,"duplication","whamg")

# delly

sens_delly = sensitivity_precision_geno_error_deletions_duplications(delly,golden,"duplication","delly")

# svaba

sens_svaba = sensitivity_precision_geno_error_deletions_duplications(svaba,golden,"duplication","svaba")

# lumpy

sens_lumpy = sensitivity_precision_geno_error_deletions_duplications(lumpy,golden,"duplication","lumpy")

# cnvnator

sens_cnvnator = sensitivity_precision_geno_error_deletions_duplications(cnvnator,golden,"duplication","cnvnator")

metrics_callers = matrix((c(t(unlist(sens_delly))[c(3,4,1,12,2,13,5,6,7,9,8,10)],
                           t(unlist(sens_lumpy))[c(3,4,1,12,2,13,5,6,7,9,8,10)],
                           t(unlist(sens_cnvnator))[c(3,4,1,12,2,13,5,6,7,9,8,10)],
                           t(unlist(sens_manta))[c(3,4,1,12,2,13,5,6,7,9,8,10)],
                           t(unlist(sens_pindel))[c(3,4,1,12,2,13,5,6,7,9,8,10)],
                           t(unlist(sens_svaba))[c(3,4,1,12,2,13,5,6,7,9,8,10)],
                           t(unlist(sens_whamg))[c(3,4,1,12,2,13,5,6,7,9,8,10)])),byrow=T,ncol = 12)

metrics_callers1 = as.data.frame(metrics_callers)

c =apply(metrics_callers, 2, as.numeric)

c = cbind(c("Delly","Lumpy","CNVnator","Manta","Pindel",
            "Svaba","Whamg"),
          round(c,2))


colnames(c) = c("Caller","TP","FP","Recall","CI Recall","Precision","CI Precision","F1-Score","Genotype error",
                              "N 0/1","Genotype error 0/1","N 1/1","Genotype error 1/1")

c = as.data.frame(c)

c$`CI Recall` = metrics_callers1$V4
c$`CI Precision` = metrics_callers1$V6

c = c %>% as.data.frame() %>% arrange(desc(`F1-Score`))

c = cbind(N=c(nrow(golden),rep("",nrow(c)-1)),
          c)

c = cbind(SV_type=c("Duplications",rep("",nrow(c)-1)),
          c)


c[,1:3] = apply(c[,1:3], 2, function(x) as.character(x));

metrics_callers = c


## prepare BBDD ####

lumpy$ID = paste0("duplication_",1:nrow(lumpy))
manta$ID = "none"
pindel$ID = "none"
delly$ID = "none"
whamg$ID = "none"
svaba$ID = "none"
cnvnator$ID = "none"
golden$ID = "none"


my_data = merge_callers_deletions_duplications(manta,lumpy,callers_ref = c("lumpy"),
                        caller_to_merge = c("manta"),svtype="duplication",repro = 0.8) 

my_data = merge_callers_deletions_duplications(pindel,my_data,callers_ref = c("lumpy","manta"),
                        caller_to_merge = c("pindel"),svtype="duplication",repro = 0.8) 

my_data = merge_callers_deletions_duplications(delly,my_data,callers_ref = c("lumpy","manta","pindel"),
                        caller_to_merge = c("delly"),svtype="duplication",repro = 0.8) 

my_data = merge_callers_deletions_duplications(whamg,my_data,callers_ref = c("lumpy","manta","pindel","delly"),
                        caller_to_merge = c("whamg"),svtype="duplication",repro = 0.8) 

my_data = merge_callers_deletions_duplications(svaba,my_data,callers_ref = c("lumpy","manta","pindel","delly",
                                                      "whamg"),
                       caller_to_merge = c("svaba"),svtype="duplication",repro = 0.8) 

my_data = merge_callers_deletions_duplications(cnvnator,my_data,callers_ref = c("lumpy","manta","pindel","delly",
                                                      "whamg","svaba"),
                        caller_to_merge = c("cnvnator"),svtype="duplication",repro = 0.8) 

my_data2 = merge_callers_deletions_duplications(golden,my_data,callers_ref = c("lumpy","manta","pindel","delly",
                                                        "whamg","svaba","cnvnator"),
                         caller_to_merge = c("golden"),svtype="duplication",repro = 0.8)

dim(my_data)
dim(my_data2)

my_data2 = unique(my_data2)


### consensus length #######

my_data2$length = consensus_length(my_data2,
                                   callers = c("lumpy","manta","pindel","delly",
                                               "whamg","svaba","cnvnator"))

summary(my_data2$length)

aux = cut(my_data2$length,
          breaks = c(30,150,500,1000,2000,3000,Inf),right = F)

my_data2$length_stretch = aux

table(my_data2$length_stretch,useNA = "always")


## add reciprocity ###

my_data2$reciprocity = reprocicity(my_data2,
                                   c("lumpy","manta","pindel","delly",
                                     "whamg","svaba","cnvnator"))
summary(my_data2$reciprocity)


## add number of callers detected #####

my_data2$callers_detected = n_callers_detected(my_data2,
                                               callers = c("lumpy","manta","pindel","delly",
                                                           "whamg","svaba","cnvnator"))

table(my_data2$callers_detected)

my_data2$callers_detected[my_data2$callers_detected=="4"] = "4-5-6-7"
my_data2$callers_detected[my_data2$callers_detected=="5"] = "4-5-6-7"
my_data2$callers_detected[my_data2$callers_detected=="6"] = "4-5-6-7"
my_data2$callers_detected[my_data2$callers_detected=="7"] = "4-5-6-7"


## add strategy ####

strategies = fread("ext/strategies.csv")

my_data2$strategy = strategy(my_data2,
                             callers = c("lumpy","manta","pindel","delly",
                                         "whamg","svaba","cnvnator"),
                             strategies)

table(my_data2$strategy)

my_data2$strategy[my_data2$strategy=="3"] = "3-4-5"
my_data2$strategy[my_data2$strategy=="4"] = "3-4-5"
my_data2$strategy[my_data2$strategy=="5"] = "3-4-5"


### colineality #####

table(my_data2$callers_detected,my_data2$strategy)
chisq.test(table(my_data2$callers_detected,my_data2$strategy))


### consensus start #####

my_data2$start = consensus_start_duplications(my_data2,
                                             callers = c("lumpy","manta","pindel","delly",
                                                         "whamg","svaba","cnvnator"))

## prepare model choose de best model####

my_data2$PASS = "YES"

my_data2$PASS[is.na(my_data2$length_golden)] = "NO"

table(my_data2$PASS)
dim(my_data2)
dim(golden)

my_data2 = my_data2[-which(duplicated(my_data2$ID)),]

my_data2 = my_data2 %>% filter(callers_detected!=0)

table(my_data2$PASS)

table(my_data2$length_stretch,my_data2$PASS)
table(my_data2$strategy,my_data2$PASS)
table(my_data2$callers_detected,my_data2$PASS) ## callers detected and strategy are not informative


my_data2$GT_manta[is.na(my_data2$GT_manta)] = "9/9"
my_data2$GT_whamg[is.na(my_data2$GT_whamg)] = "9/9"
my_data2$GT_delly[is.na(my_data2$GT_delly)] = "9/9"
my_data2$GT_lumpy[is.na(my_data2$GT_lumpy)] = "9/9"
my_data2$GT_pindel[is.na(my_data2$GT_pindel)] = "9/9"
my_data2$GT_svaba[is.na(my_data2$GT_svaba)] = "9/9"
my_data2$GT_cnvnator[is.na(my_data2$GT_cnvnator)] = "9/9"

my_data2$GT_cnvnator[my_data2$GT_cnvnator == "./1"] = "9/9"


my_data2$GT_manta2 = my_data2$GT_manta
my_data2$GT_whamg2 = my_data2$GT_whamg
my_data2$GT_delly2 = my_data2$GT_delly
my_data2$GT_lumpy2 = my_data2$GT_lumpy
my_data2$GT_pindel2 = my_data2$GT_pindel
my_data2$GT_svaba2 = my_data2$GT_svaba
my_data2$GT_cnvnator2 = my_data2$GT_cnvnator

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
my_data2$GT_cnvnator[my_data2$GT_cnvnator=="0/1"] = "0/1-1/1"
my_data2$GT_cnvnator[my_data2$GT_cnvnator=="1/1"] = "0/1-1/1"


table(my_data2$GT_manta,useNA = "always")
table(my_data2$GT_whamg,useNA = "always")
table(my_data2$GT_delly,useNA = "always")
table(my_data2$GT_svaba,useNA = "always")
table(my_data2$GT_lumpy,useNA = "always")
table(my_data2$GT_cnvnator,useNA = "always")
table(my_data2$GT_pindel,useNA = "always")


table(my_data2$PASS)

my_data2$PASS = ifelse(my_data2$PASS=="YES",0,1)

best.model = glm(formula = PASS ~ factor(GT_lumpy) + factor(GT_cnvnator) + 
                   factor(GT_manta) + factor(GT_whamg) + factor(GT_delly) + 
                   factor(GT_pindel) + factor(GT_svaba) + factor(length_stretch) + reciprocity, 
                 family = "binomial", data = my_data2)


my_data2$PASS = ifelse(my_data2$PASS==0,"YES","NO")

summary(best.model)

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
model_sens = round(model_tp/N*100,2)
#specificity
model_spec = round(model_tp/(model_fp+model_tp)*100,2)

IC_lower_sens = round(as.numeric(model_sens) - 1.96*sqrt(as.numeric(model_sens)*(100-as.numeric(model_sens))/N),digits=2)
IC_upper_sens = round(as.numeric(model_sens) + 1.96*sqrt(as.numeric(model_sens)*(100-as.numeric(model_sens))/N),digits=2)
IC_sens = paste0("[",IC_lower_sens,",",IC_upper_sens,"]")

IC_lower_spec = round(as.numeric(model_spec) - 1.96*sqrt(as.numeric(model_spec)*(100-as.numeric(model_spec))/(model_fp+model_tp)),digits=2)
IC_upper_spec = round(as.numeric(model_spec) + 1.96*sqrt(as.numeric(model_spec)*(100-as.numeric(model_spec))/(model_fp+model_tp)),digits=2)
IC_spec = paste0("[",IC_lower_spec,",",IC_upper_spec,"]")

f_score = round(((2*model_sens*model_spec)/(model_sens + model_spec))/100,2)

my_data2$prediction = my_pred
head(my_data2)


### name callers detected and prediction #####

my_data2$name_callers_detected = name_callers_detected(as.data.table(my_data2), callers = c("lumpy","manta","pindel","delly",
                                                                             "whamg","svaba","cnvnator"))

table(my_data2$name_callers_detected,my_data2$prediction) 

head(my_data2)

my_data3 = my_data2 %>% filter(name_callers_detected != "svaba") #Svaba is eliminated because is not able to classify by type the SV
my_data3$callers_detected <- gsub("4-5-6-7", 4,my_data3$callers_detected )
my_data3$strategy <- gsub("3-4-5", 3,my_data3$strategy )

table(my_data3$callers_detected)

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

det_call = my_data3 %>% filter(callers_detected %in% 3:5 & strategy %in% 2:6 & reciprocity >80)

callers3 = logical_sens_spec(det_call)

line_call3=as.data.frame(t(c("","",">=3 caller **",callers3,NA,NA,NA,NA,NA)))
colnames(line_call3) <- c("SV_type","N","Caller","TP","FP","Recall","CI Recall","Precision","CI Precision","F1-Score","Genotype error",
                          "N 0/1","Genotype error 0/1","N 1/1","Genotype error 1/1")

metrics_callers = rbind(metrics_callers,line_call3)

                
# detected by at least 4 callers and 2 strategies

det_call = my_data3 %>% filter(callers_detected %in% 4:5 & strategy %in% 2:6 & reciprocity >80)

callers4 = logical_sens_spec(det_call)

line_call4=as.data.frame(t(c("","",">=4 caller **",callers4,NA,NA,NA,NA,NA)))
colnames(line_call4) <- c("SV_type","N","Caller","TP","FP","Recall","CI Recall","Precision","CI Precision","F1-Score","Genotype error",
                          "N 0/1","Genotype error 0/1","N 1/1","Genotype error 1/1")

metrics_callers = rbind(metrics_callers,line_call4)



### CV ######

library(caret)
library(e1071)

ctrl <- trainControl(method = "repeatedcv", number = 10, 
                     savePredictions = TRUE)

## 70-30


set.seed(2589)

training = sample(1:nrow(z),2392) 

my_data2_70 = my_data2[training,]

table(my_data2_70$length_stretch,useNA = "always")
table(my_data2_70$PASS)



mod_fit <- train(factor(PASS) ~ factor(GT_lumpy)+ factor(GT_pindel) + factor(GT_cnvnator) + 
                   factor(GT_manta) + factor(GT_whamg) + factor(GT_delly) + 
                   factor(GT_svaba) + factor(GT_lumpy) + factor(length_stretch)+ reciprocity, 
                 data=my_data2_70, method="glm", family="binomial",
                 trControl = ctrl)

summary(mod_fit)

my_pred = predict(mod_fit, newdata=my_data2_70)
length(my_pred)
my_pred = ifelse(my_pred=="YES",0,1)

my_pred = ifelse(my_pred==0,"PASS","NO PASS")

table(my_data2_70$PASS)
table(my_pred,my_data2_70$PASS)
table(my_pred)

N = as.numeric(table(my_data2_70$PASS)[2])

tp_cv = as.numeric(table(my_pred,my_data2_70$PASS)[2,2])
fp_cv = as.numeric(table(my_pred)[2])-tp_cv

sens_call_cv = round(tp_cv/N*100,2)
spec_call_cv = round(tp_cv/(fp_cv+tp_cv)*100,2)

f_score_cv = round((2*sens_call_cv*spec_call_cv)/(sens_call_cv+spec_call_cv)/100,2)

IC_lower_sens = round(as.numeric(sens_call_cv) - 1.96*sqrt(as.numeric(sens_call_cv)*(100-as.numeric(sens_call_cv))/N),digits=2)
IC_upper_sens = round(as.numeric(sens_call_cv) + 1.96*sqrt(as.numeric(sens_call_cv)*(100-as.numeric(sens_call_cv))/N),digits=2)
IC_sens = paste0("[",IC_lower_sens,",",IC_upper_sens,"]")

IC_lower_spec = round(as.numeric(spec_call_cv) - 1.96*sqrt(as.numeric(spec_call_cv)*(100-as.numeric(spec_call_cv))/(tp_cv+fp_cv)),digits=2)
IC_upper_spec = round(as.numeric(spec_call_cv) + 1.96*sqrt(as.numeric(spec_call_cv)*(100-as.numeric(spec_call_cv))/(tp_cv+fp_cv)),digits=2)
IC_spec = paste0("[",IC_lower_spec,",",IC_upper_spec,"]")

line = as.data.frame(t(c("","","LR (training 70%)",tp_cv,fp_cv,sens_call_cv,IC_sens,spec_call_cv,IC_spec,f_score_cv,NA,
                         NA,NA,NA,NA)))
colnames(line) <- c("SV_type","N","Caller","TP","FP","Recall","CI Recall","Precision","CI Precision","F1-Score","Genotype error",
                    "N 0/1","Genotype error 0/1","N 1/1","Genotype error 1/1")

metrics_callers = rbind(metrics_callers,line)



# detected by at least 1 caller

det_call = my_data3_70 %>% filter(callers_detected %in% 1:9)

callers1 = logical_sens_spec(det_call)

line_call1=as.data.frame(t(c("","",">=1 caller (70%)",callers1,NA,NA,NA,NA,NA)))
colnames(line_call1) <- c("SV_type","N","Caller","TP","FP","Recall","CI Recall","Precision","CI Precision","F1-Score","Genotype error",
                          "N 0/1","Genotype error 0/1","N 1/1","Genotype error 1/1")

metrics_callers = rbind(metrics_callers,line_call1)

# detected by at least 2 callers and 2 strategies

det_call = my_data3_70 %>% filter(callers_detected %in% 2:9 & strategy %in% 2:6 & reciprocity >=80)

callers2 = logical_sens_spec(det_call)

line_call2=as.data.frame(t(c("","",">=2 caller (70%) **",callers2,NA,NA,NA,NA,NA)))
colnames(line_call2) <- c("SV_type","N","Caller","TP","FP","Recall","CI Recall","Precision","CI Precision","F1-Score","Genotype error",
                          "N 0/1","Genotype error 0/1","N 1/1","Genotype error 1/1")

metrics_callers = rbind(metrics_callers,line_call2)


# detected by at least 3 callers and 2 strategies

det_call = my_data3_70 %>% filter(callers_detected %in% 3:5 & strategy %in% 2:6 & reciprocity >80)

callers3 = logical_sens_spec(det_call)

line_call3=as.data.frame(t(c("","",">=3 caller (70%) **",callers3,NA,NA,NA,NA,NA)))
colnames(line_call3) <- c("SV_type","N","Caller","TP","FP","Recall","CI Recall","Precision","CI Precision","F1-Score","Genotype error",
                          "N 0/1","Genotype error 0/1","N 1/1","Genotype error 1/1")

metrics_callers = rbind(metrics_callers,line_call3)


# detected by at least 4 callers and 2 strategies

det_call = my_data3_70 %>% filter(callers_detected %in% 4:5 & strategy %in% 2:6 & reciprocity >80)

callers4 = logical_sens_spec(det_call)

line_call4=as.data.frame(t(c("","",">=4 caller (70%) **",callers4,NA,NA,NA,NA,NA)))
colnames(line_call4) <- c("SV_type","N","Caller","TP","FP","Recall","CI Recall","Precision","CI Precision","F1-Score","Genotype error",
                          "N 0/1","Genotype error 0/1","N 1/1","Genotype error 1/1")

metrics_callers = rbind(metrics_callers,line_call4)

                
## Validation 30%

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

sens_call_cv = round(tp_cv/N*100,2)
spec_call_cv = round(tp_cv/(fp_cv+tp_cv)*100,2)

f_score_cv = round((2*sens_call_cv*spec_call_cv)/(sens_call_cv+spec_call_cv)/100,3)

IC_lower_sens = round(as.numeric(sens_call_cv) - 1.96*sqrt(as.numeric(sens_call_cv)*(100-as.numeric(sens_call_cv))/N),digits=2)
IC_upper_sens = round(as.numeric(sens_call_cv) + 1.96*sqrt(as.numeric(sens_call_cv)*(100-as.numeric(sens_call_cv))/N),digits=2)
IC_sens = paste0("[",IC_lower_sens,",",IC_upper_sens,"]")

IC_lower_spec = round(as.numeric(spec_call_cv) - 1.96*sqrt(as.numeric(spec_call_cv)*(100-as.numeric(spec_call_cv))/(tp_cv+fp_cv)),digits=2)
IC_upper_spec = round(as.numeric(spec_call_cv) + 1.96*sqrt(as.numeric(spec_call_cv)*(100-as.numeric(spec_call_cv))/(tp_cv+fp_cv)),digits=2)
IC_spec = paste0("[",IC_lower_spec,",",IC_upper_spec,"]")

line = as.data.frame(t(c("","","LR (Validation 30%)",tp_cv,fp_cv,sens_call_cv,IC_sens,spec_call_cv,IC_spec,f_score_cv,NA,
                         NA,NA,NA,NA)))
colnames(line) <- c("SV_type","N","Caller","TP","FP","Recall","CI Recall","Precision","CI Precision","F1-Score","Genotype error",
                    "N 0/1","Genotype error 0/1","N 1/1","Genotype error 1/1")

metrics_callers = rbind(metrics_callers,line)


my_data3_30 = my_data2_30 %>% filter(name_callers_detected != "svaba") #Svaba is eliminated because is not able to classify by type the SV
my_data3_30$callers_detected <- gsub("4-5-6-7", 4,my_data3_30$callers_detected )
my_data3_30$strategy <- gsub("3-4-5", 3,my_data3_30$strategy )


# detected by at least 1 caller

det_call = my_data3_30 %>% filter(callers_detected %in% 1:9)

callers1 = logical_sens_spec(det_call)

line_call1=as.data.frame(t(c("","",">=1 caller (30%)",callers1,NA,NA,NA,NA,NA)))
colnames(line_call1) <- c("SV_type","N","Caller","TP","FP","Recall","CI Recall","Precision","CI Precision","F1-Score","Genotype error",
                          "N 0/1","Genotype error 0/1","N 1/1","Genotype error 1/1")

metrics_callers = rbind(metrics_callers,line_call1)

# detected by at least 2 callers and 2 strategies

det_call = my_data3_30 %>% filter(callers_detected %in% 2:9 & strategy %in% 2:6 & reciprocity >=80)

callers2 = logical_sens_spec(det_call)

line_call2=as.data.frame(t(c("","",">=2 caller (30%) **",callers2,NA,NA,NA,NA,NA)))
colnames(line_call2) <- c("SV_type","N","Caller","TP","FP","Recall","CI Recall","Precision","CI Precision","F1-Score","Genotype error",
                          "N 0/1","Genotype error 0/1","N 1/1","Genotype error 1/1")

metrics_callers = rbind(metrics_callers,line_call2)


# detected by at least 3 callers and 2 strategies

det_call = my_data3_30 %>% filter(callers_detected %in% 3:5 & strategy %in% 2:6 & reciprocity >80)

callers3 = logical_sens_spec(det_call)

line_call3=as.data.frame(t(c("","",">=3 caller (30%) **",callers3,NA,NA,NA,NA,NA)))
colnames(line_call3) <- c("SV_type","N","Caller","TP","FP","Recall","CI Recall","Precision","CI Precision","F1-Score","Genotype error",
                          "N 0/1","Genotype error 0/1","N 1/1","Genotype error 1/1")

metrics_callers = rbind(metrics_callers,line_call3)


# detected by at least 4 callers and 2 strategies

det_call = my_data3_30 %>% filter(callers_detected %in% 4:5 & strategy %in% 2:6 & reciprocity >80)

callers4 = logical_sens_spec(det_call)

line_call4=as.data.frame(t(c("","",">=4 caller (30%) **",callers4,NA,NA,NA,NA,NA)))
colnames(line_call4) <- c("SV_type","N","Caller","TP","FP","Recall","CI Recall","Precision","CI Precision","F1-Score","Genotype error",
                          "N 0/1","Genotype error 0/1","N 1/1","Genotype error 1/1")

metrics_callers = rbind(metrics_callers,line_call4)

