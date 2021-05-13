##### Deletions ###########

library(data.table)
library(dplyr)

source("ext/functions.R")

setwd("1_LRM/DEL/")

# read golden ######

golden_del = fread("data/Insilico_DEL")

golden_del$V1 = gsub("chr","",golden_del$V1)

golden_geno_del = NULL

for(i in names(table(golden_del$V1))){
  
  aux = golden_del[golden_del$V1==i,]
  
  which_duplicated = aux$V2[which(duplicated(aux$V2))]
  
  golden_dup = unique(aux[aux$V2 %in% which_duplicated,])
  golden_no_dup = aux[!aux$V2 %in% which_duplicated,]
  
  golden_dup$genotype = "1/1"
  golden_no_dup$genotype = "0/1"
  
  aux_geno = rbind(golden_dup,golden_no_dup)
  
  golden_geno_del = rbind(golden_geno_del,aux_geno)
}


golden = golden_geno_del
colnames(golden) = c("chr","start_golden","length_golden","GT_golden")
golden$lower_golden = golden$start_golden-0
golden$upper_golden = golden$start_golden+0
golden$chr_pos_golden = do.call(paste0,list(golden$chr,"_",golden$start_golden))

golden = golden %>% arrange(chr,start_golden) %>% as.data.table()

nrow(golden)

# read vcfs files ########

delly = fread("data/Delly_DEL")
delly$V4 = as.numeric(abs(delly$V4))
colnames(delly) = c("chr","start_delly","end","length_delly","GT_delly")
delly$lower_delly = delly$start_delly-100
delly$upper_delly = delly$start_delly+100
delly$chr_pos_delly = do.call(paste0,list(delly$chr,"_",delly$start_delly))
delly = delly %>% arrange(chr,start_delly) %>% as.data.table()


lumpy = fread("data/Lumpy_DEL")
lumpy$V4 = as.numeric(abs(lumpy$V4))
colnames(lumpy) = c("chr","start_lumpy","end","length_lumpy","GT_lumpy")
lumpy$lower_lumpy = lumpy$start_lumpy-100
lumpy$upper_lumpy = lumpy$start_lumpy+100
lumpy$chr_pos_lumpy = do.call(paste0,list(lumpy$chr,"_",lumpy$start_lumpy))
lumpy = lumpy %>% arrange(chr,start_lumpy) %>% as.data.table()


pindel = fread("data/Pindel_DEL")
colnames(pindel) = c("chr","start_pindel","end","length_pindel","GT_pindel")
pindel$lower_pindel = pindel$start_pindel-10
pindel$upper_pindel = pindel$start_pindel+10
pindel$chr_pos_pindel = do.call(paste0,list(pindel$chr,"_",pindel$start_pindel))
pindel = pindel %>% arrange(chr,start_pindel) %>% as.data.table()

whamg = fread("data/Whamg_DEL",fill = T)
whamg$V4 = as.numeric(abs(whamg$V4))
colnames(whamg) = c("chr","start_whamg","end","length_whamg","GT_whamg")
whamg$lower_whamg = whamg$start_whamg-10
whamg$upper_whamg = whamg$start_whamg+10
whamg$chr_pos_whamg = do.call(paste0,list(whamg$chr,"_",whamg$start_whamg))
whamg = whamg %>% arrange(chr,start_whamg) %>% as.data.table()

svaba = fread("data/SVaBA_DEL")
colnames(svaba) = c("chr","start_svaba","chr2","pos2","GT_svaba")

# remove translocations and change positions

remove_variants = NULL

for(i in 1:nrow(svaba)){
  
  if(svaba$chr[i]==svaba$chr2[i] & svaba$start_svaba[i]>svaba$pos2[i]){
    
    start = svaba$start_svaba[i]
    end = svaba$pos2[i]
    
    svaba$start_svaba[i] = end
    svaba$pos2[i] = start
  }
  if(svaba$chr[i]!=svaba$chr2[i]){
    
    remove_variants = c(remove_variants,i)
    
  }
}

svaba = svaba[-remove_variants,]

svaba$length_svaba = svaba$pos2-svaba$start_svaba
svaba$lower_svaba = svaba$start_svaba-100
svaba$upper_svaba = svaba$start_svaba+100
svaba$chr_pos_svaba = do.call(paste0,list(svaba$chr,"_",svaba$start_svaba))
svaba = svaba %>% filter(length_svaba>149)
svaba = svaba %>% arrange(chr,start_svaba) %>% as.data.table()

manta = fread("data/Manta_DEL",fill = T)
manta$V4 = as.numeric(abs(manta$V4))
colnames(manta) = c("chr","start_manta","end","length_manta","GT_manta")
manta$lower_manta = manta$start_manta-50
manta$upper_manta = manta$start_manta+50
manta$chr_pos_manta = do.call(paste0,list(manta$chr,"_",manta$start_manta))
manta = manta %>% arrange(chr,start_manta) %>% as.data.table()

cnvnator = fread("data/CNVnator_DEL")
cnvnator$V4 = as.numeric(abs(cnvnator$V4))
colnames(cnvnator) = c("chr","start_cnvnator","end","length_cnvnator","GT_cnvnator")
cnvnator$lower_cnvnator = cnvnator$start_cnvnator-300
cnvnator$upper_cnvnator = cnvnator$start_cnvnator+300
cnvnator$chr_pos_cnvnator = do.call(paste0,list(cnvnator$chr,"_",cnvnator$start_cnvnator))
cnvnator = cnvnator %>% arrange(chr,start_cnvnator) %>% as.data.table()


## prepare BBDD ####

lumpy$ID = paste0("deletion_",1:nrow(lumpy))
delly$ID = "none"
pindel$ID = "none"
whamg$ID = "none"
manta$ID = "none"
cnvnator$ID = "none"
svaba$ID = "none"
golden$ID = "none"


my_data = merge_callers(delly,lumpy,callers_ref = c("lumpy"),
                        caller_to_merge = c("delly"),repro = 0.7) 

my_data = merge_callers(pindel,my_data,callers_ref = c("lumpy","delly"),
                        caller_to_merge = c("pindel"),repro = 0.7) 

my_data = merge_callers(whamg,my_data,callers_ref = c("lumpy","delly","pindel"),
                        caller_to_merge = c("whamg"),repro = 0.7) 

my_data = merge_callers(manta,my_data,callers_ref = c("lumpy","delly","pindel","whamg"),
                        caller_to_merge = c("manta"),repro = 0.7) 

my_data = merge_callers(cnvnator,my_data,callers_ref = c("lumpy","delly","pindel","whamg",
                                                         "manta"),
                        caller_to_merge = c("cnvnator"),repro = 0.7) 

my_data = merge_callers(svaba,my_data,callers_ref = c("lumpy","delly","pindel","whamg",
                                                      "manta","cnvnator"),
                        caller_to_merge = c("svaba"),repro = 0.7) 

my_data2 = merge_callers(golden,my_data,callers_ref = c("lumpy","delly","pindel","whamg",
                                                        "manta","cnvnator","svaba"),
                         caller_to_merge = c("golden"),golden = TRUE,repro = 0.7) 

dim(my_data)
dim(my_data2)

my_data2 = unique(my_data2)


## add number of callers detected #####

my_data2$callers_detected = n_callers_detected(my_data2,
                                               callers = c("lumpy","delly","pindel","whamg",
                                                           "manta","cnvnator","svaba"))

table(my_data2$callers_detected)


### consensus length #######

my_data2$length = consensus_length(my_data2,
                                   callers = c("lumpy","delly","pindel","whamg",
                                               "manta","svaba"))

my_data2$length[is.na(my_data2$length)] = my_data2$length_cnvnator[is.na(my_data2$length)] 

summary(my_data2$length)

aux = cut(my_data2$length,
          breaks = c(150,500,1000,2000,3000,10000,Inf),right = F)

my_data2$length_stretch = aux


table(my_data2$length_stretch,my_data2$PASS)


## add reciprocity ###

my_data2$reciprocity = reprocicity(my_data2,
                                   callers = c("lumpy","delly","pindel","whamg",
                                               "manta","cnvnator","svaba"))
summary(my_data2$reciprocity)


## add strategy ####

strategies = fread("ext/strategies.csv")

my_data2$strategy = strategy(my_data2,
                             callers = c("lumpy","delly","pindel","whamg",
                                         "manta","cnvnator","svaba"),
                             strategies)

table(my_data2$strategy)


### colineality #####

table(my_data2$callers_detected,my_data2$strategy)
chisq.test(table(my_data2$callers_detected,my_data2$strategy))


## prepare model choose de best model####

my_data2$PASS = "YES"

my_data2$PASS[is.na(my_data2$length_golden)] = "NO"

table(my_data2$PASS)

dim(golden)

my_data2 = my_data2[-which(duplicated(my_data2$ID)),]

table(my_data2$PASS,my_data2$callers_detected)
table(my_data2$PASS,my_data2$strategy)


my_data2$GT_manta[is.na(my_data2$GT_manta)] = "9/9"
my_data2$GT_whamg[is.na(my_data2$GT_whamg)] = "9/9"
my_data2$GT_delly[is.na(my_data2$GT_delly)] = "9/9"
my_data2$GT_lumpy[is.na(my_data2$GT_lumpy)] = "9/9"
my_data2$GT_pindel[is.na(my_data2$GT_pindel)] = "9/9"
my_data2$GT_svaba[is.na(my_data2$GT_svaba)] = "9/9"
my_data2$GT_cnvnator[is.na(my_data2$GT_cnvnator)] = "9/9"

table(my_data2$GT_manta,useNA = "always")
table(my_data2$GT_whamg,useNA = "always")
table(my_data2$GT_delly,useNA = "always")
table(my_data2$GT_lumpy,useNA = "always")
table(my_data2$GT_svaba,useNA = "always")
table(my_data2$GT_cnvnator,useNA = "always")
table(my_data2$GT_pindel,useNA = "always")

my_data2$PASS = ifelse(my_data2$PASS=="YES",0,1)

my_model = glm(PASS ~ factor(GT_manta)+
                 factor(GT_whamg)+
                 factor(GT_delly)+
                 factor(GT_lumpy)+
                 factor(GT_pindel)+
                 factor(GT_cnvnator)+
                 factor(GT_svaba) + factor(length_stretch) + factor(strategy)+
                 reciprocity,
               data=my_data2,family = "binomial",
               control = list(maxit = 100))

my_data2$PASS = ifelse(my_data2$PASS==0,"YES","NO")

summary(my_model)

my_pred = predict(my_model,my_data2,type="response")

my_pred = ifelse(my_pred>0.5,1,0)

my_pred = ifelse(my_pred==1,"NO PASS","PASS")

table(my_data2$PASS)
table(my_pred,my_data2$PASS)
table(my_pred)


# without CHR and BP
#sensitivity
model_sens = (table(my_pred,my_data2$PASS)[2,2]/nrow(golden))*100
#specificity
model_spec = (table(my_pred,my_data2$PASS)[2,2])/(table(my_pred,my_data2$PASS)[2,2]+
                                                    table(my_pred,my_data2$PASS)[2,1])*100

f_score = ((2*model_sens*model_spec)/(model_sens + model_spec))/100

my_data2$prediction = my_pred

my_tab = table(my_data2$callers_detected,my_data2$PASS)
my_tab = cbind(my_tab,table(my_data2$callers_detected,my_data2$prediction))

prop.table(table(my_data2$callers_detected,my_data2$PASS),margin=1)


## genotype error in the model ######

my_data2$GT = NA

data_call = my_data2 %>% filter(prediction=="PASS") %>% filter(PASS=="YES") %>% 
  dplyr::select(GT_whamg,
                GT_svaba,
                GT_lumpy,
                GT_pindel,
                GT_delly,
                GT_manta,
                GT_cnvnator,
                GT,GT_golden,PASS,prediction)


data_call$GT_delly[data_call$GT_delly=="9/9"] = NA
data_call$GT_pindel[data_call$GT_pindel=="9/9"] = NA
data_call$GT_svaba[data_call$GT_svaba=="9/9"] = NA
data_call$GT_whamg[data_call$GT_whamg=="9/9"] = NA
data_call$GT_lumpy[data_call$GT_lumpy=="9/9"] = NA
data_call$GT_manta[data_call$GT_manta=="9/9"] = NA
data_call$GT_cnvnator[data_call$GT_cnvnator =="9/9"] = NA

table(data_call$GT_whamg)
table(data_call$GT_svaba)
table(data_call$GT_lumpy)
table(data_call$GT_pindel)
table(data_call$GT_delly)
table(data_call$GT_manta)
table(data_call$GT_cnvnator)


data_call = as.matrix(data_call)

for(i in 1:nrow(data_call)){
  
  my_tab = table(as.character(data_call[i,c(1:7)]))
  
  if(length(which(my_tab==max(my_tab)))==1){
    data_call[i,8] = names(which.max(my_tab))
  }
  
  if(length(which(my_tab==max(my_tab)))>1){
    
    data_call[i,8] = NA
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


### choose best model ####

metrics_callers$Caller = as.character(metrics_callers$Caller)

metrics_callers = metrics_callers[1:8,]

my_data2$PASS = ifelse(my_data2$PASS=="YES",0,1)

full.model =  glm(PASS ~ factor(GT_manta)+
                    factor(GT_whamg)+
                    factor(GT_delly)+
                    factor(GT_lumpy)+
                    factor(GT_pindel)+
                    factor(GT_cnvnator)+
                    factor(GT_svaba) + factor(length_stretch) + 
                    factor(strategy) +
                    reciprocity,
                  data=my_data2,family = "binomial")

base.model = glm(PASS ~ factor(GT_manta)+
                   factor(GT_whamg)+
                   factor(GT_delly)+
                   factor(GT_lumpy)+
                   factor(GT_pindel)+
                   factor(GT_cnvnator)+
                   factor(GT_svaba),
                 data=my_data2,family = "binomial")

step(full.model,
     scope = list(lower = formula(base.model),
                  upper = formula(full.model)),
     direction = "backward")


best.model = glm(formula = PASS ~ factor(GT_manta) + factor(GT_whamg) + factor(GT_delly) + 
                   factor(GT_lumpy) + factor(GT_pindel) + factor(GT_cnvnator) + 
                   factor(GT_svaba) + factor(length_stretch) + factor(strategy), 
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
model_sens = round(model_tp/N*100,3)
#specificity
model_spec = round(model_tp/(model_fp+model_tp)*100,3)

model_tp;model_fp;model_sens;model_spec

f_score = round(((2*model_sens*model_spec)/(model_sens + model_spec))/100,3)

my_data2$prediction = my_pred


# detected 1 caller

det_call = my_data2 %>% filter(callers_detected %in% 1:9)

callers1 = logical_sens_spec(det_call)

# detected 2 callers and 2 strategies

det_call = my_data2 %>% filter(callers_detected %in% 2:9 & strategy %in% 2:6 & reciprocity >80)

callers2 = logical_sens_spec(det_call)

# detected 3 callers and 2 strategies

det_call = my_data2 %>% filter(callers_detected %in% 3:9 & strategy %in% 2:6 & reciprocity >80)

callers3 = logical_sens_spec(det_call)

# detected 4 callers and 2 strategies

det_call = my_data2 %>% filter(callers_detected %in% 4:9 & strategy %in% 2:6 & reciprocity >80)

callers4 = logical_sens_spec(det_call)


metrics_callers = rbind(metrics_callers,
                        c("","","LR",model_tp,model_fp,model_sens,model_spec,f_score,
                          model_g_error))

metrics_callers = rbind(metrics_callers,
                        c("","",">=1 caller",callers1,""))

metrics_callers = rbind(metrics_callers,
                        c("","",">=2 caller **",callers2,""))

metrics_callers = rbind(metrics_callers,
                        c("","",">=3 caller **",callers3,""))

metrics_callers = rbind(metrics_callers,
                        c("","",">=4 caller **",callers4,""))



### CV ######


my_data2$GT_manta[my_data2$GT_manta=="0/1"] = "0/1-1/1"
my_data2$GT_manta[my_data2$GT_manta=="1/1"] = "0/1-1/1"

table(my_data2$GT_manta,my_data2$PASS)

my_data2$GT_cnvnator[my_data2$GT_cnvnator=="0/1"] = "0/1-1/1"
my_data2$GT_cnvnator[my_data2$GT_cnvnator=="1/1"] = "0/1-1/1"

my_data2$GT_delly[my_data2$GT_delly=="0/1"] = "0/1-1/1"
my_data2$GT_delly[my_data2$GT_delly=="1/1"] = "0/1-1/1"

my_data2$GT_lumpy[my_data2$GT_lumpy=="0/1"] = "0/1-1/1"
my_data2$GT_lumpy[my_data2$GT_lumpy=="1/1"] = "0/1-1/1"

my_data2$GT_pindel[my_data2$GT_pindel=="0/1"] = "0/1-1/1"
my_data2$GT_pindel[my_data2$GT_pindel=="1/1"] = "0/1-1/1"

my_data2$GT_svaba[my_data2$GT_svaba=="0/1"] = "0/1-1/1"
my_data2$GT_svaba[my_data2$GT_svaba=="1/1"] = "0/1-1/1"

my_data2$GT_whamg[my_data2$GT_whamg=="0/1"] = "0/1-1/1"
my_data2$GT_whamg[my_data2$GT_whamg=="1/1"] = "0/1-1/1"

my_data2$strategy[my_data2$strategy=="4"] = "4-5"
my_data2$strategy[my_data2$strategy=="5"] = "4-5"

table(my_data2$strategy,my_data2$PASS)


library(caret)
library(e1071)

ctrl <- trainControl(method = "repeatedcv", number = 10, 
                     savePredictions = TRUE)

## 70-30

0.7*nrow(my_data2)

set.seed(2589)

training = sample(1:nrow(my_data2),2392)

table(my_data2$PASS,useNA = "always")

my_data2$PASS = as.factor(my_data2$PASS)

my_data2_70 = my_data2[training,]

table(my_data2_70$PASS,useNA = "always")

mod_fit <- train(factor(PASS) ~ factor(GT_manta) + factor(GT_whamg) + factor(GT_delly) + 
                   factor(GT_lumpy) + factor(GT_pindel) + factor(GT_cnvnator) + 
                   factor(GT_svaba) + factor(length_stretch) + factor(strategy), 
                 data=my_data2_70, method="glm", family="binomial",
                 trControl = ctrl)

summary(mod_fit)


mod_fit_l_num <- train(factor(PASS) ~ factor(GT_manta) + factor(GT_whamg) + factor(GT_delly) + 
                         factor(GT_lumpy) + factor(GT_pindel) + factor(GT_cnvnator) + 
                         factor(GT_svaba) + as.numeric(length) + factor(strategy), 
                       data=my_data2_70, method="glm", family="binomial",
                       trControl = ctrl)

mod_fit_no_l <- train(factor(PASS) ~ factor(GT_manta) + factor(GT_whamg) + factor(GT_delly) + 
                        factor(GT_lumpy) + factor(GT_pindel) + factor(GT_cnvnator) + 
                        factor(GT_svaba) + factor(strategy), 
                      data=my_data2_70, method="glm", family="binomial",
                      trControl = ctrl)


# 0.5

my_pred = predict(mod_fit, newdata=my_data2_70,type = "prob")

my_pred = my_pred[,1]
my_pred = ifelse(my_pred<0.5,"PASS","NO PASS")

N = as.numeric(table(my_data2_70$PASS)[2])

tp_cv = as.numeric(table(my_pred,my_data2_70$PASS)[2,2])
fp_cv = as.numeric(table(my_pred)[2])-tp_cv

sens_call_cv = round(tp_cv/N*100,3)
spec_call_cv = round(tp_cv/(fp_cv+tp_cv)*100,3)

f_score_cv = round((2*sens_call_cv*spec_call_cv)/(sens_call_cv+spec_call_cv)/100,3)
