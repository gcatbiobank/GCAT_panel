library(data.table)
library(dplyr)

source("ext/functions.R")

setwd("1_LRM/INV/")

# read golden ######

golden_inv = fread("data/Insilico_INV",fill = T)

golden_inv = golden_inv[,-c(5)]

golden_inv$V1 = gsub("chr","",golden_inv$V1)

golden_geno_inv = NULL

for(i in names(table(golden_inv$V1))){
  
  aux = golden_inv[golden_inv$V1==i,]
  
  which_duplicated = aux$V2[which(duplicated(aux$V2))]
  
  golden_dup = unique(aux[aux$V2 %in% which_duplicated,])
  golden_no_dup = aux[!aux$V2 %in% which_duplicated,]
  
  golden_dup$genotype = "1/1"
  golden_no_dup$genotype = "0/1"
  
  aux_geno = rbind(golden_dup,golden_no_dup)
  
  golden_geno_inv = rbind(golden_geno_inv,aux_geno)
}


golden = golden_geno_inv
colnames(golden) = c("chr","start_golden","length_golden","end_golden","GT_golden")

# change start and end

golden$length_aux = golden$end_golden-golden$start_golden

for(i in 1:nrow(golden)){
  
  if(golden$length_aux[i]<0){
    
    aux = golden$start_golden[i]
    aux2 = golden$end_golden[i]
    
    golden$start_golden[i] = aux2
    golden$end_golden[i] = aux
    
  }
  
  if(golden$length_aux[i]==0){
    
    golden$end_golden[i] =  golden$start_golden[i] + golden$length_golden[i]
    
  }
  
}

golden$length_aux = golden$end_golden-golden$start_golden
length(which(golden$length<=0))

golden$lower_golden_bp1 = golden$start_golden-0
golden$upper_golden_bp1 = golden$start_golden+0
golden$lower_golden_bp2 = golden$end_golden-0
golden$upper_golden_bp2 = golden$end_golden+0
golden$chr_pos_golden = do.call(paste0,list(golden$chr,"_",golden$start_golden))

golden = golden %>% as.data.table()


## read vcf ######

delly = fread("data/Delly_INV")
delly$V4 = as.numeric(abs(delly$V4))
colnames(delly) = c("chr","start_delly","end_delly","length_delly","GT_delly")
delly = unify_breakpoints(delly,"delly")
colnames(delly)[2:5] = c("chr","start_delly","end_delly","GT_delly")
delly$length_delly = delly$end-delly$start
delly$lower_delly_bp1 = delly$start_delly-20
delly$upper_delly_bp1 = delly$start_delly+20
delly$lower_delly_bp2 = delly$end_delly-20
delly$upper_delly_bp2 = delly$end_delly+20

delly = delly %>% as.data.table()

delly = unique(delly)
delly$ID = NULL

lumpy = fread("data/Lumpy_INV")
lumpy$V4 = as.numeric(abs(lumpy$V4))
colnames(lumpy) = c("chr","start_lumpy","end_lumpy","length_lumpy","GT_lumpy")
lumpy = unify_breakpoints(lumpy,"lumpy")
colnames(lumpy)[2:5] = c("chr","start_lumpy","end_lumpy","GT_lumpy")
lumpy$length_lumpy = lumpy$end_lumpy-lumpy$start_lumpy
lumpy$lower_lumpy_bp1 = lumpy$start_lumpy-50
lumpy$upper_lumpy_bp1 = lumpy$start_lumpy+50
lumpy$lower_lumpy_bp2 = lumpy$end_lumpy-50
lumpy$upper_lumpy_bp2 = lumpy$end_lumpy+50

lumpy = lumpy %>% as.data.table()

lumpy = unique(lumpy)
lumpy$ID = NULL


pindel = fread("data/Pindel_INV")
colnames(pindel) = c("chr","start_pindel","end_pindel","length_pindel","GT_pindel")
pindel = unify_breakpoints(pindel,"pindel")
colnames(pindel)[2:5] = c("chr","start_pindel","end_pindel","GT_pindel")
pindel$length_pindel = pindel$end_pindel-pindel$start_pindel
pindel$lower_pindel_bp1 = pindel$start_pindel-10
pindel$upper_pindel_bp1 = pindel$start_pindel+10
pindel$lower_pindel_bp2 = pindel$end_pindel-10
pindel$upper_pindel_bp2 = pindel$end_pindel+10

dim(pindel)
pindel = unique(pindel)
pindel = pindel %>% as.data.table()
pindel$ID = NULL


whamg = fread("data/Whamg_INV",fill = T)
colnames(whamg) = c("chr","start_whamg","end_whamg","length_whamg","GT_whamg")
whamg = unify_breakpoints(whamg,"whamg")
colnames(whamg)[2:5] = c("chr","start_whamg","end_whamg","GT_whamg")
whamg$length_whamg = whamg$end_whamg-whamg$start_whamg
whamg$lower_whamg_bp1 = whamg$start_whamg-10
whamg$upper_whamg_bp1 = whamg$start_whamg+10
whamg$lower_whamg_bp2 = whamg$end_whamg-10
whamg$upper_whamg_bp2 = whamg$end_whamg+10

whamg = unique(whamg)
whamg = whamg %>% as.data.table()
whamg$ID = NULL

manta = fread("data/Manta_INV",fill = T)
manta$V4 = as.numeric(abs(manta$V4))
colnames(manta) = c("chr","start_manta","end_manta","length_manta","GT_manta")
manta = unify_breakpoints(manta,"manta")
colnames(manta)[2:5] = c("chr","start_manta","end_manta","GT_manta")
manta$length_manta = manta$end_manta-manta$start_manta
manta$lower_manta_bp1 = manta$start_manta-10
manta$upper_manta_bp1 = manta$start_manta+10
manta$lower_manta_bp2 = manta$end_manta-10
manta$upper_manta_bp2 = manta$end_manta+10


svaba = fread("data/SVABA_INV")
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

svaba = svaba[,c(1:3,6,4)]
colnames(svaba) = c("chr","start_svaba","end_svaba","length_svaba","GT_svaba")
svaba = unique(svaba)
svaba = unify_breakpoints(svaba,"svaba")
colnames(svaba)[2:5] = c("chr","start_svaba","end_svaba","GT_svaba")
svaba$length_svaba = svaba$end_svaba-svaba$start_svaba
svaba$lower_svaba_bp1 = svaba$start_svaba-10
svaba$upper_svaba_bp1 = svaba$start_svaba+10
svaba$lower_svaba_bp2 = svaba$end_svaba-10
svaba$upper_svaba_bp2 = svaba$end_svaba+10

svaba = unique(svaba)
svaba = as.data.table(svaba)
svaba$ID = NULL

### sensitivity-specifitiy-geno error for each caller ######

# lumpy

sens_lumpy = sensitivity_precision_geno_error_inv(lumpy,golden,"inversion","lumpy")

# pindel

sens_pindel = sensitivity_precision_geno_error_inv(pindel,golden,"inversion","pindel")

# manta

sens_manta = sensitivity_precision_geno_error_inv(manta,golden,"inversion","manta")

# whamg

sens_whamg = sensitivity_precision_geno_error_inv(whamg,golden,"inversion","whamg")

# delly

sens_delly = sensitivity_precision_geno_error_inv(delly,golden,"inversion","delly")

# svaba

sens_svaba = sensitivity_precision_geno_error_inv(svaba,golden,"inversion","svaba")


metrics_callers = matrix(as.numeric(c(t(unlist(sens_delly))[c(3,4,1,2,5,6,7,9,8,10)],
                                      t(unlist(sens_lumpy))[c(3,4,1,2,5,6,7,9,8,10)],
                                      t(unlist(sens_manta))[c(3,4,1,2,5,6,7,9,8,10)],
                                      t(unlist(sens_pindel))[c(3,4,1,2,5,6,7,9,8,10)],
                                      t(unlist(sens_svaba))[c(3,4,1,2,5,6,7,9,8,10)],
                                      t(unlist(sens_whamg))[c(3,4,1,2,5,6,7,9,8,10)])),byrow=T,ncol = 10)

metrics_callers = cbind(c("Delly","Lumpy","Manta","Pindel",
                          "Svaba","Whamg"),
                        round(metrics_callers,2))

colnames(metrics_callers) = c("Caller","TP","FP","Sensitivity","Precision","F1-Score","Genotype error",
                              "N 0/1","Genotype error 0/1","N 1/1","Genotype error 1/1")

metrics_callers = metrics_callers %>% as.data.frame() %>% arrange(desc(`F1-Score`))

metrics_callers = cbind(N=c(nrow(golden),rep("",nrow(metrics_callers)-1)),
                        metrics_callers)

metrics_callers = cbind(SV_type=c("Inversions",rep("",nrow(metrics_callers)-1)),
                        metrics_callers)

metrics_callers[,1:3] = apply(metrics_callers[,1:3], 2, function(x) as.character(x));
metrics_callers[,4:9] = apply(metrics_callers[,4:9], 2, function(x) as.numeric(x));


## prepare BBDD ####

type = "inversion"

delly$ID = "none"

delly$ID = paste0(type,1:nrow(delly))

lumpy$ID = "none"
manta$ID = "none"
whamg$ID = "none"
pindel$ID = "none"
svaba$ID = "none"
golden$ID = "none"


my_data = merge_callers_inv(lumpy,delly,callers_ref = c("delly"),
                        caller_to_merge = c("lumpy"),
                        type="inversion") 

my_data = merge_callers_inv(manta,my_data,callers_ref = c("delly","lumpy"),
                           caller_to_merge = c("manta"),
                           type="inversion") 

my_data = merge_callers_inv(whamg,my_data,callers_ref = c("delly","lumpy","manta"),
                           caller_to_merge = c("whamg"),
                           type="inversion") 

my_data = merge_callers_inv(pindel,my_data,callers_ref = c("delly","lumpy","manta","whamg"),
                           caller_to_merge = c("pindel"),
                           type="inversion") 

my_data = merge_callers_inv(svaba,my_data,callers_ref = c("delly","lumpy","manta","whamg","pindel"),
                           caller_to_merge = c("svaba"),
                           type="inversion") 

my_data2 = merge_callers_inv(golden,my_data,callers_ref = c("delly","lumpy","manta","whamg","pindel","svaba"),
                            caller_to_merge = c("golden"),
                            type="inversion") 

dim(my_data2)
dim(my_data)

my_data2 = unique(my_data2)

my_data2 = as.data.table(my_data2)

### start, end, length, 

callers_ref = c("delly","lumpy","manta","whamg","pindel","svaba")

callers_ref = c(callers_ref,"golden")

colnames_bbdd = as.vector(outer(c("start_","end_","length_",
                             "GT_"),callers_ref,paste0))

final_data = my_data2[,c("ID","chr",colnames_bbdd),with=F]

final_data = final_data %>% group_by(ID) %>% 
              dplyr::summarise(chr=names(sort(table(chr,useNA = "always"),decreasing = T))[1],
                   start = names(sort(table(c(start_pindel,start_manta,start_whamg,start_lumpy,start_delly,start_svaba,start_golden)),decreasing = T))[1],
                   end = names(sort(table(c(end_pindel,end_manta,end_whamg,end_lumpy,end_delly,end_svaba,end_golden)),decreasing = T))[1],                   
                   GT_whamg = names(sort(table(GT_whamg,useNA = "always"),decreasing = T))[1],
                   GT_pindel = names(sort(table(GT_pindel,useNA = "always"),decreasing = T))[1],
                   GT_delly = names(sort(table(GT_delly,useNA = "always"),decreasing = T))[1],
                   GT_lumpy = names(sort(table(GT_lumpy,useNA = "always"),decreasing = T))[1],
                   GT_manta = names(sort(table(GT_manta,useNA = "always"),decreasing = T))[1],
                   GT_svaba = names(sort(table(GT_svaba,useNA = "always"),decreasing = T))[1],
                   GT_golden = names(sort(table(GT_golden,useNA = "always"),decreasing = T))[1],
                   length_golden = names(sort(table(length_golden,useNA = "always"),decreasing = T))[1],
                   max_length = max(c(length_pindel,length_manta,length_whamg,length_lumpy,length_delly,length_svaba),na.rm = T),
                   min_length = min(c(length_pindel,length_manta,length_whamg,length_lumpy,length_delly,length_svaba),na.rm = T),
                   length_whamg = names(sort(table(length_whamg,useNA = "always"),decreasing = T))[1],
                   length_pindel = names(sort(table(length_pindel,useNA = "always"),decreasing = T))[1],
                   length_delly = names(sort(table(length_delly,useNA = "always"),decreasing = T))[1],
                   length_lumpy = names(sort(table(length_lumpy,useNA = "always"),decreasing = T))[1],
                   length_manta = names(sort(table(length_manta,useNA = "always"),decreasing = T))[1],
                   length_svaba = names(sort(table(length_svaba,useNA = "always"),decreasing = T))[1],
                   start_whamg = names(sort(table(start_whamg,useNA = "always"),decreasing = T))[1],
                   start_pindel = names(sort(table(start_pindel,useNA = "always"),decreasing = T))[1],
                   start_delly = names(sort(table(start_delly,useNA = "always"),decreasing = T))[1],
                   start_lumpy = names(sort(table(start_lumpy,useNA = "always"),decreasing = T))[1],
                   start_manta = names(sort(table(start_manta,useNA = "always"),decreasing = T))[1],
                   start_svaba = names(sort(table(start_svaba,useNA = "always"),decreasing = T))[1],
                   start_golden = names(sort(table(start_golden,useNA = "always"),decreasing = T))[1])
                     

final_data$length = abs(as.numeric(final_data$end)-as.numeric(final_data$start))

summary(final_data$length)

final_data = final_data %>% filter(length>30) 

final_data$reciprocity = as.numeric(final_data$min_length)/as.numeric(final_data$max_length)
final_data$length_golden = as.numeric(final_data$length_golden)

dim(final_data)

length(unique(my_data2$ID))

final_data = as.data.table(final_data)

### length stretch #####

aux = cut(final_data$length,
          breaks = c(30,150,500,1000,2000,3000,Inf),right = F)

table(aux)

final_data$length_stretch = aux


## add number of callers detected #####

final_data$callers_detected = n_callers_detected(final_data,
                                                 callers = c("whamg","lumpy","manta",
                                                             "pindel","delly","svaba"))

final_data$callers_detected_ok = n_callers_detected(final_data,
                                                 callers = c("whamg","lumpy","manta",
                                                             "pindel","delly","svaba"))


table(final_data$callers_detected)

final_data$callers_detected[final_data$callers_detected=="3"] = "3-4-5-6"
final_data$callers_detected[final_data$callers_detected=="4"] = "3-4-5-6"
final_data$callers_detected[final_data$callers_detected=="5"] = "3-4-5-6"
final_data$callers_detected[final_data$callers_detected=="6"] = "3-4-5-6"

table(final_data$callers_detected)
table(final_data$callers_detected_ok)


## add strategy ####

final_data = final_data %>% as.data.table()

strategies = fread("ext/strategies.csv")

final_data$strategy = 0

final_data$strategy = strategy(final_data,
                               callers = c("whamg","lumpy","manta",
                                           "pindel","delly","svaba"),
                               strategies)

final_data$strategy_ok = strategy(final_data,
                               callers = c("whamg","lumpy","manta",
                                           "pindel","delly","svaba"),
                               strategies)


table(final_data$strategy)

final_data$strategy[final_data$strategy=="3"] = "3-4"
final_data$strategy[final_data$strategy=="4"] = "3-4"

table(final_data$strategy)
table(final_data$strategy_ok)


### colineality #####

table(final_data$callers_detected,final_data$strategy)
chisq.test(table(final_data$callers_detected,final_data$strategy))


## prepare model ####

final_data$PASS = "YES"

final_data$PASS[is.na(final_data$GT_golden)] = "NO"

final_data = final_data %>% filter(callers_detected!=0)

table(final_data$PASS)
table(final_data$callers_detected,final_data$PASS)
table(final_data$strategy,final_data$PASS)
table(final_data$length_stretch,final_data$PASS)

final_data$GT_manta[is.na(final_data$GT_manta)] = "9/9"
final_data$GT_whamg[is.na(final_data$GT_whamg)] = "9/9"
final_data$GT_delly[is.na(final_data$GT_delly)] = "9/9"
final_data$GT_lumpy[is.na(final_data$GT_lumpy)] = "9/9"
final_data$GT_pindel[is.na(final_data$GT_pindel)] = "9/9"
final_data$GT_svaba[is.na(final_data$GT_svaba)] = "9/9"

final_data$GT_manta2 = final_data$GT_manta
final_data$GT_whamg2 = final_data$GT_whamg
final_data$GT_delly2 = final_data$GT_delly
final_data$GT_lumpy2 = final_data$GT_lumpy
final_data$GT_pindel2 = final_data$GT_pindel
final_data$GT_svaba2 = final_data$GT_svaba


final_data$GT_manta[final_data$GT_manta=="0/1"] = "0/1-1/1"
final_data$GT_manta[final_data$GT_manta=="1/1"] = "0/1-1/1"
final_data$GT_whamg[final_data$GT_whamg=="0/1"] = "0/1-1/1"
final_data$GT_whamg[final_data$GT_whamg=="1/1"] = "0/1-1/1"
final_data$GT_delly[final_data$GT_delly=="0/1"] = "0/1-1/1"
final_data$GT_delly[final_data$GT_delly=="1/1"] = "0/1-1/1"
final_data$GT_lumpy[final_data$GT_lumpy=="0/1"] = "0/1-1/1"
final_data$GT_lumpy[final_data$GT_lumpy=="1/1"] = "0/1-1/1"
final_data$GT_pindel[final_data$GT_pindel=="0/1"] = "0/1-1/1"
final_data$GT_pindel[final_data$GT_pindel=="1/1"] = "0/1-1/1"
final_data$GT_svaba[final_data$GT_svaba=="0/1"] = "0/1-1/1"
final_data$GT_svaba[final_data$GT_svaba=="1/1"] = "0/1-1/1"


table(final_data$GT_manta,useNA = "always")
table(final_data$GT_whamg,useNA = "always")
table(final_data$GT_delly,useNA = "always")
table(final_data$GT_lumpy,useNA = "always")
table(final_data$GT_svaba,useNA = "always")
table(final_data$GT_pindel,useNA = "always")



final_data$PASS = ifelse(final_data$PASS=="YES",0,1)

my_model =  glm(PASS ~ factor(GT_manta)+
                  factor(GT_whamg)+
                  factor(GT_delly)+
                  factor(GT_lumpy)+
                  factor(GT_pindel)+
                  factor(GT_svaba) + factor(length_stretch),
                data=final_data,family = "binomial",
                control = list(maxit = 200))

final_data$PASS = ifelse(final_data$PASS==0,"YES","NO")

summary(my_model)

my_pred = predict(my_model,final_data,type="response")

my_pred = ifelse(my_pred>0.5,1,0)

my_pred = ifelse(my_pred==1,"NO PASS","PASS")

table(final_data$PASS)
table(my_pred,final_data$PASS)
table(my_pred)

length(unique(golden$chr_pos))

# without CHR and BP
#sensitivity
model_sens = (table(my_pred,final_data$PASS)[2,2]/nrow(golden))*100
#specificity
model_spec = (table(my_pred,final_data$PASS)[2,2])/(table(my_pred,final_data$PASS)[2,2]+
                                                    table(my_pred,final_data$PASS)[2,1])*100

f_score = ((2*model_sens*model_spec)/(model_sens + model_spec))/100

final_data$prediction = my_pred


### name callers detected and prediction #####

final_data$name_callers_detected = name_callers_detected(as.data.table(final_data), 
                                                         callers = c("whamg","lumpy","manta",
                                                                     "pindel","delly","svaba"))

table(final_data$name_callers_detected,final_data$prediction) 


### genotype error in the model ######

# 1) Geno error: most common genotype

data_call = final_data %>% filter(prediction=="PASS") %>% filter(PASS=="YES") %>% 
  dplyr::select(GT_whamg2,
                GT_lumpy2,
                GT_pindel2,
                GT_delly2,
                GT_manta2,
                GT_golden,PASS,prediction)

data_call$GT_manta = data_call$GT_manta2
data_call$GT_whamg = data_call$GT_whamg2
data_call$GT_delly = data_call$GT_delly2
data_call$GT_lumpy = data_call$GT_lumpy2
data_call$GT_pindel = data_call$GT_pindel2
data_call$GT_svaba = data_call$GT_svaba2

data_call$GT_delly[data_call$GT_delly=="9/9"] = NA
data_call$GT_pindel[data_call$GT_pindel=="9/9"] = NA
data_call$GT_whamg[data_call$GT_whamg=="9/9"] = NA
data_call$GT_lumpy[data_call$GT_lumpy=="9/9"] = NA
data_call$GT_manta[data_call$GT_manta=="9/9"] = NA

table(data_call$GT_whamg)
table(data_call$GT_lumpy)
table(data_call$GT_pindel)
table(data_call$GT_delly)
table(data_call$GT_manta)


data_call$GT = NA

data_call = as.matrix(data_call)

for(i in 1:nrow(data_call)){
  
  my_tab = table(as.character(data_call[i,c(1:5)]))
  
  if(length(which(my_tab==max(my_tab)))==1){
    data_call[i,14] = names(which.max(my_tab))
  }
  
  if(length(which(my_tab==max(my_tab)))>1){
    
    data_call[i,14] = NA
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


# 2)  Order: Lumpy-Pindel-Whamg-Delly-Manta 

data_call = final_data %>% filter(prediction=="PASS") %>% filter(PASS=="YES") %>% 
  dplyr::select(GT_whamg2,
                GT_lumpy2,
                GT_pindel2,
                GT_delly2,
                GT_manta2,
                GT_golden,PASS,prediction)

data_call$GT_manta = data_call$GT_manta2
data_call$GT_whamg = data_call$GT_whamg2
data_call$GT_delly = data_call$GT_delly2
data_call$GT_lumpy = data_call$GT_lumpy2
data_call$GT_pindel = data_call$GT_pindel2
data_call$GT_svaba = data_call$GT_svaba2

data_call$GT_delly[data_call$GT_delly=="9/9"] = NA
data_call$GT_pindel[data_call$GT_pindel=="9/9"] = NA
data_call$GT_whamg[data_call$GT_whamg=="9/9"] = NA
data_call$GT_lumpy[data_call$GT_lumpy=="9/9"] = NA
data_call$GT_manta[data_call$GT_manta=="9/9"] = NA

table(data_call$GT_whamg)
table(data_call$GT_lumpy)
table(data_call$GT_pindel)
table(data_call$GT_delly)
table(data_call$GT_manta)

caller = c("manta","whamg","delly","lumpy","pindel")

data_call$GT = NA

data_call = as.matrix(data_call)

i=1

for(i in 1:nrow(data_call)){
  
  names_call = caller[which(!is.na(data_call[i,c(9:13)]))]
  
  if(sum(names_call %in% "lumpy")==1){
    
    data_call[i,14] = data_call[i,12]
  }
  
  
  if(sum(names_call %in% "lumpy")==0 & sum(names_call %in% "pindel")==1){
    
    data_call[i,14] = data_call[i,13]
  }
  
  if(sum(names_call %in% "lumpy")==0 & sum(names_call %in% "pindel")==0 & sum(names_call %in% "whamg")==1){
    
    data_call[i,14] = data_call[i,10]
  }
  
  if(sum(names_call %in% "lumpy")==0 & sum(names_call %in% "pindel")==0 & sum(names_call %in% "whamg")==0 &
     sum(names_call %in% "delly")==1){
    
    data_call[i,14] = data_call[i,11]
  }
  
  if(sum(names_call %in% "lumpy")==0 & sum(names_call %in% "pindel")==0 & sum(names_call %in% "whamg")==0 &
     sum(names_call %in% "delly")==0 & sum(names_call %in% "manta")==1){
    
    data_call[i,14] = data_call[i,9]
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


data_call_het = data_call[data_call$GT_golden=="0/1",]
data_call_hom = data_call[data_call$GT_golden=="1/1",]

model_g_error_het = 100 -(sum(as.character(data_call_het$GT)==as.character(data_call_het$GT_golden))/
                        nrow(data_call_het))*100
model_g_error_het


model_g_error_hom = 100 -(sum(as.character(data_call_hom$GT)==as.character(data_call_hom$GT_golden))/
                            nrow(data_call_hom))*100
model_g_error_hom


### choose best model (reprocicity and number of callers are discarded because produce large standard deviations of the betas) ####

final_data$PASS = ifelse(final_data$PASS=="YES",0,1)

full.model =  glm(PASS ~ factor(GT_manta) + factor(GT_whamg) + factor(GT_delly) + 
                    factor(GT_lumpy) + factor(GT_pindel) + factor(GT_svaba) + 
                    factor(length_stretch),
                  data=final_data,family = "binomial")

base.model = glm(PASS ~ factor(GT_manta) + factor(GT_whamg) + factor(GT_delly) + 
                   factor(GT_lumpy) + factor(GT_pindel) + factor(GT_svaba),
                 data=final_data,family = "binomial")

step(full.model,
     scope = list(lower = formula(base.model),
                  upper = formula(full.model)),
     direction = "backward")


best.model = glm(formula = PASS ~ factor(GT_manta) + factor(GT_whamg) + factor(GT_delly) + 
                   factor(GT_lumpy) + factor(GT_pindel) + factor(GT_svaba) + 
                   factor(length_stretch), family = "binomial", 
                 data = final_data)

final_data$PASS = ifelse(final_data$PASS==0,"YES","NO")

summary(best.model)

my_pred = predict(best.model,final_data,type="response")

my_pred = ifelse(my_pred>0.5,1,0)

my_pred = ifelse(my_pred==1,"NO PASS","PASS")

table(final_data$PASS)
table(my_pred,final_data$PASS)
table(my_pred)

N = nrow(golden)

model_tp = as.numeric(table(my_pred,final_data$PASS)[2,2])
model_fp = as.numeric(table(my_pred)[2]-model_tp)  
#sensitivity
model_sens = round(model_tp/N*100,3)
#specificity
model_spec = round(model_tp/(model_fp+model_tp)*100,3)

model_tp;model_fp;model_sens;model_spec

f_score = round(((2*model_sens*model_spec)/(model_sens + model_spec))/100,3)

final_data$prediction = my_pred

final_data = as.data.table(final_data)

### name callers detected by prediction #####

final_data$name_callers_detected = name_callers_detected(as.data.table(final_data), 
                                              callers =  c("pindel","manta","whamg","lumpy","svaba","delly"))

table(final_data$name_callers_detected,final_data$prediction) 


# detected by at least 1 caller

det_call = final_data 

callers1 = logical_sens_spec(det_call)

# detected by at least 2 callers and 2 strategies

det_call = final_data %>% filter(callers_detected_ok %in% c(2:6) &  strategy_ok %in% c(2:4) & reciprocity >0.80)

callers2 = logical_sens_spec(det_call)

# detected by at least 3 callers and 2 strategies

det_call = final_data %>% filter(callers_detected_ok %in% c(3:6) & strategy_ok %in% c(2:4) & reciprocity >0.80)

callers3 = logical_sens_spec(det_call)

# detected by at least 4 callers and 2 strategies

det_call = final_data %>% filter(callers_detected_ok %in% c(4:6) & strategy_ok %in% c(2:4) & reciprocity >0.80)

callers4 = logical_sens_spec(det_call)

metrics_callers = rbind(as.data.frame(metrics_callers),
                        c("","","LR",model_tp,model_fp,model_sens,model_spec,f_score,
                          round(model_g_error,3)))

metrics_callers = rbind(metrics_callers,
                        c("","",">=1 caller",callers1,""))

metrics_callers = rbind(metrics_callers,
                        c("","",">=2 caller **",callers2,""))
metrics_callers = rbind(metrics_callers,
                        c("","",">=3 caller **",callers3,""))
metrics_callers = rbind(metrics_callers,
                        c("","",">=4 caller **",callers4,""))


### CV ######

library(caret)
library(e1071)

ctrl <- trainControl(method = "repeatedcv", number = 10, 
                     savePredictions = TRUE)


## 90-10 (small sample size)

0.7*nrow(final_data)

set.seed(2589)

training = sample(1:nrow(final_data),2818)

table(final_data$PASS,useNA = "always")

final_data$PASS = as.factor(final_data$PASS)

final_data_70 = final_data[training,]

table(final_data_70$PASS,final_data_70$length_stretch)
table(final_data_70$GT_whamg,final_data_70$length_stretch)
table(final_data_70$GT_pindel,final_data_70$length_stretch)
table(final_data_70$GT_delly,final_data_70$length_stretch)
table(final_data_70$GT_lumpy,final_data_70$length_stretch)
table(final_data_70$GT_manta,final_data_70$length_stretch)
table(final_data_70$GT_svaba,final_data_70$length_stretch)



table(final_data_70$PASS,useNA = "always")

mod_fit <- train(factor(PASS) ~ factor(GT_manta) + factor(GT_whamg) + factor(GT_delly) + 
                   factor(GT_lumpy) + factor(GT_pindel) + factor(GT_svaba)+factor(length_stretch), 
                 data=final_data_70, method="glm", family="binomial",
                 trControl = ctrl)

summary(mod_fit)

my_pred = predict(mod_fit, newdata=final_data_70)

my_pred = ifelse(my_pred=="YES",0,1)

my_pred = ifelse(my_pred==0,"PASS","NO PASS")

table(final_data_70$PASS)
table(my_pred,final_data_70$PASS)
table(my_pred)

N = 283

tp_cv = as.numeric(table(my_pred,final_data_70$PASS)[2,2])
fp_cv = as.numeric(table(my_pred)[2])-tp_cv

sens_call_cv = round(tp_cv/N*100,3)
spec_call_cv = round(tp_cv/(fp_cv+tp_cv)*100,3)

f_score_cv = round((2*sens_call_cv*spec_call_cv)/(sens_call_cv+spec_call_cv)/100,3)


# detected by at least 1 caller

det_call = final_data_70 %>% filter(callers_detected_ok %in% 1:9)

callers1 = logical_sens_spec(det_call)

# detected by at least 2 callers and 2 strategies

det_call = final_data_70 %>% filter(callers_detected_ok %in% 2:6 & strategy_ok %in% 2:4)

det_call = det_call[which(det_call$reciprocity>0.80),]

callers2 = logical_sens_spec(det_call)

# detected by at least 3 callers and 2 strategies

det_call = final_data_70 %>% filter(callers_detected_ok %in% 3:9 & strategy_ok %in% 2:6)

det_call = det_call[which(det_call$reciprocity>0.80),]

callers3 = logical_sens_spec(det_call)

# detected by at least 4 callers and 2 strategies

det_call = final_data_70 %>% filter(callers_detected_ok %in% 4:9 & strategy_ok %in% 2:6)

det_call = det_call[which(det_call$reciprocity>0.80),]

callers4 = logical_sens_spec(det_call)


metrics_callers = rbind(metrics_callers,
                        c("",N,"LR (100% and CV)",
                          c(tp_cv,fp_cv,sens_call_cv,spec_call_cv,f_score_cv),""))

metrics_callers = rbind(metrics_callers,
                        c("","",">=1 caller",callers1,""))

metrics_callers = rbind(metrics_callers,
                        c("","",">=2 caller **",callers2,""))

metrics_callers = rbind(metrics_callers,
                        c("","",">=3 caller **",callers3,""))

metrics_callers = rbind(metrics_callers,
                        c("","",">=4 caller **",callers4,""))


## test set

final_data_30 = final_data[-training,]

final_data_30$GT_manta[final_data_30$GT_manta=="0/1"] = "0/1-1/1"
final_data_30$GT_manta[final_data_30$GT_manta=="1/1"] = "0/1-1/1"
final_data_30$GT_whamg[final_data_30$GT_whamg=="0/1"] = "0/1-1/1"
final_data_30$GT_whamg[final_data_30$GT_whamg=="1/1"] = "0/1-1/1"
final_data_30$GT_delly[final_data_30$GT_delly=="0/1"] = "0/1-1/1"
final_data_30$GT_delly[final_data_30$GT_delly=="1/1"] = "0/1-1/1"
final_data_30$GT_lumpy[final_data_30$GT_lumpy=="0/1"] = "0/1-1/1"
final_data_30$GT_lumpy[final_data_30$GT_lumpy=="1/1"] = "0/1-1/1"
final_data_30$GT_pindel[final_data_30$GT_pindel=="0/1"] = "0/1-1/1"
final_data_30$GT_pindel[final_data_30$GT_pindel=="1/1"] = "0/1-1/1"
final_data_30$GT_svaba[final_data_30$GT_svaba=="0/1"] = "0/1-1/1"
final_data_30$GT_svaba[final_data_30$GT_svaba=="1/1"] = "0/1-1/1"

final_data_30$GT_manta[is.na(final_data_30$GT_manta)] = "9/9"
final_data_30$GT_whamg[is.na(final_data_30$GT_whamg)] = "9/9"
final_data_30$GT_delly[is.na(final_data_30$GT_delly)] = "9/9"
final_data_30$GT_lumpy[is.na(final_data_30$GT_lumpy)] = "9/9"
final_data_30$GT_pindel[is.na(final_data_30$GT_pindel)] = "9/9"
final_data_30$GT_svaba[is.na(final_data_30$GT_svaba)] = "9/9"

my_pred = predict(mod_fit,final_data_30)

my_pred = ifelse(my_pred=="YES",0,1)

my_pred = ifelse(my_pred==0,"PASS","NO PASS")

table(final_data_30$PASS)
table(my_pred,final_data_30$PASS)
table(my_pred)

N = as.numeric(table(final_data_30$PASS)[2])

tp_cv = as.numeric(table(my_pred,final_data_30$PASS)[2,2])
fp_cv = as.numeric(table(my_pred)[2])-tp_cv

sens_call_cv = round(tp_cv/N*100,3)
spec_call_cv = round(tp_cv/(fp_cv+tp_cv)*100,3)

f_score_cv = round((2*sens_call_cv*spec_call_cv)/(sens_call_cv+spec_call_cv)/100,3)


# detected by at least 1 caller

det_call = final_data_30 %>% filter(callers_detected %in% 1:9)

callers1 = logical_sens_spec(det_call)

# detected by at least 2 callers and 2 strategies

det_call = final_data_30 %>% filter(callers_detected %in% 2:9 & strategy %in% 2:6 & reprocicity>0.8)

callers2 = logical_sens_spec(det_call)

# detected by at least 3 callers and 2 strategies

det_call = final_data_30 %>% filter(callers_detected %in% 3:9 & strategy %in% 2:6 & reprocicity>0.8)

callers3 = logical_sens_spec(det_call)

# detected by at least 4 callers and 2 strategies

det_call = final_data_30 %>% filter(callers_detected %in% 4:9 & strategy %in% 2:6 & reprocicity>0.8)

callers4 = logical_sens_spec(det_call)


metrics_callers = rbind(metrics_callers,
                        c("",N,"LR (test 10%)",
                          c(tp_cv,fp_cv,sens_call_cv,spec_call_cv,f_score_cv),""))

metrics_callers = rbind(metrics_callers,
                        c("","",">=1 caller",callers1,""))

metrics_callers = rbind(metrics_callers,
                        c("","",">=2 caller **",callers2,""))

metrics_callers = rbind(metrics_callers,
                        c("","",">=3 caller **",callers3,""))

metrics_callers = rbind(metrics_callers,
                        c("","",">=4 caller **",callers4,""))

write.csv(metrics_callers,"METRICS/inversions.csv",row.names = F)

metrics_callers = read.csv("METRICS/inversions.csv")

metrics_callers$F1.Score = ((2*metrics_callers$Sensitivity*metrics_callers$Precision)/
                              (metrics_callers$Sensitivity + metrics_callers$Precision))

