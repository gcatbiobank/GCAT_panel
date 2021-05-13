library(data.table)
library(dplyr)

setwd("/home/jordi/Mare_nostrum1/GCAT_project_all_samples/merge_all_calling_GCAT/")
setwd("/home/igalvan/Mare_folder/")

zzz = readRDS("LRM_final_data/Duplications/model_duplications.rds")
summary(zzz)



source("Deleciones/functions_JV_IC.R")

# read golden ######

golden_duplications = fread("Duplications/Golden/insilico_data/Duplis_intra_all")
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

# input Montse script

golden %>% select(chr,start_golden,end_golden,length_golden,diff,window_bp,GT_golden) %>%
  write.table("Duplications/Golden/Genotyping_Montse/golden_duplications_all1.txt",row.names = F,quote = F)



# window size sensitivity and precision #######

win_size = c(10,20,50,100,200,300)

my_callers = expand.grid(callers = c("delly","cnvnator","pindel","whamg","svaba","manta",
                                     "lumpy"),
                         win_size = win_size)

my_callers$sensitivity = 0
my_callers$precision = 0
my_callers$f.score = 0
my_callers$g.error = 0


for(i in 1:6){
  
  # read vcfs files 
  delly = fread("Duplications/Golden/insilico_data/Duplicaciones_Delly")
  colnames(delly) = c("chr","start_delly","end","length_delly","GT_delly")
  delly$lower_delly = delly$start_delly-win_size[i]
  delly$upper_delly = delly$start_delly+win_size[i]
  delly$chr_pos_delly = do.call(paste0,list(delly$chr,"_",delly$start_delly))
  summary(delly$length_delly)
  
  pindel = fread("Duplications/Golden/insilico_data/Duplicaciones_Pindel",fill = T)
  colnames(pindel) = c("chr","start_pindel","end","length_pindel","GT_pindel")
  pindel$lower_pindel = pindel$start_pindel-win_size[i]
  pindel$upper_pindel = pindel$start_pindel+win_size[i]
  pindel$chr_pos_pindel = do.call(paste0,list(pindel$chr,"_",pindel$start_pindel))
  summary(pindel$length_pindel)
  
  whamg = fread("Duplications/Golden/insilico_data/Duplicaciones_whamg")
  colnames(whamg) = c("chr","start_whamg","end","length_whamg","GT_whamg")
  whamg$lower_whamg = whamg$start_whamg-win_size[i]
  whamg$upper_whamg = whamg$start_whamg+win_size[i]
  whamg$chr_pos_whamg = do.call(paste0,list(whamg$chr,"_",whamg$start_whamg))
  summary(whamg$length_whamg)
  
  svaba = fread("Inversions/Golden/insilico_data/SVABA_total_ivan")
  svaba$length = svaba$V2-svaba$V3
  
  # change start and end
  
  svaba$length_aux = svaba$V3-svaba$V2
  
  for(j in 1:nrow(svaba)){
    
    if(svaba$length_aux[j]<0){
      
      aux = svaba$V2[j]
      aux2 = svaba$V3[j]
      
      svaba$V2[j] = aux2
      svaba$V3[j] = aux
      
    }
    
  }
  
  svaba$length_aux = svaba$V3-svaba$V2
  length(which(golden$length<=0))
  svaba = svaba[,c(1:4,6)]
  colnames(svaba) = c("chr","start_svaba","end_svaba","GT_svaba","length_svaba")
  svaba$lower_svaba = svaba$start_svaba-win_size[i]
  svaba$upper_svaba = svaba$start_svaba+win_size[i]
  svaba = unique(svaba)
  svaba$chr_pos_svaba = do.call(paste0,list(svaba$chr,"_",svaba$start_svaba))
  
  
  manta = fread("Duplications/Golden/insilico_data/Duplicaciones_manta")
  colnames(manta) = c("chr","start_manta","end","length_manta","GT_manta")
  manta$lower_manta = manta$start_manta-win_size[i]
  manta$upper_manta = manta$start_manta+win_size[i]
  manta$chr_pos_manta = do.call(paste0,list(manta$chr,"_",manta$start_manta))
  
  lumpy = fread("Duplications/Golden/insilico_data/Duplicaciones_Lumpy")
  colnames(lumpy) = c("chr","start_lumpy","end","length_lumpy","GT_lumpy")
  lumpy$lower_lumpy = lumpy$start_lumpy-win_size[i]
  lumpy$upper_lumpy = lumpy$start_lumpy+win_size[i]
  lumpy$chr_pos_lumpy = do.call(paste0,list(lumpy$chr,"_",lumpy$start_lumpy))
  
  cnvnator = fread("Duplications/Golden/insilico_data/Duplicaciones_CNVator")
  colnames(cnvnator) = c("chr","start_cnvnator","end","length_cnvnator","GT_cnvnator")
  cnvnator$lower_cnvnator = cnvnator$start_cnvnator-win_size[i]
  cnvnator$upper_cnvnator = cnvnator$start_cnvnator+win_size[i]
  cnvnator$chr_pos_cnvnator = do.call(paste0,list(cnvnator$chr,"_",cnvnator$start_cnvnator))
  
  
  # cnvnator
  
  sens_cnvnator = sensitivity_precision_geno_error_deletions_duplications(cnvnator,golden,"duplication","cnvnator")
  
  my_callers[(my_callers$callers=="cnvnator" & my_callers$win_size==win_size[i]),]$sensitivity = sens_cnvnator$sensitivity
  my_callers[(my_callers$callers=="cnvnator" & my_callers$win_size==win_size[i]),]$precision = sens_cnvnator$precision
  my_callers[(my_callers$callers=="cnvnator" & my_callers$win_size==win_size[i]),]$f.score = sens_cnvnator$f1_score
  my_callers[(my_callers$callers=="cnvnator" & my_callers$win_size==win_size[i]),]$g.error = sens_cnvnator$g_error
  
  
  # pindel
  
  sens_pindel = sensitivity_precision_geno_error_deletions_duplications(pindel,golden,"duplication","pindel")
  
  my_callers[(my_callers$callers=="pindel" & my_callers$win_size==win_size[i]),]$sensitivity = sens_pindel$sensitivity
  my_callers[(my_callers$callers=="pindel" & my_callers$win_size==win_size[i]),]$precision = sens_pindel$precision
  my_callers[(my_callers$callers=="pindel" & my_callers$win_size==win_size[i]),]$f.score = sens_pindel$f1_score
  my_callers[(my_callers$callers=="pindel" & my_callers$win_size==win_size[i]),]$g.error = sens_pindel$g_error
  
  # manta
  
  sens_manta = sensitivity_precision_geno_error_deletions_duplications(manta,golden,"duplication","manta")
  
  my_callers[(my_callers$callers=="manta" & my_callers$win_size==win_size[i]),]$sensitivity = sens_manta$sensitivity
  my_callers[(my_callers$callers=="manta" & my_callers$win_size==win_size[i]),]$precision = sens_manta$precision
  my_callers[(my_callers$callers=="manta" & my_callers$win_size==win_size[i]),]$f.score = sens_manta$f1_score
  my_callers[(my_callers$callers=="manta" & my_callers$win_size==win_size[i]),]$g.error = sens_manta$g_error
  
  # svaba
  
  sens_svaba = sensitivity_precision_geno_error_deletions_duplications(svaba,golden,"duplication","svaba")
  
  my_callers[(my_callers$callers=="svaba" & my_callers$win_size==win_size[i]),]$sensitivity = sens_svaba$sensitivity
  my_callers[(my_callers$callers=="svaba" & my_callers$win_size==win_size[i]),]$precision = sens_svaba$precision
  my_callers[(my_callers$callers=="svaba" & my_callers$win_size==win_size[i]),]$f.score = sens_svaba$f1_score
  my_callers[(my_callers$callers=="svaba" & my_callers$win_size==win_size[i]),]$g.error = sens_svaba$g_error
  
  # whamg
  
  sens_whamg = sensitivity_precision_geno_error_deletions_duplications(whamg,golden,"duplication","whamg")
  
  my_callers[(my_callers$callers=="whamg" & my_callers$win_size==win_size[i]),]$sensitivity = sens_whamg$sensitivity
  my_callers[(my_callers$callers=="whamg" & my_callers$win_size==win_size[i]),]$precision = sens_whamg$precision
  my_callers[(my_callers$callers=="whamg" & my_callers$win_size==win_size[i]),]$f.score = sens_whamg$f1_score
  my_callers[(my_callers$callers=="whamg" & my_callers$win_size==win_size[i]),]$g.error = sens_whamg$g_error
  
  # delly
  
  sens_delly = sensitivity_precision_geno_error_deletions_duplications(delly,golden,"duplication","delly")
  
  my_callers[(my_callers$callers=="delly" & my_callers$win_size==win_size[i]),]$sensitivity = sens_delly$sensitivity
  my_callers[(my_callers$callers=="delly" & my_callers$win_size==win_size[i]),]$precision = sens_delly$precision
  my_callers[(my_callers$callers=="delly" & my_callers$win_size==win_size[i]),]$f.score = sens_delly$f1_score
  my_callers[(my_callers$callers=="delly" & my_callers$win_size==win_size[i]),]$g.error = sens_delly$g_error
  
  # lumpy
  
  sens_lumpy = sensitivity_precision_geno_error_deletions_duplications(lumpy,golden,"duplication","lumpy")
  
  my_callers[(my_callers$callers=="lumpy" & my_callers$win_size==win_size[i]),]$sensitivity = sens_lumpy$sensitivity
  my_callers[(my_callers$callers=="lumpy" & my_callers$win_size==win_size[i]),]$precision = sens_lumpy$precision
  my_callers[(my_callers$callers=="lumpy" & my_callers$win_size==win_size[i]),]$f.score = sens_lumpy$f1_score
  my_callers[(my_callers$callers=="lumpy" & my_callers$win_size==win_size[i]),]$g.error = sens_lumpy$g_error
  
}

dir.create("Duplications/Golden/outputs/")
write.csv(my_callers,"Duplications/Golden/outputs/window_size_callers.csv",row.names = F)

## plot results ####

my_callers = read.csv("Duplications/Golden/outputs/window_size_callers.csv")

my_callers$f.score = my_callers$f.score*100 

my_callers$sensitivity = round(my_callers$sensitivity,2)
my_callers$precision = round(my_callers$precision,2)
my_callers$f.score = round(my_callers$f.score,2)
my_callers$g.error = round(my_callers$g.error,2)

my_callers %>% arrange(callers,win_size) %>%
  write.csv("Duplications/Golden/outputs/window_size_duplications.csv",row.names=F)


delly = my_callers %>% filter(callers=="delly")
lumpy = my_callers %>% filter(callers=="lumpy") # 50
cnvnator = my_callers %>% filter(callers=="cnvnator") # 100
whamg = my_callers %>% filter(callers=="whamg")
svaba = my_callers %>% filter(callers=="svaba")
manta = my_callers %>% filter(callers=="manta") # 20
pindel = my_callers %>% filter(callers=="pindel")

# read vcfs files ########

delly = fread("Duplications/Golden/insilico_data/Duplicaciones_Delly")
colnames(delly) = c("chr","start_delly","end","length_delly","GT_delly")
delly$lower_delly = delly$start_delly-10
delly$upper_delly = delly$start_delly+10
delly$chr_pos_delly = do.call(paste0,list(delly$chr,"_",delly$start_delly))
summary(delly$length_delly)

pindel = fread("Duplications/Golden/insilico_data/Duplicaciones_Pindel",fill = T)
colnames(pindel) = c("chr","start_pindel","end","length_pindel","GT_pindel")
pindel$lower_pindel = pindel$start_pindel-10
pindel$upper_pindel = pindel$start_pindel+10
pindel$chr_pos_pindel = do.call(paste0,list(pindel$chr,"_",pindel$start_pindel))
summary(pindel$length_pindel)

whamg = fread("Duplications/Golden/insilico_data/Duplicaciones_whamg")
colnames(whamg) = c("chr","start_whamg","end","length_whamg","GT_whamg")
whamg$lower_whamg = whamg$start_whamg-10
whamg$upper_whamg = whamg$start_whamg+10
whamg$chr_pos_whamg = do.call(paste0,list(whamg$chr,"_",whamg$start_whamg))
summary(whamg$length_whamg)

#whamg = fread("Duplications/Golden/insilico_data/Wham_bayestyper_dup")
#colnames(whamg) = c("chr","start_whamg","end","GT_whamg")
#whamg$lower_whamg = whamg$start_whamg-10
#whamg$upper_whamg = whamg$start_whamg+10
#whamg$chr_pos_whamg = do.call(paste0,list(whamg$chr,"_",whamg$start_whamg))
#whamg$length_whamg = 0
#summary(whamg$length_whamg)

svaba = fread("Inversions/Golden/insilico_data/SVABA_total_ivan")
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

manta = fread("Duplications/Golden/insilico_data/Duplicaciones_manta")
colnames(manta) = c("chr","start_manta","end","length_manta","GT_manta")
manta$lower_manta = manta$start_manta-20
manta$upper_manta = manta$start_manta+20
manta$chr_pos_manta = do.call(paste0,list(manta$chr,"_",manta$start_manta))

lumpy = fread("Duplications/Golden/insilico_data/Duplicaciones_Lumpy")
colnames(lumpy) = c("chr","start_lumpy","end","length_lumpy","GT_lumpy")
lumpy$lower_lumpy = lumpy$start_lumpy-50
lumpy$upper_lumpy = lumpy$start_lumpy+50
lumpy$chr_pos_lumpy = do.call(paste0,list(lumpy$chr,"_",lumpy$start_lumpy))

cnvnator = fread("Duplications/Golden/insilico_data/Duplicaciones_CNVator")
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

write.csv(metrics_callers,"LRM_final_data/Duplications/table_metrics_duplications_new_model.csv",row.names = F)

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

saveRDS(my_data2,"Duplications/Golden/outputs/merge_callers_JV_new.rds")

my_data2 = readRDS("Duplications/Golden/outputs/merge_callers_JV_new.rds") #new model

#my_data2 = readRDS("Duplications/Golden/outputs/merge_callers.rds") # old model
dim(my_data2)
head(my_data2)
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

strategies = fread("Deleciones/strategies.csv")

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
table(my_data2$callers_detected,my_data2$PASS) ## callers detected, strategy are not informative


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


#saveRDS(my_data2,"Duplications/Golden/outputs/merge_data_duplications.rds")

### Genotyping Montse+Jordi #####
all_model_ok = my_data2 %>% filter(prediction == "PASS")
all_model_ok = all_model_ok %>% filter(PASS == "YES")
dim(all_model_ok)

geno_montse = fread("Duplications/Golden/Genotyping_Montse/insilico3_log_regression_gt_dup_mon_all.vcf")


geno_montse$chr_pos_golden = do.call(paste0,list(geno_montse$V1,"_",geno_montse$V2))


geno_montse_final = geno_montse[geno_montse$chr_pos_golden %in% all_model_ok$chr_pos_golden,]

geno_montse_final= distinct(geno_montse_final, chr_pos_golden,.keep_all = TRUE)

head(geno_montse_final)

homo_golden = nrow(geno_montse_final[geno_montse_final$V7=="1/1",])
het_golden = nrow(geno_montse_final[geno_montse_final$V7=="0/1",])

g.error = round(100-sum(geno_montse_final$V7==geno_montse_final$V9)/nrow(geno_montse_final)*100,2)


homo_golden_het_caller = nrow(geno_montse_final[geno_montse_final$V7=="1/1" & V9=="0/1",])
homo_golden_hom_caller = nrow(geno_montse_final[geno_montse_final$V7=="1/1" & V9=="1/1",])
het_golden_het_caller = nrow(geno_montse_final[geno_montse_final$V7=="0/1" & V9=="0/1",])
het_golden_hom_caller = nrow(geno_montse_final[geno_montse_final$V7=="0/1" & V9=="1/1",])


g_error_het = round(100-nrow(geno_montse_final[geno_montse_final$V7=="0/1" & geno_montse_final$V9=="0/1",])/het_golden*100,digits=2)
g_error_hom = round(100-nrow(geno_montse_final[geno_montse_final$V7=="1/1" & geno_montse_final$V9=="1/1",])/homo_golden*100,digits=2)

line = as.data.frame(t(c("","","LR",model_tp,model_fp,model_sens,IC_sens,model_spec,IC_spec,f_score,g.error,
                         het_golden,g_error_het,homo_golden,g_error_hom)))
colnames(line) <- c("SV_type","N","Caller","TP","FP","Recall","CI Recall","Precision","CI Precision","F1-Score","Genotype error",
                    "N 0/1","Genotype error 0/1","N 1/1","Genotype error 1/1")

metrics_callers = rbind(metrics_callers,line)



###########################


my_data2 = as.data.table(my_data2)

my_data2$start = consensus_start_duplications(my_data2,
                             callers = c("lumpy","manta","pindel","delly",
                                         "whamg","svaba","cnvnator"))

my_data2$end = my_data2$start + my_data2$length

model_data = my_data2 %>% filter(prediction=="PASS")

model_data = model_data[model_data$start_golden %in% golden$start_golden,] 

model_data = model_data %>% arrange(chr,start)

model_data$diff = model_data$start

model_data$diff = ifelse(model_data$diff<0,0,model_data$diff)

model_data$window_bp = 50

model_data %>% select(chr,start,end,length,diff,window_bp,GT_golden) %>%
  write.table("Duplications/Golden/Genotyping_Montse/logistic_model_duplications_insilico_all_JV.txt",row.names = F,quote = F)



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

# que lo detecte al menos 1 caller

det_call =  my_data3 %>% filter(callers_detected %in% 1:9 & reciprocity >=80)

callers1 = logical_sens_spec(det_call)

line_call1=as.data.frame(t(c("","",">=1 caller",callers1,NA,NA,NA,NA,NA)))
colnames(line_call1) <- c("SV_type","N","Caller","TP","FP","Recall","CI Recall","Precision","CI Precision","F1-Score","Genotype error",
                          "N 0/1","Genotype error 0/1","N 1/1","Genotype error 1/1")

metrics_callers = rbind(metrics_callers,line_call1)


# que lo detecte al menos 2 caller y 2 estrategias a la vez

det_call = my_data3 %>% filter(callers_detected %in% 2:9 & strategy %in% 2:6 & reciprocity >=80)

callers2 = logical_sens_spec(det_call)

line_call2=as.data.frame(t(c("","",">=2 caller **",callers2,NA,NA,NA,NA,NA)))
colnames(line_call2) <- c("SV_type","N","Caller","TP","FP","Recall","CI Recall","Precision","CI Precision","F1-Score","Genotype error",
                          "N 0/1","Genotype error 0/1","N 1/1","Genotype error 1/1")

metrics_callers = rbind(metrics_callers,line_call2)


# que lo detecte al menos 3 caller

det_call = my_data3 %>% filter(callers_detected %in% 3:5 & strategy %in% 2:6 & reciprocity >80)

callers3 = logical_sens_spec(det_call)

line_call3=as.data.frame(t(c("","",">=3 caller **",callers3,NA,NA,NA,NA,NA)))
colnames(line_call3) <- c("SV_type","N","Caller","TP","FP","Recall","CI Recall","Precision","CI Precision","F1-Score","Genotype error",
                          "N 0/1","Genotype error 0/1","N 1/1","Genotype error 1/1")

metrics_callers = rbind(metrics_callers,line_call3)



# que lo detecte al menos 4 caller

det_call = my_data3 %>% filter(callers_detected %in% 4:5 & strategy %in% 2:6 & reciprocity >80)

callers4 = logical_sens_spec(det_call)

line_call4=as.data.frame(t(c("","",">=4 caller **",callers4,NA,NA,NA,NA,NA)))
colnames(line_call4) <- c("SV_type","N","Caller","TP","FP","Recall","CI Recall","Precision","CI Precision","F1-Score","Genotype error",
                          "N 0/1","Genotype error 0/1","N 1/1","Genotype error 1/1")

metrics_callers = rbind(metrics_callers,line_call4)


write.csv(metrics_callers,"LRM_final_data/Duplications/table_metrics_duplications_new_model.csv",row.names = F)


### CV ######

library(caret)
library(e1071)

ctrl <- trainControl(method = "repeatedcv", number = 10, 
                     savePredictions = TRUE)

## 70-30

z = my_data2[!is.na(my_data2$length_stretch),] #new model

#z = my_data2 # old model

0.7*nrow(z)
dim(z)

set.seed(2589)

training = sample(1:nrow(z),2392) #new model
#training = sample(1:nrow(z),2449) # old model


table(my_data2_70$length_stretch,useNA = "always")

z$PASS = as.factor(z$PASS)

my_data2_70 = z[training,]

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




#saveRDS(mod_fit,"Duplications/Golden/model_duplications.rds")

my_data3_70 = my_data2_70 %>% filter(name_callers_detected != "svaba") #Svaba is eliminated because is not able to classify by type the SV
my_data3_70$callers_detected <- gsub("4-5-6-7", 4,my_data3_70$callers_detected )
my_data3_70$strategy <- gsub("3-4-5", 3,my_data3_70$strategy )


# que lo detecte al menos 1 caller

det_call = my_data3_70 %>% filter(callers_detected %in% 1:9)

callers1 = logical_sens_spec(det_call)

line_call1=as.data.frame(t(c("","",">=1 caller (70%)",callers1,NA,NA,NA,NA,NA)))
colnames(line_call1) <- c("SV_type","N","Caller","TP","FP","Recall","CI Recall","Precision","CI Precision","F1-Score","Genotype error",
                          "N 0/1","Genotype error 0/1","N 1/1","Genotype error 1/1")

metrics_callers = rbind(metrics_callers,line_call1)

# que lo detecte al menos 2 caller y 2 estrategias a la vez

det_call = my_data3_70 %>% filter(callers_detected %in% 2:9 & strategy %in% 2:6 & reciprocity >=80)

callers2 = logical_sens_spec(det_call)

line_call2=as.data.frame(t(c("","",">=2 caller (70%) **",callers2,NA,NA,NA,NA,NA)))
colnames(line_call2) <- c("SV_type","N","Caller","TP","FP","Recall","CI Recall","Precision","CI Precision","F1-Score","Genotype error",
                          "N 0/1","Genotype error 0/1","N 1/1","Genotype error 1/1")

metrics_callers = rbind(metrics_callers,line_call2)


# que lo detecte al menos 3 caller

det_call = my_data3_70 %>% filter(callers_detected %in% 3:5 & strategy %in% 2:6 & reciprocity >80)

callers3 = logical_sens_spec(det_call)

line_call3=as.data.frame(t(c("","",">=3 caller (70%) **",callers3,NA,NA,NA,NA,NA)))
colnames(line_call3) <- c("SV_type","N","Caller","TP","FP","Recall","CI Recall","Precision","CI Precision","F1-Score","Genotype error",
                          "N 0/1","Genotype error 0/1","N 1/1","Genotype error 1/1")

metrics_callers = rbind(metrics_callers,line_call3)


# que lo detecte al menos 4 caller

det_call = my_data3_70 %>% filter(callers_detected %in% 4:5 & strategy %in% 2:6 & reciprocity >80)

callers4 = logical_sens_spec(det_call)

line_call4=as.data.frame(t(c("","",">=4 caller (70%) **",callers4,NA,NA,NA,NA,NA)))
colnames(line_call4) <- c("SV_type","N","Caller","TP","FP","Recall","CI Recall","Precision","CI Precision","F1-Score","Genotype error",
                          "N 0/1","Genotype error 0/1","N 1/1","Genotype error 1/1")

metrics_callers = rbind(metrics_callers,line_call4)

## Validation 30%

my_data2_30 = z[-training,]

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


# que lo detecte al menos 1 caller

det_call = my_data3_30 %>% filter(callers_detected %in% 1:9)

callers1 = logical_sens_spec(det_call)

line_call1=as.data.frame(t(c("","",">=1 caller (30%)",callers1,NA,NA,NA,NA,NA)))
colnames(line_call1) <- c("SV_type","N","Caller","TP","FP","Recall","CI Recall","Precision","CI Precision","F1-Score","Genotype error",
                          "N 0/1","Genotype error 0/1","N 1/1","Genotype error 1/1")

metrics_callers = rbind(metrics_callers,line_call1)

# que lo detecte al menos 2 caller y 2 estrategias a la vez

det_call = my_data3_30 %>% filter(callers_detected %in% 2:9 & strategy %in% 2:6 & reciprocity >=80)

callers2 = logical_sens_spec(det_call)

line_call2=as.data.frame(t(c("","",">=2 caller (30%) **",callers2,NA,NA,NA,NA,NA)))
colnames(line_call2) <- c("SV_type","N","Caller","TP","FP","Recall","CI Recall","Precision","CI Precision","F1-Score","Genotype error",
                          "N 0/1","Genotype error 0/1","N 1/1","Genotype error 1/1")

metrics_callers = rbind(metrics_callers,line_call2)


# que lo detecte al menos 3 caller

det_call = my_data3_30 %>% filter(callers_detected %in% 3:5 & strategy %in% 2:6 & reciprocity >80)

callers3 = logical_sens_spec(det_call)

line_call3=as.data.frame(t(c("","",">=3 caller (30%) **",callers3,NA,NA,NA,NA,NA)))
colnames(line_call3) <- c("SV_type","N","Caller","TP","FP","Recall","CI Recall","Precision","CI Precision","F1-Score","Genotype error",
                          "N 0/1","Genotype error 0/1","N 1/1","Genotype error 1/1")

metrics_callers = rbind(metrics_callers,line_call3)


# que lo detecte al menos 4 caller

det_call = my_data3_30 %>% filter(callers_detected %in% 4:5 & strategy %in% 2:6 & reciprocity >80)

callers4 = logical_sens_spec(det_call)

line_call4=as.data.frame(t(c("","",">=4 caller (30%) **",callers4,NA,NA,NA,NA,NA)))
colnames(line_call4) <- c("SV_type","N","Caller","TP","FP","Recall","CI Recall","Precision","CI Precision","F1-Score","Genotype error",
                          "N 0/1","Genotype error 0/1","N 1/1","Genotype error 1/1")

metrics_callers = rbind(metrics_callers,line_call4)

write.csv(metrics_callers,"LRM_final_data/Duplications/table_metrics_duplications_new_model.csv",row.names = F)

###########################################################################################################################################


head(metrics_callers)
library(xtable)
print(xtable(metrics_callers[,-1],digits = 2),include.rownames = F)


#### Plot graphics ##########

library(data.table)
library(ggplot2)
library(plyr)

aux = cut(golden$length_golden,
          breaks = c(30,75,150,300,500,1000,2000,3000,Inf),right = F)
sum(table(aux))
golden$length_stretch = aux

golden_length = golden %>% group_by(length_stretch) %>% 
  dplyr::summarise(golden = n()) %>% arrange(length_stretch) %>% as.data.frame()


### model #####

my_data2 = readRDS("Duplications/outputs/merge_callers_predicted.rds")

aux = cut(my_data2$length_golden,
          breaks = c(30,75,150,300,500,1000,2000,3000,Inf),right = F)

aux2 = cut(my_data2$length,
           breaks = c(30,75,150,300,500,1000,2000,3000,Inf),right = F)

aux[is.na(aux)] = aux2[is.na(aux)]

table(aux)
my_data2$length_stretch = aux

sensitivity_model = my_data2 %>% filter(!is.na(GT_golden) & prediction=="PASS") %>% group_by(length_stretch) %>% 
  dplyr::summarise(sensitivity = n()) %>% arrange(length_stretch) %>% as.data.frame()

sensitivity_model$sensitivity = sensitivity_model$sensitivity/golden_length$golden*100 

precision_model = my_data2 %>% filter(!is.na(GT_golden) & prediction=="PASS") %>% group_by(length_stretch) %>% 
  dplyr::summarise(truepos = n()) %>% arrange(length_stretch) %>% as.data.frame()

true_falsepositive_model = my_data2 %>% filter(prediction=="PASS") %>% group_by(length_stretch) %>% 
  dplyr::summarise(true_falsepos = n()) %>% arrange(length_stretch) %>% as.data.frame()

precision_model$precision = precision_model$truepos/(true_falsepositive_model$true_falsepos[-9])*100 

sens_precision = left_join(sensitivity_model,precision_model %>% dplyr::select(-truepos))

sens_precision$caller = "Logistic regression model"

sens_precision_final = sens_precision

sens_precision_final


### callers #####

sensitivity_precision_caller <- function(my_data,golden,caller_name){
  
  ### my_data es el dataset con el merge de todos los callers
  ### golden es el numero de variantes para cada tamaño
  ### el nombre del caller
  
  # my_data = my_data2
  # golden = golden_length
  # caller_name = "gatk"
  
  template_sens = golden
  colnames(template_sens)[2] = "sensitivity"
  template_sens$sensitivity = 0
  
  template_prec = golden
  colnames(template_prec)[2] = "precision"
  template_prec$precision = 0
  
  
  aux = my_data[as.character(t(my_data[,paste0("GT_",caller_name),with=F]))!="9/9",]
  
  sensitivity = aux %>% filter(!is.na(GT_golden)) %>% group_by(length_stretch) %>% 
    dplyr::summarise(sensitivity = n()) %>% arrange(length_stretch) %>% as.data.frame()
  
  template_sens = left_join(template_sens %>% dplyr::select(length_stretch),sensitivity)
  
  if(sum(is.na(template_sens$sensitivity))>0){
    template_sens[is.na(template_sens$sensitivity),]$sensitivity=0}
  
  template_sens$sensitivity = template_sens$sensitivity/golden$golden*100 
  
  precision = aux %>% filter(!is.na(GT_golden)) %>% group_by(length_stretch) %>% 
    dplyr::summarise(truepos = n()) %>% arrange(length_stretch) %>% as.data.frame()
  
  template_prec = left_join(template_prec %>% dplyr::select(length_stretch),precision)
  
  if(sum(is.na(template_prec$truepos))>0){
    template_prec[is.na(template_prec$truepos),]$truepos=0}
  
  true_falsepositive = aux %>% group_by(length_stretch) %>% 
    dplyr::summarise(true_falsepos = n()) %>% arrange(length_stretch) %>% as.data.frame()
  
  template_prec = left_join(template_prec %>% dplyr::select(length_stretch,truepos),true_falsepositive)
  
  if(sum(is.na(template_prec$true_falsepos))>0){
    template_prec[is.na(template_prec$true_falsepos),]$true_falsepos=1}
  
  template_prec$precision = template_prec$truepos/(template_prec$true_falsepos)*100 
  
  sens_precision = left_join(template_sens,template_prec %>% dplyr::select(-truepos,-true_falsepos))
  
  sens_precision$caller = caller_name
  
  return(sens_precision)
}


sens_precision = sensitivity_precision_caller(my_data2,golden_length,caller_name ="lumpy")
sens_precision_final = rbind(sens_precision_final,sens_precision)

sens_precision = sensitivity_precision_caller(my_data2,golden_length,caller_name ="cnvnator")
sens_precision_final = rbind(sens_precision_final,sens_precision)

sens_precision = sensitivity_precision_caller(my_data2,golden_length,caller_name ="delly")
sens_precision_final = rbind(sens_precision_final,sens_precision)

sens_precision = sensitivity_precision_caller(my_data2,golden_length,caller_name ="manta")
sens_precision_final = rbind(sens_precision_final,sens_precision)

sens_precision = sensitivity_precision_caller(my_data2,golden_length,caller_name ="svaba")
sens_precision_final = rbind(sens_precision_final,sens_precision)

sens_precision = sensitivity_precision_caller(my_data2,golden_length,caller_name ="whamg")
sens_precision_final = rbind(sens_precision_final,sens_precision)

sens_precision = sensitivity_precision_caller(my_data2,golden_length,caller_name ="pindel")
sens_precision_final = rbind(sens_precision_final,sens_precision)



### plot #####

sens_precision_final$sensitivity = -sens_precision_final$sensitivity

sens_precision_final$caller = as.factor(sens_precision_final$caller)

sens_precision_final$color = factor(sens_precision_final$caller)

sens_precision_final$color = mapvalues(sens_precision_final$color , 
                                       from = levels(sens_precision_final$caller), 
                                       to = c("blue3","firebrick2","black", "gray",
                                              "darkorange2","darkorchid2",
                                              "deeppink2","darkgoldenrod2"))


# [1] "deepvariant" dodgerblue             
# [2] "delly" firebrick2                   
# [3] "gatk"  forestgreen                   
# [4] "Logistic regression model" black
# [5] "lumpy"         gray           
# [6] "manta"          darkorange2          
# [7] "pindel"         darkorchid2          
# [8] "strelka"        cadetblue3          
# [9] "svaba"          deeppink2          
# [10] "whamg"         darkgoldenrod2
# [11] "pamir" brown
# [12] "popins" darkolivegreen2
# [13] "cnvnator" blue3

breaks = c(150,300,500,1000,2000,3000)
breaks_labels = c("150-300","300-500","500-1K","1K-2K","2K-3K","3K->3K")

c(30,75,150,300,500,1000,2000,3000,Inf)

breaks = c(30,75,150,300,500,1000,2000,3000)
breaks_labels = c("31-75","75-150","150-300","300-500","500-1K","1K-2K","2K-3K","3K->3K")

sens_precision_final$breaks = breaks

sens_precision_final[sens_precision_final==0]=NA

model_in = sens_precision_final %>% filter(caller %in% "Logistic regression model")

sens_precision_final[sens_precision_final$caller  %in% "Logistic regression model", ]$sensitivity = NA
sens_precision_final[sens_precision_final$caller  %in% "Logistic regression model", ]$precision = NA

p1 = ggplot(sens_precision_final, aes(x = breaks, y = precision,  color = caller, group = caller)) + 
  geom_point() + geom_line() +
  scale_color_manual(values = levels(factor(sens_precision_final$color))) + theme_classic() +
  ggtitle(c("505 Duplications"))+
  scale_x_log10(breaks = breaks,
                labels = breaks_labels,
                limits = c(30,3000))+ 
  scale_y_continuous(breaks = seq(-100,100,10),labels=c(rev(seq(10,100,10)),seq(0,100,10)),
                     limits = c(-100,100)) + 
  xlab("Size (bp)") + guides(color=guide_legend(title="Variant caller")) +
  geom_point( aes(x = breaks, y = precision,  color = caller, group = caller),model_in,size=2) +
  geom_line( aes(x = breaks, y = precision,  color = caller, group = caller),model_in,size=1.2)

p1 +  geom_line(sens_precision_final,mapping = aes(x = breaks, y = sensitivity,  color = caller)) + 
  geom_point(sens_precision_final,mapping = aes(x = breaks, y = sensitivity,  color = caller)) +
  geom_hline(yintercept = 0,linetype="dashed") + ylab("Sensitivity            (%)            Precision") +
  geom_point( aes(x = breaks, y = sensitivity,  color = caller, group = caller),model_in,size=2) +
  geom_line( aes(x = breaks, y = sensitivity,  color = caller, group = caller),model_in,size=1.2)



pdf("Duplications/outputs/Duplications_.pdf",width = 11,height = 5)
p1 +  geom_line(sens_precision_final,mapping = aes(x = breaks, y = sensitivity,  color = caller)) + 
  geom_point(sens_precision_final,mapping = aes(x = breaks, y = sensitivity,  color = caller)) +
  geom_hline(yintercept = 0,linetype="dashed") + ylab("Sensitivity            (%)            Precision") +
  geom_point( aes(x = breaks, y = sensitivity,  color = caller, group = caller),model_in,size=2) +
  geom_line( aes(x = breaks, y = sensitivity,  color = caller, group = caller),model_in,size=1.2)
dev.off()









