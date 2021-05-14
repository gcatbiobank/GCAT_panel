library(data.table)
library(dplyr)
library(caret)
library(e1071)
              
args <- commandArgs(trailingOnly = TRUE)
i= as.numeric(args[2])
              
source("ext/functions.R")
              
strategies = fread("ext/strategies.csv")
              
mod_fit = readRDS("1_LRM/models/LRM_model_mid_DEL.rds")
              
ids = fread("ext/all_samplesok",header = F)
              
ids = ids$V1

print(ids[i])
              
dir.create(paste0("/2_merge_callers/mid_DEL/outputs/",ids[i]))
              
j= as.numeric(args[1])
 
# read VCF files ########
              
delly = fread(paste0("/2_merge_callers/mid_DEL/data/Delly/",ids[i],"_mid_DEL"))
delly$V4 = as.numeric(abs(delly$V4))
colnames(delly) = c("chr","start_delly","end","length_delly","GT_delly")
delly$lower_delly = delly$start_delly-5
delly$upper_delly = delly$start_delly+5
delly$chr_pos_delly = do.call(paste0,list(delly$chr,"_",delly$start_delly))
delly$chr[delly$chr=="X"] = 23
delly$chr[delly$chr=="Y"] = 24
delly$chr = as.character(delly$chr)
delly = delly %>% filter(length_delly>30 & length_delly<151)
              
                      
lumpy = fread(paste0("/2_merge_callers/mid_DEL/data/Lumpy/",ids[i],"_mid_DEL"))
lumpy$V4 = as.numeric(abs(lumpy$V4))
colnames(lumpy) = c("chr","start_lumpy","end","length_lumpy","GT_lumpy")
lumpy$lower_lumpy = lumpy$start_lumpy-10
lumpy$upper_lumpy = lumpy$start_lumpy+10
lumpy$chr_pos_lumpy = do.call(paste0,list(lumpy$chr,"_",lumpy$start_lumpy))
lumpy$chr[lumpy$chr=="X"] = 23
lumpy$chr[lumpy$chr=="Y"] = 24
lumpy$chr = as.character(lumpy$chr)
lumpy = lumpy %>% filter(length_lumpy>30 & length_lumpy<151)
              
pindel = fread(paste0("/2_merge_callers/mid_DEL/data/Pindel/",ids[i],"_mid_DEL"))
colnames(pindel) = c("chr","start_pindel","end","length_pindel","GT_pindel")
pindel$lower_pindel = pindel$start_pindel-5
pindel$upper_pindel = pindel$start_pindel+5
pindel$chr_pos_pindel = do.call(paste0,list(pindel$chr,"_",pindel$start_pindel))
pindel$chr[pindel$chr=="X"] = 23
pindel$chr[pindel$chr=="Y"] = 24
pindel$chr = as.character(pindel$chr)
pindel = pindel %>% filter(length_pindel>30 & length_pindel<151)
              
whamg = fread(paste0("/2_merge_callers/mid_DEL/data/Whamg/",ids[i],"_mid_DEL"))
whamg$V4 = as.numeric(abs(whamg$V4))
colnames(whamg) = c("chr","start_whamg","end","length_whamg","GT_whamg")
whamg$lower_whamg = whamg$start_whamg-10
whamg$upper_whamg = whamg$start_whamg+10
whamg$chr_pos_whamg = do.call(paste0,list(whamg$chr,"_",whamg$start_whamg))
whamg$chr[whamg$chr=="X"] = 23
whamg$chr[whamg$chr=="Y"] = 24
whamg$chr = as.character(whamg$chr)
whamg = whamg %>% filter(length_whamg>30 & length_whamg<151)
              
svaba = fread(paste0("/2_merge_callers/mid_DEL/data/SVaBA/",ids[i],"_mid_DEL"))
svaba$length_svaba = abs(nchar(svaba$V3)-nchar(svaba$V4))
svaba = svaba %>% filter(length_svaba>30 & length_svaba<151)
colnames(svaba)[c(1,2,5)] = c("chr","start_svaba","GT_svaba")
svaba = svaba[,c(1,2,5,6)]
svaba$lower_svaba = svaba$start_svaba-5
svaba$upper_svaba = svaba$start_svaba+5
svaba$chr_pos_svaba = do.call(paste0,list(svaba$chr,"_",svaba$start_svaba))
svaba$chr[svaba$chr=="X"] = 23
svaba$chr[svaba$chr=="Y"] = 24
svaba$chr = as.character(svaba$chr)
svaba = svaba %>% filter(length_svaba>30 & length_svaba<151)
              
              
manta = fread(paste0("/2_merge_callers/mid_DEL/data/Manta/",ids[i],"_mid_DEL"))
manta$V4 = as.numeric(abs(manta$V4))
colnames(manta) = c("chr","start_manta","end","length_manta","GT_manta")
manta$lower_manta = manta$start_manta-5
manta$upper_manta = manta$start_manta+5
manta$chr_pos_manta = do.call(paste0,list(manta$chr,"_",manta$start_manta))
manta$chr[manta$chr=="X"] = 23
manta$chr[manta$chr=="Y"] = 24
manta$chr = as.character(manta$chr)
              
            
manta2 = fread(paste0("/2_merge_callers/mid_DEL/data/Manta/",ids[i],"_mid_DEL_30_50"))
manta2$length_manta = abs(nchar(manta2$V3)-nchar(manta2$V4))
manta2$GT_manta = "2/2"
manta2 = manta2 %>% filter(length_manta>30 & length_manta<51)
colnames(manta2)[c(1,2)] = c("chr","start_manta")
manta2$end = manta2$length_manta+manta2$start_manta
manta2$lower_manta = manta2$start_manta-5
manta2$upper_manta = manta2$start_manta+5
manta2$chr_pos_manta = do.call(paste0,list(manta2$chr,"_",manta2$start_manta))
manta2$chr[manta2$chr=="X"] = 23
manta2$chr[manta2$chr=="Y"] = 24
manta2$chr =as.character(manta2$chr)
manta2 = manta2[,c("chr","start_manta","end","length_manta","GT_manta","lower_manta","upper_manta","chr_pos_manta")]
manta = rbind(manta,manta2)
manta = manta %>% filter(length_manta>30 & length_manta<151)
              
gatk = fread(paste0("/2_merge_callers/mid_DEL/data/Gatk/",ids[i],"_mid_DEL"))
gatk$length_gatk = abs(nchar(gatk$V3)-nchar(gatk$V4))
colnames(gatk)[c(1,2,5)] = c("chr","start_gatk","GT_gatk")
gatk = gatk[,c(1,2,5,6)]
gatk$end = gatk$length_gatk+gatk$start_gatk
gatk$lower_gatk = gatk$start_gatk-5
gatk$upper_gatk = gatk$start_gatk+5
gatk$chr_pos_gatk = do.call(paste0,list(gatk$chr,"_",gatk$start_gatk))
gatk$chr[gatk$chr=="X"] = 23
gatk$chr[gatk$chr=="Y"] = 24
gatk$chr = as.character(gatk$chr)
gatk = gatk %>% filter(length_gatk>30 & length_gatk<151)
              
deepvariant = fread(paste0("/2_merge_callers/mid_DEL/data/Deepvariant/",ids[i],"_mid_DEL"))
deepvariant$length_deepvariant = abs(nchar(deepvariant$V3)-nchar(deepvariant$V4))
colnames(deepvariant)[c(1,2,5)] = c("chr","start_deepvariant","GT_deepvariant")
deepvariant = deepvariant[,c(1,2,5,6)]
deepvariant$end = deepvariant$length_deepvariant+deepvariant$start_deepvariant
deepvariant$lower_deepvariant = deepvariant$start_deepvariant-5
deepvariant$upper_deepvariant = deepvariant$start_deepvariant+5
deepvariant$chr_pos_deepvariant = do.call(paste0,list(deepvariant$chr,"_",deepvariant$start_deepvariant))
deepvariant$chr[deepvariant$chr=="X"] = 23
deepvariant$chr[deepvariant$chr=="Y"] = 24
deepvariant$chr = as.character(deepvariant$chr)
deepvariant = deepvariant %>% filter(length_deepvariant>30 & length_deepvariant<151)
              
strelka = fread(paste0("/2_merge_callers/mid_DEL/data/Strelka/",ids[i],"_mid_DEL"))
strelka$length_strelka = abs(nchar(strelka$V3)-nchar(strelka$V4))
colnames(strelka)[c(1,2,5)] = c("chr","start_strelka","GT_strelka")
strelka = strelka[,c(1,2,5,6)]
strelka$end = strelka$length_strelka+strelka$start_strelka
strelka$lower_strelka = strelka$start_strelka-5
strelka$upper_strelka = strelka$start_strelka+5
strelka$chr_pos_strelka = do.call(paste0,list(strelka$chr,"_",strelka$start_strelka))
strelka$chr[strelka$chr=="X"] = 23
strelka$chr[strelka$chr=="Y"] = 24
strelka$chr = as.character(strelka$chr)
strelka = strelka %>% filter(length_strelka>30 & length_strelka<151)
              
              
              
call_windows = data.frame(caller = c("delly","lumpy","pindel","whamg","svaba",
                                                   "manta","gatk","deepvariant","strelka"),
                                        window = c(5,10,5,10,5,5,5,5,5))
              
call_windows$caller = as.character(call_windows$caller)
              
                             
# filter by chromosome ########
                
delly_chr = delly %>% filter(chr==j) %>% arrange(start_delly) %>% as.data.table() %>% unique()
lumpy_chr = lumpy %>% filter(chr==j) %>% arrange(start_lumpy) %>% as.data.table() %>% unique()
pindel_chr = pindel %>% filter(chr==j) %>% arrange(start_pindel) %>% as.data.table() %>% unique()
whamg_chr = whamg %>% filter(chr==j) %>% arrange(start_whamg) %>% as.data.table() %>% unique()
svaba_chr = svaba %>% filter(chr==j) %>% arrange(start_svaba) %>% as.data.table() %>% unique()
manta_chr = manta %>% filter(chr==j) %>% arrange(start_manta) %>% as.data.table() %>% unique()
gatk_chr = gatk %>% filter (chr==j) %>% arrange(start_gatk) %>% as.data.table() %>% unique()
deepvariant_chr = deepvariant %>% filter (chr==j) %>% arrange(start_deepvariant) %>% as.data.table() %>% unique()
strelka_chr = strelka %>% filter (chr==j) %>% arrange(start_strelka) %>% as.data.table() %>% unique()
              
               
n_rows = nrow(delly_chr) + nrow(lumpy_chr) + nrow(pindel_chr) + nrow(whamg_chr) + nrow(svaba_chr) + nrow(manta_chr) + nrow(gatk_chr) + nrow(deepvariant_chr) + nrow(strelka_chr)
                
if(n_rows!=0){
                
## prepare BBDD ordered by F-score ####
                
f_score_order = c("delly","manta","svaba","lumpy",
                                  "whamg","gatk","pindel","deepvariant","strelka")
                
delly_chr$ID = paste0("deletion_",1:nrow(delly_chr))
manta_chr$ID = "none"
svaba_chr$ID = "none"
lumpy_chr$ID = "none"
whamg_chr$ID = "none"
gatk_chr$ID = "none"
pindel_chr$ID = "none"
deepvariant_chr$ID = "none"
strelka_chr$ID = "none"
                
                
caller_first = 1
                
while(nrow(get(paste0(f_score_order[caller_first],"_chr")))==0){
                  
                  caller_first = caller_first + 1
  }
                
caller_second = caller_first + 1
                
while(nrow(get(paste0(f_score_order[caller_second],"_chr")))==0){
                  
                  caller_second = caller_second + 1
  }
                
                
first_call = f_score_order[caller_first]
second_call = f_score_order[caller_second]
                
my_data = merge_callers_deletions_duplications(get(paste0(second_call,"_chr")),get(paste0(first_call,"_chr")),callers_ref = first_call,
                                        caller_to_merge = second_call,repro=0.8,svtype = "deletion") 
                
for(caller_merge in (caller_second+1):length(f_score_order)){
                  
  n_row_call = nrow(get(paste0(f_score_order[caller_merge],"_chr")))
                  
    while(n_row_call!=0){
              
                    callers_ref = tstrsplit(colnames(my_data)[grep("GT_",colnames(my_data))],split="GT_")[[2]]
                    
                    my_data = merge_callers_deletions_duplications(get(paste0(f_score_order[caller_merge],"_chr")),my_data,callers_ref = callers_ref,
                                            caller_to_merge = f_score_order[caller_merge],repro=0.8, svtype = "deletion") 
                    
                    n_row_call = 0
                          
                  }
                  
                }
                
                
callers_all = tstrsplit(colnames(my_data)[grep("GT_",colnames(my_data))],split="GT_")[[2]]
                
                
## consensus genotype #####
                
my_data$GT_manta_ok =  my_data$GT_manta
                
# remove 2/2 and replace for NA
                
my_data$GT_manta[my_data$GT_manta=="2/2"] = NA    
                
my_data$GT = consensus_genotype(my_data,callers = callers_all)
        
## geno callers detected ####
                
my_data$geno_callers = geno_callers(my_data,callers = callers_all)
                
my_data$GT_manta = my_data$GT_manta_ok
                
                
# change GT
                
my_data$GT_manta[my_data$GT_manta=="0/1"] = "0/1-1/1"
my_data$GT_manta[my_data$GT_manta=="1/1"] = "0/1-1/1"
my_data$GT_manta[my_data$GT_manta=="2/2"] = "0/1-1/1"
                
my_data$GT_delly[my_data$GT_delly=="0/1"] = "0/1-1/1"
my_data$GT_delly[my_data$GT_delly=="1/1"] = "0/1-1/1"
                
my_data$GT_gatk[my_data$GT_gatk=="0/1"] = "0/1-1/1"
my_data$GT_gatk[my_data$GT_gatk=="1/1"] = "0/1-1/1"
                
my_data$GT_lumpy[my_data$GT_lumpy=="0/1"] = "0/1-1/1"
my_data$GT_lumpy[my_data$GT_lumpy=="1/1"] = "0/1-1/1"
                
my_data$GT_pindel[my_data$GT_pindel=="0/1"] = "0/1-1/1"
my_data$GT_pindel[my_data$GT_pindel=="1/1"] = "0/1-1/1"
                
my_data$GT_svaba[my_data$GT_svaba=="0/1"] = "0/1-1/1"
my_data$GT_svaba[my_data$GT_svaba=="1/1"] = "0/1-1/1"
                
my_data$GT_whamg[my_data$GT_whamg=="0/1"] = "0/1-1/1"
my_data$GT_whamg[my_data$GT_whamg=="1/1"] = "0/1-1/1"
                
my_data$GT_deepvariant[my_data$GT_deepvariant=="0/1"] = "0/1-1/1"
my_data$GT_deepvariant[my_data$GT_deepvariant=="1/1"] = "0/1-1/1"
              
my_data$GT_strelka[my_data$GT_strelka=="0/1"] = "0/1-1/1"
my_data$GT_strelka[my_data$GT_strelka=="1/1"] = "0/1-1/1"
              
## add number of callers detected #####
                
my_data$callers_detected = n_callers_detected(my_data,callers = callers_all)
                
my_data$callers_detected[my_data$callers_detected %in% 6:9] = "6-7-8-9"
                
my_data$callers_detected_ok = n_callers_detected(my_data,callers = callers_all)
                
### consensus length #######
                
my_data$length = consensus_length(my_data,callers = callers_all)
                                             
aux = cut(my_data$length,breaks = c(30,50,75,100,125,151),right = F)
                
my_data$length_stretch = aux
                
                
### consensus start #######
                
my_data$start = consensus_start(my_data,callers = callers_all)
                
                
## add reciprocity ###
                
my_data$reciprocity = reprocicity(my_data,callers = callers_all)
                
                
## add strategy ####
                
my_data$strategy = strategy(my_data,callers = callers_all,strategies)
                
## predict using LR model #####
                
my_data$GT_manta[is.na(my_data$GT_manta)] = "9/9"
my_data$GT_whamg[is.na(my_data$GT_whamg)] = "9/9"
my_data$GT_delly[is.na(my_data$GT_delly)] = "9/9"
my_data$GT_lumpy[is.na(my_data$GT_lumpy)] = "9/9"
my_data$GT_pindel[is.na(my_data$GT_pindel)] = "9/9"
my_data$GT_svaba[is.na(my_data$GT_svaba)] = "9/9"
my_data$GT_gatk[is.na(my_data$GT_gatk)] = "9/9"
my_data$GT_strelka[is.na(my_data$GT_strelka)] = "9/9"
my_data$GT_deepvariant[is.na(my_data$GT_deepvariant)] = "9/9"
                
                
data_predict = my_data %>% select(ID,chr,start,GT_manta,GT_whamg,
                                                  GT_delly,GT_lumpy,GT_pindel,
                                                  GT_svaba,GT_gatk,GT_strelka,GT_deepvariant,length,
                                                  length_stretch,strategy,
                                                  reciprocity,callers_detected,callers_detected_ok,GT,geno_callers) %>%
                  as.data.frame(row.names=FALSE)
                
data_predict$GT_manta = as.factor(data_predict$GT_manta)
data_predict$GT_whamg = as.factor(data_predict$GT_whamg)
data_predict$GT_delly = as.factor(data_predict$GT_delly)
data_predict$GT_lumpy = as.factor(data_predict$GT_lumpy)
data_predict$GT_pindel = as.factor(data_predict$GT_pindel)
data_predict$GT_gatk = as.factor(data_predict$GT_gatk)
data_predict$GT_deepvariant = as.factor(data_predict$GT_deepvariant)
data_predict$GT_strelka = as.factor(data_predict$GT_strelka)
data_predict$GT_svaba = as.factor(data_predict$GT_svaba)
data_predict$length_stretch = as.factor(data_predict$length_stretch)
data_predict$strategy = as.factor(data_predict$strategy)
                
# predict YES/NO threshold 0.5
                
my_pred = caret::predict.train(mod_fit,as.data.frame(data_predict)) 
                
my_pred = ifelse(my_pred=="YES",0,1)
                
my_pred = ifelse(my_pred==0,"PASS","NO_PASS")
                
data_predict$PASS = my_pred
        
# predict probability 
                
my_pred = predict.train(mod_fit,as.data.frame(data_predict),type = "prob")
                
data_predict$PASS_num = my_pred[,2]
                
data_predict = as.data.table(data_predict)
                
## name callers detected ####
                
data_predict$name_callers = name_callers_detected(data_predict,callers = callers_all)
                
## minimum window size ####
                
data_predict$min_window_size = min_window_size(data_predict,callers = callers_all,call_windows)
                
## window size for merging between samples ####
                
data_predict$lower = data_predict$start-data_predict$min_window_size
data_predict$upper = data_predict$start+data_predict$min_window_size
                
## group by ID to remove duplicates ####
                
final_data = data_predict %>% select(ID,chr,start,GT,length,
                                                     strategy,reciprocity,
                                                     callers_detected_ok,name_callers,
                                                     min_window_size,lower,upper,PASS,PASS_num,geno_callers) %>%
                  arrange(start)
                
final_data = final_data %>% group_by(ID) %>% 
                  summarise(chr = unique(chr),
                            start = floor(median(start)),
                            GT = paste(unique(na.omit(GT)), collapse = ','),
                            length = median(length),
                            strategy = paste(unique(na.omit(strategy)), collapse = ','),
                            reciprocity = median(reciprocity),
                            callers_detected = median(callers_detected_ok),
                            name_callers = paste(unique(na.omit(name_callers)), collapse = ','),
                            GT_callers = paste(unique(na.omit(geno_callers)), collapse = ','),
                            min_window_size = median(min_window_size),
                            lower = median(lower),
                            upper = median(upper),
                            PASS = paste(unique(na.omit(PASS)), collapse = ','),
                            PASS_num = median(PASS_num)) %>%
                  select(-ID)
                
final_data = final_data %>% arrange(start)
                
final_data$GT[!final_data$GT %in% c("0/1","1/1")] = "./."
                
final_data = final_data %>% filter(PASS %in% c("PASS","NO_PASS"))
                
missing_geno = which(final_data$GT=="./.")
   
final_data$GT_best = NA
                
for(k in missing_geno){
                  
                  ranking_missing = c("deepvariant","strelka","delly","manta","gatk","pindel","svaba","whamg","lumpy")
                  
                  callers = unlist(tstrsplit(final_data$name_callers[k],","))
                  
                  if(final_data$length[k]<51 & sum(callers %in% "manta")>0){
                    
                    ranking_missing = c("deepvariant","strelka","delly","gatk","pindel","svaba","whamg","lumpy")
                    callers = callers[-which(callers %in% "manta")]
                  }
                  
                  d=1
                  
                  if(length(callers)>0){
                    while(!ranking_missing[d] %in% callers){
                      
                      d = d + 1
                      
                    }
                    
                    final_data$GT_best[k] = unlist(tstrsplit(final_data$GT_callers[k],","))[which(callers %in% ranking_missing[d])]
                    print(k)}
                }
                
final_data[which(final_data[,9]==""),9] = "./."
                
manta_missing = which((as.character(t(final_data[,8])) %in% "manta") & as.character(t(final_data[,3]))=="./.") 
                
final_data[manta_missing,15] = "0/0"
                
colnames(final_data) = paste0(colnames(final_data),"_",ids[i])
                
fwrite(final_data,paste0("/2_merge_callers/mid_DEL/outputs/",ids[i],"/",ids[i],"_mid_DEL_chr_",j),
                       sep = " ",row.names = F,quote = F)
                          
