library(data.table)
library(dplyr)
library(caret)
library(e1071)

args <- commandArgs(trailingOnly = TRUE)
i= as.numeric(args[2])
              
source("ext/functions.R")
              
strategies = fread("ext/strategies.csv")
              
mod_fit = readRDS("1_LRM/models/LRM_model_INS.rds")
              
ids = fread("ext/all_samplesok",header = F)
              
ids = ids$V1

print(ids[i])
              
dir.create(paste0("/2_merge_callers/INS/outputs/",ids[i]))
              

# read VCF files ########

delly = fread(paste0("/2_merge_callers/INS/data/Delly/",ids[i],"_INS"))
delly$V3 = as.numeric(abs(delly$V3))
delly = delly %>% filter(V3>30)
colnames(delly) = c("chr","start_delly","length_delly","GT_delly")
delly$lower_delly = delly$start_delly-10
delly$upper_delly = delly$start_delly+10
delly$chr_pos_delly = do.call(paste0,list(delly$chr,"_",delly$start_delly))
delly = delly %>% filter(GT_delly!="0/0")
delly$chr[delly$chr=="X"] = 23
delly$chr[delly$chr=="Y"] = 24
delly$chr = as.character(delly$chr)

popins = fread(paste0("/2_merge_callers/INS/data/Popins/",ids[i],"_INS"))
colnames(popins) = c("chr","start_popins","GT_popins")
popins$lower_popins = popins$start_popins-10
popins$upper_popins = popins$start_popins+10
popins$chr_pos_popins = do.call(paste0,list(popins$chr,"_",popins$start_popins))
popins$length_popins = 0
popins = popins %>% arrange(chr,start_popins)
popins$start_popins = as.numeric(popins$start_popins)
#popins$diff = c(0,diff(popins$start_popins))
#popins = popins %>% filter(diff<0 | diff>200)
popins = popins %>% filter(GT_popins!="0/0")
popins$chr[popins$chr=="X"] = 23
popins$chr[popins$chr=="Y"] = 24
popins$chr = as.character(popins$chr)

pindel = fread(paste0("/2_merge_callers/INS/data/Pindel/",ids[i],"_INS"))
pindel$V4 = as.numeric(nchar(pindel$V4))
pindel$V3 <- NULL
pindel$V6 = pindel$V4
pindel$V4 <- NULL
colnames(pindel) = c("chr","start_pindel","GT_pindel","length_pindel")
pindel2=fread(paste0("/gpfs/projects/bsc05/jordivalls/GCAT_project_all_samples/merge_all_calling_GCAT/Insertions/Pindel/",ids[i],"_Pindel_insertions_big"))
colnames(pindel2) = c("chr","start_pindel","GT_pindel")
pindel2$length_pindel = 0
pindel = rbind(pindel,pindel2)
pindel$lower_pindel = pindel$start_pindel-10
pindel$upper_pindel = pindel$start_pindel+10
pindel$chr_pos_pindel = do.call(paste0,list(pindel$chr,"_",pindel$start_pindel))
pindel$chr = as.character(pindel$chr)
pindel = pindel %>% filter(GT_pindel!="0/0")
pindel$chr[pindel$chr=="X"] = 23
pindel$chr[pindel$chr=="Y"] = 24
pindel$chr = as.character(pindel$chr)

whamg = fread(paste0("/2_merge_callers/INS/data/Whamg/",ids[i],"_INS"))
whamg = whamg[,-2,with=F]
colnames(whamg) = c("chr","start_whamg","length_whamg")
whamg$lower_whamg = whamg$start_whamg-10
whamg$upper_whamg = whamg$start_whamg+10
whamg$chr_pos_whamg = do.call(paste0,list(whamg$chr,"_",whamg$start_whamg))
whamg$GT_whamg = "0/1"
whamg = whamg %>% arrange(chr,chr_pos_whamg)
whamg$start_whamg = as.numeric(whamg$start_whamg)
whamg$diff = c(0,diff(whamg$start_whamg))
whamg = whamg %>% filter(diff<0 | diff>10)
whamg = whamg %>% filter(GT_whamg!="0/0")
whamg$chr[whamg$chr=="X"] = 23
whamg$chr[whamg$chr=="Y"] = 24
whamg$chr = as.character(whamg$chr)

svaba = fread(paste0("/2_merge_callers/INS/data/SVaBA/",ids[i],"_INS"))
colnames(svaba) = c("chr","start_svaba","chr2","pos2","GT_svaba")

# remove translocations and change positions

remove_variants = NULL

for(variant in 1:nrow(svaba)){
  
  if(svaba$chr[variant]==svaba$chr2[variant] & svaba$start_svaba[variant]>svaba$pos2[variant]){
    
    start = svaba$start_svaba[variant]
    end = svaba$pos2[variant]
    
    svaba$start_svaba[variant] = end
    svaba$pos2[variant] = start
  }
  if(svaba$chr[variant]!=svaba$chr2[variant]){
    
    remove_variants = c(remove_variants,variant)
    
  }
}

svaba = svaba[-remove_variants,]

svaba$chr2 <- NULL
svaba$pos2 <- NULL
svaba$length_svaba = 0
svaba1 = fread(paste0("/2_merge_callers/INS/data/SVaBA/",ids[i],"_INS_small"))
svaba1$length_svaba = abs(nchar(svaba1$V3)-nchar(svaba1$V4))
colnames(svaba1)[c(1,2,5)] = c("chr","start_svaba","GT_svaba")
svaba1 = svaba1[,c(1,2,5,6)]
svaba = rbind(svaba,svaba1)
svaba$lower_svaba = svaba$start_svaba-10
svaba$upper_svaba = svaba$start_svaba+10
svaba$chr_pos_svaba = do.call(paste0,list(svaba$chr,"_",svaba$start_svaba))
svaba = svaba %>% arrange(chr,start_svaba)
svaba$start_svaba = as.numeric(svaba$start_svaba)
svaba$diff = c(0,diff(svaba$start_svaba))
svaba = svaba %>% filter(diff<0 | diff>10)
svaba = svaba %>% filter(GT_svaba!="0/0")
svaba$chr[svaba$chr=="X"] = 23
svaba$chr[svaba$chr=="Y"] = 24
svaba$chr = as.character(svaba$chr)



manta = fread(paste0("/2_merge_callers/INS/data/Manta/",ids[i],"_INS_small"))
manta$V4 = as.numeric(abs(manta$V4))
manta$V3 = as.character(manta$V3)
manta$V3 = "0/1"
colnames(manta) = c("chr","start_manta","GT_manta","length_manta")
manta2 = fread(paste0("/2_merge_callers/INS/data/Manta/",ids[i],"_INS_large"))
colnames(manta2) = c("chr","start_manta","GT_manta")
manta2$length_manta = 0
manta = rbind(manta,manta2)
manta$lower_manta = manta$start_manta-10
manta$upper_manta = manta$start_manta+10
manta$chr_pos_manta = do.call(paste0,list(manta$chr,"_",manta$start_manta))
manta = manta %>% arrange(chr,start_manta)
manta$start_manta = as.numeric(manta$start_manta)
manta$diff = c(0,diff(manta$start_manta))
manta = manta %>% filter(diff<0 | diff>10)
manta = manta %>% filter(GT_manta!="0/0")
manta$chr[manta$chr=="X"] = 23
manta$chr[manta$chr=="Y"] = 24
manta$chr = as.character(manta$chr)

gatk = fread(paste0("/2_merge_callers/INS/data/Gatk/",ids[i],"_INS"))
gatk$length_gatk = abs(nchar(gatk$V3)-nchar(gatk$V4))
#gatk = gatk %>% filter(length_gatk>30 & length_gatk<151)
colnames(gatk)[c(1,2,5)] = c("chr","start_gatk","GT_gatk")
gatk = gatk[,c(1,2,5,6)]
gatk$end = gatk$length_gatk+gatk$start_gatk
gatk$lower_gatk = gatk$start_gatk-10
gatk$upper_gatk = gatk$start_gatk+10
gatk$chr_pos_gatk = do.call(paste0,list(gatk$chr,"_",gatk$start_gatk))
gatk = gatk %>% filter(GT_gatk!="0/0")
gatk$chr[gatk$chr=="X"] = 23
gatk$chr[gatk$chr=="Y"] = 24
gatk$chr = as.character(gatk$chr)

strelka = fread(paste0("/2_merge_callers/INS/data/Strelka/",ids[i],"_INS"))
strelka$length_strelka = abs(nchar(strelka$V3)-nchar(strelka$V4))
strelka = strelka %>% filter(length_strelka>30 & length_strelka<151)
colnames(strelka)[c(1,2,5)] = c("chr","start_strelka","GT_strelka")
strelka = strelka[,c(1,2,5,6)]
strelka$end = strelka$length_strelka+strelka$start_strelka
strelka$lower_strelka = strelka$start_strelka-10
strelka$upper_strelka = strelka$start_strelka+10
strelka$chr_pos_strelka = do.call(paste0,list(strelka$chr,"_",strelka$start_strelka))
strelka = strelka %>% filter(GT_strelka!="0/0")
strelka$chr[strelka$chr=="X"] = 23
strelka$chr[strelka$chr=="Y"] = 24
strelka$chr = as.character(strelka$chr)


call_windows = data.frame(caller = c("popins","delly","pindel","whamg",
                                     "manta","gatk","svaba", "strelka"),
                          window = c(10,10,10,10,10,10,10,10))

call_windows$caller = as.character(call_windows$caller)



# filter by chromosome ########
 
j= as.numeric(args[1])

delly_chr = delly %>% filter(chr==j) %>% arrange(start_delly) %>% as.data.table() %>% unique()
popins_chr = popins %>% filter(chr==j) %>% arrange(start_popins) %>% as.data.table() %>% unique()
pindel_chr = pindel %>% filter(chr==j) %>% arrange(start_pindel) %>% as.data.table() %>% unique()
whamg_chr = whamg %>% filter(chr==j) %>% arrange(start_whamg) %>% as.data.table() %>% unique()
svaba_chr = svaba %>% filter(chr==j) %>% arrange(start_svaba) %>% as.data.table() %>% unique()
manta_chr = manta %>% filter(chr==j) %>% arrange(start_manta) %>% as.data.table() %>% unique()
gatk_chr = gatk %>% filter(chr==j) %>% arrange(start_gatk) %>% as.data.table() %>% unique()
strelka_chr = strelka %>% filter(chr==j) %>% arrange(start_strelka) %>% as.data.table() %>% unique()

n_rows = nrow(delly_chr) + nrow(popins_chr) + nrow(pindel_chr) + nrow(whamg_chr) + nrow(svaba_chr) + nrow(manta_chr) + nrow(gatk_chr) + nrow(strelka_chr)

if(n_rows!=0){

## prepare BBDD ordered by F-score ####

f_score_order = c("manta","pindel","popins","whamg",
                  "gatk","strelka","delly","svaba")

manta_chr$ID = paste0("insertion_",1:nrow(manta_chr)) # el merge empieza por popins y asignamos IDs
delly_chr$ID = "none"
pindel_chr$ID = "none"
whamg_chr$ID = "none"
popins_chr$ID = "none"
gatk_chr$ID = "none"
svaba_chr$ID = "none"
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

my_data = merge_callers_insertions(get(paste0(second_call,"_chr")),get(paste0(first_call,"_chr")),callers_ref = first_call,
                        caller_to_merge = second_call) 

for(caller_merge in (caller_second+1):length(f_score_order)){
  
  n_row_call = nrow(get(paste0(f_score_order[caller_merge],"_chr")))
  
  while(n_row_call!=0){
      callers_ref = tstrsplit(colnames(my_data)[grep("GT_",colnames(my_data))],split="GT_")[[2]]
    
      my_data = merge_callers_insertions(get(paste0(f_score_order[caller_merge],"_chr")),my_data,callers_ref = callers_ref,
                                            caller_to_merge = f_score_order[caller_merge]) 
      
      n_row_call = 0
            
    }
    
  }
  
  
callers_all = tstrsplit(colnames(my_data)[grep("GT_",colnames(my_data))],split="GT_")[[2]]
  
caller = c("delly","popins","manta","strelka","pindel","gatk","svaba","whamg")

if(sum(caller %in% callers_all)!=8){
  
  callers_missing = caller[-which(caller %in% callers_all)]
  
  for(call_miss in 1:length(callers_missing)){
    
    my_data$start_miss = NA
    my_data$GT_miss = NA
    my_data$length_miss = NA
    
    colnames(my_data)[(ncol(my_data)-2):ncol(my_data)] = paste0(c("start_","GT_","length_"),callers_missing[call_miss])
    
  }
  
  
}
## geno callers detected ####
  
my_data$geno_callers = geno_callers(my_data,
                                      callers = callers_all)
  
# change GT
  
my_data$GT_manta2 = my_data$GT_manta
my_data$GT_whamg2 = my_data$GT_whamg
my_data$GT_delly2 = my_data$GT_delly
my_data$GT_popins2 = my_data$GT_popins
my_data$GT_pindel2 = my_data$GT_pindel
my_data$GT_svaba2 = my_data$GT_svaba
my_data$GT_gatk2 = my_data$GT_gatk
my_data$GT_strelka2 = my_data$GT_strelka

my_data$GT_manta[my_data$GT_manta=="0/1"] = "0/1-1/1"
my_data$GT_manta[my_data$GT_manta=="1/1"] = "0/1-1/1"

my_data$GT_gatk[my_data$GT_gatk=="0/1"] = "0/1-1/1"
my_data$GT_gatk[my_data$GT_gatk=="1/1"] = "0/1-1/1"

my_data$GT_delly[my_data$GT_delly=="0/1"] = "0/1-1/1"
my_data$GT_delly[my_data$GT_delly=="1/1"] = "0/1-1/1"

my_data$GT_popins[my_data$GT_popins=="0/1"] = "0/1-1/1"
my_data$GT_popins[my_data$GT_popins=="1/1"] = "0/1-1/1"

my_data$GT_pindel[my_data$GT_pindel=="0/1"] = "0/1-1/1"
my_data$GT_pindel[my_data$GT_pindel=="1/1"] = "0/1-1/1"

my_data$GT_svaba[my_data$GT_svaba=="0/1"] = "0/1-1/1"
my_data$GT_svaba[my_data$GT_svaba=="1/1"] = "0/1-1/1"
  
my_data$GT_whamg[my_data$GT_whamg=="0/1"] = "0/1-1/1"
my_data$GT_whamg[my_data$GT_whamg=="1/1"] = "0/1-1/1"
  
my_data$GT_strelka[my_data$GT_strelka=="0/1"] = "0/1-1/1"
my_data$GT_strelka[my_data$GT_strelka=="1/1"] = "0/1-1/1"
  
  
  ## add number of callers detected #####
  
my_data$callers_detected = n_callers_detected(my_data,
                                                callers = callers_all)
  
my_data$callers_detected_ok = n_callers_detected(my_data,
                                                callers = callers_all)

my_data$callers_detected[my_data$callers_detected=="4"] = "4-5-6-7"
my_data$callers_detected[my_data$callers_detected=="5"] = "4-5-6-7"
my_data$callers_detected[my_data$callers_detected=="6"] = "4-5-6-7"
my_data$callers_detected[my_data$callers_detected=="7"] = "4-5-6-7"
my_data$callers_detected[my_data$callers_detected=="8"] = "4-5-6-7"
  
  ### consensus start #######
  
my_data$start = consensus_start_insertions(my_data,
                                  callers = callers_all)
  
  ## add strategy ####
  
my_data$strategy_ok = strategy(my_data,
                                 callers = callers_all,
                                 strategies)
  
table(my_data$callers_detected)

  
  ## predict using LR model #####
  
my_data$GT_manta[is.na(my_data$GT_manta)] = "9/9"
my_data$GT_whamg[is.na(my_data$GT_whamg)] = "9/9"
my_data$GT_delly[is.na(my_data$GT_delly)] = "9/9"
my_data$GT_popins[is.na(my_data$GT_popins)] = "9/9"
my_data$GT_pindel[is.na(my_data$GT_pindel)] = "9/9"
my_data$GT_svaba[is.na(my_data$GT_svaba)] = "9/9"
my_data$GT_gatk[is.na(my_data$GT_gatk)] = "9/9"
my_data$GT_strelka[is.na(my_data$GT_strelka)] = "9/9"
my_data$GT_whamg[my_data$GT_whamg==""] = "9/9"
  
data_predict = my_data %>% select(ID,chr,start,GT_manta,GT_whamg,
                                    GT_delly,GT_popins,GT_pindel,
                                    GT_gatk,GT_svaba,GT_strelka,GT_manta2,GT_whamg2,
                                    GT_delly2,GT_popins2,GT_pindel2,
                                    GT_gatk2,GT_svaba2,GT_strelka2,
                                    length_manta,
                                    strategy_ok,callers_detected,callers_detected_ok,geno_callers) %>%
  as.data.frame(row.names=FALSE)
  
data_predict$GT_manta = as.factor(data_predict$GT_manta)
data_predict$GT_whamg = as.factor(data_predict$GT_whamg)
data_predict$GT_delly = as.factor(data_predict$GT_delly)
data_predict$GT_popins = as.factor(data_predict$GT_popins)
data_predict$GT_pindel = as.factor(data_predict$GT_pindel)
data_predict$GT_gatk = as.factor(data_predict$GT_gatk)
data_predict$GT_svaba = as.factor(data_predict$GT_svaba)
data_predict$GT_strelka = as.factor(data_predict$GT_strelka)
data_predict$callers_detected= as.factor(data_predict$callers_detected)
data_predict$geno_callers=as.factor(data_predict$geno_callers)
  
  # predict YES/NO threshold 0.5
  
my_pred = caret::predict.train(mod_fit,as.data.frame(data_predict)) 
my_pred = ifelse(my_pred=="YES",0,1)
  
my_pred = ifelse(my_pred==0,"PASS","NO_PASS")
  
data_predict$PASS = my_pred
  
  # predict probability (valores entre 0 y 1)
  
my_pred = predict.train(mod_fit,as.data.frame(data_predict),type = "prob")
  
data_predict$PASS_num = my_pred[,2]
  
data_predict = as.data.table(data_predict)
  
  ## name callers detected ####
  
data_predict$name_callers = name_callers_detected(data_predict,
                                                    callers = callers_all)
  
  ## minimum window size ####
  
data_predict$min_window_size = min_window_size(data_predict,
                                                 callers = callers_all,
                                                 call_windows)
  
  ## window size for merging between samples ####
  
data_predict$lower = data_predict$start-data_predict$min_window_size
data_predict$upper = data_predict$start+data_predict$min_window_size
  
  ## consensus genotype #####
  
data_predict$GT_manta = data_predict$GT_manta2
data_predict$GT_whamg = data_predict$GT_whamg2
data_predict$GT_delly = data_predict$GT_delly2
data_predict$GT_popins = data_predict$GT_popins2
data_predict$GT_pindel = data_predict$GT_pindel2
data_predict$GT_svaba = data_predict$GT_svaba2
data_predict$GT_gatk = data_predict$GT_gatk2
data_predict$GT_strelka = data_predict$GT_strelka2
  
data_predict$GT_whamg[data_predict$GT_whamg=="0/1"] = NA # remove whamg
  
data_predict$GT_manta[which(data_predict$length_manta!=0)] = NA
  
data_predict$GT = consensus_genotype(data_predict,
                                  callers = callers_all)
  
  
data_predict$geno_callers = geno_callers(data_predict,
                                            callers = callers_all)
  
  ## group by ID to remove duplicates ####
  
final_data = data_predict %>% select(ID,chr,start,GT,
                                       strategy_ok,
                                       callers_detected_ok,name_callers,
                                       min_window_size,lower,upper,PASS,PASS_num,geno_callers,length_manta) %>%
  arrange(start)
  
final_data = final_data %>% group_by(ID) %>% 
  summarise(chr = unique(chr),
              start = median(start),
              GT = paste(unique(na.omit(GT)), collapse = ','),
              strategy = paste(unique(na.omit(strategy_ok)), collapse = ','),
              callers_detected = median(callers_detected_ok),
              name_callers = paste(unique(na.omit(name_callers)), collapse = ','),
              GT_callers = paste(unique(na.omit(geno_callers)), collapse = ','),
              min_window_size = median(min_window_size),
              lower = median(lower),
              upper = median(upper),
              PASS = paste(unique(na.omit(PASS)), collapse = ','),
              PASS_num = median(PASS_num),
              length_manta = median(length_manta)) %>%
  select(-ID)
  
final_data = final_data %>% arrange(start)

final_data$GT[!final_data$GT %in% c("0/1","1/1")] = "./."
final_data = final_data %>% filter(PASS %in% c("PASS","NO_PASS")) # treure rares "PASS,NO_PASS"
missing_geno = which(final_data$GT=="./.")
final_data$GT_best = NA
  
final_data[missing_geno,]
  
if(length(missing_geno)>0){
for(k in missing_geno){
    
  ranking_missing = c("delly","popins","manta","strelka","pindel","gatk","svaba")
    
  if(!is.na(final_data$length_manta[k]) & final_data$length_manta[k]>0){ranking_missing = c("delly","popins","strelka","pindel","gatk","svaba")}
    
  callers = unlist(tstrsplit(final_data$name_callers[k],","))
  
  d=1
  while(!ranking_missing[d] %in% callers){
    d = d + 1
    if(is.na(ranking_missing[d]) == TRUE){
	break
    }
  }
  if(is.na(ranking_missing[d]) == FALSE){
  final_data$GT_best[k] = unlist(tstrsplit(final_data$GT_callers[k],","))[which(callers %in% ranking_missing[d])]
  }
  else{
    final_data$GT_best[k] = "0/0"
    print (k)
    print (callers)
    print(final_data$GT_best[k])
     }

}
}
  

final_data$length_manta = NULL
  
colnames(final_data) = paste0(colnames(final_data),"_",ids[i])
  
fwrite(final_data,paste0("/2_merge_callers/INS/outputs/",ids[i],"/",ids[i],"_INS_chr_",j),
                       sep = " ",row.names = F,quote = F)
