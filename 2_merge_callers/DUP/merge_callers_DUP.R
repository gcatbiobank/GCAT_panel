library(data.table)
library(dplyr)
library(caret)
library(e1071)

args <- commandArgs(trailingOnly = TRUE)
i= as.numeric(args[2])
              
source("ext/functions.R")
              
strategies = fread("ext/strategies.csv")
              
mod_fit = readRDS("1_LRM/models/LRM_model_DUP.rds")
              
ids = fread("ext/all_samplesok",header = F)
              
ids = ids$V1

print(ids[i])
              
dir.create(paste0("/2_merge_callers/DUP/outputs/",ids[i]))
              

# read VCF files ########

delly = fread(paste0("/2_merge_callers/DUP/data/Delly/",ids[i],"_DUP"))
delly = delly[,1:5]
delly$V4 = as.numeric(abs(delly$V4))
colnames(delly) = c("chr","start_delly","end","length_delly","GT_delly")
delly$lower_delly = delly$start_delly-10
delly$upper_delly = delly$start_delly+10
delly$chr_pos_delly = do.call(paste0,list(delly$chr,"_",delly$start_delly))
delly$chr[delly$chr=="X"] = 23
delly$chr[delly$chr=="Y"] = 24
delly$chr = as.character(delly$chr)
delly = delly %>% filter(chr %in% 1:23)
delly$GT2_delly = delly$GT_delly
delly = delly %>% filter(length_delly>30)


table(delly$chr)
dim(delly)

lumpy = fread(paste0("/2_merge_callers/DUP/data/Lumpy/",ids[i],"_DUP"))
lumpy$V4 = as.numeric(abs(lumpy$V4))
lumpy = lumpy [,1:5]
colnames(lumpy) = c("chr","start_lumpy","end","length_lumpy","GT_lumpy")
lumpy$lower_lumpy = lumpy$start_lumpy-50
lumpy$upper_lumpy = lumpy$start_lumpy+50
lumpy$chr_pos_lumpy = do.call(paste0,list(lumpy$chr,"_",lumpy$start_lumpy))
lumpy$chr[lumpy$chr=="X"] = 23
lumpy$chr[lumpy$chr=="Y"] = 24
lumpy$chr = as.character(lumpy$chr)
lumpy = lumpy %>% filter(chr %in% 1:23)
lumpy$GT2_lumpy= lumpy$GT_lumpy
lumpy = lumpy %>% filter(length_lumpy>30)


pindel = fread(paste0("/2_merge_callers/DUP/data/Pindel/",ids[i],"_DUP"))
pindel = pindel [,1:5]
colnames(pindel) = c("chr","start_pindel","end","length_pindel","GT_pindel")
pindel$lower_pindel = pindel$start_pindel-10
pindel$upper_pindel = pindel$start_pindel+10
pindel$chr_pos_pindel = do.call(paste0,list(pindel$chr,"_",pindel$start_pindel))
pindel$chr[pindel$chr=="X"] = 23
pindel$chr[pindel$chr=="Y"] = 24
pindel$chr = as.character(pindel$chr)
pindel = pindel %>% filter(chr %in% 1:23)
pindel$GT2_pindel = pindel$GT_pindel
table(pindel$chr)
dim(pindel)
pindel = pindel %>% filter(length_pindel>30)

whamg = fread(paste0("/2_merge_callers/DUP/data/Whamg/",ids[i],"_DUP"))
whamg = whamg [,1:5]
whamg$V4 = as.numeric(abs(whamg$V4))
colnames(whamg) = c("chr","start_whamg","end","length_whamg","GT_whamg")
whamg$lower_whamg = whamg$start_whamg-10
whamg$upper_whamg = whamg$start_whamg+10
whamg$chr_pos_whamg = do.call(paste0,list(whamg$chr,"_",whamg$start_whamg))
whamg$chr[whamg$chr=="X"] = 23
whamg$chr[whamg$chr=="Y"] = 24
whamg$chr = as.character(whamg$chr)
whamg = whamg %>% filter(chr %in% 1:23)
whamg$GT2_whamg = whamg$GT_whamg
whamg = whamg %>% filter(length_whamg>30)


svaba = fread(paste0("/2_merge_callers/DUP/data/SVaBA/",ids[i],"_DUP"))
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

svaba$lower_svaba = svaba$start_svaba-10
svaba$upper_svaba = svaba$start_svaba+10
svaba$chr_pos_svaba = do.call(paste0,list(svaba$chr,"_",svaba$start_svaba))
svaba$pos2 = as.numeric(svaba$pos2)
svaba$length_svaba = svaba$pos2-svaba$start_svaba
svaba$chr[svaba$chr=="X"] = 23
svaba$chr[svaba$chr=="Y"] = 24
svaba$chr = as.character(svaba$chr)
svaba = svaba %>% filter(chr %in% 1:23)
svaba$GT2_svaba = svaba$GT_svaba
svaba = svaba %>% filter(length_svaba>30)


manta = fread(paste0("/2_merge_callers/DUP/data/Manta/",ids[i],"_DUP"))
manta$V4 = as.numeric(abs(manta$V4))
manta= manta [,1:5]
colnames(manta) = c("chr","start_manta","end","length_manta","GT_manta")
manta$lower_manta = manta$start_manta-20
manta$upper_manta = manta$start_manta+20
manta$chr_pos_manta = do.call(paste0,list(manta$chr,"_",manta$start_manta))
manta$chr[manta$chr=="X"] = 23
manta$chr[manta$chr=="Y"] = 24
manta$chr = as.character(manta$chr)
manta = manta %>% filter(chr %in% 1:23)
manta$GT2_manta = manta$GT_manta
dim (manta)
table(manta$chr)
manta = manta %>% filter(length_manta>30)

cnvnator = fread(paste0("/2_merge_callers/DUP/data/CNVnator/",ids[i],"_DUP"))
cnvnator$V4 = as.numeric(abs(cnvnator$V4))
colnames(cnvnator) = c("chr","start_cnvnator","end","length_cnvnator","GT_cnvnator")
cnvnator$lower_cnvnator = cnvnator$start_cnvnator-100
cnvnator$upper_cnvnator = cnvnator$start_cnvnator+100
cnvnator$chr_pos_cnvnator = do.call(paste0,list(cnvnator$chr,"_",cnvnator$start_cnvnator))
cnvnator$chr[cnvnator$chr=="X"] = 23 
cnvnator$chr[cnvnator$chr=="Y"] = 24
cnvnator$chr = as.character(cnvnator$chr)
cnvnator$GT2_cnvnator = cnvnator$GT_cnvnator
dim (cnvnator)
cnvnator = cnvnator %>% filter(length_cnvnator>30)

call_windows = data.frame(caller = c("delly","lumpy","pindel","whamg",
                                     "svaba","manta","cnvnator"),
                          window = c(10,50,10,10,10,20,100))

call_windows$caller = as.character(call_windows$caller)


# filter by chromosome ########
        
j= as.numeric(args[1])

        delly_chr = delly %>% filter(chr==j) %>% arrange(start_delly) %>% as.data.table() %>% unique()
        lumpy_chr = lumpy %>% filter(chr==j) %>% arrange(start_lumpy) %>% as.data.table() %>% unique()
        pindel_chr = pindel %>% filter(chr==j) %>% arrange(start_pindel) %>% as.data.table() %>% unique()
        whamg_chr = whamg %>% filter(chr==j) %>% arrange(start_whamg) %>% as.data.table() %>% unique()
        svaba_chr = svaba %>% filter(chr==j) %>% arrange(start_svaba) %>% as.data.table() %>% unique()
        manta_chr = manta %>% filter(chr==j) %>% arrange(start_manta) %>% as.data.table() %>% unique()
        cnvnator_chr = cnvnator %>% filter(chr==j) %>% arrange(start_cnvnator) %>% as.data.table() %>% unique()
        
        n_rows = nrow(delly_chr) + nrow(lumpy_chr) + nrow(pindel_chr) + nrow(whamg_chr) + nrow(svaba_chr) + nrow(manta_chr) + nrow(cnvnator_chr)
        
        if(n_rows!=0){
        
        ## prepare BBDD ordered by F-score ####
        
        f_score_order = c("lumpy","manta","pindel","delly",
                          "whamg","svaba","cnvnator")
        
        lumpy_chr$ID = paste0("duplication_",1:nrow(lumpy_chr)) # el merge empieza por lumpy y asignamos IDs
        manta_chr$ID = "none"
        pindel_chr$ID = "none"
        delly_chr$ID = "none"
        whamg_chr$ID = "none"
        svaba_chr$ID = "none"
        cnvnator_chr$ID = "none"
        
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
                                caller_to_merge = second_call,repro=0.8,svtype="duplication") 
        
        for(caller_merge in (caller_second+1):length(f_score_order)){
          
          n_row_call = nrow(get(paste0(f_score_order[caller_merge],"_chr")))
          
          while(n_row_call!=0){
      
            callers_ref = tstrsplit(colnames(my_data)[grep("GT_",colnames(my_data))],split="GT_")[[2]]
            
            my_data = merge_callers_deletions_duplications(get(paste0(f_score_order[caller_merge],"_chr")),my_data,callers_ref = callers_ref,
                                                  caller_to_merge = f_score_order[caller_merge],repro=0.8,svtype="duplication") 
            
            n_row_call = 0
                  
          }
          
        }
        
        
        callers_all = tstrsplit(colnames(my_data)[grep("GT_",colnames(my_data))],split="GT_")[[2]]
        
        
	 callers = c("manta","whamg","delly","lumpy","pindel","svaba","cnvnator")

        if(sum(callers %in% callers_all)!=7){

          callers_missing = callers[-which(callers %in% callers_all)]

          for(call_miss in 1:length(callers_missing)){

            my_data$start_miss = NA
            my_data$end_miss = NA
            my_data$GT_miss = NA
            my_data$length_miss = NA

            colnames(my_data)[(ncol(my_data)-3):ncol(my_data)] = paste0(c("start_","end_","GT_","length_"),callers_missing[call_miss])

          }


        }



        ## geno callers detected ####
        
        my_data$geno_callers = geno_callers(my_data,
                                            callers = callers_all)
        
        
        # change GT
        
        my_data$GT2_manta = my_data$GT_manta
        my_data$GT2_cnvnator = my_data$GT_cnvnator
        my_data$GT2_delly = my_data$GT_delly
        my_data$GT2_lumpy = my_data$GT_lumpy
        my_data$GT2_pindel = my_data$GT_pindel
        my_data$GT2_svaba = my_data$GT_svaba
        my_data$GT2_whamg = my_data$GT_whamg
        
        my_data$GT_manta[my_data$GT_manta=="0/1"] = "0/1-1/1"
        my_data$GT_manta[my_data$GT_manta=="1/1"] = "0/1-1/1"
        
        my_data$GT_cnvnator[my_data$GT_cnvnator=="0/1"] = "0/1-1/1"
        my_data$GT_cnvnator[my_data$GT_cnvnator=="1/1"] = "0/1-1/1"
        
        my_data$GT_delly[my_data$GT_delly=="0/1"] = "0/1-1/1"
        my_data$GT_delly[my_data$GT_delly=="1/1"] = "0/1-1/1"
        
        my_data$GT_lumpy[my_data$GT_lumpy=="0/1"] = "0/1-1/1"
        my_data$GT_lumpy[my_data$GT_lumpy=="1/1"] = "0/1-1/1"
        
        my_data$GT_pindel[my_data$GT_pindel=="0/1"] = "0/1-1/1"
        my_data$GT_pindel[my_data$GT_pindel=="1/1"] = "0/1-1/1"
        
        my_data$GT_svaba[my_data$GT_svaba=="0/1"] = "0/1-1/1"
        my_data$GT_svaba[my_data$GT_svaba=="1/1"] = "0/1-1/1"
        
        my_data$GT_whamg[my_data$GT_whamg=="0/1"] = "0/1-1/1"
        my_data$GT_whamg[my_data$GT_whamg=="1/1"] = "0/1-1/1"
        
        
        ## add number of callers detected #####
        
        my_data$callers_detected = n_callers_detected(my_data,
                                                      callers = callers_all)
        
        my_data$callers_detected2 = n_callers_detected(my_data,
                                                      callers = callers_all)
        
        my_data$callers_detected[my_data$callers_detected=="4"] = "4-5-6-7"
        my_data$callers_detected[my_data$callers_detected=="5"] = "4-5-6-7"
        my_data$callers_detected[my_data$callers_detected=="6"] = "4-5-6-7"
        my_data$callers_detected[my_data$callers_detected=="7"] = "4-5-6-7"
        
        ### consensus length #######
        
        my_data$length = consensus_length(my_data,
                                          callers = callers_all)
        
        
        my_data$length[is.na(my_data$length)] = my_data$length_cnvnator[is.na(my_data$length)]
        
        aux = cut(my_data$length,
                  breaks = c(30,150,500,1000,2000,3000,Inf),right = F)
        
        my_data$length_stretch = aux
        
        
        ### consensus start we remove CNVator because it reports bad the start position #######
        
        callers_all_no_cnvnator = callers_all
        
        if(sum(callers_all %in% "cnvnator")){
                callers_all_no_cnvnator = callers_all[-which(callers_all %in% "cnvnator")]
        }
        
        my_data$start = consensus_start_duplications(my_data,
                                                    callers = callers_all_no_cnvnator)
        
	my_data$start[which(is.na(my_data$start))] = my_data$start_cnvnator[which(is.na(my_data$start))]

        ## add reciprocity ###
        
        my_data$reciprocity = reprocicity(my_data,
                                          callers = callers_all)
        
        
        ## add strategy ####
        
        my_data$strategy_ok = strategy(my_data,
                                       callers = callers_all,
                                       strategies)
        
        my_data$strategy = my_data$strategy_ok
        
        my_data$strategy[my_data$strategy=="3"] = "3-4-5"
        my_data$strategy[my_data$strategy=="4"] = "3-4-5"
        my_data$strategy[my_data$strategy=="5"] = "3-4-5"
       
        
        ## predict using LR model #####
        
        my_data$GT_manta[is.na(my_data$GT_manta)] = "9/9"
        my_data$GT_whamg[is.na(my_data$GT_whamg)] = "9/9"
        my_data$GT_delly[is.na(my_data$GT_delly)] = "9/9"
        my_data$GT_lumpy[is.na(my_data$GT_lumpy)] = "9/9"
        my_data$GT_pindel[is.na(my_data$GT_pindel)] = "9/9"
        my_data$GT_svaba[is.na(my_data$GT_svaba)] = "9/9"
        my_data$GT_cnvnator[is.na(my_data$GT_cnvnator)] = "9/9"
        
        
        data_predict = my_data %>% select(ID,chr,start,GT_manta,GT_whamg,
                                          GT_delly,GT_lumpy,GT_pindel,
                                          GT_cnvnator,GT_svaba,length,
                                          length_stretch,strategy,strategy_ok,
                                          reciprocity,callers_detected,callers_detected2,geno_callers,GT2_manta,GT2_whamg,
                                          GT2_delly,GT2_lumpy,GT2_pindel,
                                          GT2_cnvnator,GT2_svaba) %>%
          as.data.frame(row.names=FALSE)
        
        data_predict$GT_manta = as.factor(data_predict$GT_manta)
        data_predict$GT_whamg = as.factor(data_predict$GT_whamg)
        data_predict$GT_delly = as.factor(data_predict$GT_delly)
        data_predict$GT_lumpy = as.factor(data_predict$GT_lumpy)
        data_predict$GT_pindel = as.factor(data_predict$GT_pindel)
        data_predict$GT_cnvnator = as.factor(data_predict$GT_cnvnator)
        data_predict$GT_svaba = as.factor(data_predict$GT_svaba)
        data_predict$length_stretch = as.factor(data_predict$length_stretch)
        data_predict$strategy = as.factor(data_predict$strategy)
        data_predict$callers_detected = as.factor(data_predict$callers_detected)
        
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
        
        ## consensus Genotype ####
 
	data_predict$GT = "0/1"


      
  ## group by ID to remove duplicates ####
  
  final_data = data_predict %>% select(ID,chr,start,GT,length,
                                       strategy_ok,reciprocity,
                                       callers_detected2,name_callers,
                                       min_window_size,lower,upper,PASS,PASS_num,geno_callers) %>%
    arrange(start)
  
  final_data = final_data %>% group_by(ID) %>% 
    summarise(chr = unique(chr),
              start = median(start),
              GT = paste(unique(na.omit(GT)), collapse = ','),
              length = median(length),
              strategy = paste(unique(na.omit(strategy_ok)), collapse = ','),
              reciprocity = median(reciprocity),
              callers_detected = median(callers_detected2),
              name_callers = paste(unique(na.omit(name_callers)), collapse = ','),
              GT_callers = paste(unique(na.omit(geno_callers)), collapse = ','),
              min_window_size = median(min_window_size),
              lower = median(lower),
              upper = median(upper),
              PASS = paste(unique(na.omit(PASS)), collapse = ','),
              PASS_num = median(PASS_num)) %>%
    select(-ID)
  
  final_data = final_data %>% arrange(start)
  final_data = final_data %>% filter(PASS %in% c("PASS","NO_PASS"))
 
  colnames(final_data) = paste0(colnames(final_data),"_",ids[i])
  
fwrite(final_data,paste0("/2_merge_callers/DUP/outputs/",ids[i],"/",ids[i],"_DUP_chr_",j),
                       sep = " ",row.names = F,quote = F)
  
 }
        
  


