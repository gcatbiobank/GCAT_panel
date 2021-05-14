library(data.table)
library(dplyr)
library(caret)
library(e1071)

args <- commandArgs(trailingOnly = TRUE)
i= as.numeric(args[2])
              
source("ext/functions.R")
              
strategies = fread("ext/strategies.csv")
              
mod_fit = readRDS("1_LRM/models/LRM_model_INV.rds")
              
ids = fread("ext/all_samplesok",header = F)
              
ids = ids$V1

print(ids[i])
              
dir.create(paste0("/2_merge_callers/INV/outputs/",ids[i]))
              

# read VCF files ########

delly = fread(paste0("/2_merge_callers/INV/data/Delly/",ids[i],"_INV"))
delly$V4 = as.numeric(abs(delly$V4))
colnames(delly) = c("chr","start_delly","end_delly","length_delly","GT_delly")

delly = unify_breakpoints(delly,"delly")
colnames(delly)[2:5] = c("chr","start_delly","end_delly","GT_delly")
delly$ID = NULL
delly$length_delly = delly$end-delly$start

delly$lower_delly_bp1 = delly$start_delly-20
delly$upper_delly_bp1 = delly$start_delly+20
delly$lower_delly_bp2 = delly$end_delly-20
delly$upper_delly_bp2 = delly$end_delly+20
delly = delly %>% as.data.table()
delly$chr[delly$chr=="X"] = 23
delly$chr[delly$chr=="Y"] = 24
delly$chr = as.character(delly$chr)
delly = delly %>% filter(chr %in% 1:23)  %>% as.data.table()
delly = delly %>% filter(length_delly>30)
delly = unique(delly)
delly$GT2_delly = delly$GT_delly

table(delly$chr)
dim(delly)
summary(delly$end_delly-delly$start_delly)


lumpy = fread(paste0("/2_merge_callers/INV/data/Lumpy/",ids[i],"_INV"))
lumpy$V4 = as.numeric(abs(lumpy$V4))
colnames(lumpy) = c("chr","start_lumpy","end_lumpy","length_lumpy","GT_lumpy")

lumpy = unify_breakpoints(lumpy,"lumpy")
colnames(lumpy)[2:5] = c("chr","start_lumpy","end_lumpy","GT_lumpy")
lumpy$ID = NULL
lumpy$length_lumpy = lumpy$end-lumpy$start

lumpy$lower_lumpy_bp1 = lumpy$start_lumpy-50
lumpy$upper_lumpy_bp1 = lumpy$start_lumpy+50
lumpy$lower_lumpy_bp2 = lumpy$end_lumpy-50
lumpy$upper_lumpy_bp2 = lumpy$end_lumpy+50
lumpy$chr[lumpy$chr=="X"] = 23
lumpy$chr[lumpy$chr=="Y"] = 24
lumpy$chr = as.character(lumpy$chr)
lumpy = unique(lumpy)
lumpy = lumpy %>% filter(chr %in% 1:23) %>% as.data.table()
lumpy = lumpy %>% filter(length_lumpy>30)
lumpy$GT2_lumpy= lumpy$GT_lumpy
table(lumpy$chr)
dim(lumpy)
summary(lumpy$end_lumpy-lumpy$start_lumpy)


pindel = fread(paste0("/2_merge_callers/INV/data/Pindel/",ids[i],"_INV"))
colnames(pindel) = c("chr","start_pindel","end_pindel","length_pindel","GT_pindel")

pindel = unify_breakpoints(pindel,"pindel")
colnames(pindel)[2:5] = c("chr","start_pindel","end_pindel","GT_pindel")
pindel$ID = NULL
pindel$length_pindel = pindel$end-pindel$start

pindel$lower_pindel_bp1 = pindel$start_pindel-10
pindel$upper_pindel_bp1 = pindel$start_pindel+10
pindel$lower_pindel_bp2 = pindel$end_pindel-10
pindel$upper_pindel_bp2 = pindel$end_pindel+10
pindel$chr[pindel$chr=="X"] = 23
pindel$chr[pindel$chr=="Y"] = 24
pindel$chr = as.character(pindel$chr)
pindel = unique(pindel)
pindel = pindel %>% filter(chr %in% 1:23) %>% as.data.table()
pindel = pindel %>% filter(length_pindel>30)
pindel$GT2_pindel = pindel$GT_pindel
table(pindel$chr)
dim(pindel)
summary(pindel$end_pindel-pindel$start_pindel)


whamg = fread(paste0("/2_merge_callers/INV/data/Whamg/",ids[i],"_INV"))
colnames(whamg) = c("chr","start_whamg","end_whamg","length_whamg","GT_whamg")

whamg = unify_breakpoints(whamg,"whamg")
colnames(whamg)[2:5] = c("chr","start_whamg","end_whamg","GT_whamg")
whamg$ID = NULL
whamg$length_whamg = whamg$end-whamg$start

whamg$lower_whamg_bp1 = whamg$start_whamg-10
whamg$upper_whamg_bp1 = whamg$start_whamg+10
whamg$lower_whamg_bp2 = whamg$end_whamg-10
whamg$upper_whamg_bp2 = whamg$end_whamg+10
whamg$chr[whamg$chr=="X"] = 23
whamg$chr[whamg$chr=="Y"] = 24
whamg$chr = as.character(whamg$chr)
whamg = unique(whamg)
whamg = whamg %>% filter(chr %in% 1:23) %>% as.data.table()
whamg = whamg %>% filter(length_whamg>30)
whamg$GT2_whamg = whamg$GT_whamg
dim (whamg)
table(whamg$chr)
summary(whamg$end_whamg-whamg$start_whamg)

manta = fread(paste0("/2_merge_callers/INV/data/Manta/",ids[i],"_INV"))
manta$V4 = as.numeric(abs(manta$V4))
colnames(manta) = c("chr","start_manta","end_manta","length_manta","GT_manta")

manta = unify_breakpoints(manta,"manta")
colnames(manta)[2:5] = c("chr","start_manta","end_manta","GT_manta")
manta$ID = NULL
manta$length_manta = manta$end-manta$start

manta$lower_manta_bp1 = manta$start_manta-10
manta$upper_manta_bp1 = manta$start_manta+10
manta$lower_manta_bp2 = manta$end_manta-10
manta$upper_manta_bp2 = manta$end_manta+10
manta$chr[manta$chr=="X"] = 23
manta$chr[manta$chr=="Y"] = 24
manta$chr = as.character(manta$chr)
manta = unique(manta)
manta = manta %>% filter(chr %in% 1:23) %>% as.data.table()
manta = manta %>% filter(length_manta>30)
manta$GT2_manta = manta$GT_manta
dim (manta)
table(manta$chr)


svaba = fread(paste0("/2_merge_callers/INV/data/SVaBA/",ids[i],"_INV"))
svaba$length_aux = svaba$V4-svaba$V2

for(j in 1:nrow(svaba)){
  
  if(svaba$length_aux[j]<0){
    
    aux = svaba$V2[j]
    aux2 = svaba$V4[j]
    
    svaba$V2[j] = aux2
    svaba$V4[j] = aux
    
  }
  
}

svaba$length_aux = svaba$V4-svaba$V2
svaba = svaba[,c(1,2,3,4,5,6)]
colnames(svaba) = c("chr","start_svaba","chr2","end_svaba","GT_svaba","length_svaba")


remove_variants = NULL

for(variant in 1:nrow(svaba)){
  
  if(svaba$chr[variant]==svaba$chr2[variant] & svaba$start_svaba[variant]>svaba$end_svaba[variant]){
    
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
svaba = unique(svaba)
svaba = svaba[,c(1,2,4,6,5)]

svaba = unify_breakpoints(svaba,"svaba")
colnames(svaba)[2:5] = c("chr","start_svaba","end_svaba","GT_svaba")
svaba$ID = NULL
svaba$length_svaba = svaba$end-svaba$start

svaba$chr[svaba$chr=="X"] = 23
svaba$chr[svaba$chr=="Y"] = 24
svaba$chr = as.character(svaba$chr)
svaba = svaba %>% filter(chr %in% 1:23) %>% as.data.table()
svaba = svaba %>% filter(length_svaba>30)
svaba$GT2_svaba = svaba$GT_svaba
dim(svaba)
table(svaba$chr)

svaba$lower_svaba_bp1 = svaba$start_svaba-10
svaba$upper_svaba_bp1 = svaba$start_svaba+10
svaba$lower_svaba_bp2 = svaba$end_svaba-10
svaba$upper_svaba_bp2 = svaba$end_svaba+10

call_windows = data.frame(callers = c("delly","lumpy","pindel","whamg",
                                     "svaba","manta"),
                          window = c(20,50,10,10,10,10))

call_windows$callers = as.character(call_windows$caller)

      

# filter by chromosome ########

j= as.numeric(args[1])

        
        delly_chr = delly %>% filter(chr==j) %>% arrange(start_delly) %>% as.data.table() %>% unique()
        lumpy_chr = lumpy %>% filter(chr==j) %>% arrange(start_lumpy) %>% as.data.table() %>% unique()
        pindel_chr = pindel %>% filter(chr==j) %>% arrange(start_pindel) %>% as.data.table() %>% unique()
        whamg_chr = whamg %>% filter(chr==j) %>% arrange(start_whamg) %>% as.data.table() %>% unique()
        svaba_chr = svaba %>% filter(chr==j) %>% arrange(start_svaba) %>% as.data.table() %>% unique()
        manta_chr = manta %>% filter(chr==j) %>% arrange(start_manta) %>% as.data.table() %>% unique()
        
        n_rows = nrow(delly_chr) + nrow(lumpy_chr) + nrow(pindel_chr) + nrow(whamg_chr) + nrow(svaba_chr) + nrow(manta_chr)
        
                
                
        type = "inversion"
                
        delly_chr$ID = "none"
                
        delly_chr$ID[1] = paste0(type,"_",1)
                
        for(id in 2:nrow(delly_chr)){
                  
          aux = delly_chr[(id-1):id,c(as.vector(outer(c("start_","end_","length_"),"delly",paste0)),"chr"),with=F]
                  
          reprocicity = as.numeric(aux[1,3])/as.numeric(aux[2,3])
                  
          reprocicity = ifelse(reprocicity>1,as.numeric(aux[2,3])/as.numeric(aux[1,3]),reprocicity)
                  
           if(delly_chr$chr[id]==delly_chr$chr[id-1] & reprocicity>0.80 & (abs(delly_chr$start_delly[id]-delly_chr$start_delly[id-1])<=50 | # maxim error que tolerem
                                                                        abs(delly_chr$end_delly[id]-delly_chr$end_delly[id-1])<=50)){
                    
           delly_chr$ID[id] = delly_chr$ID[id-1]
                    
            }
            else{delly_chr$ID[id] = paste0(type,"_",id)}
          
        }
                  
                
        lumpy_chr$ID = "none"
        manta_chr$ID = "none"
        whamg_chr$ID = "none"
        pindel_chr$ID = "none"
        svaba_chr$ID = "none"
        
        f_score_order = c("delly","lumpy","manta","whamg",
                          "pindel","svaba")
        
      

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
        
        my_data = merge_callers_inv(get(paste0(second_call,"_chr")),get(paste0(first_call,"_chr")),callers_ref = first_call,
                                                       caller_to_merge = second_call,type="inversion") 
        
        for(caller_merge in (caller_second+1):length(f_score_order)){
          
          n_row_call = nrow(get(paste0(f_score_order[caller_merge],"_chr")))
          
          while(n_row_call!=0){
            
            callers_ref = tstrsplit(colnames(my_data)[grep("GT_",colnames(my_data))],split="GT_")[[2]]
            
            my_data = merge_callers_inv(get(paste0(f_score_order[caller_merge],"_chr")),my_data,callers_ref = callers_ref,
                                                           caller_to_merge = f_score_order[caller_merge],type="inversion") 
            
            n_row_call = 0
            
          }
          
        }
        
        
        callers_all = tstrsplit(colnames(my_data)[grep("GT_",colnames(my_data))],split="GT_")[[2]]
        
        callers = c("manta","whamg","delly","lumpy","pindel","svaba")
        
        if(sum(callers %in% callers_all)!=6){
          
          callers_missing = callers[-which(callers %in% callers_all)]
          
          for(call_miss in 1:length(callers_missing)){
            
            my_data$start_miss = NA
            my_data$end_miss = NA
            my_data$GT_miss = NA
            my_data$length_miss = NA
            
            colnames(my_data)[(ncol(my_data)-3):ncol(my_data)] = paste0(c("start_","end_","GT_","length_"),callers_missing[call_miss])
            
          }
          
          
        }
        
        
        
        my_data2 = my_data       
        dim(my_data)
        dim(my_data2)
        
        my_data2 = unique(my_data2)
                
                
              
                
        colnames_bbdd = as.vector(outer(c("start_","end_","length_",
                                          "GT_"),callers,paste0))
        
        final_data = my_data2[,c("ID","chr",colnames_bbdd),with=F]        
                
        final_data = final_data %>% group_by(ID) %>% 
          dplyr::summarise(chr=names(sort(table(chr,useNA = "always"),decreasing = T))[1],
                           start = names(sort(table(c(start_pindel,start_manta,start_whamg,start_lumpy,start_delly,start_svaba)),decreasing = T))[1],
                           end = names(sort(table(c(end_pindel,end_manta,end_whamg,end_lumpy,end_delly,end_svaba)),decreasing = T))[1],                   
                           GT_whamg = names(sort(table(GT_whamg,useNA = "always"),decreasing = T))[1],
                           GT_pindel = names(sort(table(GT_pindel,useNA = "always"),decreasing = T))[1],
                           GT_delly = names(sort(table(GT_delly,useNA = "always"),decreasing = T))[1],
                           GT_lumpy = names(sort(table(GT_lumpy,useNA = "always"),decreasing = T))[1],
                           GT_manta = names(sort(table(GT_manta,useNA = "always"),decreasing = T))[1],
                           GT_svaba = names(sort(table(GT_svaba,useNA = "always"),decreasing = T))[1],
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
                           start_svaba = names(sort(table(start_svaba,useNA = "always"),decreasing = T))[1]) %>%
          as.data.table()
        
        
        
        final_data$start = consensus_start_inversions(final_data,callers_all)
        
        final_data$length = consensus_length(final_data,callers_all)
        
        final_data$end = final_data$start + final_data$length
        
        final_data$length = abs(as.numeric(final_data$end)-as.numeric(final_data$start))
        
        summary(final_data$length)
        
        final_data = final_data %>% filter(length>30)     
                
        final_data$reciprocity = as.numeric(final_data$min_length)/as.numeric(final_data$max_length)
        
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
                                                         callers = callers_all)
        
        final_data$callers_detected_ok = n_callers_detected(final_data,
                                                            callers = callers_all)
        
        
        table(final_data$callers_detected)
        
        final_data$callers_detected[final_data$callers_detected=="3"] = "3-4-5-6"
        final_data$callers_detected[final_data$callers_detected=="4"] = "3-4-5-6"
        final_data$callers_detected[final_data$callers_detected=="5"] = "3-4-5-6"
        final_data$callers_detected[final_data$callers_detected=="6"] = "3-4-5-6"
        
        table(final_data$callers_detected)
        table(final_data$callers_detected_ok)        
                
        ## add strategy ####
        
        final_data = final_data %>% as.data.table()
        
        final_data$strategy = 0
        
        final_data$strategy = strategy(final_data,
                                       callers = callers_all,
                                       strategies)
        
        final_data$strategy_ok = strategy(final_data,
                                          callers = callers_all,
                                          strategies)
        
        
        table(final_data$strategy)
        
        final_data$strategy[final_data$strategy=="3"] = "3-4"
        final_data$strategy[final_data$strategy=="4"] = "3-4"
        
        table(final_data$strategy)
        table(final_data$strategy_ok)
        
        
        ## geno callers detected ####
        
        final_data$geno_callers = geno_callers(final_data,
                                            callers = callers_all)
        
        ## prepare model ####
        
        
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
 
        
        # GT ->  Order: Lumpy-Pindel-Whamg-Delly-Manta
        
        data_call = final_data %>%  
          dplyr::select(GT_whamg2,
                        GT_lumpy2,
                        GT_pindel2,
                        GT_delly2,
                        GT_manta2,GT_svaba2)
        
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
        data_call$GT_svaba[data_call$GT_svaba=="9/9"] = NA
        
        table(data_call$GT_whamg)
        table(data_call$GT_lumpy)
        table(data_call$GT_pindel)
        table(data_call$GT_delly)
        table(data_call$GT_manta)
        table(data_call$GT_svaba)
        
       
        #GT Consensus######
        
        caller = c("manta","whamg","delly","lumpy","pindel","svaba")
        

        data_call$GT = NA
        
        data_call = as.matrix(data_call)
        
        
        for(geno in 1:nrow(data_call)){
          
          names_call = caller[which(!is.na(data_call[geno,c(7:12)]))]
          
          if(sum(names_call %in% "lumpy")==1){
            
            data_call[geno,13] = data_call[geno,10]
          }
          
          
          if(sum(names_call %in% "lumpy")==0 & sum(names_call %in% "pindel")==1){
            
            data_call[geno,13] = data_call[geno,11]
          }
          
          if(sum(names_call %in% "lumpy")==0 & sum(names_call %in% "pindel")==0 & sum(names_call %in% "whamg")==1){
            
            data_call[geno,13] = data_call[geno,8]
          }
          
          if(sum(names_call %in% "lumpy")==0 & sum(names_call %in% "pindel")==0 & sum(names_call %in% "whamg")==0 &
             sum(names_call %in% "delly")==1){
            
            data_call[geno,13] = data_call[geno,9]
          }
          
          if(sum(names_call %in% "lumpy")==0 & sum(names_call %in% "pindel")==0 & sum(names_call %in% "whamg")==0 &
             sum(names_call %in% "delly")==0 & sum(names_call %in% "manta")==1){
            
            data_call[geno,13] = data_call[geno,7]
          }
          
          if(sum(names_call %in% "lumpy")==0 & sum(names_call %in% "pindel")==0 & sum(names_call %in% "whamg")==0 &
             sum(names_call %in% "delly")==0 & sum(names_call %in% "manta")==0 & sum(names_call %in% "svaba")==1){
            
            data_call[geno,13] = data_call[geno,12]
          }
        }
        
        
        head(data_call)
        dim(data_call)
        
        data_call = as.data.frame(data_call)
        
        table(data_call$GT,useNA = "always")
        
  
      final_data$GT =data_call$GT
        
      #Predict Final_data
        
        data_predict = final_data %>% select(ID,chr,start,end,GT_manta,GT_whamg,
                                          GT_delly,GT_lumpy,GT_pindel,
                                          GT_svaba,length,
                                          length_stretch,strategy,strategy_ok,
                                          reciprocity,callers_detected,callers_detected_ok,geno_callers,GT_manta2,GT_whamg2,
                                          GT_delly2,GT_lumpy2,GT_pindel2,GT_svaba2,GT) %>%
          as.data.frame(row.names=FALSE)
        
        data_predict$GT_manta = as.factor(data_predict$GT_manta)
        data_predict$GT_whamg = as.factor(data_predict$GT_whamg)
        data_predict$GT_delly = as.factor(data_predict$GT_delly)
        data_predict$GT_lumpy = as.factor(data_predict$GT_lumpy)
        data_predict$GT_pindel = as.factor(data_predict$GT_pindel)
        data_predict$GT_svaba = as.factor(data_predict$GT_svaba)
        data_predict$length_stretch = as.factor(data_predict$length_stretch)
        data_predict$strategy = as.factor(data_predict$strategy)
        data_predict$callers_detected = as.factor(data_predict$callers_detected)
        
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
        
        data_predict$name_callers = name_callers_detected(data_predict,
                                                          callers = callers_all)
        
        ## minimum window size ####
        
        data_predict$min_window_size = min_window_size(data_predict,
                                                       callers = callers_all,
                                                       call_windows)
        
        ## window size for merging between samples ####
        
        data_predict$lower = data_predict$start-data_predict$min_window_size
        data_predict$upper = data_predict$start+data_predict$min_window_size
        
      

  ## group by ID to remove duplicates ####
  
  final_predict = data_predict %>% select(ID,chr,start,GT,length,
                                       strategy_ok,reciprocity,
                                       callers_detected_ok,name_callers,
                                       min_window_size,lower,upper,PASS,PASS_num,geno_callers) %>%
    arrange(start)
  
        final_predict = final_predict %>% group_by(ID) %>% 
    summarise(chr = unique(chr),
              start = median(start),
              GT = paste(unique(na.omit(GT)), collapse = ','),
              length = median(length),
              strategy = paste(unique(na.omit(strategy_ok)), collapse = ','),
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
  
        final_predict = final_predict %>% arrange(start)
        final_predict = final_predict %>% filter(PASS %in% c("PASS","NO_PASS")) ##We delete PASS/NOPASS events in the same variant
  
 
colnames(final_predict) = paste0(colnames(final_predict),"_",ids[i])
  
fwrite(final_data,paste0("/2_merge_callers/INV/outputs/",ids[i],"/",ids[i],"_INV_chr_",j),
                       sep = " ",row.names = F,quote = F)


