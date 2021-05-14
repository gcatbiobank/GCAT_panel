library(data.table)
library(dplyr)
library(caret)
library(e1071)

args <- commandArgs(trailingOnly = TRUE)
i= as.numeric(args[2])
source("/gpfs/projects/bsc05/jordivalls/GCAT_project_all_samples/merge_all_calling_GCAT/Translocaciones/functions.R")

mod_fit = readRDS("/gpfs/projects/bsc05/jordivalls/GCAT_project_all_samples/merge_all_calling_GCAT/Translocaciones/model_translocations.rds")

ids = fread("/gpfs/projects/bsc05/jordivalls/GCAT_project_all_samples/merge_all_calling_GCAT/Deleciones/all_samplesok",header = F)

ids = ids$V1

strategies = fread("/gpfs/projects/bsc05/jordivalls/GCAT_project_all_samples/merge_all_calling_GCAT/Deleciones/strategies.csv")

print(ids[i])

dir.create(paste0("/gpfs/projects/bsc05/jordivalls/GCAT_project_all_samples/merge_all_calling_GCAT/Translocaciones/merge_callers_new/",ids[i]))

j= as.numeric(args[1])


# read vcfs Dani files ########

delly = fread(paste0("/gpfs/projects/bsc05/jordivalls/GCAT_project_all_samples/merge_all_calling_GCAT/Translocations/Delly/",ids[i],"_Delly_Translocations"))
colnames(delly) = c("chr_1_delly","start_1_delly","chr_2_delly","start_2_delly","GT_delly")

delly$chr_1_delly[delly$chr_1_delly=="X"] = "23"
delly$chr_1_delly[delly$chr_1_delly=="Y"] = "24"
delly$chr_2_delly[delly$chr_2_delly=="X"] = "23"
delly$chr_2_delly[delly$chr_2_delly=="Y"] = "24"

delly$chr_1_delly = as.numeric(delly$chr_1_delly)
delly$chr_2_delly = as.numeric(delly$chr_2_delly)
delly$start_1_delly = as.numeric(delly$start_1_delly)
delly$start_2_delly = as.numeric(delly$start_2_delly)

nas = which(apply(delly,1,function(x)sum(is.na(x)))>0) # hay chromosomas desconocidos N3

if(length(nas)>0){
  delly = delly[-nas,]}

for(delly_var in 1:nrow(delly)){

  if(delly$chr_1_delly[delly_var]>delly$chr_2_delly[delly_var]){

    chr1 = delly$chr_1_delly[delly_var]
    chr2 = delly$chr_2_delly[delly_var]
    start1 = delly$start_1_delly[delly_var]
    start2 = delly$start_2_delly[delly_var]

    delly$chr_1_delly[delly_var] = chr2
    delly$chr_2_delly[delly_var] = chr1
    delly$start_1_delly[delly_var] = start2
    delly$start_2_delly[delly_var] = start1
  }

}

delly = unique(delly)

delly = delly %>% arrange(chr_1_delly,chr_2_delly,
                          start_1_delly,start_2_delly) %>% as.data.table()

delly$length_delly = 0

delly = delly %>% filter(chr_1_delly==j)

delly_ok = delly

delly_ok$end_1_delly = delly_ok$start_1_delly
delly_ok$end_2_delly = delly_ok$start_2_delly
delly_ok$lower_delly_start_1 = delly_ok$start_1_delly-300
delly_ok$upper_delly_start_1 = delly_ok$start_1_delly+300
delly_ok$lower_delly_start_2 = delly_ok$start_2_delly-300
delly_ok$upper_delly_start_2 = delly_ok$start_2_delly+300

delly_ok$lower_delly_end_1 = delly_ok$end_1_delly-300
delly_ok$upper_delly_end_1 = delly_ok$end_1_delly+300
delly_ok$lower_delly_end_2 = delly_ok$end_2_delly-300
delly_ok$upper_delly_end_2 = delly_ok$end_2_delly+300

delly_ok = delly_ok %>% as.data.table()

if(nrow(delly)>1){

  type = "translocation"

  delly$ID = "none"

  delly$ID[1] = paste0(type,"_",1)

  for(delly_var in 2:nrow(delly)){

    if(delly$chr_1_delly[delly_var]==delly$chr_1_delly[delly_var-1] &
       delly$chr_2_delly[delly_var]==delly$chr_2_delly[delly_var-1] &
       (abs(delly$start_1_delly[delly_var]-delly$start_1_delly[delly_var-1])<=150 |
       abs(delly$start_2_delly[delly_var]-delly$start_2_delly[delly_var-1])<=150) &
       delly$GT_delly[delly_var]==delly$GT_delly[delly_var-1]){

      delly$ID[delly_var] = paste0(delly$ID[delly_var-1])
    }
    else{delly$ID[delly_var] = paste0(type,"_",delly_var)}

  }

  delly_ok = data.frame(ID = unique(delly$ID))

  delly_ok$chr_1_delly = 0
  delly_ok$start_1_delly = 0
  delly_ok$end_1_delly = 0

  delly_ok$chr_2_delly = 0
  delly_ok$start_2_delly = 0
  delly_ok$end_2_delly = 0

  delly_ok$GT_delly = 0
  delly_ok$length_delly = 0

  for(delly_var in 1:nrow(delly_ok)){

    aux = delly %>% filter(ID==delly_ok$ID[delly_var])

    delly_ok$chr_1_delly[delly_var] = unique(aux$chr_1_delly)
    delly_ok$start_1_delly[delly_var] = min(aux$start_1_delly)
    delly_ok$end_1_delly[delly_var] = max(aux$start_1_delly)

    delly_ok$chr_2_delly[delly_var] = unique(aux$chr_2_delly)
    delly_ok$start_2_delly[delly_var] = min(aux$start_2_delly)
    delly_ok$end_2_delly[delly_var] = max(aux$start_2_delly)

    delly_ok$GT_delly[delly_var] = paste0(unique(aux$GT_delly),collapse = "")
    delly_ok$length_delly[delly_var] = max(c(delly_ok$end_2_delly[delly_var]-delly_ok$start_2_delly[delly_var],
                                     delly_ok$end_1_delly[delly_var]-delly_ok$start_1_delly[delly_var]))

  }

  delly_ok$GT_delly[-which(delly_ok$GT_delly %in% c("0/1","1/1"))] = "./."

  delly_ok$lower_delly_start_1 = delly_ok$start_1_delly-300
  delly_ok$upper_delly_start_1 = delly_ok$start_1_delly+300
  delly_ok$lower_delly_start_2 = delly_ok$start_2_delly-300
  delly_ok$upper_delly_start_2 = delly_ok$start_2_delly+300

  delly_ok$lower_delly_end_1 = delly_ok$end_1_delly-300
  delly_ok$upper_delly_end_1 = delly_ok$end_1_delly+300
  delly_ok$lower_delly_end_2 = delly_ok$end_2_delly-300
  delly_ok$upper_delly_end_2 = delly_ok$end_2_delly+300

  delly_ok = delly_ok %>% as.data.table()

}

lumpy = fread(paste0("/gpfs/projects/bsc05/jordivalls/GCAT_project_all_samples/merge_all_calling_GCAT/Translocations/lumpy/",ids[i],"_lumpy_translocation"))
colnames(lumpy) = c("chr_1_lumpy","start_1_lumpy","chr_2_lumpy","start_2_lumpy","GT_lumpy")

lumpy = lumpy %>% filter(chr_1_lumpy!="hs37d5")
lumpy = lumpy %>% filter(chr_2_lumpy!="hs37d5")

lumpy$chr_1_lumpy[lumpy$chr_1_lumpy=="X"] = "23"
lumpy$chr_1_lumpy[lumpy$chr_1_lumpy=="Y"] = "24"
lumpy$chr_1_lumpy[lumpy$chr_1_lumpy=="M"] = "25"

lumpy$chr_2_lumpy[lumpy$chr_2_lumpy=="X"] = "23"
lumpy$chr_2_lumpy[lumpy$chr_2_lumpy=="Y"] = "24"
lumpy$chr_2_lumpy[lumpy$chr_2_lumpy=="M"] = "25"


lumpy$chr_1_lumpy = as.numeric(lumpy$chr_1_lumpy)
lumpy$chr_2_lumpy = as.numeric(lumpy$chr_2_lumpy)
lumpy$start_1_lumpy = as.numeric(lumpy$start_1_lumpy)
lumpy$start_2_lumpy = as.numeric(lumpy$start_2_lumpy)

lumpy = lumpy %>% filter(!is.na(chr_1_lumpy))
lumpy = lumpy %>% filter(!is.na(chr_2_lumpy))

which_intra = NULL

for(lumpy_var in 1:nrow(lumpy)){
  
  if(lumpy$chr_1_lumpy[lumpy_var]>lumpy$chr_2_lumpy[lumpy_var]){
    
    chr1 = lumpy$chr_1_lumpy[lumpy_var]
    chr2 = lumpy$chr_2_lumpy[lumpy_var]
    start1 = lumpy$start_1_lumpy[lumpy_var]
    start2 = lumpy$start_2_lumpy[lumpy_var]
    
    lumpy$chr_1_lumpy[lumpy_var] = chr2
    lumpy$chr_2_lumpy[lumpy_var] = chr1
    lumpy$start_1_lumpy[lumpy_var] = start2
    lumpy$start_2_lumpy[lumpy_var] = start1
  }
  
  if(lumpy$chr_1_lumpy[lumpy_var]==lumpy$chr_2_lumpy[lumpy_var]){
    
    which_intra = c(which_intra,lumpy_var)
    }
}

lumpy = lumpy[-which_intra,]

lumpy = unique(lumpy)

lumpy = lumpy %>% arrange(chr_1_lumpy,chr_2_lumpy,
                          start_1_lumpy,start_2_lumpy) %>% as.data.table()

lumpy$length_lumpy = 0

lumpy = lumpy %>% filter(chr_1_lumpy==j)

lumpy_ok = lumpy

lumpy_ok$end_1_lumpy = lumpy_ok$start_1_lumpy 
lumpy_ok$end_2_lumpy = lumpy_ok$start_2_lumpy 
lumpy_ok$lower_lumpy_start_1 = lumpy_ok$start_1_lumpy-200
lumpy_ok$upper_lumpy_start_1 = lumpy_ok$start_1_lumpy+200
lumpy_ok$lower_lumpy_start_2 = lumpy_ok$start_2_lumpy-200
lumpy_ok$upper_lumpy_start_2 = lumpy_ok$start_2_lumpy+200

lumpy_ok$lower_lumpy_end_1 = lumpy_ok$end_1_lumpy-200
lumpy_ok$upper_lumpy_end_1 = lumpy_ok$end_1_lumpy+200
lumpy_ok$lower_lumpy_end_2 = lumpy_ok$end_2_lumpy-200
lumpy_ok$upper_lumpy_end_2 = lumpy_ok$end_2_lumpy+200

lumpy_ok = lumpy_ok %>% as.data.table()

if(nrow(lumpy)>1){
  
  type = "translocation"
  
  lumpy$ID = "none"
  
  lumpy$ID[1] = paste0(type,"_",1)
  
  for(lumpy_var in 2:nrow(lumpy)){
    
    if(lumpy$chr_1_lumpy[lumpy_var]==lumpy$chr_1_lumpy[lumpy_var-1] & 
       lumpy$chr_2_lumpy[lumpy_var]==lumpy$chr_2_lumpy[lumpy_var-1] & 
       (abs(lumpy$start_1_lumpy[lumpy_var]-lumpy$start_1_lumpy[lumpy_var-1])<=150 | 
         abs(lumpy$start_2_lumpy[lumpy_var]-lumpy$start_2_lumpy[lumpy_var-1])<=150) &
       lumpy$GT_lumpy[lumpy_var]==lumpy$GT_lumpy[lumpy_var-1]){
      
      lumpy$ID[lumpy_var] = paste0(lumpy$ID[lumpy_var-1])
    }
    else{lumpy$ID[lumpy_var] = paste0(type,"_",lumpy_var)}
    
  }
  
  lumpy_ok = data.frame(ID = unique(lumpy$ID))
  
  lumpy_ok$chr_1_lumpy = 0
  lumpy_ok$start_1_lumpy = 0
  lumpy_ok$end_1_lumpy = 0
  
  lumpy_ok$chr_2_lumpy = 0
  lumpy_ok$start_2_lumpy = 0
  lumpy_ok$end_2_lumpy = 0
  
  lumpy_ok$GT_lumpy = 0
  lumpy_ok$length_lumpy = 0
  
  
  for(lumpy_var in 1:nrow(lumpy_ok)){
    
    aux = lumpy %>% filter(ID==lumpy_ok$ID[lumpy_var])
    
    lumpy_ok$chr_1_lumpy[lumpy_var] = unique(aux$chr_1_lumpy)
    lumpy_ok$start_1_lumpy[lumpy_var] = min(aux$start_1_lumpy)
    lumpy_ok$end_1_lumpy[lumpy_var] = max(aux$start_1_lumpy)
    
    lumpy_ok$chr_2_lumpy[lumpy_var] = unique(aux$chr_2_lumpy)
    lumpy_ok$start_2_lumpy[lumpy_var] = min(aux$start_2_lumpy)
    lumpy_ok$end_2_lumpy[lumpy_var] = max(aux$start_2_lumpy)
    
    lumpy_ok$GT_lumpy[lumpy_var] = paste0(unique(aux$GT_lumpy),collapse = "")
    lumpy_ok$length_lumpy[lumpy_var] = max(c(lumpy_ok$end_2_lumpy[lumpy_var]-lumpy_ok$start_2_lumpy[lumpy_var],
                                     lumpy_ok$end_1_lumpy[lumpy_var]-lumpy_ok$start_1_lumpy[lumpy_var]))
    
  }
  
  lumpy_ok$GT_lumpy[-which(lumpy_ok$GT_lumpy %in% c("0/1","1/1"))] = "./."
  
  
  lumpy_ok$lower_lumpy_start_1 = lumpy_ok$start_1_lumpy-200
  lumpy_ok$upper_lumpy_start_1 = lumpy_ok$start_1_lumpy+200
  lumpy_ok$lower_lumpy_start_2 = lumpy_ok$start_2_lumpy-200
  lumpy_ok$upper_lumpy_start_2 = lumpy_ok$start_2_lumpy+200
  
  
  lumpy_ok$lower_lumpy_end_1 = lumpy_ok$end_1_lumpy-200
  lumpy_ok$upper_lumpy_end_1 = lumpy_ok$end_1_lumpy+200
  lumpy_ok$lower_lumpy_end_2 = lumpy_ok$end_2_lumpy-200
  lumpy_ok$upper_lumpy_end_2 = lumpy_ok$end_2_lumpy+200
  
  lumpy_ok = lumpy_ok %>% as.data.table()

}


pindel = fread(paste0("/gpfs/projects/bsc05/jordivalls/GCAT_project_all_samples/merge_all_calling_GCAT/Translocations/Pindel/",ids[i],"_Pindel_translocations"))
pindel$GT_pindel = "0/1"
colnames(pindel) = c("chr_1_pindel","start_1_pindel","chr_2_pindel","start_2_pindel","GT_pindel")

pindel$chr_1_pindel[pindel$chr_1_pindel=="X"] = "23"
pindel$chr_1_pindel[pindel$chr_1_pindel=="Y"] = "24"
pindel$chr_2_pindel[pindel$chr_2_pindel=="X"] = "23"
pindel$chr_2_pindel[pindel$chr_2_pindel=="Y"] = "24"

pindel$chr_1_pindel = as.numeric(pindel$chr_1_pindel)
pindel$chr_2_pindel = as.numeric(pindel$chr_2_pindel)
pindel$start_1_pindel = as.numeric(pindel$start_1_pindel)
pindel$start_2_pindel = as.numeric(pindel$start_2_pindel)

pindel = pindel %>% filter(chr_1_pindel %in% 1:24)
pindel = pindel %>% filter(chr_2_pindel %in% 1:24)

nas = which(apply(pindel,1,function(x)sum(is.na(x)))>0) # hay chromosomas desconocidos N3

if(length(nas)>0){
  pindel = pindel[-nas,]}


for(pindel_var in 1:nrow(pindel)){
  
  if(pindel$chr_1_pindel[pindel_var]>pindel$chr_2_pindel[pindel_var]){
    
    chr1 = pindel$chr_1_pindel[pindel_var]
    chr2 = pindel$chr_2_pindel[pindel_var]
    start1 = pindel$start_1_pindel[pindel_var]
    start2 = pindel$start_2_pindel[pindel_var]
    
    pindel$chr_1_pindel[pindel_var] = chr2
    pindel$chr_2_pindel[pindel_var] = chr1
    pindel$start_1_pindel[pindel_var] = start2
    pindel$start_2_pindel[pindel_var] = start1
  }
  
}

pindel = unique(pindel)

pindel$chr_1_pindel = as.numeric(pindel$chr_1_pindel)
pindel$chr_2_pindel = as.numeric(pindel$chr_2_pindel)
pindel$start_1_pindel = as.numeric(pindel$start_1_pindel)
pindel$start_2_pindel = as.numeric(pindel$start_2_pindel)

pindel = pindel %>% arrange(chr_1_pindel,chr_2_pindel,
                          start_1_pindel,start_2_pindel) %>% as.data.table()

pindel$length_pindel = 0

pindel = pindel %>% filter(chr_1_pindel==j)

pindel_ok = pindel

pindel_ok$end_1_pindel = pindel_ok$start_1_pindel 
pindel_ok$end_2_pindel = pindel_ok$start_2_pindel 
pindel_ok$lower_pindel_start_1 = pindel_ok$start_1_pindel-10
pindel_ok$upper_pindel_start_1 = pindel_ok$start_1_pindel+10
pindel_ok$lower_pindel_start_2 = pindel_ok$start_2_pindel-10
pindel_ok$upper_pindel_start_2 = pindel_ok$start_2_pindel+10

pindel_ok$lower_pindel_end_1 = pindel_ok$end_1_pindel-10
pindel_ok$upper_pindel_end_1 = pindel_ok$end_1_pindel+10
pindel_ok$lower_pindel_end_2 = pindel_ok$end_2_pindel-10
pindel_ok$upper_pindel_end_2 = pindel_ok$end_2_pindel+10

pindel_ok = pindel_ok %>% as.data.table()

if(nrow(pindel)>1){
  
  type = "translocation"
  
  pindel$ID = "none"
  
  pindel$ID[1] = paste0(type,"_",1)
  
  for(pindel_var in 2:nrow(pindel)){
    
    if(pindel$chr_1_pindel[pindel_var]==pindel$chr_1_pindel[pindel_var-1] & 
       pindel$chr_2_pindel[pindel_var]==pindel$chr_2_pindel[pindel_var-1] & 
       (abs(pindel$start_1_pindel[pindel_var]-pindel$start_1_pindel[pindel_var-1])<=150 | 
         abs(pindel$start_2_pindel[pindel_var]-pindel$start_2_pindel[pindel_var-1])<=150)){
      
      pindel$ID[pindel_var] = paste0(pindel$ID[pindel_var-1])
    }
    else{pindel$ID[pindel_var] = paste0(type,"_",pindel_var)}
    
  }
  
  pindel_ok = data.frame(ID = unique(pindel$ID))
  
  pindel_ok$chr_1_pindel = 0
  pindel_ok$start_1_pindel = 0
  pindel_ok$end_1_pindel = 0
  
  pindel_ok$chr_2_pindel = 0
  pindel_ok$start_2_pindel = 0
  pindel_ok$end_2_pindel = 0
  
  pindel_ok$GT_pindel = "0/1"
  pindel_ok$length_pindel = 0
  
  for(pindel_var in 1:nrow(pindel_ok)){
    
    aux = pindel %>% filter(ID==pindel_ok$ID[pindel_var])
    
    pindel_ok$chr_1_pindel[pindel_var] = unique(aux$chr_1_pindel)
    pindel_ok$start_1_pindel[pindel_var] = min(aux$start_1_pindel)
    pindel_ok$end_1_pindel[pindel_var] = max(aux$start_1_pindel)
    
    pindel_ok$chr_2_pindel[pindel_var] = unique(aux$chr_2_pindel)
    pindel_ok$start_2_pindel[pindel_var] = min(aux$start_2_pindel)
    pindel_ok$end_2_pindel[pindel_var] = max(aux$start_2_pindel)
    
    pindel_ok$GT_pindel[pindel_var] = paste0(unique(aux$GT_pindel),collapse = "")
    pindel_ok$length_pindel[pindel_var] = max(c(pindel_ok$end_2_pindel[pindel_var]-pindel_ok$start_2_pindel[pindel_var],
                                             pindel_ok$end_1_pindel[pindel_var]-pindel_ok$start_1_pindel[pindel_var]))
    
  }
  
  pindel_ok$GT_pindel[-which(pindel_ok$GT_pindel %in% c("0/1","1/1"))] = "./."
  
  pindel_ok$lower_pindel_start_1 = pindel_ok$start_1_pindel-10
  pindel_ok$upper_pindel_start_1 = pindel_ok$start_1_pindel+10
  pindel_ok$lower_pindel_start_2 = pindel_ok$start_2_pindel-10
  pindel_ok$upper_pindel_start_2 = pindel_ok$start_2_pindel+10
  
  pindel_ok$lower_pindel_end_1 = pindel_ok$end_1_pindel-10
  pindel_ok$upper_pindel_end_1 = pindel_ok$end_1_pindel+10
  pindel_ok$lower_pindel_end_2 = pindel_ok$end_2_pindel-10
  pindel_ok$upper_pindel_end_2 = pindel_ok$end_2_pindel+10
  
  pindel_ok = pindel_ok %>% as.data.table()
  
  head(pindel_ok)

}


manta = fread(paste0("/gpfs/projects/bsc05/jordivalls/GCAT_project_all_samples/merge_all_calling_GCAT/Translocations/manta/",ids[i],"_manta_translocaciones"))
colnames(manta) = c("chr_1_manta","start_1_manta","chr_2_manta","start_2_manta","GT_manta")

manta$chr_1_manta[manta$chr_1_manta=="X"] = "23"
manta$chr_1_manta[manta$chr_1_manta=="Y"] = "24"
manta$chr_1_manta[manta$chr_1_manta=="M"] = "25"
manta$chr_2_manta[manta$chr_2_manta=="X"] = "23"
manta$chr_2_manta[manta$chr_2_manta=="Y"] = "24"
manta$chr_2_manta[manta$chr_2_manta=="M"] = "25"

manta$chr_1_manta = as.numeric(manta$chr_1_manta)
manta$chr_2_manta = as.numeric(manta$chr_2_manta)
manta$start_1_manta = as.numeric(manta$start_1_manta)
manta$start_2_manta = as.numeric(manta$start_2_manta)

manta = manta %>% filter(chr_1_manta %in% 1:24)
manta = manta %>% filter(chr_2_manta %in% 1:24)

nas = which(apply(manta,1,function(x)sum(is.na(x)))>0) # hay chromosomas desconocidos N3

if(length(nas)>0){
  manta = manta[-nas,]}

for(manta_var in 1:nrow(manta)){
  
  if(manta$chr_1_manta[manta_var]>manta$chr_2_manta[manta_var]){
    
    chr1 = manta$chr_1_manta[manta_var]
    chr2 = manta$chr_2_manta[manta_var]
    start1 = manta$start_1_manta[manta_var]
    start2 = manta$start_2_manta[manta_var]
    
    manta$chr_1_manta[manta_var] = chr2
    manta$chr_2_manta[manta_var] = chr1
    manta$start_1_manta[manta_var] = start2
    manta$start_2_manta[manta_var] = start1
  }
  
}

manta = unique(manta)

manta = manta %>% arrange(chr_1_manta,chr_2_manta,
                            start_1_manta,start_2_manta) %>% as.data.table()

manta$length_manta = 0 

manta = manta %>% filter(chr_1_manta==j)

manta_ok = manta

manta_ok$end_1_manta = manta_ok$start_1_manta 
manta_ok$end_2_manta = manta_ok$start_2_manta 
manta_ok$lower_manta_start_1 = manta_ok$start_1_manta-200
manta_ok$upper_manta_start_1 = manta_ok$start_1_manta+200
manta_ok$lower_manta_start_2 = manta_ok$start_2_manta-200
manta_ok$upper_manta_start_2 = manta_ok$start_2_manta+200

manta_ok$lower_manta_end_1 = manta_ok$end_1_manta-200
manta_ok$upper_manta_end_1 = manta_ok$end_1_manta+200
manta_ok$lower_manta_end_2 = manta_ok$end_2_manta-200
manta_ok$upper_manta_end_2 = manta_ok$end_2_manta+200

manta_ok = manta_ok %>% as.data.table()

if(nrow(manta)>1){
  
  type = "translocation"
  
  manta$ID = "none"
  
  manta$ID[1] = paste0(type,"_",1)
  
  for(manta_var in 2:nrow(manta)){
    
    if(manta$chr_1_manta[manta_var]==manta$chr_1_manta[manta_var-1] & 
       manta$chr_2_manta[manta_var]==manta$chr_2_manta[manta_var-1] & 
       (abs(manta$start_1_manta[manta_var]-manta$start_1_manta[manta_var-1])<=150 | 
        abs(manta$start_2_manta[manta_var]-manta$start_2_manta[manta_var-1])<=150) & 
       manta$GT_manta[manta_var]==manta$GT_manta[manta_var-1]){
      
      manta$ID[manta_var] = paste0(manta$ID[manta_var-1])
    }
    else{manta$ID[manta_var] = paste0(type,"_",manta_var)}
    
  }
  
  manta_ok = data.frame(ID = unique(manta$ID))
  
  manta_ok$chr_1_manta = 0
  manta_ok$start_1_manta = 0
  manta_ok$end_1_manta = 0
  
  manta_ok$chr_2_manta = 0
  manta_ok$start_2_manta = 0
  manta_ok$end_2_manta = 0
  
  manta_ok$GT_manta = 0
  manta_ok$length_manta = 0
  
  for(manta_var in 1:nrow(manta_ok)){
    
    aux = manta %>% filter(ID==manta_ok$ID[manta_var])
    
    manta_ok$chr_1_manta[manta_var] = unique(aux$chr_1_manta)
    manta_ok$start_1_manta[manta_var] = min(aux$start_1_manta)
    manta_ok$end_1_manta[manta_var] = max(aux$start_1_manta)
    
    manta_ok$chr_2_manta[manta_var] = unique(aux$chr_2_manta)
    manta_ok$start_2_manta[manta_var] = min(aux$start_2_manta)
    manta_ok$end_2_manta[manta_var] = max(aux$start_2_manta)
    
    manta_ok$GT_manta[manta_var] = paste0(unique(aux$GT_manta),collapse = "")
    manta_ok$length_manta[manta_var] = max(c(manta_ok$end_2_manta[manta_var]-manta_ok$start_2_manta[manta_var],
                                                manta_ok$end_1_manta[manta_var]-manta_ok$start_1_manta[manta_var]))
    
  }
  
  manta_ok$GT_manta[-which(manta_ok$GT_manta %in% c("0/1","1/1"))] = "./."
  
  manta_ok$lower_manta_start_1 = manta_ok$start_1_manta-200
  manta_ok$upper_manta_start_1 = manta_ok$start_1_manta+200
  manta_ok$lower_manta_start_2 = manta_ok$start_2_manta-200
  manta_ok$upper_manta_start_2 = manta_ok$start_2_manta+200
  
  manta_ok$lower_manta_end_1 = manta_ok$end_1_manta-200
  manta_ok$upper_manta_end_1 = manta_ok$end_1_manta+200
  manta_ok$lower_manta_end_2 = manta_ok$end_2_manta-200
  manta_ok$upper_manta_end_2 = manta_ok$end_2_manta+200
  
  manta_ok = manta_ok %>% as.data.table()
  
  head(manta_ok)

}

svaba = fread(paste0("/gpfs/projects/bsc05/jordivalls/GCAT_project_all_samples/merge_all_calling_GCAT/Deleciones/SVABA/",ids[i],"_SVABA"))
colnames(svaba) = c("chr_1_svaba","start_1_svaba","chr_2_svaba","start_2_svaba","GT_svaba")

svaba = svaba %>% filter(chr_1_svaba!="hs37d5")
svaba = svaba %>% filter(chr_2_svaba!="hs37d5")

svaba$chr_1_svaba[svaba$chr_1_svaba=="X"] = "23"
svaba$chr_1_svaba[svaba$chr_1_svaba=="Y"] = "24"
svaba$chr_1_svaba[svaba$chr_1_svaba=="M"] = "25"

svaba$chr_2_svaba[svaba$chr_2_svaba=="X"] = "23"
svaba$chr_2_svaba[svaba$chr_2_svaba=="Y"] = "24"
svaba$chr_2_svaba[svaba$chr_2_svaba=="M"] = "25"


svaba$chr_1_svaba = as.numeric(svaba$chr_1_svaba)
svaba$chr_2_svaba = as.numeric(svaba$chr_2_svaba)
svaba$start_1_svaba = as.numeric(svaba$start_1_svaba)
svaba$start_2_svaba = as.numeric(svaba$start_2_svaba)

which_intra = NULL

for(svaba_var in 1:nrow(svaba)){
  
  if(svaba$chr_1_svaba[svaba_var]>svaba$chr_2_svaba[svaba_var]){
    
    chr1 = svaba$chr_1_svaba[svaba_var]
    chr2 = svaba$chr_2_svaba[svaba_var]
    start1 = svaba$start_1_svaba[svaba_var]
    start2 = svaba$start_2_svaba[svaba_var]
    
    svaba$chr_1_svaba[svaba_var] = chr2
    svaba$chr_2_svaba[svaba_var] = chr1
    svaba$start_1_svaba[svaba_var] = start2
    svaba$start_2_svaba[svaba_var] = start1
  }
  
  if(svaba$chr_1_svaba[svaba_var]==svaba$chr_2_svaba[svaba_var]){
    
    which_intra = c(which_intra,svaba_var)
  }
}

svaba = svaba[-which_intra,]

svaba = unique(svaba)

svaba = svaba %>% arrange(chr_1_svaba,chr_2_svaba,
                          start_1_svaba,start_2_svaba) %>% as.data.table()

svaba$length_svaba = 0 

svaba = svaba %>% filter(chr_1_svaba==j)

svaba_ok = svaba

svaba_ok$end_1_svaba = svaba_ok$start_1_svaba 
svaba_ok$end_2_svaba = svaba_ok$start_2_svaba 
svaba_ok$lower_svaba_start_1 = svaba_ok$start_1_svaba-300
svaba_ok$upper_svaba_start_1 = svaba_ok$start_1_svaba+300
svaba_ok$lower_svaba_start_2 = svaba_ok$start_2_svaba-300
svaba_ok$upper_svaba_start_2 = svaba_ok$start_2_svaba+300

svaba_ok$lower_svaba_end_1 = svaba_ok$end_1_svaba-300
svaba_ok$upper_svaba_end_1 = svaba_ok$end_1_svaba+300
svaba_ok$lower_svaba_end_2 = svaba_ok$end_2_svaba-300
svaba_ok$upper_svaba_end_2 = svaba_ok$end_2_svaba+300

svaba_ok = svaba_ok %>% as.data.table()

if(nrow(svaba)>1){

  type = "translocation"
  
  svaba$ID = "none"
  
  svaba$ID[1] = paste0(type,"_",1)
  
  for(svaba_var in 2:nrow(svaba)){
    
    if(svaba$chr_1_svaba[svaba_var]==svaba$chr_1_svaba[svaba_var-1] & 
       svaba$chr_2_svaba[svaba_var]==svaba$chr_2_svaba[svaba_var-1] & 
       (abs(svaba$start_1_svaba[svaba_var]-svaba$start_1_svaba[svaba_var-1])<=150 | 
        abs(svaba$start_2_svaba[svaba_var]-svaba$start_2_svaba[svaba_var-1])<=150) &
       svaba$GT_svaba[svaba_var]==svaba$GT_svaba[svaba_var-1]){
      
      svaba$ID[svaba_var] = paste0(svaba$ID[svaba_var-1])
    }
    else{svaba$ID[svaba_var] = paste0(type,"_",svaba_var)}
    
  }
  
  svaba_ok = data.frame(ID = unique(svaba$ID))
  
  svaba_ok$chr_1_svaba = 0
  svaba_ok$start_1_svaba = 0
  svaba_ok$end_1_svaba = 0
  
  svaba_ok$chr_2_svaba = 0
  svaba_ok$start_2_svaba = 0
  svaba_ok$end_2_svaba = 0
  
  svaba_ok$GT_svaba = 0
  svaba_ok$length_svaba = 0
  
  
  for(svaba_var in 1:nrow(svaba_ok)){
    
    aux = svaba %>% filter(ID==svaba_ok$ID[svaba_var])
    
    svaba_ok$chr_1_svaba[svaba_var] = unique(aux$chr_1_svaba)
    svaba_ok$start_1_svaba[svaba_var] = min(aux$start_1_svaba)
    svaba_ok$end_1_svaba[svaba_var] = max(aux$start_1_svaba)
    
    svaba_ok$chr_2_svaba[svaba_var] = unique(aux$chr_2_svaba)
    svaba_ok$start_2_svaba[svaba_var] = min(aux$start_2_svaba)
    svaba_ok$end_2_svaba[svaba_var] = max(aux$start_2_svaba)
    
    svaba_ok$GT_svaba[svaba_var] = paste0(unique(aux$GT_svaba),collapse = "")
    svaba_ok$length_svaba[svaba_var] = max(c(svaba_ok$end_2_svaba[svaba_var]-svaba_ok$start_2_svaba[svaba_var],
                                             svaba_ok$end_1_svaba[svaba_var]-svaba_ok$start_1_svaba[svaba_var]))
    
  }
  
  svaba_ok$GT_svaba[-which(svaba_ok$GT_svaba %in% c("0/1","1/1"))] = "./."
  
  
  svaba_ok$lower_svaba_start_1 = svaba_ok$start_1_svaba-300
  svaba_ok$upper_svaba_start_1 = svaba_ok$start_1_svaba+300
  svaba_ok$lower_svaba_start_2 = svaba_ok$start_2_svaba-300
  svaba_ok$upper_svaba_start_2 = svaba_ok$start_2_svaba+300
  
  
  svaba_ok$lower_svaba_end_1 = svaba_ok$end_1_svaba-300
  svaba_ok$upper_svaba_end_1 = svaba_ok$end_1_svaba+300
  svaba_ok$lower_svaba_end_2 = svaba_ok$end_2_svaba-300
  svaba_ok$upper_svaba_end_2 = svaba_ok$end_2_svaba+300
  
  svaba_ok = svaba_ok %>% as.data.table()
  
  head(svaba_ok)

}  
  
call_windows = data.frame(caller = c("delly","lumpy","pindel",
                                     "svaba","manta"),
                          window = c(300,200,10,300,200))

call_windows$caller = as.character(call_windows$caller)



n_rows = nrow(delly_ok) + nrow(lumpy_ok) + nrow(pindel_ok) + nrow(svaba_ok) + nrow(manta_ok)

num_callers = c(nrow(delly_ok)>0, nrow(lumpy_ok)>0,  nrow(pindel_ok)>0, nrow(svaba_ok)>0,nrow(manta_ok)>0)

if(n_rows!=0 & sum(num_callers)>1){
  
  ## prepare BBDD ordered by F-score ####
  
  f_score_order = c("manta","lumpy","delly","svaba","pindel")
  
  if(nrow(manta_ok)>0){manta_ok$ID = "none"}
  if(nrow(lumpy_ok)>0){lumpy_ok$ID = "none"}
  if(nrow(delly_ok)>0){delly_ok$ID = "none"}
  if(nrow(svaba_ok)>0){svaba_ok$ID = "none"}
  if(nrow(pindel_ok)>0){pindel_ok$ID = "none"}
  
  
  caller_first = 1
  
  while(nrow(get(paste0(f_score_order[caller_first],"_ok")))==0){
    
    caller_first = caller_first + 1
  }
  
  caller_second = caller_first + 1
  
  while(nrow(get(paste0(f_score_order[caller_second],"_ok")))==0){
    
    caller_second = caller_second + 1
  }
  
  
  first_call = f_score_order[caller_first]
  second_call = f_score_order[caller_second]
  
  if(first_call=="manta"){manta_ok$ID = paste0("translocation_",1:nrow(manta_ok))}
  if(first_call=="lumpy"){lumpy_ok$ID = paste0("translocation_",1:nrow(lumpy_ok))}
  if(first_call=="delly"){delly_ok$ID = paste0("translocation_",1:nrow(delly_ok))}
  if(first_call=="svaba"){svaba_ok$ID = paste0("translocation_",1:nrow(svaba_ok))}
  if(first_call=="pindel"){pindel_ok$ID = paste0("translocation_",1:nrow(pindel_ok))}
  
  my_data = merge_callers_bp_trans_bis(get(paste0(second_call,"_ok")),get(paste0(first_call,"_ok")),callers_ref = first_call,
                                                 caller_to_merge = second_call,merge = T,golden = F,
                                       type="translocation") 
  
  for(caller_merge in (caller_second+1):length(f_score_order)){
    
    n_row_call = nrow(get(paste0(f_score_order[caller_merge],"_ok")))
    
    while(n_row_call!=0){
      
      callers_ref = tstrsplit(colnames(my_data)[grep("GT_",colnames(my_data))],split="GT_")[[2]]
      
      my_data = merge_callers_bp_trans_bis(get(paste0(f_score_order[caller_merge],"_ok")),my_data,callers_ref = callers_ref,
                                                     caller_to_merge = f_score_order[caller_merge],merge = T,golden = F,
                                           type="translocation")  
      
      n_row_call = 0
      
    }
    
  }
  
  
  callers_all = tstrsplit(colnames(my_data)[grep("GT_",colnames(my_data))],split="GT_")[[2]]
  
  
  ### consensus chr and start  #######
  
  if(nrow(manta_ok)==0){my_data$start_1_manta = NA}
  if(nrow(lumpy_ok)==0){my_data$start_1_lumpy = NA}
  if(nrow(delly_ok)==0){my_data$start_1_delly = NA}
  if(nrow(svaba_ok)==0){my_data$start_1_svaba = NA}
  if(nrow(pindel_ok)==0){my_data$start_1_pindel = NA}
  
  if(nrow(manta_ok)==0){my_data$start_2_manta = NA}
  if(nrow(lumpy_ok)==0){my_data$start_2_lumpy = NA}
  if(nrow(delly_ok)==0){my_data$start_2_delly = NA}
  if(nrow(svaba_ok)==0){my_data$start_2_svaba = NA}
  if(nrow(pindel_ok)==0){my_data$start_2_pindel = NA}
  
  if(nrow(manta_ok)==0){my_data$end_1_manta = NA}
  if(nrow(lumpy_ok)==0){my_data$end_1_lumpy = NA}
  if(nrow(delly_ok)==0){my_data$end_1_delly = NA}
  if(nrow(svaba_ok)==0){my_data$end_1_svaba = NA}
  if(nrow(pindel_ok)==0){my_data$end_1_pindel = NA}
  
  if(nrow(manta_ok)==0){my_data$end_2_manta = NA}
  if(nrow(lumpy_ok)==0){my_data$end_2_lumpy = NA}
  if(nrow(delly_ok)==0){my_data$end_2_delly = NA}
  if(nrow(svaba_ok)==0){my_data$end_2_svaba = NA}
  if(nrow(pindel_ok)==0){my_data$end_2_pindel = NA}
  
  my_data$start_1 = consensus_start_1_translocations(my_data,
                                                     callers = callers_all)
  
  my_data$start_2 = consensus_start_2_translocations(my_data,
                                                     callers = callers_all)
  
  my_data$end_1 = consensus_end_1_translocations(my_data,
                                                     callers = callers_all)
  
  my_data$end_2 = consensus_end_2_translocations(my_data,
                                                     callers = callers_all)
  
  my_data$chr_1 = j
  
  aux_chr = my_data[,paste0("chr_2_",callers_all),with=F]
  
  my_data$chr_2 = apply(aux_chr,1,function(x){unique(x[!is.na(x)])})
  
  my_data$length_1 = my_data$end_1-my_data$start_1
  my_data$length_2 = my_data$end_2-my_data$start_2
  
  if(nrow(manta_ok)==0){my_data$GT_manta = NA}
  if(nrow(lumpy_ok)==0){my_data$GT_lumpy = NA}
  if(nrow(delly_ok)==0){my_data$GT_delly = NA}
  if(nrow(svaba_ok)==0){my_data$GT_svaba = NA}
  if(nrow(pindel_ok)==0){my_data$GT_pindel = NA}
  
  ## neteja IDs
  
  final_data = my_data %>% group_by(ID) %>% 
    dplyr::summarise(chr_1 = unique(chr_1,na.rm = T),
                     start_1 = median(c(start_1),na.rm = T),
                     end_1 = median(c(end_1),na.rm = T),
                     chr_2 = unique(chr_2,na.rm = T),
                     start_2 = median(c(start_2),na.rm = T),
                     end_2 = median(c(end_2),na.rm = T),
                     GT_pindel = names(sort(table(GT_pindel,useNA = "always"),decreasing = T))[1],
                     GT_delly = names(sort(table(GT_delly,useNA = "always"),decreasing = T))[1],
                     GT_lumpy = names(sort(table(GT_lumpy,useNA = "always"),decreasing = T))[1],
                     GT_manta = names(sort(table(GT_manta,useNA = "always"),decreasing = T))[1],
                     GT_svaba = names(sort(table(GT_svaba,useNA = "always"),decreasing = T))[1],
                     length = max(c(length_1,
                                    length_2))) %>% as.data.frame()
  

  ## consensus Genotype 
  
  final_data$GT = NA
  
  data_call = final_data %>%
    dplyr::select(GT_lumpy,
                  GT_delly,
                  GT_manta,
                  GT)
  
  caller = c("lumpy","delly","manta")
  
  data_call = as.matrix(data_call)
  
  
  for(gt in 1:nrow(data_call)){
    
    names_call = caller[which(!is.na(data_call[gt,c(1:4)]))]
    
    if(sum(names_call %in% "manta")==1){
      
      data_call[gt,4] = data_call[gt,3]
    }
    
    
    if(sum(names_call %in% "manta")==0 & sum(names_call %in% "delly")==1){
      
      data_call[gt,4] = data_call[gt,2]
    }
    
    
    if(sum(names_call %in% "manta")==0 & sum(names_call %in% "delly")==0){
      
      data_call[gt,4] = data_call[gt,1]
    }
    
  }
  
  head(data_call)
  dim(data_call)
  
  data_call = as.data.frame(data_call)
  
  table(data_call$GT,useNA = "always")
  
  final_data$GT = as.character(data_call[,4])
  
  final_data$GT[-which(final_data$GT %in% c("0/1","1/1"))] = "./."
  
  
  final_data = final_data %>% arrange(chr_1,chr_2,start_1)
  
  final_data_ok = final_data
  
  # change GT
  
  
  final_data_ok$GT2_manta = final_data_ok$GT_manta
  final_data_ok$GT2_delly = final_data_ok$GT_delly
  final_data_ok$GT2_lumpy = final_data_ok$GT_lumpy
  final_data_ok$GT2_pindel = final_data_ok$GT_pindel
  final_data_ok$GT2_svaba = final_data_ok$GT_svaba

  final_data_ok$GT_manta[final_data_ok$GT_manta=="./."] = "0/1-1/1"
  final_data_ok$GT_manta[final_data_ok$GT_manta=="0/1"] = "0/1-1/1"
  final_data_ok$GT_manta[final_data_ok$GT_manta=="1/1"] = "0/1-1/1"

  final_data_ok$GT_delly[final_data_ok$GT_delly=="./."] = "0/1-1/1"
  final_data_ok$GT_delly[final_data_ok$GT_delly=="0/1"] = "0/1-1/1"
  final_data_ok$GT_delly[final_data_ok$GT_delly=="1/1"] = "0/1-1/1"
  
  final_data_ok$GT_lumpy[final_data_ok$GT_lumpy=="./."] = "0/1-1/1"
  final_data_ok$GT_lumpy[final_data_ok$GT_lumpy=="0/1"] = "0/1-1/1"
  final_data_ok$GT_lumpy[final_data_ok$GT_lumpy=="1/1"] = "0/1-1/1"
  
  final_data_ok$GT_pindel[final_data_ok$GT_pindel=="./."] = "0/1-1/1"
  final_data_ok$GT_pindel[final_data_ok$GT_pindel=="0/1"] = "0/1-1/1"
  final_data_ok$GT_pindel[final_data_ok$GT_pindel=="1/1"] = "0/1-1/1"
  
  final_data_ok$GT_svaba[final_data_ok$GT_svaba=="./."] = "0/1-1/1"
  final_data_ok$GT_svaba[final_data_ok$GT_svaba=="0/1"] = "0/1-1/1"
  final_data_ok$GT_svaba[final_data_ok$GT_svaba=="1/1"] = "0/1-1/1"
  
  
  final_data_ok = final_data_ok %>% as.data.table()
  
  ## add number of callers detected #####
  
  final_data_ok$callers_detected = n_callers_detected(final_data_ok,
                                                callers = callers_all)
  
  final_data_ok$callers_detected2 = n_callers_detected(final_data_ok,
                                                 callers = callers_all)
  
  final_data_ok$callers_detected[final_data_ok$callers_detected=="4"] = "4-5"
  final_data_ok$callers_detected[final_data_ok$callers_detected=="5"] = "4-5"

  
  ## add reciprocity ###
  
  final_data_ok$reciprocity = 1
  
  
  ## add strategy ####
  
  final_data_ok$strategy = strategy(final_data_ok,
                                 callers = callers_all,
                                 strategies)
  

  
  ## predict using LR model #####
  
  final_data_ok$GT_manta[is.na(final_data_ok$GT_manta)] = "9/9"
  final_data_ok$GT_delly[is.na(final_data_ok$GT_delly)] = "9/9"
  final_data_ok$GT_lumpy[is.na(final_data_ok$GT_lumpy)] = "9/9"
  final_data_ok$GT_pindel[is.na(final_data_ok$GT_pindel)] = "9/9"
  final_data_ok$GT_svaba[is.na(final_data_ok$GT_svaba)] = "9/9"

  
  data_predict = final_data_ok %>% select(ID,chr_1,start_1,end_1,chr_2,start_2,end_2,
                                    GT_manta,GT_delly,GT_lumpy,GT_pindel,
                                    GT_svaba,GT,length,
                                    strategy,
                                    reciprocity,callers_detected,callers_detected2,GT2_manta,
                                    GT2_delly,GT2_lumpy,GT2_pindel,
                                    GT2_svaba) %>%
    as.data.frame(row.names=FALSE)
  
  data_predict$GT_manta = as.factor(data_predict$GT_manta)
  data_predict$GT_delly = as.factor(data_predict$GT_delly)
  data_predict$GT_lumpy = as.factor(data_predict$GT_lumpy)
  data_predict$GT_pindel = as.factor(data_predict$GT_pindel)
  data_predict$GT_svaba = as.factor(data_predict$GT_svaba)
  data_predict$strategy = as.factor(data_predict$strategy)
  data_predict$callers_detected = as.factor(data_predict$callers_detected)
  
  # predict YES/NO threshold 0.5
  
  data_predict = as.data.frame(data_predict)
  
  my_pred = caret::predict.train(mod_fit,data_predict) 
  
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
  
  ## write data ####
  
  data_predict$callers_detected = data_predict$callers_detected2
  
  final_data_out = data_predict %>% dplyr::select(chr_1,start_1,end_1,chr_2,start_2,end_2,
                                                  GT,length,strategy,reciprocity,callers_detected,name_callers,
                                                  min_window_size,lower,upper,PASS,PASS_num) %>%
    arrange(chr_1,chr_2,start_1)
  
  colnames(final_data_out) = paste0(colnames(final_data_out),"_",ids[i])
  
  fwrite(final_data_out,
         paste0("/gpfs/projects/bsc05/jordivalls/GCAT_project_all_samples/merge_all_calling_GCAT/Translocaciones/merge_callers_new/",ids[i],"/",ids[i],"_Trans_chr_",j),
         sep = " ",row.names = F,quote = F)
  
  
        




