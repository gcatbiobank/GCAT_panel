library(data.table)
library(dplyr)

source("ext/functions.R")

setwd("1_LRM/mid_DEL/")

# read golden (order by chr1 and chr2) ######

golden_trans = fread("data/Insilico_TRA",fill = T)

golden_trans$V1 = gsub("chr","",golden_trans$V1)
golden_trans$V3 = gsub("chr","",golden_trans$V3)

golden_trans = golden_trans[,c(1,2,8,3,4,6)]

golden_trans_geno = NULL

for(i in names(table(golden_trans$V1))){
  
  aux = golden_trans[golden_trans$V1==i,]
  
  which_duplicated = aux$V2[which(duplicated(aux$V2))]
  
  golden_dup = unique(aux[aux$V2 %in% which_duplicated,])
  golden_no_dup = aux[!aux$V2 %in% which_duplicated,]
  
  golden_dup$genotype = "1/1"
  golden_no_dup$genotype = "0/1"
  
  aux_geno = rbind(golden_dup,golden_no_dup)
  
  golden_trans_geno = rbind(golden_trans_geno,aux_geno)
}


golden = golden_trans_geno

colnames(golden) = c("chr_1_golden","start_1_golden","end_1_golden",
                     "chr_2_golden","start_2_golden","end_2_golden",
                     "GT_golden")

golden$chr_1_golden[golden$chr_1_golden=="X"] = "23"
golden$chr_1_golden[golden$chr_1_golden=="Y"] = "24"
golden$chr_2_golden[golden$chr_2_golden=="X"] = "23"
golden$chr_2_golden[golden$chr_2_golden=="Y"] = "24"

golden$chr_1_golden = as.numeric(golden$chr_1_golden)
golden$chr_2_golden = as.numeric(golden$chr_2_golden)

remove_same_chr = NULL

for(i in 1:nrow(golden)){
  
  if(golden$chr_1_golden[i]>golden$chr_2_golden[i]){
    
    chr1 = golden$chr_1_golden[i]
    chr2 = golden$chr_2_golden[i]
    start1 = golden$start_1_golden[i]
    start2 = golden$start_2_golden[i]
    end1 = golden$end_1_golden[i]
    end2 = golden$end_2_golden[i]
    
    golden$chr_1_golden[i] = chr2
    golden$chr_2_golden[i] = chr1
    golden$start_1_golden[i] = start2
    golden$start_2_golden[i] = start1
    golden$end_1_golden[i] = end2
    golden$end_2_golden[i] = end1
    
  }
  
  if(golden$chr_1_golden[i]==golden$chr_2_golden[i]){
    
    remove_same_chr = c(remove_same_chr,i)
    
  }
  
}

# remove intra-chromosomal variants 

golden[remove_same_chr,] %>% write.table("data/Insilico_TRA_intrachromosomal",row.names=F,
                                         quote=F)

golden = golden[-remove_same_chr,]

golden$length_golden = 0

golden$start_1_golden = as.numeric(golden$start_1_golden)
golden$start_2_golden = as.numeric(golden$start_2_golden)
golden$end_1_golden = as.numeric(golden$end_1_golden)
golden$end_2_golden = as.numeric(golden$end_2_golden)

golden = golden %>% arrange(chr_1_golden,chr_2_golden,
                            start_1_golden,start_2_golden) %>% as.data.table()

for(i in 1:nrow(golden)){
  
  if(!is.na(golden$end_1_golden[i])){
  
  if(golden$end_1_golden[i]<golden$start_1_golden[i]){
    
    aux  = golden$start_1_golden[i]
    aux2 = golden$end_1_golden[i]
    
    golden$start_1_golden[i] = aux2
    golden$end_1_golden[i] = aux
  }
  
  if(golden$end_2_golden[i]<golden$start_2_golden[i]){
    
    aux  = golden$start_2_golden[i]
    aux2 = golden$end_2_golden[i]
    
    golden$start_2_golden[i] = aux2
    golden$end_2_golden[i] = aux
  }
  
  
  golden$length_golden[i] = max(c(golden$end_2_golden[i]-golden$start_2_golden[i],
                                   golden$end_1_golden[i]-golden$start_1_golden[i]))
  
  }
  
}

golden = golden %>% filter(length_golden>30 | length_golden==0) %>% as.data.table()

golden$end_1_golden[252] = golden$start_1_golden[252] + 1
golden$end_2_golden[252] = golden$start_2_golden[252] + 1
golden$length_golden[252] = 1


golden_pseudo = fread("data/Insilico_TRA_pseudogenes")

golden_pseudo$V4 = nchar(golden_pseudo$V4)

head(golden_pseudo)

golden_pseudo$V1 = gsub("chr","",golden_pseudo$V1)

colnames(golden_pseudo)[c(1,3,4)] = c("chr_1_golden","start_1_golden","length_golden")

golden_pseudo$end_1_golden = golden_pseudo$start_1_golden + 1
golden_pseudo$GT_golden = "0/1"
golden_pseudo$chr_2_golden = 0
golden_pseudo$start_2_golden = 0
golden_pseudo$end_2_golden = 0

golden = rbind(golden,
               golden_pseudo[,c(1,3,6,8,9,10,7,4)])



golden$lower_golden_start_1 = golden$start_1_golden-0
golden$upper_golden_start_1 = golden$start_1_golden+0
golden$lower_golden_start_2 = golden$start_2_golden-0
golden$upper_golden_start_2 = golden$start_2_golden+0

golden$lower_golden_end_1 = golden$end_1_golden-0
golden$upper_golden_end_1 = golden$end_1_golden+0
golden$lower_golden_end_2 = golden$end_2_golden-0
golden$upper_golden_end_2 = golden$end_2_golden+0

tail(golden)

golden = golden %>% as.data.table()

head(golden)


## read vcf ######

window_between_callers = 150

### delly #####

delly = fread("data/Delly_TRA")
colnames(delly) = c("chr_1_delly","start_1_delly","chr_2_delly","start_2_delly","GT_delly")

delly$chr_1_delly[delly$chr_1_delly=="X"] = "23"
delly$chr_1_delly[delly$chr_1_delly=="Y"] = "24"
delly$chr_2_delly[delly$chr_2_delly=="X"] = "23"
delly$chr_2_delly[delly$chr_2_delly=="Y"] = "24"

delly$chr_1_delly = as.numeric(delly$chr_1_delly)
delly$chr_2_delly = as.numeric(delly$chr_2_delly)

for(i in 1:nrow(delly)){
  
  if(delly$chr_1_delly[i]>delly$chr_2_delly[i]){
    
    chr1 = delly$chr_1_delly[i]
    chr2 = delly$chr_2_delly[i]
    start1 = delly$start_1_delly[i]
    start2 = delly$start_2_delly[i]
    
    delly$chr_1_delly[i] = chr2
    delly$chr_2_delly[i] = chr1
    delly$start_1_delly[i] = start2
    delly$start_2_delly[i] = start1
  }
  
}

delly = unique(delly)

delly = delly %>% arrange(chr_1_delly,chr_2_delly,
                          start_1_delly,start_2_delly) %>% as.data.table()

type = "translocation"

delly$ID = "none"

delly$ID[1] = paste0(type,"_",1)

for(i in 2:nrow(delly)){
  
  if(delly$chr_1_delly[i]==delly$chr_1_delly[i-1] & delly$chr_2_delly[i]==delly$chr_2_delly[i-1] & 
      ((abs(delly$start_1_delly[i]-delly$start_1_delly[i-1])<=window_between_callers | 
        abs(delly$start_2_delly[i]-delly$start_2_delly[i-1])<=window_between_callers)) &
        delly$GT_delly[i]==delly$GT_delly[i-1]){
    
    delly$ID[i] = paste0(delly$ID[i-1])
  }
  else{delly$ID[i] = paste0(type,"_",i)}
  
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

i=1

for(i in 1:nrow(delly_ok)){
  
  aux = delly %>% filter(ID==delly_ok$ID[i])
  
  delly_ok$chr_1_delly[i] = unique(aux$chr_1_delly)
  delly_ok$start_1_delly[i] = min(aux$start_1_delly)
  delly_ok$end_1_delly[i] = max(aux$start_1_delly)
  
  delly_ok$chr_2_delly[i] = unique(aux$chr_2_delly)
  print(paste0(i," ",unique(aux$chr_2_delly)))
  delly_ok$start_2_delly[i] = min(aux$start_2_delly)
  delly_ok$end_2_delly[i] = max(aux$start_2_delly)
  
  delly_ok$GT_delly[i] = paste0(unique(aux$GT_delly),collapse = "")
  delly_ok$length_delly[i] = max(c(delly_ok$end_2_delly[i]-delly_ok$start_2_delly[i],
                                   delly_ok$end_1_delly[i]-delly_ok$start_1_delly[i]))
  
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


### lumpy #####

lumpy = fread("data/Lumpy_TRA")
colnames(lumpy) = c("chr_1_lumpy","start_1_lumpy","chr_2_lumpy","start_2_lumpy","GT_lumpy")

lumpy = lumpy %>% filter(chr_1_lumpy!="hs37d5")
lumpy = lumpy %>% filter(chr_2_lumpy!="hs37d5")

lumpy$chr_1_lumpy[lumpy$chr_1_lumpy=="X"] = "23"
lumpy$chr_1_lumpy[lumpy$chr_1_lumpy=="Y"] = "24"
lumpy$chr_2_lumpy[lumpy$chr_2_lumpy=="X"] = "23"
lumpy$chr_2_lumpy[lumpy$chr_2_lumpy=="Y"] = "24"

lumpy$chr_1_lumpy = as.numeric(lumpy$chr_1_lumpy)
lumpy$chr_2_lumpy = as.numeric(lumpy$chr_2_lumpy)

for(i in 1:nrow(lumpy)){
  
  if(lumpy$chr_1_lumpy[i]>lumpy$chr_2_lumpy[i]){
    
    chr1 = lumpy$chr_1_lumpy[i]
    chr2 = lumpy$chr_2_lumpy[i]
    start1 = lumpy$start_1_lumpy[i]
    start2 = lumpy$start_2_lumpy[i]
    
    lumpy$chr_1_lumpy[i] = chr2
    lumpy$chr_2_lumpy[i] = chr1
    lumpy$start_1_lumpy[i] = start2
    lumpy$start_2_lumpy[i] = start1
  }
  
}

lumpy = unique(lumpy)


lumpy = lumpy %>% arrange(chr_1_lumpy,chr_2_lumpy,
                          start_1_lumpy,start_2_lumpy) %>% as.data.table()


type = "translocation"

lumpy$ID = "none"

lumpy$ID[1] = paste0(type,"_",1)

for(i in 2:nrow(lumpy)){
  
  if(lumpy$chr_1_lumpy[i]==lumpy$chr_1_lumpy[i-1] & lumpy$chr_2_lumpy[i]==lumpy$chr_2_lumpy[i-1] & 
     ((abs(lumpy$start_1_lumpy[i]-lumpy$start_1_lumpy[i-1])<=window_between_callers | 
       abs(lumpy$start_2_lumpy[i]-lumpy$start_2_lumpy[i-1])<=window_between_callers)) &
           lumpy$GT_lumpy[i]==lumpy$GT_lumpy[i-1]){
    
    lumpy$ID[i] = paste0(lumpy$ID[i-1])
  }
  else{lumpy$ID[i] = paste0(type,"_",i)}
  
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

i=1

for(i in 1:nrow(lumpy_ok)){
  
  aux = lumpy %>% filter(ID==lumpy_ok$ID[i])
  
  lumpy_ok$chr_1_lumpy[i] = unique(aux$chr_1_lumpy)
  lumpy_ok$start_1_lumpy[i] = min(aux$start_1_lumpy)
  lumpy_ok$end_1_lumpy[i] = max(aux$start_1_lumpy)
  
  lumpy_ok$chr_2_lumpy[i] = unique(aux$chr_2_lumpy)
  lumpy_ok$start_2_lumpy[i] = min(aux$start_2_lumpy)
  lumpy_ok$end_2_lumpy[i] = max(aux$start_2_lumpy)
  
  lumpy_ok$GT_lumpy[i] = paste0(unique(aux$GT_lumpy),collapse = "")
  lumpy_ok$length_lumpy[i] = max(c(lumpy_ok$end_2_lumpy[i]-lumpy_ok$start_2_lumpy[i],
                                     lumpy_ok$end_1_lumpy[i]-lumpy_ok$start_1_lumpy[i]))
  
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


### pindel #####

pindel = fread("data/Pindel_TRA")
colnames(pindel) = c("chr_1_pindel","start_1_pindel","chr_2_pindel","start_2_pindel")

pindel$chr_1_pindel[pindel$chr_1_pindel=="X"] = "23"
pindel$chr_1_pindel[pindel$chr_1_pindel=="Y"] = "24"
pindel$chr_2_pindel[pindel$chr_2_pindel=="X"] = "23"
pindel$chr_2_pindel[pindel$chr_2_pindel=="Y"] = "24"

pindel$chr_1_pindel = as.numeric(pindel$chr_1_pindel)
pindel$chr_2_pindel = as.numeric(pindel$chr_2_pindel)

for(i in 1:nrow(pindel)){
  
  if(pindel$chr_1_pindel[i]>pindel$chr_2_pindel[i]){
    
    chr1 = pindel$chr_1_pindel[i]
    chr2 = pindel$chr_2_pindel[i]
    start1 = pindel$start_1_pindel[i]
    start2 = pindel$start_2_pindel[i]
    
    pindel$chr_1_pindel[i] = chr2
    pindel$chr_2_pindel[i] = chr1
    pindel$start_1_pindel[i] = start2
    pindel$start_2_pindel[i] = start1
  }
  
}

pindel = unique(pindel)

pindel = pindel %>% arrange(chr_1_pindel,chr_2_pindel,
                            start_1_pindel,start_2_pindel) %>% as.data.table()

type = "translocation"

pindel$ID = "none"

pindel$ID[1] = paste0(type,"_",1)

for(i in 2:nrow(pindel)){
  
  if(pindel$chr_1_pindel[i]==pindel$chr_1_pindel[i-1] & pindel$chr_2_pindel[i]==pindel$chr_2_pindel[i-1] & 
     ((abs(pindel$start_1_pindel[i]-pindel$start_1_pindel[i-1])<=window_between_callers | 
       abs(pindel$start_2_pindel[i]-pindel$start_2_pindel[i-1])<=window_between_callers))){
    
    pindel$ID[i] = paste0(pindel$ID[i-1])
  }
  else{pindel$ID[i] = paste0(type,"_",i)}
  
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

i=1

for(i in 1:nrow(pindel_ok)){
  
  aux = pindel %>% filter(ID==pindel_ok$ID[i])
  
  pindel_ok$chr_1_pindel[i] = unique(aux$chr_1_pindel)
  pindel_ok$start_1_pindel[i] = min(aux$start_1_pindel)
  pindel_ok$end_1_pindel[i] = max(aux$start_1_pindel)
  
  pindel_ok$chr_2_pindel[i] = unique(aux$chr_2_pindel)
  pindel_ok$start_2_pindel[i] = min(aux$start_2_pindel)
  pindel_ok$end_2_pindel[i] = max(aux$start_2_pindel)
  
  pindel_ok$length_pindel[i] = max(c(pindel_ok$end_2_pindel[i]-pindel_ok$start_2_pindel[i],
                                   pindel_ok$end_1_pindel[i]-pindel_ok$start_1_pindel[i]))
  
}

pindel_ok$lower_pindel_start_1 = pindel_ok$start_1_pindel-10
pindel_ok$upper_pindel_start_1 = pindel_ok$start_1_pindel+10
pindel_ok$lower_pindel_start_2 = pindel_ok$start_2_pindel-10
pindel_ok$upper_pindel_start_2 = pindel_ok$start_2_pindel+10

pindel_ok$lower_pindel_end_1 = pindel_ok$end_1_pindel-10
pindel_ok$upper_pindel_end_1 = pindel_ok$end_1_pindel+10
pindel_ok$lower_pindel_end_2 = pindel_ok$end_2_pindel-10
pindel_ok$upper_pindel_end_2 = pindel_ok$end_2_pindel+10


pindel_ok = pindel_ok %>% as.data.table()


## manta #####

manta = fread("data/Manta_TRA",fill = T)
colnames(manta) = c("chr_1_manta","start_1_manta","chr_2_manta","start_2_manta","GT_manta")

manta$chr_1_manta[manta$chr_1_manta=="X"] = "23"
manta$chr_1_manta[manta$chr_1_manta=="Y"] = "24"
manta$chr_2_manta[manta$chr_2_manta=="X"] = "23"
manta$chr_2_manta[manta$chr_2_manta=="Y"] = "24"

manta$chr_1_manta = as.numeric(manta$chr_1_manta)
manta$chr_2_manta = as.numeric(manta$chr_2_manta)

for(i in 1:nrow(manta)){
  
  if(manta$chr_1_manta[i]>manta$chr_2_manta[i]){
    
    chr1 = manta$chr_1_manta[i]
    chr2 = manta$chr_2_manta[i]
    start1 = manta$start_1_manta[i]
    start2 = manta$start_2_manta[i]
    
    manta$chr_1_manta[i] = chr2
    manta$chr_2_manta[i] = chr1
    manta$start_1_manta[i] = start2
    manta$start_2_manta[i] = start1
  }
  
}

manta = unique(manta)

manta = manta %>% arrange(chr_1_manta,chr_2_manta,
                          start_1_manta,start_2_manta) %>% as.data.table()

type = "translocation"

manta$ID = "none"

manta$ID[1] = paste0(type,"_",1)

for(i in 2:nrow(manta)){
  
  if(manta$chr_1_manta[i]==manta$chr_1_manta[i-1] & manta$chr_2_manta[i]==manta$chr_2_manta[i-1] & 
     ((abs(manta$start_1_manta[i]-manta$start_1_manta[i-1])<=window_between_callers | 
       abs(manta$start_2_manta[i]-manta$start_2_manta[i-1])<=window_between_callers)) &
     manta$GT_manta[i]==manta$GT_manta[i-1]){
    manta$ID[i] = paste0(manta$ID[i-1])
  }
  else{manta$ID[i] = paste0(type,"_",i)}
  
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

i=1

for(i in 1:nrow(manta_ok)){
  
  aux = manta %>% filter(ID==manta_ok$ID[i])
  
  manta_ok$chr_1_manta[i] = unique(aux$chr_1_manta)
  manta_ok$start_1_manta[i] = min(aux$start_1_manta)
  manta_ok$end_1_manta[i] = max(aux$start_1_manta)
  
  manta_ok$chr_2_manta[i] = unique(aux$chr_2_manta)
  manta_ok$start_2_manta[i] = min(aux$start_2_manta)
  manta_ok$end_2_manta[i] = max(aux$start_2_manta)
  
  manta_ok$GT_manta[i] = paste0(unique(aux$GT_manta),collapse = "")
  manta_ok$length_manta[i] = max(c(manta_ok$end_2_manta[i]-manta_ok$start_2_manta[i],
                                   manta_ok$end_1_manta[i]-manta_ok$start_1_manta[i]))
  
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



### svaba #####

svaba = fread("data/SVaBA_TRA")
colnames(svaba) = c("chr_1_svaba","start_1_svaba","chr_2_svaba","start_2_svaba","GT_svaba")

svaba$chr_1_svaba[svaba$chr_1_svaba=="X"] = "23"
svaba$chr_1_svaba[svaba$chr_1_svaba=="Y"] = "24"
svaba$chr_2_svaba[svaba$chr_2_svaba=="X"] = "23"
svaba$chr_2_svaba[svaba$chr_2_svaba=="Y"] = "24"

svaba$chr_1_svaba = as.numeric(svaba$chr_1_svaba)
svaba$chr_2_svaba = as.numeric(svaba$chr_2_svaba)

remove_variants = NULL

for(i in 1:nrow(svaba)){
  
  if(svaba$chr_1_svaba[i]>svaba$chr_2_svaba[i]){
    
    chr1 = svaba$chr_1_svaba[i]
    chr2 = svaba$chr_2_svaba[i]
    start1 = svaba$start_1_svaba[i]
    start2 = svaba$start_2_svaba[i]
    
    svaba$chr_1_svaba[i] = chr2
    svaba$chr_2_svaba[i] = chr1
    svaba$start_1_svaba[i] = start2
    svaba$start_2_svaba[i] = start1
  }
  if(svaba$chr_1_svaba[i]==svaba$chr_2_svaba[i]){
    
    remove_variants = c(remove_variants,i)
    
  }
}

svaba = svaba[-remove_variants,]
svaba = unique(svaba)

svaba = svaba %>% arrange(chr_1_svaba,chr_2_svaba,
                          start_1_svaba,start_2_svaba) %>% as.data.table()

type = "translocation"

svaba$ID = "none"

svaba$ID[1] = paste0(type,"_",1)

for(i in 2:nrow(svaba)){
  
  if(svaba$chr_1_svaba[i]==svaba$chr_1_svaba[i-1] & svaba$chr_2_svaba[i]==svaba$chr_2_svaba[i-1] & 
     ((abs(svaba$start_1_svaba[i]-svaba$start_1_svaba[i-1])<=window_between_callers | 
       abs(svaba$start_2_svaba[i]-svaba$start_2_svaba[i-1])<=window_between_callers)) &
     svaba$GT_svaba[i]==svaba$GT_svaba[i-1]){
    svaba$ID[i] = paste0(svaba$ID[i-1])
  }
  else{svaba$ID[i] = paste0(type,"_",i)}
  
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

i=1

for(i in 1:nrow(svaba_ok)){
  
  aux = svaba %>% filter(ID==svaba_ok$ID[i])
  
  svaba_ok$chr_1_svaba[i] = unique(aux$chr_1_svaba)
  svaba_ok$start_1_svaba[i] = min(aux$start_1_svaba)
  svaba_ok$end_1_svaba[i] = max(aux$start_1_svaba)
  
  svaba_ok$chr_2_svaba[i] = unique(aux$chr_2_svaba)
  svaba_ok$start_2_svaba[i] = min(aux$start_2_svaba)
  svaba_ok$end_2_svaba[i] = max(aux$start_2_svaba)
  
  svaba_ok$GT_svaba[i] = paste0(unique(aux$GT_svaba),collapse = "")
  
  print(paste0(i," ",paste0(unique(aux$GT_svaba))),collapse = "")
  
  svaba_ok$length_svaba[i] = max(c(svaba_ok$end_2_svaba[i]-svaba_ok$start_2_svaba[i],
                                   svaba_ok$end_1_svaba[i]-svaba_ok$start_1_svaba[i]))
  
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



### sensitivity-specifitiy-geno error for each caller ######

# lumpy

sens_lumpy = sensitivity_precision_geno_error_bp_trans_bis(lumpy_ok,golden,"translocation","lumpy")

aux = sens_lumpy$all

aux = aux %>% filter(is.na(chr_1_golden))

head(aux)

#write.table(aux[,c(12,14)],"Translocations/lumpy_falsos_positivos",row.names = F,quote = F,sep = "\t")


# pindel

sens_pindel = sensitivity_precision_geno_error_bp_trans_bis(pindel_ok,golden,"translocation","pindel")

aux = sens_pindel$all

aux = aux %>% filter(is.na(chr_1_golden))

head(aux)

#write.table(aux[,c(12,14)],"Translocations/pindel_falsos_positivos",row.names = F,quote = F,sep = "\t")

# manta

sens_manta = sensitivity_precision_geno_error_bp_trans_bis(manta_ok,golden,"translocation","manta")

aux = sens_manta$all

aux = aux %>% filter(is.na(start_1_golden))

head(aux)

#write.table(aux[,c(12,14)],"Translocations/manta_falsos_positivos",row.names = F,quote = F,sep = "\t")

# delly

sens_delly = sensitivity_precision_geno_error_bp_trans_bis(delly_ok,golden,"translocation","delly")

aux = sens_delly$all

aux = aux %>% filter(is.na(chr_1_golden))

head(aux)

#write.table(aux[,c(12,14)],"Translocations/delly_falsos_positivos",row.names = F,quote = F,sep = "\t")


# svaba

sens_svaba = sensitivity_precision_geno_error_bp_trans_bis(svaba_ok,golden,"translocation","svaba")


metrics_callers = matrix(c(as.numeric(t(unlist(sens_delly))[c(3,4,1,2,5,6,7,9,8,10)]),
                           as.numeric(t(unlist(sens_lumpy))[c(3,4,1,2,5,6,7,9,8,10)]),
                           as.numeric(t(unlist(sens_svaba))[c(3,4,1,2,5,6,7,9,8,10)]),
                           as.numeric(t(unlist(sens_manta))[c(3,4,1,2,5,6,7,9,8,10)]),
                           as.numeric(t(unlist(sens_pindel))[c(3,4,1,2,5,6,7,9,8,10)])),byrow=T,ncol = 10)

metrics_callers = cbind(c("Delly","Lumpy","Svaba","Manta","Pindel"),
                        round(metrics_callers,3))

colnames(metrics_callers) = c("Caller","TP","FP","Sensitivity","Precision","F1-Score","Genotype error",
                              "N 0/1","Genotype error 0/1","N 1/1","Genotype error 1/1")

metrics_callers = metrics_callers %>% as.data.frame() %>% arrange(desc(`F1-Score`))

metrics_callers = cbind(N=c(nrow(golden),rep("",nrow(metrics_callers)-1)),
                        metrics_callers)

metrics_callers = cbind(SV_type=c("Translocations",rep("",nrow(metrics_callers)-1)),
                        metrics_callers)

metrics_callers[,1:3] = apply(metrics_callers[,1:3], 2, function(x) as.character(x));
metrics_callers[,4:9] = apply(metrics_callers[,4:9], 2, function(x) as.numeric(x));



## prepare BBDD ####

# manta

full_merge = NULL

for(chromosome in 1:22){

  manta_chr = manta_ok %>% filter(chr_1_manta==chromosome) %>% as.data.table()
  delly_chr = delly_ok %>% filter(chr_1_delly==chromosome) %>% as.data.table()
  lumpy_chr = lumpy_ok %>% filter(chr_1_lumpy==chromosome) %>% as.data.table()
  pindel_chr = pindel_ok %>% filter(chr_1_pindel==chromosome) %>% as.data.table()
  svaba_chr = svaba_ok %>% filter(chr_1_svaba==chromosome) %>% as.data.table()
  golden_chr = golden %>% filter(chr_1_golden==chromosome) %>% as.data.table()
  
  manta_chr$ID = paste0("translocation_",1:nrow(manta_chr))
  delly_chr$ID = "none"
  lumpy_chr$ID = "none"
  svaba_chr$ID = "none"
  pindel_chr$ID = "none"
  golden_chr$ID = "none"
  
  my_data = merge_callers_bp_trans_bis(delly_chr,manta_chr,callers_ref = c("manta"),
                             caller_to_merge = c("delly"),merge = T,golden = F,
                             type="translocation") 
  
  my_data = merge_callers_bp_trans_bis(lumpy_chr,my_data,callers_ref = c("manta","delly"),
                             caller_to_merge = c("lumpy"),merge = T,golden = F,
                             type="translocation") 
  
  my_data = merge_callers_bp_trans_bis(svaba_chr,my_data,callers_ref = c("manta","delly","lumpy"),
                             caller_to_merge = c("svaba"),merge = T,golden = F,
                             type="translocation") 
   
  my_data = merge_callers_bp_trans_bis(pindel_chr,my_data,callers_ref = c("manta","delly","lumpy","svaba"),
                                   caller_to_merge = c("pindel"),merge = T,golden = F,
                                   type="translocation")
  
  my_data = merge_callers_bp_trans_bis(golden_chr,my_data,callers_ref = c("manta","lumpy","delly","svaba","pindel"),
                              caller_to_merge = c("golden"),golden = TRUE,merge=TRUE,
                              type="translocation") 
  
  full_merge = rbind(full_merge,my_data)
  
  print(paste0("chromosome: ",chromosome))

}

for(chromosome in 23){
  
  manta_chr = manta_ok %>% filter(chr_1_manta==chromosome) %>% as.data.table()
  delly_chr = delly_ok %>% filter(chr_1_delly==chromosome) %>% as.data.table()
  lumpy_chr = lumpy_ok %>% filter(chr_1_lumpy==chromosome) %>% as.data.table()
  pindel_chr = pindel_ok %>% filter(chr_1_pindel==chromosome) %>% as.data.table()
  svaba_chr = svaba_ok %>% filter(chr_1_svaba==chromosome) %>% as.data.table()
  golden_chr = golden %>% filter(chr_1_golden==chromosome) %>% as.data.table()
  
  manta_chr$ID = paste0("translocation_",1:nrow(manta_chr))
  delly_chr$ID = "none"
  lumpy_chr$ID = "none"
  svaba_chr$ID = "none"
  pindel_chr$ID = "none"
  golden_chr$ID = "none"
  
  my_data = merge_callers_bp_trans_bis(delly_chr,manta_chr,callers_ref = c("manta"),
                                       caller_to_merge = c("delly"),merge = T,golden = F,
                                       type="translocation") 
  
  my_data = merge_callers_bp_trans_bis(lumpy_chr,my_data,callers_ref = c("manta","delly"),
                                       caller_to_merge = c("lumpy"),merge = T,golden = F,
                                       type="translocation") 

  my_data = cbind(my_data,rbind(svaba_chr,rep(NA,nrow(my_data)),fill=T)[,-c(1,18)])  
  
  my_data = cbind(my_data,rbind(pindel_chr,rep(NA,nrow(my_data)),fill=T)[,-c(1,18)])  

  my_data = merge_callers_bp_trans_bis(golden_chr,my_data,callers_ref = c("manta","lumpy","delly","svaba","pindel"),
                                       caller_to_merge = c("golden"),golden = TRUE,merge=TRUE,
                                       type="translocation") 
  
  full_merge = rbind(full_merge,my_data)
  
  print(paste0("chromosome: ",chromosome))
  
}

table(full_merge$chr_1_golden)
table(full_merge$chr_2_golden)


my_data2 = full_merge

callers_ref = c("manta","lumpy","delly","pindel","svaba","golden")

my_data2$chr_1 = NA

for(i in 1:nrow(my_data2)){
  
  chromosomes = as.numeric(my_data2[i,paste0("chr_1_",callers_ref),with=F])
  
  chromosomes = chromosomes[!is.na(chromosomes)]
  
  my_data2$chr_1[i] = unique(as.numeric(chromosomes))[1]
  
}

table(my_data2$chr_1)


my_data2$chr_1_golden = as.numeric(my_data2$chr_1_golden)

# neteja IDs 

final_data = NULL

chromosomes = 1

for(chromosomes in 1:23){
  
  aux = my_data2 %>% filter(chr_1 == chromosomes)
  
  aux$length_manta[which(aux$length_manta==0)] = NA
  aux$length_lumpy[which(aux$length_lumpy==0)] = NA
  aux$length_pindel[which(aux$length_pindel==0)] = NA
  aux$length_delly[which(aux$length_delly==0)] = NA
  aux$length_svaba[which(aux$length_svaba==0)] = NA
  
  aux = aux %>% group_by(ID) %>% 
    dplyr::summarise(chr_1 = median(chr_1,na.rm = T),
                     start_1 = median(c(start_1_manta,
                                        start_1_lumpy,
                                        start_1_pindel,
                                        start_1_delly,
                                        start_1_svaba),na.rm = T),
                     end_1 = median(c(end_1_manta,
                                      end_1_lumpy,
                                      end_1_pindel,
                                      end_1_delly,
                                      end_1_svaba),na.rm = T),
                     chr_2 = median(c(chr_2_manta,
                                      chr_2_lumpy,
                                      chr_2_pindel,
                                      chr_2_delly,
                                      chr_2_svaba),na.rm = T),
                     start_2 = median(c(start_2_manta,
                                        start_2_lumpy,
                                        start_2_pindel,
                                        start_2_delly,
                                        start_2_svaba),na.rm = T),
                     end_2 = median(c(end_2_manta,
                                      end_2_lumpy,
                                      end_2_pindel,
                                      end_2_delly,
                                      end_2_svaba),na.rm = T),
                     start_1_manta = median(c(start_1_manta),na.rm = T),
                     start_1_lumpy = median(c(start_1_lumpy),na.rm = T),
                     start_1_pindel = median(c(start_1_pindel),na.rm = T),
                     start_1_delly = median(c(start_1_delly),na.rm = T),
                     start_1_svaba = median(c(start_1_svaba),na.rm = T),
                     end_1_manta = median(c(end_1_manta),na.rm = T),
                     end_1_lumpy = median(c(end_1_lumpy),na.rm = T),
                     end_1_pindel = median(c(end_1_pindel),na.rm = T),
                     end_1_delly = median(c(end_1_delly),na.rm = T),
                     end_1_svaba = median(c(end_1_svaba),na.rm = T),
                     start_2_manta = median(c(start_2_manta),na.rm = T),
                     start_2_lumpy = median(c(start_2_lumpy),na.rm = T),
                     start_2_pindel = median(c(start_2_pindel),na.rm = T),
                     start_2_delly = median(c(start_2_delly),na.rm = T),
                     start_2_svaba = median(c(start_2_svaba),na.rm = T),
                     end_2_manta = median(c(end_2_manta),na.rm = T),
                     end_2_lumpy = median(c(end_2_lumpy),na.rm = T),
                     end_2_pindel = median(c(end_2_pindel),na.rm = T),
                     end_2_delly = median(c(end_2_delly),na.rm = T),
                     end_2_svaba = median(c(end_2_svaba),na.rm = T),
                     GT_pindel = names(sort(table(GT_pindel,useNA = "always"),decreasing = T))[1],
                     GT_delly = names(sort(table(GT_delly,useNA = "always"),decreasing = T))[1],
                     GT_lumpy = names(sort(table(GT_lumpy,useNA = "always"),decreasing = T))[1],
                     GT_manta = names(sort(table(GT_manta,useNA = "always"),decreasing = T))[1],
                     GT_svaba = names(sort(table(GT_svaba,useNA = "always"),decreasing = T))[1],
                     length_manta = median(c(length_manta),na.rm = T),
                     length_lumpy = median(c(length_lumpy),na.rm = T),
                     length_pindel = median(c(length_pindel),na.rm = T),
                     length_delly = median(c(length_delly),na.rm = T),
                     length_svaba = median(c(length_svaba),na.rm = T),
                     length = median(c(length_manta,
                                       length_lumpy,
                                       length_pindel,
                                       length_delly,
                                       length_svaba),na.rm = T),
                     chr_1_golden = median(c(chr_1_golden),na.rm = T),
                     start_1_golden = median(c(start_1_golden),na.rm = T),
                     end_1_golden = median(c(end_1_golden),na.rm = T),
                     chr_2_golden = median(c(chr_2_golden),na.rm = T),
                     start_2_golden = median(c(start_2_golden),na.rm = T),
                     end_2_golden = median(c(end_2_golden),na.rm = T),
                     GT_golden = names(sort(table(GT_golden,useNA = "always"),decreasing = T))[1],
                     length_golden = median(length_golden,na.rm = T))
  
  final_data = rbind(final_data,aux)
}

table(final_data$chr_1)
table(final_data$chr_2)


### add length #######

summary(final_data$length)
                              
aux = as.character(cut(final_data$length,
                    breaks = c(30,150,500,1000,2000,3000,Inf),right = F))

aux[which(is.na(aux))] ="unknown"

table(aux)

final_data$length_stretch = aux


## add strategy ####

final_data = final_data %>% as.data.table()

strategies = fread("ext/strategies.csv")

final_data$strategy = strategy(final_data,
                             callers = c("lumpy","manta","pindel","delly","svaba"),
                             strategies)

table(final_data$strategy)


## add number of callers detected #####

final_data$callers_detected = n_callers_detected(final_data,
                                                 callers = c("lumpy","manta","pindel","delly","svaba"))

final_data$callers_detected2 = n_callers_detected(final_data,
                                                 callers = c("lumpy","manta","pindel","delly","svaba"))

final_data$callers_detected[final_data$callers_detected=="4"] = "4-5"
final_data$callers_detected[final_data$callers_detected=="5"] = "4-5"

table(final_data$callers_detected)
table(final_data$callers_detected2)

### colineality #####

table(final_data$callers_detected,final_data$strategy)
chisq.test(table(final_data$callers_detected,final_data$strategy))


### GT ######

final_data$GT_final = NA

data_call = final_data %>%
  dplyr::select(GT_lumpy,
                GT_delly,
                GT_manta,
                GT_final)

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

table(data_call$GT_final,useNA = "always")

final_data$GT_final = as.character(data_call[,4])


## prepare model choose the best model ####

final_data = final_data %>% filter(length>30 | is.na(length))

final_data$PASS = "YES"

final_data$PASS[is.na(final_data$GT_golden)] = "NO"

table(final_data$PASS)
table(final_data$PASS,final_data$callers_detected)
table(final_data$PASS,final_data$length_stretch)

final_data %>% filter(PASS=="NO" & callers_detected==5)


dim(golden)

final_data$GT_manta[is.na(final_data$GT_manta)] = "9/9"
final_data$GT_delly[is.na(final_data$GT_delly)] = "9/9"
final_data$GT_lumpy[is.na(final_data$GT_lumpy)] = "9/9"
final_data$GT_pindel[is.na(final_data$GT_pindel)] = "9/9"
final_data$GT_svaba[is.na(final_data$GT_svaba)] = "9/9"


final_data$GT_manta2 = final_data$GT_manta
final_data$GT_delly2 = final_data$GT_delly
final_data$GT_lumpy2 = final_data$GT_lumpy
final_data$GT_pindel2 = final_data$GT_pindel
final_data$GT_svaba2 = final_data$GT_svaba

final_data$GT_manta[final_data$GT_manta=="./."] = "0/1-1/1"
final_data$GT_manta[final_data$GT_manta=="0/1"] = "0/1-1/1"
final_data$GT_manta[final_data$GT_manta=="1/1"] = "0/1-1/1"

final_data$GT_svaba[final_data$GT_svaba=="./."] = "0/1-1/1"
final_data$GT_svaba[final_data$GT_svaba=="0/1"] = "0/1-1/1"
final_data$GT_svaba[final_data$GT_svaba=="1/1"] = "0/1-1/1"

final_data$GT_delly[final_data$GT_delly=="./."] = "0/1-1/1"
final_data$GT_delly[final_data$GT_delly=="0/1"] = "0/1-1/1"
final_data$GT_delly[final_data$GT_delly=="1/1"] = "0/1-1/1"

final_data$GT_lumpy[final_data$GT_lumpy=="./."] = "0/1-1/1"
final_data$GT_lumpy[final_data$GT_lumpy=="0/1"] = "0/1-1/1"
final_data$GT_lumpy[final_data$GT_lumpy=="1/1"] = "0/1-1/1"

final_data$GT_pindel[final_data$GT_pindel=="./."] = "0/1-1/1"
final_data$GT_pindel[final_data$GT_pindel=="0/1"] = "0/1-1/1"
final_data$GT_pindel[final_data$GT_pindel=="1/1"] = "0/1-1/1"


table(final_data$GT_manta,useNA = "always")
table(final_data$GT_delly,useNA = "always")
table(final_data$GT_lumpy,useNA = "always")
table(final_data$GT_svaba,useNA = "always")
table(final_data$GT_pindel,useNA = "always")

final_data$PASS = ifelse(final_data$PASS=="YES",0,1)

my_model = glm(PASS ~ factor(GT_manta)+
                 factor(GT_delly)+
                 factor(GT_lumpy)+
                 factor(GT_pindel)+
                 factor(GT_svaba) + 
                 factor(callers_detected),
               data=final_data,family = "binomial",
               control = list(maxit = 100))

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

model_sens;model_spec

final_data$prediction = my_pred

## genotype error in the model ######

# 1) Geno error: most common genotype


final_data$GT = NA

data_call = final_data %>% filter(prediction=="PASS") %>% filter(PASS=="YES") %>% 
  dplyr::select(GT_lumpy2,
                GT_delly2,
                GT_manta2,
                GT,GT_golden,PASS,prediction)


data_call$GT_delly2[data_call$GT_delly2=="9/9"] = NA
data_call$GT_lumpy2[data_call$GT_lumpy2=="9/9"] = NA
data_call$GT_manta2[data_call$GT_manta2=="9/9"] = NA

table(data_call$GT_lumpy2)
table(data_call$GT_delly2)
table(data_call$GT_manta2)

data_call = as.matrix(data_call)

for(i in 1:nrow(data_call)){
  
  my_tab = table(as.character(data_call[i,c(1:3)]))
  
  if(length(which(my_tab==max(my_tab)))==1){
    data_call[i,4] = names(which.max(my_tab))
  }
  
  if(length(which(my_tab==max(my_tab)))>1){
    
    data_call[i,4] = NA
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


# 2) order manta-delly-lumpy

caller = c("lumpy","delly","manta")

data_call = as.matrix(data_call)

i=1

for(i in 1:nrow(data_call)){
  
  names_call = caller[which(!is.na(data_call[i,c(1:7)]))]
  
  if(sum(names_call %in% "manta")==1){
    
    data_call[i,4] = data_call[i,3]
  }
  
  
  if(sum(names_call %in% "manta")==0 & sum(names_call %in% "delly")==1){
    
    data_call[i,4] = data_call[i,2]
  }
  
  
  if(sum(names_call %in% "manta")==0 & sum(names_call %in% "delly")==0){
    
    data_call[i,4] = data_call[i,1]
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



### choose best model ####

final_data$PASS = ifelse(final_data$PASS=="YES",0,1)

full.model = glm(PASS ~ factor(GT_manta)+
                   factor(GT_delly)+
                   factor(GT_lumpy)+
                   factor(GT_pindel)+
                   factor(GT_svaba) + 
                   factor(callers_detected),
                 data=final_data,family = "binomial")

base.model = glm(PASS ~ factor(GT_manta)+
                   factor(GT_delly)+
                   factor(GT_lumpy)+
                   factor(GT_pindel)+
                   factor(GT_svaba),
                 data=final_data,family = "binomial")


step(full.model,
     scope = list(lower = formula(base.model),
                  upper = formula(full.model)),
     direction = "backward")

best.model = glm(formula = PASS ~ factor(GT_manta) + factor(GT_delly) + factor(GT_lumpy) + 
                   factor(GT_pindel) + factor(GT_svaba), family = "binomial", 
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


# detected by at least 1 caller

det_call = final_data %>% filter(callers_detected %in% 1:9)

callers1 = logical_sens_spec(det_call)

# detected by at least 2 callers and 2 strategies

det_call = final_data %>% filter(callers_detected2 %in% 2:9 & strategy %in% 2:6)

callers2 = logical_sens_spec(det_call)

# detected by at least 3 callers and 2 strategies

det_call = final_data %>% filter(callers_detected2 %in% 3:9 & strategy %in% 2:6)

callers3 = logical_sens_spec(det_call)

# detected by at least 4 callers and 2 strategies

det_call = final_data %>% filter(callers_detected2 %in% 4:9 & strategy %in% 2:6)

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

## 70-30

0.7*nrow(final_data)

set.seed(2589)

training = sample(1:nrow(final_data),345)

table(final_data$PASS,useNA = "always")

final_data$PASS = as.factor(final_data$PASS)

final_data_70 = final_data[training,]

table(final_data_70$PASS,useNA = "always")

mod_fit <- train(factor(PASS) ~ factor(GT_manta) + factor(GT_delly) + factor(GT_lumpy) + 
                   factor(GT_pindel) + factor(GT_svaba), 
                 data=final_data_70, method="glm", family="binomial",
                 trControl = ctrl)

my_pred = predict(mod_fit, newdata=final_data_70)

my_pred = ifelse(my_pred=="YES",0,1)

my_pred = ifelse(my_pred==0,"PASS","NO PASS")

table(final_data_70$PASS)
table(my_pred,final_data_70$PASS)
table(my_pred)

N = as.numeric(table(final_data_70$PASS)[2])

tp_cv = as.numeric(table(my_pred,final_data_70$PASS)[2,2])
fp_cv = as.numeric(table(my_pred)[2])-tp_cv

sens_call_cv = round(tp_cv/N*100,3)
spec_call_cv = round(tp_cv/(fp_cv+tp_cv)*100,3)

f_score_cv = round((2*sens_call_cv*spec_call_cv)/(sens_call_cv+spec_call_cv)/100,3)

                              
# detected by at least 1 callers and 2 strategies

det_call = final_data_70 %>% filter(callers_detected2 %in% 1:9)

callers1 = logical_sens_spec(det_call)

# detected by at least 2 callers and 2 strategies

det_call = final_data_70 %>% filter(callers_detected2 %in% 2:9 & strategy %in% 2:6)

callers2 = logical_sens_spec(det_call)

# detected by at least 3 callers and 2 strategies

det_call = final_data_70 %>% filter(callers_detected2 %in% 3:9 & strategy %in% 2:6)

callers3 = logical_sens_spec(det_call)

# detected by at least 4 callers and 2 strategies

det_call = final_data_70 %>% filter(callers_detected2 %in% 4:9 & strategy %in% 2:6)

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


final_data_30 = final_data[-training,]

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

det_call = final_data_30 %>% filter(callers_detected2 %in% 1:9)

callers1 = logical_sens_spec(det_call)

# detected by at least 2 callers and 2 strategies

det_call = final_data_30 %>% filter(callers_detected2 %in% 2:9 & strategy %in% 2:6)

callers2 = logical_sens_spec(det_call)

# detected by at least 3 callers and 2 strategies

det_call = final_data_30 %>% filter(callers_detected2 %in% 3:9 & strategy %in% 2:6)

callers3 = logical_sens_spec(det_call)

# detected by at least 4 callers and 2 strategies

det_call = final_data_30 %>% filter(callers_detected2 %in% 4:9 & strategy %in% 2:6)

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


metrics_callers$F1.Score = ((2*metrics_callers$Sensitivity*metrics_callers$Precision)/
                              (metrics_callers$Sensitivity + metrics_callers$Precision))


