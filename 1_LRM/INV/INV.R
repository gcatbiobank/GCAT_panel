library(data.table)
library(dplyr)

setwd("/home/jvall2/Mare_folder1/merge_all_calling_GCAT/")
setwd("/home/igalvan/Mare_folder/")

#setwd("D:/BSC/Machine_vector/SV/")

source("Deleciones/functions.R")

# read golden ######

golden_inv = fread("Inversions/Golden/insilico_data/Insilico3_inversiones_final",fill = T)

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


# window size sensitivity and precision #######


# start and end criteria

win_size = c(10,20,50,100,200,300)

my_callers = expand.grid(callers = c("delly","lumpy","pindel","whamg","svaba","manta"),
                         win_size = win_size)

my_callers$sensitivity = 0
my_callers$precision = 0
my_callers$f.score = 0
my_callers$g.error = 0

i=1

for(i in 1:6){
  
  # read vcfs files 
  
  delly = fread("Inversions/Golden/insilico_data/Delly_inversiones")
  delly$V4 = as.numeric(abs(delly$V4))
  colnames(delly) = c("chr","start_delly","end_delly","length_delly","GT_delly")
  delly$lower_delly_bp1 = delly$start_delly-win_size[i]
  delly$upper_delly_bp1 = delly$start_delly+win_size[i]
  delly$lower_delly_bp2 = delly$end_delly-win_size[i]
  delly$upper_delly_bp2 = delly$end_delly+win_size[i]
  
  
  lumpy = fread("Inversions/Golden/insilico_data/Lumpy_inversiones")
  lumpy$V4 = as.numeric(abs(lumpy$V4))
  colnames(lumpy) = c("chr","start_lumpy","end_lumpy","length_lumpy","GT_lumpy")
  lumpy$lower_lumpy_bp1 = lumpy$start_lumpy-win_size[i]
  lumpy$upper_lumpy_bp1 = lumpy$start_lumpy+win_size[i]
  lumpy$lower_lumpy_bp2 = lumpy$end_lumpy-win_size[i]
  lumpy$upper_lumpy_bp2 = lumpy$end_lumpy+win_size[i]
  
  pindel = fread("Inversions/Golden/insilico_data/Pindel_inversiones")
  colnames(pindel) = c("chr","start_pindel","end_pindel","length_pindel","GT_pindel")
  pindel$lower_pindel_bp1 = pindel$start_pindel-win_size[i]
  pindel$upper_pindel_bp1 = pindel$start_pindel+win_size[i]
  pindel$lower_pindel_bp2 = pindel$end_pindel-win_size[i]
  pindel$upper_pindel_bp2 = pindel$end_pindel+win_size[i]
  
  
  whamg = fread("Inversions/Golden/insilico_data/Wham_inversiones",fill = T)
  colnames(whamg) = c("chr","start_whamg","end_whamg","length_whamg","GT_whamg")
  whamg$lower_whamg_bp1 = whamg$start_whamg-win_size[i]
  whamg$upper_whamg_bp1 = whamg$start_whamg+win_size[i]
  whamg$lower_whamg_bp2 = whamg$end_whamg-win_size[i]
  whamg$upper_whamg_bp2 = whamg$end_whamg+win_size[i]
  
  manta = fread("Inversions/Golden/insilico_data/Manta_inversiones",fill = T)
  manta$V4 = as.numeric(abs(manta$V4))
  colnames(manta) = c("chr","start_manta","end_manta","length_manta","GT_manta")
  manta$lower_manta_bp1 = manta$start_manta-win_size[i]
  manta$upper_manta_bp1 = manta$start_manta+win_size[i]
  manta$lower_manta_bp2 = manta$end_manta-win_size[i]
  manta$upper_manta_bp2 = manta$end_manta+win_size[i]
  
  
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
  
  svaba = svaba[,c(1:4,6)]
  colnames(svaba) = c("chr","start_svaba","end_svaba","GT_svaba","length_svaba")
  svaba$lower_svaba_bp1 = svaba$start_svaba-win_size[i]
  svaba$upper_svaba_bp1 = svaba$start_svaba+win_size[i]
  svaba$lower_svaba_bp2 = svaba$end_svaba-win_size[i]
  svaba$upper_svaba_bp2 = svaba$end_svaba+win_size[i]
  
  svaba = unique(svaba)
  
  # lumpy
  
  sens_lumpy = sensitivity_precision_geno_error_inv(lumpy,golden,"inversion","lumpy")
  
  my_callers[(my_callers$callers=="lumpy" & my_callers$win_size==win_size[i]),]$sensitivity = sens_lumpy$sensitivity
  my_callers[(my_callers$callers=="lumpy" & my_callers$win_size==win_size[i]),]$precision = sens_lumpy$precision
  my_callers[(my_callers$callers=="lumpy" & my_callers$win_size==win_size[i]),]$f.score = sens_lumpy$f1_score
  my_callers[(my_callers$callers=="lumpy" & my_callers$win_size==win_size[i]),]$g.error = sens_lumpy$g_error
  
  
  # pindel
  
  sens_pindel = sensitivity_precision_geno_error_inv(pindel,golden,"inversion","pindel")
  
  my_callers[(my_callers$callers=="pindel" & my_callers$win_size==win_size[i]),]$sensitivity = sens_pindel$sensitivity
  my_callers[(my_callers$callers=="pindel" & my_callers$win_size==win_size[i]),]$precision = sens_pindel$precision
  my_callers[(my_callers$callers=="pindel" & my_callers$win_size==win_size[i]),]$f.score = sens_pindel$f1_score
  my_callers[(my_callers$callers=="pindel" & my_callers$win_size==win_size[i]),]$g.error = sens_pindel$g_error
  
  # manta
  
  sens_manta = sensitivity_precision_geno_error_inv(manta,golden,"inversion","manta")
  
  my_callers[(my_callers$callers=="manta" & my_callers$win_size==win_size[i]),]$sensitivity = sens_manta$sensitivity
  my_callers[(my_callers$callers=="manta" & my_callers$win_size==win_size[i]),]$precision = sens_manta$precision
  my_callers[(my_callers$callers=="manta" & my_callers$win_size==win_size[i]),]$f.score = sens_manta$f1_score
  my_callers[(my_callers$callers=="manta" & my_callers$win_size==win_size[i]),]$g.error = sens_manta$g_error
  
  # svaba
  
  sens_svaba = sensitivity_precision_geno_error_inv(svaba,golden,"inversion","svaba")
  
  my_callers[(my_callers$callers=="svaba" & my_callers$win_size==win_size[i]),]$sensitivity = sens_svaba$sensitivity
  my_callers[(my_callers$callers=="svaba" & my_callers$win_size==win_size[i]),]$precision = sens_svaba$precision
  my_callers[(my_callers$callers=="svaba" & my_callers$win_size==win_size[i]),]$f.score = sens_svaba$f1_score
  my_callers[(my_callers$callers=="svaba" & my_callers$win_size==win_size[i]),]$g.error = sens_svaba$g_error
  
  # whamg
  
  sens_whamg = sensitivity_precision_geno_error_inv(whamg,golden,"inversion","whamg")
  
  my_callers[(my_callers$callers=="whamg" & my_callers$win_size==win_size[i]),]$sensitivity = sens_whamg$sensitivity
  my_callers[(my_callers$callers=="whamg" & my_callers$win_size==win_size[i]),]$precision = sens_whamg$precision
  my_callers[(my_callers$callers=="whamg" & my_callers$win_size==win_size[i]),]$f.score = sens_whamg$f1_score
  my_callers[(my_callers$callers=="whamg" & my_callers$win_size==win_size[i]),]$g.error = sens_whamg$g_error
  
  # delly
  
  sens_delly = sensitivity_precision_geno_error_inv(delly,golden,"inversion","delly")
  
  my_callers[(my_callers$callers=="delly" & my_callers$win_size==win_size[i]),]$sensitivity = sens_delly$sensitivity
  my_callers[(my_callers$callers=="delly" & my_callers$win_size==win_size[i]),]$precision = sens_delly$precision
  my_callers[(my_callers$callers=="delly" & my_callers$win_size==win_size[i]),]$f.score = sens_delly$f1_score
  my_callers[(my_callers$callers=="delly" & my_callers$win_size==win_size[i]),]$g.error = sens_delly$g_error

  
}

write.csv(my_callers,"Inversions/Golden/outputs/window_size_callers_start_or_end_inversions_criteria.csv",row.names = F)



## plot results ####

my_callers = fread("Inversions/Golden/outputs/window_size_callers_start_or_end_inversions_criteria.csv")

delly = my_callers %>% filter(callers=="delly") #20
lumpy = my_callers %>% filter(callers=="lumpy") #50
whamg = my_callers %>% filter(callers=="whamg") #10
svaba = my_callers %>% filter(callers=="svaba") #10
manta = my_callers %>% filter(callers=="manta") #10
pindel = my_callers %>% filter(callers=="pindel") #10



## read vcf ######

delly = fread("Inversions/Golden/insilico_data/Delly_inversiones")
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

lumpy = fread("Inversions/Golden/insilico_data/Lumpy_inversiones")
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


pindel = fread("Inversions/Golden/insilico_data/Pindel_inversiones")
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


whamg = fread("Inversions/Golden/insilico_data/Wham_inversiones",fill = T)
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

manta = fread("Inversions/Golden/insilico_data/Manta_inversiones",fill = T)
manta$V4 = as.numeric(abs(manta$V4))
colnames(manta) = c("chr","start_manta","end_manta","length_manta","GT_manta")
manta = unify_breakpoints(manta,"manta")
colnames(manta)[2:5] = c("chr","start_manta","end_manta","GT_manta")
manta$length_manta = manta$end_manta-manta$start_manta
manta$lower_manta_bp1 = manta$start_manta-10
manta$upper_manta_bp1 = manta$start_manta+10
manta$lower_manta_bp2 = manta$end_manta-10
manta$upper_manta_bp2 = manta$end_manta+10

dim(manta)
manta = unique(manta)
manta = manta %>% as.data.table()
manta$ID = NULL

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

#write.csv(metrics_callers,"Inversions/Golden/outputs/inversions.csv",row.names = F)

metrics_callers = read.csv("Inversions/Golden/outputs/inversions.csv")

## prepare BBDD ####

# (neteja de cada vcf, el primer vcf es fa fora de la funció, la resta es fa dins) IDs unic per posicions properes, si te dos genotips posem dos linies, si te un unic genotip comu,
# posem una linea

#delly = delly %>% filter(chr==1) %>% as.data.table()
#lumpy = lumpy %>% filter(chr==1) %>% as.data.table()
#manta = manta %>% filter(chr==1) %>% as.data.table()
#whamg = whamg %>% filter(chr==1) %>% as.data.table()
#pindel = pindel %>% filter(chr==1) %>% as.data.table()
#svaba = svaba %>% filter(chr==1) %>% as.data.table()
#golden = golden %>% filter(chr==1) %>% as.data.table()


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

#saveRDS(my_data2,"Inversions/Golden/outputs/data_merge.rds")

my_data2 = readRDS("Inversions/Golden/outputs/data_merge.rds")

my_data2 = as.data.table(my_data2)

### start, end, length, 

# neteja bis_bis abans d'entrar al model i posar minim o maxim 

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

strategies = fread("Deleciones/strategies.csv")

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

#saveRDS(final_data,"Inversions/Golden/outputs/data_merge_inversiones.rds")

final_data = readRDS("Inversions/Golden/outputs/data_merge_inversiones.rds")


### name callers detected and prediction #####

final_data$name_callers_detected = name_callers_detected(as.data.table(final_data), 
                                                         callers = c("whamg","lumpy","manta",
                                                                     "pindel","delly","svaba"))

table(final_data$name_callers_detected,final_data$prediction) 


### genotype error in the model ######

# 1) classic geno error

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


# 2)  Escala Lumpy-Pindel-Whamg-Delly-Manta (La mejor opcion!)

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


# que lo detecte al menos 1 caller

det_call = final_data 

callers1 = logical_sens_spec(det_call)

# que lo detecte al menos 2 caller y 2 estrategias a la vez

det_call = final_data %>% filter(callers_detected_ok %in% c(2:6) &  strategy_ok %in% c(2:4) & reciprocity >0.80)

callers2 = logical_sens_spec(det_call)

# que lo detecte al menos 3 caller

det_call = final_data %>% filter(callers_detected_ok %in% c(3:6) & strategy_ok %in% c(2:4) & reciprocity >0.80)

callers3 = logical_sens_spec(det_call)

# que lo detecte al menos 4 caller

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


#write.csv(metrics_callers,"Inversions/Golden/outputs/inversions.csv",row.names = F)

### CV ######

metrics_callers = read.csv("Inversions/Golden/outputs/inversions.csv")
metrics_callers$Caller = as.character(metrics_callers$Caller)

library(caret)
library(e1071)

ctrl <- trainControl(method = "repeatedcv", number = 10, 
                     savePredictions = TRUE)


## 90-10 (no hi ha prou dades, fem 100%)

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


# que lo detecte al menos 1 caller

det_call = final_data_70 %>% filter(callers_detected_ok %in% 1:9)

callers1 = logical_sens_spec(det_call)

# que lo detecte al menos 2 caller y 2 estrategias a la vez

det_call = final_data_70 %>% filter(callers_detected_ok %in% 2:6 & strategy_ok %in% 2:4)

det_call = det_call[which(det_call$reciprocity>0.80),]

callers2 = logical_sens_spec(det_call)

# que lo detecte al menos 3 caller

det_call = final_data_70 %>% filter(callers_detected_ok %in% 3:9 & strategy_ok %in% 2:6)

det_call = det_call[which(det_call$reciprocity>0.80),]

callers3 = logical_sens_spec(det_call)

# que lo detecte al menos 4 caller

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


#write.csv(metrics_callers,"Inversions/Golden/outputs/inversions.csv",row.names = F)

metrics_callers = read.csv("Inversions/Golden/outputs/inversions.csv")

#saveRDS(mod_fit,"Inversions/Golden/model_inversions.rds")

mod_fit = readRDS("Inversions/Golden/model_inversions.rds")


## no hi ha test set

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


# que lo detecte al menos 1 caller

det_call = final_data_30 %>% filter(callers_detected %in% 1:9)

callers1 = logical_sens_spec(det_call)

# que lo detecte al menos 2 caller y 2 estrategias a la vez

det_call = final_data_30 %>% filter(callers_detected %in% 2:9 & strategy %in% 2:6 & reprocicity>0.8)

callers2 = logical_sens_spec(det_call)

# que lo detecte al menos 3 caller

det_call = final_data_30 %>% filter(callers_detected %in% 3:9 & strategy %in% 2:6 & reprocicity>0.8)

callers3 = logical_sens_spec(det_call)

# que lo detecte al menos 4 caller

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

head(metrics_callers)
library(xtable)
print(xtable(metrics_callers[,-1],digits = 2),include.rownames = F)



### PRUEBAS ########

### glmnet #####

final_data$reciprocity[is.na(final_data$reciprocity)] = 1

#install.packages("glmnet", repos = "http://cran.us.r-project.org")

library(glmnet)

final_data$PASS = ifelse(final_data$PASS=="YES",0,1)

x = final_data %>% select(GT_whamg,GT_pindel,GT_delly,GT_lumpy,GT_manta,GT_svaba,length_stretch,strategy,reciprocity)

y = final_data$PASS

x_train <- model.matrix( ~ .-1, x)

fit = glmnet(x_train, y, family = "binomial")

plot(fit, xvar = "dev", label = TRUE)

cvfit = cv.glmnet(x_train, y, family = "binomial", type.measure = "class")

plot(cvfit)

cvfit$lambda.min
cvfit$lambda.1se
cvfit$lambda

my_pred = predict(cvfit, newx = x_train, s = "lambda.min", type = "class")

my_pred = ifelse(my_pred==1,"NO PASS","PASS")

y = ifelse(y==0,"YES","NO")

table(my_pred,y)

table(my_pred)

sens = 236/283*100

spec = 236/243*100

((2*sens*spec)/(sens + spec))/100

foldid=sample(1:10,size=length(y),replace=TRUE)
cv1=cv.glmnet(x_train,y,foldid=foldid,alpha=1, family = "binomial")
cv.5=cv.glmnet(x_train,y,foldid=foldid,alpha=.5, family = "binomial")
cv0=cv.glmnet(x_train,y,foldid=foldid,alpha=0, family = "binomial")

par(mfrow=c(2,2))
plot(cv1);plot(cv.5);plot(cv0)
plot(log(cv1$lambda),cv1$cvm,pch=19,col="red",xlab="log(Lambda)",ylab=cv1$name)
points(log(cv.5$lambda),cv.5$cvm,pch=19,col="grey")
points(log(cv0$lambda),cv0$cvm,pch=19,col="blue")
legend("topleft",legend=c("alpha= 1","alpha= .5","alpha 0"),pch=19,col=c("red","grey","blue"))

final_data$prediction = my_pred

final_data$PASS = ifelse(final_data$PASS==0,"YES","NO")

# geno error

data_call = final_data %>% filter(prediction=="PASS") %>% filter(PASS=="YES") %>% 
  dplyr::select(GT_whamg,
                GT_lumpy,
                GT_pindel,
                GT_delly,
                GT_manta,
                GT_golden,PASS,prediction)

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
    data_call[i,9] = names(which.max(my_tab))
  }
  
  if(length(which(my_tab==max(my_tab)))>1){
    
    data_call[i,9] = NA
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


## LASSO 70%

training$PASS = ifelse(training$PASS=="YES",0,1)

x = training %>% select(GT_whamg,GT_pindel,GT_delly,GT_lumpy,GT_manta,GT_svaba,length_stretch,strategy,reciprocity)

y = training$PASS

x_train <- model.matrix( ~ .-1, x)

fit = glmnet(x_train, y, family = "binomial")

cvfit = cv.glmnet(x_train, y, family = "binomial", type.measure = "class")

my_pred = predict(cvfit, newx = x_train, s = "lambda.min", type = "class")

my_pred = ifelse(my_pred==1,"NO PASS","PASS")

y = ifelse(y==0,"YES","NO")

table(my_pred,y)

table(my_pred)

# without CHR and BP
#sensitivity
model_sens = (table(my_pred,training$PASS)[2,1]/table(training$PASS)[1])*100
#specificity
model_spec = (table(my_pred,training$PASS)[2,1])/(table(my_pred,training$PASS)[2,2]+
                                                    table(my_pred,training$PASS)[2,1])*100

f_score = ((2*model_sens*model_spec)/(model_sens + model_spec))/100


c(model_sens,model_spec,f_score)

training$prediction = my_pred

training$PASS = ifelse(training$PASS==0,"YES","NO")


# geno error

data_call = training %>% filter(prediction=="PASS") %>% filter(PASS=="YES") %>% 
  dplyr::select(GT_whamg,
                GT_lumpy,
                GT_pindel,
                GT_delly,
                GT_manta,
                GT_golden,PASS,prediction)

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
    data_call[i,9] = names(which.max(my_tab))
  }
  
  if(length(which(my_tab==max(my_tab)))>1){
    
    data_call[i,9] = NA
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


## LASSO 30%

validation$PASS = ifelse(validation$PASS=="YES",0,1)

x = validation %>% select(GT_whamg,GT_pindel,GT_delly,GT_lumpy,GT_manta,GT_svaba,length_stretch,strategy,reciprocity)

y = validation$PASS

x_train <- model.matrix( ~ .-1, x)

fit = glmnet(x_train, y, family = "binomial")

cvfit = cv.glmnet(x_train, y, family = "binomial", type.measure = "class")

my_pred = predict(cvfit, newx = x_train, s = "lambda.min", type = "class")

my_pred = ifelse(my_pred==1,"NO PASS","PASS")

y = ifelse(y==0,"YES","NO")

table(my_pred,y)

table(my_pred)

# without CHR and BP
#sensitivity
model_sens = (table(my_pred,validation$PASS)[2,1]/table(validation$PASS)[1])*100
#specificity
model_spec = (table(my_pred,validation$PASS)[2,1])/(table(my_pred,validation$PASS)[2,2]+
                                                    table(my_pred,validation$PASS)[2,1])*100

f_score = ((2*model_sens*model_spec)/(model_sens + model_spec))/100


c(model_sens,model_spec,f_score)

validation$prediction = my_pred

validation$PASS = ifelse(validation$PASS==0,"YES","NO")


# geno error

data_call = validation %>% filter(prediction=="PASS") %>% filter(PASS=="YES") %>% 
  dplyr::select(GT_whamg,
                GT_lumpy,
                GT_pindel,
                GT_delly,
                GT_manta,
                GT_golden,PASS,prediction)

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
    data_call[i,9] = names(which.max(my_tab))
  }
  
  if(length(which(my_tab==max(my_tab)))>1){
    
    data_call[i,9] = NA
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



#### SVM ##########

saveRDS(training,"Inversiones/outputs/training.rds")
saveRDS(validation,"Inversiones/outputs/validation.rds")
saveRDS(final_data,"Inversiones/outputs/final_data.rds")

install.packages("penalizedSVM")
library(penalizedSVM)


final_data$PASS = ifelse(final_data$PASS=="YES",-1,1)

x = final_data %>% dplyr::select(GT_whamg,GT_pindel,GT_delly,GT_lumpy,GT_manta,GT_svaba,length_stretch,strategy,reciprocity)

y = final_data$PASS

x_train <- model.matrix( ~ .-1, x)

Lambda.scad <- c(0.01 ,0.05)
Lambda.scad <- seq(0.01 ,0.05, 0.01)


## 1 Norm "1norm" LASSO

fit.1norm <- svmfs(x=x_train,y=y, fs.method="1norm",
                   maxIter = 700,cross.inner =50,
                  lambda1.set=Lambda.scad,grid.search = "interval")

test.error.1norm<-predict(fit.1norm, newdata=x_train,newdata.labels=y )
(test.error.1norm$tab)



## "scad+L2" for Elastic SCAD

fit.scad.l2 <- svmfs(x=x_train,y=y, fs.method="scad+L2", 
                     maxIter = 10,
                     lambda1.set=Lambda.scad)

test.error.scad.l2<-predict(fit.scad.l2, newdata=x_train,newdata.labels=y )
(test.error.scad.l2$tab)




#### Plot graphics ##########

library(data.table)
library(ggplot2)
library(plyr)

# read golden ######

golden_inv = fread("Inversiones/Insilico3_inversiones_final",fill = T)

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

aux = cut(golden$length_golden,
          breaks = c(0,150,300,500,1000,2000,3000,Inf),right = F)
table(aux)
golden$length_stretch = aux

golden_length = golden %>% group_by(length_stretch) %>% 
  dplyr::summarise(golden = n()) %>% arrange(length_stretch) %>% as.data.frame()


### model #####

my_data2 = readRDS("Inversiones/outputs/data_merge_inversiones.rds")

aux = cut(my_data2$length_golden,
          breaks = c(0,150,300,500,1000,2000,3000,Inf),right = F)

aux2 = cut(my_data2$length,
           breaks = c(0,150,300,500,1000,2000,3000,Inf),right = F)

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

precision_model$precision = precision_model$truepos/(true_falsepositive_model$true_falsepos)*100 

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
  
  
  aux = my_data[as.character(t(my_data[,paste0("GT_",caller_name)]))!="9/9",]
  
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

sens_precision = sensitivity_precision_caller(my_data2,golden_length,caller_name ="lumpy")
sens_precision_final = rbind(sens_precision_final,sens_precision)


### plot #####

sens_precision_final$sensitivity = -sens_precision_final$sensitivity

sens_precision_final$caller = as.factor(sens_precision_final$caller)

sens_precision_final$color = factor(sens_precision_final$caller)

sens_precision_final$color = mapvalues(sens_precision_final$color , 
                                       from = levels(sens_precision_final$caller), 
                                       to = c("firebrick2","black",
                                              "gray","darkorange2","darkorchid2","deeppink2","darkgoldenrod2"))


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
# [12] "cnvnator" darkolivegreen2
# [13] "genomestrip" blue3


breaks = c(50,150,300,500,1000,2000,3000)
breaks_labels = c("0-150","150-300","300-500","500-1K","1K-2K","2K-3K","3K->3K")

sens_precision_final$breaks = breaks

sens_precision_final[sens_precision_final==0]=NA

model_in = sens_precision_final %>% filter(caller %in% "Logistic regression model")

sens_precision_final[sens_precision_final$caller  %in% "Logistic regression model", ]$sensitivity = NA
sens_precision_final[sens_precision_final$caller  %in% "Logistic regression model", ]$precision = NA

p1 = ggplot(sens_precision_final, aes(x = breaks, y = precision,  color = caller, group = caller)) + 
  geom_point() + geom_line() +
  scale_color_manual(values = levels(factor(sens_precision_final$color))) + theme_classic() +
  ggtitle(c("283 Inversions"))+
  scale_x_log10(breaks = breaks,
                labels = breaks_labels,
                limits = c(50,3000))+ 
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



pdf("Inversiones//outputs/plots/Inversions.pdf",width = 11,height = 5)
p1 +  geom_line(sens_precision_final,mapping = aes(x = breaks, y = sensitivity,  color = caller)) + 
  geom_point(sens_precision_final,mapping = aes(x = breaks, y = sensitivity,  color = caller)) +
  geom_hline(yintercept = 0,linetype="dashed") + ylab("Sensitivity            (%)            Precision") +
  geom_point( aes(x = breaks, y = sensitivity,  color = caller, group = caller),model_in,size=2) +
  geom_line( aes(x = breaks, y = sensitivity,  color = caller, group = caller),model_in,size=1.2)
dev.off()



