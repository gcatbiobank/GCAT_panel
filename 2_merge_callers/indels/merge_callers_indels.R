library(data.table)
library(dplyr)
library(rlang)
library(caret)
library(e1071)
library(R.utils)

args <- commandArgs(trailingOnly = TRUE)
chrm <- args[1]
options(digits = 16,scipen = 16)
strategies = fread("./inputs/strategies.csv")
mod_fit = readRDS("/home/jordi/web/indels/glm_model.rds")
ids  = fread("./inputs/samples",header = F)
#ids= fread("/gpfs/projects/bsc05/jordivalls/GCAT_project_all_samples/merge_all_calling_GCAT/male_samples",header = F) #activate only for male
ids = ids$V1
#q <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
q<- as.numeric(args[2])
print(ids[q])
dir.create(paste0("./outputs/",ids[q]))

## read vcfs Dani files ########

deepvariant = fread(paste0("./inputs/Deepvariant/",ids[q],"/",ids[q],"_",chrm,"_indels.vcf.gz"))
colnames(deepvariant) = c("chr","position","REF","ALT","GT_deepvariant","AD_deep","DP_deep","GQ_deep")
aux = tstrsplit(deepvariant$AD_deep,split=",")
deepvariant$AD_REF_deep = aux[[1]]
deepvariant$AD_ALT_deep = aux[[2]]
deepvariant$ID = do.call(paste0,list(deepvariant$chr,"_",deepvariant$position,"_",deepvariant$REF,"_",deepvariant$ALT))
deepvariant$chr[deepvariant$chr=="X"] = 23     
deepvariant$chr[deepvariant$chr=="Y"] = 24

gatk = fread(paste0("./inputs/GATK/",ids[q],"/",ids[q],"_indel_",chrm,"_indels.vcf.gz"))
colnames(gatk) = c("chr","position","REF","ALT","GT_gatk","AD_gatk","DP_gatk","GQ_gatk")
aux = tstrsplit(gatk$AD_gatk,split=",")
gatk$AD_REF_gatk = aux[[1]]
gatk$AD_ALT_gatk = aux[[2]]
gatk$ID = do.call(paste0,list(gatk$chr,"_",gatk$position,"_",gatk$REF,"_",gatk$ALT))
gatk$chr[gatk$chr=="X"] = 23
gatk$chr[gatk$chr=="Y"] = 24

strelka = fread(paste0("./inputs/Strelka2/",ids[q],"/",ids[q],"_",chrm,"_indels.vcf.gz"))
colnames(strelka) = c("chr","position","REF","ALT","GT_strelka","AD_strelka","DP_strelka","GQ_strelka")
aux = tstrsplit(strelka$AD_strelka,split=",")
strelka$AD_REF_strelka = aux[[1]]
strelka$AD_ALT_strelka = aux[[2]]
strelka$ID = do.call(paste0,list(strelka$chr,"_",strelka$position,"_",strelka$REF,"_",strelka$ALT))
strelka$chr[strelka$chr=="X"] = 23
strelka$chr[strelka$chr=="Y"] = 24

## merge_callers #######

my_data = full_join(deepvariant,gatk)
my_data = full_join(my_data,strelka)


## length indel ###

my_data$length = abs(nchar(my_data$REF)-nchar(my_data$ALT))

## consensus genotype

callers=c("deepvariant","gatk","strelka")
n_callers = length(callers)
data_length= my_data %>% select("GT_deepvariant","GT_gatk", "GT_strelka")
data_length[data_length=="9/9"] = NA

data_length$consensus = NA

data_length = as.matrix(data_length)

for(i in 1:nrow(data_length)){
  
  my_tab = table(as.character(data_length[i,c(1:n_callers)]))
  
  if(length(which(my_tab==max(my_tab)))==1){
    data_length[i,n_callers+1] = names(which.max(my_tab))
  }
  
  if(length(which(my_tab==max(my_tab)))>1){
    
    data_length[i,n_callers+1] = "./."
  }
}

data_length = as.data.frame(data_length)
my_data$GT = as.character(data_length$consensus)


## num callers_detected

callers=c("deepvariant","gatk","strelka")
n_callers = length(callers)
data_length= my_data %>% select("GT_deepvariant","GT_gatk", "GT_strelka")
data_length[data_length=="9/9"] = NA
data_length$consensus = NA

data_length = as.matrix(data_length)
for(i in 1:nrow(data_length)){
  
  my_tab = sum(is.na(data_length[i,1:n_callers]))
  
  data_length[i,n_callers+1] = n_callers-my_tab
  
}

data_length = as.data.frame(data_length)

my_data$callers_detected = as.numeric(as.character.numeric_version(data_length$consensus))

##num strategies

callers=c("deepvariant","gatk","strelka")
n_callers = length(callers)
data_length= my_data %>% select("GT_deepvariant","GT_gatk", "GT_strelka")
data_length[data_length=="9/9"] = NA
data_length$consensus = NA

data_length = as.matrix(data_length)

for(i in 1:nrow(data_length)){
  
  which_na = !is.na(data_length[i,1:n_callers])
  
  my_tab = names(data_length[i,which_na])
  
  my_tab = gsub("GT_","",my_tab)
  
  my_strat = strategies[strategies$CALLER %in% my_tab,2]$STRATEGY
  
  data_length[i,n_callers+1] = length(unique(my_strat))
  
  
}

data_length = as.data.frame(data_length)
my_data$strategy = as.numeric(as.character.numeric_version(data_length$consensus))


## AD_REF_consensus

callers=c("deepvariant","gatk","strelka")
n_callers = length(callers)
data_length= my_data %>% select("AD_REF_deep","AD_REF_gatk", "AD_REF_strelka")
data_length[data_length=="9/9"] = NA
data_length$consensus = NA

data_length = as.matrix(data_length)

for(i in 1:nrow(data_length)){
  my_tab = as.numeric(data_length[i,1:n_callers])
  data_length[i,n_callers+1] = median(my_tab,na.rm = T)
}

data_length = as.data.frame(data_length)
my_data$AD_median_REF_consensus=ceiling(as.numeric(as.character.numeric_version(data_length$consensus)))


##AD_ALT_consensus

callers=c("deepvariant","gatk","strelka")
n_callers = length(callers)
data_length= my_data %>% select("AD_ALT_deep","AD_ALT_gatk", "AD_ALT_strelka")
data_length[data_length=="9/9"] = NA
data_length$consensus = NA


data_length=as.matrix(data_length)
for(i in 1:nrow(data_length)){
  my_tab = as.numeric(data_length[i,1:n_callers])
  data_length[i,n_callers+1] = median(my_tab,na.rm = T)
  
  
}

data_length = as.data.frame(data_length)
my_data$AD_median_ALT_consensus=ceiling(as.numeric(as.character.numeric_version(data_length$consensus)))


## AD_consensus

my_data$AD_consensus= paste0(my_data$AD_median_REF_consensus,",",my_data$AD_median_ALT_consensus) 


## GQ_each_caller

callers=c("deepvariant","gatk","strelka")
n_callers = length(callers)
data_length= my_data %>% select("GQ_deep","GQ_gatk", "GQ_strelka")
data_length[data_length=="9/9"] = NA
data_length$consensus = NA

data_length = as.matrix(data_length)
for(i in 1:nrow(data_length)){
  aux = data_length[i,1:n_callers]
  
  my_GQ = aux[which(data_length[i,1:n_callers]!="9/9")]
  
  data_length[i,n_callers+1] = paste(my_GQ,collapse = ",")
  
}

data_length = as.data.frame(data_length)
my_data$GQ_deep_gatk_strelka=as.character(data_length$consensus)


##DP_consensus

my_data$DP_consensus = my_data$AD_median_REF_consensus + my_data$AD_median_ALT_consensus



## predict using LR model #####

my_data$GT_deepvariant[my_data$GT_deepvariant!="./." & !is.na(my_data$GT_deepvariant)] = "0/1-1/1-2/1"

my_data$GT_gatk[my_data$GT_gatk!="./." & !is.na(my_data$GT_gatk)] = "0/1-1/1-2/1"

my_data$GT_strelka[my_data$GT_strelka!="./." & !is.na(my_data$GT_strelka)] = "0/1-1/1-2/1"

my_data$GT_deepvariant[is.na(my_data$GT_deepvariant)] = "9/9"
my_data$GT_gatk[is.na(my_data$GT_gatk)] = "9/9"
my_data$GT_strelka[is.na(my_data$GT_strelka)] = "9/9"

data_predict = my_data %>% select(ID,chr,position,REF,ALT,GT_deepvariant,GT_gatk,
                                  GT_strelka,GT,callers_detected,strategy,GQ_deep_gatk_strelka,
                                  AD_median_REF_consensus,AD_median_ALT_consensus,DP_consensus,length,AD_consensus) 

data_predict$GT_deepvariant = as.factor(data_predict$GT_deepvariant)
data_predict$GT_gatk = as.factor(data_predict$GT_gatk)
data_predict$GT_strelka = as.factor(data_predict$GT_strelka)


## Filtrar todas las variantes que tengan 0 en AD_median_ALT_consensus ####

data_predict =data_predict[data_predict$AD_median_ALT_consensus!=0, ]


## predict YES/NO threshold 0.5 #######

my_pred = predict(mod_fit,data_predict,type="response")

my_pred = ifelse(my_pred>0.5,1,0)
my_pred = ifelse(my_pred==0,"PASS","NO_PASS")

data_predict$PASS = my_pred


## predict probability (valores entre 0 y 1)

my_pred = predict(mod_fit,data_predict,type = "response") 

data_predict$PASS_num = 1-my_pred 
data_predict = as.data.table(data_predict)

## name callers detected ####

callers=c("deepvariant","gatk","strelka")
n_callers = length(callers)
data_length= data_predict %>% select("GT_deepvariant","GT_gatk", "GT_strelka")
data_length[data_length=="9/9"] = NA
data_length$consensus = NA

data_length = as.matrix(data_length)

for(i in 1:nrow(data_length)){
  
  aux = data_length[i,1:n_callers]
  
  my_names = names(aux)[which(data_length[i,1:n_callers]!="9/9")]
  
  my_names = gsub("GT_","",my_names)
  
  data_length[i,n_callers+1] = paste(my_names,collapse = ",")
  
}

data_length = as.data.frame(data_length)
data_predict$name_callers= as.character(data_length$consensus)


## write final_data ####

final_data = data_predict %>% select(chr,position,REF,ALT,GT,strategy,callers_detected,name_callers,GQ_deep_gatk_strelka,
                                     AD_consensus,DP_consensus,PASS,PASS_num,length) %>%
  arrange(position)


final_data$GT[!final_data$GT %in% c("0/1","1/1","2/1","1/2","1/0")] = "./."

colnames(final_data) = paste0(colnames(final_data),"_",ids[q])

if(chrm =="X"){
  fwrite(final_data,
         paste0("./outputs/",ids[q],"/",ids[q],"_indels_chr_","23"),
         sep = " ",row.names = F,quote = F)
}else if(chrm =="Y"){
  fwrite(final_data,
         paste0("./outputs/",ids[q],"/",ids[q],"_indels_chr_","24"),
         sep = " ",row.names = F,quote = F)
  
}else{
  fwrite(final_data,
         paste0("./outputs/",ids[q],"/",ids[q],"_indels_chr_",chrm),
         sep = " ",row.names = F,quote = F)
  
}
