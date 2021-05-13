# hacemos insertions juntas porque manta, pindel te dan 2 archivos hasta 50bp, 
# no podemos saber el tamaño hasta 150 para separar 30-150 i >150

library(data.table)
library(dplyr)

setwd("/run/media/igalvan/TOSHIBA EXT/BSC/Machine_vector/SV")
setwd("D:/BSC/Machine_vector/SV")

source("functions/functions.R")

# read golden ######

golden_ins = fread("Insertions/data/Insilico3_inserciones_grandes_virus",fill=T)

golden_ins = golden_ins[,c(1,3,4),with=F]

golden_ins$V4 = nchar(golden_ins$V4)

summary(golden_ins$V4)

golden_ins$V1 = gsub("chr","",golden_ins$V1)

golden_ins2 = fread("Insertions/data/Insilico3_inserciones_grandes_virus",fill=T,skip = 180)

golden_ins2 = golden_ins2[,c(1,3,4),with=F]

golden_ins2$V4 = nchar(golden_ins2$V4)

summary(golden_ins2$V4)

golden_ins2$V1 = gsub("chr","",golden_ins2$V1)

golden_ins = rbind(golden_ins,golden_ins2)

golden_geno_ins = NULL

for(i in names(table(golden_ins$V1))){
  
  aux = golden_ins[golden_ins$V1==i,]
  
  which_duplicated = aux$V3[which(duplicated(aux$V3))]
  
  golden_dup = unique(aux[aux$V3 %in% which_duplicated,])
  golden_no_dup = aux[!aux$V3 %in% which_duplicated,]
  
  golden_dup$genotype = "1/1"
  golden_no_dup$genotype = "0/1"
  
  aux_geno = rbind(golden_dup,golden_no_dup)
  
  golden_geno_ins = rbind(golden_geno_ins,aux_geno)
}

golden = golden_geno_ins
colnames(golden) = c("chr","start_golden","length_golden","GT_golden")
golden$lower_golden = golden$start_golden-0
golden$upper_golden = golden$start_golden+0
golden$chr_pos_golden = do.call(paste0,list(golden$chr,"_",golden$start_golden))





# window size sensitivity and precision #######

win_size = c(10,20,50,100,200,300)

my_callers = expand.grid(callers = c("delly","popins","pindel","whamg","svaba","manta",
                                     "strelka","gatk"),
                         win_size = win_size)

my_callers$sensitivity = 0
my_callers$precision = 0
my_callers$f.score = 0
my_callers$g.error = 0


for(i in 1:6){
  
  # read vcfs files 
  
  delly = fread("Insertions/data/Delly_inserciones_30_50")
  delly$V4 = as.numeric(abs(delly$V4))
  delly = delly %>% filter(V4>30)
  colnames(delly) = c("chr","start_delly","end","length_delly","GT_delly")
  delly$lower_delly = delly$start_delly-win_size[i]
  delly$upper_delly = delly$start_delly+win_size[i]
  delly$chr_pos_delly = do.call(paste0,list(delly$chr,"_",delly$start_delly))
  summary(delly$length_delly)
  
  # pindel no da genotipo ni length para SV
  
  pindel = fread("Insertions/data/Pindel_inserciones_totales",fill = T)
  pindel$V4 = as.numeric(nchar(pindel$V4))
  colnames(pindel) = c("chr","start_pindel","end","length_pindel","GT_pindel")
  #pindel = pindel %>% filter(length_pindel>30 & length_pindel<151)
  pindel$lower_pindel = pindel$start_pindel-win_size[i]
  pindel$upper_pindel = pindel$start_pindel+win_size[i]
  pindel$chr_pos_pindel = do.call(paste0,list(pindel$chr,"_",pindel$start_pindel))
  summary(pindel$length_pindel)
  pindel$chr = as.character(pindel$chr)
  pindel$GT_pindel[pindel$GT_pindel==""] = "0/1"
  
  
  whamg = fread("Insertions/data/Wham_inserciones_totales")
  #whamg = whamg %>% filter(V3>30 & V3<151)
  whamg = whamg[,-2,with=F]
  colnames(whamg) = c("chr","start_whamg","length_whamg")
  whamg$lower_whamg = whamg$start_whamg-win_size[i]
  whamg$upper_whamg = whamg$start_whamg+win_size[i]
  whamg$chr_pos_whamg = do.call(paste0,list(whamg$chr,"_",whamg$start_whamg))
  whamg$GT_whamg = "0/1"
  summary(whamg$length_whamg)
  
  whamg = whamg %>% arrange(chr,chr_pos_whamg)
  head(whamg,40)
  whamg$start_whamg = as.numeric(whamg$start_whamg)
  whamg$diff = c(0,diff(whamg$start_whamg))
  whamg = whamg %>% filter(diff<0 | diff>win_size[i])
  
  
  #svaba no da length para SV
  svaba = fread("Insertions/data/Svaba_insercionesmedianas_totalSV",fill=T,skip = 5289)
  svaba$length_svaba = abs(nchar(svaba$V3)-nchar(svaba$V4))
  colnames(svaba)[c(1,2,5)] = c("chr","start_svaba","GT_svaba")
  svaba = svaba[,c(1,2,5,6)]
  svaba2 = fread("Insertions/data/Svaba_insercionesmedianas_totalSV",fill=T)
  colnames(svaba2) = c("chr","start_svaba","GT_svaba")
  svaba2$length_svaba = 0
  svaba = rbind(svaba,svaba2)
  svaba$lower_svaba = svaba$start_svaba-win_size[i]
  svaba$upper_svaba = svaba$start_svaba+win_size[i]
  svaba$chr_pos_svaba = do.call(paste0,list(svaba$chr,"_",svaba$start_svaba))
  summary(svaba$length_svaba)
  
  svaba = svaba %>% arrange(chr,start_svaba)
  head(svaba,40)
  svaba$start_svaba = as.numeric(svaba$start_svaba)
  svaba$diff = c(0,diff(svaba$start_svaba))
  svaba = svaba %>% filter(diff<0 | diff>win_size[i])
  
  gatk = fread("Insertions/data/Gatk_inserciones_30_150")
  gatk$length_gatk = abs(nchar(gatk$V3)-nchar(gatk$V4))
  #gatk = gatk %>% filter(length_gatk>30 & length_gatk<151)
  colnames(gatk)[c(1,2,5)] = c("chr","start_gatk","GT_gatk")
  gatk = gatk[,c(1,2,5,6)]
  gatk$end = gatk$length_gatk+gatk$start_gatk
  gatk$lower_gatk = gatk$start_gatk-win_size[i]
  gatk$upper_gatk = gatk$start_gatk+win_size[i]
  gatk$chr_pos_gatk = do.call(paste0,list(gatk$chr,"_",gatk$start_gatk))
  summary(gatk$length_gatk)
  
  strelka = fread("Insertions/data/Strelka_inserciones_30_50")
  strelka$length_strelka = abs(nchar(strelka$V3)-nchar(strelka$V4))
  strelka = strelka %>% filter(length_strelka>30 & length_strelka<151)
  colnames(strelka)[c(1,2,5)] = c("chr","start_strelka","GT_strelka")
  strelka = strelka[,c(1,2,5,6)]
  strelka$end = strelka$length_strelka+strelka$start_strelka
  strelka$lower_strelka = strelka$start_strelka-win_size[i]
  strelka$upper_strelka = strelka$start_strelka+win_size[i]
  strelka$chr_pos_strelka = do.call(paste0,list(strelka$chr,"_",strelka$start_strelka))
  summary(strelka$length_strelka)
  
  
  #manta no da length para SV y no da genotipo para pequeñas
  
  manta = fread("Insertions/data/Manta_inserciones_totales")
  manta$V4 = as.numeric(abs(manta$V4))
  #manta = manta %>% filter(V4>30 & V4<151)
  manta$V3 = as.character(manta$V3)
  manta$V3 = "0/1"
  colnames(manta) = c("chr","start_manta","GT_manta","length_manta")
  manta2 = fread("Insertions/data/Manta_inserciones_totales",skip=147)
  head(manta2)
  colnames(manta2) = c("chr","start_manta","GT_manta")
  manta2$length_manta = 0
  manta = rbind(manta,manta2)
  manta$lower_manta = manta$start_manta-win_size[i]
  manta$upper_manta = manta$start_manta+win_size[i]
  manta$chr_pos_manta = do.call(paste0,list(manta$chr,"_",manta$start_manta))
  
  manta = manta %>% arrange(chr,start_manta)
  head(manta,40)
  manta$start_manta = as.numeric(manta$start_manta)
  manta$diff = c(0,diff(manta$start_manta))
  manta = manta %>% filter(diff<0 | diff>win_size[i])
  
  
  
  #popins no da length
  popins = fread("Insertions/data/Popins_inserciones_totales")
  colnames(popins) = c("chr","start_popins","GT_popins")
  popins$lower_popins = popins$start_popins-win_size[i]
  popins$upper_popins = popins$start_popins+win_size[i]
  popins$chr_pos_popins = do.call(paste0,list(popins$chr,"_",popins$start_popins))
  popins$length_popins = 0
  summary(popins$length_popins)
  
  popins = popins %>% arrange(chr,start_popins)
  head(popins,40)
  popins$start_popins = as.numeric(popins$start_popins)
  popins$diff = c(0,diff(popins$start_popins))
  popins = popins %>% filter(diff<0 | diff>300)
  
  
  
  
  # popins
  
  sens_popins = sensitivity_precision_geno_error(popins,golden,"deletion","popins")
  
  my_callers[(my_callers$callers=="popins" & my_callers$win_size==win_size[i]),]$sensitivity = sens_popins$sensitivity
  my_callers[(my_callers$callers=="popins" & my_callers$win_size==win_size[i]),]$precision = sens_popins$precision
  my_callers[(my_callers$callers=="popins" & my_callers$win_size==win_size[i]),]$f.score = sens_popins$f1_score
  my_callers[(my_callers$callers=="popins" & my_callers$win_size==win_size[i]),]$g.error = sens_popins$g_error
  
  
  # pindel
  
  sens_pindel = sensitivity_precision_geno_error(pindel,golden,"deletion","pindel")
  
  my_callers[(my_callers$callers=="pindel" & my_callers$win_size==win_size[i]),]$sensitivity = sens_pindel$sensitivity
  my_callers[(my_callers$callers=="pindel" & my_callers$win_size==win_size[i]),]$precision = sens_pindel$precision
  my_callers[(my_callers$callers=="pindel" & my_callers$win_size==win_size[i]),]$f.score = sens_pindel$f1_score
  my_callers[(my_callers$callers=="pindel" & my_callers$win_size==win_size[i]),]$g.error = sens_pindel$g_error
  
  # manta
  
  sens_manta = sensitivity_precision_geno_error(manta,golden,"deletion","manta")
  
  my_callers[(my_callers$callers=="manta" & my_callers$win_size==win_size[i]),]$sensitivity = sens_manta$sensitivity
  my_callers[(my_callers$callers=="manta" & my_callers$win_size==win_size[i]),]$precision = sens_manta$precision
  my_callers[(my_callers$callers=="manta" & my_callers$win_size==win_size[i]),]$f.score = sens_manta$f1_score
  my_callers[(my_callers$callers=="manta" & my_callers$win_size==win_size[i]),]$g.error = sens_manta$g_error
  
  # svaba
  
  sens_svaba = sensitivity_precision_geno_error(svaba,golden,"deletion","svaba")
  
  my_callers[(my_callers$callers=="svaba" & my_callers$win_size==win_size[i]),]$sensitivity = sens_svaba$sensitivity
  my_callers[(my_callers$callers=="svaba" & my_callers$win_size==win_size[i]),]$precision = sens_svaba$precision
  my_callers[(my_callers$callers=="svaba" & my_callers$win_size==win_size[i]),]$f.score = sens_svaba$f1_score
  my_callers[(my_callers$callers=="svaba" & my_callers$win_size==win_size[i]),]$g.error = sens_svaba$g_error
  
  # whamg
  
  sens_whamg = sensitivity_precision_geno_error(whamg,golden,"deletion","whamg")
  
  my_callers[(my_callers$callers=="whamg" & my_callers$win_size==win_size[i]),]$sensitivity = sens_whamg$sensitivity
  my_callers[(my_callers$callers=="whamg" & my_callers$win_size==win_size[i]),]$precision = sens_whamg$precision
  my_callers[(my_callers$callers=="whamg" & my_callers$win_size==win_size[i]),]$f.score = sens_whamg$f1_score
  my_callers[(my_callers$callers=="whamg" & my_callers$win_size==win_size[i]),]$g.error = sens_whamg$g_error
  
  # delly
  
  sens_delly = sensitivity_precision_geno_error(delly,golden,"deletion","delly")
  
  my_callers[(my_callers$callers=="delly" & my_callers$win_size==win_size[i]),]$sensitivity = sens_delly$sensitivity
  my_callers[(my_callers$callers=="delly" & my_callers$win_size==win_size[i]),]$precision = sens_delly$precision
  my_callers[(my_callers$callers=="delly" & my_callers$win_size==win_size[i]),]$f.score = sens_delly$f1_score
  my_callers[(my_callers$callers=="delly" & my_callers$win_size==win_size[i]),]$g.error = sens_delly$g_error
  
  # gatk
  
  sens_gatk = sensitivity_precision_geno_error(gatk,golden,"deletion","gatk")
  
  my_callers[(my_callers$callers=="gatk" & my_callers$win_size==win_size[i]),]$sensitivity = sens_gatk$sensitivity
  my_callers[(my_callers$callers=="gatk" & my_callers$win_size==win_size[i]),]$precision = sens_gatk$precision
  my_callers[(my_callers$callers=="gatk" & my_callers$win_size==win_size[i]),]$f.score = sens_gatk$f1_score
  my_callers[(my_callers$callers=="gatk" & my_callers$win_size==win_size[i]),]$g.error = sens_gatk$g_error
  
  # strelka
  
  sens_strelka = sensitivity_precision_geno_error(strelka,golden,"deletion","strelka")
  
  my_callers[(my_callers$callers=="strelka" & my_callers$win_size==win_size[i]),]$sensitivity = sens_strelka$sensitivity
  my_callers[(my_callers$callers=="strelka" & my_callers$win_size==win_size[i]),]$precision = sens_strelka$precision
  my_callers[(my_callers$callers=="strelka" & my_callers$win_size==win_size[i]),]$f.score = sens_strelka$f1_score
  my_callers[(my_callers$callers=="strelka" & my_callers$win_size==win_size[i]),]$g.error = sens_strelka$g_error
  
  
}

dir.create("Insertions/outputs/")
write.csv(my_callers,"Insertions/outputs/window_size_callers.csv",row.names = F)


## plot results ####

my_callers = read.csv("Insertions/outputs/window_size_callers.csv")

my_callers$f.score = my_callers$f.score*100 
my_callers$g.error = 100 - my_callers$g.error

my_callers$sensitivity = round(my_callers$sensitivity,2)
my_callers$precision = round(my_callers$precision,2)
my_callers$f.score = round(my_callers$f.score,2)
my_callers$g.error = round(my_callers$g.error,2)

my_callers %>% arrange(callers,win_size) %>%
  write.csv("METRICS/window_size_insertions.csv",row.names=F)


delly = my_callers %>% filter(callers=="delly")
popins = my_callers %>% filter(callers=="popins")
gatk = my_callers %>% filter(callers=="gatk")
strelka = my_callers %>% filter(callers=="strelka")
whamg = my_callers %>% filter(callers=="whamg")
svaba = my_callers %>% filter(callers=="svaba")
manta = my_callers %>% filter(callers=="manta")
pindel = my_callers %>% filter(callers=="pindel")


png("Deleciones_150_tobig/outputs/window_size_callers.png",height =12,width = 12,res = 300,units="in")
par(mfrow=c(2,2))
plot(delly$win_size,delly$sensitivity,col=1,type="b",xlab="Window size (bp)", ylab = "Sensitivity (%)",
     ylim=c(min(c(my_callers$sensitivity,my_callers$precision,my_callers$f.score,my_callers$g.error)),
            max(c(my_callers$sensitivity,my_callers$precision,my_callers$f.score,my_callers$g.error))),lwd=2)
points(lumpy$win_size,lumpy$sensitivity,col=2,type="b",lwd=2)
points(genomestrip$win_size,genomestrip$sensitivity,col=3,type="b",lwd=2)
points(cnvnator$win_size,cnvnator$sensitivity,col=4,type="b",lwd=2)
points(whamg$win_size,whamg$sensitivity,col=5,type="b",lwd=2)
points(svaba$win_size,svaba$sensitivity,col=6,type="b",lwd=2)
points(manta$win_size,manta$sensitivity,col="darkgoldenrod2",type="b",lwd=2)
points(pindel$win_size,pindel$sensitivity,col=8,type="b",lwd=2)

legend("bottomright",c("delly","lumpy","genomestrip","cnvnator","whamg","svaba","manta",
                       "pindel"),col=c(1:6,"darkgoldenrod2",8),pch=16)


plot(delly$win_size,delly$precision,col=1,type="b",xlab="Window size (bp)", ylab = "Precision (%)",
     ylim=c(min(c(my_callers$sensitivity,my_callers$precision,my_callers$f.score,my_callers$g.error)),
            max(c(my_callers$sensitivity,my_callers$precision,my_callers$f.score,my_callers$g.error))))
points(lumpy$win_size,lumpy$precision,col=2,type="b",lwd=2)
points(genomestrip$win_size,genomestrip$precision,col=3,type="b",lwd=2)
points(cnvnator$win_size,cnvnator$precision,col=4,type="b",lwd=2)
points(whamg$win_size,whamg$precision,col=5,type="b",lwd=2)
points(svaba$win_size,svaba$precision,col=6,type="b",lwd=2)
points(manta$win_size,manta$precision,col="darkgoldenrod2",type="b",lwd=2)
points(pindel$win_size,pindel$precision,col=8,type="b",lwd=2)

plot(delly$win_size,delly$f.score,col=1,type="b",lwd=2,xlab="Window size (bp)", ylab = "F1-Score (%)",
     ylim=c(min(c(my_callers$sensitivity,my_callers$precision,my_callers$f.score,my_callers$g.error)),
            max(c(my_callers$sensitivity,my_callers$precision,my_callers$f.score,my_callers$g.error))))
points(lumpy$win_size,lumpy$f.score,col=2,type="b",lwd=2)
points(genomestrip$win_size,genomestrip$f.score,col=3,type="b",lwd=2)
points(cnvnator$win_size,cnvnator$f.score,col=4,type="b",lwd=2)
points(whamg$win_size,whamg$f.score,col=5,type="b",lwd=2)
points(svaba$win_size,svaba$f.score,col=6,type="b",lwd=2)
points(manta$win_size,manta$f.score,col="darkgoldenrod2",type="b",lwd=2)
points(pindel$win_size,pindel$f.score,col=8,type="b",lwd=2)

plot(delly$win_size,delly$g.error,col=1,type="b",lwd=2,xlab="Window size (bp)", ylab = "Genotype concordance (%)",
     ylim=c(min(c(my_callers$sensitivity,my_callers$precision,my_callers$f.score,my_callers$g.error)),
            max(c(my_callers$sensitivity,my_callers$precision,my_callers$f.score,my_callers$g.error))))
points(lumpy$win_size,lumpy$g.error,col=2,type="b",lwd=2)
points(genomestrip$win_size,genomestrip$g.error,col=3,type="b",lwd=2)
points(cnvnator$win_size,cnvnator$g.error,col=4,type="b",lwd=2)
points(whamg$win_size,whamg$g.error,col=5,type="b",lwd=2)
points(svaba$win_size,svaba$g.error,col=6,type="b",lwd=2)
points(manta$win_size,manta$g.error,col="darkgoldenrod2",type="b",lwd=2)
points(pindel$win_size,pindel$g.error,col=8,type="b",lwd=2)

dev.off()




# read vcfs files ########

delly = fread("Insertions/data/Delly_inserciones_30_50")
delly$V4 = as.numeric(abs(delly$V4))
delly = delly %>% filter(V4>30)
colnames(delly) = c("chr","start_delly","end","length_delly","GT_delly")
delly$lower_delly = delly$start_delly-10
delly$upper_delly = delly$start_delly+10
delly$chr_pos_delly = do.call(paste0,list(delly$chr,"_",delly$start_delly))
summary(delly$length_delly)

# pindel no da genotipo ni length para SV

pindel = fread("Insertions/data/Pindel_inserciones_totales",fill = T)
pindel = pindel[1:60,c(1,2,5,4)]
pindel$V4 = as.numeric(nchar(pindel$V4))
colnames(pindel) = c("chr","start_pindel","GT_pindel","length_pindel")

pindel2 = fread("Insertions/data/Inserciones_largas_genotipadas_pindel_ivan")
colnames(pindel2) = c("chr","start_pindel","GT_pindel")
pindel2$length_pindel = 0

pindel = rbind(pindel,pindel2)

#pindel = pindel %>% filter(length_pindel>30 & length_pindel<151)
pindel$lower_pindel = pindel$start_pindel-10
pindel$upper_pindel = pindel$start_pindel+10
pindel$chr_pos_pindel = do.call(paste0,list(pindel$chr,"_",pindel$start_pindel))
summary(pindel$length_pindel)
pindel$chr = as.character(pindel$chr)
#pindel$GT_pindel[pindel$GT_pindel==""] = "0/1"


whamg = fread("Insertions/data/wham_inserciones_con_genotipo")
#whamg = whamg %>% filter(V3>30 & V3<151)
#whamg = whamg[,-2,with=F]
colnames(whamg) = c("chr","start_whamg","length_whamg","GT_whamg")
whamg$lower_whamg = whamg$start_whamg-10
whamg$upper_whamg = whamg$start_whamg+10
whamg$chr_pos_whamg = do.call(paste0,list(whamg$chr,"_",whamg$start_whamg))
#whamg$GT_whamg = "0/1"
summary(whamg$length_whamg)

whamg = whamg %>% arrange(chr,chr_pos_whamg)
head(whamg,40)
whamg$start_whamg = as.numeric(whamg$start_whamg)
whamg$diff = c(0,diff(whamg$start_whamg))
whamg = whamg %>% filter(diff<0 | diff>10)


#svaba no da length para SV
svaba = fread("Insertions/data/Svaba_insercionesmedianas_totalSV",fill=T,skip = 5289)
svaba$length_svaba = abs(nchar(svaba$V3)-nchar(svaba$V4))
colnames(svaba)[c(1,2,5)] = c("chr","start_svaba","GT_svaba")
svaba = svaba[,c(1,2,5,6)]
#svaba2 = fread("Insertions/data/Svaba_insercionesmedianas_totalSV",fill=T)
#colnames(svaba2) = c("chr","start_svaba","GT_svaba")
#svaba2$length_svaba = 0
#svaba = rbind(svaba,svaba2)
svaba$lower_svaba = svaba$start_svaba-10
svaba$upper_svaba = svaba$start_svaba+10
svaba$chr_pos_svaba = do.call(paste0,list(svaba$chr,"_",svaba$start_svaba))
summary(svaba$length_svaba)

svaba = svaba %>% arrange(chr,start_svaba)
head(svaba,40)
svaba$start_svaba = as.numeric(svaba$start_svaba)
svaba$diff = c(0,diff(svaba$start_svaba))
svaba = svaba %>% filter(diff<0 | diff>10)



gatk = fread("Insertions/data/Gatk_inserciones_30_150")
gatk$length_gatk = abs(nchar(gatk$V3)-nchar(gatk$V4))
#gatk = gatk %>% filter(length_gatk>30 & length_gatk<151)
colnames(gatk)[c(1,2,5)] = c("chr","start_gatk","GT_gatk")
gatk = gatk[,c(1,2,5,6)]
gatk$end = gatk$length_gatk+gatk$start_gatk
gatk$lower_gatk = gatk$start_gatk-10
gatk$upper_gatk = gatk$start_gatk+10
gatk$chr_pos_gatk = do.call(paste0,list(gatk$chr,"_",gatk$start_gatk))
summary(gatk$length_gatk)

strelka = fread("Insertions/data/Strelka_inserciones_30_50")
strelka$length_strelka = abs(nchar(strelka$V3)-nchar(strelka$V4))
strelka = strelka %>% filter(length_strelka>30 & length_strelka<151)
colnames(strelka)[c(1,2,5)] = c("chr","start_strelka","GT_strelka")
strelka = strelka[,c(1,2,5,6)]
strelka$end = strelka$length_strelka+strelka$start_strelka
strelka$lower_strelka = strelka$start_strelka-10
strelka$upper_strelka = strelka$start_strelka+10
strelka$chr_pos_strelka = do.call(paste0,list(strelka$chr,"_",strelka$start_strelka))
summary(strelka$length_strelka)


#manta no da length para SV y no da genotipo para pequeñas

manta = fread("Insertions/data/Manta_inserciones_totales")
manta$V4 = as.numeric(abs(manta$V4))
#manta = manta %>% filter(V4>30 & V4<151)
manta$V3 = as.character(manta$V3)
manta$V3 = "0/1"
#manta$V3 = "9/9" #only for g.error
colnames(manta) = c("chr","start_manta","GT_manta","length_manta")
manta2 = fread("Insertions/data/Manta_inserciones_totales",skip=147)
head(manta2)
colnames(manta2) = c("chr","start_manta","GT_manta")
manta2$length_manta = 0
manta = rbind(manta,manta2)
manta$lower_manta = manta$start_manta-10
manta$upper_manta = manta$start_manta+10
manta$chr_pos_manta = do.call(paste0,list(manta$chr,"_",manta$start_manta))

manta = manta %>% arrange(chr,start_manta)
head(manta,40)
manta$start_manta = as.numeric(manta$start_manta)
manta$diff = c(0,diff(manta$start_manta))
manta = manta %>% filter(diff<0 | diff>10)


#popins no da length
popins = fread("Insertions/data/Popins_inserciones_totales")
colnames(popins) = c("chr","start_popins","GT_popins")
popins$lower_popins = popins$start_popins-200
popins$upper_popins = popins$start_popins+200
popins$chr_pos_popins = do.call(paste0,list(popins$chr,"_",popins$start_popins))
popins$length_popins = 0
summary(popins$length_popins)

popins = popins %>% arrange(chr,start_popins)
head(popins,40)
popins$start_popins = as.numeric(popins$start_popins)
popins$diff = c(0,diff(popins$start_popins))
popins = popins %>% filter(diff<0 | diff>200)



### sensitivity-specifitiy-geno error for each caller ######

# gatk 

sens_gatk = sensitivity_precision_geno_error(gatk,golden,"insertion","gatk")

# strelka

strelka = strelka %>% as.data.table()

sens_strelka = sensitivity_precision_geno_error(strelka,golden,"insertion","strelka")

# popins

popins = popins %>% as.data.table()

sens_popins = sensitivity_precision_geno_error(popins,golden,"insertion","popins")

# pindel

pindel = pindel %>% as.data.table()

sens_pindel = sensitivity_precision_geno_error(pindel,golden,"insertion","pindel")

# manta

manta = manta %>% as.data.table()

sens_manta = sensitivity_precision_geno_error(manta,golden,"insertion","manta")

# whamg

whamg = whamg %>% as.data.table()

sens_whamg = sensitivity_precision_geno_error(whamg,golden,"insertion","whamg")

# delly

delly = delly %>% as.data.table()

sens_delly = sensitivity_precision_geno_error(delly,golden,"insertion","delly")

# svaba

svaba = svaba %>% as.data.table()

sens_svaba = sensitivity_precision_geno_error(svaba,golden,"insertion","svaba")


metrics_callers = matrix(c(t(unlist(sens_delly))[c(3,4,1,2,5,6)],
                           t(unlist(sens_gatk))[c(3,4,1,2,5,6)],
                           t(unlist(sens_popins))[c(3,4,1,2,5,6)],
                           t(unlist(sens_manta))[c(3,4,1,2,5,6)],
                           t(unlist(sens_pindel))[c(3,4,1,2,5,6)],
                           t(unlist(sens_strelka))[c(3,4,1,2,5,6)],
                           t(unlist(sens_svaba))[c(3,4,1,2,5,6)],
                           t(unlist(sens_whamg))[c(3,4,1,2,5,6)]),byrow=T,ncol = 6)

metrics_callers = cbind(c("Delly","Gatk","Popins","Manta","Pindel","Strelka","Svaba","Whamg"),
                        round(metrics_callers,2))

colnames(metrics_callers) = c("Caller","TP","FP","Sensitivity","Precision","F1-Score","Genotype error")

metrics_callers = metrics_callers %>% as.data.frame() %>% arrange(desc(`F1-Score`))

metrics_callers = cbind(N=c(nrow(golden),rep("",nrow(metrics_callers)-1)),
                        metrics_callers)

metrics_callers = cbind(SV_type=c("Insertions",rep("",nrow(metrics_callers)-1)),
                        metrics_callers)

metrics_callers[,1:3] = apply(metrics_callers[,1:3], 2, function(x) as.character(x));
metrics_callers[,4:9] = apply(metrics_callers[,4:9], 2, function(x) as.numeric(x));

write.csv(metrics_callers,"METRICS/insertions.csv",row.names = F)


## prepare BBDD ####

popins$ID = paste0("insertion_",1:nrow(popins))
pindel$ID = "none"
manta$ID = "none"
whamg$ID = "none"
gatk$ID = "none"
svaba$ID = "none"
strelka$ID = "none"
delly$ID = "none"
golden$ID = "none"


my_data = merge_callers(pindel,popins,callers_ref = c("popins"),
                        caller_to_merge = c("pindel")) 

my_data = merge_callers(manta,my_data,callers_ref = c("popins","pindel"),
                        caller_to_merge = c("manta")) 

my_data = merge_callers(whamg,my_data,callers_ref = c("popins","pindel","manta"),
                        caller_to_merge = c("whamg")) 

my_data = merge_callers(gatk,my_data,callers_ref = c("popins","pindel","manta","whamg"),
                        caller_to_merge = c("gatk")) 

my_data = merge_callers(svaba,my_data,callers_ref = c("popins","pindel","manta","whamg",
                                                      "gatk"),
                        caller_to_merge = c("svaba")) 

my_data = merge_callers(strelka,my_data,callers_ref = c("popins","pindel","manta","whamg",
                                                      "gatk","svaba"),
                        caller_to_merge = c("strelka")) 

my_data = merge_callers(delly,my_data,callers_ref = c("popins","pindel","manta","whamg",
                                                      "gatk","svaba","strelka"),
                        caller_to_merge = c("delly")) 


my_data2 = merge_callers(golden,my_data,callers_ref = c("popins","pindel","manta","whamg",
                                                        "gatk","svaba","strelka","delly"),
                         caller_to_merge = c("golden"))


dim(my_data)
dim(my_data2)

my_data2 = unique(my_data2)

#saveRDS(my_data2,"Insertions/outputs/merge_callers.rds")

my_data2 = readRDS("Insertions/outputs/merge_callers.rds")


## add number of callers detected #####

my_data2$length_svaba = 0

my_data2$callers_detected = n_callers_detected(my_data2,
                                               c("popins","pindel","manta","whamg",
                                                 "gatk","svaba","strelka","delly"))

table(my_data2$callers_detected)

## add strategy ####

strategies = fread("strategies.csv")

my_data2$strategy = strategy(my_data2,
                             c("popins","pindel","manta","whamg",
                               "gatk","svaba","strelka","delly"),strategies)

table(my_data2$strategy)

### colineality #####

table(my_data2$callers_detected,my_data2$strategy)
chisq.test(table(my_data2$callers_detected,my_data2$strategy))

## prepare model ####

my_data2$PASS = "YES"

my_data2$PASS[is.na(my_data2$length_golden)] = "NO"

dim(golden)

my_data2 = my_data2[-which(duplicated(my_data2$ID)),]

my_data2 = my_data2 %>% filter(callers_detected!=0)

#my_data2 = my_data2 %>% filter(callers_detected!=1)


table(my_data2$PASS)

my_data2$GT_gatk[is.na(my_data2$GT_gatk)] = "9/9"
my_data2$GT_strelka[is.na(my_data2$GT_strelka)] = "9/9"
my_data2$GT_manta[is.na(my_data2$GT_manta)] = "9/9"
my_data2$GT_whamg[is.na(my_data2$GT_whamg)] = "9/9"
my_data2$GT_delly[is.na(my_data2$GT_delly)] = "9/9"
my_data2$GT_popins[is.na(my_data2$GT_popins)] = "9/9"
my_data2$GT_pindel[is.na(my_data2$GT_pindel)] = "9/9"
my_data2$GT_svaba[is.na(my_data2$GT_svaba)] = "9/9"

my_data2$GT_whamg[my_data2$GT_whamg==""] = "9/9"
my_data2$GT_popins[my_data2$GT_popins=="0/0"] = "9/9"

table(my_data2$GT_gatk,useNA = "always")
table(my_data2$GT_strelka,useNA = "always")
table(my_data2$GT_manta,useNA = "always")
table(my_data2$GT_whamg,useNA = "always")
table(my_data2$GT_delly,useNA = "always")
table(my_data2$GT_popins,useNA = "always")
table(my_data2$GT_svaba,useNA = "always")

table(my_data2$PASS)

my_data2$PASS = ifelse(my_data2$PASS=="YES",0,1)

my_model = glm(PASS ~ factor(GT_gatk) + 
                 factor(GT_strelka)+
                 factor(GT_manta)+
                 factor(GT_whamg)+
                 factor(GT_delly)+
                 factor(GT_popins)+
                 factor(GT_pindel)+
                 factor(GT_svaba) + 
                 factor(callers_detected), 
                 #factor(strategy)+
               data=my_data2,family = "binomial")

my_data2$PASS = ifelse(my_data2$PASS==0,"YES","NO")

summary(my_model)

my_pred = predict(my_model,my_data2,type="response")

my_pred = ifelse(my_pred>0.5,1,0)

my_pred = ifelse(my_pred==1,"NO PASS","PASS")

table(my_data2$PASS)
table(my_pred,my_data2$PASS)
table(my_pred)

length(unique(golden$chr_pos))

# without CHR and BP
#sensitivity
model_sens = (table(my_pred,my_data2$PASS)[2,2]/nrow(golden))*100
#specificity
model_spec = (table(my_pred,my_data2$PASS)[2,2])/(table(my_pred,my_data2$PASS)[2,2]+
                                                    table(my_pred,my_data2$PASS)[2,1])*100

f_score = ((2*model_sens*model_spec)/(model_sens + model_spec))/100


table(my_pred,my_data2$callers_detected)


### choose best model ####

my_data2$PASS = ifelse(my_data2$PASS=="YES",0,1)

full.model =  glm(PASS ~ factor(GT_gatk) + 
                    factor(GT_strelka)+
                    factor(GT_manta)+
                    factor(GT_whamg)+
                    factor(GT_delly)+
                    factor(GT_popins)+
                    factor(GT_pindel)+
                    factor(GT_svaba) + 
                    factor(callers_detected)+ factor(strategy),
                  data=my_data2,family = "binomial")

base.model = glm(PASS ~ factor(GT_gatk) + 
                   factor(GT_strelka)+
                   factor(GT_manta)+
                   factor(GT_whamg)+
                   factor(GT_delly)+
                   factor(GT_popins)+
                   factor(GT_pindel)+
                   factor(GT_svaba),
                 data=my_data2,family = "binomial")

step(full.model,
     scope = list(lower = formula(base.model),
                  upper = formula(full.model)),
     direction = "backward")

best.model = glm(formula = PASS ~ factor(GT_gatk) + factor(GT_strelka) + factor(GT_manta) + 
                   factor(GT_whamg) + factor(GT_delly) + factor(GT_popins) + 
                   factor(GT_pindel) + factor(GT_svaba) + factor(callers_detected), 
                 family = "binomial", data = my_data2)

my_data2$PASS = ifelse(my_data2$PASS==0,"YES","NO")

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
model_sens = round(model_tp/N*100,3)
#specificity
model_spec = round(model_tp/(model_fp+model_tp)*100,3)

model_tp;model_fp;model_sens;model_spec

f_score = round(((2*model_sens*model_spec)/(model_sens + model_spec))/100,3)

my_data2$prediction = my_pred


# que lo detecte al menos 1 caller

det_call = my_data2 %>% filter(callers_detected %in% 1:9)

callers1 = logical_sens_spec(det_call)

# que lo detecte al menos 2 caller y 2 estrategias a la vez

det_call = my_data2 %>% filter(callers_detected %in% 2:9 & strategy %in% 2:6)

callers2 = logical_sens_spec(det_call)

# que lo detecte al menos 3 caller

det_call = my_data2 %>% filter(callers_detected %in% 3:9 & strategy %in% 2:6)

callers3 = logical_sens_spec(det_call)

# que lo detecte al menos 4 caller

det_call = my_data2 %>% filter(callers_detected %in% 4:9 & strategy %in% 2:6)

callers4 = logical_sens_spec(det_call)

metrics_callers = rbind(metrics_callers,
                        c("","","LR",model_tp,model_fp,model_sens,model_spec,f_score,
                          ""))

metrics_callers = rbind(metrics_callers,
                        c("","",">=1 caller",callers1,""))

metrics_callers = rbind(metrics_callers,
                        c("","",">=2 caller **",callers2,""))

metrics_callers = rbind(metrics_callers,
                        c("","",">=3 caller **",callers3,""))

metrics_callers = rbind(metrics_callers,
                        c("","",">=4 caller **",callers4,""))


write.csv(metrics_callers,"METRICS/insertions.csv",row.names = F)


## update Manta, Pindel, Whamg genotypes #####

my_data2 = my_data2 %>% select(-GT_manta,-GT_pindel,-GT_whamg)

my_data2 = left_join(my_data2,manta %>% select(chr_pos_manta,GT_manta))
my_data2 = left_join(my_data2,pindel %>% select(chr_pos_pindel,GT_pindel))
my_data2 = left_join(my_data2,whamg %>% select(chr_pos_whamg,GT_whamg))

## genotype error in the model ######

data_call = my_data2 %>% filter(prediction=="PASS") %>% filter(PASS=="YES") %>% 
  dplyr::select(GT_whamg,
                GT_strelka,
                GT_pindel,
                GT_delly,
                GT_manta,
                GT_popins,
                GT_gatk,
                GT_svaba,
                GT_golden,PASS,prediction)

data_call$GT_delly[data_call$GT_delly=="9/9"] = NA
data_call$GT_pindel[data_call$GT_pindel=="9/9"] = NA
data_call$GT_whamg[data_call$GT_whamg=="9/9"] = NA
data_call$GT_strelka[data_call$GT_strelka=="9/9"] = NA
data_call$GT_manta[data_call$GT_manta=="9/9"] = NA
data_call$GT_svaba[data_call$GT_svaba=="9/9"] = NA
data_call$GT_gatk[data_call$GT_gatk=="9/9"] = NA
data_call$GT_popins[data_call$GT_popins=="9/9"] = NA

table(data_call$GT_whamg)
table(data_call$GT_strelka)
table(data_call$GT_pindel)
table(data_call$GT_delly)
table(data_call$GT_manta)
table(data_call$GT_svaba)
table(data_call$GT_gatk)
table(data_call$GT_popins)


data_call$GT = NA

data_call = as.matrix(data_call)

for(i in 1:nrow(data_call)){
  
  my_tab = table(as.character(data_call[i,c(1:8)]))
  
  if(length(which(my_tab==max(my_tab)))==1){
    data_call[i,12] = names(which.max(my_tab))
  }
  
  if(length(which(my_tab==max(my_tab)))>1){
    
    data_call[i,12] = NA
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



# caret package cv k-fold ####

metrics_callers = read.csv("METRICS/insertions.csv")
metrics_callers$Caller = as.character(metrics_callers$Caller)

library(caret)
library(e1071)

ctrl <- trainControl(method = "repeatedcv", number = 10, 
                     savePredictions = TRUE)

## 70-30

0.7*nrow(my_data2)

set.seed(2589)

training = sample(1:nrow(my_data2),1223)

table(my_data2$PASS,useNA = "always")

my_data2$PASS = as.factor(my_data2$PASS)

my_data2_70 = my_data2[training,]

table(my_data2_70$PASS,useNA = "always")

mod_fit <- train(factor(PASS) ~ factor(GT_gatk) + factor(GT_strelka) + factor(GT_manta) + 
                   factor(GT_whamg) + factor(GT_delly) + factor(GT_popins) + 
                   factor(GT_pindel) + factor(GT_svaba) + factor(callers_detected), 
                 data=my_data2_70, method="glm", family="binomial",
                 trControl = ctrl)

my_pred = predict(mod_fit, newdata=my_data2_70)

my_pred = ifelse(my_pred=="YES",0,1)

my_pred = ifelse(my_pred==0,"PASS","NO PASS")

table(my_data2_70$PASS)
table(my_pred,my_data2_70$PASS)
table(my_pred)

N = as.numeric(table(my_data2_70$PASS)[2])

tp_cv = as.numeric(table(my_pred,my_data2_70$PASS)[2,2])
fp_cv = as.numeric(table(my_pred)[2])-tp_cv

sens_call_cv = round(tp_cv/N*100,3)
spec_call_cv = round(tp_cv/(fp_cv+tp_cv)*100,3)

f_score_cv = round((2*sens_call_cv*spec_call_cv)/(sens_call_cv+spec_call_cv)/100,3)

saveRDS(mod_fit,"Inversiones/outputs/model.rds")

mod_fit = readRDS("Inversiones/outputs/model.rds")


# que lo detecte al menos 1 caller

det_call = my_data2_70 %>% filter(callers_detected %in% 1:9)

callers1 = logical_sens_spec(det_call)

# que lo detecte al menos 2 caller y 2 estrategias a la vez

det_call = my_data2_70 %>% filter(callers_detected %in% 2:9 & strategy %in% 2:6)

callers2 = logical_sens_spec(det_call)

# que lo detecte al menos 3 caller

det_call = my_data2_70 %>% filter(callers_detected %in% 3:9 & strategy %in% 2:6)

callers3 = logical_sens_spec(det_call)

# que lo detecte al menos 4 caller

det_call = my_data2_70 %>% filter(callers_detected %in% 4:9 & strategy %in% 2:6)

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


my_data2_30 = my_data2[-training,]

my_pred = predict(mod_fit,my_data2_30)

my_pred = ifelse(my_pred=="YES",0,1)

my_pred = ifelse(my_pred==0,"PASS","NO PASS")

table(my_data2_30$PASS)
table(my_pred,my_data2_30$PASS)
table(my_pred)

N = as.numeric(table(my_data2_30$PASS)[2])

tp_cv = as.numeric(table(my_pred,my_data2_30$PASS)[2,2])
fp_cv = as.numeric(table(my_pred)[2])-tp_cv

sens_call_cv = round(tp_cv/N*100,3)
spec_call_cv = round(tp_cv/(fp_cv+tp_cv)*100,3)

f_score_cv = round((2*sens_call_cv*spec_call_cv)/(sens_call_cv+spec_call_cv)/100,3)


# que lo detecte al menos 1 caller

det_call = my_data2_30 %>% filter(callers_detected %in% 1:9)

callers1 = logical_sens_spec(det_call)

# que lo detecte al menos 2 caller y 2 estrategias a la vez

det_call = my_data2_30 %>% filter(callers_detected %in% 2:9 & strategy %in% 2:6)

callers2 = logical_sens_spec(det_call)

# que lo detecte al menos 3 caller

det_call = my_data2_30 %>% filter(callers_detected %in% 3:9 & strategy %in% 2:6)

callers3 = logical_sens_spec(det_call)

# que lo detecte al menos 4 caller

det_call = my_data2_30 %>% filter(callers_detected %in% 4:9 & strategy %in% 2:6)

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

metrics_callers$Sensitivity = as.numeric(metrics_callers$Sensitivity)
metrics_callers$Precision = as.numeric(metrics_callers$Precision)

metrics_callers$F1.Score = ((2*metrics_callers$Sensitivity*metrics_callers$Precision)/
                              (metrics_callers$Sensitivity + metrics_callers$Precision))

head(metrics_callers)
library(xtable)
print(xtable(metrics_callers[,-1],digits = 2),include.rownames = F)



## Random forest #####

x_train = cbind(PASS = my_data2_70$PASS,as.matrix(my_data2_70[,which(colnames(my_data2_70) %in% 
                       c("GT_popins","GT_pindel","GT_manta","GT_whamg",
                         "GT_svaba","GT_delly","GT_strelka","GT_gatk","strategy")),with=F ]))

x_train = as.data.frame(x_train)
head(x_train)
str(x_train)
rf <- randomForest(PASS ~ ., x_train)
rf

table(my_data2_70$PASS)
table(rf$predicted)

337/402
337/378

x_train_30 = cbind(PASS = my_data2_30$PASS,as.matrix(my_data2_30[,which(colnames(my_data2_30) %in% 
                                                                       c("GT_popins","GT_pindel","GT_manta","GT_whamg",
                                                                         "GT_svaba","GT_delly","GT_strelka","GT_gatk","strategy")),with=F ]))

my_pred = predict(rf,x_train_30)
table(my_pred)
table(my_data2_30$PASS)
table(my_pred,my_data2_30$PASS)

153/184
153/161

rf$confusion



#saveRDS(my_model,"Indels_30_150/outputs/model_insertion_30_150.rds")

#my_model = readRDS("Indels_30_150/outputs/model_insertion_30_150.rds")



# que lo detecte al menos 1 caller

det_call = my_data2 %>% filter(callers_detected %in% 1:9)

sens_call1 = table(det_call$PASS)[2]/nrow(golden)*100
spec_call1 = table(det_call$PASS)[2]/nrow(det_call)*100

f_score1 = (2*sens_call1*spec_call1)/(sens_call1+spec_call1)/100


# que lo detecte al menos 2 caller y 2 estrategias a la vez

det_call = my_data2 %>% filter(callers_detected %in% 2:9 & strategy %in% 2:6 & reciprocity >80)

sens_call2 = table(det_call$PASS)[2]/nrow(golden)*100
spec_call2 = table(det_call$PASS)[2]/nrow(det_call)*100

f_score2 = (2*sens_call2*spec_call2)/(sens_call2+spec_call2)/100


# que lo detecte al menos 3 caller

det_call = my_data2 %>% filter(callers_detected %in% 3:9 & strategy %in% 2:6 & reciprocity >80)

sens_call3 = table(det_call$PASS)[2]/nrow(golden)*100
spec_call3 = table(det_call$PASS)[2]/nrow(det_call)*100

f_score3 = (2*sens_call3*spec_call3)/(sens_call3+spec_call3)/100


# que lo detecte al menos 4 caller

det_call = my_data2 %>% filter(callers_detected %in% 4:9 & strategy %in% 2:6 & reciprocity >80)

sens_call4 = table(det_call$PASS)[2]/nrow(golden)*100
spec_call4 = table(det_call$PASS)[2]/nrow(det_call)*100

f_score4 = (2*sens_call4*spec_call4)/(sens_call4+spec_call4)/100


my_tab_metrics = matrix(c(sens_call1,spec_call1,f_score1,
                          sens_call2,spec_call2,f_score2,
                          sens_call3,spec_call3,f_score3,
                          sens_call4,spec_call4,f_score4),byrow=T,ncol=3)

print(xtable(my_tab_metrics,digits = 3),include.rownames = F)




### genotype error in the model ######

my_data2$GT = NA

data_call = my_data2 %>% filter(prediction=="PASS") %>% filter(PASS=="YES") %>% 
  dplyr::select(GT_gatk,
                GT_whamg,
                GT_svaba,
                GT_pamir,
                GT_pindel,
                GT_delly,
                GT_strelka,
                GT_deepvariant,
                GT_manta,
                GT,GT_golden,PASS,prediction)



data_call$GT_deepvariant[data_call$GT_deepvariant=="9/9"] = NA
data_call$GT_gatk[data_call$GT_gatk=="9/9"] = NA
data_call$GT_strelka[data_call$GT_strelka=="9/9"] = NA
data_call$GT_delly[data_call$GT_delly=="9/9"] = NA
data_call$GT_pindel[data_call$GT_pindel=="9/9"] = NA
data_call$GT_svaba[data_call$GT_svaba=="9/9"] = NA
data_call$GT_whamg[data_call$GT_whamg=="9/9"] = NA
data_call$GT_pamir[data_call$GT_pamir=="9/9"] = NA
data_call$GT_manta[data_call$GT_manta=="9/9"] = NA


table(data_call$GT_gatk)
table(data_call$GT_whamg)
table(data_call$GT_svaba)
table(data_call$GT_pamir)
table(data_call$GT_pindel)
table(data_call$GT_delly)
table(data_call$GT_strelka)
table(data_call$GT_deepvariant)
table(data_call$GT_manta)


data_call = as.matrix(data_call)

for(i in 1:nrow(data_call)){
  
  my_tab = table(as.character(data_call[i,c(1:9)]))
  
  if(length(which(my_tab==max(my_tab)))==1){
    data_call[i,10] = names(which.max(my_tab))
  }
  
  if(length(which(my_tab==max(my_tab)))>1){
    
    data_call[i,10] = NA
  }
}


head(data_call)
dim(data_call)


# use genotype from GATK because SVABA is bad

data_call[is.na(data_call[,10]),10] = data_call[is.na(data_call[,10]),1] 


data_call = as.data.frame(data_call)

table(data_call$GT,useNA = "always")

# concordance genotype of the model with golden ######

model_g_error = 100 -(sum(as.character(data_call$GT)==as.character(data_call$GT_golden))/
                        nrow(data_call))*100



#### Plot graphics ##########

library(data.table)
library(ggplot2)
library(plyr)

# read golden ######

golden_ins = fread("Insertions/data/Insilico3_inserciones_grandes_virus",fill=T)

golden_ins = golden_ins[,c(1,3,4),with=F]

golden_ins$V4 = nchar(golden_ins$V4)

summary(golden_ins$V4)

golden_ins$V1 = gsub("chr","",golden_ins$V1)

golden_ins2 = fread("Insertions/data/Insilico3_inserciones_grandes_virus",fill=T,skip = 180)

golden_ins2 = golden_ins2[,c(1,3,4),with=F]

golden_ins2$V4 = nchar(golden_ins2$V4)

summary(golden_ins2$V4)

golden_ins2$V1 = gsub("chr","",golden_ins2$V1)

golden_ins = rbind(golden_ins,golden_ins2)

golden_geno_ins = NULL

for(i in names(table(golden_ins$V1))){
  
  aux = golden_ins[golden_ins$V1==i,]
  
  which_duplicated = aux$V3[which(duplicated(aux$V3))]
  
  golden_dup = unique(aux[aux$V3 %in% which_duplicated,])
  golden_no_dup = aux[!aux$V3 %in% which_duplicated,]
  
  golden_dup$genotype = "1/1"
  golden_no_dup$genotype = "0/1"
  
  aux_geno = rbind(golden_dup,golden_no_dup)
  
  golden_geno_ins = rbind(golden_geno_ins,aux_geno)
}

golden = golden_geno_ins
colnames(golden) = c("chr","start_golden","length_golden","GT_golden")
golden$lower_golden = golden$start_golden-0
golden$upper_golden = golden$start_golden+0
golden$chr_pos_golden = do.call(paste0,list(golden$chr,"_",golden$start_golden))


aux = cut(golden$length_golden,
          breaks = c(30,75,150,300,500,1000,2000,3000,Inf),right = F)
sum(table(aux))
golden$length_stretch = aux

golden_length = golden %>% group_by(length_stretch) %>% 
  dplyr::summarise(golden = n()) %>% arrange(length_stretch) %>% as.data.frame()


### model #####

my_data2 = readRDS("Insertions/outputs/merge_callers_predicted.rds")

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


sens_precision = sensitivity_precision_caller(my_data2,golden_length,caller_name ="gatk")
sens_precision_final = rbind(sens_precision_final,sens_precision)

sens_precision = sensitivity_precision_caller(my_data2,golden_length,caller_name ="strelka")
sens_precision_final = rbind(sens_precision_final,sens_precision)

sens_precision = sensitivity_precision_caller(my_data2,golden_length,caller_name ="delly")
sens_precision_final = rbind(sens_precision_final,sens_precision)

sens_precision = sensitivity_precision_caller(my_data2,golden_length,caller_name ="manta")
sens_precision_final = rbind(sens_precision_final,sens_precision)

sens_precision = sensitivity_precision_caller(my_data2,golden_length,caller_name ="svaba")
sens_precision_final = rbind(sens_precision_final,sens_precision)

sens_precision = sensitivity_precision_caller(my_data2,golden_length,caller_name ="popins")
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
                                       to = c("firebrick2","forestgreen","black",
                                              "darkorange2","darkorchid2","darkolivegreen2",
                                              "cadetblue3","deeppink2","darkgoldenrod2"))


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

# extra colors "blue3"

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
  ggtitle(c("119 Insertions"))+
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



pdf("Insertions/outputs/Insertions_large_size_svaba_150.pdf",width = 11,height = 5)
p1 +  geom_line(sens_precision_final,mapping = aes(x = breaks, y = sensitivity,  color = caller)) + 
  geom_point(sens_precision_final,mapping = aes(x = breaks, y = sensitivity,  color = caller)) +
  geom_hline(yintercept = 0,linetype="dashed") + ylab("Sensitivity            (%)            Precision") +
  geom_point( aes(x = breaks, y = sensitivity,  color = caller, group = caller),model_in,size=2) +
  geom_line( aes(x = breaks, y = sensitivity,  color = caller, group = caller),model_in,size=1.2)
dev.off()






############ without PAMIR ###################

library(data.table)
library(dplyr)

setwd("/run/media/igalvan/TOSHIBA EXT/BSC/Machine_vector/SV")


source("functions/functions.R")

# read golden ######

golden_ins = fread("Indels_30_150/Inserciones/Insilico3_inserciones_30_150")

golden_ins$V3 = nchar(golden_ins$V3)

golden_ins = golden_ins %>% filter(V3>30 & V3<151)

golden_ins$V1 = gsub("chr","",golden_ins$V1)

golden_geno_ins = NULL

for(i in names(table(golden_ins$V1))){
  
  aux = golden_ins[golden_ins$V1==i,]
  
  which_duplicated = aux$V2[which(duplicated(aux$V2))]
  
  golden_dup = unique(aux[aux$V2 %in% which_duplicated,])
  golden_no_dup = aux[!aux$V2 %in% which_duplicated,]
  
  golden_dup$genotype = "1/1"
  golden_no_dup$genotype = "0/1"
  
  aux_geno = rbind(golden_dup,golden_no_dup)
  
  golden_geno_ins = rbind(golden_geno_ins,aux_geno)
}

golden = golden_geno_ins
colnames(golden) = c("chr","start_golden","length_golden","GT_golden")
golden$lower_golden = golden$start_golden-0
golden$upper_golden = golden$start_golden+0
golden$chr_pos_golden = do.call(paste0,list(golden$chr,"_",golden$start_golden))


# read vcfs files ########

delly = fread("Indels_30_150/Inserciones/Delly_inserciones_30_50")
delly$V4 = as.numeric(abs(delly$V4))
delly = delly %>% filter(V4>30 & V4<151)
colnames(delly) = c("chr","start_delly","end","length_delly","GT_delly")
delly$lower_delly = delly$start_delly-10
delly$upper_delly = delly$start_delly+10
delly$chr_pos_delly = do.call(paste0,list(delly$chr,"_",delly$start_delly))
summary(delly$length_delly)

pindel = fread("Indels_30_150/Inserciones/Pindels_inserciones_30_50")
pindel$V4 = as.numeric(nchar(pindel$V4))
colnames(pindel) = c("chr","start_pindel","end","length_pindel","GT_pindel")
pindel = pindel %>% filter(length_pindel>30 & length_pindel<151)
pindel$lower_pindel = pindel$start_pindel-10
pindel$upper_pindel = pindel$start_pindel+10
pindel$chr_pos_pindel = do.call(paste0,list(pindel$chr,"_",pindel$start_pindel))
summary(pindel$length_pindel)
pindel$chr = as.character(pindel$chr)

whamg = fread("Indels_30_150/Inserciones/Whamg_inserciones_30_150",fill = T)
whamg = whamg %>% filter(V3>30 & V3<151)
colnames(whamg) = c("chr","start_whamg","length_whamg")
whamg$lower_whamg = whamg$start_whamg-10
whamg$upper_whamg = whamg$start_whamg+10
whamg$chr_pos_whamg = do.call(paste0,list(whamg$chr,"_",whamg$start_whamg))
whamg$GT_whamg = "0/1"
summary(whamg$length_whamg)

svaba = fread("Indels_30_150/Inserciones/Svaba_inserciones_30_150")
svaba$length_svaba = abs(nchar(svaba$V3)-nchar(svaba$V4))
svaba = svaba %>% filter(length_svaba>30 & length_svaba<151)
colnames(svaba)[c(1,2,5)] = c("chr","start_svaba","GT_svaba")
svaba = svaba[,c(1,2,5,6)]
svaba$lower_svaba = svaba$start_svaba-10
svaba$upper_svaba = svaba$start_svaba+10
svaba$chr_pos_svaba = do.call(paste0,list(svaba$chr,"_",svaba$start_svaba))
summary(svaba$length_svaba)


gatk = fread("Indels_30_150/Inserciones/Gatk_inserciones_30_150")
gatk$length_gatk = abs(nchar(gatk$V3)-nchar(gatk$V4))
gatk = gatk %>% filter(length_gatk>30 & length_gatk<151)
colnames(gatk)[c(1,2,5)] = c("chr","start_gatk","GT_gatk")
gatk = gatk[,c(1,2,5,6)]
gatk$end = gatk$length_gatk+gatk$start_gatk
gatk$lower_gatk = gatk$start_gatk-10
gatk$upper_gatk = gatk$start_gatk+10
gatk$chr_pos_gatk = do.call(paste0,list(gatk$chr,"_",gatk$start_gatk))
summary(gatk$length_gatk)


deepvariant = fread("Indels_30_150/Inserciones/Deepvariant_inserciones_30_50")
deepvariant$length_deepvariant = abs(nchar(deepvariant$V3)-nchar(deepvariant$V4))
deepvariant = deepvariant %>% filter(length_deepvariant>30 & length_deepvariant<151)
colnames(deepvariant)[c(1,2,5)] = c("chr","start_deepvariant","GT_deepvariant")
deepvariant = deepvariant[,c(1,2,5,6)]
deepvariant$end = deepvariant$length_deepvariant+deepvariant$start_deepvariant
deepvariant$lower_deepvariant = deepvariant$start_deepvariant-10
deepvariant$upper_deepvariant = deepvariant$start_deepvariant+10
deepvariant$chr_pos_deepvariant = do.call(paste0,list(deepvariant$chr,"_",deepvariant$start_deepvariant))
summary(deepvariant$length_deepvariant)


strelka = fread("Indels_30_150/Inserciones/Strelka_inserciones_30_50")
strelka$length_strelka = abs(nchar(strelka$V3)-nchar(strelka$V4))
strelka = strelka %>% filter(length_strelka>30 & length_strelka<151)
colnames(strelka)[c(1,2,5)] = c("chr","start_strelka","GT_strelka")
strelka = strelka[,c(1,2,5,6)]
strelka$end = strelka$length_strelka+strelka$start_strelka
strelka$lower_strelka = strelka$start_strelka-10
strelka$upper_strelka = strelka$start_strelka+10
strelka$chr_pos_strelka = do.call(paste0,list(strelka$chr,"_",strelka$start_strelka))
summary(strelka$length_strelka)


manta = fread("Indels_30_150/Inserciones/Manta_inserciones_30_50",fill = T)
manta$V4 = as.numeric(abs(manta$V4))
manta = manta %>% filter(V4>30 & V4<151)
manta$V5 = "0/1"
colnames(manta) = c("chr","start_manta","end","length_manta","GT_manta")
manta$lower_manta = manta$start_manta-10
manta$upper_manta = manta$start_manta+10
manta$chr_pos_manta = do.call(paste0,list(manta$chr,"_",manta$start_manta))


## prepare BBDD ####


svaba$ID = paste0("insertion_",1:nrow(svaba))
manta$ID = "none"
gatk$ID = "none"
strelka$ID = "none"
pindel$ID = "none"
delly$ID = "none"
whamg$ID = "none"
deepvariant$ID = "none"
golden$ID = "none"


my_data = merge_callers(manta,svaba,callers_ref = c("svaba"),
                        caller_to_merge = c("manta")) 

my_data = merge_callers(gatk,my_data,callers_ref = c("svaba","manta"),
                        caller_to_merge = c("gatk")) 

my_data = merge_callers(strelka,my_data,callers_ref = c("svaba","manta","gatk"),
                        caller_to_merge = c("strelka")) 

my_data = merge_callers(pindel,my_data,callers_ref = c("svaba","manta","gatk","strelka"),
                        caller_to_merge = c("pindel")) 

my_data = merge_callers(delly,my_data,callers_ref = c("svaba","manta","gatk","strelka",
                                                      "pindel"),
                        caller_to_merge = c("delly")) 

my_data = merge_callers(whamg,my_data,callers_ref = c("svaba","manta","gatk","strelka",
                                                      "pindel","delly"),
                        caller_to_merge = c("whamg")) 

my_data = merge_callers(deepvariant,my_data,callers_ref = c("svaba","manta","gatk","strelka",
                                                            "pindel","delly","whamg"),
                        caller_to_merge = c("deepvariant")) 


my_data2 = merge_callers(golden,my_data,callers_ref = c("svaba","manta","gatk","strelka",
                                                        "pindel","delly","whamg","deepvariant"),
                         caller_to_merge = c("golden"))


dim(my_data)
dim(my_data2)

my_data2 = unique(my_data2)

### consensus length #######

data_length = my_data2 %>% dplyr::select(length_gatk,
                                         length_whamg,
                                         length_svaba,
                                         length_pindel,
                                         length_delly,
                                         length_strelka,
                                         length_deepvariant,
                                         length_manta)

data_length$consensus = NA

data_length = as.matrix(data_length)

for(i in 1:nrow(data_length)){
  
  my_tab = table(data_length[i,1:8])
  
  if(length(which(my_tab==max(my_tab)))==1){
    data_length[i,9] = names(which.max(my_tab))
  }
  
  if(length(which(my_tab==max(my_tab)))>1){
    
    data_length[i,9] = names(my_tab)[1]
  }
}

data_length = as.data.frame(data_length)

table(data_length$consensus,useNA = "always")

my_data2$length = as.numeric(as.character.numeric_version(data_length$consensus))

my_data2 = my_data2 %>% filter(!is.na(length))

aux = cut(my_data2$length,
          breaks = c(30,50,75,125,150))

my_data2$length_stretch = aux


## add numero de callers detectados ##########

data_call = my_data2 %>% dplyr::select(length_gatk,
                                       length_whamg,
                                       length_svaba,
                                       length_pindel,
                                       length_delly,
                                       length_strelka,
                                       length_deepvariant,
                                       length_manta)

data_call$consensus = NA

data_call = as.matrix(data_call)

for(i in 1:nrow(data_call)){
  
  my_tab = sum(is.na(data_call[i,1:8]))
  
  
  data_call[i,9] = 8-my_tab
  
}

data_call = as.data.frame(data_call)

table(data_call$consensus,useNA = "always")

my_data2$callers_detected = as.numeric(as.character.numeric_version(data_call$consensus))


## add reciprocity ###

data_call = my_data2 %>% dplyr::select(length_gatk,
                                       length_whamg,
                                       length_svaba,
                                       length_pindel,
                                       length_delly,
                                       length_strelka,
                                       length_deepvariant,
                                       length_manta)

data_call$reciprocity = NA

data_call = as.matrix(data_call)

for(i in 1:nrow(data_call)){
  
  if(sum(is.na(data_call[i,1:8]))<7){
    
    end = max(data_call[i,1:8],na.rm = T)
    start = min(data_call[i,1:8],na.rm = T)
    
    data_call[i,9] = start/end*100
    
    
  }
  
  if(sum(is.na(data_call[i,1:8]))==7){
    
    data_call[i,9] = 100
    
    
  }
  
}

data_call = as.data.frame(data_call)

table(data_call$reciprocity,useNA = "always")

my_data2$reciprocity = as.numeric(as.character.numeric_version(data_call$reciprocity))



## add strategy ####

strategies = fread("strategies.csv")

data_call = my_data2 %>% dplyr::select(length_gatk,
                                       length_whamg,
                                       length_svaba,
                                       length_pindel,
                                       length_delly,
                                       length_strelka,
                                       length_deepvariant,
                                       length_manta)

data_call$strategy = NA

data_call = as.matrix(data_call)

for(i in 1:nrow(data_call)){
  
  which_na = !is.na(data_call[i,1:8])
  
  my_tab = names(data_call[i,which_na])
  
  my_tab = gsub("length_","",my_tab)
  
  my_strat = strategies[strategies$CALLER %in% my_tab,2]$STRATEGY
  
  data_call[i,9] = length(unique(my_strat))
  
  
}

data_call = as.data.frame(data_call)

table(data_call$strategy,useNA = "always")

my_data2$strategy = as.numeric(as.character.numeric_version(data_call$strategy))


## prepare model ####

my_data2$PASS = "YES"

my_data2$PASS[is.na(my_data2$length_golden)] = "NO"

table(my_data2$PASS)

dim(golden)

my_data2 = my_data2[-which(duplicated(my_data2$ID)),]

my_data2$GT_deepvariant[is.na(my_data2$GT_deepvariant)] = "9/9"
my_data2$GT_gatk[is.na(my_data2$GT_gatk)] = "9/9"
my_data2$GT_strelka[is.na(my_data2$GT_strelka)] = "9/9"
my_data2$GT_manta[is.na(my_data2$GT_manta)] = "9/9"
my_data2$GT_whamg[is.na(my_data2$GT_whamg)] = "9/9"
my_data2$GT_delly[is.na(my_data2$GT_delly)] = "9/9"
my_data2$GT_pindel[is.na(my_data2$GT_pindel)] = "9/9"
my_data2$GT_svaba[is.na(my_data2$GT_svaba)] = "9/9"

my_data2$GT_deepvariant[my_data2$GT_deepvariant=="0/0"] = "9/9"
my_data2$GT_whamg[my_data2$GT_whamg==""] = "9/9"

table(my_data2$GT_deepvariant,useNA = "always")
table(my_data2$GT_gatk,useNA = "always")
table(my_data2$GT_strelka,useNA = "always")
table(my_data2$GT_manta,useNA = "always")
table(my_data2$GT_whamg,useNA = "always")
table(my_data2$GT_delly,useNA = "always")
table(my_data2$GT_svaba,useNA = "always")

table(my_data2$PASS)

my_data2$PASS = ifelse(my_data2$PASS=="YES",0,1)

my_model = glm(PASS ~ factor(GT_deepvariant) + 
                 factor(GT_gatk) + 
                 factor(GT_strelka)+
                 factor(GT_manta)+
                 factor(GT_whamg)+
                 factor(GT_delly)+
                 factor(GT_pindel)+
                 factor(GT_svaba) + factor(length)+ 
                 factor(callers_detected) + 
                 #factor(strategy)+
                 reciprocity,
               data=my_data2,family = "binomial")

my_data2$PASS = ifelse(my_data2$PASS==0,"YES","NO")

summary(my_model)

my_pred = predict(my_model,my_data2,type="response")

my_pred = ifelse(my_pred>0.5,1,0)

my_pred = ifelse(my_pred==1,"NO PASS","PASS")

table(my_data2$PASS)
table(my_pred,my_data2$PASS)
table(my_pred)

length(unique(golden$chr_pos))

# without CHR and BP
#sensitivity
model_sens = (table(my_pred,my_data2$PASS)[2,2]/nrow(golden))*100
#specificity
model_spec = (table(my_pred,my_data2$PASS)[2,2])/(table(my_pred,my_data2$PASS)[2,2]+
                                                    table(my_pred,my_data2$PASS)[2,1])*100

f_score = ((2*model_sens*model_spec)/(model_sens + model_spec))/100

my_data2$prediction = my_pred


### genotype error in the model ######

my_data2$GT = NA

data_call = my_data2 %>% filter(prediction=="PASS") %>% filter(PASS=="YES") %>% 
  dplyr::select(GT_gatk,
                GT_whamg,
                GT_svaba,
                GT_pindel,
                GT_delly,
                GT_strelka,
                GT_deepvariant,
                GT_manta,
                GT,GT_golden,PASS,prediction)



data_call$GT_deepvariant[data_call$GT_deepvariant=="9/9"] = NA
data_call$GT_gatk[data_call$GT_gatk=="9/9"] = NA
data_call$GT_strelka[data_call$GT_strelka=="9/9"] = NA
data_call$GT_delly[data_call$GT_delly=="9/9"] = NA
data_call$GT_pindel[data_call$GT_pindel=="9/9"] = NA
data_call$GT_svaba[data_call$GT_svaba=="9/9"] = NA
data_call$GT_whamg[data_call$GT_whamg=="9/9"] = NA
data_call$GT_manta[data_call$GT_manta=="9/9"] = NA


table(data_call$GT_gatk)
table(data_call$GT_whamg)
table(data_call$GT_svaba)
table(data_call$GT_pindel)
table(data_call$GT_delly)
table(data_call$GT_strelka)
table(data_call$GT_deepvariant)
table(data_call$GT_manta)


data_call = as.matrix(data_call)

for(i in 1:nrow(data_call)){
  
  my_tab = table(as.character(data_call[i,c(1:8)]))
  
  if(length(which(my_tab==max(my_tab)))==1){
    data_call[i,9] = names(which.max(my_tab))
  }
  
  if(length(which(my_tab==max(my_tab)))>1){
    
    data_call[i,9] = NA
  }
}


head(data_call)
dim(data_call)


# use genotype from GATK because SVABA is bad

data_call[is.na(data_call[,9]),9] = data_call[is.na(data_call[,9]),1] 

data_call = as.data.frame(data_call)

table(data_call$GT,useNA = "always")

# concordance genotype of the model with golden ######

data_call = data_call[!is.na(data_call$GT),]

model_g_error = 100 -(sum(as.character(data_call$GT)==as.character(data_call$GT_golden))/
                        nrow(data_call))*100




