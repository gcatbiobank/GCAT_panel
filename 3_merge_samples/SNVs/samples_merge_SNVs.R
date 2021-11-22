library(data.table)
library(dplyr)
library(R.utils)

args <- commandArgs(trailingOnly = TRUE)
batch_job <- args[1] 
num_samples <- args[2]
real_num_samples <- as.integer(num_samples)
num_samples = as.integer(num_samples) +1
#batch_job = 57606
options(digits = 16,scipen = 16)
# chr22 batch 56586, 57611
path_inputs <- args[3]
path_inputs = "/home/jordi/merge_samples/SNVs/inputs/"
#source("/home/jordi/merge_samples/functions.R")


## create batches #####

regions = fread("./length_chromosomes_for_merge.txt")

batch = 5e4

region_out = NULL

for(batches_chr in 1:nrow(regions)){
  
  if(as.numeric(regions$pos[batches_chr])>batch){
    
    aux = seq(1,regions$pos[batches_chr],by = batch )
    
    aux2 = cbind(aux[-length(aux)],as.numeric(aux[-1]-1))
    
    aux2 = cbind(rep(regions$chr[batches_chr],nrow(aux2)),as.data.frame(aux2))
    
    colnames(aux2) = c("CHR","start","end")
    
    region_out = rbind(region_out,aux2)
    
    colnames(region_out) = c("CHR","start","end")
    
  }
  if(as.numeric(regions$pos[batches_chr])<=batch){
    
    aux3 = data.frame(CHR=regions$chr[batches_chr],start=1,end=regions$pos[batches_chr])
    print(aux3)
    region_out = rbind(region_out,aux3)
    
  }
  
}

ids = fread("./samplesok",header = F)

ids = ids$V1



print(ids)

chunk = region_out[batch_job,]

chr = as.numeric(chunk$CHR)
start = as.numeric(chunk$start)
end = as.numeric(chunk$end)

#start = 18720682
#end = 18720682+50000

n_row = 0

empty_samples = 0

sample_initial=1 # start sample 1

my_samples = data.frame()

while(n_row==0 & sample_initial<num_samples){#2
  
  sample1 = fread(paste0(path_inputs,"/",
                         ids[sample_initial],"/",ids[sample_initial],"_SNPs_chr_",chr),
                  header = T)
  
  sample1 = sample1[order(sample1[,2]),]
  
  pos_sample = as.numeric(t(sample1[,2]))
  
  my_samples1 = sample1[which(pos_sample>=start & pos_sample<=end),]
  
  n_row = nrow(my_samples1)
  
  if(n_row==0){
    
    my_samples = cbind(my_samples,my_samples1)
    
    empty_samples = empty_samples + 1
    
  }
  
  if(n_row>0){
    
    my_samples1$ID = do.call(paste0,list(t(my_samples1[,1]),":",t(my_samples1[,2]),":",t(my_samples1[,3]),":",t(my_samples1[,4])))
    
    names_col = colnames(my_samples)
    
    if(sample_initial==1){my_samples = my_samples1}
    
    if(sample_initial!=1){
      
      my_samples = rbind(my_samples,matrix(NA,nrow = nrow(my_samples1),ncol=ncol(my_samples)))
      
      colnames(my_samples) = names_col
      
      my_samples = cbind(my_samples,my_samples1)
    }
    
  }
  
  sample_initial = sample_initial + 1
}

my_samples = as.data.table(my_samples)


if(sample_initial==num_samples){write.table(region_out[batch_job,],paste0("./outputs/empty_batches/chr_",chr,"_start_",start,"_end_",end),row.names=F,quote=F)
}

if(sample_initial<num_samples){
  
  for(sample_merge in sample_initial:real_num_samples){
    
    sample1 = fread(paste0(path_inputs,"/",
                           ids[sample_merge],"/",ids[sample_merge],"_SNPs_chr_",chr),header = T)
    
    
    sample1 = sample1[order(sample1[,2]),]
    
    pos_sample = as.numeric(t(sample1[,2]))
    
    sample1 = sample1[which(pos_sample>=start & pos_sample<=end),]
    
    n_row = nrow(sample1)
    
    if(n_row>0){
      
      sample1$ID = do.call(paste0,list(t(sample1[,1]),":",t(sample1[,2]),":",t(sample1[,3]),":",t(sample1[,4])))
      
      my_samples= full_join(my_samples,sample1) ### merge samples
      
      if(sum(duplicated(my_samples$ID))>0){
        
        my_samples = my_samples[order(my_samples$ID),]
        
        ids_merge = my_samples$ID
        
        unique_ids = unique(ids_merge)
        
        for(id_loop in 1:length(unique_ids)){
          
          length_id = length(ids_merge[ids_merge %in% unique_ids[id_loop]])
          
          if(length_id>1){
            ids_merge[ids_merge %in% unique_ids[id_loop]] = paste0(ids_merge[ids_merge %in% unique_ids[id_loop]],".",1:length_id)
          }}
        
        my_samples$ID = ids_merge
        
      }
      
    }
    
    if(n_row==0){
      
      empty_samples = empty_samples + 1
      
      my_samples$V1 =  NA
      my_samples$V2 =  NA
      my_samples$V3 =  NA
      my_samples$V4 =  NA
      my_samples$V5 =  NA
      my_samples$V6 =  NA
      my_samples$V7 =  NA
      my_samples$V8 =  NA
      my_samples$V9 =  NA
      my_samples$V10 =  NA
      my_samples$V11 =  NA
      my_samples$V12 =  NA
      my_samples$V13 =  NA

      
      colnames(my_samples)[(ncol(my_samples)-12):ncol(my_samples)] = paste0(c("chr_", # 12, per N-1
                                                                              "position_",
                                                                              "REF_",
                                                                              "ALT_",
                                                                              "GT_",
                                                                              "strategy_",
                                                                              "callers_detected_",
                                                                              "name_callers_",
                                                                              "GQ_deep_gatk_strelka_",
                                                                              "AD_consensus_",
                                                                              "DP_consensus_",
                                                                              "PASS_",
                                                                              "PASS_num_"),ids[sample_merge])
      
    }
    print(paste0("The number of variants of the merging is: ",nrow(my_samples)))
    print(paste0("The number of samples of the merging is: ",sample_merge))
    print(paste0("The number of empty samples is: ",empty_samples))
    
  }
  
  my_samples$ID = tstrsplit(my_samples$ID,split="\\.")[[1]]
  
  # consensus #####
  
  final_merge = data.frame(ID = unique(my_samples$ID))
  
  final_merge$CHR = chr
  final_merge$start = tstrsplit(unique(my_samples$ID),split=":")[[2]]
  
  
  genos = my_samples %>%
    group_by(ID) %>% 
    summarize_at(paste0("GT_",ids),~paste(unique(na.omit(.)), collapse = ',')) %>%
    as.matrix()
  
  ncalls = my_samples %>%
    group_by(ID) %>% 
    summarize_at(paste0("callers_detected_",ids),~median(.,na.rm=T)) %>%
    as.data.frame()
  
  name_calls = my_samples %>%
    group_by(ID) %>% 
    summarize_at(paste0("name_callers_",ids),~paste(unique(na.omit(.)), collapse = ',')) %>%
    as.data.frame()
  
  strategy_var = my_samples %>%
    group_by(ID) %>% 
    summarize_at(paste0("strategy_",ids),~paste(unique(na.omit(.)), collapse = ',')) %>%
    as.data.frame()
  
  pass = my_samples %>%
    group_by(ID) %>% 
    summarize_at(paste0("PASS_",ids),~paste(unique(na.omit(.)), collapse = ',')) %>%
    as.data.frame()
  
  pass_num = my_samples %>%
    group_by(ID) %>% 
    summarize_at(paste0("PASS_num_",ids),~median(.,na.rm=T)) %>%
    as.data.frame()
  
  GQ = my_samples %>% group_by(ID) %>% summarize_at(paste0("GQ_deep_gatk_strelka_",ids),~paste(unique(na.omit(.)), collapse = ',')) %>%
    as.data.frame()
  
  AD = my_samples %>% group_by(ID) %>% summarize_at(paste0("AD_consensus_",ids),~paste(unique(na.omit(.)), collapse = ',')) %>%
    as.data.frame()
  
  DP = my_samples %>% group_by(ID) %>% summarize_at(paste0("DP_consensus_",ids),~median(.,na.rm=T)) %>%
    as.data.frame()
  
  # add 0/0 genotype 
  
  if(nrow(genos)==1){
    
    genos[which(!(genos %in% c("./.","0/1","1/1")))[-1]] = "0/0"
    
    genos = as.data.frame(genos)
    
    genos = genos %>% mutate_at(c(paste0("GT_",ids)), funs(as.character(.)))
  }
  
  if(nrow(genos)>1){
    
    id = genos[,1]
    genos = genos[,-1]
    
    genos[!genos %in% c("./.","0/1","1/1")] = ""
    genos[genos %in% c("")] = "0/0"
    
    genos = as.data.frame(cbind(ID=id,genos))
  }
  
  # left join final merge
  
  final_merge = left_join(final_merge,genos) 
  final_merge = left_join(final_merge,ncalls) 
  final_merge = left_join(final_merge,name_calls) 
  final_merge = left_join(final_merge,strategy_var) 
  final_merge = left_join(final_merge,pass) 
  final_merge = left_join(final_merge,pass_num)
  final_merge = left_join(final_merge,GQ)
  final_merge = left_join(final_merge,AD)
  final_merge = left_join(final_merge,DP)
  
  final_merge = final_merge %>% arrange(start)
  

  
  print(dim(final_merge))
  
  final_merge %>% 
    fwrite(paste0("./outputs/chr_",chr,
                  "/SNPs_chr_",chr,"_",as.character(start),"_",as.character(end)),row.names=F,sep=" ",na = "NA")
  
  print("Now compress please")
  
  system(paste0("gzip ./outputs/chr_",chr,
              "/SNPs_chr_",chr,"_",as.character(start),"_",as.character(end)))
  
  print("Now generate VCF")
  
  final_merge = fread(paste0("zcat ./outputs/chr_",chr,
                             "/SNPs_chr_",chr,"_",as.character(start),"_",as.character(end),".gz"))
  
  final_merge = final_merge %>% as.data.frame()
  dim(final_merge)
  
  ## make VCF #####
  
  # AC ####
  
  aux = final_merge[,paste0("GT_",ids)]
  
  final_merge$AC = apply(aux,1,function(x){
    
    alleles = c(tstrsplit(x,split="/")[[1]],
                tstrsplit(x,split="/")[[2]])
    
    ac = sum(alleles==1)
    
    return(paste0("AC=",ac,";"))
    
  })
  
  
  # MAC ####
  
  final_merge$MAC = apply(aux,1,function(x){
    
    alleles = c(tstrsplit(x,split="/")[[1]],tstrsplit(x,split="/")[[2]])
    
    ac1 = sum(alleles==1)
    ac2 = sum(alleles==0)
    
    return(paste0("MAC=",min(ac1,ac2),";"))
    
  })
  
  
  # AF ####
  
  final_merge$AF = apply(aux,1,function(x){
    
    alleles = c(tstrsplit(x,split="/")[[1]],tstrsplit(x,split="/")[[2]])
    
    af = sum(alleles==1)
    
    total = sum(alleles==1) + sum(alleles==0)
    
    af = af/total
    
    return(paste0("AF=",round(af,4),";"))
    
  })
  
  
  # AN ####
  
  final_merge$AN = apply(aux,1,function(x){
    
    alleles = c(tstrsplit(x,split="/")[[1]],tstrsplit(x,split="/")[[2]])
    
    total = sum(alleles==1) + sum(alleles==0)
    
    return(paste0("AN=",total,";"))
    
  })
  
  
  # MAF ####
  
  final_merge$MAF = apply(aux,1,function(x){
    
    alleles = c(tstrsplit(x,split="/")[[1]],tstrsplit(x,split="/")[[2]])
    
    af = sum(alleles==1)
    
    total = sum(alleles==1) + sum(alleles==0)
    
    af = af/total
    
    maf = ifelse(af>0.5,1-af,af)
    
    return(paste0("MAF=",round(maf,4),";"))  
  })
  
  
  # SVTYPE ####
  
  final_merge$SVTYPE = paste0("SVTYPE=SNP;")
  
  # POPVAR 
  
  maf = apply(aux,1,function(x){
    
    alleles = c(tstrsplit(x,split="/")[[1]],tstrsplit(x,split="/")[[2]])
    
    af = sum(alleles==1)
    
    total = sum(alleles==1) + sum(alleles==0)
    
    af = af/total
    
    maf = ifelse(af>0.5,1-af,af)
    
    return(round(maf,4))  
  })
  
  
  final_merge$POPVAR = ""
  
  final_merge$POPVAR[which(maf>=0.05)] = "POPVAR=COMMON"
  final_merge$POPVAR[which(maf<0.05 & maf>=0.01)] = "POPVAR=LOW_FREQ"
  final_merge$POPVAR[which(maf<0.01)] = "POPVAR=RARE"
  final_merge$POPVAR[which(maf==0)] = "POPVAR=MONOMORPHIC"
  
  
  genotype_numbers = apply(aux,1,function(x){
    
    geno_num = sum(!x %in% c("./.","0/0"))
    
    return(geno_num)
    
  })
  
  final_merge$POPVAR[which(genotype_numbers==1)] = "POPVAR=SINGLETON"
  final_merge$POPVAR[which(genotype_numbers==2)] = "POPVAR=DOUBLETON"
  
  
  
  # QUAL #####
  
  aux = final_merge[,paste0("PASS_",ids)]
  
  final_merge$qual = apply(aux,1,function(x){
    
    pass = sum(x=="PASS",na.rm=T)
    no_pass = sum(x=="NO_PASS",na.rm=T)
    
    return(round(pass/(pass+no_pass),3))
    
  })
  
  # PASS #####
  
  final_merge$pass = apply(aux,1,function(x){
    
    pass = sum(x=="PASS",na.rm=T)

    return(pass)
    
  })
  
  # NO_PASS #####
  
  final_merge$no_pass = apply(aux,1,function(x){
    
    no_pass = sum(x=="NO_PASS",na.rm=T)
    
    return(no_pass)
    
  })

  final_merge$no_var = real_num_samples - (final_merge$pass + final_merge$no_pass)
  
  final_merge$pass = paste0(";N_PASS=",final_merge$pass,";") 
  final_merge$no_pass = paste0("N_NOPASS=",final_merge$no_pass,";") 
  
  # Genotype counts ####
  
  aux = final_merge[,paste0("GT_",ids)]
  
  final_merge$Geno_00 = apply(aux,1,function(x){
    
    return(sum(x=="0/0",na.rm=T))
    
  })
  
  final_merge$Geno_01 = apply(aux,1,function(x){
    
    return(sum(x=="0/1",na.rm=T))
    
  })
  
  
  final_merge$Geno_11 = apply(aux,1,function(x){
    
    return(sum(x=="1/1",na.rm=T))
    
  })
  final_merge$Geno_.. = apply(aux,1,function(x){
    return(sum(x=="./.",na.rm=T))
  })
  
  final_merge$Geno_00 = paste0("GENO_00=",final_merge$Geno_00,";")
  final_merge$Geno_01 = paste0("GENO_01=",final_merge$Geno_01,";")
  final_merge$Geno_11 = paste0("GENO_11=",final_merge$Geno_11,";")  
  final_merge$Geno_.. = paste0("GENO_..=",final_merge$Geno_..)
  
  # CREATE VCF
  
  my_vcf = data.frame("CHR"=final_merge$CHR)
  
  colnames(my_vcf) = "#CHROM"
  
  my_vcf$POS = final_merge$start
  my_vcf$ID = final_merge$ID
  my_vcf$REF = tstrsplit(final_merge$ID,split=":")[[3]]
  my_vcf$ALT = tstrsplit(final_merge$ID,split=":")[[4]]
  my_vcf$QUAL = final_merge$qual
  my_vcf$FILTER = ifelse(final_merge$qual>0.5,"PASS","NO_PASS")
  
  my_vcf$INFO = do.call(paste0,list(final_merge$AC,
                                    final_merge$MAC,
                                    final_merge$AF,
                                    final_merge$AN,
                                    final_merge$MAF,
                                    final_merge$SVTYPE,
                                    final_merge$POPVAR,
                                    final_merge$pass,
                                    final_merge$no_pass,
                                    final_merge$Geno_00,
                                    final_merge$Geno_01,
                                    final_merge$Geno_11,
                                    final_merge$Geno_..))
  
  my_vcf$FORMAT = "GT:GQ:AD:DP:LM:PVM:SVMETHOD:SVCALL:SVSTRAT"
  
  
  for(sample_vcf in 1:real_num_samples){
    
    sample1 = final_merge[,paste0(c("GT_",
                                    "GQ_deep_gatk_strelka_",
                                    "AD_consensus_",
                                    "DP_consensus_",
                                    "PASS_num_",
                                    "PASS_",
                                    "name_callers_",
                                    "callers_detected_",
                                    "strategy_"),ids[sample_vcf])]
    
    sample1[,5] = round(sample1[,5],4)
    
    sample1[is.na(sample1)] = ""
    
    sample1 = do.call(paste0,list(sample1[,1],":",
                                  sample1[,2],":",
                                  sample1[,3],":",
                                  sample1[,4],":",
                                  sample1[,5],":",
                                  sample1[,6],":",
                                  sample1[,7],":",
                                  sample1[,8],":",
                                  sample1[,9]))
    
    
    my_vcf = cbind(my_vcf,sample1)
    
    colnames(my_vcf)[ncol(my_vcf)] = ids[sample_vcf]
    
  }
  
  
  fwrite(my_vcf,paste0("./outputs/chr_",chr,
                       "/SNPs_chr_",chr,"_",as.character(start),"_",as.character(end),".vcf"),row.names=F,sep="\t")
  
  print("Now compress again please")
  
  system(paste0("cat ./header ./outputs/chr_",chr,
                "/SNPs_chr_",chr,"_",as.character(start),"_",as.character(end),".vcf > ./outputs/chr_",chr,
                "/SNPs_chr_",chr,"_",as.character(start),"_",as.character(end),".vcf1" ))
  
  system(paste0("mv ./outputs/chr_",chr,
                "/SNPs_chr_",chr,"_",as.character(start),"_",as.character(end),".vcf1 ./outputs/chr_",chr,
                "/SNPs_chr_",chr,"_",as.character(start),"_",as.character(end),".vcf" ))
  
  
  system(paste0("gzip ./outputs/chr_",chr,
              "/SNPs_chr_",chr,"_",as.character(start),"_",as.character(end),".vcf"))
  

  system(paste0("rm ./outputs/chr_",chr,
                "/SNPs_chr_",chr,"_",as.character(start),"_",as.character(end),".gz"))

  
  
  
  
}
