library(data.table)
library(dplyr)
library(R.utils)

args <- commandArgs(trailingOnly = TRUE)
batch_job <- args[1] 

#batch_job = 1000

options(digits = 16,scipen = 16)

source("/gpfs/projects/bsc05/jordivalls/GCAT_project_all_samples/merge_all_calling_GCAT//Deleciones/functions.R")

regions = fread("/gpfs/projects/bsc05/jordivalls/GCAT_project_all_samples/merge_all_calling_GCAT/length_chromosomes_for_merge.txt")

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

ids = fread("/gpfs/projects/bsc05/jordivalls/GCAT_project_all_samples/merge_all_calling_GCAT//Deleciones/all_samplesok",header = F)

ids = ids$V1

#ids = ids[1:785]

chunk = region_out[batch_job,]

chr = as.numeric(chunk$CHR)
start = as.numeric(chunk$start)
end = as.numeric(chunk$end)
        
# chr = 1
# start = 14400001
# end = 14450000

n_row = 0

empty_samples = 0

sample_initial=1 # start sample 1

my_samples = data.frame()

while(n_row==0 & sample_initial<786){
  
  sample1 = fread(paste0("/gpfs/projects/bsc05/jordivalls/GCAT_project_all_samples/merge_all_calling_GCAT/Inversion/merge_callers/",
                         ids[sample_initial],"/",ids[sample_initial],"_Inv_chr_",chr),
                  header = T,fill = T)
  
  sample1$PL = "0,255,255"
  sample1$PL[sample1[,3]=="0/1"] = "255,0,255"
  sample1$PL[sample1[,3]=="1/1"] = "255,255,0"
  sample1$PL[sample1[,3]=="./."] = "0,255,255"
  
  colnames(sample1)[15] = paste0("PL_",ids[sample_initial])
  
  sample1 = sample1[order(sample1[,2]),]
  
  pos_sample = as.numeric(t(sample1[,2]))
  
  my_samples1 = sample1[which(pos_sample>=start & pos_sample<=end),]
  
  n_row = nrow(my_samples1)
  
  if(n_row==0){
    
    my_samples = cbind(my_samples,my_samples1) # merge null samples
    
    empty_samples = empty_samples + 1
    
  }
  
  if(n_row>0){
    
    my_samples1$ID = paste0("inversion_",1:nrow(my_samples1))
    
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


dir.create("/gpfs/projects/bsc05/jordivalls/GCAT_project_all_samples/merge_all_calling_GCAT/Inversion/merge_samples_new/empty_batches/")

if(sample_initial==786){write.table(region_out[batch_job,],paste0("/gpfs/projects/bsc05/jordivalls/GCAT_project_all_samples/merge_all_calling_GCAT/Inversion/merge_samples_new/empty_batches/chr_",chr,"_start_",start,"_end_",end),row.names=F,quote=F)
}

if(sample_initial<786){
  
  for(sample_merge in sample_initial:785){
    
    sample1 = fread(paste0("/gpfs/projects/bsc05/jordivalls/GCAT_project_all_samples/merge_all_calling_GCAT/Inversion/merge_callers/",
                           ids[sample_merge],"/",ids[sample_merge],"_Inv_chr_",chr),header = T,fill = T)
    
    sample1$PL = "0,255,255"
    sample1$PL[sample1[,3]=="0/1"] = "255,0,255"
    sample1$PL[sample1[,3]=="1/1"] = "255,255,0"
    sample1$PL[sample1[,3]=="./."] = "0,255,255"
    
    colnames(sample1)[15] = paste0("PL_",ids[sample_merge])
    
    sample1 = sample1[order(sample1[,2]),]
    
    pos_sample = as.numeric(t(sample1[,2]))
    
    sample1 = sample1[which(pos_sample>=start & pos_sample<=end),]
    
    n_row = nrow(sample1)
    
    if(n_row>0){
      
      sample1$ID = "none"
      
      my_samples = merge_samples(sample1,my_samples,ids[sample_merge],ids[1:(sample_merge-1)])
      
      my_samples = unique(my_samples)
      
      
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
      my_samples$V14 =  NA
      my_samples$V15 =  NA
      #my_samples$V16 =  NA
  
      
      
      colnames(my_samples)[(ncol(my_samples)-14):ncol(my_samples)] = paste0(c("chr_", # N-1 ncol
                                                                              "start_",
                                                                              "GT_",
                                                                              "length_",
                                                                              "strategy_",
                                                                              "reciprocity_",
                                                                              "callers_detected_",
                                                                              "name_callers_",
                                                                              "GT_callers_",
                                                                              "min_window_size_",
                                                                              "lower_",
                                                                              "upper_",
                                                                              "PASS_",
                                                                              "PASS_num_",
                                                                              "PL_"),ids[sample_merge])
      
    }
    print(paste0("The number of variants of the merging is: ",nrow(my_samples)))
    print(paste0("The number of samples of the merging is: ",sample_merge))
    print(paste0("The number of empty samples is: ",empty_samples))
    
  }
  
  dim(my_samples[,paste0("chr_",ids),with=F])
  dim(my_samples[,paste0("start_",ids),with=F])
  dim(my_samples[,paste0("GT_",ids),with=F])
  dim(my_samples[,paste0("length_",ids),with=F])
  dim(my_samples[,paste0("strategy_",ids),with=F])
  dim(my_samples[,paste0("reciprocity_",ids),with=F])
  dim(my_samples[,paste0("callers_detected_",ids),with=F])
  dim(my_samples[,paste0("name_callers_",ids),with=F])
  dim(my_samples[,paste0("GT_callers_",ids),with=F])
  dim(my_samples[,paste0("min_window_size_",ids),with=F])
  dim(my_samples[,paste0("lower_",ids),with=F])
  dim(my_samples[,paste0("upper_",ids),with=F])
  dim(my_samples[,paste0("PASS_",ids),with=F])
  dim(my_samples[,paste0("PASS_num_",ids),with=F])
  dim(my_samples[,paste0("PL_",ids),with=F])
  
  my_samples[,-c(paste0("chr_",ids),
                    paste0("start_",ids),
                    paste0("GT_",ids),
                    paste0("length_",ids),
                    paste0("strategy_",ids),
                    paste0("reciprocity_",ids),
                    paste0("callers_detected_",ids),
                    paste0("name_callers_",ids),
                    paste0("GT_callers_",ids),
                    paste0("min_window_size_",ids),
                    paste0("lower_",ids),
                    paste0("upper_",ids),
                    paste0("PASS_",ids),
                    paste0("PASS_num_",ids),
                    paste0("PL_",ids)),with=F]
  
  
  my_samples$ID = tstrsplit(my_samples$ID,split="\\.")[[1]]
  
  # consensus #####
  
  final_merge = data.frame(ID = unique(my_samples$ID))
  
  final_merge$CHR = chr
  
  final_merge$start = consensus_merge_samples(my_samples,unique(my_samples$ID),"start_")
  final_merge$start_error = consensus_merge_samples(my_samples,unique(my_samples$ID),"min_window_size_")
  final_merge$length = consensus_merge_samples(my_samples,unique(my_samples$ID),"length_")
  
  
  genos = my_samples %>%
    group_by(ID) %>% 
    summarize_at(paste0("GT_",ids),~paste(unique(na.omit(.)), collapse = ',')) %>%
    as.matrix()
  
  pls = my_samples %>%
    group_by(ID) %>% 
    summarize_at(paste0("PL_",ids),~paste(unique(na.omit(.)), collapse = ',')) %>%
    as.data.frame()
  
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
  
  repro = my_samples %>%
    group_by(ID) %>% 
    summarize_at(paste0("reciprocity_",ids),~median(.,na.rm=T)) %>%
    as.data.frame()
  
  #genos_best = my_samples %>%
    #group_by(ID) %>%
    #summarize_at(paste0("GT_best_",ids),~paste(unique(na.omit(.)), collapse = ',')) %>%
    #as.data.frame()
  
  
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
  final_merge = left_join(final_merge,pls) 
  final_merge = left_join(final_merge,ncalls) 
  final_merge = left_join(final_merge,name_calls) 
  final_merge = left_join(final_merge,strategy_var) 
  final_merge = left_join(final_merge,pass) 
  final_merge = left_join(final_merge,pass_num) 
  final_merge = left_join(final_merge,repro) 
  #final_merge = left_join(final_merge,genos_best)
  
  final_merge$start = as.numeric(floor(final_merge$start))
  
  final_merge = final_merge %>% arrange(start)
  
  print(dim(final_merge))
  
  dir.create(paste0("/gpfs/projects/bsc05/jordivalls/GCAT_project_all_samples/merge_all_calling_GCAT/Inversion/merge_samples_new/chr_",chr))
  
  final_merge %>% 
    fwrite(paste0("/gpfs/projects/bsc05/jordivalls/GCAT_project_all_samples/merge_all_calling_GCAT/Inversion/merge_samples_new/chr_",chr,
                  "/Inv_chr_",chr,"_",as.character(start),"_",as.character(end)),row.names=F,sep=" ",na = "NA")
  
  print("Now compress pleaseeeeeeeeeeeeeeee")
  
  system(paste0("gzip /gpfs/projects/bsc05/jordivalls/GCAT_project_all_samples/merge_all_calling_GCAT/Inversion/merge_samples_new/chr_",chr,
                "/Inv_chr_",chr,"_",as.character(start),"_",as.character(end)))
  
  print("Now generate VCF")
  
  final_merge = fread(paste0("zcat /gpfs/projects/bsc05/jordivalls/GCAT_project_all_samples/merge_all_calling_GCAT/Inversion/merge_samples_new/chr_",chr,
                             "/Inv_chr_",chr,"_",as.character(start),"_",as.character(end),".gz"))
  
  # warning, some variants repeated in the VCF
  
  ## make VCF #####
  
  final_merge2 = as.data.frame(final_merge)
  
  # replace "./." for the best genotype for variants with <5% missing
  
  # miss_calculation = final_merge2[,paste0(c("GT_"),ids)]
  # 
  # miss_rate = NULL
  # 
  # for(i in 1:nrow(miss_calculation)){
  #   
  #   aux = sum(miss_calculation[i,]=="./.",na.rm = T)/ncol(miss_calculation)
  #   
  #   miss_rate = c(miss_rate,aux)
  # }
  # 
  # which_miss_05 = which(miss_rate<0.05)
  # 
  # for(sample_vcf in 1:785){
  #   
  #   aux2 = final_merge2[which_miss_05,
  #                       paste0(c("GT_"),ids[sample_vcf])]
  #   
  #   aux = final_merge2[which_miss_05,
  #                      paste0(c("GT_best_"),ids[sample_vcf])]
  #   
  #   which_best = which(t(aux) %in% c("0/0","0/1","1/1"))
  #   
  #   if(length(which_best)>0){
  #     
  #     final_merge2[which_miss_05,
  #                  paste0(c("GT_"),ids[sample_vcf])][which_best] = as.character(aux[which_best])
  #     
  #   }
  #   
  #   print(ids[sample_vcf])
  # }
  # 
  
  # AC allele counts ALT ####
  
  aux = final_merge2[,paste0("GT_",ids)]
  
  final_merge2$AC = apply(aux,1,function(x){
    
    alleles = c(tstrsplit(x,split="/")[[1]],
                tstrsplit(x,split="/")[[2]])
    
    ac = sum(alleles==1)
    
    return(paste0("AC=",ac,";"))
    
  })
  
  
  # MAC minor allele count ####
  
  final_merge2$MAC = apply(aux,1,function(x){
    
    alleles = c(tstrsplit(x,split="/")[[1]],tstrsplit(x,split="/")[[2]])
    
    ac1 = sum(alleles==1)
    ac2 = sum(alleles==0)
    
    return(paste0("MAC=",min(ac1,ac2),";"))
    
  })
  
  
  # AF allele frequency ALT ####
  
  final_merge2$AF = apply(aux,1,function(x){
    
    alleles = c(tstrsplit(x,split="/")[[1]],tstrsplit(x,split="/")[[2]])
    
    af = sum(alleles==1)
    
    total = sum(alleles==1) + sum(alleles==0)
    
    af = af/total
    
    return(paste0("AF=",round(af,4),";"))
    
  })
  
  
  # AN allele number ####
  
  final_merge2$AN = apply(aux,1,function(x){
    
    alleles = c(tstrsplit(x,split="/")[[1]],tstrsplit(x,split="/")[[2]])
    
    total = sum(alleles==1) + sum(alleles==0)
    
    return(paste0("AN=",total,";"))
    
  })
  
  
  # MAF ####
  
  final_merge2$MAF = apply(aux,1,function(x){
    
    alleles = c(tstrsplit(x,split="/")[[1]],tstrsplit(x,split="/")[[2]])
    
    af = sum(alleles==1)
    
    total = sum(alleles==1) + sum(alleles==0)
    
    af = af/total
    
    maf = ifelse(af>0.5,1-af,af)
    
    return(paste0("MAF=",round(maf,4),";"))  
  })
  
  
  # END #####
  
  final_merge2$END = paste0("END=",final_merge2$start + final_merge2$length,";") 
  
  # SVLEN #####
  
  final_merge2$SVLEN = paste0("SVLEN=",final_merge2$length,";") 
  
  
  # ERRORBKP ####
  
  final_merge2$ERRBKP = paste0("ERRBKP=+-",final_merge2$start_error,";")
  
  
  # SVTYPE ####
  
  final_merge2$SVTYPE = paste0("SVTYPE=INV;")
  
  # POPVAR 
  
  maf = apply(aux,1,function(x){
    
    alleles = c(tstrsplit(x,split="/")[[1]],tstrsplit(x,split="/")[[2]])
    
    af = sum(alleles==1)
    
    total = sum(alleles==1) + sum(alleles==0)
    
    af = af/total
    
    maf = ifelse(af>0.5,1-af,af)
    
    return(round(maf,4))  
  })
  
  
  final_merge2$POPVAR = ""
  
  final_merge2$POPVAR[which(maf>=0.05)] = "POPVAR=COMMON"
  final_merge2$POPVAR[which(maf<0.05 & maf>=0.01)] = "POPVAR=LOW_FREQ"
  final_merge2$POPVAR[which(maf<0.01)] = "POPVAR=RARE"
  final_merge2$POPVAR[which(maf==0)] = "POPVAR=MONOMORPHIC"
  
  
  genotype_numbers = apply(aux,1,function(x){
    
    geno_num = sum(!x %in% c("./.","0/0"))
    
    return(geno_num)
    
  })
  
  final_merge2$POPVAR[which(genotype_numbers==1)] = "POPVAR=SINGLETON"
  final_merge2$POPVAR[which(genotype_numbers==2)] = "POPVAR=DOUBLETON"
  
  
  
  # QUAL #####
  
  aux = final_merge2[,paste0("PASS_",ids)]
  
  final_merge2$qual = apply(aux,1,function(x){
    
    pass = sum(x=="PASS",na.rm=T)
    no_pass = sum(x=="NO_PASS",na.rm=T)
    
    return(round(pass/(pass+no_pass),3))
    
  })
  
  # PASS #####
  
  final_merge2$pass = apply(aux,1,function(x){
    
    pass = sum(x=="PASS",na.rm=T)
    
    return(pass)
    
  })
  
  # NO_PASS #####
  
  final_merge2$no_pass = apply(aux,1,function(x){
    
    no_pass = sum(x=="NO_PASS",na.rm=T)
    
    return(no_pass)
    
  })
  
  final_merge2$pass = paste0(";N_PASS=",final_merge2$pass,";") 
  final_merge2$no_pass = paste0("N_NOPASS=",final_merge2$no_pass,";") 
  
  
  # Genotype counts ####
  
  aux = final_merge2[,paste0("GT_",ids)]
  
  final_merge2$Geno_00 = apply(aux,1,function(x){
    
    return(sum(x=="0/0",na.rm=T))
    
  })
  
  final_merge2$Geno_01 = apply(aux,1,function(x){
    
    return(sum(x=="0/1",na.rm=T))
    
  })
  
  
  final_merge2$Geno_11 = apply(aux,1,function(x){
    
    return(sum(x=="1/1",na.rm=T))
    
  })
  
  final_merge2$Geno_.. = apply(aux,1,function(x){
    return(sum(x=="./.",na.rm=T))
  })
  
  
  final_merge2$Geno_00 = paste0("GENO_00=",final_merge2$Geno_00,";") 
  final_merge2$Geno_01 = paste0("GENO_01=",final_merge2$Geno_01,";") 
  final_merge2$Geno_11 = paste0("GENO_11=",final_merge2$Geno_11,";") 
  final_merge2$Geno_.. = paste0("GENO_..=",final_merge2$Geno_..,";")
  
  
  # Caller stats #######
  
  aux = final_merge2[,paste0("callers_detected_",ids)]
  
  final_merge2$callers_tab = apply(aux,1,function(x){
    
    aux = factor(as.numeric(x))
    
    my_c = as.numeric(levels(aux))
    
    no_calls = which(!1:7 %in% my_c)
    
    levels(aux) = c(levels(aux),no_calls)
    
    return(paste0(table(aux)[order(names(table(aux)))],collapse = ","))
    
  })
  
  final_merge2$min_calls = apply(aux,1,function(x){
    
    return(min(as.numeric(x),na.rm=T))
    
  })
  
  final_merge2$max_calls = apply(aux,1,function(x){
    
    return(max(as.numeric(x),na.rm=T))
    
  })
  
  final_merge2$median_calls = apply(aux,1,function(x){
    
    return(median(as.numeric(x),na.rm=T))
    
  })
  
  final_merge2$callers_tab = paste0("TAB_CALLERS=",final_merge2$callers_tab,";") 
  final_merge2$min_calls = paste0("MIN_CALLERS=",final_merge2$min_calls,";") 
  final_merge2$max_calls = paste0("MAX_CALLERS=",final_merge2$max_calls,";") 
  final_merge2$median_calls = paste0("MEDIAN_CALLERS=",final_merge2$median_calls) 
  
  
  # CREATE VCF
  
  my_vcf = data.frame("CHR"=final_merge2$CHR)
  
  colnames(my_vcf) = "#CHROM"
  
  my_vcf$POS = final_merge2$start
  my_vcf$ID = paste0(chr,":",my_vcf$POS,":A:<INV>")
  my_vcf$REF = "A"
  my_vcf$ALT = "<INV>"
  my_vcf$QUAL = final_merge2$qual
  my_vcf$FILTER = ifelse(final_merge2$qual>0.5,"PASS","NO_PASS")
  
  my_vcf$INFO = do.call(paste0,list(final_merge2$AC,
                                    final_merge2$MAC,
                                    final_merge2$AF,
                                    final_merge2$AN,
                                    final_merge2$MAF,
                                    final_merge2$END,
                                    final_merge2$SVLEN,
                                    final_merge2$ERRBKP,
                                    final_merge2$SVTYPE,
                                    final_merge2$POPVAR,
                                    final_merge2$pass,
                                    final_merge2$no_pass,
                                    final_merge2$Geno_00,
                                    final_merge2$Geno_01,
                                    final_merge2$Geno_11,
                                    final_merge2$Geno_..,
                                    final_merge2$callers_tab,
                                    final_merge2$min_calls,
                                    final_merge2$max_calls,
                                    final_merge2$median_calls))
  
  my_vcf$FORMAT = "GT:PL:LM:PVM:SVMETHOD:SVCALL:SVSTRAT:SVREPRO"
  
  for(sample_vcf in 1:785){
    
    sample1 = final_merge2[,paste0(c("GT_",
                                     "PL_",
                                     "PASS_num_",
                                     "PASS_",
                                     "name_callers_",
                                     "callers_detected_",
                                     "strategy_",
                                     "reciprocity_"),ids[sample_vcf])]
    
    sample1[,3] = round(sample1[,3],4)
    
    sample1[is.na(sample1)] = ""
    
    sample1[,2][sample1[,1]=="0/0"] = "0,255,255"
    sample1[,2][sample1[,1]=="0/1"] = "255,0,255"
    sample1[,2][sample1[,1]=="1/1"] = "255,255,0"
    
    sample1 = do.call(paste0,list(sample1[,1],":",
                                  sample1[,2],":",
                                  sample1[,3],":",
                                  sample1[,4],":",
                                  sample1[,5],":",
                                  sample1[,6],":",
                                  sample1[,7],":",
                                  sample1[,8]))
    
    
    my_vcf = cbind(my_vcf,sample1)
    
    colnames(my_vcf)[ncol(my_vcf)] = ids[sample_vcf]
    
  }
  
  
  fwrite(my_vcf,paste0("/gpfs/projects/bsc05/jordivalls/GCAT_project_all_samples/merge_all_calling_GCAT/Inversion/merge_samples_new/chr_",chr,
                       "/Inv_chr_",chr,"_",as.character(start),"_",as.character(end),".vcf"),row.names=F,sep="\t")
  
  print("Now compress again pleaseeeeeeeeeeeeeeee")
  
  
  system(paste0("gzip /gpfs/projects/bsc05/jordivalls/GCAT_project_all_samples/merge_all_calling_GCAT/Inversion/merge_samples_new/chr_",chr,
                "/Inv_chr_",chr,"_",as.character(start),"_",as.character(end),".vcf"))
  
  
  
}



  