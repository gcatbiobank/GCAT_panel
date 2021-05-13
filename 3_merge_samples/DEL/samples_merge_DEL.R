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

chunk = region_out[batch_job,]

chr = as.numeric(chunk$CHR)
start = as.numeric(chunk$start)
end = as.numeric(chunk$end)

n_row = 0

empty_samples = 0

sample_initial=1 # start sample 1

my_samples = data.frame()

while(n_row==0 & sample_initial<786){
  
  sample1 = fread(paste0("/gpfs/projects/bsc05/jordivalls/GCAT_project_all_samples/merge_all_calling_GCAT//Deleciones/merge_callers/",
                         ids[sample_initial],"/",ids[sample_initial],"_Del_chr_",chr),
                  header = T,fill = T)
  
  sample1$PL = "0,255,255"
  sample1$PL[sample1[,3]=="0/1"] = "255,0,255"
  sample1$PL[sample1[,3]=="1/1"] = "255,255,0"
  sample1$PL[sample1[,3]=="./."] = "0,255,255"
  
  colnames(sample1)[16] = paste0("PL_",ids[sample_initial])
  
  sample1 = sample1[order(sample1[,2]),]
  
  pos_sample = as.numeric(t(sample1[,2]))
  
  my_samples1 = sample1[which(pos_sample>=start & pos_sample<=end),]
  
  n_row = nrow(my_samples1)
  
  if(n_row==0){
    
    my_samples = cbind(my_samples,my_samples1) # merge null samples
    
    empty_samples = empty_samples + 1
    
  }
  
  if(n_row>0){
    
    my_samples1$ID = paste0("deletion_",1:nrow(my_samples1))
    
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


if(sample_initial==786){write.table(region_out[batch_job,],paste0("/gpfs/projects/bsc05/jordivalls/GCAT_project_all_samples/merge_all_calling_GCAT//Deleciones/merge_samples/empty_batches/chr_",chr,"_start_",start,"_end_",end),row.names=F,quote=F)
}

if(sample_initial<786){
  
  for(sample_merge in sample_initial:785){
    
    sample1 = fread(paste0("/gpfs/projects/bsc05/jordivalls/GCAT_project_all_samples/merge_all_calling_GCAT//Deleciones/merge_callers/",
                           ids[sample_merge],"/",ids[sample_merge],"_Del_chr_",chr),header = T,fill = T)
    
    sample1$PL = "0,255,255"
    sample1$PL[sample1[,3]=="0/1"] = "255,0,255"
    sample1$PL[sample1[,3]=="1/1"] = "255,255,0"
    sample1$PL[sample1[,3]=="./."] = "0,255,255"
    
    colnames(sample1)[16] = paste0("PL_",ids[sample_merge])
    
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
      my_samples$V16 =  NA
      my_samples$V17 =  NA
      
      
      colnames(my_samples)[(ncol(my_samples)-15):ncol(my_samples)] = paste0(c("chr_", # N-1 ncol
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
                                                                              "GT_best_",
                                                                              "PL_"),ids[sample_merge])
      
    }
    print(paste0("The number of variants of the merging is: ",nrow(my_samples)))
    print(paste0("The number of samples of the merging is: ",sample_merge))
    print(paste0("The number of empty samples is: ",empty_samples))
    
  }
  
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
  genos_best = my_samples %>%
    group_by(ID) %>%
    summarize_at(paste0("GT_best_",ids),~paste(unique(na.omit(.)), collapse = ',')) %>%
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
  final_merge = left_join(final_merge,pls) 
  final_merge = left_join(final_merge,ncalls) 
  final_merge = left_join(final_merge,name_calls) 
  final_merge = left_join(final_merge,strategy_var) 
  final_merge = left_join(final_merge,pass) 
  final_merge = left_join(final_merge,pass_num) 
  final_merge = left_join(final_merge,repro) 
  final_merge = left_join(final_merge,genos_best)
  
  final_merge = final_merge %>% arrange(start)
  
  final_merge$start = floor(final_merge$start)
  
  print(dim(final_merge))
  
  final_merge %>% 
    fwrite(paste0("/gpfs/projects/bsc05/jordivalls/GCAT_project_all_samples/merge_all_calling_GCAT//Deleciones/merge_samples/chr_",chr,
                  "/Del_chr_",chr,"_",as.character(start),"_",as.character(end)),row.names=F,sep=" ",na = "NA")
  
  print("Now compress pleaseeeeeeeeeeeeeeee")
  
  system(paste0("gzip /gpfs/projects/bsc05/jordivalls/GCAT_project_all_samples/merge_all_calling_GCAT//Deleciones/merge_samples/chr_",chr,
              "/Del_chr_",chr,"_",as.character(start),"_",as.character(end)))
  
  print("Now generate VCF")
  
  final_merge = fread(paste0("zcat /gpfs/projects/bsc05/jordivalls/GCAT_project_all_samples/merge_all_calling_GCAT//Deleciones/merge_samples/chr_",chr,
                             "/Del_chr_",chr,"_",as.character(start),"_",as.character(end),".gz"))
  
  final_merge = final_merge %>% as.data.frame()
  
  
  # warning, some variants repeated in the VCF
  
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
  
  
  # END #####
  
  final_merge$END = paste0("END=",final_merge$start + final_merge$length,";") 
  
  # SVLEN #####
  
  final_merge$SVLEN = paste0("SVLEN=",final_merge$length,";") 
  
  
  # ERRORBKP ####
  
  final_merge$ERRBKP = paste0("ERRBKP=+-",final_merge$start_error,";")
  
  
  # SVTYPE ####
  
  final_merge$SVTYPE = paste0("SVTYPE=DEL;")
  
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
  
  final_merge$no_var = 785 - (final_merge$pass + final_merge$no_pass)
  
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
  
  final_merge$Geno_00 = paste0("GENO_00=",final_merge$Geno_00,";") 
  final_merge$Geno_01 = paste0("GENO_01=",final_merge$Geno_01,";") 
  final_merge$Geno_11 = paste0("GENO_11=",final_merge$Geno_11,";") 
  
  
  # Caller stats #######
  
  aux = final_merge[,paste0("callers_detected_",ids)]
  
  final_merge$callers_tab = apply(aux,1,function(x){
    
    aux = factor(as.numeric(x))
    
    my_c = as.numeric(levels(aux))
    
    no_calls = which(!1:7 %in% my_c)
    
    levels(aux) = c(levels(aux),no_calls)
    
    return(paste0(table(aux)[order(names(table(aux)))],collapse = ","))
    
  })
  
  final_merge$min_calls = apply(aux,1,function(x){

    return(min(as.numeric(x),na.rm=T))
    
  })
  
  final_merge$max_calls = apply(aux,1,function(x){
    
    return(max(as.numeric(x),na.rm=T))
    
  })
  
  final_merge$median_calls = apply(aux,1,function(x){
    
    return(median(as.numeric(x),na.rm=T))
    
  })
  
  final_merge$callers_tab = paste0("TAB_CALLERS=",final_merge$callers_tab,";") 
  final_merge$min_calls = paste0("MIN_CALLERS=",final_merge$min_calls,";") 
  final_merge$max_calls = paste0("MAX_CALLERS=",final_merge$max_calls,";") 
  final_merge$median_calls = paste0("MEDIAN_CALLERS=",final_merge$median_calls) 
  
  
  
  
  # CREATE VCF
  
  my_vcf = data.frame("CHR"=final_merge$CHR)
  
  colnames(my_vcf) = "#CHROM"
  
  my_vcf$POS = final_merge$start
  my_vcf$ID = paste0(chr,":",my_vcf$POS,":A:<DEL>")
  my_vcf$REF = "A"
  my_vcf$ALT = "<DEL>"
  my_vcf$QUAL = final_merge$qual
  my_vcf$FILTER = ifelse(final_merge$qual>0.5,"PASS","NO_PASS")
  
  my_vcf$INFO = do.call(paste0,list(final_merge$AC,
                                    final_merge$MAC,
                                    final_merge$AF,
                                    final_merge$AN,
                                    final_merge$MAF,
                                    final_merge$END,
                                    final_merge$SVLEN,
                                    final_merge$ERRBKP,
                                    final_merge$SVTYPE,
                                    final_merge$POPVAR,
                                    final_merge$pass,
                                    final_merge$no_pass,
                                    final_merge$Geno_00,
                                    final_merge$Geno_01,
                                    final_merge$Geno_11,
                                    final_merge$callers_tab,
                                    final_merge$min_calls,
                                    final_merge$max_calls,
                                    final_merge$median_calls))
  
  my_vcf$FORMAT = "GT:PL:LM:PVM:SVMETHOD:SVCALL:SVSTRAT:SVREPRO"
  
  
  miss_calculation = final_merge[,paste0(c("GT_"),ids)]
  
  miss_rate = NULL
  
  for(i in 1:nrow(miss_calculation)){
    
    aux = sum(miss_calculation[i,]=="./.")/ncol(miss_calculation)
    
    miss_rate = c(miss_rate,aux)
  }
  
  which_miss_05 = which(miss_rate<0.05)
  
  
  for(sample_vcf in 1:785){
    
    sample1 = final_merge[,paste0(c("GT_",
                                    "PL_",
                                    "PASS_num_",
                                    "PASS_",
                                    "name_callers_",
                                    "callers_detected_",
                                    "strategy_",
                                    "reciprocity_",
                                    "GT_best_"),ids[sample_vcf])]
    
    sample1[,2][sample1[,2]==""] = "0,255,255"
    
    sample1[,3] = round(sample1[,3],4)
    
    sample1[is.na(sample1)] = ""
    
    sample1[,1] = as.character(sample1[,1])
    sample1[,9] = as.character(sample1[,9])
    
    
    # replace "./." for the best genotype for variants with <5% missing
    
    if(sum(sample1[which_miss_05,1]=="./.")>0){
      
      if(sum(sample1[which_miss_05,9]!="")>0){
        
        sample1[which_miss_05,1][which(sample1[which_miss_05,9]!="")] = sample1[which_miss_05,9][which(sample1[which_miss_05,9]!="")]
        
      }
      
    }
    
    
    sample1[which_miss_05,1][sample1[which_miss_05,1]=="./."] = "0/0"
    
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
  
  
  fwrite(my_vcf[which_miss_05,],paste0("/gpfs/projects/bsc05/jordivalls/GCAT_project_all_samples/merge_all_calling_GCAT//Deleciones/merge_samples/chr_",chr,
                       "/Del_chr_",chr,"_",as.character(start),"_",as.character(end),".vcf"),row.names=F,sep="\t")
  
  print("Now compress again pleaseeeeeeeeeeeeeeee")
  
  
  system(paste0("gzip /gpfs/projects/bsc05/jordivalls/GCAT_project_all_samples/merge_all_calling_GCAT//Deleciones/merge_samples/chr_",chr,
              "/Del_chr_",chr,"_",as.character(start),"_",as.character(end),".vcf"))
  
  fwrite(my_vcf[-which_miss_05,],paste0("/gpfs/projects/bsc05/jordivalls/GCAT_project_all_samples/merge_all_calling_GCAT//Deleciones/merge_samples/more_5_percent_missings/",
                                        "/Del_large_size_missings_chr_",chr,"_",as.character(start),"_",as.character(end),".vcf"),row.names=F,sep="\t")
  
  print("Now compress again pleaseeeeeeeeeeeeeeee")
  
  
  system(paste0("gzip /gpfs/projects/bsc05/jordivalls/GCAT_project_all_samples/merge_all_calling_GCAT//Deleciones/merge_samples/more_5_percent_missings/",
                "/Del_large_size_missings_chr_",chr,"_",as.character(start),"_",as.character(end),".vcf"))
  
  
}



