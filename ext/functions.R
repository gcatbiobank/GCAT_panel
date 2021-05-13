##########################################################################################
##  Outline                                                                           ####
##                                                                                      ##
## 1 - Deletions and duplications (one breakpoint and length resolution)                ##
## 2 - Insertions (one breakpoint and no length resolution)                             ##
## 3 - Inversions (two breakpoints at same chromosome and length resolution)            ##
## 4 - Translocations (two breakpoints in both chromosomes and no length resolution)    ##
##                                                                                      ##
##########################################################################################


### merge callers deletions and duplications #####

merge_callers_deletions_duplications = function(call1,call2,callers_ref,caller_to_merge,
                                                golden=FALSE,repro,svtype){
  
  # call2: son las variantes de referencia para la intersección
  # call1: al que vas a buscar las posiciones
  # callers_ref: los callers de referencia que hay en call2
  # caller_to_merge: el caller que vamos a hacer el merge
  # golden: si el merge es con la golden
  
  #call2 = lumpy
  #call1 = delly
  #callers_ref = c("lumpy")
  #caller_to_merge = c("delly")
  
  n_callers = length(callers_ref)
  
  for(i in 1:nrow(call1)){
    
    k = 1
    
    while(k<n_callers+1){
      
      aux = call1[i,c("chr",paste0(c("lower_","upper_","length_"),caller_to_merge)),with=F]
      
      call1_interval = (as.numeric(aux[1,2])):(as.numeric(aux[1,3]))
      call1_chr = aux$chr
      
      call2_find = call2[call2$chr==call1_chr,
                         c("ID",paste0(c("lower_","upper_","length_"),callers_ref[k])),
                         with=F]
      
      if(sum(is.na(call2_find[,2]))>0){
        call2_find = call2_find[which(!is.na(call2_find[,2])),]} # remove NA of some variants when the merge between more than 1 caller is done
      
      if(nrow(call2_find)>0){
        
        for(j in 1:nrow(call2_find)){
          
          my_find = findInterval(call1_interval,
                                 as.numeric(call2_find[j,2:3]))
          
          reprocicity = as.numeric(aux[1,4])/as.numeric(call2_find[j,4])
          
          reprocicity = ifelse(reprocicity>1,1/reprocicity,reprocicity)
          
          if(sum(my_find==1)>0 & reprocicity>repro){
            
            call1[i,]$ID = as.character(call2_find[j,]$ID)
            
            k = n_callers + 1
            
          }
        }
      }
      
      k = k + 1
      
    }
    
    print(i)
    
    
  }
  
  
  
  if(golden==FALSE){
    
    new_sv = nrow(call1[call1$ID=="none",])
    
    call1[call1$ID=="none",]$ID = paste0(svtype,"_",(nrow(call2)+1):(nrow(call2)+new_sv))
    
    
    my_data = full_join(call2[,c("ID","chr",
                                 as.vector(outer(c("chr_pos_","start_","length_",
                                                   "GT_","lower_","upper_"),callers_ref,paste0))),with=F],
                        call1[,c("ID","chr",
                                 paste0(c("chr_pos_","start_","length_",
                                          "GT_","lower_","upper_"),caller_to_merge)),with=F])
  }
  if(golden==TRUE){
    
    my_data = left_join(call2[,c("ID","chr",
                                 as.vector(outer(c("chr_pos_","start_","length_",
                                                   "GT_","lower_","upper_"),callers_ref,paste0))),with=F],
                        call1[,c("ID","chr",
                                 paste0(c("chr_pos_","start_","length_",
                                          "GT_","lower_","upper_"),caller_to_merge)),with=F])
  }
  return(as.data.table(my_data))
  
}


# sensitivity precision geno length error Deletions and Duplications ######

sensitivity_precision_geno_error_deletions_duplications <- function(call,golden,type,call_name){
  
  golden$ID = "none"
  
  call$ID = paste0(type,1:nrow(call))
  
  call = merge_callers_deletions_duplications(call1=golden,call2=call,callers_ref = call_name,
                                              caller_to_merge = c("golden"),golden = TRUE, repro=0.8,svtype=type)
  all = call
  
  call$PASS = "YES"
  call$PASS[is.na(call$length_golden)] = "NO"
  table(call$PASS)
  
  merge_table = call
  
  sens = (table(call$PASS)[2]/nrow(golden))*100
  
  spec = (table(call$PASS)[2]/nrow(call))*100
  
  IC_lower_sens = round(as.numeric(sens) - 1.96*sqrt(as.numeric(sens)*(100-as.numeric(sens))/nrow(golden)),digits=2)
  IC_upper_sens = round(as.numeric(sens) + 1.96*sqrt(as.numeric(sens)*(100-as.numeric(sens))/nrow(golden)),digits=2)
  IC_sens = paste0("[",IC_lower_sens,",",IC_upper_sens,"]")
  
  IC_lower_spec = round(as.numeric(spec) - 1.96*sqrt(as.numeric(spec)*(100-as.numeric(spec))/nrow(call)),digits=2)
  IC_upper_spec = round(as.numeric(spec) + 1.96*sqrt(as.numeric(spec)*(100-as.numeric(spec))/nrow(call)),digits=2)
  IC_spec = paste0("[",IC_lower_spec,",",IC_upper_spec,"]")
  
  
  true_positive = table(call$PASS)[2]
  false_positive = nrow(call) - table(call$PASS)[2]
  
  call = call[!is.na(call$GT_golden),]
  
  geno1 = as.character(t(call[,paste0("GT_",call_name),with=F]))
  geno2 = as.character(t(call$GT_golden))
  
  
  if(sum(geno1 %in% "9/9")>0){  # replace non-genotype (Manta)
    which_99 = which(geno1 %in% "9/9")
    geno1 = geno1[-which_99]
    geno2 = geno2[-which_99]}
  
  g.error = 100-sum(geno1==geno2)/length(geno1)*100
  
  
  length1 = as.numeric(t(call[,paste0("length_",call_name),with=F]))
  length2 = as.numeric(t(call$length_golden))
  
  l.error = median(length1-length2,na.rm = T)
  
  homo_golden = nrow(call[call$GT_golden=="1/1",])
  het_golden = nrow(call[call$GT_golden=="0/1",])
  
  
  homo_golden_het_caller = nrow(call[call$GT_golden=="1/1" & geno1=="0/1",])
  homo_golden_hom_caller = nrow(call[call$GT_golden=="1/1" & geno1=="1/1",])
  het_golden_het_caller = nrow(call[call$GT_golden=="0/1" & geno1=="0/1",])
  het_golden_hom_caller = nrow(call[call$GT_golden=="0/1" & geno1=="1/1",])
  
  g_error_het = 100-nrow(call[call$GT_golden=="0/1" & geno1=="0/1",])/het_golden*100
  g_error_hom = 100-nrow(call[call$GT_golden=="1/1" & geno1=="1/1",])/homo_golden*100
  
  return(list(sensitivity = sens, precision = spec,
              true_positive = true_positive,
              false_positive = false_positive,
              f1_score = ((2*sens*spec)/(sens+spec))/100,
              g_error = g.error,
              n_het = het_golden,
              n_homo = homo_golden,
              g_error_het = g_error_het,
              g_error_hom = g_error_hom,
              l_error = l.error,IC_sens = IC_sens,IC_spec = IC_spec,
              all = all))
  
}


### merge callers insertions #####

merge_callers_insertions = function(call1,call2,callers_ref,caller_to_merge,
                                    golden=FALSE){
  
  # call2: son las variantes de referencia para la intersección
  # call1: al que vas a buscar las posiciones
  # callers_ref: los callers de referencia que hay en call2
  # caller_to_merge: el caller que vamos a hacer el merge
  # golden: si el merge es con la golden
  
  #call2 = lumpy
  #call1 = delly
  #callers_ref = c("lumpy")
  #caller_to_merge = c("delly")
  
  n_callers = length(callers_ref)
  
  for(i in 1:nrow(call1)){
    
    k = 1
    
    while(k<n_callers+1){
      
      aux = call1[i,c("chr",paste0(c("lower_","upper_","length_"),caller_to_merge)),with=F]
      
      call1_interval = (as.numeric(aux[1,2])):(as.numeric(aux[1,3]))
      call1_chr = aux$chr
      
      call2_find = call2[call2$chr==call1_chr,
                         c("ID",paste0(c("lower_","upper_","length_"),callers_ref[k])),
                         with=F]
      
      if(sum(is.na(call2_find[,2]))>0){
        call2_find = call2_find[which(!is.na(call2_find[,2])),]} # remove NA of some variants when the merge between more than 1 caller is done
      
      if(nrow(call2_find)>0){
        
        for(j in 1:nrow(call2_find)){
          
          my_find = findInterval(call1_interval,
                                 as.numeric(call2_find[j,2:3]))
          
          if(sum(my_find==1)>0){
            
            call1[i,]$ID = as.character(call2_find[j,]$ID)
            
            k = n_callers + 1
            
          }
        }
      }
      
      k = k + 1
      
    }
    
    print(i)
    
    
  }
  
  
  
  if(golden==FALSE){
    
    new_insertions = nrow(call1[call1$ID=="none",])
    
    call1[call1$ID=="none",]$ID = paste0("insertion_",(nrow(call2)+1):(nrow(call2)+new_insertions))
    
    
    my_data = full_join(call2[,c("ID","chr",
                                 as.vector(outer(c("chr_pos_","start_","length_",
                                                   "GT_","lower_","upper_"),callers_ref,paste0))),with=F],
                        call1[,c("ID","chr",
                                 paste0(c("chr_pos_","start_","length_",
                                          "GT_","lower_","upper_"),caller_to_merge)),with=F])
  }
  if(golden==TRUE){
    
    my_data = left_join(call2[,c("ID","chr",
                                 as.vector(outer(c("chr_pos_","start_","length_",
                                                   "GT_","lower_","upper_"),callers_ref,paste0))),with=F],
                        call1[,c("ID","chr",
                                 paste0(c("chr_pos_","start_","length_",
                                          "GT_","lower_","upper_"),caller_to_merge)),with=F])
  }
  return(as.data.table(my_data))
  
}


# sensitivity precision geno length error Insertions ######

sensitivity_precision_geno_error_insertions <- function(call,golden,type,call_name){
  
  golden$ID = "none"
  
  call$ID = paste0(type,1:nrow(call))
  
  call = merge_callers_insertions(call1=golden,call2=call,callers_ref = call_name,
                                  caller_to_merge = c("golden"),golden = TRUE)
  all = call
  
  call$PASS = "YES"
  call$PASS[is.na(call$length_golden)] = "NO"
  table(call$PASS)
  
  merge_table = call
  
  sens = (table(call$PASS)[2]/nrow(golden))*100
  
  spec = (table(call$PASS)[2]/nrow(call))*100

  IC_lower_sens = round(as.numeric(sens) - 1.96*sqrt(as.numeric(sens)*(100-as.numeric(sens))/nrow(golden)),digits=2)
  IC_upper_sens = round(as.numeric(sens) + 1.96*sqrt(as.numeric(sens)*(100-as.numeric(sens))/nrow(golden)),digits=2)
  IC_sens = paste0("[",IC_lower_sens,",",IC_upper_sens,"]")
  
  IC_lower_spec = round(as.numeric(spec) - 1.96*sqrt(as.numeric(spec)*(100-as.numeric(spec))/nrow(call)),digits=2)
  IC_upper_spec = round(as.numeric(spec) + 1.96*sqrt(as.numeric(spec)*(100-as.numeric(spec))/nrow(call)),digits=2)
  IC_spec = paste0("[",IC_lower_spec,",",IC_upper_spec,"]")
  
  
  true_positive = table(call$PASS)[2]
  false_positive = nrow(call) - table(call$PASS)[2]
  
  call = call[!is.na(call$GT_golden),]
  
  geno1 = as.character(t(call[,paste0("GT_",call_name),with=F]))
  geno2 = as.character(t(call$GT_golden))
  
  
  if(sum(geno1 %in% "9/9")>0){  # replace non-genotype (Manta)
    which_99 = which(geno1 %in% "9/9")
    geno1 = geno1[-which_99]
    geno2 = geno2[-which_99]}
  
  g.error = 100-sum(geno1==geno2)/length(geno1)*100
  
  
  length1 = as.numeric(t(call[,paste0("length_",call_name),with=F]))
  length2 = as.numeric(t(call$length_golden))
  
  homo_golden = nrow(call[call$GT_golden=="1/1",])
  het_golden = nrow(call[call$GT_golden=="0/1",])
  
  
  homo_golden_het_caller = nrow(call[call$GT_golden=="1/1" & geno1=="0/1",])
  homo_golden_hom_caller = nrow(call[call$GT_golden=="1/1" & geno1=="1/1",])
  het_golden_het_caller = nrow(call[call$GT_golden=="0/1" & geno1=="0/1",])
  het_golden_hom_caller = nrow(call[call$GT_golden=="0/1" & geno1=="1/1",])
  
  g_error_het = 100-nrow(call[call$GT_golden=="0/1" & geno1=="0/1",])/het_golden*100
  g_error_hom = 100-nrow(call[call$GT_golden=="1/1" & geno1=="1/1",])/homo_golden*100
  
  return(list(sensitivity = sens, precision = spec,
              true_positive = true_positive,
              false_positive = false_positive,
              f1_score = ((2*sens*spec)/(sens+spec))/100,
              g_error = g.error,
              n_het = het_golden,
              n_homo = homo_golden,
              g_error_het = g_error_het,
              g_error_hom = g_error_hom,
              IC_sens = IC_sens,
              IC_spec = IC_spec,
              all = all))
  
}


## unify breakpoints Inversions ######


unify_breakpoints = function(caller,name){
  
  type = "inversion"
  
  caller$ID = "none"
  
  caller$ID[1] = paste0(type,"_",1)
  
  for(i in 2:nrow(caller)){
    
    aux = caller[(i-1):i,c(as.vector(outer(c("start_","end_","length_"),name,paste0)),"chr"),with=F]
    
    reprocicity = as.numeric(aux[1,3])/as.numeric(aux[2,3])
    
    reprocicity = ifelse(reprocicity>1,as.numeric(aux[2,3])/as.numeric(aux[1,3]),reprocicity)
    
    if(caller[i,1]==caller[i-1,1] & reprocicity>0.80 & (abs(caller[i,2]-caller[i-1,2])<=50 | # maxim error que tolerem
                                                        abs(caller[i,3]-caller[i-1,3])<=50) &
       caller[i,5]==caller[i-1,5]){
      
      caller$ID[i] = caller$ID[i-1]
      
    }
    else{caller$ID[i] = paste0(type,"_",i)}
    
  }
  
  caller_ok = data.frame(ID = unique(caller$ID))
  
  caller_ok$chr = 0
  caller_ok$start = 0
  caller_ok$end = 0
  
  caller_ok$GT = 0
  
  
  for(i in 1:nrow(caller_ok)){
    
    aux = caller %>% filter(ID==caller_ok$ID[i])
    
    caller_ok$chr[i] = unique(aux[,1])
    caller_ok$start[i] = min(aux[,2])
    caller_ok$end[i] = max(aux[,3])
    
    caller_ok$GT[i] = paste0(unique(aux[,5]),collapse = "")
    
    
  }
  
  return(caller_ok)
  
}



## merge callers 2 breakpoints same chromosome (inversions) #####


merge_callers_inv = function(call1,call2,callers_ref,caller_to_merge,
                             golden=FALSE,type,merge=FALSE){
  
  # call2: son las variantes de referencia para la intersección
  # call1: al que vas a buscar las posiciones
  # callers_ref: los callers de referencia que hay en call2
  # caller_to_merge: el caller que vamos a hacer el merge
  # golden: si el merge es con la golden
  
  #call1 = call
  #call2 = golden
  #caller_to_merge = call_name
  #callers_ref = c("golden")
  #call2$ID = paste0("inversion",1:nrow(call2))
  #call1$ID = "none"
  #type = "inversion"
  #i=2
  
  n_callers = length(callers_ref)
  
  
  for(i in 1:nrow(call1)){
    
    k = 1
    
    while(k<n_callers+1){
      
      aux = call1[i,c("chr",paste0(c("lower_","upper_"),caller_to_merge,c("_bp1","_bp1")),paste0(c("length_"),caller_to_merge)),with=F]
      aux2 = call1[i,c("chr",paste0(c("lower_","upper_"),caller_to_merge,c("_bp2","_bp2"))),with=F]
      
      call1_interval = (as.numeric(aux[1,2])):(as.numeric(aux[1,3]))
      call1_chr = aux$chr
      
      call1_interval_bp = (as.numeric(aux2[1,2])):(as.numeric(aux2[1,3]))
      call1_chr_bp = aux$chr
      
      call2_find = call2[call2$chr==call1_chr,c("ID",paste0(c("lower_","upper_"),callers_ref[k],c("_bp1","_bp1")),paste0(c("length_"),callers_ref[k])),with=F]
      call2_find_bp = call2[call2$chr==call1_chr,c("ID",paste0(c("lower_","upper_"),callers_ref[k],c("_bp2","_bp2"))),with=F]
      
      if(sum(is.na(call2_find[,2]))>0){call2_find = call2_find[which(!is.na(call2_find[,2])),]} # remove NA of some variants
      if(sum(is.na(call2_find_bp[,2]))>0){call2_find_bp = call2_find_bp[which(!is.na(call2_find_bp[,2])),]} # remove NA of some variants
      
      if(nrow(call2_find)>0){
        
        for(j in 1:nrow(call2_find)){
          
          #j=3
          
          if(golden==TRUE){
            
            call_interval = as.numeric(call2_find[j,2]):as.numeric(call2_find[j,3])
            call_interval_bp = as.numeric(call2_find_bp[j,2]):as.numeric(call2_find_bp[j,3])
            
            golden_start = call1_interval
            golden_end = call1_interval_bp
            
            my_find = sum(golden_start %in% call_interval)
            my_find_bp = sum(golden_end %in% call_interval_bp)
            
            
            if((my_find==1 | my_find_bp==1)){ # the first or second bp
              
              call1[i,]$ID = as.character(call2_find[j,]$ID)
              
              k = n_callers + 1
              
            }
            
          }
          
          if(golden==FALSE){
            
            my_find = sum(call1_interval %in% as.numeric(call2_find[j,2]):as.numeric(call2_find[j,3]))
            my_find_bp = sum(call1_interval_bp %in% as.numeric(call2_find_bp[j,2]):as.numeric(call2_find_bp[j,3]))
            
            reprocicity = as.numeric(aux[1,4])/as.numeric(call2_find[j,4])
            
            reprocicity = ifelse(reprocicity>1,as.numeric(call2_find[j,4])/as.numeric(aux[1,4]),reprocicity)
            
            if((my_find>=1 | my_find_bp>=1) & reprocicity>0.8){
              
              call1[i,]$ID = as.character(call2_find[j,]$ID)
              
              k = n_callers + 1
              
            }
            
          }
          
          
          
        }
      }
      
      k = k + 1
      
    }
    
    print(i)
    
    
  }
  
  # IDs unic per posicions properes (loop by chromosome)
  
  new_variants = call1[call1$ID=="none",]
  
  if(nrow(new_variants)==1){
    id = nrow(call2)+1
    new_variants$ID[1] = paste0(type,"_",id)
    call1[call1$ID=="none",]$ID = new_variants$ID}
  
  if(nrow(new_variants)>1){
    
    id = nrow(call2)+1
    
    new_variants$ID[1] = paste0(type,"_",id)
    
    id = id + 1
    
    for(i in 2:nrow(new_variants)){
      
      aux = new_variants[(i-1):i,c(as.vector(outer(c("start_","end_","length_"),caller_to_merge,paste0)),"chr"),with=F]
      
      reprocicity = as.numeric(aux[1,3])/as.numeric(aux[2,3])
      
      reprocicity = ifelse(reprocicity>1,as.numeric(aux[2,3])/as.numeric(aux[1,3]),reprocicity)
      
      if((aux$chr[1]==aux$chr[2]) & reprocicity>0.80 & (abs(aux[1,1]-aux[2,1])<=50 | (abs(aux[1,2]-aux[2,2])<=50))){ # remove also rows with close bp
        new_variants$ID[i] = new_variants$ID[i-1]}
      
      else{
        
        new_variants$ID[i] = paste0(type,"_",id)
        
        id = id + 1
        
      }
      
      
    }
    
  }
  
  
  call1[call1$ID=="none",]$ID = new_variants$ID
  
  bp_names = as.vector(outer(c("lower_","upper_"),callers_ref,paste0))
  
  bp_names1 = as.vector(outer(bp_names,c("_bp1"),paste0))
  bp_names2 = as.vector(outer(bp_names,c("_bp2"),paste0))
  
  
  if(golden==FALSE){
    
    my_data = full_join(call2[,c("ID","chr",
                                 as.vector(outer(c("start_","end_","length_",
                                                   "GT_"),callers_ref,paste0)),
                                 bp_names1,bp_names2),with=F],
                        call1[,c("ID","chr",
                                 paste0(c("start_","end_","length_",
                                          "GT_"),caller_to_merge),
                                 paste0(c("lower_","upper_"),caller_to_merge,c("_bp1","_bp1")),
                                 paste0(c("lower_","upper_"),caller_to_merge,c("_bp2","_bp2"))),with=F])
  }
  if(golden==TRUE & merge==FALSE){
    
    
    my_data = right_join(call2[,c("ID","chr",
                                  as.vector(outer(c("start_","end_","length_",
                                                    "GT_"),callers_ref,paste0)),
                                  bp_names1,bp_names2),with=F],
                         call1[,c("ID","chr",
                                  paste0(c("start_","end_","length_",
                                           "GT_"),caller_to_merge),
                                  paste0(c("lower_","upper_"),caller_to_merge,c("_bp1","_bp1")),
                                  paste0(c("lower_","upper_"),caller_to_merge,c("_bp2","_bp2"))),with=F])
  }
  if(golden==TRUE & merge==TRUE){
    
    
    my_data = left_join(call2[,c("ID","chr",
                                 as.vector(outer(c("start_","end_","length_",
                                                   "GT_"),callers_ref,paste0)),
                                 bp_names1,bp_names2),with=F],
                        call1[,c("ID","chr",
                                 paste0(c("start_","end_","length_",
                                          "GT_"),caller_to_merge),
                                 paste0(c("lower_","upper_"),caller_to_merge,c("_bp1","_bp1")),
                                 paste0(c("lower_","upper_"),caller_to_merge,c("_bp2","_bp2"))),with=F])
  }
  
  
  return(as.data.table(my_data))
  
}


## sensitivity-precision 2 breakpoints same chromosome (inversions) #####


sensitivity_precision_geno_error_inv <- function(call,golden,type,call_name){
  
  #call = lumpy
  #golden = golden
  #type = "inversion"
  #call_name = "lumpy"
  
  golden$ID = paste0(type,1:nrow(golden))
  
  call$ID = "none"
  call = merge_callers_inv(call1=call,call2=golden,callers_ref =  c("golden"),
                           caller_to_merge = call_name,golden = TRUE,type)
  
  all = call
  
  # sensitivity precision
  
  if(sum(duplicated(call$ID))>1){
    
    call_sens = call[-which(duplicated(call$ID)),]}
  
  if(sum(duplicated(call$ID))==0){
    
    call_sens = call}
  
  
  call_sens$PASS = "YES"
  call_sens$PASS[is.na(call_sens$length_golden)] = "NO"
  table(call_sens$PASS)
  
  sens = (table(call_sens$PASS)[2]/nrow(golden))*100
  
  spec = (table(call_sens$PASS)[2]/nrow(call_sens))*100
  
  IC_lower_sens = round(as.numeric(sens) - 1.96*sqrt(as.numeric(sens)*(100-as.numeric(sens))/nrow(golden)),digits=2)
  IC_upper_sens = round(as.numeric(sens) + 1.96*sqrt(as.numeric(sens)*(100-as.numeric(sens))/nrow(golden)),digits=2)
  IC_sens = paste0("[",IC_lower_sens,",",IC_upper_sens,"]")
  
  IC_lower_spec = round(as.numeric(spec) - 1.96*sqrt(as.numeric(spec)*(100-as.numeric(spec))/nrow(call)),digits=2)
  IC_upper_spec = round(as.numeric(spec) + 1.96*sqrt(as.numeric(spec)*(100-as.numeric(spec))/nrow(call)),digits=2)
  IC_spec = paste0("[",IC_lower_spec,",",IC_upper_spec,"]")
  
  true_positive = table(call_sens$PASS)[2]
  false_positive = nrow(call_sens) - table(call_sens$PASS)[2]
  
  
  # genotype error
  
  call = call[!is.na(call$GT_golden),c("ID",paste0("GT_",call_name),"GT_golden"),with=F]
  call = as.data.frame(unique(call))
  
  if(sum(duplicated(call$ID))>0){
    
    var_dup = call[which(duplicated(call$ID)),]$ID
    
    call[call$ID %in% var_dup,paste0("GT_",call_name)] = "./."
    call = unique(call)
    
  }
  
  geno1 = as.character(t(call[,paste0("GT_",call_name)]))
  geno2 = as.character(t(call$GT_golden))
  
  
  if(sum(geno1 %in% "9/9")>0){  
    which_99 = which(geno1 %in% "9/9")
    geno1 = geno1[-which_99]
    geno2 = geno2[-which_99]}
  
  
  g.error = 100-sum(geno1==geno2)/length(geno1)*100
  
  homo_golden = nrow(call[call$GT_golden=="1/1",])
  het_golden = nrow(call[call$GT_golden=="0/1",])
  
  
  homo_golden_het_caller = nrow(call[call$GT_golden=="1/1" & geno1=="0/1",])
  homo_golden_hom_caller = nrow(call[call$GT_golden=="1/1" & geno1=="1/1",])
  het_golden_het_caller = nrow(call[call$GT_golden=="0/1" & geno1=="0/1",])
  het_golden_hom_caller = nrow(call[call$GT_golden=="0/1" & geno1=="1/1",])
  
  g_error_het = 100-nrow(call[call$GT_golden=="0/1" & geno1=="0/1",])/het_golden*100
  g_error_hom = 100-nrow(call[call$GT_golden=="1/1" & geno1=="1/1",])/homo_golden*100
  
  return(list(sensitivity = sens, precision = spec,
              true_positive = true_positive,
              false_positive = false_positive,
              f1_score = ((2*sens*spec)/(sens+spec))/100,
              g_error = g.error,
              n_het = het_golden,
              n_homo = homo_golden,
              g_error_het = g_error_het,
              g_error_hom = g_error_hom,
              IC_sens = IC_sens,
              IC_spec = IC_spec,
              all = all))
}







#### consensus length ######

#merge = my_data2
#callers = c("lumpy","delly","pindel","whamg",
#           "manta","cnvnator","svaba")


consensus_length <- function(merge,callers){
  
  n_callers = length(callers)
  
  data_length = merge[,paste0("length_",callers),with=F]
  
  data_length$consensus = NA
  
  data_length = as.matrix(data_length)
  
  for(i in 1:nrow(data_length)){
    
    my_tab = as.numeric(data_length[i,1:n_callers])
    
    
    data_length[i,n_callers+1] = floor(median(my_tab,na.rm = T))
    
  }
  
  data_length = as.data.frame(data_length)
  
  out = as.numeric(as.character.numeric_version(data_length$consensus))
  
  return(out)
  
  
}


#### consensus start Deletions #####


consensus_start <- function(merge,callers){
  
  n_callers = length(callers)
  
  data_length = merge[,paste0("start_",callers),with=F]
  
  data_length$consensus = NA
  
  data_length = as.matrix(data_length)
  
  for(i in 1:nrow(data_length)){
    
    my_tab = as.numeric(data_length[i,1:n_callers])
    
    which_na = which(is.na(my_tab))
    
    if(sum(callers[-which_na] %in% "pindel")==1){
      
      data_length[i,n_callers+1] = data_length[i,"start_pindel"]
      
    }
    
    if(sum(callers[-which_na] %in% "pindel")==0 & sum(callers[-which_na] %in% "whamg")==1){
      
      data_length[i,n_callers+1] = data_length[i,"start_whamg"]
      
    }
    if(sum(callers[-which_na] %in% "pindel")==0 & sum(callers[-which_na] %in% "whamg")==0 & 
       sum(callers[-which_na] %in% "delly")==1){
      
      data_length[i,n_callers+1] = data_length[i,"start_delly"]
      
    }
    if(sum(callers[-which_na] %in% "pindel")==0 & sum(callers[-which_na] %in% "whamg")==0 & 
       sum(callers[-which_na] %in% "delly")==0 & sum(callers[-which_na] %in% "manta")==1){
      
      data_length[i,n_callers+1] = data_length[i,"start_manta"]
      
    }
    
    if(sum(callers[-which_na] %in% "pindel")==0 & sum(callers[-which_na] %in% "whamg")==0 & 
       sum(callers[-which_na] %in% "delly")==0 & sum(callers[-which_na] %in% "manta")==0){
      
      data_length[i,n_callers+1] = floor(median(my_tab,na.rm = T))
      
    }
    
  }
  
  data_length = as.data.frame(data_length)
  
  out = as.numeric(as.character.numeric_version(data_length$consensus))
  
  return(out)
}



#### consensus start Duplications #####


consensus_start_duplications <- function(merge,callers){
  
  n_callers = length(callers)
  
  data_length = merge[,paste0("start_",callers),with=F]
  
  data_length$consensus = NA
  
  data_length = as.matrix(data_length)
  
  for(i in 1:nrow(data_length)){
    
    my_tab = as.numeric(data_length[i,1:n_callers])
    
    which_na = which(is.na(my_tab))
    
    if(!is.na(merge$start_svaba[i])){
      
      data_length[i,n_callers+1] = merge$start_svaba[i]
      
    }
    
    if(is.na(merge$start_svaba[i]) & !is.na(merge$start_pindel[i])){
      
      data_length[i,n_callers+1] = merge$start_pindel[i]
      
    }
    
    if(is.na(merge$start_svaba[i]) & is.na(merge$start_pindel[i]) &
       !is.na(merge$start_delly[i])){
      
      data_length[i,n_callers+1] = merge$start_delly[i]
      
    }
    
    if(is.na(merge$start_svaba[i]) & is.na(merge$start_pindel[i]) &
       is.na(merge$start_delly[i]) & !is.na(merge$start_whamg[i])){      
      data_length[i,n_callers+1] = merge$start_whamg[i]
      
    }
    
    
    if(is.na(merge$start_svaba[i]) & is.na(merge$start_pindel[i]) &
       is.na(merge$start_delly[i]) & is.na(merge$start_whamg[i])){        
      data_length[i,n_callers+1] = floor(median(my_tab,na.rm = T))
      
    }
    
  }
  
  data_length = as.data.frame(data_length)
  
  out = as.numeric(as.character.numeric_version(data_length$consensus))
  
  return(out)
}



#### consensus start Insertions #####


consensus_start_insertions<- function(merge,callers){
  
  n_callers = length(callers)
  
  data_length = merge[,paste0("start_",callers),with=F]
  
  data_length$consensus = NA
  
  data_length = as.matrix(data_length)
  
  for(i in 1:nrow(data_length)){
    
    my_tab = as.numeric(data_length[i,1:n_callers])
    
    if(!is.na(merge$start_pindel[i])){
      
      data_length[i,n_callers+1] = merge$start_pindel[i]
      
    }
    
    if(is.na(merge$start_pindel[i]) & !is.na(merge$start_delly[i])){
      
      data_length[i,n_callers+1] = merge$start_delly[i]
      
    }
    
    if(is.na(merge$start_pindel[i]) & is.na(merge$start_delly[i]) & 
       !is.na(merge$start_strelka[i])){
      
      data_length[i,n_callers+1] = merge$start_strelka[i]
      
    }
    
    if(is.na(merge$start_pindel[i]) & is.na(merge$start_delly[i]) & 
       is.na(merge$start_strelka[i]) & !is.na(merge$start_svaba[i]) ){
      
      data_length[i,n_callers+1] = merge$start_svaba[i]
      
    }
    
    if(is.na(merge$start_pindel[i]) & is.na(merge$start_delly[i]) & 
       is.na(merge$start_strelka[i]) & is.na(merge$start_svaba[i]) & !is.na(merge$start_manta[i])){
      
      data_length[i,n_callers+1] = merge$start_manta[i]
      
    }
    
    if(is.na(merge$start_pindel[i]) & is.na(merge$start_delly[i]) & 
       is.na(merge$start_strelka[i]) & is.na(merge$start_svaba[i]) & is.na(merge$start_manta[i])){
      
      data_length[i,n_callers+1] = floor(median(my_tab,na.rm = T))
      
    }
    
  }
  
  data_length = as.data.frame(data_length)
  
  out = as.numeric(as.character.numeric_version(data_length$consensus))
  
  return(out)
}


#### consensus start Inversions #####


consensus_start_inversions<- function(merge,callers){
  
  n_callers = length(callers)
  
  data_length = merge[,paste0("start_",callers),with=F]
  
  data_length$consensus = NA
  
  data_length = as.matrix(data_length)
  
  for(i in 1:nrow(data_length)){
    
    my_tab = as.numeric(data_length[i,1:n_callers])
    
    if(!is.na(merge$start_lumpy[i])){
      
      data_length[i,n_callers+1] = merge$start_lumpy[i]
      
    }
    
    if(is.na(merge$start_lumpy[i]) & !is.na(merge$start_pindel[i])){
      
      data_length[i,n_callers+1] = merge$start_pindel[i]
      
    }
    
    if(is.na(merge$start_pindel[i]) & is.na(merge$start_pindel[i]) & 
       !is.na(merge$start_delly[i])){
      
      data_length[i,n_callers+1] = merge$start_delly[i]
      
    }
    
    
    if(is.na(merge$start_pindel[i]) & is.na(merge$start_pindel[i]) & 
       is.na(merge$start_delly[i])){
      
      data_length[i,n_callers+1] = floor(median(my_tab,na.rm = T))
      
    }
    
  }
  
  data_length = as.data.frame(data_length)
  
  out = as.numeric(as.character.numeric_version(data_length$consensus))
  
  return(out)
}



#### consensus start Translocations #####


consensus_start_1_translocations <- function(merge,callers){
  
  n_callers = length(callers)
  
  data_length = merge[,paste0("start_1_",callers),with=F]
  
  data_length$consensus = NA
  
  data_length = as.matrix(data_length)
  
  if(!callers %in% "manta"){merge$start_1_manta = NA}
  
  for(i in 1:nrow(data_length)){
    
    my_tab = as.numeric(data_length[i,1:n_callers])
    
    if(!is.na(merge$start_1_manta[i])){
      
      data_length[i,n_callers+1] = merge$start_1_manta[i]
      
    }
    if(is.na(merge$start_1_manta[i])){
      
      data_length[i,n_callers+1] = floor(median(my_tab,na.rm = T))
    }
  }
  
  data_length = as.data.frame(data_length)
  
  out = as.numeric(as.character.numeric_version(data_length$consensus))
  
  return(out)
}

consensus_start_2_translocations <- function(merge,callers){
  
  n_callers = length(callers)
  
  data_length = merge[,paste0("start_2_",callers),with=F]
  
  data_length$consensus = NA
  
  data_length = as.matrix(data_length)
  
  if(!callers %in% "manta"){merge$start_2_manta = NA}
  
  for(i in 1:nrow(data_length)){
    
    my_tab = as.numeric(data_length[i,1:n_callers])
    
    if(!is.na(merge$start_2_manta[i])){
      
      data_length[i,n_callers+1] = merge$start_2_manta[i]
      
    }
    if(is.na(merge$start_2_manta[i])){
      
      data_length[i,n_callers+1] = floor(median(my_tab,na.rm = T))
    }
    
  }
  
  data_length = as.data.frame(data_length)
  
  out = as.numeric(as.character.numeric_version(data_length$consensus))
  
  return(out)
}

consensus_end_1_translocations <- function(merge,callers){
  
  n_callers = length(callers)
  
  data_length = merge[,paste0("end_1_",callers),with=F]
  
  data_length$consensus = NA
  
  data_length = as.matrix(data_length)
  
  if(!callers %in% "manta"){merge$end_1_manta = NA}
  
  for(i in 1:nrow(data_length)){
    
    my_tab = as.numeric(data_length[i,1:n_callers])
    
    if(!is.na(merge$end_1_manta[i])){
      
      data_length[i,n_callers+1] = merge$end_1_manta[i]
      
    }
    if(is.na(merge$end_1_manta[i])){
      
      data_length[i,n_callers+1] = floor(median(my_tab,na.rm = T))
    }
    
  }
  
  data_length = as.data.frame(data_length)
  
  out = as.numeric(as.character.numeric_version(data_length$consensus))
  
  return(out)
}

consensus_end_2_translocations <- function(merge,callers){
  
  n_callers = length(callers)
  
  data_length = merge[,paste0("end_2_",callers),with=F]
  
  data_length$consensus = NA
  
  data_length = as.matrix(data_length)
  
  if(!callers %in% "manta"){merge$end_2_manta = NA}
  
  for(i in 1:nrow(data_length)){
    
    my_tab = as.numeric(data_length[i,1:n_callers])
    
    if(!is.na(merge$end_2_manta[i])){
      
      data_length[i,n_callers+1] = merge$end_2_manta[i]
      
    }
    if(is.na(merge$end_2_manta[i])){
      
      data_length[i,n_callers+1] = floor(median(my_tab,na.rm = T))
    }
    
  }
  
  data_length = as.data.frame(data_length)
  
  out = as.numeric(as.character.numeric_version(data_length$consensus))
  
  return(out)
}





#### consensus genotype ######

consensus_genotype <- function(merge,callers){
  
  n_callers = length(callers)
  
  data_length = merge[,paste0("GT_",callers),with=F]
  
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
  
  out = as.character(data_length$consensus)
  
  return(out)
  
  
}

#### genotypes from callers ######

geno_callers <- function(merge,callers){
  
  n_callers = length(callers)
  
  data_length = merge[,paste0("GT_",callers),with=F]
  
  data_length$consensus = NA
  
  data_length = as.matrix(data_length)
  
  for(i in 1:nrow(data_length)){
    
    aux = data_length[i,1:n_callers]
    
    my_names = aux[which(data_length[i,1:n_callers]!="9/9")]
    
    data_length[i,n_callers+1] = paste(my_names,collapse = ",")
    
  }
  
  data_length = as.data.frame(data_length)
  
  out = as.character(data_length$consensus)
  
  return(out)
  
  
}



#### number of callers detecting the variant ######

n_callers_detected <- function(merge,callers){
  
  n_callers = length(callers)
  
  data_length = merge[,paste0("GT_",callers),with=F]
  
  data_length$consensus = NA
  
  data_length = as.matrix(data_length)
  
  for(i in 1:nrow(data_length)){
    
    my_tab = sum(is.na(data_length[i,1:n_callers]))
    
    data_length[i,n_callers+1] = n_callers-my_tab
    
  }
  
  
  data_length = as.data.frame(data_length)
  
  out = as.numeric(as.character.numeric_version(data_length$consensus))
  
  return(out)
  
  
}

#### name of callers detecting the variant ######

#merge = data_predict
#callers = c("lumpy","delly","pindel","whamg",
#            "manta","cnvnator","svaba")


name_callers_detected <- function(merge,callers){
  
  n_callers = length(callers)
  
  data_length = merge[,paste0("GT_",callers),with=F]
  
  data_length$consensus = NA
  
  data_length = as.matrix(data_length)
  
  for(i in 1:nrow(data_length)){
    
    aux = data_length[i,1:n_callers]
    
    my_names = names(aux)[which(data_length[i,1:n_callers]!="9/9")]
    
    my_names = gsub("GT_","",my_names)
    
    data_length[i,n_callers+1] = paste(my_names,collapse = ",")
    
  }
  
  data_length = as.data.frame(data_length)
  
  out = as.character(data_length$consensus)
  
  return(out)
  
  
}


#### minimum window size of callers detecting the variant ######

#merge = data_predict
#callers = c("lumpy","delly","pindel","whamg",
#            "manta","cnvnator","svaba")


min_window_size <- function(merge,callers,call_windows){
  
  n_callers = length(callers)
  
  data_length = merge[,paste0("GT_",callers),with=F]
  
  data_length$consensus = NA
  
  data_length = as.matrix(data_length)
  
  for(i in 1:nrow(data_length)){
    
    aux = data_length[i,1:n_callers]
    
    my_names = names(aux)[which(data_length[i,1:n_callers]!="9/9")]
    
    my_names = gsub("GT_","",my_names)
    
    out = min(call_windows$window[call_windows$caller %in% my_names])
    
    data_length[i,n_callers+1] = out
    
  }
  
  data_length = as.data.frame(data_length)
  
  out = as.character(data_length$consensus)
  
  return(as.numeric(out))
  
  
}


### reprocicity ######

reprocicity <- function(merge,callers){
  
  n_callers = length(callers)
  
  data_length = merge[,paste0("length_",callers),with=F]
  
  data_length$consensus = NA
  
  data_length = as.matrix(data_length)
  
  for(i in 1:nrow(data_length)){
    
    if(sum(is.na(data_length[i,1:n_callers]))<(n_callers-1)){
      
      end = max(data_length[i,1:n_callers],na.rm = T)
      start = min(data_length[i,1:n_callers],na.rm = T)
      
      data_length[i,n_callers+1] = round(start/end*100,2)
      
      
    }
    
    if(sum(is.na(data_length[i,1:n_callers]))==(n_callers-1)){
      
      data_length[i,n_callers+1] = 100
      
      
    }
    
  }
  
  data_length = as.data.frame(data_length)
  
  out = as.numeric(as.character.numeric_version(data_length$consensus))
  
  return(out)
  
  
}


## strategy ####

strategy <- function(merge,callers,strategies){
  
  n_callers = length(callers)
  
  data_length = merge[,paste0("GT_",callers),with=F]
  
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
  
  out = as.numeric(as.character.numeric_version(data_length$consensus))
  
  return(out)
  
  
}


### logical sensitivity specificity > 1 caller, >2 callers etc  ####

logical_sens_spec <- function(x){
  
  tp = as.numeric(table(x$PASS)[2])
  fp = nrow(x)-tp
  
  sens_call = round(tp/N*100,2)
  spec_call = round(tp/(fp+tp)*100,2)
  
  f_score = round((2*sens_call*spec_call)/(sens_call+spec_call)/100,2)
  
  IC_lower_sens = round(as.numeric(sens_call) - 1.96*sqrt(as.numeric(sens_call)*(100-as.numeric(sens_call))/N),digits=2)
  IC_upper_sens = round(as.numeric(sens_call) + 1.96*sqrt(as.numeric(sens_call)*(100-as.numeric(sens_call))/N),digits=2)
  IC_sens = paste0("[",IC_lower_sens,",",IC_upper_sens,"]")
  
  IC_lower_spec = round(as.numeric(spec_call) - 1.96*sqrt(as.numeric(spec_call)*(100-as.numeric(spec_call))/(fp+tp)),digits=2)
  IC_upper_spec = round(as.numeric(spec_call) + 1.96*sqrt(as.numeric(spec_call)*(100-as.numeric(spec_call))/(fp+tp)),digits=2)
  IC_spec = paste0("[",IC_lower_spec,",",IC_upper_spec,"]")
  
  
  return(c(tp,fp,sens_call,IC_sens,spec_call,IC_spec,f_score))
  
  
}


### consensus merge samples #####

consensus_merge_samples <- function(x,ids_variants,var){
  
  var_consensus = NULL
  
  for(i in ids_variants){
    
    aux = x[x$ID %in% i,]
    
    aux = floor(median(as.matrix(aux[,paste0(var,ids),with=F]),na.rm=T))
    
    var_consensus = c(var_consensus,aux)
    
  }
  
  
  return(var_consensus)
  
  
}


## merge callers 2 breakpoints (Translocations) #####


merge_callers_bp_trans = function(call1,call2,callers_ref,caller_to_merge,
                                  golden=FALSE,type,merge=FALSE){
  
  # call2: son las variantes de referencia para la intersección
  # call1: al que vas a buscar las posiciones
  # callers_ref: los callers de referencia que hay en call2
  # caller_to_merge: el caller que vamos a hacer el merge
  # golden: si el merge es con la golden
  
  #call1 = pindel_chr
  #call2 = my_data
  #caller_to_merge = "pindel"
  #callers_ref = c("manta","lumpy","delly")
  #call2$ID = paste0("translocation",1:nrow(call2))
  #call1$ID = "none"
  #golden=FALSE 
  
  n_callers = length(callers_ref)
  
  for(i in 1:nrow(call1)){
    
    k = 1
    
    while(k<n_callers+1){
      
      aux = call1[i,c(paste0("chr_1_",caller_to_merge),paste0(c("lower_","upper_"),caller_to_merge,c("_bp1","_bp1"))),with=F]
      aux2 = call1[i,c(paste0("chr_2_",caller_to_merge),paste0(c("lower_","upper_"),caller_to_merge,c("_bp2","_bp2"))),with=F]
      
      call1_interval = (as.numeric(aux[1,2])):(as.numeric(aux[1,3]))
      call1_chr = as.character(aux[1,1])
      
      call1_interval_bp = (as.numeric(aux2[1,2])):(as.numeric(aux2[1,3]))
      call1_chr_bp = as.character(aux2[1,1])
      
      call2_find = call2[as.character(t(call2[,paste0("chr_1_",callers_ref[k]),with=F]))==call1_chr,
                         c("ID",paste0(c("lower_","upper_"),callers_ref[k],c("_bp1","_bp1"))),with=F]
      call2_find_bis = call2[as.character(t(call2[,paste0("chr_2_",callers_ref[k]),with=F]))==call1_chr,
                             c("ID",paste0(c("lower_","upper_"),callers_ref[k],c("_bp2","_bp2"))),with=F]
      
      colnames(call2_find_bis) = colnames(call2_find)
      
      call2_find = rbind(call2_find,call2_find_bis)
      
      call2_find_bp = call2[as.character(t(call2[,paste0("chr_1_",callers_ref[k]),with=F]))==call1_chr_bp,
                            c("ID",paste0(c("lower_","upper_"),callers_ref[k],c("_bp1","_bp1"))),with=F]
      call2_find_bp_bis = call2[as.character(t(call2[,paste0("chr_2_",callers_ref[k]),with=F]))==call1_chr_bp,
                                c("ID",paste0(c("lower_","upper_"),callers_ref[k],c("_bp2","_bp2"))),with=F]
      
      colnames(call2_find_bp_bis) = colnames(call2_find_bp)
      
      call2_find_bp = rbind(call2_find_bp,call2_find_bp_bis)
      
      
      if(sum(is.na(call2_find[,2]))>0){call2_find = call2_find[which(!is.na(call2_find[,2])),]} # remove NA of some variants
      if(sum(is.na(call2_find_bp[,2]))>0){call2_find_bp = call2_find_bp[which(!is.na(call2_find_bp[,2])),]} # remove NA of some variants
      
      if(nrow(call2_find)>0){
        
        for(j in 1:nrow(call2_find)){
          
          call_interval = as.numeric(call2_find[j,2]):as.numeric(call2_find[j,3])
          
          golden_start = c(call1_interval,call1_interval_bp)
          
          my_find = sum(golden_start %in% call_interval)
          
          
          if((my_find>0)){ # the first bp
            
            call1[i,]$ID = as.character(call2_find[j,]$ID)
            
            k = n_callers + 1
            
          }
        }
      }
      
      # Si es troba en un dels dos chromosomes ja es la mateixa variant (els 2 chromosomes si que coincideixen obligatoriament)
      
      if(nrow(call2_find_bp)>0){
        
        for(j in 1:nrow(call2_find_bp)){
          
          call_interval_bp = as.numeric(call2_find_bp[j,2]):as.numeric(call2_find_bp[j,3])
          
          golden_start = c(call1_interval,call1_interval_bp)
          
          my_find_bp = sum(golden_start %in% call_interval_bp)
          
          if((my_find_bp>0)){ # the second bp
            
            call1[i,]$ID = as.character(call2_find_bp[j,]$ID)
            
            k = n_callers + 1
            
          }
        }
      }
      
      
      
      k = k + 1
      
    }
    
    print(i)
    
    
  }
  
  # IDs unic per posicions properes (loop by chromosome)
  
  new_variants = call1[call1$ID=="none",]
  
  if(nrow(new_variants)==1){
    id = nrow(call2)+1
    new_variants$ID[1] = paste0(type,"_",id)
    call1[call1$ID=="none",]$ID = new_variants$ID}
  
  if(nrow(new_variants)>1){
    id = nrow(call2)+1
    
    new_variants$ID[1] = paste0(type,"_",id)
    
    id = id + 1
    
    for(i in 2:nrow(new_variants)){
      
      aux = new_variants[(i-1):i,c(as.vector(outer(c("start_1_","start_2_"),caller_to_merge,paste0))),with=F]
      
      if(abs(aux$start_1[1]-aux$start_1[2])<10 | abs(aux$start_2[1]-aux$start_2[2])<10){ #variantes proximas
        
        new_variants$ID[i] = paste0(new_variants$ID[i-1])}
      
      
      else{
        
        new_variants$ID[i] = paste0(type,"_",id)
        
        id = id + 1
        
      }
      
      
    }
    
    call1[call1$ID=="none",]$ID = new_variants$ID
    
  }
  
  
  
  bp_names = as.vector(outer(c("lower_","upper_"),callers_ref,paste0))
  
  bp_names1 = as.vector(outer(bp_names,c("_bp1"),paste0))
  bp_names2 = as.vector(outer(bp_names,c("_bp2"),paste0))
  
  
  if(golden==FALSE){
    
    my_data = full_join(call2[,c("ID",
                                 as.vector(outer(c("chr_1_","chr_2_","start_1_","start_2_","length_",
                                                   "GT_"),callers_ref,paste0)),
                                 bp_names1,bp_names2),with=F],
                        call1[,c("ID",
                                 paste0(c("chr_1_","chr_2_","start_1_","start_2_","length_",
                                          "GT_"),caller_to_merge),
                                 paste0(c("lower_","upper_"),caller_to_merge,c("_bp1","_bp1")),
                                 paste0(c("lower_","upper_"),caller_to_merge,c("_bp2","_bp2"))),with=F])
  }
  if(golden==TRUE & merge==FALSE){
    
    
    my_data = right_join(call2[,c("ID",
                                  as.vector(outer(c("chr_1_","chr_2_","start_1_","start_2_","length_",
                                                    "GT_"),callers_ref,paste0)),
                                  bp_names1,bp_names2),with=F],
                         call1[,c("ID",
                                  paste0(c("chr_1_","chr_2_","start_1_","start_2_","length_",
                                           "GT_"),caller_to_merge),
                                  paste0(c("lower_","upper_"),caller_to_merge,c("_bp1","_bp1")),
                                  paste0(c("lower_","upper_"),caller_to_merge,c("_bp2","_bp2"))),with=F])
  }
  if(golden==TRUE & merge==TRUE){
    
    
    my_data = left_join(call2[,c("ID",
                                 as.vector(outer(c("chr_1_","chr_2_","start_1_","start_2_","length_",
                                                   "GT_"),callers_ref,paste0)),
                                 bp_names1,bp_names2),with=F],
                        call1[,c("ID",
                                 paste0(c("chr_1_","chr_2_","start_1_","start_2_","length_",
                                          "GT_"),caller_to_merge),
                                 paste0(c("lower_","upper_"),caller_to_merge,c("_bp1","_bp1")),
                                 paste0(c("lower_","upper_"),caller_to_merge,c("_bp2","_bp2"))),with=F])
  }
  
  
  return(as.data.table(my_data))
  
}


## sensitivity-precision 2 breakpoints (Translocations) #####


sensitivity_precision_geno_error_bp_trans <- function(call,golden,type,call_name){
  
  # call = delly
  # golden = golden
  # type = "translocation"
  # call_name = "delly"
  
  golden$ID = paste0(type,1:nrow(golden))
  
  call$ID = "none"
  call = merge_callers_bp_trans(call1=call,call2=golden,callers_ref =  c("golden"),
                                caller_to_merge = call_name,golden = TRUE,type)
  
  all = call
  
  # sensitivity precision
  
  if(sum(duplicated(call$ID))>1){
    
    call_sens = call[-which(duplicated(call$ID)),]}
  
  if(sum(duplicated(call$ID))==0){
    
    call_sens = call}
  
  
  call_sens$PASS = "YES"
  call_sens$PASS[is.na(call_sens$length_golden)] = "NO"
  table(call_sens$PASS)
  
  sens = (table(call_sens$PASS)[2]/nrow(golden))*100
  
  spec = (table(call_sens$PASS)[2]/nrow(call_sens))*100
  
  true_positive = table(call_sens$PASS)[2]
  false_positive = nrow(call_sens) - table(call_sens$PASS)[2]
  
  
  # genotype error
  
  call = call[!is.na(call$GT_golden),c("ID",paste0("GT_",call_name),"GT_golden"),with=F]
  call = as.data.frame(unique(call))
  
  if(sum(duplicated(call$ID))>0){
    
    var_dup = call[which(duplicated(call$ID)),]$ID
    
    call[call$ID %in% var_dup,paste0("GT_",call_name)] = "./."
    call = unique(call)
    
  }
  
  geno1 = as.character(t(call[,paste0("GT_",call_name)]))
  geno2 = as.character(t(call$GT_golden))
  
  
  if(sum(geno1 %in% "9/9")>0){  
    which_99 = which(geno1 %in% "9/9")
    geno1 = geno1[-which_99]
    geno2 = geno2[-which_99]}
  
  
  g.error = 100-sum(geno1==geno2)/length(geno1)*100
  
  
  return(list(sensitivity = sens, precision = spec,
              true_positive = true_positive,
              false_positive = false_positive,
              f1_score = ((2*sens*spec)/(sens+spec))/100,
              g_error = g.error,
              all = all))
  
}




### merge samples mid Del, Deletions > 150 and Inversions #####

merge_samples = function(sample1,sample2,sample_new,samples_old){
  
  # sample2: variantes de las muestras con el merge hecho
  # sample1: la nueva muestra
  # sample_new: ID nueva muestra
  # samples_old: ID muestras con el merge hecho
  
  #sample_new = "CWGS026"
  #samples_old = "CWGS025"
  #sample1 = fread(paste0("merge_callers/",ids[1],"_Del_chr_20"),header = T)
  #sample2 = my_samples
  #sample1$ID = paste0("deletion_",1:nrow(sample1))
  #sample2$ID = "none"
  
  n_samples = length(samples_old)
  
  for(i in 1:nrow(sample1)){
    
    k=1 # the number of sample to compare 
    
    while(k<n_samples+1){
      
      aux = sample1[i,c(paste0(c("lower_","upper_","length_"),sample_new)),with=F]
      
      sample1_interval = (as.numeric(aux[1,1])):(as.numeric(aux[1,2]))
      
      sample2_find = sample2[,c("ID",paste0(c("lower_","upper_","length_"),samples_old[k])),with=F]
      
      if(sum(is.na(sample2_find[,2]))>0){sample2_find = sample2_find[which(!is.na(sample2_find[,2])),]} # remove NA of some variants
      
      if(nrow(sample2_find)>0){
        
        for(j in 1:nrow(sample2_find)){
          
          my_find = findInterval(sample1_interval,
                                 as.numeric(sample2_find[j,2:3]))
          
          reprocicity = aux[1,3]/as.numeric(sample2_find[j,4])
          
          if(reprocicity>1){reprocicity = 1/reprocicity}
          
          if(sum(my_find==1)>0 & reprocicity>=0.80){
            
            sample1[i,]$ID = as.character(sample2_find[j,]$ID)
            
            k = n_samples + 1
            
          }
        }
      }  
      
      k = k + 1
      
    }
    
    print(i)
    
    
  }
  
  new_deletions = nrow(sample1[sample1$ID=="none",])
  
  sample1[sample1$ID=="none",]$ID = paste0("deletion_",(nrow(sample2)+1):(nrow(sample2)+new_deletions))
  
  my_data = full_join(sample1,sample2)
  
  return(as.data.table(my_data))
  
}



### merge samples Insertions #####

merge_samples_insertions = function(sample1,sample2,sample_new,samples_old){
  
  # sample2: variantes de las muestras con el merge hecho
  # sample1: la nueva muestra
  # sample_new: ID nueva muestra
  # samples_old: ID muestras con el merge hecho
  
  #sample_new = "CWGS026"
  #samples_old = "CWGS025"
  #sample1 = fread(paste0("merge_callers/",ids[1],"_Del_chr_20"),header = T)
  #sample2 = my_samples
  #sample1$ID = paste0("deletion_",1:nrow(sample1))
  #sample2$ID = "none"
  
  n_samples = length(samples_old)
  
  for(i in 1:nrow(sample1)){
    
    k=1 # the number of sample to compare 
    
    while(k<n_samples+1){
      
      aux = sample1[i,c(paste0(c("lower_","upper_"),sample_new)),with=F]
      
      sample1_interval = (as.numeric(aux[1,1])):(as.numeric(aux[1,2]))
      
      sample2_find = sample2[,c("ID",paste0(c("lower_","upper_"),samples_old[k])),with=F]
      
      if(sum(is.na(sample2_find[,2]))>0){sample2_find = sample2_find[which(!is.na(sample2_find[,2])),]} # remove NA of some variants
      
      if(nrow(sample2_find)>0){
        
        for(j in 1:nrow(sample2_find)){
          
          my_find = findInterval(sample1_interval,
                                 as.numeric(sample2_find[j,2:3]))
          
          #reprocicity = aux[1,3]/as.numeric(sample2_find[j,4])
          
          #if(reprocicity>1){reprocicity = 1/reprocicity}
          
          if(sum(my_find==1)>0){
            
            sample1[i,]$ID = as.character(sample2_find[j,]$ID)
            
            k = n_samples + 1
            
          }
        }
      }  
      
      k = k + 1
      
    }
    
    print(i)
    
    
  }
  
  new_deletions = nrow(sample1[sample1$ID=="none",])
  
  sample1[sample1$ID=="none",]$ID = paste0("insertion_",(nrow(sample2)+1):(nrow(sample2)+new_deletions))
  
  my_data = full_join(sample1,sample2)
  
  return(as.data.table(my_data))
  
}




## merge callers 2 breakpoints (Translocations) new #####


merge_callers_bp_trans_bis = function(call1,call2,callers_ref,caller_to_merge,
                                      golden=FALSE,type,merge=FALSE){
  
  # call2: son las variantes de referencia para la intersección
  # call1: al que vas a buscar las posiciones
  # callers_ref: los callers de referencia que hay en call2
  # caller_to_merge: el caller que vamos a hacer el merge
  # golden: si el merge es con la golden
  
  # call1 = manta_ok
  # call2 = golden
  # caller_to_merge = "manta"
  # callers_ref = c("golden")
  # call2$ID = paste0("translocation_",1:nrow(call2))
  # call1$ID = "none"
  # golden=FALSE
  # 
  # i = 1
  
  n_callers = length(callers_ref)
  
  for(i in 1:nrow(call1)){
    
    k = 1
    
    while(k<n_callers+1){
      
      aux = call1[i,c(paste0("chr_1_",caller_to_merge),
                      paste0(c("lower_","upper_"),caller_to_merge,c("_start_1","_start_1")),
                      paste0(c("lower_","upper_"),caller_to_merge,c("_end_1","_end_1")),
                      paste0("length_",caller_to_merge)),with=F]
      
      aux2 = call1[i,c(paste0("chr_2_",caller_to_merge),
                       paste0(c("lower_","upper_"),caller_to_merge,c("_start_2","_start_2")),
                       paste0(c("lower_","upper_"),caller_to_merge,c("_end_2","_end_2"))),with=F]
      
      call1_interval = c((as.numeric(aux[1,2])):(as.numeric(aux[1,3])),(as.numeric(aux[1,4])):(as.numeric(aux[1,5])))
      call1_chr = as.character(aux[1,1])
      
      call1_interval_bp = c((as.numeric(aux2[1,2])):(as.numeric(aux2[1,3])),(as.numeric(aux2[1,4])):(as.numeric(aux2[1,5])))
      call1_chr_bp = as.character(aux2[1,1])
      
      call2_find = call2[as.character(t(call2[,paste0("chr_1_",callers_ref[k]),with=F]))==call1_chr & 
                           as.character(t(call2[,paste0("chr_2_",callers_ref[k]),with=F]))==call1_chr_bp,
                         c("ID",paste0(c("lower_","upper_"),callers_ref[k],c("_start_1","_start_1")),
                           paste0(c("length_","chr_1_","chr_2_"),callers_ref[k])),with=F]
      
      call2_find_bis = call2[as.character(t(call2[,paste0("chr_1_",callers_ref[k]),with=F]))==call1_chr & 
                               as.character(t(call2[,paste0("chr_2_",callers_ref[k]),with=F]))==call1_chr_bp,
                             c("ID",paste0(c("lower_","upper_"),callers_ref[k],c("_end_1","_end_1")),
                               paste0(c("length_","chr_1_","chr_2_"),callers_ref[k])),with=F]
      
      colnames(call2_find_bis) = colnames(call2_find)
      
      call2_find = rbind(call2_find,call2_find_bis)
      
      call2_find_bp = call2[as.character(t(call2[,paste0("chr_1_",callers_ref[k]),with=F]))==call1_chr & 
                              as.character(t(call2[,paste0("chr_2_",callers_ref[k]),with=F]))==call1_chr_bp,
                            c("ID",paste0(c("lower_","upper_"),callers_ref[k],c("_start_2","_start_2")),
                              paste0(c("length_","chr_1_","chr_2_"),callers_ref[k])),with=F]
      
      call2_find_bp_bis = call2[as.character(t(call2[,paste0("chr_1_",callers_ref[k]),with=F]))==call1_chr & 
                                  as.character(t(call2[,paste0("chr_2_",callers_ref[k]),with=F]))==call1_chr_bp,
                                c("ID",paste0(c("lower_","upper_"),callers_ref[k],c("_end_2","_end_2")),
                                  paste0(c("length_","chr_1_","chr_2_"),callers_ref[k])),with=F]
      
      colnames(call2_find_bp_bis) = colnames(call2_find_bp)
      
      call2_find_bp = rbind(call2_find_bp,call2_find_bp_bis)
      
      
      if(sum(is.na(call2_find[,2]))>0){call2_find = call2_find[which(!is.na(call2_find[,2])),]} # remove NA of some variants
      if(sum(is.na(call2_find_bp[,2]))>0){call2_find_bp = call2_find_bp[which(!is.na(call2_find_bp[,2])),]} # remove NA of some variants
      
      
      merge_bp1 = 0 
      id_bp1 = NULL
      
      if(nrow(call2_find)>0){
        
        for(j in 1:nrow(call2_find)){
          
          call_interval = as.numeric(call2_find[j,2]):as.numeric(call2_find[j,3])
          
          golden_start = c(call1_interval,call1_interval_bp)
          
          my_find = sum(golden_start %in% call_interval)
          
          reprocicity = as.numeric(aux[,4])/as.numeric(call2_find[j,4])
          
          reprocicity = ifelse(reprocicity>1,as.numeric(call2_find[j,4])/as.numeric(aux[,4]),reprocicity)
          
          if((my_find>0)){ # the first bp
            
            # call1[i,]$ID = as.character(call2_find[j,]$ID)
            # 
            # k = n_callers + 1
            # 
            # print(as.character(call2_find[j,]$ID))
            
            merge_bp1 = merge_bp1 + 1
            id_bp1 = c(id_bp1,as.character(call2_find_bp[j,]$ID))
            
          }
        }
      }
      
      # Si es troba en un dels dos chromosomes ja es la mateixa variant (els 2 chromosomes si que coincideixen obligatoriament)
      
      merge_bp2 = 0 
      
      id_bp2 = NULL
      
      if(nrow(call2_find_bp)>0){
        
        for(j in 1:nrow(call2_find_bp)){
          
          call_interval_bp = as.numeric(call2_find_bp[j,2]):as.numeric(call2_find_bp[j,3])
          
          golden_start = c(call1_interval,call1_interval_bp)
          
          my_find_bp = sum(golden_start %in% call_interval_bp)
          
          reprocicity = as.numeric(aux[,4])/as.numeric(call2_find_bp[j,4])
          
          reprocicity = ifelse(reprocicity>1,as.numeric(call2_find_bp[j,4])/as.numeric(aux[,4]),reprocicity)
          
          if((my_find_bp>0)){ # the second bp
            
            # call1[i,]$ID = as.character(call2_find_bp[j,]$ID)
            # 
            # k = n_callers + 1
            # print(as.character(call2_find_bp[j,]$ID))
            
            merge_bp2 = merge_bp2 + 1
            
            id_bp2 = c(id_bp2,as.character(call2_find_bp[j,]$ID))
            
          }
        }
      }
      
      
      if(merge_bp1>0 & merge_bp2>0){
        
        id = names(sort(table(c(id_bp1,id_bp2)),decreasing = T))[1]
        
        call1[i,]$ID = id
        
        k = n_callers + 1
      }
      
      
      k = k + 1
      
    }
    
    print(i)
    
    
  }
  
  # IDs unic per posicions properes (loop by chromosome)
  
  new_variants = call1[call1$ID=="none",]
  
  new_variants$ID = as.character(new_variants$ID)
  
  if(nrow(new_variants)==1){
    id = nrow(call2)+1
    new_variants$ID[1] = paste0(type,"_",id)
    call1[call1$ID=="none",]$ID = new_variants$ID}
  
  if(nrow(new_variants)>1){
    
    id = nrow(call2)+1
    
    new_variants$ID[1] = paste0(type,"_",id)
    
    id = id + 1
    
    for(i in 2:nrow(new_variants)){
      
      aux = new_variants[(i-1):i,c(as.vector(outer(c("start_1_","start_2_","end_1_","end_2_","chr_1_","chr_2_"),caller_to_merge,paste0))),with=F]
      
      if((aux[1,5]==aux[2,5]) & (aux[1,6]==aux[2,6])  & 
         (abs(aux[1,1]-aux[2,1])<10 | abs(aux[1,2]-aux[2,2])<10 | 
          abs(aux[1,3]-aux[2,3])<10 | abs(aux[1,4]-aux[2,4])<10)){ #variantes proximas
        
        new_variants$ID[i] = paste0(new_variants$ID[i-1])}
      
      
      else{
        
        new_variants$ID[i] = paste0(type,"_",id)
        
        id = id + 1
        
      }
      
      
    }
    
    call1[call1$ID=="none",]$ID = new_variants$ID
    
  }
  
  call1$ID = as.character(call1$ID)
  call2$ID = as.character(call2$ID)
  
  bp_names = as.vector(outer(c("lower_","upper_"),callers_ref,paste0))
  
  bp_names1 = as.vector(outer(bp_names,c("_start_1"),paste0))
  bp_names2 = as.vector(outer(bp_names,c("_start_2"),paste0))
  bp_names3 = as.vector(outer(bp_names,c("_end_1"),paste0))
  bp_names4 = as.vector(outer(bp_names,c("_end_2"),paste0))
  
  if(golden==FALSE){
    
    my_data = full_join(call2[,c("ID",
                                 as.vector(outer(c("chr_1_","chr_2_","start_1_","start_2_","end_1_","end_2_","length_",
                                                   "GT_"),callers_ref,paste0)),
                                 bp_names1,bp_names2,bp_names3,bp_names4),with=F],
                        call1[,c("ID",
                                 paste0(c("chr_1_","chr_2_","start_1_","start_2_","end_1_","end_2_","length_",
                                          "GT_"),caller_to_merge),
                                 paste0(c("lower_","upper_"),caller_to_merge,c("_start_1","_start_1")),
                                 paste0(c("lower_","upper_"),caller_to_merge,c("_start_2","_start_2")),
                                 paste0(c("lower_","upper_"),caller_to_merge,c("_end_1","_end_1")),
                                 paste0(c("lower_","upper_"),caller_to_merge,c("_end_2","_end_2"))),with=F])
  }
  if(golden==TRUE & merge==FALSE){
    
    
    my_data = right_join(call2[,c("ID",
                                  as.vector(outer(c("chr_1_","chr_2_","start_1_","start_2_","end_1_","end_2_","length_",
                                                    "GT_"),callers_ref,paste0)),
                                  bp_names1,bp_names2,bp_names3,bp_names4),with=F],
                         call1[,c("ID",
                                  paste0(c("chr_1_","chr_2_","start_1_","start_2_","end_1_","end_2_","length_",
                                           "GT_"),caller_to_merge),
                                  paste0(c("lower_","upper_"),caller_to_merge,c("_start_1","_start_1")),
                                  paste0(c("lower_","upper_"),caller_to_merge,c("_start_2","_start_2")),
                                  paste0(c("lower_","upper_"),caller_to_merge,c("_end_1","_end_1")),
                                  paste0(c("lower_","upper_"),caller_to_merge,c("_end_2","_end_2"))),with=F])
  }
  if(golden==TRUE & merge==TRUE){
    
    
    my_data = left_join(call2[,c("ID",
                                 as.vector(outer(c("chr_1_","chr_2_","start_1_","start_2_","end_1_","end_2_","length_",
                                                   "GT_"),callers_ref,paste0)),
                                 bp_names1,bp_names2,bp_names3,bp_names4),with=F],
                        call1[,c("ID",
                                 paste0(c("chr_1_","chr_2_","start_1_","start_2_","end_1_","end_2_","length_",
                                          "GT_"),caller_to_merge),
                                 paste0(c("lower_","upper_"),caller_to_merge,c("_start_1","_start_1")),
                                 paste0(c("lower_","upper_"),caller_to_merge,c("_start_2","_start_2")),
                                 paste0(c("lower_","upper_"),caller_to_merge,c("_end_1","_end_1")),
                                 paste0(c("lower_","upper_"),caller_to_merge,c("_end_2","_end_2"))),with=F])
  }
  
  
  return(as.data.table(my_data))
  
}


## sensitivity-precision 2 breakpoints (Translocations) #####


sensitivity_precision_geno_error_bp_trans_bis <- function(call,golden,type,call_name){
  
  #call = pindel_ok
  #golden = golden
  #type = "translocation_"
  #call_name = "pindel"
  
  golden$ID = paste0(type,1:nrow(golden))
  
  call$ID = "none"
  call = merge_callers_bp_trans_bis(call1=call,call2=golden,callers_ref =  c("golden"),
                                    caller_to_merge = call_name,golden = TRUE,type)
  
  all = call
  
  # sensitivity precision
  
  if(sum(duplicated(call$ID))>=1){
    
    call_sens = call[-which(duplicated(call$ID)),]}
  
  if(sum(duplicated(call$ID))==0){
    
    call_sens = call}
  
  
  call_sens$PASS = "YES"
  call_sens$PASS[is.na(call_sens$length_golden)] = "NO"
  table(call_sens$PASS)
  
  sens = (table(call_sens$PASS)[2]/nrow(golden))*100
  
  spec = (table(call_sens$PASS)[2]/nrow(call_sens))*100
  
  IC_lower_sens = round(as.numeric(sens) - 1.96*sqrt(as.numeric(sens)*(100-as.numeric(sens))/nrow(golden)),digits=2)
  IC_upper_sens = round(as.numeric(sens) + 1.96*sqrt(as.numeric(sens)*(100-as.numeric(sens))/nrow(golden)),digits=2)
  IC_sens = paste0("[",IC_lower_sens,",",IC_upper_sens,"]")
  
  IC_lower_spec = round(as.numeric(spec) - 1.96*sqrt(as.numeric(spec)*(100-as.numeric(spec))/nrow(call)),digits=2)
  IC_upper_spec = round(as.numeric(spec) + 1.96*sqrt(as.numeric(spec)*(100-as.numeric(spec))/nrow(call)),digits=2)
  IC_spec = paste0("[",IC_lower_spec,",",IC_upper_spec,"]")
  
  
  true_positive = table(call_sens$PASS)[2]
  false_positive = nrow(call_sens) - table(call_sens$PASS)[2]
  
  
  # genotype error
  
  call = call[!is.na(call$GT_golden),c("ID",paste0("GT_",call_name),"GT_golden"),with=F]
  call = as.data.frame(unique(call))
  
  if(sum(duplicated(call$ID))>0){
    
    var_dup = call[which(duplicated(call$ID)),]$ID
    
    call[call$ID %in% var_dup,paste0("GT_",call_name)] = "./."
    call = unique(call)
    
  }
  
  geno1 = as.character(t(call[,paste0("GT_",call_name)]))
  geno2 = as.character(t(call$GT_golden))
  
  
  if(sum(geno1 %in% "9/9")>0){  
    which_99 = which(geno1 %in% "9/9")
    geno1 = geno1[-which_99]
    geno2 = geno2[-which_99]}
  
  
  g.error = 100-sum(geno1==geno2)/length(geno1)*100
  
  homo_golden = nrow(call[call$GT_golden=="1/1",])
  het_golden = nrow(call[call$GT_golden=="0/1",])
  
  
  homo_golden_het_caller = nrow(call[call$GT_golden=="1/1" & geno1=="0/1",])
  homo_golden_hom_caller = nrow(call[call$GT_golden=="1/1" & geno1=="1/1",])
  het_golden_het_caller = nrow(call[call$GT_golden=="0/1" & geno1=="0/1",])
  het_golden_hom_caller = nrow(call[call$GT_golden=="0/1" & geno1=="1/1",])
  
  g_error_het = 100-nrow(call[call$GT_golden=="0/1" & geno1=="0/1",])/het_golden*100
  g_error_hom = 100-nrow(call[call$GT_golden=="1/1" & geno1=="1/1",])/homo_golden*100
  
  return(list(sensitivity = sens, precision = spec,
              true_positive = true_positive,
              false_positive = false_positive,
              f1_score = ((2*sens*spec)/(sens+spec))/100,
              g_error = g.error,
              n_het = het_golden,
              n_homo = homo_golden,
              g_error_het = g_error_het,
              g_error_hom = g_error_hom,
              IC_sens = IC_sens,
              IC_spec = IC_spec,
              all = all))
  
}



