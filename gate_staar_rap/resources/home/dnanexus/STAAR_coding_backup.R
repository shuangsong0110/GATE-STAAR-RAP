print(commandArgs(TRUE))
args<- commandArgs(TRUE)
arrayid<- as.numeric(args[1])
region <- as.character(args[2])
chr <- arrayid
##################################  User specify (required) ##############################################

token <- '' ### 1. get from DNAnexus
if(T){
  system("pip3 install setuptools-rust")
  system("pip3 install --upgrade --ignore-installed pip")
  system("pip3 install dxpy")
  system(paste0("dx login --token ", token))
}

system(paste0("dx download -f [path_to_nullmodel].rda")) ### 2. null model
system(paste0("dx download -f [path_to_nullmodel].varianceRatio.txt")) ### 3. varRatio
system(paste0("dx download -f [path_to_Annotation_name_catalog.RData]")) ### 4. annotation catalog
output_path <- '[path_to_destination]'
system(paste0("dx mkdir -p ",output_path)) ### 5. output path in DNAnexus

system(paste0("dx download -r ",output_path,"/* --overwrite")) ### 6. for resumable file transfer
gds.path <- paste0("ukb.500k.wgs.chr",chr,".pass.annotated.gds") ### 7. agds path
if(!file.exists(paste0(gds.path))){
  system(paste0("dx download -f [path_to_gds_files]",gds.path))
  print("GDS File Downloaded from UKB")
  print(list.files())
}
load('./[nullmodel].rda')
library(data.table)
ratio <- fread(paste0('[nullmodel].varianceRatio.txt'))$V1


##################################  User specify (if necessary) ##############################################
rare_maf_cutoff=0.01;
rv_num_cutoff=2;
varRatio=ratio;
QC_label= 'annotation/info/QC_label2'
#QC_label= 'annotation/filter'
variant_type=c("SNV");
geno_missing_imputation=c("mean");
Annotation_dir <- "annotation/info/FunctionalAnnotation"
Use_annotation_weights <- TRUE
## Annotation name
Annotation_name <- c("CADD","LINSIGHT","FATHMM.XF","aPC.EpigeneticActive","aPC.EpigeneticRepressed","aPC.EpigeneticTranscription",
                     "aPC.Conservation","aPC.LocalDiversity","aPC.Mappability","aPC.TF","aPC.Protein")

if(region=='coding'){
output_path1 <- ".//step3_coding_output//"
dir.create(output_path1)
## output file name
output_file_name <- paste0("Coding_")
}
if(region=='noncoding'){
  output_path1 <- ".//step3_noncoding_output//"
  dir.create(output_path1)
  ## output file name
  output_file_name <- paste0("Noncoding_")
}

########################################################################################

if(T){
  coding_survival <- function(chr,gene_name,genofile,obj_nullmodel,genes,
                              rare_maf_cutoff=0.01,rv_num_cutoff=2,varRatio=1,
                              QC_label="annotation/filter",variant_type=c("SNV","Indel","variant"),geno_missing_imputation=c("mean","minor"),
                              Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,
                              Use_annotation_weights=c(TRUE,FALSE),Annotation_name=NULL,
                              SPA_p_filter=FALSE,p_filter_cutoff=0.05,silent=FALSE){
    
    ## evaluate choices
    #t1 <- proc.time()
    variant_type <- match.arg(variant_type)
    geno_missing_imputation <- match.arg(geno_missing_imputation)
    
    phenotype.id <- as.character(obj_nullmodel$id_include)
    n_pheno <- obj_nullmodel$n.pheno
    
    ## SPA status
    use_SPA <- TRUE
    use_SPA0 <- F
    # if(!is.null(obj_nullmodel$use_SPA))
    # {
    #   use_SPA <- obj_nullmodel$use_SPA
    # }else
    # {
    #   use_SPA <- FALSE
    # }
    
    ## get SNV id, position, REF, ALT (whole genome)
    filter <- seqGetData(genofile, QC_label)
    if(variant_type=="variant")
    {
      SNVlist <- filter == "PASS"
    }
    
    if(variant_type=="SNV")
    {
      SNVlist <- (filter == "PASS") & isSNV(genofile)
    }
    
    if(variant_type=="Indel")
    {
      SNVlist <- (filter == "PASS") & (!isSNV(genofile))
    }
    
    position <- as.numeric(seqGetData(genofile, "position"))
    variant.id <- seqGetData(genofile, "variant.id")
    
    rm(filter)
    #gc()
    
    ### Gene
    kk <- which(genes[,1]==gene_name)
    
    sub_start_loc <- genes[kk,3]
    sub_end_loc <- genes[kk,4]
    
    is.in <- (SNVlist)&(position>=sub_start_loc)&(position<=sub_end_loc)
    variant.id.gene <- variant.id[is.in]
    rm(SNVlist)
    rm(is.in)
    rm(position)
    #gc()
    
    seqSetFilter(genofile,variant.id=variant.id.gene,sample.id=phenotype.id)
    ##### HWE
    # print('HWE')
    # print(length(variant.id.gene))
    # print(length(hwe_all))
    # AF <- seqGetData(genofile, 'annotation/info/AF')$data
    # print(summary(AF))
    # print(length(AF))
    # AC <- as.numeric(AF)*length(phenotype.id)*2
    # #mac0 <- minorAlleleCount(genofile)
    # variant.id.gene <- variant.id.gene[which(AC<20 | hwe_all>1e-9)]
    # rm(AC)
    # rm(AF)
    # print(length(variant.id.gene))
    # seqSetFilter(genofile,variant.id=variant.id.gene)
    # print('HWE done')
    
    ## Gencode_Exonic
    GENCODE.EXONIC.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.EXONIC.Category")]))
    ## Gencode
    GENCODE.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Category")]))
    ## Meta.SVM.Pred
    MetaSVM_pred <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="MetaSVM")]))
    
    ################################################
    #           Coding
    ################################################
    variant.id.gene <- seqGetData(genofile, "variant.id")
    lof.in.coding <- (GENCODE.EXONIC.Category=="stopgain")|(GENCODE.EXONIC.Category=="stoploss")|(GENCODE.Category=="splicing")|(GENCODE.Category=="exonic;splicing")|(GENCODE.Category=="ncRNA_splicing")|(GENCODE.Category=="ncRNA_exonic;splicing")|(GENCODE.EXONIC.Category=="nonsynonymous SNV")|(GENCODE.EXONIC.Category=="synonymous SNV")
    variant.id.gene <- variant.id.gene[lof.in.coding]
    #hwe.gene <- hwe_all[lof.in.coding]
    
    seqSetFilter(genofile,variant.id=variant.id.gene)
    
    ## Gencode_Exonic
    GENCODE.EXONIC.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.EXONIC.Category")]))
    ## Gencode
    GENCODE.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Category")]))
    ## Meta.SVM.Pred
    MetaSVM_pred <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="MetaSVM")]))
    
    ## Annotation
    Anno.Int.PHRED.sub <- NULL
    Anno.Int.PHRED.sub.name <- NULL
    
    if(variant_type=="SNV")
    {
      if(Use_annotation_weights)
      {
        for(k in 1:length(Annotation_name))
        {
          if(Annotation_name[k]%in%Annotation_name_catalog$name)
          {
            Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,Annotation_name[k])
            Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name==Annotation_name[k])]))
            
            if(Annotation_name[k]=="CADD")
            {
              Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
            }
            
            if(Annotation_name[k]=="aPC.LocalDiversity")
            {
              Annotation.PHRED.2 <- -10*log10(1-10^(-Annotation.PHRED/10))
              Annotation.PHRED <- cbind(Annotation.PHRED,Annotation.PHRED.2)
              Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,paste0(Annotation_name[k],"(-)"))
            }
            Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub,Annotation.PHRED)
          }
        }
        
        Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
        colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
      }
    }
    
    ## genotype id
    id.genotype <- seqGetData(genofile,"sample.id")
    
    ################################################
    #                  plof_ds
    ################################################
    #t2 <- proc.time()
    variant.id.gene <- seqGetData(genofile, "variant.id")
    lof.in.plof <- (GENCODE.EXONIC.Category=="stopgain")|(GENCODE.EXONIC.Category=="stoploss")|(GENCODE.Category=="splicing")|(GENCODE.Category=="exonic;splicing")|(GENCODE.Category=="ncRNA_splicing")|(GENCODE.Category=="ncRNA_exonic;splicing")|((GENCODE.EXONIC.Category=="nonsynonymous SNV")&(MetaSVM_pred=="D"))
    variant.id.gene.category <- variant.id.gene[lof.in.plof]
    
    seqSetFilter(genofile,variant.id=variant.id.gene.category)
    
    
    # id.genotype.match <- rep(0,length(id.genotype))
    
    id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
    phenotype.id.merge <- data.frame(phenotype.id)
    phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
    id.genotype.match <- phenotype.id.merge$index
    
    ## Genotype
    Geno <- seqGetData(genofile, "$dosage")
    Geno <- Geno[id.genotype.match,,drop=FALSE]
    ## impute missing
    if(!is.null(dim(Geno)))
    {
      if(dim(Geno)[2]>0)
      {
        if(geno_missing_imputation=="mean")
        {
          tmp <- matrix_flip_mean(Geno)
          Geno <- tmp$Geno
          MAF <- tmp$MAF
          rm(tmp)
          # gc()
        }
        if(geno_missing_imputation=="minor")
        {
          Geno <- matrix_flip_minor(Geno)$Geno
        }
      }
    }
    gc()
    ## Annotation
    Anno.Int.PHRED.sub.category <- Anno.Int.PHRED.sub[lof.in.plof,]
    #t3 <- proc.time()
    pvalues <- 1
    if(n_pheno == 1)
    {
      if(F)
      {
        try(pvalues <- STAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,
                             rv_num_cutoff=rv_num_cutoff),silent=silent)
      }else{
        #print('ooooo')
        #t1 <- proc.time()
        try(pvalues <- STAAR_Binary_SPA_survival(Geno,obj_nullmodel,Anno.Int.PHRED.sub.category,
                                                 rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
                                                 varRatio=varRatio,
                                                 SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff,
                                                 MAF=MAF),silent=silent)
        # t2 <- proc.time()
        # t2-t1
        #gc()
      }
    }else
    {
      try(pvalues <- MultiSTAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff),silent=silent)
    }
    
    results_plof_ds <- c()
    if(inherits(pvalues, "list"))
    {
      results_temp <- as.vector(genes[kk,])
      results_temp[3] <- "plof_ds"
      results_temp[2] <- chr
      results_temp[1] <- as.character(genes[kk,1])
      results_temp[4] <- pvalues$num_variant
      
      if(!use_SPA0)
      {
        results_temp <- c(results_temp,pvalues$cMAC,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
                          pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
                          pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)
      }else
      {
        results_temp <- c(results_temp,pvalues$cMAC,
                          pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_B)
      }
      
      results_plof_ds <- rbind(results_plof_ds,results_temp)
    }
    
    if(!is.null(results_plof_ds))
    {
      if(!use_SPA0)
      {
        colnames(results_plof_ds) <- colnames(results_plof_ds, do.NULL = FALSE, prefix = "col")
        colnames(results_plof_ds)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
        colnames(results_plof_ds)[(dim(results_plof_ds)[2]-1):dim(results_plof_ds)[2]] <- c("ACAT-O","STAAR-O")
      }else
      {
        colnames(results_plof_ds) <- colnames(results_plof_ds, do.NULL = FALSE, prefix = "col")
        colnames(results_plof_ds)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
        colnames(results_plof_ds)[dim(results_plof_ds)[2]] <- c("STAAR-B")
        
      }
    }
    #t4 <- proc.time()
    # print
    # print(t2-t1)
    # print(t3-t2)
    # print(t4-t3)
    #####################################################
    #                      plof
    #####################################################
    lof.in.plof <- (GENCODE.EXONIC.Category=="stopgain")|(GENCODE.EXONIC.Category=="stoploss")|(GENCODE.Category=="splicing")|(GENCODE.Category=="exonic;splicing")|(GENCODE.Category=="ncRNA_splicing")|(GENCODE.Category=="ncRNA_exonic;splicing")
    variant.id.gene.category <- variant.id.gene[lof.in.plof]
    
    seqSetFilter(genofile,variant.id=variant.id.gene.category)
    
    ## genotype id
    #id.genotype <- seqGetData(genofile,"sample.id")
    #id.genotype <- id.genotype0
    # id.genotype.match <- rep(0,length(id.genotype))
    
    # id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
    # phenotype.id.merge <- data.frame(phenotype.id)
    # phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
    # id.genotype.match <- phenotype.id.merge$index
    
    ## Genotype
    Geno <- seqGetData(genofile, "$dosage")
    Geno <- Geno[id.genotype.match,,drop=FALSE]
    #gc()
    ## impute missing
    if(!is.null(dim(Geno)))
    {
      if(dim(Geno)[2]>0)
      {
        if(geno_missing_imputation=="mean")
        {
          #Geno <- matrix_flip_mean(Geno)$Geno
          tmp <- matrix_flip_mean(Geno)
          Geno <- tmp$Geno
          MAF <- tmp$MAF
          rm(tmp)
          
        }
        if(geno_missing_imputation=="minor")
        {
          Geno <- matrix_flip_minor(Geno)$Geno
        }
      }
    }
    gc()
    ## Annotation
    Anno.Int.PHRED.sub.category <- Anno.Int.PHRED.sub[lof.in.plof,]
    
    pvalues <- 1
    if(n_pheno == 1)
    {
      if(!use_SPA)
      {
        try(pvalues <- STAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff),silent=silent)
      }else{
        try(pvalues <- STAAR_Binary_SPA_survival(Geno,obj_nullmodel,Anno.Int.PHRED.sub.category,
                                                 rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
                                                 SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff,
                                                 varRatio=varRatio, MAF=MAF),silent=silent)
        #gc()
      }
    }else
    {
      try(pvalues <- MultiSTAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff),silent=silent)
    }
    
    results_plof <- c()
    if(inherits(pvalues, "list"))
    {
      results_temp <- as.vector(genes[kk,])
      results_temp[3] <- "plof"
      results_temp[2] <- chr
      results_temp[1] <- as.character(genes[kk,1])
      results_temp[4] <- pvalues$num_variant
      
      if(!use_SPA0)
      {
        results_temp <- c(results_temp,pvalues$cMAC,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
                          pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
                          pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)
      }else
      {
        results_temp <- c(results_temp,pvalues$cMAC,
                          pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_B)
      }
      
      results_plof <- rbind(results_plof,results_temp)
    }
    
    if(!is.null(results_plof))
    {
      if(!use_SPA0)
      {
        colnames(results_plof) <- colnames(results_plof, do.NULL = FALSE, prefix = "col")
        colnames(results_plof)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
        colnames(results_plof)[(dim(results_plof)[2]-1):dim(results_plof)[2]] <- c("ACAT-O","STAAR-O")
      }else
      {
        colnames(results_plof) <- colnames(results_plof, do.NULL = FALSE, prefix = "col")
        colnames(results_plof)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
        colnames(results_plof)[dim(results_plof)[2]] <- c("STAAR-B")
      }
    }
    
    #############################################
    #             synonymous
    #############################################
    lof.in.synonymous <- (GENCODE.EXONIC.Category=="synonymous SNV")
    variant.id.gene.category <- variant.id.gene[lof.in.synonymous]
    
    seqSetFilter(genofile,variant.id=variant.id.gene.category)
    
    ## genotype id
    #id.genotype <- seqGetData(genofile,"sample.id")
    # id.genotype.match <- rep(0,length(id.genotype))
    
    # id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
    # phenotype.id.merge <- data.frame(phenotype.id)
    # phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
    # id.genotype.match <- phenotype.id.merge$index
    
    ## Genotype
    Geno <- seqGetData(genofile, "$dosage")
    Geno <- Geno[id.genotype.match,,drop=FALSE]
    #gc()
    ## impute missing
    if(!is.null(dim(Geno)))
    {
      if(dim(Geno)[2]>0)
      {
        if(geno_missing_imputation=="mean")
        {
          #Geno <- matrix_flip_mean(Geno)$Geno
          tmp <- matrix_flip_mean(Geno)
          Geno <- tmp$Geno
          MAF <- tmp$MAF
          rm(tmp)
          
        }
        if(geno_missing_imputation=="minor")
        {
          Geno <- matrix_flip_minor(Geno)$Geno
        }
      }
    }
    gc()
    ## Annotation
    Anno.Int.PHRED.sub.category <- Anno.Int.PHRED.sub[lof.in.synonymous,]
    
    pvalues <- 1
    if(n_pheno == 1)
    {
      if(!use_SPA)
      {
        try(pvalues <- STAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff),silent=silent)
      }else{
        try(pvalues <- STAAR_Binary_SPA_survival(Geno,obj_nullmodel,Anno.Int.PHRED.sub.category,
                                                 rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
                                                 SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff,
                                                 varRatio=varRatio,MAF=MAF),silent=silent)
        #gc()
      }
    }else
    {
      try(pvalues <- MultiSTAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff),silent=silent)
    }
    
    results_synonymous <- c()
    if(inherits(pvalues, "list"))
    {
      results_temp <- as.vector(genes[kk,])
      results_temp[3] <- "synonymous"
      results_temp[2] <- chr
      results_temp[1] <- as.character(genes[kk,1])
      results_temp[4] <- pvalues$num_variant
      
      if(!use_SPA0)
      {
        results_temp <- c(results_temp,pvalues$cMAC,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
                          pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
                          pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)
      }else
      {
        results_temp <- c(results_temp,pvalues$cMAC,
                          pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_B)
      }
      
      results_synonymous <- rbind(results_synonymous,results_temp)
    }
    
    if(!is.null(results_synonymous))
    {
      if(!use_SPA0)
      {
        colnames(results_synonymous) <- colnames(results_synonymous, do.NULL = FALSE, prefix = "col")
        colnames(results_synonymous)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
        colnames(results_synonymous)[(dim(results_synonymous)[2]-1):dim(results_synonymous)[2]] <- c("ACAT-O","STAAR-O")
      }else
      {
        colnames(results_synonymous) <- colnames(results_synonymous, do.NULL = FALSE, prefix = "col")
        colnames(results_synonymous)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
        colnames(results_synonymous)[dim(results_synonymous)[2]] <- c("STAAR-B")
      }
      
    }
    
    #################################################
    #        missense
    #################################################
    lof.in.missense <- (GENCODE.EXONIC.Category=="nonsynonymous SNV")
    variant.id.gene.category <- variant.id.gene[lof.in.missense]
    
    seqSetFilter(genofile,variant.id=variant.id.gene.category)
    
    ## genotype id
    #id.genotype <- seqGetData(genofile,"sample.id")
    # id.genotype.match <- rep(0,length(id.genotype))
    
    # id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
    # phenotype.id.merge <- data.frame(phenotype.id)
    # phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
    # id.genotype.match <- phenotype.id.merge$index
    
    ## Genotype
    Geno <- seqGetData(genofile, "$dosage")
    Geno <- Geno[id.genotype.match,,drop=FALSE]
    #()
    ## impute missing
    if(!is.null(dim(Geno)))
    {
      if(dim(Geno)[2]>0)
      {
        if(geno_missing_imputation=="mean")
        {
          tmp <- matrix_flip_mean(Geno)
          Geno <- tmp$Geno
          MAF <- tmp$MAF
          rm(tmp)
          #gc()
          
        }
        if(geno_missing_imputation=="minor")
        {
          Geno <- matrix_flip_minor(Geno)$Geno
        }
      }
    }
    gc()
    ## Annotation
    Anno.Int.PHRED.sub.category <- Anno.Int.PHRED.sub[lof.in.missense,]
    
    pvalues <- 1
    if(n_pheno == 1)
    {
      if(!use_SPA)
      {
        try(pvalues <- STAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff),silent=silent)
      }else{
        try(pvalues <- STAAR_Binary_SPA_survival(Geno,obj_nullmodel,Anno.Int.PHRED.sub.category,
                                                 rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
                                                 SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff,
                                                 varRatio=varRatio, MAF=MAF),silent=silent)
        #gc()
      }
    }else
    {
      try(pvalues <- MultiSTAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff),silent=silent)
    }
    
    results <- c()
    if(inherits(pvalues, "list"))
    {
      results_temp <- as.vector(genes[kk,])
      results_temp[3] <- "missense"
      results_temp[2] <- chr
      results_temp[1] <- as.character(genes[kk,1])
      results_temp[4] <- pvalues$num_variant
      
      if(!use_SPA0)
      {
        results_temp <- c(results_temp,pvalues$cMAC,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
                          pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
                          pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)
      }else
      {
        results_temp <- c(results_temp,pvalues$cMAC,
                          pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_B)
      }
      
      results <- rbind(results,results_temp)
    }
    
    #################################################
    #         disruptive missense
    #################################################
    lof.in.dmissense <- (GENCODE.EXONIC.Category=="nonsynonymous SNV")&(MetaSVM_pred=="D")
    variant.id.gene.category <- variant.id.gene[lof.in.dmissense]
    
    seqSetFilter(genofile,variant.id=variant.id.gene.category)
    
    ## genotype id
    #id.genotype <- seqGetData(genofile,"sample.id")
    # id.genotype.match <- rep(0,length(id.genotype))
    
    # id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
    # phenotype.id.merge <- data.frame(phenotype.id)
    # phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
    # id.genotype.match <- phenotype.id.merge$index
    
    ## Genotype
    Geno <- seqGetData(genofile, "$dosage")
    Geno <- Geno[id.genotype.match,,drop=FALSE]
    #gc()
    ## impute missing
    if(!is.null(dim(Geno)))
    {
      if(dim(Geno)[2]>0)
      {
        if(geno_missing_imputation=="mean")
        {
          tmp <- matrix_flip_mean(Geno)
          Geno <- tmp$Geno
          MAF <- tmp$MAF
          rm(tmp)
          
          
        }
        if(geno_missing_imputation=="minor")
        {
          Geno <- matrix_flip_minor(Geno)$Geno
        }
      }
    }
    gc()
    ## Annotation
    Anno.Int.PHRED.sub.category <- Anno.Int.PHRED.sub[lof.in.dmissense,]
    
    pvalues <- 1
    if(n_pheno == 1)
    {
      if(!use_SPA)
      {
        try(pvalues <- STAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff),silent=silent)
      }else{
        try(pvalues <- STAAR_Binary_SPA_survival(Geno,obj_nullmodel,Anno.Int.PHRED.sub.category,
                                                 rare_maf_cutoff=rare_maf_cutoff,
                                                 rv_num_cutoff=rv_num_cutoff,SPA_p_filter=SPA_p_filter,
                                                 p_filter_cutoff=p_filter_cutoff,
                                                 varRatio=varRatio, MAF=MAF),silent=silent)
        #gc()
      }
    }else
    {
      try(pvalues <- MultiSTAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff),silent=silent)
    }
    
    if(inherits(pvalues, "list"))
    {
      results_temp <- as.vector(genes[kk,])
      results_temp[3] <- "disruptive_missense"
      results_temp[2] <- chr
      results_temp[1] <- as.character(genes[kk,1])
      results_temp[4] <- pvalues$num_variant
      
      if(!use_SPA0)
      {
        results_temp <- c(results_temp,pvalues$cMAC,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
                          pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
                          pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)
      }else
      {
        results_temp <- c(results_temp,pvalues$cMAC,
                          pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_B)
      }
      
      results <- rbind(results,results_temp)
    }
    
    if(!is.null(results))
    {
      if(!use_SPA0)
      {
        colnames(results) <- colnames(results, do.NULL = FALSE, prefix = "col")
        colnames(results)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
        colnames(results)[(dim(results)[2]-1):dim(results)[2]] <- c("ACAT-O","STAAR-O")
      }else
      {
        colnames(results) <- colnames(results, do.NULL = FALSE, prefix = "col")
        colnames(results)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
        colnames(results)[dim(results)[2]] <- c("STAAR-B")
      }
      
      if(dim(results)[1]==1)
      {
        if(results[3]!="disruptive_missense")
        {
          if(!use_SPA0)
          {
            results <- cbind(results,matrix(1,1,6))
            colnames(results)[(dim(results)[2]-5):dim(results)[2]] <- c("SKAT(1,25)-Disruptive","SKAT(1,1)-Disruptive","Burden(1,25)-Disruptive","Burden(1,1)-Disruptive","ACAT-V(1,25)-Disruptive","ACAT-V(1,1)-Disruptive")
            results_missense <- results
            results_ds <- c()
          }else{
            results <- cbind(results,matrix(1,1,2))
            colnames(results)[(dim(results)[2]-1):dim(results)[2]] <- c("Burden(1,25)-Disruptive","Burden(1,1)-Disruptive")
            results_missense <- results
            results_ds <- c()
          }
        }else
        {
          results_missense <- c()
          results_ds <- results
          results <- c()
        }
      }
      
      if(!is.null(results))
      {
        if(dim(results)[1]==2)
        {
          if(!use_SPA0)
          {
            results_m <- c(results[1,],rep(0,6))
            names(results_m)[(length(results_m)-5):length(results_m)] <- c("SKAT(1,25)-Disruptive","SKAT(1,1)-Disruptive","Burden(1,25)-Disruptive","Burden(1,1)-Disruptive","ACAT-V(1,25)-Disruptive","ACAT-V(1,1)-Disruptive")
            results_m[(length(results_m)-5):length(results_m)] <- results[2,c("SKAT(1,25)","SKAT(1,1)","Burden(1,25)","Burden(1,1)","ACAT-V(1,25)","ACAT-V(1,1)")]
            apc_num <- (length(results_m)-19)/6
            p_seq <- c(1:apc_num,1:apc_num+(apc_num+1),1:apc_num+2*(apc_num+1),1:apc_num+3*(apc_num+1),1:apc_num+4*(apc_num+1),1:apc_num+5*(apc_num+1),(6*apc_num+9):(6*apc_num+14))
            results_m["STAAR-O"] <- CCT(as.numeric(results_m[6:length(results_m)][p_seq]))
            results_m["STAAR-S(1,25)"] <- CCT(as.numeric(results_m[6:length(results_m)][c(1:apc_num,6*apc_num+9)]))
            results_m["STAAR-S(1,1)"] <- CCT(as.numeric(results_m[6:length(results_m)][c(1:apc_num+(apc_num+1),6*apc_num+10)]))
            results_m["STAAR-B(1,25)"] <- CCT(as.numeric(results_m[6:length(results_m)][c(1:apc_num+2*(apc_num+1),6*apc_num+11)]))
            results_m["STAAR-B(1,1)"] <- CCT(as.numeric(results_m[6:length(results_m)][c(1:apc_num+3*(apc_num+1),6*apc_num+12)]))
            results_m["STAAR-A(1,25)"] <- CCT(as.numeric(results_m[6:length(results_m)][c(1:apc_num+4*(apc_num+1),6*apc_num+13)]))
            results_m["STAAR-A(1,1)"] <- CCT(as.numeric(results_m[6:length(results_m)][c(1:apc_num+5*(apc_num+1),6*apc_num+14)]))
            
            results_ds <- c()
            results_ds <- rbind(results_ds,results[2,])
            
            results <- c()
            results <- rbind(results,results_m)
          }else
          {
            results_m <- c(results[1,],rep(0,2))
            names(results_m)[(length(results_m)-1):length(results_m)] <- c("Burden(1,25)-Disruptive","Burden(1,1)-Disruptive")
            results_m[(length(results_m)-1):length(results_m)] <- results[2,c("Burden(1,25)","Burden(1,1)")]
            
            ## check whether the p-values is NA. If so, set NA equals 1.
            if(is.na(results_m[(length(results_m)-1)]))
            {
              results_m[(length(results_m)-1)] <- 1
            }
            
            if(is.na(results_m[length(results_m)]))
            {
              results_m[length(results_m)] <- 1
            }
            
            apc_num <- (length(results_m)-10)/2
            p_seq <- c(1:apc_num,1:apc_num+(apc_num+1),(length(results_m)-6):(length(results_m)-5))
            
            ## calculate STAAR-B
            pvalues_sub <- as.numeric(results_m[6:length(results_m)][p_seq])
            if(sum(is.na(pvalues_sub))>0)
            {
              if(sum(is.na(pvalues_sub))==length(pvalues_sub))
              {
                results_m["STAAR-B"] <- 1
              }else
              {
                ## not all NAs
                pvalues_sub <- pvalues_sub[!is.na(pvalues_sub)]
                if(sum(pvalues_sub[pvalues_sub<1])>0)
                {
                  ## not all ones
                  results_m["STAAR-B"] <- CCT(pvalues_sub[pvalues_sub<1])
                  
                }else
                {
                  results_m["STAAR-B"] <- 1
                  
                }
              }
            }else
            {
              if(sum(pvalues_sub[pvalues_sub<1])>0)
              {
                results_m["STAAR-B"] <- CCT(pvalues_sub[pvalues_sub<1])
              }else
              {
                results_m["STAAR-B"] <- 1
              }
            }
            
            ## calculate STAAR-B(1,25)
            pvalues_sub <- as.numeric(results_m[6:length(results_m)][c(1:apc_num,(length(results_m)-6))])
            if(sum(is.na(pvalues_sub))>0)
            {
              if(sum(is.na(pvalues_sub))==length(pvalues_sub))
              {
                results_m["STAAR-B(1,25)"] <- 1
              }else
              {
                ## not all NAs
                pvalues_sub <- pvalues_sub[!is.na(pvalues_sub)]
                if(sum(pvalues_sub[pvalues_sub<1])>0)
                {
                  ## not all ones
                  results_m["STAAR-B(1,25)"] <- CCT(pvalues_sub[pvalues_sub<1])
                  
                }else
                {
                  results_m["STAAR-B(1,25)"] <- 1
                  
                }
              }
            }else
            {
              if(sum(pvalues_sub[pvalues_sub<1])>0)
              {
                results_m["STAAR-B(1,25)"] <- CCT(pvalues_sub[pvalues_sub<1])
              }else
              {
                results_m["STAAR-B(1,25)"] <- 1
              }
            }
            
            ## calculate STAAR-B(1,1)
            pvalues_sub <- as.numeric(results_m[6:length(results_m)][c(1:apc_num+(apc_num+1),(length(results_m)-5))])
            if(sum(is.na(pvalues_sub))>0)
            {
              if(sum(is.na(pvalues_sub))==length(pvalues_sub))
              {
                results_m["STAAR-B(1,1)"] <- 1
              }else
              {
                ## not all NAs
                pvalues_sub <- pvalues_sub[!is.na(pvalues_sub)]
                if(sum(pvalues_sub[pvalues_sub<1])>0)
                {
                  ## not all ones
                  results_m["STAAR-B(1,1)"] <- CCT(pvalues_sub[pvalues_sub<1])
                  
                }else
                {
                  results_m["STAAR-B(1,1)"] <- 1
                  
                }
              }
            }else
            {
              if(sum(pvalues_sub[pvalues_sub<1])>0)
              {
                results_m["STAAR-B(1,1)"] <- CCT(pvalues_sub[pvalues_sub<1])
              }else
              {
                results_m["STAAR-B(1,1)"] <- 1
              }
            }
            
            results_ds <- c()
            results_ds <- rbind(results_ds,results[2,])
            
            results <- c()
            results <- rbind(results,results_m)
          }
        }
      }
    }else
    {
      results <- c()
      results_ds <- c()
    }
    
    results_coding <- list(plof=results_plof,plof_ds=results_plof_ds,missense=results,disruptive_missense=results_ds,synonymous=results_synonymous)
    
    seqResetFilter(genofile)
    
    return(results_coding)
  }
  
  
  # mac0 <- apply(Geno,2,sum)
  # genotype <- Geno[,which(mac0<10)]
  # annotation_phred=Anno.Int.PHRED.sub.category[which(mac0<10),]
  # max_iter=1000
  # SPA_p_filter=FALSE
  # p_filter_cutoff=0.05
  
  # STAAR_Binary_SPA_survival(Geno,obj_nullmodel,Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff
  # genotype <- Geno
  # max_iter=1000
  # SPA_p_filter=FALSE
  # p_filter_cutoff=0.05
  # annotation_phred=Anno.Int.PHRED.sub.category
  
  
  STAAR_Binary_SPA_survival <- function(genotype,obj_nullmodel,annotation_phred=NULL,
                                        rare_maf_cutoff=0.01,rv_num_cutoff=2,varRatio=1,
                                        tol=.Machine$double.eps^0.25,max_iter=1000,
                                        SPA_p_filter=FALSE,p_filter_cutoff=0.05, MAF, minpv=1e-45){
    #print('aaaaa')
    if(class(genotype)[1] != "matrix" && !(!is.null(attr(class(genotype), "package")) && attr(class(genotype), "package") == "Matrix")){
      stop("genotype is not a matrix!")
    }
    
    if(dim(genotype)[2] == 1){
      stop(paste0("Number of rare variant in the set is less than 2!"))
    }
    
    annotation_phred <- as.data.frame(annotation_phred)
    if(dim(annotation_phred)[1] != 0 & dim(genotype)[2] != dim(annotation_phred)[1]){
      stop(paste0("Dimensions don't match for genotype and annotation!"))
    }
    
    if(!is.null(attr(class(genotype), "package")) && attr(class(genotype), "package") == "Matrix"){
      genotype <- as.matrix(genotype)
    }
    #genotype <- matrix_flip(genotype)
    #MAF <- genotype$MAF
    RV_label <- as.vector((MAF<rare_maf_cutoff)&(MAF>0))
    # Geno_rare <- genotype$Geno[,RV_label]
    G <- genotype[,RV_label]
    rm(genotype)
    #gc()
    MAC <- round(MAF[RV_label]*2*nrow(G))
    
    # rm(genotype)
    # gc()
    annotation_phred <- annotation_phred[RV_label,,drop=FALSE]
    
    if(sum(RV_label) >= rv_num_cutoff){
      # G <- as(Geno_rare,"dgCMatrix")
      MAF <- MAF[RV_label]
      # rm(Geno_rare)
      # gc()
      
      annotation_rank <- 1 - 10^(-annotation_phred/10)
      
      ## beta(1,25)
      w_1 <- dbeta(MAF,1,25)
      ## beta(1,1)
      w_2 <- dbeta(MAF,1,1)
      if(dim(annotation_phred)[2] == 0){
        ## Burden, SKAT, ACAT-V
        w_B <- w_S <- as.matrix(cbind(w_1,w_2))
        w_A <- as.matrix(cbind(w_1^2/dbeta(MAF,0.5,0.5)^2,w_2^2/dbeta(MAF,0.5,0.5)^2))
        
      }else{
        ## Burden
        w_B_1 <- annotation_rank*w_1
        w_B_1 <- cbind(w_1,w_B_1)
        w_B_2 <- annotation_rank*w_2
        w_B_2 <- cbind(w_2,w_B_2)
        w_B <- cbind(w_B_1,w_B_2)
        w_B <- as.matrix(w_B)
        
        ## SKAT
        w_S_1 <- sqrt(annotation_rank)*w_1
        w_S_1 <- cbind(w_1,w_S_1)
        w_S_2 <- sqrt(annotation_rank)*w_2
        w_S_2 <- cbind(w_2,w_S_2)
        w_S <- cbind(w_S_1,w_S_2)
        w_S <- as.matrix(w_S)
        
        ## ACAT-V
        w_A_1 <- annotation_rank*w_1^2/dbeta(MAF,0.5,0.5)^2
        w_A_1 <- cbind(w_1^2/dbeta(MAF,0.5,0.5)^2,w_A_1)
        w_A_2 <- annotation_rank*w_2^2/dbeta(MAF,0.5,0.5)^2
        w_A_2 <- cbind(w_2^2/dbeta(MAF,0.5,0.5)^2,w_A_2)
        w_A <- cbind(w_A_1,w_A_2)
        w_A <- as.matrix(w_A)
      }
      
      # residuals.phenotype <- obj_nullmodel$residuals
      # muhat <- c(obj_nullmodel$mu)
      #print('aaaaa')
      #t1 <- proc.time()
      pvalues <- STAAR_B_Binary_SPA_survival(G=G,mac=MAC,
                                             #residuals=residuals.phenotype,
                                             varRatio=varRatio,
                                             #muhat=muhat,
                                             weights_B=w_B,weights_S=w_S,weights_A=w_A,
                                             tol=tol,max_iter=max_iter,obj_nullmodel=obj_nullmodel)
      #  t2 <- proc.time()
      #t2-t1
      #print(pvalues)
      pvalues[pvalues==0] <- minpv
      num_variant <- sum(RV_label) #dim(G)[2]
      cMAC <- sum(G)
      num_annotation <- dim(annotation_phred)[2]+1
      if(T){
        #CCT <- ACAT
        results_STAAR_O <- CCT(pvalues[1:(6*num_annotation)])
        results_ACAT_O <- CCT(pvalues[c(1,num_annotation+1,2*num_annotation+1,3*num_annotation+1,4*num_annotation+1,5*num_annotation+1)])
        pvalues_STAAR_S_1_25 <- CCT(pvalues[1:num_annotation])
        pvalues_STAAR_S_1_1 <- CCT(pvalues[(num_annotation+1):(2*num_annotation)])
        pvalues_STAAR_B_1_25 <- CCT(pvalues[(2*num_annotation+1):(3*num_annotation)])
        pvalues_STAAR_B_1_1 <- CCT(pvalues[(3*num_annotation+1):(4*num_annotation)])
        pvalues_STAAR_A_1_25 <- CCT(pvalues[(4*num_annotation+1):(5*num_annotation)])
        pvalues_STAAR_A_1_1 <- CCT(pvalues[(5*num_annotation+1):(6*num_annotation)])
        
        results_STAAR_S_1_25 <- c(pvalues[1:num_annotation],pvalues_STAAR_S_1_25)
        results_STAAR_S_1_25 <- data.frame(t(results_STAAR_S_1_25))
        
        results_STAAR_S_1_1 <- c(pvalues[(num_annotation+1):(2*num_annotation)],pvalues_STAAR_S_1_1)
        results_STAAR_S_1_1 <- data.frame(t(results_STAAR_S_1_1))
        
        results_STAAR_B_1_25 <- c(pvalues[(2*num_annotation+1):(3*num_annotation)],pvalues_STAAR_B_1_25)
        results_STAAR_B_1_25 <- data.frame(t(results_STAAR_B_1_25))
        
        results_STAAR_B_1_1 <- c(pvalues[(3*num_annotation+1):(4*num_annotation)],pvalues_STAAR_B_1_1)
        results_STAAR_B_1_1 <- data.frame(t(results_STAAR_B_1_1))
        
        results_STAAR_A_1_25 <- c(pvalues[(4*num_annotation+1):(5*num_annotation)],pvalues_STAAR_A_1_25)
        results_STAAR_A_1_25 <- data.frame(t(results_STAAR_A_1_25))
        
        results_STAAR_A_1_1 <- c(pvalues[(5*num_annotation+1):(6*num_annotation)],pvalues_STAAR_A_1_1)
        results_STAAR_A_1_1 <- data.frame(t(results_STAAR_A_1_1))
        
        if(dim(annotation_phred)[2] == 0){
          colnames(results_STAAR_S_1_25) <- c("SKAT(1,25)","STAAR-S(1,25)")
          colnames(results_STAAR_S_1_1) <- c("SKAT(1,1)","STAAR-S(1,1)")
          colnames(results_STAAR_B_1_25) <- c("Burden(1,25)","STAAR-B(1,25)")
          colnames(results_STAAR_B_1_1) <- c("Burden(1,1)","STAAR-B(1,1)")
          colnames(results_STAAR_A_1_25) <- c("ACAT-V(1,25)","STAAR-A(1,25)")
          colnames(results_STAAR_A_1_1) <- c("ACAT-V(1,1)","STAAR-A(1,1)")
        }else{
          colnames(results_STAAR_S_1_25) <- c("SKAT(1,25)",
                                              paste0("SKAT(1,25)-",colnames(annotation_phred)),
                                              "STAAR-S(1,25)")
          colnames(results_STAAR_S_1_1) <- c("SKAT(1,1)",
                                             paste0("SKAT(1,1)-",colnames(annotation_phred)),
                                             "STAAR-S(1,1)")
          colnames(results_STAAR_B_1_25) <- c("Burden(1,25)",
                                              paste0("Burden(1,25)-",colnames(annotation_phred)),
                                              "STAAR-B(1,25)")
          colnames(results_STAAR_B_1_1) <- c("Burden(1,1)",
                                             paste0("Burden(1,1)-",colnames(annotation_phred)),
                                             "STAAR-B(1,1)")
          colnames(results_STAAR_A_1_25) <- c("ACAT-V(1,25)",
                                              paste0("ACAT-V(1,25)-",colnames(annotation_phred)),
                                              "STAAR-A(1,25)")
          colnames(results_STAAR_A_1_1) <- c("ACAT-V(1,1)",
                                             paste0("ACAT-V(1,1)-",colnames(annotation_phred)),
                                             "STAAR-A(1,1)")
        }
        
        return(list(num_variant = num_variant,
                    cMAC = cMAC,
                    RV_label = RV_label,
                    results_STAAR_O = results_STAAR_O,
                    results_ACAT_O = results_ACAT_O,
                    results_STAAR_S_1_25 = results_STAAR_S_1_25,
                    results_STAAR_S_1_1 = results_STAAR_S_1_1,
                    results_STAAR_B_1_25 = results_STAAR_B_1_25,
                    results_STAAR_B_1_1 = results_STAAR_B_1_1,
                    results_STAAR_A_1_25 = results_STAAR_A_1_25,
                    results_STAAR_A_1_1 = results_STAAR_A_1_1))
      }
    }else{
      stop(paste0("Number of rare variant in the set is less than ",rv_num_cutoff,"!"))
    }
    
  }
  
  #library(ACAT)
  # mac=MAC
  # weights_B=w_B
  # weights_A=w_A
  # weights_S=w_S
  # mac_thres=10
  # 
  STAAR_B_Binary_SPA_survival <- function(G,mac,
                                          varRatio,
                                          weights_B,weights_S,weights_A,
                                          tol,max_iter, mac_thres=20, obj_nullmodel ){
    #G <- matrix_flip(G)$Geno
    
    nsnp00 <- ncol(G)
    muhat <- c(obj_nullmodel$mu)
    ind <- which(obj_nullmodel$y==1)
    ## mu0
    if(T){
      mu0 <- colSums(G * muhat)
      #gc()
      score_skat0_unweight <- colSums(G[ind, ])-mu0
      ii <- which(mu0==0)
      if(length(ii)>0){
        #ii <- tmp
        G <- G[,-ii]
        mu0 <- mu0[-ii]
        weights_B <- weights_B[-ii, ]
        weights_A <- weights_A[-ii, ]
        weights_S <- weights_S[-ii, ]
        score_skat0_unweight <- score_skat0_unweight[-ii]
        mac <- mac[-ii]
      }
    }
    #gc()
    #R_pois <- t(G) %*% (G * (muhat+muhat^2))
    library(Matrix)
    R_pois_direct <- function(G, muhat) {
      # G: a 'dgCMatrix' (or other sparse format), dimension n x p
      # muhat: numeric vector of length n
      # Returns: p x p matrix = t(G) %*% [G * (muhat + muhat^2)]
      #
      # We compute t(G) %*% (W %*% G), where W is diagonal with w_i = muhat_i + muhat_i^2.
      
      # w <- muhat + muhat^2
      W <- Diagonal(x = muhat)         # sparse diagonal matrix
      R <- crossprod(G, W %*% G)   # crossprod(G, X) = t(G) %*% X
      return(R)
    }
    R_pois<- R_pois_direct(G, muhat)
    var0 <- diag(R_pois)
    # library(propagate)
    R00 <- cov2cor(as.matrix(R_pois))
    ### LD
    if(T){
      R00_tmp <- abs(R00)
      R00_tmp[lower.tri(R00_tmp, diag=T)] <- 0
      #tmp <- abs(R00-diag(ncol(R00)))
      rr <- which(R00_tmp>sqrt(.95), arr.ind=T)
      ind_ld <- unique(rr[,2])
      if(length(ind_ld)>0){
        G <- G[,-ind_ld]
        weights_B <- weights_B[-ind_ld, ]
        weights_A <- weights_A[-ind_ld, ]
        weights_S <- weights_S[-ind_ld, ]
        mu0=mu0[-ind_ld]
        R00 <- R00[-ind_ld, -ind_ld]
        var0 <- var0[-ind_ld]
        mac=mac[-ind_ld]
        score_skat0_unweight <- score_skat0_unweight[-ind_ld]
      }
    }
    cs <- sqrt(colSums(abs(R00)))
    gc()
    # ww_mu0 <- (length(mu0)/mu0)/sum(1/mu0)
    # t2 <- proc.time()
    # t2-t1
    if(nsnp00-length(ind_ld)==1){
      wn <- length(weights_B)
      #res <- 1
      res <- rep(scoreTest_SAIGE_survivalTrait_cond_sparseSigma_fast(
        G0=G, obj.noK=obj_nullmodel, varRatio=varRatio)$p.value,5 * wn)
      return(res)
      
    }else{
      residuals <- obj_nullmodel$residuals
      vn <- nrow(G)
      un <- ncol(G)
      wn <- ncol(weights_B)
      
      G_cumu <- G %*% weights_B
      # Bisection initialization
      # xmin <- -100.0
      # xmax <- 100.0
      res <- rep(1,5 * wn)
      # print(res)
      # print(wn)
      # ACAT
      #mac00 <- apply(G[which(c(obj_nullmodel$y)==0),],2,sum)
      # mac11 <- apply(G[ind,],2,sum)
      # mac00 <- mac-mac11
      # mac_min <- pmin(mac00, mac11)
      # print(mac00)
      # print(mac11)
      # print(mac_min)
      # print(length(mac))
      # print(dim(G))
      id_veryrare <- which(mac <= mac_thres)
      id_common <- which(mac > mac_thres)
      # print(id_veryrare)
      # print(id_common)
      n0 <- length(id_veryrare)
      n1 <- length(id_common)
      # print(n0)
      # print(n1)
      pseq <- rep(1, un)
      wseq <- rep(0, un)
      G_sub <- G[,id_common]
      #t1 <- proc.time()
      if(n1>0){
        if(n1==1){
          pseq[1] <- scoreTest_SAIGE_survivalTrait_cond_sparseSigma_fast(
            G0=G_sub, obj.noK=obj_nullmodel, varRatio=varRatio)$p.value
        }else{
          #start_time <- proc.time()
          for (k in 1:n1) {
            #print('ccccc')
            # print(dim(G[,id_common[k]]))
            # print(dim(G_cumu[,1]))
            # print(dim(G))
            #print(k)
            #start_time <- proc.time()
            pseq[k] <- scoreTest_SAIGE_survivalTrait_cond_sparseSigma_fast(
              G0=G_sub[,k], obj.noK=obj_nullmodel, varRatio=varRatio)$p.value
            #end_time <- proc.time()
            #print(end_time-start_time)
            # if(k%%20==0){
            #   gc()
            # }
            #print('dddd')
          }
          #end_time <- proc.time()
          # print(end_time-start_time)
          
        }
      }  
      # t2 <- proc.time()
      # t2-t1
      rm(G_sub)
      if(n0>0){
        #print('ccccc')
        # print(dim(G))
        if(n0==1){
          G_cumu_ultra <- matrix(G[,id_veryrare], nrow=nrow(G), ncol=ncol(weights_B))
        }else{
          G_cumu_ultra <- G[,id_veryrare] %*% weights_B[id_veryrare, ]
        }
        #print('dddd')
      }
      gc()
      #print(res)
      ### calculate LD ratio for SKAT
      #RR <- cor(G)
      if(T){
        #score_skat0 <- weights*apply(g*(y-mu), 2, sum)^2
        #mu0 <- apply(G*c(obj_nullmodel$mu), 2, sum)
        ### gamma  skew & mean
        para <- compute_k_theta_v(mu0)
        # kk <- 4*mu0*(2*mu0+1)^3/(8*mu0^2+22*mu0+1)^2
        # theta1 <- mu0/kk
        # theta2 <-  sqrt((3*mu0^2+mu0)/(kk+kk^2))
        # theta <- (theta1+theta2)/2
        kk <- para[1,]
        theta <- para[2,]
        #R_pois <- cor2cov(R00^2, var=var0)
        
        # R2 <- R00*(2+R00)/3
        # R_pois <- cor2cov(R2, var=var0)
        #denominator <- outer(sqrt(G*kk*theta), sqrt(G*kk*theta))
      }
      
      # print('iter')
      # t1 <- proc.time()
      for (i in 1:wn) {
        
        #print(i)
        # sum0 <- 0.0
        # for (k in 1:n) {
        #   sum0 <- sum0 + x[k] * weights_B[k, i]
        # }
        
        # Calculate p-values
        #respart1 <- Saddle_Binary_SPA_survival(abs(sum0), muhat, G_cumu[, i], tol, max_iter, lower1 = FALSE, varRatio=varRatio)
        if(T){
          
          ### Burden
          # print(G_cumu[,i])
          # t1 <- proc.time()
          respart1 <- scoreTest_SAIGE_survivalTrait_cond_sparseSigma_fast(
            G0=G_cumu[,i]/sum(weights_B[,i]), obj.noK=obj_nullmodel, varRatio=varRatio)$p.value
          res[wn+i] <- respart1 
          
          ### ACAT
          if(n1>0){
            wseq[1:n1] <- weights_A[id_common,i]
          }
          if(n0==0){
            res[2*wn+i] <- CCT(pseq, weights = weights_A[,i]) 
          }else{
            #t1 <- proc.time()
            # print(summary(G_cumu_ultra[,i]))
            pseq[n1 + 1] <-  scoreTest_SAIGE_survivalTrait_cond_sparseSigma_fast(
              G0=G_cumu_ultra[,i]/sum(weights_B[id_veryrare, i]), obj.noK=obj_nullmodel, varRatio=varRatio)$p.value
            # t2 <- proc.time()
            # t2-t1
            wseq[n1 + 1] <- sum(weights_A[id_veryrare, i])/n0 
            res[2 * wn + i] <- CCT(pseq[1:(n1 + 1)], wseq[1:(n1 + 1)])
            
          }
          tmp0 <- pv_skat_ld5(#mu=c(obj_nullmodel$mu), 
            weights=weights_S[,i]^2/sum(weights_S[,i]^2), 
            score_skat0_unweight=score_skat0_unweight,cs=cs,
            #y=c(obj_nullmodel$y), g=G, 
            #R_pois=abs(R00),
            kk=kk, theta=theta, mu0=mu0, varRatio=varRatio, useLD=F, R_pois=abs(R00))
          aa <- try( tmp <- pv_skat_ld5(#mu=c(obj_nullmodel$mu), 
            weights=weights_S[,i]^2/sum(weights_S[,i]^2), 
            score_skat0_unweight=score_skat0_unweight,cs=cs,
            #y=c(obj_nullmodel$y), g=G, 
            #R_pois=abs(R00),
            kk=kk, theta=theta, mu0=mu0, varRatio=varRatio, useLD=T, R_pois=abs(R00)))
          if(inherits(aa, "try-error")){
            
            tmp <- pv_skat_ld4(#mu=c(obj_nullmodel$mu), 
              weights=weights_S[,i]^2/sum(weights_S[,i]^2), 
              score_skat0_unweight=score_skat0_unweight,cs=cs,
              #y=c(obj_nullmodel$y), g=G, 
              #R_pois=abs(R00),
              kk=kk, theta=theta, mu0=mu0, varRatio=varRatio)
          }
          res[i] <- max(tmp, tmp0)
          
          
          if(is.nan(res[i])){
            res[i] <- 0.999
          }
          #t2 <- proc.time()
          #t3=proc.time()
          #res[i] <- respart1 
          #respart1
          # }else{
          #  res[i] <- res[ wn + i]
          # }
          
          
          #   res[i] <-  res[wn+i]
          
          #  }
          # print('aaa')
          # print(t2-t1)
          # print(t3-t2)
          
        }
        
        # If saddle point approximation fails, set the p-value to 1
        if (res[i] > .9999) {
          res[i] <- .9999
        }
        if (res[wn+i] > .9999) {
          res[wn+i] <- .9999
        }
        if (res[2*wn+i] > .9999) {
          res[2*wn+i] <- .9999
        }
        if (res[3*wn+i] > .9999) {
          res[3*wn+i] <- .9999
        }
        if (res[4*wn+i] > .9999) {
          res[4*wn+i] <- .9999
        }
        # if(i==10){
        #   gc()
        # }
        if(i%%14==0){
          gc()
        }
        #gc()
        # t2 <- proc.time()
        # print(i)
        # print(t2-t1)
      }
      # t2 <- proc.time()
      #  t2-t1
      #gc()
      print('iter done')
      gc()
      return(res)
      
      
    }
    
    
  }
  
  
  
  Gene_Centric_Coding_survival <- function(chr,gene_name,category=c("all_categories","plof","plof_ds","missense","disruptive_missense","synonymous","ptv","ptv_ds","all_categories_incl_ptv"),
                                           genofile,obj_nullmodel,rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                           QC_label="annotation/filter",variant_type=c("SNV","Indel","variant"),geno_missing_imputation=c("mean","minor"),
                                           Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,
                                           Use_annotation_weights=c(TRUE,FALSE),Annotation_name=NULL,
                                           SPA_p_filter=TRUE,p_filter_cutoff=0.05,silent=FALSE,varRatio=1){
    
    ## evaluate choices
    category <- match.arg(category)
    variant_type <- match.arg(variant_type)
    geno_missing_imputation <- match.arg(geno_missing_imputation)
    
    genes <- genes_info[genes_info[,2]==chr,]
    
    #if(category=="all_categories")
    # {
    results <- coding_survival(chr,gene_name,genofile,obj_nullmodel,genes,
                               rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
                               QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                               Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                               Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                               SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff,silent=silent,varRatio=varRatio)
    # }
    
    return(results)
  }
  
  
  
  
  
  
  
  if(T){
    
    
    K_wgamma <- function(alpha, beta, weights, tt, R_pois=NA){
      if(!is.na(R_pois[1])){
        return(-sum(alpha*log(1-weights/beta*tt))+weights%*%(R_pois-diag(diag(R_pois)))%*%weights/2*tt^2)
      }else{
        return(-sum(alpha*log(1-weights/beta*tt)))
      }
    }
    K_wgamma_prime <- function(alpha, beta, weights, tt, R_pois=NA){
      if(!is.na(R_pois[1])){
        return(sum(alpha*weights/(beta-weights*tt))+weights%*%(R_pois-diag(diag(R_pois)))%*%weights*tt)
      }else{
        return(sum(alpha*weights/(beta-weights*tt)))
      }
    }
    K_wgamma_prime2 <- function(alpha, beta, weights, tt, R_pois=NA){
      if(!is.na(R_pois[1])){
        return(sum(alpha*weights^2/(beta-weights*tt)^2)+weights%*%(R_pois-diag(diag(R_pois)))%*%weights)
      }else{
        return(sum(alpha*weights^2/(beta-weights*tt)^2))
      }
    }
    K_wgamma_prime_adj <- function(alpha, beta, weights, tt, x, R_pois=NA){
      if(!is.na(R_pois[1])){
        return(sum(alpha*weights/(beta-weights*tt))+weights%*%(R_pois-diag(diag(R_pois)))%*%weights*tt-x)
      }else{  
        return(sum(alpha*weights/(beta-weights*tt))-x)
      }
    }
    
    pv_wgamma <- function(alpha, beta, weights, cut, R_pois=NA, pv.unadj=NULL){
      #mean0 <- sum(alpha*weights/beta)
      # sd0 <- sqrt(sum(alpha*weights^2/beta^2))
      #if(is.null(pv.unadj)){
      #pv.unadj0 <- 1-pnorm(cut, mean=mean0, sd=sd0)
      #}
      #if(is.null(pv.unadj)){
      #  pv.unadj <- pv.unadj0
      # }
      if(pv.unadj<0.001){
        t_hat <- uniroot(K_wgamma_prime_adj, 
                         lower = -100, upper = min(beta/weights) - 1e-6,
                         alpha=alpha, beta=beta, weights=weights, x=cut, R_pois=R_pois)$root
        #print(t_hat)
        w <- sign(t_hat)*(2*(t_hat*cut-K_wgamma(alpha=alpha, beta=beta, weights=weights,tt=t_hat)))^0.5
        v <- t_hat*(K_wgamma_prime2(alpha=alpha, beta=beta, weights=weights,tt=t_hat,R_pois=R_pois))^0.5
        pdf_approx <- 1-pnorm(w+1/w*log(v/w))
        SPA=T
      }else{
        pdf_approx <- pv.unadj
        SPA=F
      }
      if(pdf_approx==0){
        pdf_approx <- pv.unadj
        SPA=F
      }
      return(list(pv=pdf_approx, pv_unadj=pv.unadj, SPA=SPA))
    }
    
    
    # 
    # mu=c(obj_nullmodel$mu)
    # weights=weights_S[,i]/sum(weights_S[,i])
    # y=c(obj_nullmodel$y)
    # g=G
    pv_skat_ld4 <- function(weights,score_skat0_unweight,  cs, kk, theta, mu0, varRatio=1){
      # mu0 <- apply(g*mu, 2, sum)
      # mu02 <- apply(g*(mu^2), 2, sum)
      # denominator <- outer(sqrt(mu0+mu02), sqrt(mu0+mu02))
      # R_pois <- t(g) %*% (g * (mu+mu^2))/denominator
      #### score_skat0 <- weights*apply(g*(y-mu), 2, sum)^2
      score_skat0 <- weights*score_skat0_unweight^2
      #cs <- sqrt(colSums(R_pois))
      score_skat0 <- score_skat0/cs/sqrt(varRatio)
      
      #score_skat0 <- weights*apply(g*(y-mu), 2, sum)^2/sqrt(theta)
      #1-pgamma(cut, shape=kk, scale=theta)
      # res <- pv_wgamma(alpha=kk, beta=1, weights=w1, cut=sum(score_skat0/(ldratio)))
      # res$score <- sum(score_skat0)
      # res$ldratio <- ldratio
      #w1 <- weights*theta*sqrt(varRatio)
      #w1 <-  weights*theta
      #library(mgcv)
      pv.unadj <- psum.chisq(sum(score_skat0)/varRatio, df=rep(1,length(theta)),lb=weights*mu0)
      
      #w1 <- weights*sqrt(theta)
      # pv_wgamma(alpha=abs(shape_params_w), beta=1/abs(scale_params_w_new), 
      #           w=sign(shape_params_w*scale_params_w_new ), cut=sum(score_skat0))$pv
      return(pv_wgamma(alpha=kk, beta=1, weights=weights*theta, cut=sum(score_skat0), R_pois=NA, pv.unadj=pv.unadj)$pv)
    }
    gamma_weighted_decomp <- function(alpha, beta, w, rho, tol = 1e-3) {
      # ---- basic checks --------------------------------------------------------
      K <- length(alpha)
      if (length(beta) != K || length(w) != K)
        stop("alpha, beta, and w must have the same length (K)")
      if (!all(dim(rho) == c(K, K)))
        stop("rho must be a K×K matrix")
      
      # force diagonal to 1 for safety
      diag(rho) <- 1
      
      # ---- compute alpha_0 -----------------------------------------------------
      idx_lower <- lower.tri(rho)
      pair_sqrt <- sqrt(outer(alpha, alpha))
      alpha0 <- min(rho[idx_lower] * pair_sqrt[idx_lower])
      
      # ---- compute alpha_ij matrix --------------------------------------------
      alpha_ij <- matrix(0, K, K)
      for (i in 1:(K - 1)) {
        for (j in (i + 1):K) {
          alpha_ij[i, j] <- rho[i, j] * sqrt(alpha[i] * alpha[j]) - alpha0
          alpha_ij[j, i] <- alpha_ij[i, j]  # make symmetric for convenience
        }
      }
      
      # ---- compute alpha_z vector ---------------------------------------------
      alpha_z <- alpha - (alpha0 + rowSums(alpha_ij))
      
      # ---- feasibility check ---------------------------------------------------
      if (any(c(alpha0, alpha_ij[idx_lower], alpha_z) < -tol))
        stop("Derived shape parameters contain negative values; correlation matrix is infeasible under the Gamma additive model.")
      
      # Treat tiny negatives as zero
      alpha0          <- if (alpha0 < tol) 0 else alpha0
      alpha_ij[abs(alpha_ij) < tol] <- 0
      alpha_z[abs(alpha_z) < tol]   <- 0
      
      # ---- assemble independent components ------------------------------------
      alpha_list <- numeric(0)
      beta_list  <- numeric(0)
      name_list  <- character(0)
      
      # global component
      if (alpha0 > 0) {
        alpha_list <- c(alpha_list, alpha0)
        beta_list  <- c(beta_list, sum(w * beta))
        name_list  <- c(name_list, "U_global")
      }
      
      # pairwise components
      for (i in 1:(K - 1)) {
        for (j in (i + 1):K) {
          if (alpha_ij[i, j] > 0) {
            alpha_list <- c(alpha_list, alpha_ij[i, j])
            beta_list  <- c(beta_list, w[i] * beta[i] + w[j] * beta[j])
            name_list  <- c(name_list, sprintf("Y_%d_%d", i, j))
          }
        }
      }
      
      # private components
      for (i in 1:K) {
        if (alpha_z[i] > 0) {
          alpha_list <- c(alpha_list, alpha_z[i])
          beta_list  <- c(beta_list, w[i] * beta[i])
          name_list  <- c(name_list, sprintf("Z_%d", i))
        }
      }
      
      names(alpha_list) <- name_list
      names(beta_list)  <- name_list
      
      return(list(alpha_new = alpha_list, beta_new = beta_list))
    }
    
    pv_skat_ld5 <- function(weights,score_skat0_unweight,  cs, kk, theta, mu0, varRatio=1, useLD=T, R_pois){
      # mu0 <- apply(g*mu, 2, sum)
      # mu02 <- apply(g*(mu^2), 2, sum)
      # denominator <- outer(sqrt(mu0+mu02), sqrt(mu0+mu02))
      # R_pois <- t(g) %*% (g * (mu+mu^2))/denominator
      #### score_skat0 <- weights*apply(g*(y-mu), 2, sum)^2
      score_skat0 <- weights*score_skat0_unweight^2
      #cs <- sqrt(colSums(R_pois))
      ii <- which(cs>1.01)
      if(length(ii)>1 & useLD){
        tmp <-   gamma_weighted_decomp(kk[ii], theta[ii], weights[ii], R_pois[ii,ii])
        weights_new <- c(weights[-ii], rep(1,length(tmp$alpha_new)))
        kk_new <- c(kk[-ii], tmp$alpha_new)
        theta_new <- c(theta[-ii], tmp$beta_new)
      }else{
        theta_new <- theta
        kk_new <- kk
        weights_new <- weights
      }
      
      score_skat0 <- score_skat0/sqrt(varRatio)
      
      #score_skat0 <- weights*apply(g*(y-mu), 2, sum)^2/sqrt(theta)
      #1-pgamma(cut, shape=kk, scale=theta)
      # res <- pv_wgamma(alpha=kk, beta=1, weights=w1, cut=sum(score_skat0/(ldratio)))
      # res$score <- sum(score_skat0)
      # res$ldratio <- ldratio
      #w1 <- weights*theta*sqrt(varRatio)
      #w1 <-  weights*theta
      #library(mgcv)
      pv.unadj <- psum.chisq(sum(score_skat0)/varRatio, df=rep(1,length(theta)),lb=weights*mu0)
      
      #w1 <- weights*sqrt(theta)
      # pv_wgamma(alpha=abs(shape_params_w), beta=1/abs(scale_params_w_new), 
      #           w=sign(shape_params_w*scale_params_w_new ), cut=sum(score_skat0))$pv
      return(pv_wgamma(alpha=kk_new, beta=1, weights=weights_new*theta_new, cut=sum(score_skat0), R_pois=NA, pv.unadj=pv.unadj)$pv)
    }
    
    
  }
  
  if(T){
    compute_k_theta <- function(mu, d1=5, d2=7){
      cut1 <- floor(sqrt(uniroot(ff1_pv,interval=c(0,max(5000,500*mu)), mu=mu, pv=10^(-d1))$root))^2
      cut2 <- floor(sqrt(uniroot(ff1_pv,interval=c(0,max(5000,500*mu)), mu=mu, pv=10^(-d2))$root))^2
      p1 <- ff1(mu, cut1)
      p2 <- ff1(mu, cut2)
      kk <- 4*mu*(2*mu+1)^3/(8*mu^2+22*mu+1)^2
      theta <- mu/kk
      return(optim(c(kk, theta), cut1=cut1, cut2=cut2,p1=p1, p2=p2, fn=ff_gamma_optim)$par)
      
      
    }
    compute_k_theta_v <- Vectorize(compute_k_theta, 'mu')
    ff1 <- function(mu,cut){
      return( 1-ppois(sqrt(cut)+mu, lambda = mu)+ppois(-sqrt(cut)-mu, lambda=mu)
      )
    }
    
    
    ff1_pv <- function(mu,cut, pv){
      return( 1-ppois(sqrt(cut)+mu, lambda = mu)+ppois(-sqrt(cut)-mu, lambda=mu)-pv
      )
    }
    
    ff_gamma_optim <- function(xx, cut1, cut2, p1, p2){
      kk <- xx[1]
      theta <- xx[2]
      return((-log10(1-pgamma(cut1, shape=kk, scale=theta))+log10(p1))^2+
               (-log10(1-pgamma(cut2, shape=kk, scale=theta))+log10(p2))^2)
    }
    
    
  }
  
  
  
}
if(T){
  ### line 36: score change
  
  scoreTest_SAIGE_survivalTrait_cond_sparseSigma_fast=function(G0,  obj.noK, 
                                                               varRatio=1, Cutoff=2, 
                                                               sparseSigma=NULL, isCondition=FALSE, 
                                                               OUT_cond=NULL, G1tilde_P_G2tilde = NULL, 
                                                               G2tilde_P_G2tilde_inv=NULL, aa=NULL, type='origin', 
                                                               score_skat=NULL, pv_acat=NULL,wvar=0, wvar2=0){
    
    #print(summary(G0))
    #G0[is.na(G0)] <- mean(G0, na.rm=T)
    # if(is.na(G0[1])){
    #   print(length(G0))
    #   print(G0[1:5])
    # }
    #t1 <- proc.time()
    N = length(G0)
    AC <- sum(G0)
    AF <- AC/N/2
    MAF <- min(AF, 1-AF)
    mu.a <- obj.noK$mu
    y <- obj.noK$y
    mu2.a=mu.a ##poisson
    #print(table(G0))
    #print('aaaaa')
    if(AF > 0.5){
      G0 = 2-G0
      AC2 = 2*N - AC
    }else{
      AC2 = AC
    }
    Run1=TRUE
    idx_no0<-which(G0>0)
    G0 = matrix(G0, ncol = 1)
    
    #tp0 = proc.time()
    #mug = mean(G0)
    #tp1 = proc.time()
    #print("tp1-tp0")
    #print(tp1-tp0)
    ##########old way to get g_mc
    #G0_mc = matrix(G0-mean(G0), ncol = 1)
    #XVG0_mc = eigenMapMatMult(obj.noK$XV, G0_mc)
    #g_mc = G0_mc - eigenMapMatMult(obj.noK$XXVX_inv, XVG0_mc)
    #########
    #if(!isCondition){
    #  if(IsSparse==TRUE){
    if(type!='skat'){
      if(MAF < 0.05){
        out.score<- Score_Test_Sparse_Survival(obj.noK, G0, mu.a, mu2.a, varRatio);
        #tp2a = proc.time()
        #print("tp2a-tp1")
        #print(tp2a-tp1)
      }else{
        out.score<- Score_Test_Survival(obj.noK, G0, mu.a, mu2.a, varRatio);
        #tp2b = proc.time()
        #print("tp2b-tp1")
        #print(tp2b-tp1)
      }
    }
    #t2 <- proc.time()
    #cat("out.score$Tstat: ", out.score$Tstat, "\n")
    #cat("out.score1$Tstat: ", out.score1$Tstat, "\n")
    #cat("out.score$var2: ", out.score$var2, "\n")
    #cat("out.score1$var2: ", out.score1$var2, "\n")
    ##if(out.score["pval.noadj"] > 0.05){
    if(type=='skat'){
      # XVG0 = (obj.noK$XV_fg%*%G0)
      # g = G0 - (obj.noK$XXVX_inv_fg%*% XVG0)
      g <- G0
      
    }else{
      if(!isCondition){
        if(abs(as.numeric(unlist(out.score["Tstat"])[1])/sqrt(as.numeric(unlist(out.score["var1"])[1]))) < Cutoff){
          if(AF > 0.5){
            out.score$BETA = (-1)*out.score$BETA
            out.score$Tstat = (-1)*out.score$Tstat
          }
          outVec = list(BETA = out.score$BETA, SE = out.score$SE, Tstat = out.score$Tstat, p.value = out.score$pval.noadj,
                        p.value.NA = out.score$pval.noadj, Is.converge = 0, var1 = out.score$var1, var2 = out.score$var2)
          Run1=FALSE
          #print('PASS SPA')
        }else{
          if(MAF < 0.05){
            # XVG0 = eigenMapMatMult(obj.noK$XV_fg, G0)
            # g = G0 - eigenMapMatMult(obj.noK$XXVX_inv_fg, XVG0)
            XVG0 = (obj.noK$XV_fg%*%G0)
            g = G0 - (obj.noK$XXVX_inv_fg%*% XVG0)
            #tp3a = proc.time()
            #print("tp3a-tp2a")
            #print(tp3a-tp2a)
            #print(g[1:20])
            #print(g_mc[1:20])
          }else{
            g = out.score$g_tilde	
            #tp3b = proc.time()
            #print("tp3b-tp2b")
            #print(tp3b-tp2b)
          }
        }
      }else{ #if(!isCondition){
        if(MAF < 0.05){
          # XVG0 = eigenMapMatMult(obj.noK$XV_fg, G0)
          # g = G0 - eigenMapMatMult(obj.noK$XXVX_inv_fg, XVG0)
          XVG0 = (obj.noK$XV_fg%*% G0)
          g = G0 - (obj.noK$XXVX_inv_fg%*% XVG0)
        }else{
          g = out.score$g_tilde
          #tp3b = proc.time()
          #print("tp3b-tp2b")
          #print(tp3b-tp2b)
        }
        
      }
    }
    # t3 <- proc.time()
    # }
    #}else{ #if(!isCondition){
    
    #}
    
    #cat("Run1: ", Run1, "\n")
    if(Run1){
      if(type=='skat'){
        #NAset = which(G0==0)
        #tp4 = proc.time()
        
        # print('qqqqq')
        out1 = scoreTest_SPAGMMAT_survivalTrait_cond_sparseSigma_fast(g, Score = score_skat, pval.noadj = pv_acat, 
                                                                      var1_a = 1, var2_a = 1, AC2, AC,NAset, 
                                                                      y, mu.a, varRatio, Cutoff, sparseSigma=sparseSigma, isCondition=isCondition, 
                                                                      OUT_cond=OUT_cond, G1tilde_P_G2tilde = G1tilde_P_G2tilde, 
                                                                      G2tilde_P_G2tilde_inv=G2tilde_P_G2tilde_inv,
                                                                      aa=aa, type=type,Score_skat=score_skat,wvar=wvar,wvar2=wvar2)
        # print(out1)
        outVec = list(p.value = out1["p.value"]$p.value, 
                      p.value.NA = out1["p.value.NA"]$p.value.NA,
                      Is.converge=out1["Is.converge"]$Is.converge)
        
        return(outVec)
        
      }else{
        #G0 = matrix(G0, ncol = 1)
        #XVG0 = eigenMapMatMult(obj.noK$XV, G0)
        #G = G0  -  eigenMapMatMult(obj.noK$XXVX_inv, XVG0) # G is X adjusted
        #g = G
        NAset = which(G0==0)
        #tp4 = proc.time()
        out1 = scoreTest_SPAGMMAT_survivalTrait_cond_sparseSigma_fast(g, Score = out.score$Tstat, pval.noadj = out.score$pval.noadj, 
                                                                      var1_a = out.score$var1, var2_a = out.score$var2, AC2, AC,NAset, 
                                                                      y, mu.a, varRatio, Cutoff, sparseSigma=sparseSigma, isCondition=isCondition, 
                                                                      OUT_cond=OUT_cond, G1tilde_P_G2tilde = G1tilde_P_G2tilde, 
                                                                      G2tilde_P_G2tilde_inv=G2tilde_P_G2tilde_inv, aa=aa, type=type)
        # print(out1)
        
        if(AF > 0.5){
          out1$BETA = (-1)*out1$BETA
          out1$Tstat = (-1)*out1$Tstat
          if(isCondition){
            out1$BETA_c = (-1) * out1$BETA_c
            out1$Tstat_c = (-1) * out1$Tstat_c
          }
        }
        #tp5 = proc.time()
        #print("tp5-tp4")
        #print(tp5-tp4)	
        if(isCondition){
          outVec = list(BETA = out1["BETA"], SE = out1["SE"], Tstat = out1["Tstat"],p.value = out1["p.value"], 
                        p.value.NA = out1["p.value.NA"], Is.converge=out1["Is.converge"], var1 = out1["var1"], 
                        var2 = out1["var2"], Tstat_c = out1["Tstat_c"], p.value.c = out1["p.value.c"], 
                        var1_c = out1["var1_c"], BETA_c = out1["BETA_c"], SE_c = out1["SE_c"])
          
        }else{
          outVec = list(BETA = out1["BETA"]$BETA, SE = out1["SE"]$SE, 
                        Tstat = out1["Tstat"]$Tstat,p.value = out1["p.value"]$p.value, 
                        p.value.NA = out1["p.value.NA"]$p.value.NA,
                        Is.converge=out1["Is.converge"]$Is.converge, 
                        var1 = out1["var1"]$var1, var2 = out1["var2"]$var2)
          #outVec = list(BETA = BETA, SE = SE, Tstat = Tstat,p.value = p.value, var1 = var1, var2 = var2)
        }
      }
      #cat("p.value: ", as.numeric(outVec$p.value), "\n")
      #cat("p.value.NA: ", as.numeric(outVec$p.value.NA), "\n")
      # t4 <- proc.time()
      # print(t2-t1)
      # print(t3-t2)
      # print(t4-t3)
      #return(outVec)
    }#else{
    # outVec = list(p.value = 1,
    #                p.value.NA = 1)
    return(outVec)
    # }
  }
  
  
  # Score = score_skat
  # pval.noadj = 0.5
  # var1_a = 1
  # var2_a = 1
  # AC=AC2
  # AC_true=AC
  # mu <-  mu.a
  # Score_skat=score_skat
  scoreTest_SPAGMMAT_survivalTrait_cond_sparseSigma_fast=function(g, Score, pval.noadj, var1_a, var2_a, AC, AC_true, NAset, y, 
                                                                  mu, varRatio, Cutoff, sparseSigma=NULL, isCondition=FALSE, OUT_cond=NULL, 
                                                                  G1tilde_P_G2tilde = NULL, G2tilde_P_G2tilde_inv=NULL,
                                                                  type='origin', aa=NULL,Score_skat=NULL, wvar=0,wvar2=0){
    
    #g = G/sqrt(AC)
    #q = innerProduct(g, y)
    #print(type)
    #t1 <- proc.time()
    m1 = innerProduct(g, mu)
    #Tstat = q-m1
    #Tstat = Score
    #var2 = innerProduct(mu, g*g)
    #var2c_old = innerProduct(mu, g_mc*g_mc)
    #var2c = var2 - 2*mug*innerProduct(g,Wq) + mug^2*qW1
    
    #cat("var2c_old ", var2c_old, "\n")
    #cat("var2c ", var2c, "\n")
    #var1 = var2c * varRatio
    var1 = var1_a
    var2 = var2_a
    #t2 <- proc.time()
    if(!is.null(sparseSigma)){
      #pcginvSigma<-pcg(sparseSigma, g)
      pcginvSigma<-solve(sparseSigma, g, sparse=T)
      var2b = as.matrix(t(g) %*% pcginvSigma)
      var1 = var2b * varRatio
    }
    
    if(isCondition){
      T2stat = OUT_cond[,2]
      G1tilde_P_G2tilde = matrix(G1tilde_P_G2tilde,nrow=1)
      Tstat_c = Score - G1tilde_P_G2tilde %*% G2tilde_P_G2tilde_inv %*% T2stat
      var1_c = var1 - G1tilde_P_G2tilde %*% G2tilde_P_G2tilde_inv %*% t(G1tilde_P_G2tilde)
    }
    
    AF = AC_true/(2*length(y))
    #if(AF > 0.5){
    #  Tstat = (-1)*Tstat
    #  if(isCondition){
    #    Tstat_c = (-1)*Tstat_c
    #  }
    #}
    qq <- g
    #t3 <- proc.time()
    # t1 <- proc.time()
    if(type=='origin'){
      aa <- -mu*g
      
      #print(head(qq[-NAset]))
      #qtilde = Tstat/sqrt(var1) * sqrt(var2) + m1
      #Score2 = Score/sqrt(var1) * sqrt(var2)
      #Score2 = Score/sqrt(varRatio)
      #   Score2 = Score/sqrt(varRatio)
      # print('NAset')
      # print(length(NAset))
      Score2 = Score/sqrt(varRatio)
      #Score2 = Score
      if(length(NAset)/length(g) < 0.5){
        # print("Saddle_Prob_Poisson")
        out1 = Saddle_Prob_Poisson(Score=Score2, pval.noadj=pval.noadj, mu = mu, qq=qq, aa=aa, Cutoff = Cutoff, alpha=5*10^-8, m1=m1, var1=var2)
      }else{
        #  print("Saddle_Prob_Poisson_fast")
        out1 = Saddle_Prob_Poisson_fast(Score=Score2, pval.noadj=pval.noadj, qq=qq, aa=aa, mu = mu, qqNB=qq[-NAset], aaNB=aa[-NAset], muNA = mu[NAset], muNB = mu[-NAset], Cutoff = Cutoff, alpha = 5*10^-8, m1=m1, var1=var2)
        #print(out1)
        
      }
      #t3 <- proc.time()
      out1$var1 = var1
      out1$var2 = var2
      
      logOR = Score/var1
      SE = abs(logOR/qnorm(out1$p.value/2))
      out1$BETA=logOR
      out1$SE=SE
      out1$Tstat = Score
    }else if(type=='skat'){
      # <- sum(mu)^2-sum(mu^2)
      out1 = Saddle_Prob_Poisson_fast_skat(Score=Score_skat, pval.noadj=pval.noadj, qq=qq, aa=aa, mu = mu, wvar=wvar,wvar2=wvar2)
      
    }
    #print(t3-t2)
    # t4 <- proc.time()
    # t3 <- proc.time()
    if(isCondition){
      if(var1_c <= (.Machine$double.xmin)^2){
        out1 = c(out1, var1_c = var1_c,BETA_c = NA, SE_c = NA, Tstat_c = Tstat_c, p.value.c = 1, p.value.NA.c = 1)
      }else{
        
        #qtilde_c = Tstat_c/sqrt(var1_c) * sqrt(var2) + m1
        pval.noadj_c<-pchisq((Tstat_c)^2/(var1_c), lower.tail = FALSE, df=1)
        if(length(NAset)/length(g) < 0.5){
          ######To improve
          out1_c = Saddle_Prob_Poisson(Score=Tstat_c, pval.noadj=pval.noadj_c, mu = mu, qq=qq, aa=aa, Cutoff = Cutoff, alpha=5*10^-8, m1=m1, var1=var1_c)
        }else{
          ######To improve
          out1_c = Saddle_Prob_Poisson_fast(Score=Tstat_c, pval.noadj=pval.noadj_c, qq=qq, aa=aa, mu = mu,qqNB=qq[-NAset], aaNB=aa[-NAset], muNA = mu[NAset], muNB = mu[-NAset], Cutoff = Cutoff, alpha = 5*10^-8, m1=m1, var1=var1_c)
        }
        logOR_c = Tstat_c/var1_c
        SE_c = abs(logOR_c/qnorm(out1_c$p.value/2))
        out1 = c(out1, var1_c = var1_c,BETA_c = logOR_c, SE_c = SE_c, Tstat_c = Tstat_c, p.value.c = out1_c$p.value, p.value.NA.c = out1_c$p.value.NA)
      }
      
    }
    return(out1)
  }
  
  
  Saddle_Prob_Poisson_fast=function (Score, pval.noadj, mu, qq,aa, qqNB, aaNB, muNA,muNB, Cutoff = 2, alpha = 5*10^-8, m1, var1){
    #m1 <- sum(mu * g)
    #var1 <- sum(mu * g^2)
    p1 = NULL
    p2 = NULL
    
    #NAmu= m1-sum(gNB*muNB)
    if(F){
      NAsigma=var1-sum(muNB*gNB^2)
    }else{
      NAsigma=var1-sum(muNB*qqNB^2)
    }
    #cat("Score is ", Score, "\n")
    #cat("NAsigma is ", NAsigma, "\n")
    #print(mu[1:20])
    #print(g[1:20])
    
    #NAsigma = sum(muNA*gNA^2)
    #Score <- q - m1
    #qinv = -sign(q - m1) * abs(q - m1) + m1
    #pval.noadj <- pchisq((q - m1)^2/var1, lower.tail = FALSE,
    #    df = 1)
    #Is.converge = TRUE
    
    #if (abs(q - m1)/sqrt(var1) < Cutoff) {
    #    pval = pval.noadj
    #}else {
    #print("Saddle_Prob_Poisson_fast >= Cutoff")
    
    #Korg_Poi_result = Korg_Poi(t=0.1, mu, g)
    #Korg_Poi_fast_result = Korg_Poi_fast(t=0.1, mu, g, gNA,gNB,muNA,muNB,NAmu,NAsigma)	
    #cat("Korg_Poi_result: ", Korg_Poi_result, "\n")
    #cat("Korg_Poi_fast_result: ", Korg_Poi_fast_result, "\n")
    
    
    out.uni1 <- getroot_K1_Poi_fast(0, mu = mu, qq=qq, aa=aa,  q = Score, qqNB=qqNB, aaNB=aaNB, muNA=muNA,muNB=muNB,NAsigma=NAsigma)
    # print(out.uni1)
    out.uni2 <- getroot_K1_Poi_fast(0, mu = mu, qq=qq, aa=aa,  q = (-1)*Score,qqNB=qqNB, aaNB=aaNB, muNA=muNA,muNB=muNB,NAsigma=NAsigma)
    #cat("out.uni1 out.uni2: ", out.uni1$root, " ", out.uni2$root, "\n")
    #print(out.uni2)
    # print('uni1')
    # print(out.uni1)
    # print('uni2')
    # print(out.uni2)
    #print('dddd')
    #print(Score)
    if (out.uni1$Is.converge == TRUE && out.uni2$Is.converge == TRUE) {
      #   if(T){
      p1<-tryCatch(Get_Saddle_Prob_Poi_fast(out.uni1$root, mu, qq=qq, aa=aa, q=Score,qqNB=qqNB, aaNB=aaNB,muNA,muNB,NAsigma),error=function(e) {return(pval.noadj/2)})
      #print(p1)
      p2<-tryCatch(Get_Saddle_Prob_Poi_fast(out.uni2$root, mu,qq=qq, aa=aa,  q=(-1)*Score,qqNB=qqNB, aaNB=aaNB, muNA,muNB,NAsigma),error=function(e) {return(pval.noadj/2)})	
      #	cat("p1 p2: ", p1, " ", p2, "\n")
      
      pval = abs(p1) + abs(p2)
      Is.converge = TRUE
      # }else if(out.uni1$Is.converge == TRUE){
      #   p1<-tryCatch(Get_Saddle_Prob_Poi_fast(out.uni1$root, mu, qq=qq, aa=aa, q=Score,qqNB=qqNB, aaNB=aaNB,muNA,muNB,NAsigma),error=function(e) {return(pval.noadj/2)})
      #   #print(p1)
      #   pval = abs(p1) 
      #   Is.converge = TRUE
      # }else if(out.uni2$Is.converge == TRUE){
      #   p2<-tryCatch(Get_Saddle_Prob_Poi_fast(out.uni2$root, mu,qq=qq, aa=aa,  q=(-1)*Score,qqNB=qqNB, aaNB=aaNB, muNA,muNB,NAsigma),error=function(e) {return(pval.noadj/2)})	
      #   pval = abs(p1) + abs(p2)
      #   Is.converge = TRUE
      #   
    }else{
      print("Error_Converge")
      pval <- 1
      Is.converge = FALSE
    }
    #}
    return(list(p.value = pval, p.value.NA = pval.noadj,
                Is.converge = Is.converge, Score = Score))
  }
  
  #(obj.noK, G0, mu.a, mu2.a, varRatio)
  
  Score_Test_Sparse_Survival <-function(obj.null, G, mu, mu2, varRatio){
    # mu=mu.a; mu2= mu2.a; G=G0; obj.null=obj.noK
    #tp2a0 = proc.time()
    idx_no0<-which(G>0)
    g1<-G[idx_no0]
    noCov = FALSE
    if(dim(obj.null$X1_fg)[2] == 1){
      noCov = TRUE 
    }
    
    A1<-obj.null$XVX_inv_XV_fg[idx_no0,,drop=F]
    X1_fg<-obj.null$X1_fg[idx_no0,,drop=F]
    mu21<-mu2[idx_no0]
    mu1<-mu[idx_no0]
    y1<-obj.null$y[idx_no0]
    if(length(idx_no0) > 1){
      #    cat("idx_no0 ", idx_no0, "\n")
      #cat("dim(X1) ", dim(X1), "\n")
      #cat("dim(A1) ", dim(A1), "\n")
      Z = t(A1) %*% g1
      B<-X1_fg %*% Z
      #cat("dim(Z) ", dim(Z), "\n")
      #cat("dim(B) ", dim(B), "\n")
      g_tilde1 = g1 - B
      #print(g_tilde1[1:100])
      var2 = t(Z) %*% obj.null$XVX_fg %*% Z - t(B^2) %*% mu21 + t(g_tilde1^2) %*% mu21
      var1 = var2 * varRatio
      S1 = crossprod(y1-mu1, g_tilde1)
      
      if(!noCov){
        S_a2 = obj.null$S_a - colSums(X1_fg * (y1 - mu1))
      }else{
        S_a2 = obj.null$S_a - crossprod(X1_fg, y1 - mu1)
      }
      
      S2 = -S_a2 %*% Z
    }else{
      Z = A1 * g1
      B<- sum(X1_fg * Z)
      g_tilde1 = g1 - B
      var2 = (Z) %*% obj.null$XVX_fg %*% t(Z) - t(B^2) %*% mu21 + t(g_tilde1^2) %*% mu21
      var1 = var2 * varRatio
      S1 = crossprod(y1-mu1, g_tilde1)
      S_a2 = obj.null$S_a - X1_fg * (y1 - mu1)
      S2 = sum(-S_a2 * Z)
    }
    
    S<- S1+S2
    ##meanG = mean(G)
    #S = S - meanG*resq
    ##cat("dim(B): ", dim(B), "\n")
    ##cat("dim(Z): ", dim(Z), "\n")
    ##cat("dim(XWq): ", dim(XWq), "\n")
    ##cat("dim(Wq): ", dim(Wq), "\n")
    #tp2a1 = proc.time()
    #print("tp2a1-tp2a0")
    #print(tp2a1-tp2a0)  
    
    #var1centered1 = t(Z) %*% XWq - t(B) %*% (Wq[idx_no0,])
    #var1centered = var2 - 2*meanG*var1centered1 + meanG^2*qW1  
    #var1 = var1centered * varRatio
    
    #tp2a2 = proc.time()
    #print("tp2a2-tp2a1")
    #print(tp2a2-tp2a1)
    
    
    pval.noadj<-pchisq((S)^2/(var1), lower.tail = FALSE, df=1)
    ##add on 10-25-2017
    BETA = S/var1
    SE = abs(BETA/qnorm(pval.noadj/2))
    Tstat = S
    #tp2a3 = proc.time()
    #print("tp2a3-tp2a2")
    #print(tp2a3-tp2a2)
    #return(c(BETA, SE, Tstat, pval.noadj, pval.noadj, 1, var1, var2))
    return(list(BETA=BETA, SE=SE, Tstat=Tstat, pval.noadj=pval.noadj, pval.noadj=pval.noadj, is.converge=TRUE, var1=var1, var2=var2, B=B, Z=Z, g_tilde1=g_tilde1))	
  }
  
  
  
  
  
  Score_Test_Survival<-function(obj.null, G, mu, mu2, varRatio){
    #G = G - meanG
    g <- G  -  obj.null$XXVX_inv_fg %*%  (obj.null$XV_fg %*% G)
    q<- crossprod(g, obj.null$y) 
    m1<-crossprod(mu, g)
    var2<- crossprod(mu2, g^2)
    #var1 = var2 * varRatio
    #S = (q-m1) + meanG*resq
    S = q-m1
    #meanG = mean(G)
    #S = S - meanG*resq
    #var1centered1 = t(Z) %*% XWq - t(B) %*% Wq
    #var1centered = var2 - 2*meanG*var1centered1 + meanG^2*qW1
    var1 = var2 * varRatio
    
    pval.noadj<-pchisq((S)^2/var1, lower.tail = FALSE, df=1)
    
    ##add on 10-25-2017
    BETA = S/var1
    SE = abs(BETA/qnorm(pval.noadj/2))
    #Tstat = S^2
    Tstat = S
    
    #return(c(BETA, SE, Tstat, pval.noadj, pval.noadj, NA, var1, var2))
    #return(c(pval.noadj, pval.noadj, TRUE, var1, var2))
    return(list(BETA=BETA, SE=SE, Tstat=Tstat, pval.noadj=pval.noadj, pval.noadj=pval.noadj, is.converge=TRUE, var1=var1, var2=var2,  g_tilde=g))
  }
  
  
  
  
  
  
  innerProduct <- function(x, y) {
    if (length(x) != length(y)) {
      stop("Vectors must be of the same length")
    }
    return(sum(x * y))
  }
  
  
  
  
  
  Saddle_Prob_Poisson=function (Score, pval.noadj, mu,qq, aa, Cutoff = 2, alpha = 5*10^-8, m1, var1){
    #m1 <- sum(mu * g)
    #var1 <- sum(mu * g^2)
    p1 = NULL
    p2 = NULL
    #cat("Score is ", Score, "\n")
    #print(g[1:20])
    
    
    #Score <- q - m1
    #qinv = -sign(q - m1) * abs(q - m1) + m1
    #pval.noadj <- pchisq((q - m1)^2/var1, lower.tail = FALSE,
    #    df = 1)
    Is.converge = TRUE
    
    #if (abs(q - m1)/sqrt(var1) < Cutoff) {
    #    pval = pval.noadj
    #}else {
    #        print("Saddle_Prob_Poisson")
    #t1 <- proc.time()
    out.uni1 <- getroot_K1_Poi(0, mu = mu, aa=aa, qq=qq, q = Score)
    out.uni2 <- getroot_K1_Poi(0, mu = mu, aa=aa, qq=qq, q = (-1)*Score)
    # print('uni1')
    # print(out.uni1)
    # print('uni2')
    # print(out.uni2)
    #t2 <- proc.time()
    if (out.uni1$Is.converge == TRUE && out.uni2$Is.converge == TRUE) {
      p1 <- tryCatch(Get_Saddle_Prob_Poi(out.uni1$root, mu,qq=qq,aa=aa,  q=Score), error=function(e) {return(pval.noadj/2)})	
      p2 <- tryCatch(Get_Saddle_Prob_Poi(out.uni2$root, mu,qq=qq,aa=aa, q = (-1)*Score), error=function(e) {return(pval.noadj/2)})	
      #p1 <- Get_Saddle_Prob_Poi(out.uni1$root, mu, g, q)
      #p2 <- Get_Saddle_Prob_Poi(out.uni2$root, mu, g, qinv)
      #	    cat("p1 p2: ", p1, " ", p2, "\n")	
      
      pval = abs(p1) + abs(p2)
      Is.converge = TRUE
    }
    else {
      print("Error_Converge")
      #pval <- pval.noadj
      pval <- 1
      Is.converge = FALSE
    }
    # t3 <- proc.time()
    # print(t2-t1)
    # print(t3-t2)
    
    #}
    return(list(p.value = pval, p.value.NA = pval.noadj,
                Is.converge = Is.converge, Score = Score))
  }
  
  
  if(T){
    ##saddlepoint approxmation for sum of weighted Poisson distribution
    Korg_Poi<-function(t, mu, qq, aa)
    {
      n.t<-length(t)
      out<-rep(0,n.t)
      
      for(i in 1:n.t){
        t1<-t[i]
        temp<-mu*(exp(qq*t1)  - 1) + aa*t1
        out[i]<-sum(temp)
      }
      return(out)
    }
    
    K1_Poi<-function(t, mu, qq, aa)
    {
      n.t<-length(t)
      out<-rep(0,n.t)
      
      for(i in 1:n.t){
        t1<-t[i]
        temp<-mu * qq * exp(qq*t1) + aa
        out[i]<-sum(temp)
      }
      return(out)
    }
    
    
    K1_adj_Poi<-function(t, mu, qq, aa, q)
    {
      n.t<-length(t)	
      out<-rep(0,n.t)
      
      for(i in 1:n.t){
        t1<-t[i]
        temp<-mu * qq * exp(qq*t1)  + aa
        out[i]<-sum(temp)-q
      }
      return(out)
    }
    
    
    K2_Poi<-function(t, mu, qq,aa)
    {
      n.t<-length(t)
      out<-rep(0,n.t)
      
      for(i in 1:n.t){
        t1<-t[i]
        temp<-mu * qq^2 * exp(qq*t1)
        out[i]<-sum(temp, na.rm=TRUE)
      }
      return(out)
    }
    
    
    getroot_K1_Poi<-function(init,mu,qq,aa,q,m1,tol=.Machine$double.eps^0.25,maxiter=1000)
    {
      t<-init
      K1_eval<-K1_adj_Poi(t,mu,qq,aa,q)
      #cat("K1_eval ", K1_eval, "\n")
      prevJump<- Inf
      rep<-1
      repeat
      {
        t1 <- proc.time()
        K2_eval<-K2_Poi(t,mu,qq,aa)
        tnew<-t-K1_eval/K2_eval
        if(is.na(tnew))
        {
          conv=FALSE
          break
        }
        if(abs(tnew-t)<tol)
        {
          conv<-TRUE
          break
        }
        if(rep==maxiter)
        {
          conv<-FALSE
          break
        }
        # t2 <- proc.time()
        newK1<-K1_adj_Poi(tnew,mu,qq,aa,q)
        if(sign(K1_eval)!=sign(newK1))
        {
          if(abs(tnew-t)>prevJump-tol)
          {
            tnew<-t+sign(newK1-K1_eval)*prevJump/2
            newK1<-K1_adj_Poi(tnew,mu,qq,aa,q)
            prevJump<-prevJump/2
          } else {
            prevJump<-abs(tnew-t)
          }
        }
        #t3 <- proc.time()
        # print('ttt')
        # print(t3-t2)
        # print(t2-t1)
        rep<-rep+1
        t<-tnew
        K1_eval<-newK1
      } 
      # print('repeat')
      # print(rep)
      # print('root')
      # print(t)
      return(list(root=t,n.iter=rep,Is.converge=conv))
      # }
    }
    
    
    Get_Saddle_Prob_Poi<-function(zeta, mu, qq,aa, q) 
    {
      k1<-Korg_Poi(zeta, mu, qq,aa)
      #cat("k1 is ", k1, "\n")
      k2<-K2_Poi(zeta, mu, qq,aa)
      #cat("k2 is ", k2, "\n")
      if(is.finite(k1) && is.finite(k2))
      {
        temp1<-zeta * q - k1
        
        
        w<-sign(zeta) * (2 *temp1)^{1/2}
        v<- zeta * (k2)^{1/2}
        
        Z.test<-w + 1/w * log(v/w)	
        
        if(Z.test > 0){
          pval<-pnorm(Z.test, lower.tail = FALSE)
        } else {
          pval= -pnorm(Z.test, lower.tail = TRUE)
        }	
      } else {
        pval<-0
      }
      
      return(pval)
    }
  }
  
  if(T){
    ##saddlepoint approxmation for sum of weighted Poisson distribution
    Korg_Poi_fast <- function(t, mu, qq, aa, muNA,muNB,NAsigma, qqNB, aaNB)
    {
      n.t<-length(t)
      out<-rep(0,n.t)
      
      for(i in 1:n.t){
        t1<-t[i]
        temp<-muNB*(exp(qqNB*t1)  - 1) + aaNB*t1
        #out[i]<-sum(temp)+NAmu*t1+0.5*NAsigma*t1^2
        out[i]<-sum(temp)+0.5*NAsigma*t1^2 
      }
      return(out)
    }
    
    
    
    K1_adj_Poi_fast <-function(t, mu, qq, aa, q, muNA,muNB,NAsigma, qqNB, aaNB)
    {
      n.t<-length(t)	
      out<-rep(0,n.t)
      
      for(i in 1:n.t){
        t1<-t[i]
        temp<-muNB * qqNB * exp(qqNB*t1) + aaNB
        #temp2<-NAmu+NAsigma*t1
        temp2<-NAsigma*t1
        out[i]<-sum(temp)-q + temp2
      }
      return(out)
    }
    
    
    K2_Poi_fast<-function(t, mu, qq, aa, muNA,muNB,NAsigma, qqNB, aaNB)
    {
      n.t<-length(t)
      out<-rep(0,n.t)
      
      for(i in 1:n.t){
        t1<-t[i]
        temp<-muNB * qqNB^2 * exp(qqNB*t1)
        out[i]<-sum(temp, na.rm=TRUE) + NAsigma
      }
      return(out)
    }
    
    
    
    getroot_K1_Poi_fast<-function(init,mu,qq,aa, q,m1, qqNB, aaNB,muNA,muNB,NAsigma,
                                  tol=.Machine$double.eps^0.25,maxiter=1000)
    {
      t<-init
      # print('bbbb')
      # print(q)
      # print(K1_adj_Poi_fast(t,mu,qq,aa, q, qqNB=qqNB, aaNB=aaNB,NAsigma=NAsigma, muNA=muNA,muNB=muNB))
      # print(K1_adj_Poi_fast(t,mu,qq,-aa, q, qqNB=qqNB, aaNB=-aaNB,NAsigma=NAsigma, muNA=muNA,muNB=muNB))
      # print(K1_adj_Poi_fast(t,mu,-qq,aa, q, qqNB=-qqNB, aaNB=aaNB,NAsigma=NAsigma, muNA=muNA,muNB=muNB))
      # print(K1_adj_Poi_fast(t,mu,-qq,-aa, q, qqNB=-qqNB, aaNB=-aaNB,NAsigma=NAsigma, muNA=muNA,muNB=muNB))
      K1_eval<-K1_adj_Poi_fast(t,mu,qq,aa, q, qqNB=qqNB, aaNB=aaNB,NAsigma=NAsigma, muNA=muNA,muNB=muNB)
      #cat("K1_eval: ", K1_eval, "\n")
      prevJump<- Inf
      rep<-1
      repeat
      {
        K2_eval<-K2_Poi_fast(t,mu,qq,aa, qqNB=qqNB, aaNB=aaNB,NAsigma=NAsigma, muNA=muNA,muNB=muNB)
        tnew<- t-K1_eval/K2_eval
        if(is.na(tnew))
        {
          conv=FALSE
          break
        }
        if(abs(tnew-t)<tol)
        {
          conv<-TRUE
          break
        }
        if(rep==maxiter)
        {
          conv<-FALSE
          break
        }
        
        newK1<-K1_adj_Poi_fast(tnew,mu,qq,aa,q,muNA,muNB,NAsigma, qqNB, aaNB)
        if(sign(K1_eval)!=sign(newK1))
        {
          if(abs(tnew-t)>prevJump-tol)
          {
            tnew<-t+sign(newK1-K1_eval)*prevJump/2
            newK1<-K1_adj_Poi_fast(tnew,mu,qq,aa,q,muNA,muNB,NAsigma, qqNB, aaNB)
            prevJump<-prevJump/2
          } else {
            prevJump<-abs(tnew-t)
          }
        }
        
        rep<-rep+1
        t<-tnew
        K1_eval<-newK1
      } 
      return(list(root=t,n.iter=rep,Is.converge=conv))
      # }
    }
    
    
    Get_Saddle_Prob_Poi_fast<-function(zeta, mu, qq,aa, q,muNA,muNB,NAsigma, qqNB, aaNB) 
    {
      
      
      k1<-Korg_Poi_fast(zeta, mu, qq,aa,muNA,muNB,NAsigma, qqNB, aaNB)
      #cat("k1 is ", k1, "\n")
      k2<-K2_Poi_fast(zeta, mu,qq,aa,muNA,muNB,NAsigma, qqNB, aaNB)
      #cat("k2 is ", k2, "\n")
      if(is.finite(k1) && is.finite(k2))
      {
        temp1<-zeta * q - k1
        
        
        w<-sign(zeta) * (2 *temp1)^{1/2}
        v<- zeta * (k2)^{1/2}
        
        Z.test<-w + 1/w * log(v/w)	
        
        if(Z.test > 0){
          pval<-pnorm(Z.test, lower.tail = FALSE)
        } else {
          pval= -pnorm(Z.test, lower.tail = TRUE)
        }	
      } else {
        pval<-0
      }
      
      return(pval)
    }
  }
  
  
  
  
  
  
  
  
  if(T){
    ##saddlepoint approxmation for sum of weighted Poisson distribution
    Korg_Poi_fast_skat <- function(t, mu, qq, aa, wvar=0, wvar2=0)
    {
      NBset <- which(qq!=0)
      n.t<-length(t)
      out<-rep(0,n.t)
      
      for(i in 1:n.t){
        t1<-t[i]
        temp<-mu[NBset]*(exp(qq[NBset]*t1)  - 1) + aa[NBset]*t1
        #out[i]<-sum(temp)+NAmu*t1+0.5*NAsigma*t1^2
        out[i]<-sum(temp)+0.5*wvar*t1^2+wvar2*(exp(t1)-1)-wvar2*t1
        #+0.5*NAsigma*t1^2 
      }
      return(out)
    }
    
    #K1_Poi<-function(t, mu, g)
    #{
    #  n.t<-length(t)
    #  out<-rep(0,n.t)
    
    #  for(i in 1:n.t){
    #    t1<-t[i]
    #    temp<-mu * g * exp(g*t1) - mu * g
    #    out[i]<-sum(temp)
    #  }
    #  return(out)
    #}
    
    # tnew,mu,qq,aa,q,wvar=wvar, wvar2=wvar2
    K1_adj_Poi_fast_skat <-function(t, mu, qq, aa, q, wvar=0, wvar2=0)
    {
      n.t<-length(t)	
      out<-rep(0,n.t)
      NBset <- which(qq!=0)
      
      for(i in 1:n.t){
        t1<-t[i]
        temp<- mu[NBset] * qq[NBset] * exp(qq[NBset]*t1) + aa[NBset]
        #temp2<-NAmu+NAsigma*t1
        #temp2<-NAsigma*t1
        #out[i]<-sum(temp)-q + temp2
        if(wvar2==0){
          out[i] <- sum(temp)-q+wvar*t1
        }else{
          out[i] <- sum(temp)-q+wvar*t1+wvar2*exp(t1)-wvar2
        }
      }
      return(out)
    }
    
    
    K2_Poi_fast_skat<-function(t, mu, qq, wvar=0, wvar2=0)
    {
      n.t<-length(t)
      out<-rep(0,n.t)
      NBset <- which(qq!=0)
      
      for(i in 1:n.t){
        t1<-t[i]
        temp<-mu[NBset] * qq[NBset]^2 * exp(qq[NBset]*t1)
        out[i]<-sum(temp, na.rm=TRUE) +wvar+wvar2*exp(t1)
      }
      return(out)
    }
    
    
    getroot_K1_Poi_fast_skat<-function(init,mu,qq,aa, q,
                                       tol=.Machine$double.eps^0.25,
                                       maxiter=1000, wvar=0, wvar2=0)
    {
      t<-init
      # print('bbbb')
      # print(q)
      # print(K1_adj_Poi_fast(t,mu,qq,aa, q, qqNB=qqNB, aaNB=aaNB,NAsigma=NAsigma, muNA=muNA,muNB=muNB))
      # print(K1_adj_Poi_fast(t,mu,qq,-aa, q, qqNB=qqNB, aaNB=-aaNB,NAsigma=NAsigma, muNA=muNA,muNB=muNB))
      # print(K1_adj_Poi_fast(t,mu,-qq,aa, q, qqNB=-qqNB, aaNB=aaNB,NAsigma=NAsigma, muNA=muNA,muNB=muNB))
      # print(K1_adj_Poi_fast(t,mu,-qq,-aa, q, qqNB=-qqNB, aaNB=-aaNB,NAsigma=NAsigma, muNA=muNA,muNB=muNB))
      K1_eval<-K1_adj_Poi_fast_skat(t,mu,qq,aa, q,wvar=wvar, wvar2=wvar2)
      #cat("K1_eval: ", K1_eval, "\n")
      prevJump<- Inf
      rep<-1
      repeat
      {
        K2_eval<-K2_Poi_fast_skat(t,mu,qq,wvar=wvar, wvar2=wvar2)
        tnew<- t-K1_eval/K2_eval
        if(is.na(tnew))
        {
          conv=FALSE
          break
        }
        if(abs(tnew)==Inf){
          conv=FALSE
          break
        }
        if(abs(tnew-t)<tol)
        {
          conv<-TRUE
          break
        }
        if(rep==maxiter)
        {
          conv<-FALSE
          break
        }
        # print('tnew')
        # print(tnew)
        newK1<-K1_adj_Poi_fast_skat(tnew,mu,qq,aa,q,wvar=wvar, wvar2=wvar2)
        # print('newK1')
        # print(newK1)
        # print('K1_eval')
        # print(K1_eval)
        if(sign(K1_eval)!=sign(newK1))
        {
          if(abs(tnew-t)>prevJump-tol)
          {
            tnew<-t+sign(newK1-K1_eval)*prevJump/2
            newK1<-K1_adj_Poi_fast_skat(tnew,mu,qq,aa,q,wvar=wvar, wvar2=wvar2)
            prevJump<-prevJump/2
          } else {
            prevJump<-abs(tnew-t)
          }
        }
        
        rep<-rep+1
        t<-tnew
        K1_eval<-newK1
      } 
      return(list(root=t,n.iter=rep,Is.converge=conv))
      # }
    }
    
    
    Get_Saddle_Prob_Poi_fast_skat<-function(zeta, mu, qq,aa, q, wvar=0, wvar2=0) 
    {
      
      
      k1<-Korg_Poi_fast_skat(zeta, mu, qq,aa,wvar=wvar, wvar2=wvar2)
      #cat("k1 is ", k1, "\n")
      k2<-K2_Poi_fast_skat(zeta, mu,qq,wvar=wvar, wvar2=wvar2)
      #cat("k2 is ", k2, "\n")
      if(is.finite(k1) && is.finite(k2))
      {
        temp1<-zeta * q - k1
        
        
        w<-sign(zeta) * (2 *temp1)^{1/2}
        v<- zeta * (k2)^{1/2}
        
        Z.test<- w + 1/w * log(v/w)	
        # print('Z.test')
        # print(Z.test)
        #print(Z.test)
        if(Z.test > 0){
          pval<-pnorm(Z.test, lower.tail = FALSE)
        } else {
          pval= -pnorm(Z.test, lower.tail = TRUE)
        }	
      } else {
        pval<-0
      }
      
      return(pval)
    }
  }
  
  
  
  Saddle_Prob_Poisson_fast_skat=function (Score, pval.noadj, mu, qq,aa,wvar=0, wvar2=0){
    #m1 <- sum(mu * g)
    #var1 <- sum(mu * g^2)
    p1 = NULL
    p2 = NULL
    
    out.uni1 <- getroot_K1_Poi_fast_skat(0, mu = mu, qq=qq, aa=aa,  q = Score,wvar=wvar, wvar2=wvar2)
    # print(out.uni1)
    # out.uni2 <- getroot_K1_Poi_fast(0, mu = mu, qq=qq, aa=aa,  q = (-1)*Score)
    #cat("out.uni1 out.uni2: ", out.uni1$root, " ", out.uni2$root, "\n")
    # print(out.uni2)
    
    # print('dddd')
    # print(Score)
    if(abs(out.uni1$root)>1000){
      out.uni1$Is.converge <- FALSE
    }
    if (out.uni1$Is.converge == TRUE) {
      #if(T){
      p1<-tryCatch(Get_Saddle_Prob_Poi_fast_skat(out.uni1$root, mu=mu, qq=qq, aa=aa, q=Score,wvar=wvar, wvar2=wvar2),error=function(e) {return(pval.noadj)})
      #print(p1)
      # p2<-tryCatch(Get_Saddle_Prob_Poi_fast(out.uni2$root, mu,qq=qq, aa=aa,  q=(-1)*Score),error=function(e) {return(pval.noadj/2)})	
      # cat("p1 p2: ", p1, " ", p2, "\n")
      # 
      pval = abs(p1) 
      Is.converge = TRUE
      # }else if(out.uni1$Is.converge == TRUE){
      #   p1<-tryCatch(Get_Saddle_Prob_Poi_fast(out.uni1$root, mu, qq=qq, aa=aa, q=Score,qqNB=qqNB, aaNB=aaNB,muNA,muNB,NAsigma),error=function(e) {return(pval.noadj/2)})
      #   #print(p1)
      #   pval = abs(p1)
      #   Is.converge = TRUE
      # }else if(out.uni2$Is.converge == TRUE){
      #   p2<-tryCatch(Get_Saddle_Prob_Poi_fast(out.uni2$root, mu,qq=qq, aa=aa,  q=(-1)*Score,qqNB=qqNB, aaNB=aaNB, muNA,muNB,NAsigma),error=function(e) {return(pval.noadj/2)})
      #   pval = abs(p1) + abs(p2)
      #   Is.converge = TRUE
    }else{
      print("Error_Converge")
      #pval <- pval.noadj
      pval <- 1
      Is.converge = FALSE
    }
    #}
    return(list(p.value = pval, p.value.NA = pval.noadj,
                Is.converge = Is.converge, Score = Score))
  }
  
  
  
  
}
if(T){
  quiet <- function(expr) {
    tf  <- tempfile()                # 日志文件（想保存日志就改成固定文件名）
    con <- file(tf, open = "wt")     # 显式打开 text connection
    
    ## 把两条流都重定向到同一个 connection
    sink(con)                        # 标准输出（print/cat 等）
    sink(con, type = "message")      # 消息流（message/warning/error）
    
    on.exit({
      sink(type = "message")         # 先恢复 message 连接
      sink()                         # 再恢复 stdout 连接
      close(con)                     # 关闭文件句柄
    }, add = TRUE)
    
    suppressWarnings( suppressMessages( force(expr) ) )
    invisible(NULL)
  }
  
  
  noncoding_survival <- function(chr,gene_name,genofile,obj_nullmodel,
                                 rare_maf_cutoff=0.01,rv_num_cutoff=2,rv_num_cutoff_min_prefilter=2,
                                 rv_num_cutoff_max_prefilter=1e9,
                                 QC_label="annotation/filter",variant_type=c("SNV","Indel","variant"),geno_missing_imputation=c("mean","minor"),
                                 Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,
                                 Use_annotation_weights=c(TRUE,FALSE),Annotation_name=NULL,
                                 SPA_p_filter=FALSE,p_filter_cutoff=0.05,
                                 silent=FALSE, varRatio=1){
    
    ## evaluate choices
    variant_type <- match.arg(variant_type)
    geno_missing_imputation <- match.arg(geno_missing_imputation)
    
    phenotype.id <- as.character(obj_nullmodel$id_include)
    n_pheno <- obj_nullmodel$n.pheno
    
    ## SPA status
    # if(!is.null(obj_nullmodel$use_SPA))
    # {
    #   use_SPA <- obj_nullmodel$use_SPA
    # }else
    # {
    #   use_SPA <- FALSE
    # }
    # 
    use_SPA <- TRUE
    use_SPA0 <- F
    #####################################
    #   Gene Info
    ## get SNV id
    filter <- seqGetData(genofile, QC_label)
    if(variant_type=="variant")
    {
      SNVlist <- filter == "PASS"
    }
    
    if(variant_type=="SNV")
    {
      SNVlist <- (filter == "PASS") & isSNV(genofile)
    }
    
    if(variant_type=="Indel")
    {
      SNVlist <- (filter == "PASS") & (!isSNV(genofile))
    }
    
    variant.id <- seqGetData(genofile, "variant.id")
    
    rm(filter)
    gc()
    
    # ww <- length(id_all)
    # if(ww<7){
    #   for(ll in (ww+1):7){
    #     id_all[[ll]] <- integer(0)
    #   }
    # }
    # for(ll in 1:7){
    #   if(length(id_all[[ll]])==0){
    #     id_all[[ll]] <- integer(0)  ## if it is NULL, filter will select all SNPs
    #   }
    # }
    ## genotype id
    
    # print(length(id_all[[1]])>=rv_num_cutoff_min_prefilter & length(id_all[[1]])<rv_num_cutoff_max_prefilter)
    # if(length(id_all[[1]])>=rv_num_cutoff_min_prefilter & length(id_all[[1]])<rv_num_cutoff_max_prefilter){
    
    GENCODE.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Category")]))
    is.in <- (GENCODE.Category=="downstream")&(SNVlist)
    variant.id.downstream <- variant.id[is.in]
    
    seqSetFilter(genofile,variant.id=variant.id.downstream,sample.id=phenotype.id)
    
    rm(variant.id.downstream)
    gc()
    
    GENCODE.Info <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
    GENCODE.Info.split <- strsplit(GENCODE.Info, split = "[,]")
    variant_gene_num <- sapply(GENCODE.Info.split,function(z) length(z))
    
    variant.id.SNV <- seqGetData(genofile, "variant.id")
    variant.id.SNV <- rep(variant.id.SNV,variant_gene_num)
    
    rm(GENCODE.Info)
    gc()
    
    rm(variant_gene_num)
    gc()
    
    Gene <- as.character(unlist(GENCODE.Info.split))
    
    rm(GENCODE.Info.split)
    gc()
    
    seqResetFilter(genofile)
    
    ### Gene
    is.in <- which(Gene==gene_name)
    variant.is.in <- variant.id.SNV[is.in]
    
    seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)
    
    ## genotype id
    id.genotype <- seqGetData(genofile,"sample.id")
    # id.genotype.match <- rep(0,length(id.genotype))
    
    id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
    phenotype.id.merge <- data.frame(phenotype.id)
    phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
    id.genotype.match <- phenotype.id.merge$index
    
    ## Genotype
    Geno <- NULL
    if(length(seqGetData(genofile, "variant.id"))<rv_num_cutoff_max_prefilter )
    {
      Geno <- seqGetData(genofile, "$dosage")
      Geno <- Geno[id.genotype.match,,drop=FALSE]
    }
    
    ## impute missing
    if(!is.null(dim(Geno)))
    {
      if(dim(Geno)[2]>0)
      {
        if(geno_missing_imputation=="mean")
        {
          tmp <- matrix_flip_mean(Geno)
          Geno <- tmp$Geno
          MAF <- tmp$MAF
          rm(tmp)
        }
        if(geno_missing_imputation=="minor")
        {
          Geno <- matrix_flip_minor(Geno)$Geno
        }
      }
    }
    
    
    ## Annotation
    Anno.Int.PHRED.sub <- NULL
    Anno.Int.PHRED.sub.name <- NULL
    
    if(variant_type=="SNV")
    {
      if(Use_annotation_weights)
      {
        for(k in 1:length(Annotation_name))
        {
          if(Annotation_name[k]%in%Annotation_name_catalog$name)
          {
            Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,Annotation_name[k])
            Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name==Annotation_name[k])]))
            
            if(Annotation_name[k]=="CADD")
            {
              Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
            }
            
            if(Annotation_name[k]=="aPC.LocalDiversity")
            {
              Annotation.PHRED.2 <- -10*log10(1-10^(-Annotation.PHRED/10))
              Annotation.PHRED <- cbind(Annotation.PHRED,Annotation.PHRED.2)
              Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,paste0(Annotation_name[k],"(-)"))
            }
            Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub,Annotation.PHRED)
          }
        }
        
        Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
        colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
      }
    }
    
    pvalues <- 1
    if(n_pheno == 1)
    {
      if(F)
      {
        try(pvalues <- STAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff),silent=silent)
      }else
      {
        try(pvalues <- STAAR_Binary_SPA_survival(Geno,obj_nullmodel,Anno.Int.PHRED.sub,
                                                 rare_maf_cutoff=rare_maf_cutoff,
                                                 rv_num_cutoff=rv_num_cutoff,
                                                 varRatio=varRatio,MAF=MAF,
                                                 #rv_num_cutoff_max=rv_num_cutoff_max,
                                                 SPA_p_filter=SPA_p_filter,
                                                 p_filter_cutoff=p_filter_cutoff),silent=silent)
      }
    }else
    {
      try(pvalues <- MultiSTAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff),silent=silent)
    }
    
    results_downstream <- c()
    if(inherits(pvalues, "list"))
    {
      results_temp <- rep(NA,4)
      results_temp[3] <- "downstream"
      results_temp[2] <- chr
      results_temp[1] <- as.character(gene_name)
      results_temp[4] <- pvalues$num_variant
      
      if(!use_SPA0)
      {
        results_temp <- c(results_temp,pvalues$cMAC,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
                          pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
                          pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)
      }else
      {
        results_temp <- c(results_temp,pvalues$cMAC,
                          pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_B)
      }
      
      results_downstream <- rbind(results_downstream,results_temp)
    }
    
    if(!is.null(results_downstream))
    {
      if(!use_SPA0)
      {
        colnames(results_downstream) <- colnames(results_downstream, do.NULL = FALSE, prefix = "col")
        colnames(results_downstream)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
        colnames(results_downstream)[(dim(results_downstream)[2]-1):dim(results_downstream)[2]] <- c("ACAT-O","STAAR-O")
      }else
      {
        colnames(results_downstream) <- colnames(results_downstream, do.NULL = FALSE, prefix = "col")
        colnames(results_downstream)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
        colnames(results_downstream)[dim(results_downstream)[2]] <- c("STAAR-B")
      }
    }
    # }else{
    #   seqSetFilter(genofile,sample.id=phenotype.id)
    #   id.genotype <- seqGetData(genofile,"sample.id")
    #   # id.genotype.match <- rep(0,length(id.genotype))
    #   
    #   id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
    #   phenotype.id.merge <- data.frame(phenotype.id)
    #   phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
    #   id.genotype.match <- phenotype.id.merge$index
    #   
    #   results_downstream <- c()
    # }
    
    #seqResetFilter(genofile)
    
    ########################################
    #   Upstream
    
    
    #if(length(id_all[[2]])>=rv_num_cutoff_min_prefilter & length(id_all[[2]])<rv_num_cutoff_max_prefilter){
    
    is.in <- (GENCODE.Category=="upstream")&(SNVlist)
    variant.id.upstream <- variant.id[is.in]
    
    seqSetFilter(genofile,variant.id=variant.id.upstream,sample.id=phenotype.id)
    
    rm(variant.id.upstream)
    gc()
    
    GENCODE.Info <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
    GENCODE.Info.split <- strsplit(GENCODE.Info, split = "[,]")
    variant_gene_num <- sapply(GENCODE.Info.split,function(z) length(z))
    
    variant.id.SNV <- seqGetData(genofile, "variant.id")
    variant.id.SNV <- rep(variant.id.SNV,variant_gene_num)
    
    rm(GENCODE.Info)
    gc()
    
    rm(variant_gene_num)
    gc()
    
    Gene <- as.character(unlist(GENCODE.Info.split))
    
    rm(GENCODE.Info.split)
    gc()
    
    seqResetFilter(genofile)
    
    ### Gene
    is.in <- which(Gene==gene_name)
    variant.is.in <- variant.id.SNV[is.in]
    
    seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)
    
    ## genotype id
    id.genotype <- seqGetData(genofile,"sample.id")
    # id.genotype.match <- rep(0,length(id.genotype))
    
    id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
    phenotype.id.merge <- data.frame(phenotype.id)
    phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
    id.genotype.match <- phenotype.id.merge$index
    
    ## Genotype
    Geno <- NULL
    if(length(seqGetData(genofile, "variant.id"))<rv_num_cutoff_max_prefilter)
    {
      Geno <- seqGetData(genofile, "$dosage")
      Geno <- Geno[id.genotype.match,,drop=FALSE]
    }
    
    ## impute missing
    if(!is.null(dim(Geno)))
    {
      if(dim(Geno)[2]>0)
      {
        if(geno_missing_imputation=="mean")
        {
          tmp <- matrix_flip_mean(Geno)
          Geno <- tmp$Geno
          MAF <- tmp$MAF
          rm(tmp)
        }
        if(geno_missing_imputation=="minor")
        {
          Geno <- matrix_flip_minor(Geno)$Geno
        }
      }
    }
    
    ## Annotation
    Anno.Int.PHRED.sub <- NULL
    Anno.Int.PHRED.sub.name <- NULL
    
    if(variant_type=="SNV")
    {
      if(Use_annotation_weights)
      {
        for(k in 1:length(Annotation_name))
        {
          if(Annotation_name[k]%in%Annotation_name_catalog$name)
          {
            Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,Annotation_name[k])
            Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name==Annotation_name[k])]))
            
            if(Annotation_name[k]=="CADD")
            {
              Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
            }
            
            if(Annotation_name[k]=="aPC.LocalDiversity")
            {
              Annotation.PHRED.2 <- -10*log10(1-10^(-Annotation.PHRED/10))
              Annotation.PHRED <- cbind(Annotation.PHRED,Annotation.PHRED.2)
              Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,paste0(Annotation_name[k],"(-)"))
            }
            Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub,Annotation.PHRED)
          }
        }
        
        Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
        colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
      }
    }
    
    pvalues <- 1
    if(n_pheno == 1)
    {
      if(F)
      {
        try(pvalues <- STAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,),silent=silent)
      }else
      {
        try(pvalues <- STAAR_Binary_SPA_survival(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff,varRatio=varRatio,MAF=MAF),silent=silent)
      }
    }else
    {
      try(pvalues <- MultiSTAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff),silent=silent)
    }
    
    results_upstream <- c()
    if(inherits(pvalues, "list"))
    {
      results_temp <- rep(NA,4)
      results_temp[3] <- "upstream"
      results_temp[2] <- chr
      results_temp[1] <- as.character(gene_name)
      results_temp[4] <- pvalues$num_variant
      
      if(!use_SPA0)
      {
        results_temp <- c(results_temp,pvalues$cMAC,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
                          pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
                          pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)
      }else
      {
        results_temp <- c(results_temp,pvalues$cMAC,
                          pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_B)
      }
      
      results_upstream <- rbind(results_upstream,results_temp)
    }
    
    if(!is.null(results_upstream))
    {
      if(!use_SPA0)
      {
        colnames(results_upstream) <- colnames(results_upstream, do.NULL = FALSE, prefix = "col")
        colnames(results_upstream)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
        colnames(results_upstream)[(dim(results_upstream)[2]-1):dim(results_upstream)[2]] <- c("ACAT-O","STAAR-O")
      }else
      {
        colnames(results_upstream) <- colnames(results_upstream, do.NULL = FALSE, prefix = "col")
        colnames(results_upstream)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
        colnames(results_upstream)[dim(results_upstream)[2]] <- c("STAAR-B")
      }
    }
    # }else{
    #   results_upstream <- c()
    # }
    
    #seqResetFilter(genofile)
    
    ########################################################
    #                UTR
    is.in <- ((GENCODE.Category=="UTR3")|(GENCODE.Category=="UTR5")|(GENCODE.Category=="UTR5;UTR3"))&(SNVlist)
    variant.id.UTR <- variant.id[is.in]
    
    rm(GENCODE.Category)
    gc()
    
    seqSetFilter(genofile,variant.id=variant.id.UTR,sample.id=phenotype.id)
    
    rm(variant.id.UTR)
    gc()
    
    GENCODE.Info <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
    GENCODE.Info.split <- strsplit(GENCODE.Info, split = "[(]")
    
    rm(GENCODE.Info)
    gc()
    
    Gene <- as.character(sapply(GENCODE.Info.split,function(z) z[1]))
    
    rm(GENCODE.Info.split)
    gc()
    
    variant.id.SNV <- seqGetData(genofile, "variant.id")
    
    seqResetFilter(genofile)
    
    ### Gene
    is.in <- which(Gene==gene_name)
    variant.is.in <- variant.id.SNV[is.in]
    
    seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)
    
    ## genotype id
    id.genotype <- seqGetData(genofile,"sample.id")
    # id.genotype.match <- rep(0,length(id.genotype))
    
    id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
    phenotype.id.merge <- data.frame(phenotype.id)
    phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
    id.genotype.match <- phenotype.id.merge$index
    
    ## Genotype
    Geno <- NULL
    if(length(seqGetData(genofile, "variant.id"))<rv_num_cutoff_max_prefilter)
    {
      Geno <- seqGetData(genofile, "$dosage")
      Geno <- Geno[id.genotype.match,,drop=FALSE]
    }
    
    ## impute missing
    if(!is.null(dim(Geno)))
    {
      if(dim(Geno)[2]>0)
      {
        if(geno_missing_imputation=="mean")
        {
          tmp <- matrix_flip_mean(Geno)
          Geno <- tmp$Geno
          MAF <- tmp$MAF
          rm(tmp)
        }
        if(geno_missing_imputation=="minor")
        {
          Geno <- matrix_flip_minor(Geno)$Geno
        }
      }
    }
    
    #gc()
    ## Annotation
    Anno.Int.PHRED.sub <- NULL
    Anno.Int.PHRED.sub.name <- NULL
    
    if(variant_type=="SNV")
    {
      if(Use_annotation_weights)
      {
        for(k in 1:length(Annotation_name))
        {
          if(Annotation_name[k]%in%Annotation_name_catalog$name)
          {
            Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,Annotation_name[k])
            Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name==Annotation_name[k])]))
            
            if(Annotation_name[k]=="CADD")
            {
              Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
            }
            
            if(Annotation_name[k]=="aPC.LocalDiversity")
            {
              Annotation.PHRED.2 <- -10*log10(1-10^(-Annotation.PHRED/10))
              Annotation.PHRED <- cbind(Annotation.PHRED,Annotation.PHRED.2)
              Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,paste0(Annotation_name[k],"(-)"))
            }
            Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub,Annotation.PHRED)
          }
        }
        
        Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
        colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
      }
    }
    
    pvalues <- 1
    if(n_pheno == 1)
    {
      if(F)
      {
        try(pvalues <- STAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff),silent=silent)
      }else
      {
        try(pvalues <- STAAR_Binary_SPA_survival(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff, varRatio=varRatio, MAF=MAF),silent=silent)
      }
    }else
    {
      try(pvalues <- MultiSTAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff),silent=silent)
    }
    
    results_UTR <- c()
    if(inherits(pvalues, "list"))
    {
      results_temp <- rep(NA,4)
      results_temp[3] <- "UTR"
      results_temp[2] <- chr
      results_temp[1] <- as.character(gene_name)
      results_temp[4] <- pvalues$num_variant
      
      if(!use_SPA0)
      {
        results_temp <- c(results_temp,pvalues$cMAC,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
                          pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
                          pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)
      }else
      {
        results_temp <- c(results_temp,pvalues$cMAC,
                          pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_B)
      }
      
      results_UTR <- rbind(results_UTR,results_temp)
    }
    
    if(!is.null(results_UTR))
    {
      if(!use_SPA0)
      {
        colnames(results_UTR) <- colnames(results_UTR, do.NULL = FALSE, prefix = "col")
        colnames(results_UTR)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
        colnames(results_UTR)[(dim(results_UTR)[2]-1):dim(results_UTR)[2]] <- c("ACAT-O","STAAR-O")
      }else
      {
        colnames(results_UTR) <- colnames(results_UTR, do.NULL = FALSE, prefix = "col")
        colnames(results_UTR)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
        colnames(results_UTR)[dim(results_UTR)[2]] <- c("STAAR-B")
      }
    }
    # }else{
    #   results_UTR <- c()
    # }
    
    seqResetFilter(genofile)
    
    #############################################
    #   Promoter-CAGE
    
    ## Promoter
    varid <- seqGetData(genofile, "variant.id")
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    promGobj <- promoters(genes(txdb), upstream = 3000, downstream = 3000)
    
    # Subsetting Promoters that within +/-3kb of TSS and have CAGE signals
    CAGEAnno <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="CAGE")]))
    CAGEBvt <- CAGEAnno!=""
    CAGEidx <- which(CAGEBvt,useNames=TRUE)
    seqSetFilter(genofile,variant.id=varid[CAGEidx])
    seqSetFilter(genofile,promGobj,intersect=TRUE)
    CAGEpromgene <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
    CAGEGene <- unlist(lapply(strsplit(CAGEpromgene,"\\(|\\,|;|-"),`[[`,1))
    ##obtain variants info
    CAGEvchr <- as.numeric(seqGetData(genofile,"chromosome"))
    CAGEvpos <- as.numeric(seqGetData(genofile,"position"))
    CAGEvref <- as.character(seqGetData(genofile,"$ref"))
    CAGEvalt <- as.character(seqGetData(genofile,"$alt"))
    dfPromCAGEVarGene <- data.frame(CAGEvchr,CAGEvpos,CAGEvref,CAGEvalt,CAGEGene)
    
    ## get SNV id
    filter <- seqGetData(genofile, QC_label)
    if(variant_type=="variant")
    {
      SNVlist <- filter == "PASS"
    }
    
    if(variant_type=="SNV")
    {
      SNVlist <- (filter == "PASS") & isSNV(genofile)
    }
    
    if(variant_type=="Indel")
    {
      SNVlist <- (filter == "PASS") & (!isSNV(genofile))
    }
    if(sum(SNVlist)>0){
      variant.id <- seqGetData(genofile, "variant.id")
      variant.id.SNV <- variant.id[SNVlist]
      
      dfPromCAGEVarGene.SNV <- dfPromCAGEVarGene[SNVlist,]
      dfPromCAGEVarGene.SNV$CAGEvpos <- as.character(dfPromCAGEVarGene.SNV$CAGEvpos)
      dfPromCAGEVarGene.SNV$CAGEvref <- as.character(dfPromCAGEVarGene.SNV$CAGEvref)
      dfPromCAGEVarGene.SNV$CAGEvalt <- as.character(dfPromCAGEVarGene.SNV$CAGEvalt)
      
      seqResetFilter(genofile)
      
      rm(dfPromCAGEVarGene)
      gc()
      
      ### Gene
      is.in <- which(dfPromCAGEVarGene.SNV[,5]==gene_name)
      variant.is.in <- variant.id.SNV[is.in]
      
      seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)
      
      ## genotype id
      id.genotype <- seqGetData(genofile,"sample.id")
      # id.genotype.match <- rep(0,length(id.genotype))
      
      id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
      phenotype.id.merge <- data.frame(phenotype.id)
      phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
      id.genotype.match <- phenotype.id.merge$index
      
      ## Genotype
      Geno <- NULL
      if(length(seqGetData(genofile, "variant.id"))<rv_num_cutoff_max_prefilter)
      {
        Geno <- seqGetData(genofile, "$dosage")
        Geno <- Geno[id.genotype.match,,drop=FALSE]
      }
      
      ## impute missing
      if(!is.null(dim(Geno)))
      {
        if(dim(Geno)[2]>0)
        {
          if(geno_missing_imputation=="mean")
          {
            tmp <- matrix_flip_mean(Geno)
            Geno <- tmp$Geno
            MAF <- tmp$MAF
            rm(tmp)
          }
          if(geno_missing_imputation=="minor")
          {
            Geno <- matrix_flip_minor(Geno)$Geno
          }
        }
      }
      
      ## Annotation
      Anno.Int.PHRED.sub <- NULL
      Anno.Int.PHRED.sub.name <- NULL
      
      if(variant_type=="SNV")
      {
        if(Use_annotation_weights)
        {
          for(k in 1:length(Annotation_name))
          {
            if(Annotation_name[k]%in%Annotation_name_catalog$name)
            {
              Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,Annotation_name[k])
              Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name==Annotation_name[k])]))
              
              if(Annotation_name[k]=="CADD")
              {
                Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
              }
              
              if(Annotation_name[k]=="aPC.LocalDiversity")
              {
                Annotation.PHRED.2 <- -10*log10(1-10^(-Annotation.PHRED/10))
                Annotation.PHRED <- cbind(Annotation.PHRED,Annotation.PHRED.2)
                Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,paste0(Annotation_name[k],"(-)"))
              }
              Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub,Annotation.PHRED)
            }
          }
          
          Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
          colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
        }
      }
      
      pvalues <- 1
      if(n_pheno == 1)
      {
        if(0)
        {
          try(pvalues <- STAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff),silent=silent)
        }else
        {
          try(pvalues <- STAAR_Binary_SPA_survival(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff, varRatio=varRatio, MAF=MAF),silent=silent)
        }
      }else
      {
        try(pvalues <- MultiSTAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff),silent=silent)
      }
      
      results_promoter_CAGE <- c()
      if(inherits(pvalues, "list"))
      {
        results_temp <- c()
        results_temp[3] <- "promoter_CAGE"
        results_temp[2] <- chr
        results_temp[1] <- as.character(gene_name)
        results_temp[4] <- pvalues$num_variant
        
        if(!use_SPA0)
        {
          results_temp <- c(results_temp,pvalues$cMAC,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
                            pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
                            pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)
        }else
        {
          results_temp <- c(results_temp,pvalues$cMAC,
                            pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_B)
        }
        
        results_promoter_CAGE <- rbind(results_promoter_CAGE,results_temp)
      }
      
      if(!is.null(results_promoter_CAGE))
      {
        if(!use_SPA0)
        {
          colnames(results_promoter_CAGE) <- colnames(results_promoter_CAGE, do.NULL = FALSE, prefix = "col")
          colnames(results_promoter_CAGE)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
          colnames(results_promoter_CAGE)[(dim(results_promoter_CAGE)[2]-1):dim(results_promoter_CAGE)[2]] <- c("ACAT-O","STAAR-O")
        }else
        {
          colnames(results_promoter_CAGE) <- colnames(results_promoter_CAGE, do.NULL = FALSE, prefix = "col")
          colnames(results_promoter_CAGE)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
          colnames(results_promoter_CAGE)[dim(results_promoter_CAGE)[2]] <- c("STAAR-B")
        }
      }
    }else{
      results_promoter_CAGE <- c()
    }
    
    seqResetFilter(genofile)
    
    ##################################################
    #       Promoter-DHS
    
    # Subsetting Promoters that within +/-3kb of TSS and have rOCRs signals
    
    rOCRsAnno <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="DHS")]))
    rOCRsBvt <- rOCRsAnno!=""
    rOCRsidx <- which(rOCRsBvt,useNames=TRUE)
    seqSetFilter(genofile,variant.id=varid[rOCRsidx])
    
    seqSetFilter(genofile,promGobj,intersect=TRUE)
    rOCRspromgene <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
    rOCRsGene <- unlist(lapply(strsplit(rOCRspromgene,"\\(|\\,|;|-"),`[[`,1))
    ## obtain variants info
    rOCRsvchr <- as.numeric(seqGetData(genofile,"chromosome"))
    rOCRsvpos <- as.numeric(seqGetData(genofile,"position"))
    rOCRsvref <- as.character(seqGetData(genofile,"$ref"))
    rOCRsvalt <- as.character(seqGetData(genofile,"$alt"))
    dfPromrOCRsVarGene <- data.frame(rOCRsvchr,rOCRsvpos,rOCRsvref,rOCRsvalt,rOCRsGene)
    
    ## get SNV id
    filter <- seqGetData(genofile, QC_label)
    if(variant_type=="variant")
    {
      SNVlist <- filter == "PASS"
    }
    
    if(variant_type=="SNV")
    {
      SNVlist <- (filter == "PASS") & isSNV(genofile)
    }
    
    if(variant_type=="Indel")
    {
      SNVlist <- (filter == "PASS") & (!isSNV(genofile))
    }
    
    variant.id <- seqGetData(genofile, "variant.id")
    variant.id.SNV <- variant.id[SNVlist]
    
    dfPromrOCRsVarGene.SNV <- dfPromrOCRsVarGene[SNVlist,]
    dfPromrOCRsVarGene.SNV$rOCRsvpos <- as.character(dfPromrOCRsVarGene.SNV$rOCRsvpos)
    dfPromrOCRsVarGene.SNV$rOCRsvref <- as.character(dfPromrOCRsVarGene.SNV$rOCRsvref)
    dfPromrOCRsVarGene.SNV$rOCRsvalt <- as.character(dfPromrOCRsVarGene.SNV$rOCRsvalt)
    
    seqResetFilter(genofile)
    
    rm(dfPromrOCRsVarGene)
    gc()
    if(sum(SNVlist)>0){
      ### Gene
      is.in <- which(dfPromrOCRsVarGene.SNV[,5]==gene_name)
      variant.is.in <- variant.id.SNV[is.in]
      
      seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)
      
      ## genotype id
      id.genotype <- seqGetData(genofile,"sample.id")
      # id.genotype.match <- rep(0,length(id.genotype))
      
      id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
      phenotype.id.merge <- data.frame(phenotype.id)
      phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
      id.genotype.match <- phenotype.id.merge$index
      
      ## Genotype
      Geno <- NULL
      if(length(seqGetData(genofile, "variant.id"))<rv_num_cutoff_max_prefilter)
      {
        Geno <- seqGetData(genofile, "$dosage")
        Geno <- Geno[id.genotype.match,,drop=FALSE]
      }
      
      ## impute missing
      if(!is.null(dim(Geno)))
      {
        if(dim(Geno)[2]>0)
        {
          if(geno_missing_imputation=="mean")
          {
            tmp <- matrix_flip_mean(Geno)
            Geno <- tmp$Geno
            MAF <- tmp$MAF
            rm(tmp)
          }
          if(geno_missing_imputation=="minor")
          {
            Geno <- matrix_flip_minor(Geno)$Geno
          }
        }
      }
      
      ## Annotation
      Anno.Int.PHRED.sub <- NULL
      Anno.Int.PHRED.sub.name <- NULL
      
      if(variant_type=="SNV")
      {
        if(Use_annotation_weights)
        {
          for(k in 1:length(Annotation_name))
          {
            if(Annotation_name[k]%in%Annotation_name_catalog$name)
            {
              Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,Annotation_name[k])
              Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name==Annotation_name[k])]))
              
              if(Annotation_name[k]=="CADD")
              {
                Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
              }
              
              if(Annotation_name[k]=="aPC.LocalDiversity")
              {
                Annotation.PHRED.2 <- -10*log10(1-10^(-Annotation.PHRED/10))
                Annotation.PHRED <- cbind(Annotation.PHRED,Annotation.PHRED.2)
                Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,paste0(Annotation_name[k],"(-)"))
              }
              Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub,Annotation.PHRED)
            }
          }
          
          Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
          colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
        }
      }
      
      pvalues <- 1
      if(n_pheno == 1)
      {
        if(F)
        {
          try(pvalues <- STAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff),silent=silent)
        }else
        {
          try(pvalues <- STAAR_Binary_SPA_survival(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff, varRatio=varRatio, MAF=MAF),silent=silent)
        }
      }else
      {
        try(pvalues <- MultiSTAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff),silent=silent)
      }
      
      results_promoter_DHS <- c()
      if(inherits(pvalues, "list"))
      {
        results_temp <- c()
        results_temp[3] <- "promoter_DHS"
        results_temp[2] <- chr
        results_temp[1] <- as.character(gene_name)
        results_temp[4] <- pvalues$num_variant
        
        if(!use_SPA0)
        {
          results_temp <- c(results_temp,pvalues$cMAC,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
                            pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
                            pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)
        }else
        {
          results_temp <- c(results_temp,pvalues$cMAC,
                            pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_B)
        }
        
        results_promoter_DHS <- rbind(results_promoter_DHS ,results_temp)
      }
      
      if(!is.null(results_promoter_DHS))
      {
        if(!use_SPA0)
        {
          colnames(results_promoter_DHS) <- colnames(results_promoter_DHS, do.NULL = FALSE, prefix = "col")
          colnames(results_promoter_DHS)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
          colnames(results_promoter_DHS)[(dim(results_promoter_DHS)[2]-1):dim(results_promoter_DHS)[2]] <- c("ACAT-O","STAAR-O")
        }else
        {
          colnames(results_promoter_DHS) <- colnames(results_promoter_DHS, do.NULL = FALSE, prefix = "col")
          colnames(results_promoter_DHS)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
          colnames(results_promoter_DHS)[dim(results_promoter_DHS)[2]] <- c("STAAR-B")
        }
      }
    }else{
      results_promoter_DHS <- c()
    }
    
    seqResetFilter(genofile)
    
    ###########################################
    #        Enhancer-CAGE
    
    #Now extract the GeneHancer with CAGE Signal Overlay
    genehancerAnno <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GeneHancer")]))
    genehancer <- genehancerAnno!=""
    
    CAGEAnno <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="CAGE")]))
    CAGE <- CAGEAnno!=""
    CAGEGeneHancervt <- CAGEAnno!=""&genehancerAnno!=""
    CAGEGeneHanceridx <- which(CAGEGeneHancervt,useNames=TRUE)
    seqSetFilter(genofile,variant.id=varid[CAGEGeneHanceridx])
    
    # variants that covered by whole GeneHancer without CAGE overlap.
    genehancerSet <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GeneHancer")]))
    enhancerGene <- unlist(lapply(strsplit(genehancerSet,"="),`[[`,4))
    enhancer2GENE <- unlist(lapply(strsplit(enhancerGene,";"),`[[`,1))
    enhancervchr <- as.numeric(seqGetData(genofile,"chromosome"))
    enhancervpos <- as.numeric(seqGetData(genofile,"position"))
    enhancervref <- as.character(seqGetData(genofile,"$ref"))
    enhancervalt <- as.character(seqGetData(genofile,"$alt"))
    dfHancerCAGEVarGene <- data.frame(enhancervchr,enhancervpos,enhancervref,enhancervalt,enhancer2GENE)
    
    ## get SNV id
    filter <- seqGetData(genofile, QC_label)
    if(variant_type=="variant")
    {
      SNVlist <- filter == "PASS"
    }
    
    if(variant_type=="SNV")
    {
      SNVlist <- (filter == "PASS") & isSNV(genofile)
    }
    
    if(variant_type=="Indel")
    {
      SNVlist <- (filter == "PASS") & (!isSNV(genofile))
    }
    
    variant.id <- seqGetData(genofile, "variant.id")
    variant.id.SNV <- variant.id[SNVlist]
    
    dfHancerCAGEVarGene.SNV <- dfHancerCAGEVarGene[SNVlist,]
    dfHancerCAGEVarGene.SNV$enhancervpos <- as.character(dfHancerCAGEVarGene.SNV$enhancervpos)
    dfHancerCAGEVarGene.SNV$enhancervref <- as.character(dfHancerCAGEVarGene.SNV$enhancervref)
    dfHancerCAGEVarGene.SNV$enhancervalt <- as.character(dfHancerCAGEVarGene.SNV$enhancervalt)
    
    seqResetFilter(genofile)
    
    rm(dfHancerCAGEVarGene)
    gc()
    
    ### Gene
    is.in <- which(dfHancerCAGEVarGene.SNV[,5]==gene_name)
    variant.is.in <- variant.id.SNV[is.in]
    
    seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)
    
    ## genotype id
    id.genotype <- seqGetData(genofile,"sample.id")
    # id.genotype.match <- rep(0,length(id.genotype))
    
    id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
    phenotype.id.merge <- data.frame(phenotype.id)
    phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
    id.genotype.match <- phenotype.id.merge$index
    
    ## Genotype
    Geno <- NULL
    if(length(seqGetData(genofile, "variant.id"))<rv_num_cutoff_max_prefilter)
    {
      Geno <- seqGetData(genofile, "$dosage")
      Geno <- Geno[id.genotype.match,,drop=FALSE]
    }
    
    ## impute missing
    if(!is.null(dim(Geno)))
    {
      if(dim(Geno)[2]>0)
      {
        if(geno_missing_imputation=="mean")
        {
          tmp <- matrix_flip_mean(Geno)
          Geno <- tmp$Geno
          MAF <- tmp$MAF
          rm(tmp)
        }
        if(geno_missing_imputation=="minor")
        {
          Geno <- matrix_flip_minor(Geno)$Geno
        }
      }
    }
    
    ## Annotation
    Anno.Int.PHRED.sub <- NULL
    Anno.Int.PHRED.sub.name <- NULL
    
    if(variant_type=="SNV")
    {
      if(Use_annotation_weights)
      {
        for(k in 1:length(Annotation_name))
        {
          if(Annotation_name[k]%in%Annotation_name_catalog$name)
          {
            Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,Annotation_name[k])
            Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name==Annotation_name[k])]))
            
            if(Annotation_name[k]=="CADD")
            {
              Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
            }
            
            if(Annotation_name[k]=="aPC.LocalDiversity")
            {
              Annotation.PHRED.2 <- -10*log10(1-10^(-Annotation.PHRED/10))
              Annotation.PHRED <- cbind(Annotation.PHRED,Annotation.PHRED.2)
              Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,paste0(Annotation_name[k],"(-)"))
            }
            Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub,Annotation.PHRED)
          }
        }
        
        Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
        colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
      }
    }
    
    pvalues <- 1
    if(n_pheno == 1)
    {
      if(F)
      {
        try(pvalues <- STAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff),silent=silent)
      }else
      {
        try(pvalues <- STAAR_Binary_SPA_survival(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff, varRatio=varRatio, MAF=MAF),silent=silent)
      }
    }else
    {
      try(pvalues <- MultiSTAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff),silent=silent)
    }
    
    results_enhancer_CAGE <- c()
    if(inherits(pvalues, "list"))
    {
      results_temp <- c()
      results_temp[3] <- "enhancer_CAGE"
      results_temp[2] <- chr
      results_temp[1] <- as.character(gene_name)
      results_temp[4] <- pvalues$num_variant
      
      if(!use_SPA0)
      {
        results_temp <- c(results_temp,pvalues$cMAC,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
                          pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
                          pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)
      }else
      {
        results_temp <- c(results_temp,pvalues$cMAC,
                          pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_B)
      }
      
      results_enhancer_CAGE <- rbind(results_enhancer_CAGE,results_temp)
    }
    
    if(!is.null(results_enhancer_CAGE))
    {
      if(!use_SPA0)
      {
        colnames(results_enhancer_CAGE) <- colnames(results_enhancer_CAGE, do.NULL = FALSE, prefix = "col")
        colnames(results_enhancer_CAGE)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
        colnames(results_enhancer_CAGE)[(dim(results_enhancer_CAGE)[2]-1):dim(results_enhancer_CAGE)[2]] <- c("ACAT-O","STAAR-O")
      }else
      {
        colnames(results_enhancer_CAGE) <- colnames(results_enhancer_CAGE, do.NULL = FALSE, prefix = "col")
        colnames(results_enhancer_CAGE)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
        colnames(results_enhancer_CAGE)[dim(results_enhancer_CAGE)[2]] <- c("STAAR-B")
      }
    }
    # }else{
    #   results_enhancer_CAGE <- c()
    # }
    
    seqResetFilter(genofile)
    
    ##################################################
    #       Enhancer-DHS
    
    rOCRsAnno <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="DHS")]))
    rOCRs <- rOCRsAnno!=""
    rOCRsGeneHancervt <- rOCRsAnno!=""&genehancerAnno!=""
    rOCRsGeneHanceridx <- which(rOCRsGeneHancervt,useNames=TRUE)
    seqSetFilter(genofile,variant.id=varid[rOCRsGeneHanceridx])
    # variants that covered by whole GeneHancer without rOCRs overlap.
    
    genehancerSet <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GeneHancer")]))
    enhancerGene <- unlist(lapply(strsplit(genehancerSet,"="),`[[`,4))
    enhancer2GENE <- unlist(lapply(strsplit(enhancerGene,";"),`[[`,1))
    enhancervchr <- as.numeric(seqGetData(genofile,"chromosome"))
    enhancervpos <- as.numeric(seqGetData(genofile,"position"))
    enhancervref <- as.character(seqGetData(genofile,"$ref"))
    enhancervalt <- as.character(seqGetData(genofile,"$alt"))
    dfHancerrOCRsVarGene <- data.frame(enhancervchr,enhancervpos,enhancervref,enhancervalt,enhancer2GENE)
    
    rm(varid)
    gc()
    
    ## get SNV id
    filter <- seqGetData(genofile, QC_label)
    if(variant_type=="variant")
    {
      SNVlist <- filter == "PASS"
    }
    
    if(variant_type=="SNV")
    {
      SNVlist <- (filter == "PASS") & isSNV(genofile)
    }
    
    if(variant_type=="Indel")
    {
      SNVlist <- (filter == "PASS") & (!isSNV(genofile))
    }
    
    variant.id <- seqGetData(genofile, "variant.id")
    variant.id.SNV <- variant.id[SNVlist]
    
    dfHancerrOCRsVarGene.SNV <- dfHancerrOCRsVarGene[SNVlist,]
    dfHancerrOCRsVarGene.SNV$enhancervpos <- as.character(dfHancerrOCRsVarGene.SNV$enhancervpos)
    dfHancerrOCRsVarGene.SNV$enhancervref <- as.character(dfHancerrOCRsVarGene.SNV$enhancervref)
    dfHancerrOCRsVarGene.SNV$enhancervalt <- as.character(dfHancerrOCRsVarGene.SNV$enhancervalt)
    
    seqResetFilter(genofile)
    
    rm(dfHancerrOCRsVarGene)
    gc()
    
    ### Gene
    is.in <- which(dfHancerrOCRsVarGene.SNV[,5]==gene_name)
    variant.is.in <- variant.id.SNV[is.in]
    
    seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)
    
    ## genotype id
    id.genotype <- seqGetData(genofile,"sample.id")
    # id.genotype.match <- rep(0,length(id.genotype))
    
    id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
    phenotype.id.merge <- data.frame(phenotype.id)
    phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
    id.genotype.match <- phenotype.id.merge$index
    
    ## Genotype
    Geno <- NULL
    if(length(seqGetData(genofile, "variant.id"))<rv_num_cutoff_max_prefilter)
    {
      Geno <- seqGetData(genofile, "$dosage")
      Geno <- Geno[id.genotype.match,,drop=FALSE]
    }
    ## impute missing
    if(!is.null(dim(Geno)))
    {
      if(dim(Geno)[2]>0)
      {
        if(geno_missing_imputation=="mean")
        {
          tmp <- matrix_flip_mean(Geno)
          Geno <- tmp$Geno
          MAF <- tmp$MAF
          rm(tmp)
        }
        if(geno_missing_imputation=="minor")
        {
          Geno <- matrix_flip_minor(Geno)$Geno
        }
      }
    }
    
    ## Annotation
    Anno.Int.PHRED.sub <- NULL
    Anno.Int.PHRED.sub.name <- NULL
    
    if(variant_type=="SNV")
    {
      if(Use_annotation_weights)
      {
        for(k in 1:length(Annotation_name))
        {
          if(Annotation_name[k]%in%Annotation_name_catalog$name)
          {
            Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,Annotation_name[k])
            Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name==Annotation_name[k])]))
            
            if(Annotation_name[k]=="CADD")
            {
              Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
            }
            
            if(Annotation_name[k]=="aPC.LocalDiversity")
            {
              Annotation.PHRED.2 <- -10*log10(1-10^(-Annotation.PHRED/10))
              Annotation.PHRED <- cbind(Annotation.PHRED,Annotation.PHRED.2)
              Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,paste0(Annotation_name[k],"(-)"))
            }
            Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub,Annotation.PHRED)
          }
        }
        
        Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
        print(dim(Anno.Int.PHRED.sub))
        colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
      }
    }
    
    pvalues <- 1
    if(n_pheno == 1)
    {
      if(F)
      {
        try(pvalues <- STAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff),silent=silent)
      }else
      {
        try(pvalues <- STAAR_Binary_SPA_survival(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff,varRatio=varRatio,MAF=MAF),silent=silent)
      }
    }else
    {
      try(pvalues <- MultiSTAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff, varRatio=varRatio, MAF=MAF),silent=silent)
    }
    
    results_enhancer_DHS <- c()
    if(inherits(pvalues, "list"))
    {
      results_temp <- c()
      results_temp[3] <- "enhancer_DHS"
      results_temp[2] <- chr
      results_temp[1] <- as.character(gene_name)
      results_temp[4] <- pvalues$num_variant
      
      if(!use_SPA0)
      {
        results_temp <- c(results_temp,pvalues$cMAC,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
                          pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
                          pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)
      }else
      {
        results_temp <- c(results_temp,pvalues$cMAC,
                          pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_B)
      }
      
      results_enhancer_DHS <- rbind(results_enhancer_DHS,results_temp)
    }
    
    if(!is.null(results_enhancer_DHS))
    {
      if(!use_SPA0)
      {
        colnames(results_enhancer_DHS) <- colnames(results_enhancer_DHS, do.NULL = FALSE, prefix = "col")
        colnames(results_enhancer_DHS)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
        colnames(results_enhancer_DHS)[(dim(results_enhancer_DHS)[2]-1):dim(results_enhancer_DHS)[2]] <- c("ACAT-O","STAAR-O")
      }else
      {
        colnames(results_enhancer_DHS) <- colnames(results_enhancer_DHS, do.NULL = FALSE, prefix = "col")
        colnames(results_enhancer_DHS)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
        colnames(results_enhancer_DHS)[dim(results_enhancer_DHS)[2]] <- c("STAAR-B")
      }
    }
    # }else{
    #   results_enhancer_DHS <- c()
    # }
    
    seqResetFilter(genofile)
    
    ############################################
    #           results
    
    results_noncoding <- list(upstream=results_upstream,downstream=results_downstream,UTR=results_UTR,
                              promoter_CAGE=results_promoter_CAGE,promoter_DHS=results_promoter_DHS,
                              enhancer_CAGE=results_enhancer_CAGE,enhancer_DHS=results_enhancer_DHS)
    
    return(results_noncoding)
  }
  
  
  
  
  
  STAAR_Binary_SPA_survival <- function(genotype,obj_nullmodel,annotation_phred=NULL,
                                        rare_maf_cutoff=0.01,rv_num_cutoff=2,varRatio=1,
                                        tol=.Machine$double.eps^0.25,max_iter=1000,
                                        SPA_p_filter=FALSE,p_filter_cutoff=0.05, MAF, minpv=1e-45){
    #print('aaaaa')
    if(class(genotype)[1] != "matrix" && !(!is.null(attr(class(genotype), "package")) && attr(class(genotype), "package") == "Matrix")){
      stop("genotype is not a matrix!")
    }
    
    if(dim(genotype)[2] == 1){
      stop(paste0("Number of rare variant in the set is less than 2!"))
    }
    
    annotation_phred <- as.data.frame(annotation_phred)
    if(dim(annotation_phred)[1] != 0 & dim(genotype)[2] != dim(annotation_phred)[1]){
      stop(paste0("Dimensions don't match for genotype and annotation!"))
    }
    
    if(!is.null(attr(class(genotype), "package")) && attr(class(genotype), "package") == "Matrix"){
      genotype <- as.matrix(genotype)
    }
    #genotype <- matrix_flip(genotype)
    #MAF <- genotype$MAF
    RV_label <- as.vector((MAF<rare_maf_cutoff)&(MAF>0))
    # Geno_rare <- genotype$Geno[,RV_label]
    G <- genotype[,RV_label]
    rm(genotype)
    #gc()
    MAC <- round(MAF[RV_label]*2*nrow(G))
    
    # rm(genotype)
    # gc()
    annotation_phred <- annotation_phred[RV_label,,drop=FALSE]
    
    if(sum(RV_label) >= rv_num_cutoff){
      # G <- as(Geno_rare,"dgCMatrix")
      MAF <- MAF[RV_label]
      # rm(Geno_rare)
      # gc()
      
      annotation_rank <- 1 - 10^(-annotation_phred/10)
      
      ## beta(1,25)
      w_1 <- dbeta(MAF,1,25)
      ## beta(1,1)
      w_2 <- dbeta(MAF,1,1)
      if(dim(annotation_phred)[2] == 0){
        ## Burden, SKAT, ACAT-V
        w_B <- w_S <- as.matrix(cbind(w_1,w_2))
        w_A <- as.matrix(cbind(w_1^2/dbeta(MAF,0.5,0.5)^2,w_2^2/dbeta(MAF,0.5,0.5)^2))
        
      }else{
        ## Burden
        w_B_1 <- annotation_rank*w_1
        w_B_1 <- cbind(w_1,w_B_1)
        w_B_2 <- annotation_rank*w_2
        w_B_2 <- cbind(w_2,w_B_2)
        w_B <- cbind(w_B_1,w_B_2)
        w_B <- as.matrix(w_B)
        
        ## SKAT
        w_S_1 <- sqrt(annotation_rank)*w_1
        w_S_1 <- cbind(w_1,w_S_1)
        w_S_2 <- sqrt(annotation_rank)*w_2
        w_S_2 <- cbind(w_2,w_S_2)
        w_S <- cbind(w_S_1,w_S_2)
        w_S <- as.matrix(w_S)
        
        ## ACAT-V
        w_A_1 <- annotation_rank*w_1^2/dbeta(MAF,0.5,0.5)^2
        w_A_1 <- cbind(w_1^2/dbeta(MAF,0.5,0.5)^2,w_A_1)
        w_A_2 <- annotation_rank*w_2^2/dbeta(MAF,0.5,0.5)^2
        w_A_2 <- cbind(w_2^2/dbeta(MAF,0.5,0.5)^2,w_A_2)
        w_A <- cbind(w_A_1,w_A_2)
        w_A <- as.matrix(w_A)
      }
      
      # residuals.phenotype <- obj_nullmodel$residuals
      # muhat <- c(obj_nullmodel$mu)
      #print('aaaaa')
      #t1 <- proc.time()
      pvalues <- STAAR_B_Binary_SPA_survival(G=G,mac=MAC,
                                             #residuals=residuals.phenotype,
                                             varRatio=varRatio,
                                             #muhat=muhat,
                                             weights_B=w_B,weights_S=w_S,weights_A=w_A,
                                             tol=tol,max_iter=max_iter,obj_nullmodel=obj_nullmodel)
      #  t2 <- proc.time()
      #t2-t1
      #print(pvalues)
      pvalues[pvalues==0] <- minpv
      num_variant <- sum(RV_label) #dim(G)[2]
      cMAC <- sum(G)
      num_annotation <- dim(annotation_phred)[2]+1
      if(T){
        #CCT <- ACAT
        results_STAAR_O <- CCT(pvalues[1:(6*num_annotation)])
        results_ACAT_O <- CCT(pvalues[c(1,num_annotation+1,2*num_annotation+1,3*num_annotation+1,4*num_annotation+1,5*num_annotation+1)])
        pvalues_STAAR_S_1_25 <- CCT(pvalues[1:num_annotation])
        pvalues_STAAR_S_1_1 <- CCT(pvalues[(num_annotation+1):(2*num_annotation)])
        pvalues_STAAR_B_1_25 <- CCT(pvalues[(2*num_annotation+1):(3*num_annotation)])
        pvalues_STAAR_B_1_1 <- CCT(pvalues[(3*num_annotation+1):(4*num_annotation)])
        pvalues_STAAR_A_1_25 <- CCT(pvalues[(4*num_annotation+1):(5*num_annotation)])
        pvalues_STAAR_A_1_1 <- CCT(pvalues[(5*num_annotation+1):(6*num_annotation)])
        
        results_STAAR_S_1_25 <- c(pvalues[1:num_annotation],pvalues_STAAR_S_1_25)
        results_STAAR_S_1_25 <- data.frame(t(results_STAAR_S_1_25))
        
        results_STAAR_S_1_1 <- c(pvalues[(num_annotation+1):(2*num_annotation)],pvalues_STAAR_S_1_1)
        results_STAAR_S_1_1 <- data.frame(t(results_STAAR_S_1_1))
        
        results_STAAR_B_1_25 <- c(pvalues[(2*num_annotation+1):(3*num_annotation)],pvalues_STAAR_B_1_25)
        results_STAAR_B_1_25 <- data.frame(t(results_STAAR_B_1_25))
        
        results_STAAR_B_1_1 <- c(pvalues[(3*num_annotation+1):(4*num_annotation)],pvalues_STAAR_B_1_1)
        results_STAAR_B_1_1 <- data.frame(t(results_STAAR_B_1_1))
        
        results_STAAR_A_1_25 <- c(pvalues[(4*num_annotation+1):(5*num_annotation)],pvalues_STAAR_A_1_25)
        results_STAAR_A_1_25 <- data.frame(t(results_STAAR_A_1_25))
        
        results_STAAR_A_1_1 <- c(pvalues[(5*num_annotation+1):(6*num_annotation)],pvalues_STAAR_A_1_1)
        results_STAAR_A_1_1 <- data.frame(t(results_STAAR_A_1_1))
        
        if(dim(annotation_phred)[2] == 0){
          colnames(results_STAAR_S_1_25) <- c("SKAT(1,25)","STAAR-S(1,25)")
          colnames(results_STAAR_S_1_1) <- c("SKAT(1,1)","STAAR-S(1,1)")
          colnames(results_STAAR_B_1_25) <- c("Burden(1,25)","STAAR-B(1,25)")
          colnames(results_STAAR_B_1_1) <- c("Burden(1,1)","STAAR-B(1,1)")
          colnames(results_STAAR_A_1_25) <- c("ACAT-V(1,25)","STAAR-A(1,25)")
          colnames(results_STAAR_A_1_1) <- c("ACAT-V(1,1)","STAAR-A(1,1)")
        }else{
          colnames(results_STAAR_S_1_25) <- c("SKAT(1,25)",
                                              paste0("SKAT(1,25)-",colnames(annotation_phred)),
                                              "STAAR-S(1,25)")
          colnames(results_STAAR_S_1_1) <- c("SKAT(1,1)",
                                             paste0("SKAT(1,1)-",colnames(annotation_phred)),
                                             "STAAR-S(1,1)")
          colnames(results_STAAR_B_1_25) <- c("Burden(1,25)",
                                              paste0("Burden(1,25)-",colnames(annotation_phred)),
                                              "STAAR-B(1,25)")
          colnames(results_STAAR_B_1_1) <- c("Burden(1,1)",
                                             paste0("Burden(1,1)-",colnames(annotation_phred)),
                                             "STAAR-B(1,1)")
          colnames(results_STAAR_A_1_25) <- c("ACAT-V(1,25)",
                                              paste0("ACAT-V(1,25)-",colnames(annotation_phred)),
                                              "STAAR-A(1,25)")
          colnames(results_STAAR_A_1_1) <- c("ACAT-V(1,1)",
                                             paste0("ACAT-V(1,1)-",colnames(annotation_phred)),
                                             "STAAR-A(1,1)")
        }
        
        return(list(num_variant = num_variant,
                    cMAC = cMAC,
                    RV_label = RV_label,
                    results_STAAR_O = results_STAAR_O,
                    results_ACAT_O = results_ACAT_O,
                    results_STAAR_S_1_25 = results_STAAR_S_1_25,
                    results_STAAR_S_1_1 = results_STAAR_S_1_1,
                    results_STAAR_B_1_25 = results_STAAR_B_1_25,
                    results_STAAR_B_1_1 = results_STAAR_B_1_1,
                    results_STAAR_A_1_25 = results_STAAR_A_1_25,
                    results_STAAR_A_1_1 = results_STAAR_A_1_1))
      }
    }else{
      stop(paste0("Number of rare variant in the set is less than ",rv_num_cutoff,"!"))
    }
    
  }
  
  #library(ACAT)
  # mac=MAC
  # weights_B=w_B
  # weights_A=w_A
  # weights_S=w_S
  # mac_thres=10
  # 
  STAAR_B_Binary_SPA_survival <- function(G,mac,
                                          varRatio,
                                          weights_B,weights_S,weights_A,
                                          tol,max_iter, mac_thres=20, obj_nullmodel ){
    #G <- matrix_flip(G)$Geno
    
    nsnp00 <- ncol(G)
    muhat <- c(obj_nullmodel$mu)
    ind <- which(obj_nullmodel$y==1)
    ## mu0
    if(T){
      mu0 <- colSums(G * muhat)
      #gc()
      score_skat0_unweight <- colSums(G[ind, ])-mu0
      ii <- which(mu0==0)
      if(length(ii)>0){
        #ii <- tmp
        G <- G[,-ii]
        mu0 <- mu0[-ii]
        weights_B <- weights_B[-ii, ]
        weights_A <- weights_A[-ii, ]
        weights_S <- weights_S[-ii, ]
        score_skat0_unweight <- score_skat0_unweight[-ii]
        mac <- mac[-ii]
      }
    }
    #gc()
    #R_pois <- t(G) %*% (G * (muhat+muhat^2))
    library(Matrix)
    R_pois_direct <- function(G, muhat) {
      # G: a 'dgCMatrix' (or other sparse format), dimension n x p
      # muhat: numeric vector of length n
      # Returns: p x p matrix = t(G) %*% [G * (muhat + muhat^2)]
      #
      # We compute t(G) %*% (W %*% G), where W is diagonal with w_i = muhat_i + muhat_i^2.
      
      # w <- muhat + muhat^2
      W <- Diagonal(x = muhat)         # sparse diagonal matrix
      R <- crossprod(G, W %*% G)   # crossprod(G, X) = t(G) %*% X
      return(R)
    }
    R_pois<- R_pois_direct(G, muhat)
    var0 <- diag(R_pois)
    # library(propagate)
    R00 <- cov2cor(as.matrix(R_pois))
    ### LD
    if(T){
      R00_tmp <- abs(R00)
      R00_tmp[lower.tri(R00_tmp, diag=T)] <- 0
      #tmp <- abs(R00-diag(ncol(R00)))
      rr <- which(R00_tmp>sqrt(.95), arr.ind=T)
      ind_ld <- unique(rr[,2])
      if(length(ind_ld)>0){
        G <- G[,-ind_ld]
        weights_B <- weights_B[-ind_ld, ]
        weights_A <- weights_A[-ind_ld, ]
        weights_S <- weights_S[-ind_ld, ]
        mu0=mu0[-ind_ld]
        R00 <- R00[-ind_ld, -ind_ld]
        var0 <- var0[-ind_ld]
        mac=mac[-ind_ld]
        score_skat0_unweight <- score_skat0_unweight[-ind_ld]
      }
    }
    cs <- sqrt(colSums(abs(R00)))
    gc()
    # ww_mu0 <- (length(mu0)/mu0)/sum(1/mu0)
    # t2 <- proc.time()
    # t2-t1
    if(nsnp00-length(ind_ld)==1){
      wn <- length(weights_B)
      #res <- 1
      res <- rep(scoreTest_SAIGE_survivalTrait_cond_sparseSigma_fast(
        G0=G, obj.noK=obj_nullmodel, varRatio=varRatio)$p.value,5 * wn)
      return(res)
      
    }else{
      residuals <- obj_nullmodel$residuals
      vn <- nrow(G)
      un <- ncol(G)
      wn <- ncol(weights_B)
      
      G_cumu <- G %*% weights_B
      # Bisection initialization
      # xmin <- -100.0
      # xmax <- 100.0
      res <- rep(1,5 * wn)
      # print(res)
      # print(wn)
      # ACAT
      #mac00 <- apply(G[which(c(obj_nullmodel$y)==0),],2,sum)
      # mac11 <- apply(G[ind,],2,sum)
      # mac00 <- mac-mac11
      # mac_min <- pmin(mac00, mac11)
      # print(mac00)
      # print(mac11)
      # print(mac_min)
      # print(length(mac))
      # print(dim(G))
      id_veryrare <- which(mac <= mac_thres)
      id_common <- which(mac > mac_thres)
      # print(id_veryrare)
      # print(id_common)
      n0 <- length(id_veryrare)
      n1 <- length(id_common)
      # print(n0)
      # print(n1)
      pseq <- rep(1, un)
      wseq <- rep(0, un)
      G_sub <- G[,id_common]
      #t1 <- proc.time()
      if(n1>0){
        if(n1==1){
          pseq[1] <- scoreTest_SAIGE_survivalTrait_cond_sparseSigma_fast(
            G0=G_sub, obj.noK=obj_nullmodel, varRatio=varRatio)$p.value
        }else{
          #start_time <- proc.time()
          for (k in 1:n1) {
            #print('ccccc')
            # print(dim(G[,id_common[k]]))
            # print(dim(G_cumu[,1]))
            # print(dim(G))
            #print(k)
            #start_time <- proc.time()
            pseq[k] <- scoreTest_SAIGE_survivalTrait_cond_sparseSigma_fast(
              G0=G_sub[,k], obj.noK=obj_nullmodel, varRatio=varRatio)$p.value
            #end_time <- proc.time()
            #print(end_time-start_time)
            # if(k%%20==0){
            #   gc()
            # }
            #print('dddd')
          }
          #end_time <- proc.time()
          # print(end_time-start_time)
          
        }
      }  
      # t2 <- proc.time()
      # t2-t1
      rm(G_sub)
      if(n0>0){
        #print('ccccc')
        # print(dim(G))
        if(n0==1){
          G_cumu_ultra <- matrix(G[,id_veryrare], nrow=nrow(G), ncol=ncol(weights_B))
        }else{
          G_cumu_ultra <- G[,id_veryrare] %*% weights_B[id_veryrare, ]
        }
        #print('dddd')
      }
      gc()
      #print(res)
      ### calculate LD ratio for SKAT
      #RR <- cor(G)
      if(T){
        #score_skat0 <- weights*apply(g*(y-mu), 2, sum)^2
        #mu0 <- apply(G*c(obj_nullmodel$mu), 2, sum)
        ### gamma  skew & mean
        para <- compute_k_theta_v(mu0)
        # kk <- 4*mu0*(2*mu0+1)^3/(8*mu0^2+22*mu0+1)^2
        # theta1 <- mu0/kk
        # theta2 <-  sqrt((3*mu0^2+mu0)/(kk+kk^2))
        # theta <- (theta1+theta2)/2
        kk <- para[1,]
        theta <- para[2,]
        #R_pois <- cor2cov(R00^2, var=var0)
        
        # R2 <- R00*(2+R00)/3
        # R_pois <- cor2cov(R2, var=var0)
        #denominator <- outer(sqrt(G*kk*theta), sqrt(G*kk*theta))
      }
      
      # print('iter')
      # t1 <- proc.time()
      for (i in 1:wn) {
        
        #print(i)
        # sum0 <- 0.0
        # for (k in 1:n) {
        #   sum0 <- sum0 + x[k] * weights_B[k, i]
        # }
        
        # Calculate p-values
        #respart1 <- Saddle_Binary_SPA_survival(abs(sum0), muhat, G_cumu[, i], tol, max_iter, lower1 = FALSE, varRatio=varRatio)
        if(T){
          
          ### Burden
          # print(G_cumu[,i])
          # t1 <- proc.time()
          respart1 <- scoreTest_SAIGE_survivalTrait_cond_sparseSigma_fast(
            G0=G_cumu[,i]/sum(weights_B[,i]), obj.noK=obj_nullmodel, varRatio=varRatio)$p.value
          res[wn+i] <- respart1 
          
          ### ACAT
          if(n1>0){
            wseq[1:n1] <- weights_A[id_common,i]
          }
          if(n0==0){
            res[2*wn+i] <- CCT(pseq, weights = weights_A[,i]) 
          }else{
            #t1 <- proc.time()
            # print(summary(G_cumu_ultra[,i]))
            pseq[n1 + 1] <-  scoreTest_SAIGE_survivalTrait_cond_sparseSigma_fast(
              G0=G_cumu_ultra[,i]/sum(weights_B[id_veryrare, i]), obj.noK=obj_nullmodel, varRatio=varRatio)$p.value
            # t2 <- proc.time()
            # t2-t1
            wseq[n1 + 1] <- sum(weights_A[id_veryrare, i])/n0 
            res[2 * wn + i] <- CCT(pseq[1:(n1 + 1)], wseq[1:(n1 + 1)])
            
          }
          if(T){
            tmp0 <- pv_skat_ld5(#mu=c(obj_nullmodel$mu), 
              weights=weights_S[,i]^2/sum(weights_S[,i]^2), 
              score_skat0_unweight=score_skat0_unweight,cs=cs,
              #y=c(obj_nullmodel$y), g=G, 
              #R_pois=abs(R00),
              kk=kk, theta=theta, mu0=mu0, varRatio=varRatio, useLD=F, R_pois=abs(R00))
            aa <- try( tmp <- pv_skat_ld5(#mu=c(obj_nullmodel$mu), 
              weights=weights_S[,i]^2/sum(weights_S[,i]^2), 
              score_skat0_unweight=score_skat0_unweight,cs=cs,
              #y=c(obj_nullmodel$y), g=G, 
              #R_pois=abs(R00),
              kk=kk, theta=theta, mu0=mu0, varRatio=varRatio, useLD=T, R_pois=abs(R00)),
              silent = TRUE)
            if(inherits(aa, "try-error")){
              
              tmp <- pv_skat_ld4(#mu=c(obj_nullmodel$mu), 
                weights=weights_S[,i]^2/sum(weights_S[,i]^2), 
                score_skat0_unweight=score_skat0_unweight,cs=cs,
                #y=c(obj_nullmodel$y), g=G, 
                #R_pois=abs(R00),
                kk=kk, theta=theta, mu0=mu0, varRatio=varRatio)
            }
          }
          if(F){
            
            tmp0 <- quiet(
              pv_skat_ld5(
                weights  = weights_S[, i]^2 / sum(weights_S[, i]^2),
                score_skat0_unweight = score_skat0_unweight,
                cs       = cs,
                kk       = kk,
                theta    = theta,
                mu0      = mu0,
                varRatio = varRatio,
                useLD    = FALSE,
                R_pois   = abs(R00)
              )
            )
            
            aa <- try(
              quiet(
                tmp <- pv_skat_ld5(
                  weights  = weights_S[, i]^2 / sum(weights_S[, i]^2),
                  score_skat0_unweight = score_skat0_unweight,
                  cs       = cs,
                  kk       = kk,
                  theta    = theta,
                  mu0      = mu0,
                  varRatio = varRatio,
                  useLD    = TRUE,
                  R_pois   = abs(R00)
                )
              ),
              silent = TRUE
            )
            
            if (inherits(aa, "try-error")) {
              tmp <- quiet(
                pv_skat_ld4(
                  weights  = weights_S[, i]^2 / sum(weights_S[, i]^2),
                  score_skat0_unweight = score_skat0_unweight,
                  cs       = cs,
                  kk       = kk,
                  theta    = theta,
                  mu0      = mu0,
                  varRatio = varRatio
                )
              )
            }
            
          }
          res[i] <- max(tmp, tmp0)
          # print('tmp')
          # print(tmp)
          # print(tmp0)
          # 
          # 
          # res[i] <- pv_skat_ld4(#mu=c(obj_nullmodel$mu), 
          #   weights=weights_S[,i]^2/sum(weights_S[,i]^2), 
          #   score_skat0_unweight=score_skat0_unweight,cs=cs,
          #   #y=c(obj_nullmodel$y), g=G, 
          #   #R_pois=abs(R00),
          #   kk=kk, theta=theta, mu0=mu0, varRatio=varRatio)
          if(is.nan(res[i])){
            res[i] <- 0.999
          }
          #t2 <- proc.time()
          #t3=proc.time()
          #res[i] <- respart1 
          #respart1
          # }else{
          #  res[i] <- res[ wn + i]
          # }
          
          
          #   res[i] <-  res[wn+i]
          
          #  }
          # print('aaa')
          # print(t2-t1)
          # print(t3-t2)
          
        }
        
        # If saddle point approximation fails, set the p-value to 1
        if (res[i] > .9999) {
          res[i] <- .9999
        }
        if (res[wn+i] > .9999) {
          res[wn+i] <- .9999
        }
        if (res[2*wn+i] > .9999) {
          res[2*wn+i] <- .9999
        }
        if (res[3*wn+i] > .9999) {
          res[3*wn+i] <- .9999
        }
        if (res[4*wn+i] > .9999) {
          res[4*wn+i] <- .9999
        }
        # if(i==10){
        #   gc()
        # }
        if(i%%14==0){
          gc()
        }
        #gc()
        # t2 <- proc.time()
        # print(i)
        # print(t2-t1)
      }
      # t2 <- proc.time()
      #  t2-t1
      #gc()
      # print('iter done')
      gc()
      return(res)
      
      
    }
    
    
  }
  
  
  
  Gene_Centric_Noncoding_survival <- function(chr,gene_name,
                                              #category=c("all_categories","plof","plof_ds","missense","disruptive_missense","synonymous","ptv","ptv_ds","all_categories_incl_ptv"),
                                              genofile,obj_nullmodel,rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                              rv_num_cutoff_min_prefilter=2, rv_num_cutoff_max_prefilter=2000,
                                              QC_label="annotation/filter",variant_type=c("SNV","Indel","variant"),geno_missing_imputation=c("mean","minor"),
                                              Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,
                                              Use_annotation_weights=c(TRUE,FALSE),Annotation_name=NULL,
                                              SPA_p_filter=TRUE,p_filter_cutoff=0.05,silent=FALSE,varRatio=1){
    
    ## evaluate choices
    #category <- match.arg(category)
    variant_type <- match.arg(variant_type)
    geno_missing_imputation <- match.arg(geno_missing_imputation)
    
    genes <- genes_info[genes_info[,2]==chr,]
    
    #if(category=="all_categories")
    # {
    results <- noncoding_survival(chr,gene_name,genofile,obj_nullmodel,
                                  #genes,
                                  rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
                                  rv_num_cutoff_min_prefilter=rv_num_cutoff_min_prefilter,
                                  rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
                                  QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                  Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                  Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                                  SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff,silent=silent,varRatio=varRatio)
    # }
    
    return(results)
  }
  
  
  
  
  
  
  
  if(T){
    gamma_weighted_decomp <- function(alpha, beta, w, rho, tol = 1e-3) {
      # ---- basic checks --------------------------------------------------------
      K <- length(alpha)
      if (length(beta) != K || length(w) != K)
        stop("alpha, beta, and w must have the same length (K)")
      if (!all(dim(rho) == c(K, K)))
        stop("rho must be a K×K matrix")
      
      # force diagonal to 1 for safety
      diag(rho) <- 1
      
      # ---- compute alpha_0 -----------------------------------------------------
      idx_lower <- lower.tri(rho)
      pair_sqrt <- sqrt(outer(alpha, alpha))
      alpha0 <- min(rho[idx_lower] * pair_sqrt[idx_lower])
      
      # ---- compute alpha_ij matrix --------------------------------------------
      alpha_ij <- matrix(0, K, K)
      for (i in 1:(K - 1)) {
        for (j in (i + 1):K) {
          alpha_ij[i, j] <- rho[i, j] * sqrt(alpha[i] * alpha[j]) - alpha0
          alpha_ij[j, i] <- alpha_ij[i, j]  # make symmetric for convenience
        }
      }
      
      # ---- compute alpha_z vector ---------------------------------------------
      alpha_z <- alpha - (alpha0 + rowSums(alpha_ij))
      
      # ---- feasibility check ---------------------------------------------------
      if (any(c(alpha0, alpha_ij[idx_lower], alpha_z) < -tol))
        stop("Derived shape parameters contain negative values; correlation matrix is infeasible under the Gamma additive model.")
      
      # Treat tiny negatives as zero
      alpha0          <- if (alpha0 < tol) 0 else alpha0
      alpha_ij[abs(alpha_ij) < tol] <- 0
      alpha_z[abs(alpha_z) < tol]   <- 0
      
      # ---- assemble independent components ------------------------------------
      alpha_list <- numeric(0)
      beta_list  <- numeric(0)
      name_list  <- character(0)
      
      # global component
      if (alpha0 > 0) {
        alpha_list <- c(alpha_list, alpha0)
        beta_list  <- c(beta_list, sum(w * beta))
        name_list  <- c(name_list, "U_global")
      }
      
      # pairwise components
      for (i in 1:(K - 1)) {
        for (j in (i + 1):K) {
          if (alpha_ij[i, j] > 0) {
            alpha_list <- c(alpha_list, alpha_ij[i, j])
            beta_list  <- c(beta_list, w[i] * beta[i] + w[j] * beta[j])
            name_list  <- c(name_list, sprintf("Y_%d_%d", i, j))
          }
        }
      }
      
      # private components
      for (i in 1:K) {
        if (alpha_z[i] > 0) {
          alpha_list <- c(alpha_list, alpha_z[i])
          beta_list  <- c(beta_list, w[i] * beta[i])
          name_list  <- c(name_list, sprintf("Z_%d", i))
        }
      }
      
      names(alpha_list) <- name_list
      names(beta_list)  <- name_list
      
      return(list(alpha_new = alpha_list, beta_new = beta_list))
    }
    get_blocks <- function(R, thresh = 0) {
      stopifnot(is.matrix(R), nrow(R) == ncol(R))
      
      ## 1. 把矩阵转成 0/1 邻接矩阵
      A <- abs(R) > thresh
      diag(A) <- FALSE                      # 去掉自环
      
      ## 2. 识别连通分量
      g        <- igraph::graph_from_adjacency_matrix(A, mode = "undirected")
      comp_lab <- igraph::components(g)$membership
      
      ## 3. 返回列表：每个元素是一组索引
      split(seq_len(nrow(R)), comp_lab)
    }
    
    pv_skat_ld5 <- function(weights,score_skat0_unweight,  cs, kk, theta, mu0, varRatio=1, useLD=T, R_pois, thresh=0.05){
      # mu0 <- apply(g*mu, 2, sum)
      # mu02 <- apply(g*(mu^2), 2, sum)
      # denominator <- outer(sqrt(mu0+mu02), sqrt(mu0+mu02))
      # R_pois <- t(g) %*% (g * (mu+mu^2))/denominator
      #### score_skat0 <- weights*apply(g*(y-mu), 2, sum)^2
      score_skat0 <- weights*score_skat0_unweight^2
      #cs <- sqrt(colSums(R_pois))
      ii <- which(cs>1.01)
      if(length(ii)>1 & useLD){
        blocks <- get_blocks(R_pois[ii,ii], thresh = thresh)
        weights_tmp <- kk_tmp <- theta_tmp <- c()
        for(pp in 1:length(blocks)){
          ll <- blocks[[pp]]
          if(length(ll)>1){
            tmp <-   gamma_weighted_decomp(kk[ii[ll]], theta[ii[ll]], weights[ii[ll]], R_pois[ii[ll],ii[ll]])
            weights_tmp <- c(weights_tmp, rep(1,length(tmp$alpha_new)))
            kk_tmp <- c(kk_tmp, tmp$alpha_new)
            theta_tmp <- c(theta_tmp, tmp$beta_new)
          }else{
            weights_tmp <- c(weights_tmp,weights[ii[ll]])
            kk_tmp <- c(kk_tmp, kk[ii[ll]])
            theta_tmp <- c(theta_tmp, theta[ii[ll]])
          }
        }
        weights_new <- c(weights[-ii], weights_tmp)
        kk_new <- c(kk[-ii], kk_tmp)
        theta_new <- c(theta[-ii], theta_tmp)
      }else{
        theta_new <- theta
        kk_new <- kk
        weights_new <- weights
      }
      
      score_skat0 <- score_skat0/sqrt(varRatio)
      
      #score_skat0 <- weights*apply(g*(y-mu), 2, sum)^2/sqrt(theta)
      #1-pgamma(cut, shape=kk, scale=theta)
      # res <- pv_wgamma(alpha=kk, beta=1, weights=w1, cut=sum(score_skat0/(ldratio)))
      # res$score <- sum(score_skat0)
      # res$ldratio <- ldratio
      #w1 <- weights*theta*sqrt(varRatio)
      #w1 <-  weights*theta
      #library(mgcv)
      pv.unadj <- psum.chisq(sum(score_skat0)/varRatio, df=rep(1,length(theta)),lb=weights*mu0)
      
      #w1 <- weights*sqrt(theta)
      # pv_wgamma(alpha=abs(shape_params_w), beta=1/abs(scale_params_w_new), 
      #           w=sign(shape_params_w*scale_params_w_new ), cut=sum(score_skat0))$pv
      return(pv_wgamma(alpha=kk_new, beta=1, weights=weights_new*theta_new, cut=sum(score_skat0), R_pois=NA, pv.unadj=pv.unadj)$pv)
    }
    
    
    K_wgamma <- function(alpha, beta, weights, tt, R_pois=NA){
      if(!is.na(R_pois[1])){
        return(-sum(alpha*log(1-weights/beta*tt))+weights%*%(R_pois-diag(diag(R_pois)))%*%weights/2*tt^2)
      }else{
        return(-sum(alpha*log(1-weights/beta*tt)))
      }
    }
    K_wgamma_prime <- function(alpha, beta, weights, tt, R_pois=NA){
      if(!is.na(R_pois[1])){
        return(sum(alpha*weights/(beta-weights*tt))+weights%*%(R_pois-diag(diag(R_pois)))%*%weights*tt)
      }else{
        return(sum(alpha*weights/(beta-weights*tt)))
      }
    }
    K_wgamma_prime2 <- function(alpha, beta, weights, tt, R_pois=NA){
      if(!is.na(R_pois[1])){
        return(sum(alpha*weights^2/(beta-weights*tt)^2)+weights%*%(R_pois-diag(diag(R_pois)))%*%weights)
      }else{
        return(sum(alpha*weights^2/(beta-weights*tt)^2))
      }
    }
    K_wgamma_prime_adj <- function(alpha, beta, weights, tt, x, R_pois=NA){
      if(!is.na(R_pois[1])){
        return(sum(alpha*weights/(beta-weights*tt))+weights%*%(R_pois-diag(diag(R_pois)))%*%weights*tt-x)
      }else{  
        return(sum(alpha*weights/(beta-weights*tt))-x)
      }
    }
    
    pv_wgamma <- function(alpha, beta, weights, cut, R_pois=NA, pv.unadj=NULL){
      #mean0 <- sum(alpha*weights/beta)
      # sd0 <- sqrt(sum(alpha*weights^2/beta^2))
      #if(is.null(pv.unadj)){
      #pv.unadj0 <- 1-pnorm(cut, mean=mean0, sd=sd0)
      #}
      #if(is.null(pv.unadj)){
      #  pv.unadj <- pv.unadj0
      # }
      if(pv.unadj<0.001){
        t_hat <- uniroot(K_wgamma_prime_adj, 
                         lower = -100, upper = min(beta/weights) - 1e-6,
                         alpha=alpha, beta=beta, weights=weights, x=cut, R_pois=R_pois)$root
        #print(t_hat)
        w <- sign(t_hat)*(2*(t_hat*cut-K_wgamma(alpha=alpha, beta=beta, weights=weights,tt=t_hat)))^0.5
        v <- t_hat*(K_wgamma_prime2(alpha=alpha, beta=beta, weights=weights,tt=t_hat,R_pois=R_pois))^0.5
        pdf_approx <- 1-pnorm(w+1/w*log(v/w))
        SPA=T
      }else{
        pdf_approx <- pv.unadj
        SPA=F
      }
      if(pdf_approx==0){
        pdf_approx <- pv.unadj
        SPA=F
      }
      return(list(pv=pdf_approx, pv_unadj=pv.unadj, SPA=SPA))
    }
    
    
    # 
    # mu=c(obj_nullmodel$mu)
    # weights=weights_S[,i]/sum(weights_S[,i])
    # y=c(obj_nullmodel$y)
    # g=G
    pv_skat_ld4 <- function(weights,score_skat0_unweight,  cs, kk, theta, mu0, varRatio=1){
      # mu0 <- apply(g*mu, 2, sum)
      # mu02 <- apply(g*(mu^2), 2, sum)
      # denominator <- outer(sqrt(mu0+mu02), sqrt(mu0+mu02))
      # R_pois <- t(g) %*% (g * (mu+mu^2))/denominator
      #### score_skat0 <- weights*apply(g*(y-mu), 2, sum)^2
      score_skat0 <- weights*score_skat0_unweight^2
      #cs <- sqrt(colSums(R_pois))
      score_skat0 <- score_skat0/cs/sqrt(varRatio)
      
      #score_skat0 <- weights*apply(g*(y-mu), 2, sum)^2/sqrt(theta)
      #1-pgamma(cut, shape=kk, scale=theta)
      # res <- pv_wgamma(alpha=kk, beta=1, weights=w1, cut=sum(score_skat0/(ldratio)))
      # res$score <- sum(score_skat0)
      # res$ldratio <- ldratio
      #w1 <- weights*theta*sqrt(varRatio)
      #w1 <-  weights*theta
      #library(mgcv)
      pv.unadj <- psum.chisq(sum(score_skat0)/varRatio, df=rep(1,length(theta)),lb=weights*mu0)
      
      #w1 <- weights*sqrt(theta)
      # pv_wgamma(alpha=abs(shape_params_w), beta=1/abs(scale_params_w_new), 
      #           w=sign(shape_params_w*scale_params_w_new ), cut=sum(score_skat0))$pv
      return(pv_wgamma(alpha=kk, beta=1, weights=weights*theta, cut=sum(score_skat0), R_pois=NA, pv.unadj=pv.unadj)$pv)
    }
    
    
    
  }
  
  if(T){
    compute_k_theta <- function(mu, d1=5, d2=7){
      cut1 <- floor(sqrt(uniroot(ff1_pv,interval=c(0,max(5000,500*mu)), mu=mu, pv=10^(-d1))$root))^2
      cut2 <- floor(sqrt(uniroot(ff1_pv,interval=c(0,max(5000,500*mu)), mu=mu, pv=10^(-d2))$root))^2
      p1 <- ff1(mu, cut1)
      p2 <- ff1(mu, cut2)
      kk <- 4*mu*(2*mu+1)^3/(8*mu^2+22*mu+1)^2
      theta <- mu/kk
      return(optim(c(kk, theta), cut1=cut1, cut2=cut2,p1=p1, p2=p2, fn=ff_gamma_optim)$par)
      
      
    }
    compute_k_theta_v <- Vectorize(compute_k_theta, 'mu')
    ff1 <- function(mu,cut){
      return( 1-ppois(sqrt(cut)+mu, lambda = mu)+ppois(-sqrt(cut)-mu, lambda=mu)
      )
    }
    
    
    ff1_pv <- function(mu,cut, pv){
      return( 1-ppois(sqrt(cut)+mu, lambda = mu)+ppois(-sqrt(cut)-mu, lambda=mu)-pv
      )
    }
    
    ff_gamma_optim <- function(xx, cut1, cut2, p1, p2){
      kk <- xx[1]
      theta <- xx[2]
      return((-log10(1-pgamma(cut1, shape=kk, scale=theta))+log10(p1))^2+
               (-log10(1-pgamma(cut2, shape=kk, scale=theta))+log10(p2))^2)
    }
    
    
  }
  
  
  
}
if(T){
  ### line 36: score change
  
  scoreTest_SAIGE_survivalTrait_cond_sparseSigma_fast=function(G0,  obj.noK, 
                                                               varRatio=1, Cutoff=2, 
                                                               sparseSigma=NULL, isCondition=FALSE, 
                                                               OUT_cond=NULL, G1tilde_P_G2tilde = NULL, 
                                                               G2tilde_P_G2tilde_inv=NULL, aa=NULL, type='origin', 
                                                               score_skat=NULL, pv_acat=NULL,wvar=0, wvar2=0){
    
    #print(summary(G0))
    #G0[is.na(G0)] <- mean(G0, na.rm=T)
    # if(is.na(G0[1])){
    #   print(length(G0))
    #   print(G0[1:5])
    # }
    #t1 <- proc.time()
    N = length(G0)
    AC <- sum(G0)
    AF <- AC/N/2
    MAF <- min(AF, 1-AF)
    mu.a <- obj.noK$mu
    y <- obj.noK$y
    mu2.a=mu.a ##poisson
    #print(table(G0))
    #print('aaaaa')
    if(AF > 0.5){
      G0 = 2-G0
      AC2 = 2*N - AC
    }else{
      AC2 = AC
    }
    Run1=TRUE
    idx_no0<-which(G0>0)
    G0 = matrix(G0, ncol = 1)
    
    #tp0 = proc.time()
    #mug = mean(G0)
    #tp1 = proc.time()
    #print("tp1-tp0")
    #print(tp1-tp0)
    ##########old way to get g_mc
    #G0_mc = matrix(G0-mean(G0), ncol = 1)
    #XVG0_mc = eigenMapMatMult(obj.noK$XV, G0_mc)
    #g_mc = G0_mc - eigenMapMatMult(obj.noK$XXVX_inv, XVG0_mc)
    #########
    #if(!isCondition){
    #  if(IsSparse==TRUE){
    if(type!='skat'){
      if(MAF < 0.05){
        out.score<- Score_Test_Sparse_Survival(obj.noK, G0, mu.a, mu2.a, varRatio);
        #tp2a = proc.time()
        #print("tp2a-tp1")
        #print(tp2a-tp1)
      }else{
        out.score<- Score_Test_Survival(obj.noK, G0, mu.a, mu2.a, varRatio);
        #tp2b = proc.time()
        #print("tp2b-tp1")
        #print(tp2b-tp1)
      }
    }
    #t2 <- proc.time()
    #cat("out.score$Tstat: ", out.score$Tstat, "\n")
    #cat("out.score1$Tstat: ", out.score1$Tstat, "\n")
    #cat("out.score$var2: ", out.score$var2, "\n")
    #cat("out.score1$var2: ", out.score1$var2, "\n")
    ##if(out.score["pval.noadj"] > 0.05){
    if(type=='skat'){
      # XVG0 = (obj.noK$XV_fg%*%G0)
      # g = G0 - (obj.noK$XXVX_inv_fg%*% XVG0)
      g <- G0
      
    }else{
      if(!isCondition){
        if(abs(as.numeric(unlist(out.score["Tstat"])[1])/sqrt(as.numeric(unlist(out.score["var1"])[1]))) < Cutoff){
          if(AF > 0.5){
            out.score$BETA = (-1)*out.score$BETA
            out.score$Tstat = (-1)*out.score$Tstat
          }
          outVec = list(BETA = out.score$BETA, SE = out.score$SE, Tstat = out.score$Tstat, p.value = out.score$pval.noadj,
                        p.value.NA = out.score$pval.noadj, Is.converge = 0, var1 = out.score$var1, var2 = out.score$var2)
          Run1=FALSE
          #print('PASS SPA')
        }else{
          if(MAF < 0.05){
            # XVG0 = eigenMapMatMult(obj.noK$XV_fg, G0)
            # g = G0 - eigenMapMatMult(obj.noK$XXVX_inv_fg, XVG0)
            XVG0 = (obj.noK$XV_fg%*%G0)
            g = G0 - (obj.noK$XXVX_inv_fg%*% XVG0)
            #tp3a = proc.time()
            #print("tp3a-tp2a")
            #print(tp3a-tp2a)
            #print(g[1:20])
            #print(g_mc[1:20])
          }else{
            g = out.score$g_tilde	
            #tp3b = proc.time()
            #print("tp3b-tp2b")
            #print(tp3b-tp2b)
          }
        }
      }else{ #if(!isCondition){
        if(MAF < 0.05){
          # XVG0 = eigenMapMatMult(obj.noK$XV_fg, G0)
          # g = G0 - eigenMapMatMult(obj.noK$XXVX_inv_fg, XVG0)
          XVG0 = (obj.noK$XV_fg%*% G0)
          g = G0 - (obj.noK$XXVX_inv_fg%*% XVG0)
        }else{
          g = out.score$g_tilde
          #tp3b = proc.time()
          #print("tp3b-tp2b")
          #print(tp3b-tp2b)
        }
        
      }
    }
    # t3 <- proc.time()
    # }
    #}else{ #if(!isCondition){
    
    #}
    
    #cat("Run1: ", Run1, "\n")
    if(Run1){
      if(type=='skat'){
        #NAset = which(G0==0)
        #tp4 = proc.time()
        
        # print('qqqqq')
        out1 = scoreTest_SPAGMMAT_survivalTrait_cond_sparseSigma_fast(g, Score = score_skat, pval.noadj = pv_acat, 
                                                                      var1_a = 1, var2_a = 1, AC2, AC,NAset, 
                                                                      y, mu.a, varRatio, Cutoff, sparseSigma=sparseSigma, isCondition=isCondition, 
                                                                      OUT_cond=OUT_cond, G1tilde_P_G2tilde = G1tilde_P_G2tilde, 
                                                                      G2tilde_P_G2tilde_inv=G2tilde_P_G2tilde_inv,
                                                                      aa=aa, type=type,Score_skat=score_skat,wvar=wvar,wvar2=wvar2)
        # print(out1)
        outVec = list(p.value = out1["p.value"]$p.value, 
                      p.value.NA = out1["p.value.NA"]$p.value.NA,
                      Is.converge=out1["Is.converge"]$Is.converge)
        
        return(outVec)
        
      }else{
        #G0 = matrix(G0, ncol = 1)
        #XVG0 = eigenMapMatMult(obj.noK$XV, G0)
        #G = G0  -  eigenMapMatMult(obj.noK$XXVX_inv, XVG0) # G is X adjusted
        #g = G
        NAset = which(G0==0)
        #tp4 = proc.time()
        out1 = scoreTest_SPAGMMAT_survivalTrait_cond_sparseSigma_fast(g, Score = out.score$Tstat, pval.noadj = out.score$pval.noadj, 
                                                                      var1_a = out.score$var1, var2_a = out.score$var2, AC2, AC,NAset, 
                                                                      y, mu.a, varRatio, Cutoff, sparseSigma=sparseSigma, isCondition=isCondition, 
                                                                      OUT_cond=OUT_cond, G1tilde_P_G2tilde = G1tilde_P_G2tilde, 
                                                                      G2tilde_P_G2tilde_inv=G2tilde_P_G2tilde_inv, aa=aa, type=type)
        # print(out1)
        
        if(AF > 0.5){
          out1$BETA = (-1)*out1$BETA
          out1$Tstat = (-1)*out1$Tstat
          if(isCondition){
            out1$BETA_c = (-1) * out1$BETA_c
            out1$Tstat_c = (-1) * out1$Tstat_c
          }
        }
        #tp5 = proc.time()
        #print("tp5-tp4")
        #print(tp5-tp4)	
        if(isCondition){
          outVec = list(BETA = out1["BETA"], SE = out1["SE"], Tstat = out1["Tstat"],p.value = out1["p.value"], 
                        p.value.NA = out1["p.value.NA"], Is.converge=out1["Is.converge"], var1 = out1["var1"], 
                        var2 = out1["var2"], Tstat_c = out1["Tstat_c"], p.value.c = out1["p.value.c"], 
                        var1_c = out1["var1_c"], BETA_c = out1["BETA_c"], SE_c = out1["SE_c"])
          
        }else{
          outVec = list(BETA = out1["BETA"]$BETA, SE = out1["SE"]$SE, 
                        Tstat = out1["Tstat"]$Tstat,p.value = out1["p.value"]$p.value, 
                        p.value.NA = out1["p.value.NA"]$p.value.NA,
                        Is.converge=out1["Is.converge"]$Is.converge, 
                        var1 = out1["var1"]$var1, var2 = out1["var2"]$var2)
          #outVec = list(BETA = BETA, SE = SE, Tstat = Tstat,p.value = p.value, var1 = var1, var2 = var2)
        }
      }
      #cat("p.value: ", as.numeric(outVec$p.value), "\n")
      #cat("p.value.NA: ", as.numeric(outVec$p.value.NA), "\n")
      # t4 <- proc.time()
      # print(t2-t1)
      # print(t3-t2)
      # print(t4-t3)
      #return(outVec)
    }#else{
    # outVec = list(p.value = 1,
    #                p.value.NA = 1)
    return(outVec)
    # }
  }
  
  
  # Score = score_skat
  # pval.noadj = 0.5
  # var1_a = 1
  # var2_a = 1
  # AC=AC2
  # AC_true=AC
  # mu <-  mu.a
  # Score_skat=score_skat
  scoreTest_SPAGMMAT_survivalTrait_cond_sparseSigma_fast=function(g, Score, pval.noadj, var1_a, var2_a, AC, AC_true, NAset, y, 
                                                                  mu, varRatio, Cutoff, sparseSigma=NULL, isCondition=FALSE, OUT_cond=NULL, 
                                                                  G1tilde_P_G2tilde = NULL, G2tilde_P_G2tilde_inv=NULL,
                                                                  type='origin', aa=NULL,Score_skat=NULL, wvar=0,wvar2=0){
    
    #g = G/sqrt(AC)
    #q = innerProduct(g, y)
    #print(type)
    #t1 <- proc.time()
    m1 = innerProduct(g, mu)
    #Tstat = q-m1
    #Tstat = Score
    #var2 = innerProduct(mu, g*g)
    #var2c_old = innerProduct(mu, g_mc*g_mc)
    #var2c = var2 - 2*mug*innerProduct(g,Wq) + mug^2*qW1
    
    #cat("var2c_old ", var2c_old, "\n")
    #cat("var2c ", var2c, "\n")
    #var1 = var2c * varRatio
    var1 = var1_a
    var2 = var2_a
    #t2 <- proc.time()
    if(!is.null(sparseSigma)){
      #pcginvSigma<-pcg(sparseSigma, g)
      pcginvSigma<-solve(sparseSigma, g, sparse=T)
      var2b = as.matrix(t(g) %*% pcginvSigma)
      var1 = var2b * varRatio
    }
    
    if(isCondition){
      T2stat = OUT_cond[,2]
      G1tilde_P_G2tilde = matrix(G1tilde_P_G2tilde,nrow=1)
      Tstat_c = Score - G1tilde_P_G2tilde %*% G2tilde_P_G2tilde_inv %*% T2stat
      var1_c = var1 - G1tilde_P_G2tilde %*% G2tilde_P_G2tilde_inv %*% t(G1tilde_P_G2tilde)
    }
    
    AF = AC_true/(2*length(y))
    #if(AF > 0.5){
    #  Tstat = (-1)*Tstat
    #  if(isCondition){
    #    Tstat_c = (-1)*Tstat_c
    #  }
    #}
    qq <- g
    #t3 <- proc.time()
    # t1 <- proc.time()
    if(type=='origin'){
      aa <- -mu*g
      
      #print(head(qq[-NAset]))
      #qtilde = Tstat/sqrt(var1) * sqrt(var2) + m1
      #Score2 = Score/sqrt(var1) * sqrt(var2)
      #Score2 = Score/sqrt(varRatio)
      #   Score2 = Score/sqrt(varRatio)
      # print('NAset')
      # print(length(NAset))
      Score2 = Score/sqrt(varRatio)
      #Score2 = Score
      if(length(NAset)/length(g) < 0.5){
        # print("Saddle_Prob_Poisson")
        out1 = Saddle_Prob_Poisson(Score=Score2, pval.noadj=pval.noadj, mu = mu, qq=qq, aa=aa, Cutoff = Cutoff, alpha=5*10^-8, m1=m1, var1=var2)
      }else{
        #  print("Saddle_Prob_Poisson_fast")
        out1 = Saddle_Prob_Poisson_fast(Score=Score2, pval.noadj=pval.noadj, qq=qq, aa=aa, mu = mu, qqNB=qq[-NAset], aaNB=aa[-NAset], muNA = mu[NAset], muNB = mu[-NAset], Cutoff = Cutoff, alpha = 5*10^-8, m1=m1, var1=var2)
        #print(out1)
        
      }
      #t3 <- proc.time()
      out1$var1 = var1
      out1$var2 = var2
      
      logOR = Score/var1
      SE = abs(logOR/qnorm(out1$p.value/2))
      out1$BETA=logOR
      out1$SE=SE
      out1$Tstat = Score
    }else if(type=='skat'){
      # <- sum(mu)^2-sum(mu^2)
      out1 = Saddle_Prob_Poisson_fast_skat(Score=Score_skat, pval.noadj=pval.noadj, qq=qq, aa=aa, mu = mu, wvar=wvar,wvar2=wvar2)
      
    }
    #print(t3-t2)
    # t4 <- proc.time()
    # t3 <- proc.time()
    if(isCondition){
      if(var1_c <= (.Machine$double.xmin)^2){
        out1 = c(out1, var1_c = var1_c,BETA_c = NA, SE_c = NA, Tstat_c = Tstat_c, p.value.c = 1, p.value.NA.c = 1)
      }else{
        
        #qtilde_c = Tstat_c/sqrt(var1_c) * sqrt(var2) + m1
        pval.noadj_c<-pchisq((Tstat_c)^2/(var1_c), lower.tail = FALSE, df=1)
        if(length(NAset)/length(g) < 0.5){
          ######To improve
          out1_c = Saddle_Prob_Poisson(Score=Tstat_c, pval.noadj=pval.noadj_c, mu = mu, qq=qq, aa=aa, Cutoff = Cutoff, alpha=5*10^-8, m1=m1, var1=var1_c)
        }else{
          ######To improve
          out1_c = Saddle_Prob_Poisson_fast(Score=Tstat_c, pval.noadj=pval.noadj_c, qq=qq, aa=aa, mu = mu,qqNB=qq[-NAset], aaNB=aa[-NAset], muNA = mu[NAset], muNB = mu[-NAset], Cutoff = Cutoff, alpha = 5*10^-8, m1=m1, var1=var1_c)
        }
        logOR_c = Tstat_c/var1_c
        SE_c = abs(logOR_c/qnorm(out1_c$p.value/2))
        out1 = c(out1, var1_c = var1_c,BETA_c = logOR_c, SE_c = SE_c, Tstat_c = Tstat_c, p.value.c = out1_c$p.value, p.value.NA.c = out1_c$p.value.NA)
      }
      
    }
    return(out1)
  }
  
  
  Saddle_Prob_Poisson_fast=function (Score, pval.noadj, mu, qq,aa, qqNB, aaNB, muNA,muNB, Cutoff = 2, alpha = 5*10^-8, m1, var1){
    #m1 <- sum(mu * g)
    #var1 <- sum(mu * g^2)
    p1 = NULL
    p2 = NULL
    
    #NAmu= m1-sum(gNB*muNB)
    if(F){
      NAsigma=var1-sum(muNB*gNB^2)
    }else{
      NAsigma=var1-sum(muNB*qqNB^2)
    }
    #cat("Score is ", Score, "\n")
    #cat("NAsigma is ", NAsigma, "\n")
    #print(mu[1:20])
    #print(g[1:20])
    
    #NAsigma = sum(muNA*gNA^2)
    #Score <- q - m1
    #qinv = -sign(q - m1) * abs(q - m1) + m1
    #pval.noadj <- pchisq((q - m1)^2/var1, lower.tail = FALSE,
    #    df = 1)
    #Is.converge = TRUE
    
    #if (abs(q - m1)/sqrt(var1) < Cutoff) {
    #    pval = pval.noadj
    #}else {
    #print("Saddle_Prob_Poisson_fast >= Cutoff")
    
    #Korg_Poi_result = Korg_Poi(t=0.1, mu, g)
    #Korg_Poi_fast_result = Korg_Poi_fast(t=0.1, mu, g, gNA,gNB,muNA,muNB,NAmu,NAsigma)	
    #cat("Korg_Poi_result: ", Korg_Poi_result, "\n")
    #cat("Korg_Poi_fast_result: ", Korg_Poi_fast_result, "\n")
    
    
    out.uni1 <- getroot_K1_Poi_fast(0, mu = mu, qq=qq, aa=aa,  q = Score, qqNB=qqNB, aaNB=aaNB, muNA=muNA,muNB=muNB,NAsigma=NAsigma)
    # print(out.uni1)
    out.uni2 <- getroot_K1_Poi_fast(0, mu = mu, qq=qq, aa=aa,  q = (-1)*Score,qqNB=qqNB, aaNB=aaNB, muNA=muNA,muNB=muNB,NAsigma=NAsigma)
    #cat("out.uni1 out.uni2: ", out.uni1$root, " ", out.uni2$root, "\n")
    #print(out.uni2)
    # print('uni1')
    # print(out.uni1)
    # print('uni2')
    # print(out.uni2)
    #print('dddd')
    #print(Score)
    if (out.uni1$Is.converge == TRUE && out.uni2$Is.converge == TRUE) {
      #   if(T){
      p1<-tryCatch(Get_Saddle_Prob_Poi_fast(out.uni1$root, mu, qq=qq, aa=aa, q=Score,qqNB=qqNB, aaNB=aaNB,muNA,muNB,NAsigma),error=function(e) {return(pval.noadj/2)})
      #print(p1)
      p2<-tryCatch(Get_Saddle_Prob_Poi_fast(out.uni2$root, mu,qq=qq, aa=aa,  q=(-1)*Score,qqNB=qqNB, aaNB=aaNB, muNA,muNB,NAsigma),error=function(e) {return(pval.noadj/2)})	
      #	cat("p1 p2: ", p1, " ", p2, "\n")
      
      pval = abs(p1) + abs(p2)
      Is.converge = TRUE
      # }else if(out.uni1$Is.converge == TRUE){
      #   p1<-tryCatch(Get_Saddle_Prob_Poi_fast(out.uni1$root, mu, qq=qq, aa=aa, q=Score,qqNB=qqNB, aaNB=aaNB,muNA,muNB,NAsigma),error=function(e) {return(pval.noadj/2)})
      #   #print(p1)
      #   pval = abs(p1) 
      #   Is.converge = TRUE
      # }else if(out.uni2$Is.converge == TRUE){
      #   p2<-tryCatch(Get_Saddle_Prob_Poi_fast(out.uni2$root, mu,qq=qq, aa=aa,  q=(-1)*Score,qqNB=qqNB, aaNB=aaNB, muNA,muNB,NAsigma),error=function(e) {return(pval.noadj/2)})	
      #   pval = abs(p1) + abs(p2)
      #   Is.converge = TRUE
      #   
    }else{
      print("Error_Converge")
      pval <- 1
      Is.converge = FALSE
    }
    #}
    return(list(p.value = pval, p.value.NA = pval.noadj,
                Is.converge = Is.converge, Score = Score))
  }
  
  #(obj.noK, G0, mu.a, mu2.a, varRatio)
  
  Score_Test_Sparse_Survival <-function(obj.null, G, mu, mu2, varRatio){
    # mu=mu.a; mu2= mu2.a; G=G0; obj.null=obj.noK
    #tp2a0 = proc.time()
    idx_no0<-which(G>0)
    g1<-G[idx_no0]
    noCov = FALSE
    if(dim(obj.null$X1_fg)[2] == 1){
      noCov = TRUE 
    }
    
    A1<-obj.null$XVX_inv_XV_fg[idx_no0,,drop=F]
    X1_fg<-obj.null$X1_fg[idx_no0,,drop=F]
    mu21<-mu2[idx_no0]
    mu1<-mu[idx_no0]
    y1<-obj.null$y[idx_no0]
    if(length(idx_no0) > 1){
      #    cat("idx_no0 ", idx_no0, "\n")
      #cat("dim(X1) ", dim(X1), "\n")
      #cat("dim(A1) ", dim(A1), "\n")
      Z = t(A1) %*% g1
      B<-X1_fg %*% Z
      #cat("dim(Z) ", dim(Z), "\n")
      #cat("dim(B) ", dim(B), "\n")
      g_tilde1 = g1 - B
      #print(g_tilde1[1:100])
      var2 = t(Z) %*% obj.null$XVX_fg %*% Z - t(B^2) %*% mu21 + t(g_tilde1^2) %*% mu21
      var1 = var2 * varRatio
      S1 = crossprod(y1-mu1, g_tilde1)
      
      if(!noCov){
        S_a2 = obj.null$S_a - colSums(X1_fg * (y1 - mu1))
      }else{
        S_a2 = obj.null$S_a - crossprod(X1_fg, y1 - mu1)
      }
      
      S2 = -S_a2 %*% Z
    }else{
      Z = A1 * g1
      B<- sum(X1_fg * Z)
      g_tilde1 = g1 - B
      var2 = (Z) %*% obj.null$XVX_fg %*% t(Z) - t(B^2) %*% mu21 + t(g_tilde1^2) %*% mu21
      var1 = var2 * varRatio
      S1 = crossprod(y1-mu1, g_tilde1)
      S_a2 = obj.null$S_a - X1_fg * (y1 - mu1)
      S2 = sum(-S_a2 * Z)
    }
    
    S<- S1+S2
    ##meanG = mean(G)
    #S = S - meanG*resq
    ##cat("dim(B): ", dim(B), "\n")
    ##cat("dim(Z): ", dim(Z), "\n")
    ##cat("dim(XWq): ", dim(XWq), "\n")
    ##cat("dim(Wq): ", dim(Wq), "\n")
    #tp2a1 = proc.time()
    #print("tp2a1-tp2a0")
    #print(tp2a1-tp2a0)  
    
    #var1centered1 = t(Z) %*% XWq - t(B) %*% (Wq[idx_no0,])
    #var1centered = var2 - 2*meanG*var1centered1 + meanG^2*qW1  
    #var1 = var1centered * varRatio
    
    #tp2a2 = proc.time()
    #print("tp2a2-tp2a1")
    #print(tp2a2-tp2a1)
    
    
    pval.noadj<-pchisq((S)^2/(var1), lower.tail = FALSE, df=1)
    ##add on 10-25-2017
    BETA = S/var1
    SE = abs(BETA/qnorm(pval.noadj/2))
    Tstat = S
    #tp2a3 = proc.time()
    #print("tp2a3-tp2a2")
    #print(tp2a3-tp2a2)
    #return(c(BETA, SE, Tstat, pval.noadj, pval.noadj, 1, var1, var2))
    return(list(BETA=BETA, SE=SE, Tstat=Tstat, pval.noadj=pval.noadj, pval.noadj=pval.noadj, is.converge=TRUE, var1=var1, var2=var2, B=B, Z=Z, g_tilde1=g_tilde1))	
  }
  
  
  
  
  
  Score_Test_Survival<-function(obj.null, G, mu, mu2, varRatio){
    #G = G - meanG
    g <- G  -  obj.null$XXVX_inv_fg %*%  (obj.null$XV_fg %*% G)
    q<- crossprod(g, obj.null$y) 
    m1<-crossprod(mu, g)
    var2<- crossprod(mu2, g^2)
    #var1 = var2 * varRatio
    #S = (q-m1) + meanG*resq
    S = q-m1
    #meanG = mean(G)
    #S = S - meanG*resq
    #var1centered1 = t(Z) %*% XWq - t(B) %*% Wq
    #var1centered = var2 - 2*meanG*var1centered1 + meanG^2*qW1
    var1 = var2 * varRatio
    
    pval.noadj<-pchisq((S)^2/var1, lower.tail = FALSE, df=1)
    
    ##add on 10-25-2017
    BETA = S/var1
    SE = abs(BETA/qnorm(pval.noadj/2))
    #Tstat = S^2
    Tstat = S
    
    #return(c(BETA, SE, Tstat, pval.noadj, pval.noadj, NA, var1, var2))
    #return(c(pval.noadj, pval.noadj, TRUE, var1, var2))
    return(list(BETA=BETA, SE=SE, Tstat=Tstat, pval.noadj=pval.noadj, pval.noadj=pval.noadj, is.converge=TRUE, var1=var1, var2=var2,  g_tilde=g))
  }
  
  
  
  
  
  
  innerProduct <- function(x, y) {
    if (length(x) != length(y)) {
      stop("Vectors must be of the same length")
    }
    return(sum(x * y))
  }
  
  
  
  
  
  Saddle_Prob_Poisson=function (Score, pval.noadj, mu,qq, aa, Cutoff = 2, alpha = 5*10^-8, m1, var1){
    #m1 <- sum(mu * g)
    #var1 <- sum(mu * g^2)
    p1 = NULL
    p2 = NULL
    #cat("Score is ", Score, "\n")
    #print(g[1:20])
    
    
    #Score <- q - m1
    #qinv = -sign(q - m1) * abs(q - m1) + m1
    #pval.noadj <- pchisq((q - m1)^2/var1, lower.tail = FALSE,
    #    df = 1)
    Is.converge = TRUE
    
    #if (abs(q - m1)/sqrt(var1) < Cutoff) {
    #    pval = pval.noadj
    #}else {
    #        print("Saddle_Prob_Poisson")
    #t1 <- proc.time()
    out.uni1 <- getroot_K1_Poi(0, mu = mu, aa=aa, qq=qq, q = Score)
    out.uni2 <- getroot_K1_Poi(0, mu = mu, aa=aa, qq=qq, q = (-1)*Score)
    # print('uni1')
    # print(out.uni1)
    # print('uni2')
    # print(out.uni2)
    #t2 <- proc.time()
    if (out.uni1$Is.converge == TRUE && out.uni2$Is.converge == TRUE) {
      p1 <- tryCatch(Get_Saddle_Prob_Poi(out.uni1$root, mu,qq=qq,aa=aa,  q=Score), error=function(e) {return(pval.noadj/2)})	
      p2 <- tryCatch(Get_Saddle_Prob_Poi(out.uni2$root, mu,qq=qq,aa=aa, q = (-1)*Score), error=function(e) {return(pval.noadj/2)})	
      #p1 <- Get_Saddle_Prob_Poi(out.uni1$root, mu, g, q)
      #p2 <- Get_Saddle_Prob_Poi(out.uni2$root, mu, g, qinv)
      #	    cat("p1 p2: ", p1, " ", p2, "\n")	
      
      pval = abs(p1) + abs(p2)
      Is.converge = TRUE
    }
    else {
      print("Error_Converge")
      #pval <- pval.noadj
      pval <- 1
      Is.converge = FALSE
    }
    # t3 <- proc.time()
    # print(t2-t1)
    # print(t3-t2)
    
    #}
    return(list(p.value = pval, p.value.NA = pval.noadj,
                Is.converge = Is.converge, Score = Score))
  }
  
  
  if(T){
    ##saddlepoint approxmation for sum of weighted Poisson distribution
    Korg_Poi<-function(t, mu, qq, aa)
    {
      n.t<-length(t)
      out<-rep(0,n.t)
      
      for(i in 1:n.t){
        t1<-t[i]
        temp<-mu*(exp(qq*t1)  - 1) + aa*t1
        out[i]<-sum(temp)
      }
      return(out)
    }
    
    K1_Poi<-function(t, mu, qq, aa)
    {
      n.t<-length(t)
      out<-rep(0,n.t)
      
      for(i in 1:n.t){
        t1<-t[i]
        temp<-mu * qq * exp(qq*t1) + aa
        out[i]<-sum(temp)
      }
      return(out)
    }
    
    
    K1_adj_Poi<-function(t, mu, qq, aa, q)
    {
      n.t<-length(t)	
      out<-rep(0,n.t)
      
      for(i in 1:n.t){
        t1<-t[i]
        temp<-mu * qq * exp(qq*t1)  + aa
        out[i]<-sum(temp)-q
      }
      return(out)
    }
    
    
    K2_Poi<-function(t, mu, qq,aa)
    {
      n.t<-length(t)
      out<-rep(0,n.t)
      
      for(i in 1:n.t){
        t1<-t[i]
        temp<-mu * qq^2 * exp(qq*t1)
        out[i]<-sum(temp, na.rm=TRUE)
      }
      return(out)
    }
    
    
    getroot_K1_Poi<-function(init,mu,qq,aa,q,m1,tol=.Machine$double.eps^0.25,maxiter=1000)
    {
      t<-init
      K1_eval<-K1_adj_Poi(t,mu,qq,aa,q)
      #cat("K1_eval ", K1_eval, "\n")
      prevJump<- Inf
      rep<-1
      repeat
      {
        t1 <- proc.time()
        K2_eval<-K2_Poi(t,mu,qq,aa)
        tnew<-t-K1_eval/K2_eval
        if(is.na(tnew))
        {
          conv=FALSE
          break
        }
        if(abs(tnew-t)<tol)
        {
          conv<-TRUE
          break
        }
        if(rep==maxiter)
        {
          conv<-FALSE
          break
        }
        # t2 <- proc.time()
        newK1<-K1_adj_Poi(tnew,mu,qq,aa,q)
        if(sign(K1_eval)!=sign(newK1))
        {
          if(abs(tnew-t)>prevJump-tol)
          {
            tnew<-t+sign(newK1-K1_eval)*prevJump/2
            newK1<-K1_adj_Poi(tnew,mu,qq,aa,q)
            prevJump<-prevJump/2
          } else {
            prevJump<-abs(tnew-t)
          }
        }
        #t3 <- proc.time()
        # print('ttt')
        # print(t3-t2)
        # print(t2-t1)
        rep<-rep+1
        t<-tnew
        K1_eval<-newK1
      } 
      # print('repeat')
      # print(rep)
      # print('root')
      # print(t)
      return(list(root=t,n.iter=rep,Is.converge=conv))
      # }
    }
    
    
    Get_Saddle_Prob_Poi<-function(zeta, mu, qq,aa, q) 
    {
      k1<-Korg_Poi(zeta, mu, qq,aa)
      #cat("k1 is ", k1, "\n")
      k2<-K2_Poi(zeta, mu, qq,aa)
      #cat("k2 is ", k2, "\n")
      if(is.finite(k1) && is.finite(k2))
      {
        temp1<-zeta * q - k1
        
        
        w<-sign(zeta) * (2 *temp1)^{1/2}
        v<- zeta * (k2)^{1/2}
        
        Z.test<-w + 1/w * log(v/w)	
        
        if(Z.test > 0){
          pval<-pnorm(Z.test, lower.tail = FALSE)
        } else {
          pval= -pnorm(Z.test, lower.tail = TRUE)
        }	
      } else {
        pval<-0
      }
      
      return(pval)
    }
  }
  
  if(T){
    ##saddlepoint approxmation for sum of weighted Poisson distribution
    Korg_Poi_fast <- function(t, mu, qq, aa, muNA,muNB,NAsigma, qqNB, aaNB)
    {
      n.t<-length(t)
      out<-rep(0,n.t)
      
      for(i in 1:n.t){
        t1<-t[i]
        temp<-muNB*(exp(qqNB*t1)  - 1) + aaNB*t1
        #out[i]<-sum(temp)+NAmu*t1+0.5*NAsigma*t1^2
        out[i]<-sum(temp)+0.5*NAsigma*t1^2 
      }
      return(out)
    }
    
    
    
    K1_adj_Poi_fast <-function(t, mu, qq, aa, q, muNA,muNB,NAsigma, qqNB, aaNB)
    {
      n.t<-length(t)	
      out<-rep(0,n.t)
      
      for(i in 1:n.t){
        t1<-t[i]
        temp<-muNB * qqNB * exp(qqNB*t1) + aaNB
        #temp2<-NAmu+NAsigma*t1
        temp2<-NAsigma*t1
        out[i]<-sum(temp)-q + temp2
      }
      return(out)
    }
    
    
    K2_Poi_fast<-function(t, mu, qq, aa, muNA,muNB,NAsigma, qqNB, aaNB)
    {
      n.t<-length(t)
      out<-rep(0,n.t)
      
      for(i in 1:n.t){
        t1<-t[i]
        temp<-muNB * qqNB^2 * exp(qqNB*t1)
        out[i]<-sum(temp, na.rm=TRUE) + NAsigma
      }
      return(out)
    }
    
    
    
    getroot_K1_Poi_fast<-function(init,mu,qq,aa, q,m1, qqNB, aaNB,muNA,muNB,NAsigma,
                                  tol=.Machine$double.eps^0.25,maxiter=1000)
    {
      t<-init
      # print('bbbb')
      # print(q)
      # print(K1_adj_Poi_fast(t,mu,qq,aa, q, qqNB=qqNB, aaNB=aaNB,NAsigma=NAsigma, muNA=muNA,muNB=muNB))
      # print(K1_adj_Poi_fast(t,mu,qq,-aa, q, qqNB=qqNB, aaNB=-aaNB,NAsigma=NAsigma, muNA=muNA,muNB=muNB))
      # print(K1_adj_Poi_fast(t,mu,-qq,aa, q, qqNB=-qqNB, aaNB=aaNB,NAsigma=NAsigma, muNA=muNA,muNB=muNB))
      # print(K1_adj_Poi_fast(t,mu,-qq,-aa, q, qqNB=-qqNB, aaNB=-aaNB,NAsigma=NAsigma, muNA=muNA,muNB=muNB))
      K1_eval<-K1_adj_Poi_fast(t,mu,qq,aa, q, qqNB=qqNB, aaNB=aaNB,NAsigma=NAsigma, muNA=muNA,muNB=muNB)
      #cat("K1_eval: ", K1_eval, "\n")
      prevJump<- Inf
      rep<-1
      repeat
      {
        K2_eval<-K2_Poi_fast(t,mu,qq,aa, qqNB=qqNB, aaNB=aaNB,NAsigma=NAsigma, muNA=muNA,muNB=muNB)
        tnew<- t-K1_eval/K2_eval
        if(is.na(tnew))
        {
          conv=FALSE
          break
        }
        if(abs(tnew-t)<tol)
        {
          conv<-TRUE
          break
        }
        if(rep==maxiter)
        {
          conv<-FALSE
          break
        }
        
        newK1<-K1_adj_Poi_fast(tnew,mu,qq,aa,q,muNA,muNB,NAsigma, qqNB, aaNB)
        if(sign(K1_eval)!=sign(newK1))
        {
          if(abs(tnew-t)>prevJump-tol)
          {
            tnew<-t+sign(newK1-K1_eval)*prevJump/2
            newK1<-K1_adj_Poi_fast(tnew,mu,qq,aa,q,muNA,muNB,NAsigma, qqNB, aaNB)
            prevJump<-prevJump/2
          } else {
            prevJump<-abs(tnew-t)
          }
        }
        
        rep<-rep+1
        t<-tnew
        K1_eval<-newK1
      } 
      return(list(root=t,n.iter=rep,Is.converge=conv))
      # }
    }
    
    
    Get_Saddle_Prob_Poi_fast<-function(zeta, mu, qq,aa, q,muNA,muNB,NAsigma, qqNB, aaNB) 
    {
      
      
      k1<-Korg_Poi_fast(zeta, mu, qq,aa,muNA,muNB,NAsigma, qqNB, aaNB)
      #cat("k1 is ", k1, "\n")
      k2<-K2_Poi_fast(zeta, mu,qq,aa,muNA,muNB,NAsigma, qqNB, aaNB)
      #cat("k2 is ", k2, "\n")
      if(is.finite(k1) && is.finite(k2))
      {
        temp1<-zeta * q - k1
        
        
        w<-sign(zeta) * (2 *temp1)^{1/2}
        v<- zeta * (k2)^{1/2}
        
        Z.test<-w + 1/w * log(v/w)	
        
        if(Z.test > 0){
          pval<-pnorm(Z.test, lower.tail = FALSE)
        } else {
          pval= -pnorm(Z.test, lower.tail = TRUE)
        }	
      } else {
        pval<-0
      }
      
      return(pval)
    }
  }
  
  
  
  
  
  
  
  
  if(T){
    ##saddlepoint approxmation for sum of weighted Poisson distribution
    Korg_Poi_fast_skat <- function(t, mu, qq, aa, wvar=0, wvar2=0)
    {
      NBset <- which(qq!=0)
      n.t<-length(t)
      out<-rep(0,n.t)
      
      for(i in 1:n.t){
        t1<-t[i]
        temp<-mu[NBset]*(exp(qq[NBset]*t1)  - 1) + aa[NBset]*t1
        #out[i]<-sum(temp)+NAmu*t1+0.5*NAsigma*t1^2
        out[i]<-sum(temp)+0.5*wvar*t1^2+wvar2*(exp(t1)-1)-wvar2*t1
        #+0.5*NAsigma*t1^2 
      }
      return(out)
    }
    
    #K1_Poi<-function(t, mu, g)
    #{
    #  n.t<-length(t)
    #  out<-rep(0,n.t)
    
    #  for(i in 1:n.t){
    #    t1<-t[i]
    #    temp<-mu * g * exp(g*t1) - mu * g
    #    out[i]<-sum(temp)
    #  }
    #  return(out)
    #}
    
    # tnew,mu,qq,aa,q,wvar=wvar, wvar2=wvar2
    K1_adj_Poi_fast_skat <-function(t, mu, qq, aa, q, wvar=0, wvar2=0)
    {
      n.t<-length(t)	
      out<-rep(0,n.t)
      NBset <- which(qq!=0)
      
      for(i in 1:n.t){
        t1<-t[i]
        temp<- mu[NBset] * qq[NBset] * exp(qq[NBset]*t1) + aa[NBset]
        #temp2<-NAmu+NAsigma*t1
        #temp2<-NAsigma*t1
        #out[i]<-sum(temp)-q + temp2
        if(wvar2==0){
          out[i] <- sum(temp)-q+wvar*t1
        }else{
          out[i] <- sum(temp)-q+wvar*t1+wvar2*exp(t1)-wvar2
        }
      }
      return(out)
    }
    
    
    K2_Poi_fast_skat<-function(t, mu, qq, wvar=0, wvar2=0)
    {
      n.t<-length(t)
      out<-rep(0,n.t)
      NBset <- which(qq!=0)
      
      for(i in 1:n.t){
        t1<-t[i]
        temp<-mu[NBset] * qq[NBset]^2 * exp(qq[NBset]*t1)
        out[i]<-sum(temp, na.rm=TRUE) +wvar+wvar2*exp(t1)
      }
      return(out)
    }
    
    
    getroot_K1_Poi_fast_skat<-function(init,mu,qq,aa, q,
                                       tol=.Machine$double.eps^0.25,
                                       maxiter=1000, wvar=0, wvar2=0)
    {
      t<-init
      # print('bbbb')
      # print(q)
      # print(K1_adj_Poi_fast(t,mu,qq,aa, q, qqNB=qqNB, aaNB=aaNB,NAsigma=NAsigma, muNA=muNA,muNB=muNB))
      # print(K1_adj_Poi_fast(t,mu,qq,-aa, q, qqNB=qqNB, aaNB=-aaNB,NAsigma=NAsigma, muNA=muNA,muNB=muNB))
      # print(K1_adj_Poi_fast(t,mu,-qq,aa, q, qqNB=-qqNB, aaNB=aaNB,NAsigma=NAsigma, muNA=muNA,muNB=muNB))
      # print(K1_adj_Poi_fast(t,mu,-qq,-aa, q, qqNB=-qqNB, aaNB=-aaNB,NAsigma=NAsigma, muNA=muNA,muNB=muNB))
      K1_eval<-K1_adj_Poi_fast_skat(t,mu,qq,aa, q,wvar=wvar, wvar2=wvar2)
      #cat("K1_eval: ", K1_eval, "\n")
      prevJump<- Inf
      rep<-1
      repeat
      {
        K2_eval<-K2_Poi_fast_skat(t,mu,qq,wvar=wvar, wvar2=wvar2)
        tnew<- t-K1_eval/K2_eval
        if(is.na(tnew))
        {
          conv=FALSE
          break
        }
        if(abs(tnew)==Inf){
          conv=FALSE
          break
        }
        if(abs(tnew-t)<tol)
        {
          conv<-TRUE
          break
        }
        if(rep==maxiter)
        {
          conv<-FALSE
          break
        }
        # print('tnew')
        # print(tnew)
        newK1<-K1_adj_Poi_fast_skat(tnew,mu,qq,aa,q,wvar=wvar, wvar2=wvar2)
        # print('newK1')
        # print(newK1)
        # print('K1_eval')
        # print(K1_eval)
        if(sign(K1_eval)!=sign(newK1))
        {
          if(abs(tnew-t)>prevJump-tol)
          {
            tnew<-t+sign(newK1-K1_eval)*prevJump/2
            newK1<-K1_adj_Poi_fast_skat(tnew,mu,qq,aa,q,wvar=wvar, wvar2=wvar2)
            prevJump<-prevJump/2
          } else {
            prevJump<-abs(tnew-t)
          }
        }
        
        rep<-rep+1
        t<-tnew
        K1_eval<-newK1
      } 
      return(list(root=t,n.iter=rep,Is.converge=conv))
      # }
    }
    
    
    Get_Saddle_Prob_Poi_fast_skat<-function(zeta, mu, qq,aa, q, wvar=0, wvar2=0) 
    {
      
      
      k1<-Korg_Poi_fast_skat(zeta, mu, qq,aa,wvar=wvar, wvar2=wvar2)
      #cat("k1 is ", k1, "\n")
      k2<-K2_Poi_fast_skat(zeta, mu,qq,wvar=wvar, wvar2=wvar2)
      #cat("k2 is ", k2, "\n")
      if(is.finite(k1) && is.finite(k2))
      {
        temp1<-zeta * q - k1
        
        
        w<-sign(zeta) * (2 *temp1)^{1/2}
        v<- zeta * (k2)^{1/2}
        
        Z.test<- w + 1/w * log(v/w)	
        # print('Z.test')
        # print(Z.test)
        #print(Z.test)
        if(Z.test > 0){
          pval<-pnorm(Z.test, lower.tail = FALSE)
        } else {
          pval= -pnorm(Z.test, lower.tail = TRUE)
        }	
      } else {
        pval<-0
      }
      
      return(pval)
    }
  }
  
  
  
  Saddle_Prob_Poisson_fast_skat=function (Score, pval.noadj, mu, qq,aa,wvar=0, wvar2=0){
    #m1 <- sum(mu * g)
    #var1 <- sum(mu * g^2)
    p1 = NULL
    p2 = NULL
    
    out.uni1 <- getroot_K1_Poi_fast_skat(0, mu = mu, qq=qq, aa=aa,  q = Score,wvar=wvar, wvar2=wvar2)
    # print(out.uni1)
    # out.uni2 <- getroot_K1_Poi_fast(0, mu = mu, qq=qq, aa=aa,  q = (-1)*Score)
    #cat("out.uni1 out.uni2: ", out.uni1$root, " ", out.uni2$root, "\n")
    # print(out.uni2)
    
    # print('dddd')
    # print(Score)
    if(abs(out.uni1$root)>1000){
      out.uni1$Is.converge <- FALSE
    }
    if (out.uni1$Is.converge == TRUE) {
      #if(T){
      p1<-tryCatch(Get_Saddle_Prob_Poi_fast_skat(out.uni1$root, mu=mu, qq=qq, aa=aa, q=Score,wvar=wvar, wvar2=wvar2),error=function(e) {return(pval.noadj)})
      #print(p1)
      # p2<-tryCatch(Get_Saddle_Prob_Poi_fast(out.uni2$root, mu,qq=qq, aa=aa,  q=(-1)*Score),error=function(e) {return(pval.noadj/2)})	
      # cat("p1 p2: ", p1, " ", p2, "\n")
      # 
      pval = abs(p1) 
      Is.converge = TRUE
      # }else if(out.uni1$Is.converge == TRUE){
      #   p1<-tryCatch(Get_Saddle_Prob_Poi_fast(out.uni1$root, mu, qq=qq, aa=aa, q=Score,qqNB=qqNB, aaNB=aaNB,muNA,muNB,NAsigma),error=function(e) {return(pval.noadj/2)})
      #   #print(p1)
      #   pval = abs(p1)
      #   Is.converge = TRUE
      # }else if(out.uni2$Is.converge == TRUE){
      #   p2<-tryCatch(Get_Saddle_Prob_Poi_fast(out.uni2$root, mu,qq=qq, aa=aa,  q=(-1)*Score,qqNB=qqNB, aaNB=aaNB, muNA,muNB,NAsigma),error=function(e) {return(pval.noadj/2)})
      #   pval = abs(p1) + abs(p2)
      #   Is.converge = TRUE
    }else{
      print("Error_Converge")
      #pval <- pval.noadj
      pval <- 1
      Is.converge = FALSE
    }
    #}
    return(list(p.value = pval, p.value.NA = pval.noadj,
                Is.converge = Is.converge, Score = Score))
  }
  
  
  
  
}


## load required packages
install.packages("mgcv",repos = "http://cran.us.r-project.org")
library(mgcv)
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(STAAR)
library(STAARpipeline)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

#install.packages('propagate',repos = "http://cran.us.r-project.org")
#library(propagate)
library(Matrix)
# library(parallel)
# install.packages("pbapply",repos = "http://cran.us.r-project.org")
# library(pbapply)
# install.packages("doParallel",repos = "http://cran.us.r-project.org")
# library(doParallel)
# install.packages("gtools",repos = "http://cran.us.r-project.org")
# library(gtools)
# #install.packages("dplyr",repos = "http://cran.us.r-project.org")
# library(dplyr)
# library(readr)
# install.packages("fastmatch",repos = "http://cran.us.r-project.org")
# library(fastmatch)
`%fin%` <- function(x, table) {
  #stopifnot(require(fastmatch))
  fmatch(x, table, nomatch = 0L) > 0L
}
library(tidyr)
library(tibble)

print("R Package Loading and Installation Done")

load('Annotation_name_catalog.RData')
modglmm$obj.noK$XVX_inv_XV_fg <- modglmm$obj.noK$XXVX_inv_fg*modglmm$obj.noK$V
modglmm$obj.noK$XVX_fg <- t(modglmm$obj.noK$X1_fg)%*%(modglmm$obj.noK$X1_fg*c(modglmm$obj.noK$mu))
modglmm$obj.noK$S_a <- colSums(modglmm$obj.noK$X1_fg * c(modglmm$obj.noK$y - modglmm$obj.noK$mu))
obj_nullmodel <- modglmm$obj.noK;
obj_nullmodel$id_include <- modglmm$sampleID
obj_nullmodel$residuals <- modglmm$residuals
obj_nullmodel$n.pheno <- 1
rm(modglmm)
#gc()
## output path
## input array id from batch file (Harvard FAS RC cluster)

#if(!file.exists(paste0(output_path,output_file_name,arrayid,".Rdata"))){
###########################################################
#           Main Function 
###########################################################
## gene number in job
chr <- arrayid
print(paste("Chromosome:",chr))
#gds.path <- paste0("ukb.200k.wgs.chr",chr,".pass.annotated.gds")
#seqSetFilter(genofile, variant.sel = which(filt=='PASS'))
#ind0 <-  which((filt0=='PASS'))
#seqSetFilter(genofile, variant.sel =ind0[which(miss<0.01)])

genes <- genes_info
genes_info_chr <- genes_info[genes_info[,2]==chr,]
results_coding <- c()

#genes_info_chr_sig <- genes_info_chr[which(genes_info_chr$hgnc_symbol%in%gene_all), ]
genes_info_chr_sig <- genes_info_chr
sub_seq_id <- 1:nrow(genes_info_chr_sig)
print(sub_seq_id)
if(region=='coding'){
  run_id <- c()
  if(nrow(genes_info_chr_sig)>0){
    for(i in sub_seq_id){
      if(!file.exists(paste0('chr',chr,'_Coding_',i,'.Rdata'))){
        run_id <- c(run_id, i)
      }
    }
  }
  
  if(length(run_id)>0){
    ######################################################################################
    genofile <- seqOpen(gds.path)
    #########################################################################################
    
    for(kk in run_id)
    {
      print(kk)
      
      gene_name <- genes_info_chr_sig[kk,1]
      kk0 <- which(genes_info_chr$hgnc_symbol==gene_name)
      
      results <- Gene_Centric_Coding_survival(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                              rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                              QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                              Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                              Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,varRatio=varRatio)
      #end_time <- proc.time()
      #end_time-start_time
      #results_coding <- append(results_coding,results)
      save(results,file=paste0(output_path1,output_file_name,kk,".Rdata"))
      
      # system(paste0("dx mkdir -p ",
      system(paste0("dx upload ",output_path1,output_file_name,kk,".Rdata",
                    " --path '",output_path,"/chr",arrayid,'_',
                    output_file_name,kk,".Rdata","'"))
      
      
    }
    seqClose(genofile)
  }
}
if(region=='noncoding'){
  run_id <- c()
  if(nrow(genes_info_chr_sig)>0){
    for(i in sub_seq_id){
      if(!file.exists(paste0('chr',chr,'_Noncoding_',i,'.Rdata'))){
        run_id <- c(run_id, i)
      }
    }
  }
  
  if(length(run_id)>0){
    ######################################################################################
   
    genofile <- seqOpen(gds.path)
    #########################################################################################
    
    for(kk in run_id)
    {
      print(kk)
      
      gene_name <- genes_info_chr[kk,1]
      # for(kk in 100:188){
      # print(kk)
      # print(length(which(hwe_all<1e-9)))
      # }
      #start_time <- proc.time()
      results <- Gene_Centric_Noncoding_survival(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                                     rare_maf_cutoff=0.01,rv_num_cutoff=2,rv_num_cutoff_max_prefilter=2000,
                                                     QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                                     Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                     Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,varRatio=varRatio)
      #end_time <- proc.time()
      gc()
      #end_time-start_time
      #results_coding <- append(results_coding,results)
      save(results,file=paste0(output_path1,output_file_name,kk,".Rdata"))
      system(paste0("dx mkdir -p ",
                    output_path))
      # system(paste0("dx mkdir -p ",
      system(paste0("dx upload ",output_path1,output_file_name,kk,".Rdata",
                    " --path '",output_path,"/",
                    output_file_name,kk,"_chr",arrayid,".Rdata","'"))
      
    }
    seqClose(genofile)
  }
  
}






