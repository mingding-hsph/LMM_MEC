
loop<-1000

source('./LMM_MEC.macro.r')

sample_test<-data.frame(c(2000, 2000, 2000))
colnames(sample_test)=c("Sample_n")
sample_test$lp=c(0.1, 0, -0.1)


library(mvtnorm)

for (IV_loop in 1:nrow(sample_test)){
  
  MRMV_model<-matrix(NA,nrow=loop, ncol=5*1) 
  MRMV_IVW_model<-matrix(NA,nrow=loop, ncol=5*1) 
  MMR_2stage_model<-matrix(NA,nrow=loop, ncol=5*1) 
  
  colnames(MRMV_model)<-rbind("mean_x1","se_x1", "pvalue_x1", "ll_x1", "ul_x1")
  colnames(MRMV_IVW_model)<-rbind("mean_x1","se_x1", "pvalue_x1", "ll_x1", "ul_x1")
  colnames(MMR_2stage_model)<-rbind("mean_x1","se_x1", "pvalue_x1", "ll_x1", "ul_x1")
  
  ##start simulation
  
  for (simu_loop in 1:loop)
    
  {
    
    total_n=sample_test$Sample_n[IV_loop]    ##number of participants
    exp_num=1          ##number of risk factors
    snp_num=100     ##number of genetic variants
    
    lp<-matrix(nrow=exp_num)  
    lp[1]=sample_test$lp[IV_loop]
    
    
    
    ###1. simulate snps
    
    corr_0 <- matrix(1, nrow = snp_num, ncol = snp_num)
    
    for (i in 1:(snp_num-1))
    {
      for (j in (i+1):snp_num){
        corr_0[i,j]=corr_0[j,i]= runif(1,min=0, max=0)   
      }
    }
    
    corr_1 <- matrix(1, nrow = snp_num, ncol = snp_num)
    
    for (i in 1:(snp_num-1))
    {
      for (j in (i+1):snp_num){
        corr_1[i,j]=corr_1[j,i]= runif(1,min=sqrt(0.0001), max=sqrt(0.1))
      }
    }
    
    corr_3 <- matrix(1, nrow = snp_num, ncol = snp_num)
    
    for (i in 1:(snp_num-1))
    {
      for (j in (i+1):snp_num){
        corr_3[i,j]=corr_3[j,i]= runif(1,min=sqrt(0.1), max=sqrt(0.3))
      }
    }
    
    
    isSymmetric(corr_0)
    isSymmetric(corr_1)
    isSymmetric(corr_3)
    
    mean_snp<-matrix(0, nrow=snp_num, ncol=1) 
    
    ##the two snps are simulated from the same mean p value and same ld structure
    
    snp_con_1<-rmvnorm(total_n, mean =mean_snp, sigma =corr_1, method = "svd") ##sigma is the covariance matrix.
    snp_con_2<-rmvnorm(total_n, mean =mean_snp, sigma =corr_1, method = "svd") ##sigma is the covariance matrix.
    
    snp_cat_1<-matrix(nrow=total_n, ncol=snp_num)
    snp_cat_2<-matrix(nrow=total_n, ncol=snp_num)
    
    cut_snp<-t(matrix(runif(snp_num,min=0, max=0.8)))  ##generate p values from 0.20-0.5
    
    for (i in 1:snp_num){
      snp_cat_1[,i]=ifelse(snp_con_1[,i]>=cut_snp[i], 1, 0)
      snp_cat_2[,i]=ifelse(snp_con_2[,i]>=cut_snp[i], 1, 0)
    }
    
    #summary(snp_cat_1)
    #summary(snp_cat_2)
    
    data<-data.frame(snp_cat_1+snp_cat_2)
    colnames(data)<-paste0('snp', 1:snp_num)
    #summary(data)
    
    ##check correlation between snps 
    
    corr_snps<-round(cor(data, method = "pearson"), digits=3)
    
    ld_snps<-corr_snps^2
    
    #####2. simulate coefficients of genetic variants with x and y (in pleiotropy scenario)
    
    ##lmk_x represents beta of k genetic variants with m exposures and Y by simulating correlated uniform distrubtions
    
    beta_x<-matrix(0.5, nrow=snp_num,ncol=1)
    
    #beta_x<-runif(snp_num,c(0.8,1.0))
    
    lmk_x<-rmvnorm(exp_num, mean =beta_x, sigma =0.03*corr_snps) 
    
    for (i in 1:exp_num){
      for (j in 1:snp_num){
        lmk_x[i,j]<-0.5*lmk_x[i,j]*sample(c(-1,1,1),1)  
      }
    }
    
    #check correlation by simulating 5000 rmvnorm samples, which should be similar to corr_snps
    #round(cor(lmk_x, method = "pearson"),digits=3)
    
    
    
    #####3. simulate risk factors x
    
    for (i in 1:exp_num)
    {
      for (j in 1:snp_num)
      {
        data[,paste0('a_sub_snp', i,j)]<-data[,paste0('snp', j)]*lmk_x[i,j]
      }
    }
    
    for (i in 1:exp_num)
    {
      ##sum up snps*lmk_x for exposure am
      data[,paste0('a',i)]<-rowSums(data[,c(paste0("a_sub_snp", i, 1:snp_num))])+ rnorm(total_n, mean = 0, sd = 1) 
    }
    
    
    ######4. simulate outcome y
    
    ##lp is the causal effects of exposure on outcome y
    
    
    for (i in 1:exp_num)
    {
      data[,paste0('y_sub_a', i)]<-data[,paste0('a', i)]*lp[i]
    }
    
    ####pleiotropy effects on y
    for (j in 1:snp_num)
    {
      ##scenario 1: no direct pleiotropy
      data[,paste0('y_sub_snp', j)]<-0
      
      ##scenario 2: balanced pleiotropy
      # data[,paste0('y_sub_snp', j)]<-lmk_y_bal[j]*data[,paste0('snp', j)]
      
      ##scenario 3: unbalanced pleiotropy
      
      #data[,paste0('y_sub_snp', j)]<-lmk_y_unbal[j]*data[,paste0('snp', j)]
      
    }
    
    
    data$y<-rowSums(matrix(data[,c(paste0("y_sub_a", 1:exp_num))]))+rnorm(total_n, mean = 0, sd = 1)+rowSums(data[,c(paste0("y_sub_snp", 1:snp_num))])
    
    
    ####5. Obtain summary statistics of b_gx and b_gy
    ##b_gx and b_gy are estimated from two datasets simulated from the same lmk_x, lmk_y, and lp
    
    data_x<-data[1:total_n/2,]
    
    data_y<-data[(total_n/2+1):total_n,]
    
    ##obtain summary statistics
    
    new_data_b<-data.frame(matrix(nrow = snp_num,ncol = exp_num+1))
    new_data_se<-data.frame(matrix(nrow = snp_num,ncol = exp_num+1))
    new_data_pvalue<-data.frame(matrix(nrow = snp_num,ncol = exp_num+1))
    new_data_r2<-data.frame(matrix(nrow = snp_num,ncol = exp_num+1))
    new_data_f<-data.frame(matrix(nrow = snp_num,ncol = exp_num+1))
    
    for (i in 1:exp_num)
    {
      colnames(new_data_b)[i]<-paste0('bx', i)
      colnames(new_data_se)[i]<-paste0('bx_se', i)
      colnames(new_data_pvalue)[i]<-paste0('bx_p', i)
      colnames(new_data_r2)[i]<-paste0('bx_r', i)
      colnames(new_data_f)[i]<-paste0('bx_f', i)
      
    }
    
    colnames(new_data_b)[exp_num+1]<-paste0('by')
    colnames(new_data_se)[exp_num+1]<-paste0('by_se')
    colnames(new_data_pvalue)[exp_num+1]<-paste0('by_p')
    colnames(new_data_r2)[exp_num+1]<-paste0('by_r')
    colnames(new_data_f)[exp_num+1]<-paste0('by_f')
    
    ##obtain estimated lmk_x using data_x
    
    for (j in 1:snp_num)
    {
      for (i in 1:exp_num)
      {
        new_data_b[j,i]<-summary(lm(data_x[,c(paste0('a', i))]~data_x[,c(paste0('snp', j))]))$coeff[2,1]
        
        new_data_se[j,i]<-summary(lm(data_x[,c(paste0('a', i))]~data_x[,c(paste0('snp', j))]))$coeff[2,2]
        
        new_data_pvalue[j,i]<-summary(lm(data_x[,c(paste0('a', i))]~data_x[,c(paste0('snp', j))]))$coeff[2,4]
        
        new_data_r2[j,i]<-summary(lm(data_x[,c(paste0('a', i))]~data_x[,c(paste0('snp', j))]))$adj.r.squared
        
        new_data_f[j,i]<-summary(lm(data_x[,c(paste0('a', i))]~data_x[,c(paste0('snp', j))]))$fstatistic[1]
        
      }
    }
    
    ##obtain estimated effect of genetic variants on y using data_y, which is the sum of
    ##1. direct effect: lmk_y_bal/lmk_y_unbal
    ##2. effect mediated through x1,2,x3
    
    for (j in 1:snp_num)
    {
      new_data_b[j,exp_num+1]<-summary(lm(data_y$y~data_y[,c(paste0('snp', j))]))$coeff[2,1]
      
      new_data_se[j,exp_num+1]<-summary(lm(data_y$y~data_y[,c(paste0('snp', j))]))$coeff[2,2]
      
      new_data_pvalue[j,exp_num+1]<-summary(lm(data_y$y~data_y[,c(paste0('snp', j))]))$coeff[2,4]
      
      new_data_r2[j,exp_num+1]<-summary(lm(data_y$y~data_y[,c(paste0('snp', j))]))$adj.r.squared
      
      new_data_f[j,exp_num+1]<-summary(lm(data_y$y~data_y[,c(paste0('snp', j))]))$fstatistic[1]
      
    }
    
    new_data<-cbind(new_data_b,new_data_se,new_data_pvalue,new_data_r2,new_data_f)
  
    
    corr_snps
    
    corr_x<-matrix(1)
    

    
    MMR<-LMM_MEC(
      betaX=cbind(new_data$bx1),
      betaY=new_data$by, 
      betaX_se=cbind(new_data$bx_se1),
      betaY_se=new_data$by_se,
      corr_snps_lmmmec =corr_snps,
      corr_X_lmmmec =corr_x,
      loop_rem=500,
      cutoff_rem=0.00001,
      mec_loop = 100
    )
    
    
    MMR_2stage_model[simu_loop,1]<-t(MMR[,1])
    MMR_2stage_model[simu_loop,2]<-t(MMR[,2])
    MMR_2stage_model[simu_loop,3]<-t(MMR[,3])
    MMR_2stage_model[simu_loop,4]<-t(MMR[,4])
    MMR_2stage_model[simu_loop,5]<-t(MMR[,5])
    
    ##install.packages("MendelianRandomization")
    
    library(MendelianRandomization)
    
    MRMV<-mr_mvegger(mr_mvinput(bx = cbind(new_data$bx1), 
                                bxse = cbind(new_data$bx_se1),
                                by = new_data$by, byse = new_data$by_se, correlation=ld_snps))
    
    MRMV_IVW<-mr_mvivw(mr_mvinput(bx = cbind(new_data$bx1), 
                                  bxse = cbind(new_data$bx_se1),
                                  by = new_data$by, byse = new_data$by_se, correlation=ld_snps))
    
    
    ##output
    
    MRMV_model[simu_loop,1]<-MRMV$Estimate
    MRMV_model[simu_loop,2]<-MRMV$StdError.Est
    MRMV_model[simu_loop,3]<-MRMV$Pvalue.Est
    MRMV_model[simu_loop,4]<-MRMV$CILower.Est
    MRMV_model[simu_loop,5]<-MRMV$CIUpper.Est
    
    MRMV_IVW_model[simu_loop,1]<-MRMV_IVW$Estimate
    MRMV_IVW_model[simu_loop,2]<-MRMV_IVW$StdError
    MRMV_IVW_model[simu_loop,3]<-MRMV_IVW$Pvalue
    MRMV_IVW_model[simu_loop,4]<-MRMV_IVW$CILower
    MRMV_IVW_model[simu_loop,5]<-MRMV_IVW$CIUpper
    
  }     ##end of simu_loop
  
  
  
  ##estimate mean strength of IV
  
  library(tidyverse)
  
  meanF<-new_data%>%
    summarise(
      F1=colMeans(matrix((new_data$bx1)^2/(new_data$bx_se1)^2))
    )
  
  
  mean_F123<-meanF$F1
  
  MMR_2stage_model<-data.frame(MMR_2stage_model)
  
  
  
  MMR_2stage_model$lp1=lp[1]
  MMR_2stage_model$coverage1=NA
  MMR_2stage_model$power1=NA
  
  
  for (i in 1:nrow(MMR_2stage_model))
  {  
    ##Coverage: proportion of 95% CIs including the true causal effect over all simulated datasets; 
    
    if (MMR_2stage_model$lp1[i]>=MMR_2stage_model$ll_x1[i] & MMR_2stage_model$lp1[i]<=MMR_2stage_model$ul_x1[i])
    { MMR_2stage_model$coverage1[i]=1
    } else {
      MMR_2stage_model$coverage1[i]=0
    }
    
    
    
    ##power is estimated by the proportion of 95% CIs excluding zero in settings with significant causal effect; 
    
    if (MMR_2stage_model$pvalue_x1[i]<0.05 & MMR_2stage_model$lp1[i]!=0) 
    { MMR_2stage_model$power1[i]=1
    } else if (MMR_2stage_model$pvalue_x1[i]>=0.05 & MMR_2stage_model$lp1[i]!=0)
    { MMR_2stage_model$power1[i]=0
    }
    
    
  }  
  
  MRMV_model<-data.frame(MRMV_model)
  
  MRMV_model$lp1=lp[1]
  MRMV_model$coverage1=NA
  MRMV_model$power1=NA
  
  for (i in 1:nrow(MRMV_model))
  {  
    
    ##Coverage: proportion of 95% CIs including the true causal effect over all simulated datasets; 
    
    if (MRMV_model$lp1[i]>=MRMV_model$ll_x1[i] & MRMV_model$lp1[i]<=MRMV_model$ul_x1[i])
    { MRMV_model$coverage1[i]=1
    } else {
      MRMV_model$coverage1[i]=0
    }
    
    
    ##power is estimated by the proportion of 95% CIs excluding zero in settings with significant causal effect; 
    
    if (MRMV_model$pvalue_x1[i]<0.05 & MRMV_model$lp1[i]!=0) 
    { MRMV_model$power1[i]=1
    } else if (MRMV_model$pvalue_x1[i]>=0.05 & MRMV_model$lp1[i]!=0)
    { MRMV_model$power1[i]=0
    }
    
    
  }
  
  
  
  ##ivw
  
  MRMV_IVW_model<-data.frame(MRMV_IVW_model)
  
  MRMV_IVW_model$lp1=lp[1]
  MRMV_IVW_model$coverage1=NA
  MRMV_IVW_model$power1=NA
  
  
  for (i in 1:nrow(MRMV_IVW_model))
  {  
    
    ##Coverage: proportion of 95% CIs including the true causal effect over all simulated datasets; 
    
    if (MRMV_IVW_model$lp1[i]>=MRMV_IVW_model$ll_x1[i] & MRMV_IVW_model$lp1[i]<=MRMV_IVW_model$ul_x1[i])
    { MRMV_IVW_model$coverage1[i]=1
    } else {
      MRMV_IVW_model$coverage1[i]=0
    }
    
    
    
    ##power is estimated by the proportion of 95% CIs excluding zero in settings with significant causal effect; 
    
    if (MRMV_IVW_model$pvalue_x1[i]<0.05 & MRMV_IVW_model$lp1[i]!=0) 
    { MRMV_IVW_model$power1[i]=1
    } else if (MRMV_IVW_model$pvalue_x1[i]>=0.05 & MRMV_IVW_model$lp1[i]!=0)
    { MRMV_IVW_model$power1[i]=0
    }
    
    
  }
  
  
  MMR_v1_mean=mean(MMR_2stage_model$mean_x1)
  MMR_v1_std=mean(MMR_2stage_model$se_x1)
  MMR_v1_coverage=mean(MMR_2stage_model$coverage1)
  MMR_v1_power=mean(MMR_2stage_model$power1)
  
  MRMV_v1_mean=mean(MRMV_model$mean_x1)
  MRMV_v1_std=mean(MRMV_model$se_x1)
  MRMV_v1_coverage=mean(MRMV_model$coverage1)
  MRMV_v1_power=mean(MRMV_model$power1)
  
  MRMV_IVW_v1_mean=mean(MRMV_IVW_model$mean_x1)
  MRMV_IVW_v1_std=mean(MRMV_IVW_model$se_x1)
  MRMV_IVW_v1_coverage=mean(MRMV_IVW_model$coverage1)
  MRMV_IVW_v1_power=mean(MRMV_IVW_model$power1)
  
  sample_test$meanF[IV_loop]<-mean_F123
  sample_test$tslmm_mean[IV_loop]<-MMR_v1_mean
  sample_test$tslmm_se[IV_loop]<-MMR_v1_std
  sample_test$tslmm_cover[IV_loop]<-MMR_v1_coverage
  sample_test$tslmm_power[IV_loop]<-MMR_v1_power
  
  sample_test$egger_mean[IV_loop]<-MRMV_v1_mean
  sample_test$egger_se[IV_loop]<-MRMV_v1_std
  sample_test$egger_cover[IV_loop]<-MRMV_v1_coverage
  sample_test$egger_power[IV_loop]<-MRMV_v1_power
  
  sample_test$ivw_mean[IV_loop]<-MRMV_IVW_v1_mean
  sample_test$ivw_se[IV_loop]<-MRMV_IVW_v1_std
  sample_test$ivw_cover[IV_loop]<-MRMV_IVW_v1_coverage
  sample_test$ivw_power[IV_loop]<-MRMV_IVW_v1_power
  
  write.csv(sample_test,"./results/tableS2/ld1.noP.U.sample.csv")
  
}


