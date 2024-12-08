
loop<-1000

source('./LMM_MEC.macro.r')


sample_test<-data.frame(c(2000, 10000, 20000))
colnames(sample_test)=c("Sample_n")

sample_test$meanF<-NA
sample_test$p1_mean<-NA
sample_test$p2_mean<-NA
sample_test$p3_mean<-NA

sample_test$p1_se<-NA
sample_test$p2_se<-NA
sample_test$p3_se<-NA

sample_test$p1_cover<-NA
sample_test$p2_cover<-NA
sample_test$p3_cover<-NA

sample_test$p1_power<-NA
sample_test$p2_power<-NA
sample_test$p3_power<-NA

sample_test$p1_FPR<-NA
sample_test$p2_FPR<-NA
sample_test$p3_FPR<-NA


library(mvtnorm)

for (IV_loop in 1:nrow(sample_test)){
  
  MMR_2stage_model<-matrix(NA,nrow=loop, ncol=5*3+1) 
  
  ##start simulation
  
  for (simu_loop in 1:loop)
    
  {
    
    total_n=sample_test$Sample_n[IV_loop]    ##number of participants
    exp_num=3          ##number of risk factors
    snp_num=100     ##number of genetic variants
    
    
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
    
    ##balanced pleiotropy
    #not violate INSIDE
    ##account for LD
    
    #beta_y_bal<-runif(snp_num,min=-0.001, max=0.001)
    
    ##generate balanced pleiotropy 
    
    beta_y_bal<-matrix(0, nrow=snp_num,ncol=1)
    
    lmk_y_bal<-rmvnorm(1, mean =beta_y_bal, sigma =0.03*corr_snps) 
    
    ## scale the effect of g on y
    
    lmk_y_bal[1,j]<-0.2*lmk_y_bal[1,j]
    
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
    
    lp<-matrix(nrow=exp_num)  
    
    for (i in 1:exp_num)
    {
      
      lp[1]=0.3
      lp[2]=0
      lp[3]=-0.3
      data[,paste0('y_sub_a', i)]<-data[,paste0('a', i)]*lp[i]
    }
    
    ####pleiotropy effects on y
    for (j in 1:snp_num)
    {
      ##scenario 1: no direct pleiotropy
      #data[,paste0('y_sub_snp', j)]<-0
      
      ##scenario 2: balanced pleiotropy
      data[,paste0('y_sub_snp', j)]<-lmk_y_bal[j]*data[,paste0('snp', j)]
      
    }
    
    
    data$y<-rowSums(data[,c(paste0("y_sub_a", 1:exp_num))])+rnorm(total_n, mean = 0, sd = 1)+rowSums(data[,c(paste0("y_sub_snp", 1:snp_num))])
    
    
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
    
    ####6. The new_data will be used for MVMR 
    
    corr_snps
    
    corr_x<-cor(data[,c('a1','a2','a3')])
    

    MMR<-LMM_MEC(
      betaX=cbind(new_data$bx1, new_data$bx2, new_data$bx3),
      betaY=new_data$by, 
      betaX_se=cbind(new_data$bx_se1, new_data$bx_se2, new_data$bx_se3),
      betaY_se=new_data$by_se,
      corr_snps_lmmmec =corr_snps,
      corr_X_lmmmec =corr_x,
      loop_rem=500,
      cutoff_rem=0.00001,
      mec_loop = 100
    )
    
    
    MMR_2stage_model[simu_loop,1:3]<-t(MMR[,1])
    MMR_2stage_model[simu_loop,(1:3)+3]<-t(MMR[,2])
    MMR_2stage_model[simu_loop,(1:3)+2*3]<-t(MMR[,3])
    MMR_2stage_model[simu_loop,(1:3)+3*3]<-t(MMR[,4])
    MMR_2stage_model[simu_loop,(1:3)+4*3]<-t(MMR[,5])
    
    colnames(MMR_2stage_model)<-rbind("mean_x1","mean_x2","mean_x3","se_x1","se_x2","se_x3",
                                      "pvalue_x1", "pvalue_x2","pvalue_x3",
                                      "ll_x1","ll_x2","ll_x3", "ul_x1",  "ul_x2", "ul_x3",
                                      "meanF")
    
  }     ##end of simu_loop
  
  
  
  ##estimate mean strength of IV
  
  library(tidyverse)
  
  meanF<-new_data%>%
    summarise(
      F1=mean((new_data$bx1)^2/(new_data$bx_se1)^2),
      F2=mean((new_data$bx2)^2/(new_data$bx_se2)^2),
      F3=mean((new_data$bx3)^2/(new_data$bx_se3)^2),
      
    )
  
  mean_F123=(meanF$F1+meanF$F2+meanF$F3)/3
  
  scenario=cbind(lp[1],lp[2],lp[3])  
  
  MMR_2stage_model<-data.frame(MMR_2stage_model)
  
  MMR_2stage_model$lp1=scenario[1,1]
  MMR_2stage_model$lp2=scenario[1,2]
  MMR_2stage_model$lp3=scenario[1,3]
  
  MMR_2stage_model$coverage1=NA
  MMR_2stage_model$power1=NA
  MMR_2stage_model$FPR1=NA
  
  MMR_2stage_model$coverage2=NA
  MMR_2stage_model$power2=NA
  MMR_2stage_model$FPR2=NA
  
  MMR_2stage_model$coverage3=NA
  MMR_2stage_model$power3=NA
  MMR_2stage_model$FPR3=NA
  
  for (i in 1:nrow(MMR_2stage_model))
  {  
    ##Coverage: proportion of 95% CIs including the true causal effect over all simulated datasets; 
    
    if (MMR_2stage_model$lp1[i]>=MMR_2stage_model$ll_x1[i] & MMR_2stage_model$lp1[i]<=MMR_2stage_model$ul_x1[i])
    { MMR_2stage_model$coverage1[i]=1
    } else {
      MMR_2stage_model$coverage1[i]=0
    }
    
    if (MMR_2stage_model$lp2[i]>=MMR_2stage_model$ll_x2[i] & MMR_2stage_model$lp2[i]<=MMR_2stage_model$ul_x2[i])
    { MMR_2stage_model$coverage2[i]=1
    } else {
      MMR_2stage_model$coverage2[i]=0
    }
    
    
    if (MMR_2stage_model$lp3[i]>=MMR_2stage_model$ll_x3[i] & MMR_2stage_model$lp3[i]<=MMR_2stage_model$ul_x3[i])
    { MMR_2stage_model$coverage3[i]=1
    } else {
      MMR_2stage_model$coverage3[i]=0
    }
    
    
    ##false-positive rate: the proportion of significant causal estimates in settings with no causal effect. 
    
    if (MMR_2stage_model$pvalue_x1[i]<0.05 & MMR_2stage_model$lp1[i]==0) 
    { MMR_2stage_model$FPR1[i]=1
    } else if (MMR_2stage_model$pvalue_x1[i]>=0.05 & MMR_2stage_model$lp1[i]==0) {
      MMR_2stage_model$FPR1[i]=0  
    }
    
    if (MMR_2stage_model$pvalue_x2[i]<0.05 & MMR_2stage_model$lp2[i]==0) 
    { MMR_2stage_model$FPR2[i]=1
    } else if (MMR_2stage_model$pvalue_x2[i]>=0.05 & MMR_2stage_model$lp2[i]==0) {
      MMR_2stage_model$FPR2[i]=0  
    }
    
    
    if (MMR_2stage_model$pvalue_x3[i]<0.05 & MMR_2stage_model$lp3[i]==0) 
    { MMR_2stage_model$FPR3[i]=1
    } else if (MMR_2stage_model$pvalue_x3[i]>=0.05 & MMR_2stage_model$lp3[i]==0) {
      MMR_2stage_model$FPR3[i]=0  
    }
    
    
    ##power is estimated by the proportion of 95% CIs excluding zero in settings with significant causal effect; 
    
    if (MMR_2stage_model$pvalue_x1[i]<0.05 & MMR_2stage_model$lp1[i]!=0) 
    { MMR_2stage_model$power1[i]=1
    } else if (MMR_2stage_model$pvalue_x1[i]>=0.05 & MMR_2stage_model$lp1[i]!=0)
    { MMR_2stage_model$power1[i]=0
    }
    
    if (MMR_2stage_model$pvalue_x2[i]<0.05 & MMR_2stage_model$lp2[i]!=0) 
    { MMR_2stage_model$power2[i]=1
    } else if (MMR_2stage_model$pvalue_x2[i]>=0.05 & MMR_2stage_model$lp2[i]!=0)
    { MMR_2stage_model$power2[i]=0
    }
    
    if (MMR_2stage_model$pvalue_x3[i]<0.05 & MMR_2stage_model$lp3[i]!=0) 
    { MMR_2stage_model$power3[i]=1
    } else if (MMR_2stage_model$pvalue_x3[i]>=0.05 & MMR_2stage_model$lp3[i]!=0)
    { MMR_2stage_model$power3[i]=0
    }
    
    
  }  
  
  
  MMR_v1_mean=mean(MMR_2stage_model$mean_x1)
  MMR_v1_std=mean(MMR_2stage_model$se_x1)
  MMR_v1_coverage=mean(MMR_2stage_model$coverage1)
  MMR_v1_power=mean(MMR_2stage_model$power1)
  MMR_v1_fpr=mean(MMR_2stage_model$FPR1)
  
  MMR_v2_mean=mean(MMR_2stage_model$mean_x2)
  MMR_v2_std=mean(MMR_2stage_model$se_x2)
  MMR_v2_coverage=mean(MMR_2stage_model$coverage2)
  MMR_v2_power=mean(MMR_2stage_model$power2)
  MMR_v2_fpr=mean(MMR_2stage_model$FPR2)
  
  MMR_v3_mean=mean(MMR_2stage_model$mean_x3)
  MMR_v3_std=mean(MMR_2stage_model$se_x3)
  MMR_v3_coverage=mean(MMR_2stage_model$coverage3)
  MMR_v3_power=mean(MMR_2stage_model$power3)
  MMR_v3_fpr=mean(MMR_2stage_model$FPR3)
  
  sample_test$meanF[IV_loop]<-mean_F123
  
  sample_test$p1_mean[IV_loop]<-MMR_v1_mean
  sample_test$p2_mean[IV_loop]<-MMR_v2_mean
  sample_test$p3_mean[IV_loop]<-MMR_v3_mean
  
  sample_test$p1_se[IV_loop]<-MMR_v1_std
  sample_test$p2_se[IV_loop]<-MMR_v2_std
  sample_test$p3_se[IV_loop]<-MMR_v3_std
  
  sample_test$p1_cover[IV_loop]<-MMR_v1_coverage
  sample_test$p2_cover[IV_loop]<-MMR_v2_coverage
  sample_test$p3_cover[IV_loop]<-MMR_v3_coverage
  
  sample_test$p1_power[IV_loop]<-MMR_v1_power
  sample_test$p2_power[IV_loop]<-MMR_v2_power 
  sample_test$p3_power[IV_loop]<-MMR_v3_power
  
  sample_test$p1_FPR[IV_loop]<-MMR_v1_fpr
  sample_test$p2_FPR[IV_loop]<-MMR_v2_fpr
  sample_test$p3_FPR[IV_loop]<-MMR_v3_fpr
  

}


write.csv(sample_test,"./results/tableS1/ld1.BP.sample.csv")


