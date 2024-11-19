########################################
##LMM_MEC is an R macro to run multivariable Mendelian randomization (MVMR) using summary statistics 
##R code author: Ming Ding
##Publication: Ming Ding, Fei Zou. A Linear Mixed Model with Measurement Error Correction (LMM-MEC): A Method for Summary-data-based Multivariable Mendelian Randomization
##Contact: ming_ding@med.unc.edu
########################################


##install the following packages before running the macro

library(matrixcalc)  ##check singular.matrix
library(MASS)  ##Moore-Penrose inverse
library(mvtnorm)  ##simulate multivariate normal distribution
library(numDeriv)  ##jacobian matrix


##betaX: a N*M matrix of point estimates of summary statistics of risk factors, 
##########where N is the no. of genetic variants, and M is the no. of risk factors

##betaX_se: a N*M matrix of variances of summary statistics of risk factors

##betaY: a N*1 matrix of point estimates of summary statistics of disease 

##betaY_se: a N*1 matrix of variances of summary statistics of disease 

##corr_snps_lmmmec: a N*N matrix of correlations between genetic variants, 

##corr_X_lmmmec: a M*M matrix of correlations between summary statistics of risk factors
##can be estimated based on the data

##loop_rem: loop of the iteratively re-weighted least squares algorithm in stage 1. 

##cutoff_rem: criteria for the loop to end, suggested cutoff_rem=0.00001 

##mec_loop: number of validation studies simulated for measurement error correction


LMM_MEC=function(betaX,betaY,betaX_se, betaY_se,corr_snps_lmmmec,corr_X_lmmmec,loop_rem,cutoff_rem, mec_loop)
{
  
  betaX<-as.matrix(betaX)
  betaY<-as.matrix(betaY)
  betaX_se<-as.matrix(betaX_se)
  betaY_se<-as.matrix(betaY_se)
  corr_snps_lmmmec<-as.matrix(corr_snps_lmmmec)
  corr_X_lmmmec<-as.matrix(corr_X_lmmmec)
  
  betaX_var<-matrix(NA, nrow=nrow(betaX_se), ncol=ncol(betaX_se))
  
  for (i_var_X in 1:nrow(betaX_se))
  {
    for (j_var_X in 1:ncol(betaX_se))
    {
      betaX_var[i_var_X,j_var_X]<-betaX_se[i_var_X,j_var_X]^2
    }
  }
  
  betaY_var<-matrix(NA, nrow=nrow(betaY_se), ncol=1)
  
  for (i_var_Y in 1:nrow(betaY_se))
  {
    betaY_var[i_var_Y,1]<-betaY_se[i_var_Y,1]^2
  }
  
  intercept<-matrix(1,nrow=nrow(betaX),ncol=1)
  
  betaX_new=cbind(betaX,intercept)
  
  ##############################################################################################
  ###I. obtain lambda for measurement error correction in the summary statistics of risk factors
  
  #simulate a validation study
  
  matrix_s<-matrix(NA, nrow = nrow(betaX)*ncol(betaX), ncol = nrow(betaX)*ncol(betaX))
  
  
  for (u1 in 1:ncol(betaX))
  {
    for (u2 in 1:ncol(betaX))
    {
      
      
      for (i_svar in 1:nrow(betaX))
      {
        for (j_svar in 1:nrow(betaX))
        {
          matrix_s[nrow(betaX)*(u1-1)+i_svar,nrow(betaX)*(u2-1)+j_svar]<-
            sqrt(betaX_var[i_svar,u1])*sqrt(betaX_var[j_svar,u2])*corr_snps_lmmmec[i_svar,j_svar]*corr_X_lmmmec[u1,u2]
        }
      }
      
      
    }    
  }
  
  
  mean_s<-matrix(NA, ncol=1, nrow= nrow(betaX)*ncol(betaX))
  
  for (i_s in 1:nrow(betaX))
  {    
    for (j_s in 1:ncol(betaX))
    {
      mean_s[nrow(betaX)*(j_s-1)+i_s,1]<-betaX[i_s,j_s]
    }
  }
  
  nrow(mean_s)
  ncol(mean_s)
  
  #simulate the validation dataset
  
  lam_mean<- list()  
  lam_var<- list() 
  
  for (mec_simu in 1:mec_loop){
    
    #set.seed(mec_simu)
    
    exposure_s <- t(rmvnorm(1, mean=mean_s, sigma=matrix_s, method = "svd") ) ##sigma is the covariance for each variable
    
    nrow(exposure_s)
    ncol(exposure_s)
    
    # set.seed(Sys.time()) 
    
    betaX_simu=data.frame(betaX)  
    
    for (i_exp in 1:nrow(betaX))
    {    
      for (j_exp in 1:ncol(betaX))
      {
        betaX_simu[i_exp,paste0('s', j_exp)]<- exposure_s[nrow(betaX)*(j_exp-1)+i_exp,1]
      }
    }
    
    
    ##multivariate model to obtain lambda
    
    y_sim<-data.frame((as.matrix(betaX_simu[,paste0('s', 1:ncol(betaX))])))
    
    x_sim<-as.matrix(betaX_simu[,1:ncol(betaX)]) ##choose the first 3 columns
    
    mvmod <- lm(x_sim ~ . , data=y_sim)  
    
    lam_mean[[mec_simu]]<-(as.matrix(coef(mvmod)))[2:(ncol(betaX)+1),]
    
    lam_var[[mec_simu]]<-(vcov(mvmod))[-((ncol(betaX)+1)*((1:ncol(betaX))-1)+1),-((ncol(betaX)+1)*((1:ncol(betaX))-1)+1)]
    
  }
  
  lam_mean_pool<-Reduce("+",  lam_mean)/mec_loop
  
  lam_var_pool<-Reduce("+",  lam_var)/mec_loop
 
  ##obtain inverse lam
  
  #mean
  lam_inverse<-solve(lam_mean_pool)
  
  ##variance using delta method
  
  func_lam <- function(x_lam) {
    solve(x_lam)
  }
  
  x_lam <-lam_mean_pool   ##mean value of lam
  
  lam_jaco_pool<-jacobian(func_lam, x_lam)  ##jacobian matrix
  
  lam_var_inverse<-lam_jaco_pool%*%lam_var_pool%*%t(lam_jaco_pool)
  
  ##variables to be used for the next step
  
  lambda<-as.matrix(lam_inverse)  
  
  lambda_var<-as.matrix(lam_var_inverse)
  
  
  ##############################################################################################
  ######################II. choose random/fixed-effects model
  
  ##test for heterogeneity of betaY
  
  var_mtrx_f_ld<-matrix(NA, nrow = nrow(betaY), ncol = nrow(betaY))
  
  for (i_y in 1:nrow(betaY))
  {
    for (j_y in 1:nrow(betaY))
    {
      var_mtrx_f_ld[i_y,j_y]<-betaY_se[i_y]*betaY_se[j_y]*corr_snps_lmmmec[i_y,j_y]  ##change
    }
  }
  
  
  if (is.singular.matrix(var_mtrx_f_ld)==T){
    W_f_ld<-ginv(var_mtrx_f_ld)   
  } else {
    W_f_ld<-solve(var_mtrx_f_ld)
  }
  
  Q_statistic<-t(betaY-mean(betaY))%*%W_f_ld%*%(betaY-mean(betaY))
  
  I2<-(Q_statistic-(nrow(betaY)-1))/Q_statistic
  
  if (I2<0.5) {
    
    ####fixed-effects model accounting for LD
    
    beta_f_ld_cal<-t(betaX_new)%*%W_f_ld%*%betaX_new
    
    if (is.singular.matrix(beta_f_ld_cal)==T){
      beta_f_ld<-ginv(beta_f_ld_cal)%*%t(betaX_new)%*%W_f_ld%*%betaY
      
      SSE_f_ld<-t(betaY-betaX_new%*%beta_f_ld)%*%W_f_ld%*%(betaY-betaX_new%*%beta_f_ld)
      delta2_f_ld<-SSE_f_ld/(nrow(betaX_new)-ncol(betaX_new))
      
      beta_var_f_ld<-c(delta2_f_ld)*ginv(beta_f_ld_cal)
      
    } else {
      beta_f_ld<-solve(beta_f_ld_cal)%*%t(betaX_new)%*%W_f_ld%*%betaY
      
      SSE_f_ld<-t(betaY-betaX_new%*%beta_f_ld)%*%W_f_ld%*%(betaY-betaX_new%*%beta_f_ld)
      delta2_f_ld<-SSE_f_ld/(nrow(betaX_new)-ncol(betaX_new))
      
      beta_var_f_ld<-c(delta2_f_ld)*solve(beta_f_ld_cal)
    }
    
    
    ####fixed-effects model measurement error correction
    
    b_wo_mec<-as.matrix(beta_f_ld[1:ncol(betaX),]) ##column name of x is b1, b2, and intercept
    
    cov_wo_mec<-as.matrix(beta_var_f_ld[1:ncol(betaX),1:ncol(betaX)])
    
    ##given that k=k'*lamdata, b_wt_mec=t(b_wo_mec)%*%as.matrix(lambda)
    ##we transpose b_wt_mec, as it is easier to calculate the covariance matrix
    
    b_wt_mec<-t(as.matrix(lambda))%*%b_wo_mec
    
    ##no need to transpose variance matrix, as it is sysmetric
    
    lambda_new<-t(as.matrix(lambda)) 
    
    ##covariance matrix
    
    cov_wt_mec<-matrix(NA,nrow=ncol(betaX),ncol=ncol(betaX))
    
    for (m in 1:ncol(betaX))
    {
      for (n in 1:ncol(betaX))
      {
        
        item1<-matrix(NA,nrow=ncol(betaX), ncol=ncol(betaX))
        
        item2<-matrix(NA,nrow=ncol(betaX), ncol=ncol(betaX))
        
        item3<-matrix(NA,nrow=ncol(betaX), ncol=ncol(betaX))
        
        
        for (i_cov in 1:ncol(betaX))
        {
          for (j_cov in 1:ncol(betaX))
          {
            
            item1[i_cov,j_cov]<-lambda_var[m+(i_cov-1)*ncol(betaX), n+(j_cov-1)*ncol(betaX)]*cov_wo_mec[i_cov,j_cov]  
            
            item2[i_cov,j_cov]<-b_wo_mec[i_cov]*b_wo_mec[j_cov]*lambda_var[m+(i_cov-1)*ncol(betaX),n+(j_cov-1)*ncol(betaX)]
            
            item3[i_cov,j_cov]<-lambda_new[m,i_cov]*lambda_new[n,j_cov]*cov_wo_mec[i_cov,j_cov]  
            
          }
        }
        
        item1_sum<-colSums(matrix(rowSums(item1))) ##sum up all numbers in item1
        item2_sum<-colSums(matrix(rowSums(item2)))
        item3_sum<-colSums(matrix(rowSums(item3)))
        
        cov_wt_mec[m,n]<-item1_sum+item2_sum+item3_sum
        
      }
    }
    
    ##check whether the covariance matrix is symmetric
    cov_wt_mec
    
    ###calculate 95% CI
    
    ld_wt_mec<-data.frame(NA,nrow=ncol(betaX), ncol=10)
    
    for (i_out in 1:ncol(betaX))
    {
      
      ld_wt_mec[i_out,1]<-b_wt_mec[i_out]
      ld_wt_mec[i_out,2]<-sqrt(cov_wt_mec[i_out,i_out])
      ld_wt_mec[i_out,3]<-2*pnorm(-abs(b_wt_mec[i_out]/sqrt(cov_wt_mec[i_out,i_out])),mean=0, sd=1)
      ld_wt_mec[i_out,4]<-b_wt_mec[i_out]-1.96*sqrt(cov_wt_mec[i_out,i_out])
      ld_wt_mec[i_out,5]<-b_wt_mec[i_out]+1.96*sqrt(cov_wt_mec[i_out,i_out])
      
      ld_wt_mec[i_out,6]<-b_wo_mec[i_out]
      ld_wt_mec[i_out,7]<-sqrt(cov_wo_mec[i_out,i_out])
      ld_wt_mec[i_out,8]<-2*pnorm(-abs(b_wo_mec[i_out]/sqrt(cov_wo_mec[i_out,i_out])),mean=0, sd=1)
      ld_wt_mec[i_out,9]<-b_wo_mec[i_out]-1.96*sqrt(cov_wo_mec[i_out,i_out])
      ld_wt_mec[i_out,10]<-b_wo_mec[i_out]+1.96*sqrt(cov_wo_mec[i_out,i_out])
      
    }
    
    colnames(ld_wt_mec)=rbind("Mean_wt_mec", "SE_wt_mec","Pvalue_wt_mec", "ll_wt_mec", "ul_wt_mec",
                              "Mean_wo_mec", "SE_wo_mec","Pvalue_wo_mec", "ll_wo_mec", "ul_wo_mec")
    
    print("Fixed-effects model was used")
    
    ld_wt_mec    ##output results
    
    
  }  else if (I2>=0.5) {
    
    ## random-effects model accounting for LD
    
    A_track<-matrix(NA, nrow = loop_rem, ncol = 1)
    
    A<-0
    
    A_track[1,1]<-0
    
    for (i in 2:loop_rem)
    {
      
      ##calculate  weight
      
      var_mtrx_r<-matrix(NA, nrow = nrow(betaY), ncol = nrow(betaY))
      
      for (i_random in 1:nrow(betaY))
      {
        for (j_random in 1:nrow(betaY))
        {
          var_mtrx_r[i_random,j_random]<-sqrt(betaY_var[i_random]+A)*sqrt(betaY_var[j_random]+A)*corr_snps_lmmmec[i_random,j_random]
        }
      }
      
      if (is.singular.matrix(var_mtrx_r)==T){
        W_r<-ginv(var_mtrx_r)
      } else {
        W_r<-solve(var_mtrx_r)
      }
      
      
      beta_r_ld<-solve(t(betaX_new)%*%W_r%*%betaX_new)%*%t(betaX_new)%*%W_r%*%betaY
      
      ##obtain updated A
      
      A_s1<-(nrow(betaY)/(nrow(betaY)-ncol(betaX_new)))*((betaY-betaX_new%*%beta_r_ld)^2-betaY_var)  
      
      A_s2<-matrix(NA, nrow = nrow(betaY), ncol = 1)
      
      for (u in 1:nrow(betaY))
      {
        A_s2[u,1]<-W_r[u,u]*A_s1[u,1]
      }
      
      A<-abs(sum(A_s2)/sum(diag(W_r)))
      
      A_track[i,1]<-A
      
      ###set a criteria for the loop to end
      
      if (abs(A_track[i,1]-A_track[i-1,1])<=cutoff_rem)
      {
        break
      }
      
      stop=F
      
      if (i==loop_rem & abs(A_track[i,1]-A_track[i-1,1])>cutoff_rem)
      {
        stop=T
      }
      
      
    }
    
    
    
    beta_r_ld
    
    SSE_r_ld<-t(betaY-betaX_new%*%beta_r_ld)%*%W_r%*%(betaY-betaX_new%*%beta_r_ld)
    delta2_r_ld<-SSE_r_ld/(nrow(betaX_new)-ncol(betaX_new))
    beta_var_r_ld<-c(delta2_r_ld)*solve(t(betaX_new)%*%W_r%*%betaX_new)
    
    ####random-effects model for measurement error correction
    
    b_wo_mec<-as.matrix(beta_r_ld[1:ncol(betaX),]) ##column name of x is b1, b2, and intercept
    
    cov_wo_mec<-as.matrix(beta_var_r_ld[1:ncol(betaX),1:ncol(betaX)])
    
    ##given that k=k'*lamdata, b_wt_mec=t(b_wo_mec)%*%as.matrix(lambda)
    ##we transpose b_wt_mec, as it is easier to calculate the covariance matrix
    
    b_wt_mec<-t(as.matrix(lambda))%*%b_wo_mec
    
    ##no need to transpose the covariance matrix, as it is sysmetric
    
    lambda_new<-t(as.matrix(lambda)) 
    
    ##covariance
    
    cov_wt_mec<-matrix(NA,nrow=ncol(betaX),ncol=ncol(betaX))
    
    for (m in 1:ncol(betaX))
    {
      for (n in 1:ncol(betaX))
      {
        
        item1<-matrix(NA,nrow=ncol(betaX), ncol=ncol(betaX))
        
        item2<-matrix(NA,nrow=ncol(betaX), ncol=ncol(betaX))
        
        item3<-matrix(NA,nrow=ncol(betaX), ncol=ncol(betaX))
        
        
        for (i_cov in 1:ncol(betaX))
        {
          for (j_cov in 1:ncol(betaX))
          {
            
            item1[i_cov,j_cov]<-lambda_var[m+(i_cov-1)*ncol(betaX),n+(j_cov-1)*ncol(betaX)]*cov_wo_mec[i_cov,j_cov]  
            
            item2[i_cov,j_cov]<-b_wo_mec[i_cov]*b_wo_mec[j_cov]*lambda_var[m+(i_cov-1)*ncol(betaX),n+(j_cov-1)*ncol(betaX)]
            
            item3[i_cov,j_cov]<-lambda_new[m,i_cov]*lambda_new[n,j_cov]*cov_wo_mec[i_cov,j_cov]  
            
          }
        }
        
        item1_sum<-colSums(matrix(rowSums(item1))) ##sum up all numbers in item1
        item2_sum<-colSums(matrix(rowSums(item2)))
        item3_sum<-colSums(matrix(rowSums(item3)))
        
        cov_wt_mec[m,n]<-item1_sum+item2_sum+item3_sum
        
      }
    }
    
    ##check whether the covariance matrix is symmetric
    cov_wt_mec
    
    ###calculate 95% CI
    
    ld_wt_mec<-data.frame(NA,nrow=ncol(betaX), ncol=10)
    
    for (i_out in 1:ncol(betaX))
    {
      
      ld_wt_mec[i_out,1]<-b_wt_mec[i_out]
      ld_wt_mec[i_out,2]<-sqrt(cov_wt_mec[i_out,i_out])
      ld_wt_mec[i_out,3]<-2*pnorm(-abs(b_wt_mec[i_out]/sqrt(cov_wt_mec[i_out,i_out])),mean=0, sd=1)
      ld_wt_mec[i_out,4]<-b_wt_mec[i_out]-1.96*sqrt(cov_wt_mec[i_out,i_out])
      ld_wt_mec[i_out,5]<-b_wt_mec[i_out]+1.96*sqrt(cov_wt_mec[i_out,i_out])
      
      ld_wt_mec[i_out,6]<-b_wo_mec[i_out]
      ld_wt_mec[i_out,7]<-sqrt(cov_wo_mec[i_out,i_out])
      ld_wt_mec[i_out,8]<-2*pnorm(-abs(b_wo_mec[i_out]/sqrt(cov_wo_mec[i_out,i_out])),mean=0, sd=1)
      ld_wt_mec[i_out,9]<-b_wo_mec[i_out]-1.96*sqrt(cov_wo_mec[i_out,i_out])
      ld_wt_mec[i_out,10]<-b_wo_mec[i_out]+1.96*sqrt(cov_wo_mec[i_out,i_out])
      
    }
    
    colnames(ld_wt_mec)=rbind("Mean_wt_mec", "SE_wt_mec","Pvalue_wt_mec", "ll_wt_mec", "ul_wt_mec",
                              "Mean_wo_mec", "SE_wo_mec","Pvalue_wo_mec", "ll_wo_mec", "ul_wo_mec")
    
    rownames(ld_wt_mec)=paste0("variable",1:ncol(betaX))
    
    # if   (stop)  { 
    #   print("random-effects model does not converge")}
    
    print("Random-effects model was used")
    
    ld_wt_mec
  }
  
  
}

