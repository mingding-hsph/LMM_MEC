
library(mvtnorm)

total_n=2*1000   ##number of participants
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

snp_con_1<-rmvnorm(total_n, mean =mean_snp, sigma =corr_0, method = "svd") ##sigma is the covariance matrix.
snp_con_2<-rmvnorm(total_n, mean =mean_snp, sigma =corr_0, method = "svd") ##sigma is the covariance matrix.

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

lp<-matrix(nrow=exp_num)  

for (i in 1:exp_num)
{

  lp[1]=0.1
  lp[2]=0
  lp[3]=-0.1
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

##new data contains summary statistics of gx and gy estimated from two independent datasets
#In the new_data, by is the outcome for MVMR, bx1, bx2, bx3 are the exposures for MVMR

#estimate correlation of summary statistics of Y, which will be used for MVMR 
##the lmx_x and lmx_y are simulated with covariance matrix corr_snps

corr_snps


corr_x<-cor(data[,c('a1','a2','a3')])

##check

colMeans(new_data_f)

colMeans(new_data_r2)

corr_x

summary(lm(by~bx1+bx2+bx3,data=new_data)) 

cor(t(lmk_x))

