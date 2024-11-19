setwd("/Users/mingding/upload to github")

##download and read in LMM_MEC macro
source('./LMM_MEC.macro.r')

##read in dataset
new_data<-read.csv("./data.csv", row.names = NULL)

##corr_snps is correlation between snps
#Corr(β_1, β_2) can be approximated by Pearson correlation (h_12) between genetic variants, and this can be estimated in two steps in practice. 
#First, h_12^2, which is LD between genetic variants, can be estimated using LDmatrix Tool in LDlink. 
#Second, effective allele frequency (EAF) is used to decide positive or negative square root of LD,  
#with a positive value assigned if the absolute difference in EAF of two genetic variants is smaller than the absolute difference 
#between EAF of one genetic variant and minor allele frequency (MAF) of the other genetic variant.

corr_snps_readin<-read.csv("./corr_snps.csv", row.names = NULL)
corr_snps<- corr_snps_readin[,-1]  ##exclude the first column which is a character variable


##corr_x is correlation between risk factors. 
##This can be estimated as correlation between phenotypes using individual-level data, 
#or correlation between summary statistics of risk factors as a surrogate

corr_x<-cor(new_data[, c('bx1', 'bx2', 'bx3')])
  
##apply LMM_MEC
  
  model_output<-LMM_MEC(
    betaX=cbind(new_data$bx1, new_data$bx2, new_data$bx3),
    betaY=new_data$by, 
    betaX_se=cbind(new_data$bx_se1, new_data$bx_se2, new_data$bx_se3),
    betaY_se=new_data$by_se,
    corr_snps_lmmmec =corr_snps,
    corr_X_lmmmec =corr_x,
    loop_rem=500,
    cutoff_rem=0.00001,
    mec_loop=1
  )
  
  model_output
  
