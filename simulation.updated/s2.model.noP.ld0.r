
source('./LMM_MEC.macro.r')


loop<-1000

qhet_model<-matrix(NA,nrow=loop, ncol=5*3) 
MRMV_model<-matrix(NA,nrow=loop, ncol=5*3) 
MRMV_IVW_model<-matrix(NA,nrow=loop, ncol=5*3) 
MRMV_MEDIAN_model<-matrix(NA,nrow=loop, ncol=5*3) 
MMR_2stage_model<-matrix(NA,nrow=loop, ncol=10*3) 
presso_model<-matrix(NA,nrow=loop, ncol=5*3) 
IV_strength<-matrix(NA,nrow=loop, ncol=3) 
IV_con_strength<-matrix(NA,nrow=loop, ncol=3) ##this include conditional IV
incept_pvalue<-matrix(NA,nrow=loop, ncol=1)

colnames(MRMV_model)<-rbind("mean_x1","mean_x2","mean_x3","se_x1","se_x2","se_x3",
                            "pvalue_x1", "pvalue_x2","pvalue_x3",
                            "ll_x1","ll_x2","ll_x3", "ul_x1",  "ul_x2", "ul_x3")
colnames(MRMV_IVW_model)<-rbind("mean_x1","mean_x2","mean_x3","se_x1","se_x2","se_x3",
                                "pvalue_x1", "pvalue_x2","pvalue_x3",
                                "ll_x1","ll_x2","ll_x3", "ul_x1",  "ul_x2", "ul_x3")
colnames(MRMV_MEDIAN_model)<-rbind("mean_x1","mean_x2","mean_x3","se_x1","se_x2","se_x3",
                                   "pvalue_x1", "pvalue_x2","pvalue_x3",
                                   "ll_x1","ll_x2","ll_x3", "ul_x1",  "ul_x2", "ul_x3")
colnames(MMR_2stage_model)<-rbind("mean_x1","mean_x2","mean_x3","se_x1","se_x2","se_x3",
                                  "pvalue_x1", "pvalue_x2","pvalue_x3",
                                  "ll_x1","ll_x2","ll_x3", "ul_x1",  "ul_x2", "ul_x3",
                                  "mean_x1_wo","mean_x2_wo","mean_x3_wo","se_x1_wo","se_x2_wo","se_x3_wo",
                                  "pvalue_x1_wo", "pvalue_x2_wo","pvalue_x3_wo",
                                  "ll_x1_wo","ll_x2_wo","ll_x3_wo", "ul_x1_wo",  "ul_x2_wo", "ul_x3_wo")
colnames(presso_model)<-rbind("mean_x1","mean_x2","mean_x3","se_x1","se_x2","se_x3",
                              "pvalue_x1", "pvalue_x2","pvalue_x3",
                              "ll_x1","ll_x2","ll_x3", "ul_x1",  "ul_x2", "ul_x3")

colnames(incept_pvalue)<-rbind("Pvalue")
                        
for (simu_loop in 1:loop)
  
{
 
  source('./s1.data.ld0.noP.r')
  
  ##estimate mean strength of IV
  
  library(tidyverse)
  
  F1=mean((new_data$bx1)^2/(new_data$bx_se1)^2)
  
  F2=mean((new_data$bx2)^2/(new_data$bx_se2)^2)
  
  F3=mean((new_data$bx3)^2/(new_data$bx_se3)^2)
  
  IV_strength[simu_loop,1]<-F1
  IV_strength[simu_loop,2]<-F2
  IV_strength[simu_loop,3]<-F3
  
  ##check pleiotropy
  
  incept_pvalue[simu_loop,1]=summary(lm(by~bx1+bx2+bx3,data=new_data))$coeff[1,4]  ##p value of intercept
  
  
##install.packages("MendelianRandomization")

library(MendelianRandomization)

MRMV<-mr_mvegger(mr_mvinput(bx = cbind(new_data$bx1, new_data$bx2, new_data$bx3), 
                           bxse = cbind(new_data$bx_se1, new_data$bx_se2, new_data$bx_se3),
                           by = new_data$by, byse = new_data$by_se, correlation=corr_snps))

MRMV_IVW<-mr_mvivw(mr_mvinput(bx = cbind(new_data$bx1, new_data$bx2, new_data$bx3), 
                    bxse = cbind(new_data$bx_se1, new_data$bx_se2, new_data$bx_se3),
                    by = new_data$by, byse = new_data$by_se, correlation=corr_snps))

MRMV_MEDIAN<-mr_mvmedian(mr_mvinput(bx = cbind(new_data$bx1, new_data$bx2, new_data$bx3), 
                       bxse = cbind(new_data$bx_se1, new_data$bx_se2, new_data$bx_se3),
                       by = new_data$by, byse = new_data$by_se, correlation=corr_snps), 
                       iterations = 1000)

##output

MRMV_model[simu_loop,1:3]<-MRMV$Estimate
MRMV_model[simu_loop,(1:3)+3]<-MRMV$StdError.Est
MRMV_model[simu_loop,(1:3)+2*3]<-MRMV$Pvalue.Est
MRMV_model[simu_loop,(1:3)+3*3]<-MRMV$CILower.Est
MRMV_model[simu_loop,(1:3)+4*3]<-MRMV$CIUpper.Est

MRMV_IVW_model[simu_loop,1:3]<-MRMV_IVW$Estimate
MRMV_IVW_model[simu_loop,(1:3)+3]<-MRMV_IVW$StdError
MRMV_IVW_model[simu_loop,(1:3)+2*3]<-MRMV_IVW$Pvalue
MRMV_IVW_model[simu_loop,(1:3)+3*3]<-MRMV_IVW$CILower
MRMV_IVW_model[simu_loop,(1:3)+4*3]<-MRMV_IVW$CIUpper

MRMV_MEDIAN_model[simu_loop,1:3]<-MRMV_MEDIAN$Estimate
MRMV_MEDIAN_model[simu_loop,(1:3)+3]<-MRMV_MEDIAN$StdError
MRMV_MEDIAN_model[simu_loop,(1:3)+2*3]<-MRMV_MEDIAN$Pvalue
MRMV_MEDIAN_model[simu_loop,(1:3)+3*3]<-MRMV_MEDIAN$CILower
MRMV_MEDIAN_model[simu_loop,(1:3)+4*3]<-MRMV_MEDIAN$CIUpper

###MMR-2stage

##estimate conditional IV strength
#install.packages("MVMR")

#install.packages("remotes")
#remotes::install_github("WSpiller/MVMR")

library(MVMR)

con_F <- format_mvmr(BXGs =t(rbind(new_data$bx1,new_data$bx2,new_data$bx3)),
                     BYG = new_data$by,
                     seBXGs =t(rbind(new_data$bx_se1,new_data$bx_se2,new_data$bx_se3)),
                     seBYG = new_data$by_se,
                     RSID = seq(1:nrow(new_data)))

mvmrcovmatrix<-corr_x

Xcovmat<-phenocov_mvmr(mvmrcovmatrix, con_F[,7:9])

con_F_s <- strength_mvmr(r_input = con_F, gencov = Xcovmat)

IV_con_strength[simu_loop,1]=as.numeric(con_F_s[1])
IV_con_strength[simu_loop,2]=as.numeric(con_F_s[2])
IV_con_strength[simu_loop,3]=as.numeric(con_F_s[3])


##LMM_MEC

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

MMR_2stage_model[simu_loop,1:3+5*3]<-t(MMR[,6])
MMR_2stage_model[simu_loop,(1:3)+6*3]<-t(MMR[,7])
MMR_2stage_model[simu_loop,(1:3)+7*3]<-t(MMR[,8])
MMR_2stage_model[simu_loop,(1:3)+8*3]<-t(MMR[,9])
MMR_2stage_model[simu_loop,(1:3)+9*3]<-t(MMR[,10])

##install mrpresso
#if (!require("devtools")) { install.packages("devtools") } else {}
#devtools::install_github("rondolab/MR-PRESSO")

#install.packages("remotes")
#remotes::install_github("rondolab/MR-PRESSO")

library(MRPRESSO)

# Run MR-PRESSO on a multi-variable MR (MMR) model specifying several exposures

new_data_press<-data.frame(new_data)

MRPRESSO=mr_presso(BetaOutcome = "by", BetaExposure = c("bx1", "bx2", "bx3"), 
                   SdOutcome = "by_se", SdExposure = c("bx_se1", "bx_se2", "bx_se3"), 
                   OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = new_data_press, 
                   NbDistribution = 1000,  SignifThreshold = 0.05)

##checked, the warnings() are fine

presso_model[simu_loop,1:3]<-MRPRESSO$`Main MR results`[1:3,3]  ##estimate
presso_model[simu_loop,(1:3)+3]<-MRPRESSO$`Main MR results`[1:3,4]  ##although write sd, this shoudl be se
presso_model[simu_loop,(1:3)+2*3]<-MRPRESSO$`Main MR results`[1:3,6]  ##p value
presso_model[simu_loop,(1:3)+3*3]<-MRPRESSO$`Main MR results`[1:3,3]-1.96*MRPRESSO$`Main MR results`[1:3,4]
presso_model[simu_loop,(1:3)+4*3]<-MRPRESSO$`Main MR results`[1:3,3]+1.96*MRPRESSO$`Main MR results`[1:3,4]

write.csv(MRMV_model,"./results/table1/ld0.noP.egger.csv")
write.csv(MRMV_IVW_model,"./results/table1/ld0.noP.IVW.csv")
write.csv(MRMV_MEDIAN_model,"./results/table1/ld0.noP.MEDIAN.csv")
write.csv(MMR_2stage_model,"./results/table1/ld0.noP.LMM.MEC.csv")
write.csv(presso_model,"./results/table1/ld0.noP.presso.csv")

write.csv(IV_con_strength,"./results/table1/ld0.noP.IV_con.csv")
write.csv(IV_strength,"./results/table1/ld0.noP.IV.csv")
write.csv(incept_pvalue,"./results/table1/ld0.noP.incept.pvalue.csv")

}

