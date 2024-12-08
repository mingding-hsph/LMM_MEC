# set working directory here
this_dir <- function(directory)
  setwd( file.path(getwd(), directory) )

getwd()

##ld1
MRMV_model=read.csv("./ld1.noP.egger.csv", row.names = NULL)
MRMV_IVW_model=read.csv("./ld1.noP.IVW.csv", row.names = NULL)
MMR_2stage_model=read.csv("./ld1.noP.LMM.MEC.csv", row.names = NULL)
scenario=cbind(0.1,0,-0.1)  

MRMV_model=read.csv("./ld1.BP.egger.csv", row.names = NULL)
MRMV_IVW_model=read.csv("./ld1.BP.IVW.csv", row.names = NULL)
MMR_2stage_model=read.csv("./ld1.BP.LMM.MEC.csv", row.names = NULL)
scenario=cbind(0.3,0,-0.3)  


MRMV_model=read.csv("./ld1.UBP.egger.csv", row.names = NULL)
MRMV_IVW_model=read.csv("./ld1.UBP.IVW.csv", row.names = NULL)
MMR_2stage_model=read.csv("./ld1.UBP.LMM.MEC.csv", row.names = NULL)
scenario=cbind(0.1,0,-0.1)  

##ld3
MRMV_model=read.csv("./ld3.noP.egger.csv", row.names = NULL)
MRMV_IVW_model=read.csv("./ld3.noP.IVW.csv", row.names = NULL)
MMR_2stage_model=read.csv("./ld3.noP.LMM.MEC.csv", row.names = NULL)
scenario=cbind(0.1,0,-0.1)  

MRMV_model=read.csv("./ld3.BP.egger.csv", row.names = NULL)
MRMV_IVW_model=read.csv("./ld3.BP.IVW.csv", row.names = NULL)
MMR_2stage_model=read.csv("./ld3.BP.LMM.MEC.csv", row.names = NULL)
scenario=cbind(0.3,0,-0.3)  

MRMV_model=read.csv("./ld3.UBP.egger.csv", row.names = NULL)
MRMV_IVW_model=read.csv("./ld3.UBP.IVW.csv", row.names = NULL)
MMR_2stage_model=read.csv("./ld3.UBP.LMM.MEC.csv", row.names = NULL)
scenario=cbind(0.1,0,-0.1)  






MRMV_model$lp1=scenario[1,1]
MRMV_model$lp2=scenario[1,2]
MRMV_model$lp3=scenario[1,3]

MRMV_model$coverage1=NA
MRMV_model$power1=NA
MRMV_model$FPR1=NA
MRMV_model$mean1=NA


MRMV_model$coverage2=NA
MRMV_model$power2=NA
MRMV_model$FPR2=NA
MRMV_model$mean2=NA

MRMV_model$coverage3=NA
MRMV_model$power3=NA
MRMV_model$FPR3=NA
MRMV_model$mean3=NA

for (i in 1:nrow(MRMV_model))
{  

  ##Coverage: proportion of 95% CIs including the true causal effect over all simulated datasets; 
  
  if (MRMV_model$lp1[i]>=MRMV_model$ll_x1[i] & MRMV_model$lp1[i]<=MRMV_model$ul_x1[i])
  { MRMV_model$coverage1[i]=1
  } else {
    MRMV_model$coverage1[i]=0
  }
  
  if (MRMV_model$lp2[i]>=MRMV_model$ll_x2[i] & MRMV_model$lp2[i]<=MRMV_model$ul_x2[i])
  { MRMV_model$coverage2[i]=1
  } else {
    MRMV_model$coverage2[i]=0
  }
  
  
  if (MRMV_model$lp3[i]>=MRMV_model$ll_x3[i] & MRMV_model$lp3[i]<=MRMV_model$ul_x3[i])
  { MRMV_model$coverage3[i]=1
  } else {
    MRMV_model$coverage3[i]=0
  }
  
  
  ##false-positive rate: the proportion of significant causal estimates in settings with no causal effect. 
  
  if (MRMV_model$pvalue_x1[i]<0.05 & MRMV_model$lp1[i]==0) 
  { MRMV_model$FPR1[i]=1
  } else if (MRMV_model$pvalue_x1[i]>=0.05 & MRMV_model$lp1[i]==0) {
    MRMV_model$FPR1[i]=0  
  }
  
  if (MRMV_model$pvalue_x2[i]<0.05 & MRMV_model$lp2[i]==0) 
  { MRMV_model$FPR2[i]=1
  } else if (MRMV_model$pvalue_x2[i]>=0.05 & MRMV_model$lp2[i]==0) {
    MRMV_model$FPR2[i]=0  
  }

  
  if (MRMV_model$pvalue_x3[i]<0.05 & MRMV_model$lp3[i]==0) 
  { MRMV_model$FPR3[i]=1
  } else if (MRMV_model$pvalue_x3[i]>=0.05 & MRMV_model$lp3[i]==0) {
    MRMV_model$FPR3[i]=0  
  }
  

  ##power is estimated by the proportion of 95% CIs excluding zero in settings with significant causal effect; 
  
  if (MRMV_model$pvalue_x1[i]<0.05 & MRMV_model$lp1[i]!=0) 
  { MRMV_model$power1[i]=1
    } else if (MRMV_model$pvalue_x1[i]>=0.05 & MRMV_model$lp1[i]!=0)
    { MRMV_model$power1[i]=0
    }
  
  if (MRMV_model$pvalue_x2[i]<0.05 & MRMV_model$lp2[i]!=0) 
  { MRMV_model$power2[i]=1
  } else if (MRMV_model$pvalue_x2[i]>=0.05 & MRMV_model$lp2[i]!=0)
  { MRMV_model$power2[i]=0
  }
  
  if (MRMV_model$pvalue_x3[i]<0.05 & MRMV_model$lp3[i]!=0) 
  { MRMV_model$power3[i]=1
  } else if (MRMV_model$pvalue_x3[i]>=0.05 & MRMV_model$lp3[i]!=0)
  { MRMV_model$power3[i]=0
  }
  
  
  }

mean(MRMV_model$coverage1)
mean(MRMV_model$power1)
mean(MRMV_model$FPR1)

mean(MRMV_model$coverage2)
mean(MRMV_model$power2)
mean(MRMV_model$FPR2)

mean(MRMV_model$coverage3)
mean(MRMV_model$power3)
mean(MRMV_model$FPR3)


##ivw

##change the numbers here

MRMV_IVW_model$lp1=scenario[1,1]
MRMV_IVW_model$lp2=scenario[1,2]
MRMV_IVW_model$lp3=scenario[1,3]

MRMV_IVW_model$coverage1=NA
MRMV_IVW_model$power1=NA
MRMV_IVW_model$FPR1=NA
MRMV_IVW_model$mean1=NA

MRMV_IVW_model$coverage2=NA
MRMV_IVW_model$power2=NA
MRMV_IVW_model$FPR2=NA
MRMV_IVW_model$mean2=NA

MRMV_IVW_model$coverage3=NA
MRMV_IVW_model$power3=NA
MRMV_IVW_model$FPR3=NA
MRMV_IVW_model$mean3=NA

for (i in 1:nrow(MRMV_IVW_model))
{  
  
  ##Coverage: proportion of 95% CIs including the true causal effect over all simulated datasets; 
  
  if (MRMV_IVW_model$lp1[i]>=MRMV_IVW_model$ll_x1[i] & MRMV_IVW_model$lp1[i]<=MRMV_IVW_model$ul_x1[i])
  { MRMV_IVW_model$coverage1[i]=1
  } else {
    MRMV_IVW_model$coverage1[i]=0
  }
  
  if (MRMV_IVW_model$lp2[i]>=MRMV_IVW_model$ll_x2[i] & MRMV_IVW_model$lp2[i]<=MRMV_IVW_model$ul_x2[i])
  { MRMV_IVW_model$coverage2[i]=1
  } else {
    MRMV_IVW_model$coverage2[i]=0
  }
  
  
  if (MRMV_IVW_model$lp3[i]>=MRMV_IVW_model$ll_x3[i] & MRMV_IVW_model$lp3[i]<=MRMV_IVW_model$ul_x3[i])
  { MRMV_IVW_model$coverage3[i]=1
  } else {
    MRMV_IVW_model$coverage3[i]=0
  }
  
  
  ##false-positive rate: the proportion of significant causal estimates in settings with no causal effect. 
  
  if (MRMV_IVW_model$pvalue_x1[i]<0.05 & MRMV_IVW_model$lp1[i]==0) 
  { MRMV_IVW_model$FPR1[i]=1
  } else if (MRMV_IVW_model$pvalue_x1[i]>=0.05 & MRMV_IVW_model$lp1[i]==0) {
    MRMV_IVW_model$FPR1[i]=0  
  }
  
  if (MRMV_IVW_model$pvalue_x2[i]<0.05 & MRMV_IVW_model$lp2[i]==0) 
  { MRMV_IVW_model$FPR2[i]=1
  } else if (MRMV_IVW_model$pvalue_x2[i]>=0.05 & MRMV_IVW_model$lp2[i]==0) {
    MRMV_IVW_model$FPR2[i]=0  
  }
  
  
  if (MRMV_IVW_model$pvalue_x3[i]<0.05 & MRMV_IVW_model$lp3[i]==0) 
  { MRMV_IVW_model$FPR3[i]=1
  } else if (MRMV_IVW_model$pvalue_x3[i]>=0.05 & MRMV_IVW_model$lp3[i]==0) {
    MRMV_IVW_model$FPR3[i]=0  
  }
  
  
  ##power is estimated by the proportion of 95% CIs excluding zero in settings with significant causal effect; 
  
  if (MRMV_IVW_model$pvalue_x1[i]<0.05 & MRMV_IVW_model$lp1[i]!=0) 
  { MRMV_IVW_model$power1[i]=1
  } else if (MRMV_IVW_model$pvalue_x1[i]>=0.05 & MRMV_IVW_model$lp1[i]!=0)
  { MRMV_IVW_model$power1[i]=0
  }
  
  if (MRMV_IVW_model$pvalue_x2[i]<0.05 & MRMV_IVW_model$lp2[i]!=0) 
  { MRMV_IVW_model$power2[i]=1
  } else if (MRMV_IVW_model$pvalue_x2[i]>=0.05 & MRMV_IVW_model$lp2[i]!=0)
  { MRMV_IVW_model$power2[i]=0
  }
  
  if (MRMV_IVW_model$pvalue_x3[i]<0.05 & MRMV_IVW_model$lp3[i]!=0) 
  { MRMV_IVW_model$power3[i]=1
  } else if (MRMV_IVW_model$pvalue_x3[i]>=0.05 & MRMV_IVW_model$lp3[i]!=0)
  { MRMV_IVW_model$power3[i]=0
  }
  
  
}

mean(MRMV_IVW_model$coverage1)
mean(MRMV_IVW_model$power1)
mean(MRMV_IVW_model$FPR1)

mean(MRMV_IVW_model$coverage2)
mean(MRMV_IVW_model$power2)
mean(MRMV_IVW_model$FPR2)

mean(MRMV_IVW_model$coverage3)
mean(MRMV_IVW_model$power3)
mean(MRMV_IVW_model$FPR3)



###our model

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

mean(MMR_2stage_model$coverage1)
mean(MMR_2stage_model$power1)
mean(MMR_2stage_model$FPR1)

mean(MMR_2stage_model$coverage2)
mean(MMR_2stage_model$power2)
mean(MMR_2stage_model$FPR2)

mean(MMR_2stage_model$coverage3)
mean(MMR_2stage_model$power3)
mean(MMR_2stage_model$FPR3)


##output into tables

##install.packages("plotrix")                           # Install plotrix R package
library("plotrix")                                    # to calculate standard error


MRMV_v1_mean=mean(MRMV_model$mean_x1)
MRMV_v1_std=mean(MRMV_model$se_x1)
MRMV_v1_coverage=mean(MRMV_model$coverage1)
MRMV_v1_power=mean(MRMV_model$power1)
MRMV_v1_fpr=mean(MRMV_model$FPR1)

MRMV_v2_mean=mean(MRMV_model$mean_x2)
MRMV_v2_std=mean(MRMV_model$se_x2)
MRMV_v2_coverage=mean(MRMV_model$coverage2)
MRMV_v2_power=mean(MRMV_model$power2)
MRMV_v2_fpr=mean(MRMV_model$FPR2)

MRMV_v3_mean=mean(MRMV_model$mean_x3)
MRMV_v3_std=mean(MRMV_model$se_x3)
MRMV_v3_coverage=mean(MRMV_model$coverage3)
MRMV_v3_power=mean(MRMV_model$power3)
MRMV_v3_fpr=mean(MRMV_model$FPR3)

MRMV_IVW_v1_mean=mean(MRMV_IVW_model$mean_x1)
MRMV_IVW_v1_std=mean(MRMV_IVW_model$se_x1)
MRMV_IVW_v1_coverage=mean(MRMV_IVW_model$coverage1)
MRMV_IVW_v1_power=mean(MRMV_IVW_model$power1)
MRMV_IVW_v1_fpr=mean(MRMV_IVW_model$FPR1)

MRMV_IVW_v2_mean=mean(MRMV_IVW_model$mean_x2)
MRMV_IVW_v2_std=mean(MRMV_IVW_model$se_x2)
MRMV_IVW_v2_coverage=mean(MRMV_IVW_model$coverage2)
MRMV_IVW_v2_power=mean(MRMV_IVW_model$power2)
MRMV_IVW_v2_fpr=mean(MRMV_IVW_model$FPR2)

MRMV_IVW_v3_mean=mean(MRMV_IVW_model$mean_x3)
MRMV_IVW_v3_std=mean(MRMV_IVW_model$se_x3)
MRMV_IVW_v3_coverage=mean(MRMV_IVW_model$coverage3)
MRMV_IVW_v3_power=mean(MRMV_IVW_model$power3)
MRMV_IVW_v3_fpr=mean(MRMV_IVW_model$FPR3)


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



v1_mean=round(rbind(MMR_v1_mean, MRMV_v1_mean, MRMV_IVW_v1_mean),digits=2)
v1_std=round(rbind(MMR_v1_std, MRMV_v1_std, MRMV_IVW_v1_std),digits=2)
v1_coverage=100*rbind(MMR_v1_coverage, MRMV_v1_coverage, MRMV_IVW_v1_coverage)
v1_power=100*rbind(MMR_v1_power, MRMV_v1_power, MRMV_IVW_v1_power)
v1_FPR=100*rbind(MMR_v1_fpr, MRMV_v1_fpr, MRMV_IVW_v1_fpr)

v2_mean=round(rbind(MMR_v2_mean, MRMV_v2_mean, MRMV_IVW_v2_mean),digits=2)
v2_std=round(rbind(MMR_v2_std, MRMV_v2_std, MRMV_IVW_v2_std),digits=2)
v2_coverage=100*rbind(MMR_v2_coverage, MRMV_v2_coverage, MRMV_IVW_v2_coverage)
v2_power=100*rbind(MMR_v2_power, MRMV_v2_power, MRMV_IVW_v2_power)
v2_FPR=100*rbind(MMR_v2_fpr, MRMV_v2_fpr, MRMV_IVW_v2_fpr)

v3_mean=round(rbind(MMR_v3_mean, MRMV_v3_mean, MRMV_IVW_v3_mean),digits=2)
v3_std=round(rbind(MMR_v3_std, MRMV_v3_std, MRMV_IVW_v3_std),digits=2)
v3_coverage=100*rbind(MMR_v3_coverage, MRMV_v3_coverage, MRMV_IVW_v3_coverage)
v3_power=100*rbind(MMR_v3_power, MRMV_v3_power, MRMV_IVW_v3_power)
v3_FPR=100*rbind(MMR_v3_fpr, MRMV_v3_fpr, MRMV_IVW_v3_fpr)


combined_table=cbind(rbind(v1_mean,v1_std,v1_coverage,v1_power,v1_FPR),
                     rbind(v2_mean,v2_std,v2_coverage,v2_power,v2_FPR),
                     rbind(v3_mean,v3_std,v3_coverage,v3_power,v3_FPR))

write.csv(combined_table,"./a.table1_ld3_UBP.csv")
