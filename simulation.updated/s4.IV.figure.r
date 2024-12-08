setwd("/Users/mingding/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/a2.R21.MR/Paper 1.MR method/Rcode/simulation.v4")

setwd("./results from longleaf/Figure2/")


##balanced P
IV_0<-read.csv("./BP.csv", row.names = NULL)
colnames(IV_0)

jpeg('IV_balanced_P.png', width = 7, height = 7, units = 'in', res = 600)

lo_p1_cover <- loess(IV_0$p1_cover~IV_0$meanF)
lo_p2_cover <- loess(IV_0$p2_cover~IV_0$meanF)
lo_p3_cover <- loess(IV_0$p3_cover~IV_0$meanF)
lo_p1_power <- loess(IV_0$p1_power~IV_0$meanF)
lo_p3_power <- loess(IV_0$p3_power~IV_0$meanF)
lo_p2_FPR <- loess(IV_0$p2_FPR~IV_0$meanF)

par(mfrow=c(1,1))

plot(main="Balanced pleiotropy", IV_0$meanF, 100*IV_0$p1_cover, type="p", pch=16, cex = 0.5, 
     col="black", xlab="Mean F statistic", ylab="Performances (%)",ylim=range(0, 100),  
     xlim=rev(range(0,20)))
lines(IV_0$meanF, 100*IV_0$p2_cover, type="p", pch=16, cex = 0.5, col="black")
lines(IV_0$meanF, 100*IV_0$p3_cover, type="p", pch=16, cex = 0.5, col="black")
lines(IV_0$meanF, 100*IV_0$p1_power, type="p", pch=16, cex = 0.5, col="red")
lines(IV_0$meanF, 100*IV_0$p3_power, type="p", pch=16, cex = 0.5, col="red")

lines(IV_0$meanF, 100*predict(lo_p1_cover), type="l", lty=1,lwd=2, col="black")
lines(IV_0$meanF, 100*predict(lo_p2_cover), type="l", lty=2,lwd=2, col="black")
lines(IV_0$meanF, 100*predict(lo_p3_cover), type="l", lty=3,lwd=2, col="black")
lines(IV_0$meanF, 100*predict(lo_p1_power), type="l", lty=1,lwd=2, col="red")
lines(IV_0$meanF, 100*predict(lo_p3_power), type="l", lty=2,lwd=2, col="red")

legend(20, 70, legend=c("k1 coverage rate", "k2 coverage rate", 
                         "k3 coverage rate", "k1 power", 
                         "k3 power"),
       col=c("black", "black", "black", "red", "red"), 
       lty=c(1,2,3,1,2,1),cex=1, bty="n")

dev.off()

####no P

IV_0<-read.csv("./noP.csv", row.names = NULL)
colnames(IV_0)


jpeg('IV_noP.png', width = 7, height = 7, units = 'in', res = 600)

lo_p1_cover <- loess(IV_0$p1_cover~IV_0$meanF)
lo_p2_cover <- loess(IV_0$p2_cover~IV_0$meanF)
lo_p3_cover <- loess(IV_0$p3_cover~IV_0$meanF)
lo_p1_power <- loess(IV_0$p1_power~IV_0$meanF)
lo_p3_power <- loess(IV_0$p3_power~IV_0$meanF)
lo_p2_FPR <- loess(IV_0$p2_FPR~IV_0$meanF)

par(mfrow=c(1,1))

plot(main="No pleiotropy with disease outcome", IV_0$meanF, 100*IV_0$p1_cover, type="p", pch=16, cex = 0.5, 
     col="black", xlab="Mean F statistic", ylab="Performances (%)",ylim=range(0, 100),  
     xlim=rev(range(0,20)))
lines(IV_0$meanF, 100*IV_0$p2_cover, type="p", pch=16, cex = 0.5, col="black")
lines(IV_0$meanF, 100*IV_0$p3_cover, type="p", pch=16, cex = 0.5, col="black")
lines(IV_0$meanF, 100*IV_0$p1_power, type="p", pch=16, cex = 0.5, col="red")
lines(IV_0$meanF, 100*IV_0$p3_power, type="p", pch=16, cex = 0.5, col="red")

lines(IV_0$meanF, 100*predict(lo_p1_cover), type="l", lty=1,lwd=2, col="black")
lines(IV_0$meanF, 100*predict(lo_p2_cover), type="l", lty=2,lwd=2, col="black")
lines(IV_0$meanF, 100*predict(lo_p3_cover), type="l", lty=3,lwd=2, col="black")
lines(IV_0$meanF, 100*predict(lo_p1_power), type="l", lty=1,lwd=2, col="red")
lines(IV_0$meanF, 100*predict(lo_p3_power), type="l", lty=2,lwd=2, col="red")

legend(20, 70, legend=c("k1 coverage rate", "k2 coverage rate", 
                        "k3 coverage rate", "k1 power", 
                        "k3 power"),
       col=c("black", "black", "black", "red", "red"), 
       lty=c(1,2,3,1,2,1),cex=1, bty="n")

dev.off()



##directional P

IV_0<-read.csv("./UBP.csv", row.names = NULL)
colnames(IV_0)


jpeg('IV_direct_P.png', width = 7, height = 7, units = 'in', res = 600)

lo_p1_cover <- loess(IV_0$p1_cover~IV_0$meanF)
lo_p2_cover <- loess(IV_0$p2_cover~IV_0$meanF)
lo_p3_cover <- loess(IV_0$p3_cover~IV_0$meanF)
lo_p1_power <- loess(IV_0$p1_power~IV_0$meanF)
lo_p3_power <- loess(IV_0$p3_power~IV_0$meanF)
lo_p2_FPR <- loess(IV_0$p2_FPR~IV_0$meanF)

par(mfrow=c(1,1))

plot(main="Directional pleiotropy", IV_0$meanF, 100*IV_0$p1_cover, type="p", pch=16, cex = 0.5, 
     col="black", xlab="Mean F statistic", ylab="Performances (%)",ylim=range(0, 100),  
     xlim=rev(range(0,20)))
lines(IV_0$meanF, 100*IV_0$p2_cover, type="p", pch=16, cex = 0.5, col="black")
lines(IV_0$meanF, 100*IV_0$p3_cover, type="p", pch=16, cex = 0.5, col="black")
lines(IV_0$meanF, 100*IV_0$p1_power, type="p", pch=16, cex = 0.5, col="red")
lines(IV_0$meanF, 100*IV_0$p3_power, type="p", pch=16, cex = 0.5, col="red")

lines(IV_0$meanF, 100*predict(lo_p1_cover), type="l", lty=1,lwd=2, col="black")
lines(IV_0$meanF, 100*predict(lo_p2_cover), type="l", lty=2,lwd=2, col="black")
lines(IV_0$meanF, 100*predict(lo_p3_cover), type="l", lty=3,lwd=2, col="black")
lines(IV_0$meanF, 100*predict(lo_p1_power), type="l", lty=1,lwd=2, col="red")
lines(IV_0$meanF, 100*predict(lo_p3_power), type="l", lty=2,lwd=2, col="red")

legend(20, 70, legend=c("k1 coverage rate", "k2 coverage rate", 
                        "k3 coverage rate", "k1 power", 
                        "k3 power"),
       col=c("black", "black", "black", "red", "red"), 
       lty=c(1,2,3,1,2,1),cex=1, bty="n")

dev.off()

