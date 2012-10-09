require(reshape)
work <- getwd()

pdf("Chi_square_residuals.pdf")

#### Projected ######

res.chi.5 <- vector("list", 3)
res.chi.2k <- vector("list", 3)
res.cor.5 <- vector("list", 3)
res.cor.2k <- vector("list", 3)

#### 500m ####

#### 2000 - 1994 #####

#### Chi Square ####
setwd(work)
mc.94.00 <- read.csv("MC_9400p.csv")
change.5 <- data.frame("p94" = mc.94.00$p94, "p00" = mc.94.00$p00, "change"= mc.94.00$cc_500) 
change.5$num <- 1 
change.5$num <- 1 
change.5$ID <- 0 
change.5$ID <- NA
#change.5$ID[change.5$p94==1 & change.5$p00 ==1] <- "1 1" 
change.5$ID[change.5$p94==1 & change.5$p00 ==0] <- "1 0" 
change.5$ID[change.5$p94==0 & change.5$p00 ==1] <- "0 1" 
change.5.b <- change.5[!is.na(change.5$ID),]
change.melt <-melt(change.5.b, id=c("change", "p94", "p00", "ID")) 
ch.5 <- cast(change.melt, ID ~ change ~ variable, sum)
ch.5 <- as.data.frame(ch.5)
chi.5 <- chisq.test(ch.5)
res.chi.5[[1]] <- chi.5

#### Chi square Residuals ####

chi.5$residuals 
chi.res <- c(sapply(chi.5$residuals, sum)) 
pos <- c(1,3,2,4) 
vals <- chi.res[pos] 
chi.res.df <- data.frame("cat" = rep(c("0 1", "1 0"), each=2), "gp" = rep(c("-1", "+1"), 2), "vals"=vals)
chi.res.df<-tapply(chi.res.df$vals,list(chi.res.df$gp,chi.res.df$cat),sum)
barplot(chi.res.df,beside=T,col=c("black","grey"), xlab="Group", ylab="Chi-square contribution", main="Chi-square contribution 1994-2000 at 500m", ylim=c(-3, 3)) 
abline(h=1, lty=2, col="red") 
abline(h=-1, lty=2, col="red") 
legend(7, -1.8,rownames(chi.res.df),fill=c("black","grey"))

### Correlation ###

growth <- (mc.94.00$Count_00 - mc.94.00$Count_94)/6
cor.5 <- data.frame("mc.5" = mc.94.00$mc_500, "growth"=growth)
plot(cor.5, xlab="Change in maxent score", ylab="Growth", main="1994-2000 at 500m (Projected)")
c.5 <- cor.test(cor.5$mc.5, cor.5$growth, method="spearman")
res.cor.5[[1]] <- c.5

#### 1994 - 2008 ####

rm(change.5)
rm(change.melt)
rm(ch.5)
rm(chi.5)
rm(cor.5)
rm(c.5)
rm(change.2k)
rm(change.melt)
rm(ch.2k)
rm(chi.2k)
rm(cor.2k)
rm(c.2k)
rm(growth)

#### Chi Square ####
setwd(work)
mc.94.08 <- read.csv("MC_9408p.csv")
change.5 <- data.frame("p94" = mc.94.08$p94, "p08" = mc.94.08$p08, "change"= mc.94.08$cc_500) 
change.5$num <- 1 
change.5$num <- 1 
change.5$ID <- 0 
change.5$ID <- NA
#change.5$ID[change.5$p94==1 & change.5$p00 ==1] <- "1 1" 
change.5$ID[change.5$p94==1 & change.5$p08 ==0] <- "1 0" 
change.5$ID[change.5$p94==0 & change.5$p08 ==1] <- "0 1" 
change.5.b <- change.5[!is.na(change.5$ID),]
change.melt <-melt(change.5.b, id=c("change", "p94", "p08", "ID")) 
ch.5 <- cast(change.melt, ID ~ change ~ variable, sum)
ch.5 <- as.data.frame(ch.5)
chi.5 <- chisq.test(ch.5)
res.chi.5[[2]] <- chi.5

#### Chi square residuals ####

chi.5$residuals 
chi.res <- c(sapply(chi.5$residuals, sum)) 
pos <- c(1,3,2,4) 
vals <- chi.res[pos] 
chi.res.df <- data.frame("cat" = rep(c("0 1", "1 0"), each=2), "gp" = rep(c("-1", "+1"), 2), "vals"=vals)
chi.res.df<-tapply(chi.res.df$vals,list(chi.res.df$gp,chi.res.df$cat),sum)
barplot(chi.res.df,beside=T,col=c("black","grey"), xlab="Group", ylab="Chi-square contribution", main="Chi-square contribution 1994-2008 at 500m", ylim=c(-3, 3)) 
abline(h=1, lty=2, col="red") 
abline(h=-1, lty=2, col="red") 
legend(7, -1.8,rownames(chi.res.df),fill=c("black","grey"))

#### Correlation ####

growth <- (mc.94.08$Count_08 - mc.94.08$Count_94)/14
cor.5 <- data.frame("mc.5" = mc.94.08$mc_500, "growth"=growth)
plot(cor.5, xlab="Change in maxent score", ylab="Growth", main="1994-2008 at 500m (Projected)")
c.5 <- cor.test(cor.5$mc.5, cor.5$growth, method="spearman")
res.cor.5[[2]] <- c.5

#### 2000 - 2008 ####

rm(change.5)
rm(change.melt)
rm(ch.5)
rm(chi.5)
rm(cor.5)
rm(c.5)
rm(change.2k)
rm(change.melt)
rm(ch.2k)
rm(chi.2k)
rm(cor.2k)
rm(c.2k)
rm(growth)

#### Chi-square ####

setwd(work)
mc.00.08 <- read.csv("MC_00p08p.csv")
change.5 <- data.frame("p00" = mc.00.08$p00, "p08" = mc.00.08$p08, "change"= mc.00.08$cc_500) 
change.5$num <- 1 
change.5$num <- 1 
change.5$ID <- 0 
change.5$ID <- NA
#change.5$ID[change.5$p94==1 & change.5$p00 ==1] <- "1 1" 
change.5$ID[change.5$p00==1 & change.5$p08 ==0] <- "1 0" 
change.5$ID[change.5$p00==0 & change.5$p08 ==1] <- "0 1" 
change.5.b <- change.5[!is.na(change.5$ID),]
change.melt <-melt(change.5.b, id=c("change", "p00", "p08", "ID")) 
ch.5 <- cast(change.melt, ID ~ change ~ variable, sum)
ch.5 <- as.data.frame(ch.5)
chi.5 <- chisq.test(ch.5)
res.chi.5[[3]] <- chi.5

##### Chi square residuals #####

chi.5$residuals 
chi.res <- c(sapply(chi.5$residuals, sum)) 
pos <- c(1,3,2,4) 
vals <- chi.res[pos] 
chi.res.df <- data.frame("cat" = rep(c("0 1", "1 0"), each=2), "gp" = rep(c("-1", "+1"), 2), "vals"=vals)
chi.res.df<-tapply(chi.res.df$vals,list(chi.res.df$gp,chi.res.df$cat),sum)
barplot(chi.res.df,beside=T,col=c("black","grey"), xlab="Group", ylab="Chi-square contribution", main="Chi-square contribution 2000-2008 at 500m", ylim=c(-3, 3)) 
abline(h=1, lty=2, col="red") 
abline(h=-1, lty=2, col="red") 
legend(7, -1.8,rownames(chi.res.df),fill=c("black","grey"))

#### Correlation #####

growth <- (mc.00.08$Count_08 - mc.00.08$Count_00)/8
cor.5 <- data.frame("mc.5" = mc.00.08$mc_500, "growth"=growth)
plot(cor.5, xlab="Change in maxent score", ylab="Growth", main="2000-2008 at 500m (Projected)")
c.5 <- cor.test(cor.5$mc.5, cor.5$growth, method="spearman")
res.cor.5[[3]] <- c.5

#### 2km ####

rm(change.5)
rm(change.melt)
rm(ch.5)
rm(chi.5)
rm(cor.5)
rm(c.5)
rm(change.2k)
rm(change.melt)
rm(ch.2k)
rm(chi.2k)
rm(cor.2k)
rm(c.2k)
rm(growth)

#### 1994 - 2000 ####

setwd(work)
mc.94.00 <- read.csv("MC_9400p.csv")
change.2k <- data.frame("p94" = mc.94.00$p94, "p00" = mc.94.00$p00, "change"= mc.94.00$cc_2k) 
change.2k$num <- 1 
change.2k$ID <- NA
#change.2k$ID[change.2k$p94==1 & change.2k$p00 ==1] <- "1 1" 
change.2k$ID[change.2k$p94==1 & change.2k$p00 ==0] <- "1 0" 
change.2k$ID[change.2k$p94==0 & change.2k$p00 ==1] <- "0 1" 
change.2k.b <- change.2k[!is.na(change.2k$ID),]
change.melt <-melt(change.2k.b, id=c("change", "p94", "p00", "ID")) 
ch.2k <- cast(change.melt, ID ~ change ~ variable, sum)
ch.2k <- as.data.frame(ch.2k)
chi.2k <- chisq.test(ch.2k)
res.chi.2k[[1]] <- chi.2k

#### Chi square residuals ####

chi.2k$residuals 
chi.res <- c(sapply(chi.2k$residuals, sum)) 
pos <- c(1,3,2,4) 
vals <- chi.res[pos] 
chi.res.df <- data.frame("cat" = rep(c("0 1", "1 0"), each=2), "gp" = rep(c("-1", "+1"), 2), "vals"=vals)
chi.res.df<-tapply(chi.res.df$vals,list(chi.res.df$gp,chi.res.df$cat),sum)
barplot(chi.res.df,beside=T,col=c("black","grey"), xlab="Group", ylab="Chi-square contribution", main="Chi-square contribution 1994-2000 at 2km", ylim=c(-3, 3)) 
abline(h=1, lty=2, col="red") 
abline(h=-1, lty=2, col="red") 
legend(7, -1.8,rownames(chi.res.df),fill=c("black","grey"))

#### Correlation ####

growth <- (mc.94.00$Count_00 - mc.94.00$Count_94)/6
cor.2k <- data.frame("mc.2k" = mc.94.00$mc_2k, "growth"=growth)
plot(cor.2k, xlab="Change in maxent score", ylab="Growth", main="1994 - 2000 at 2km (Projected)")
c.2k <- cor.test(cor.2k$mc.2k, cor.2k$growth, method="spearman")
res.cor.2k[[1]] <- c.2k

#### 1994 - 2008 ####

rm(change.5)
rm(change.melt)
rm(ch.5)
rm(chi.5)
rm(cor.5)
rm(c.5)
rm(change.2k)
rm(change.melt)
rm(ch.2k)
rm(chi.2k)
rm(cor.2k)
rm(c.2k)
rm(growth)

#### Chi square ####

setwd(work)
mc.94.08 <- read.csv("MC_9408p.csv")
change.2k <- data.frame("p94" = mc.94.08$p94, "p08" = mc.94.08$p08, "change"= mc.94.08$cc_2k) 
change.2k$num <- 1 
change.2k$num <- 1 
change.2k$ID <- NA
#change.2k$ID[change.2k$p94==1 & change.2k$p00 ==1] <- "1 1" 
change.2k$ID[change.2k$p94==1 & change.2k$p08 ==0] <- "1 0" 
change.2k$ID[change.2k$p94==0 & change.2k$p08 ==1] <- "0 1" 
change.2k.b <- change.2k[!is.na(change.2k$ID),]
change.melt <-melt(change.2k.b, id=c("change", "p94", "p08", "ID")) 
ch.2k <- cast(change.melt, ID ~ change ~ variable, sum)
ch.2k <- as.data.frame(ch.2k)
chi.2k <- chisq.test(ch.2k)
res.chi.2k[[2]] <- chi.2k

#### Chi square residuals ####

chi.2k$residuals 
chi.res <- c(sapply(chi.2k$residuals, sum)) 
pos <- c(1,3,2,4) 
vals <- chi.res[pos] 
chi.res.df <- data.frame("cat" = rep(c("0 1", "1 0"), each=2), "gp" = rep(c("-1", "+1"), 2), "vals"=vals)
chi.res.df<-tapply(chi.res.df$vals,list(chi.res.df$gp,chi.res.df$cat),sum)
barplot(chi.res.df,beside=T,col=c("black","grey"), xlab="Group", ylab="Chi-square contribution", main="Chi-square contribution 1994-2008 at 2km", ylim=c(-3, 3)) 
abline(h=1, lty=2, col="red") 
abline(h=-1, lty=2, col="red") 
legend(7, -1.8,rownames(chi.res.df),fill=c("black","grey"))

#### Correlation ####

growth <- (mc.94.08$Count_08 - mc.94.08$Count_94)/14
cor.2k <- data.frame("mc.2k" = mc.94.08$mc_2k, "growth"=growth)
plot(cor.2k, xlab="Change in maxent score", ylab="Growth", main="1994 - 2008 at 2km (Projected)")
c.2k <- cor.test(cor.2k$mc.2k, cor.2k$growth, method="spearman")
res.cor.2k[[2]] <- c.2k

#### 2000 - 2008 ####

rm(change.5)
rm(change.melt)
rm(ch.5)
rm(chi.5)
rm(cor.5)
rm(c.5)
rm(change.2k)
rm(change.melt)
rm(ch.2k)
rm(chi.2k)
rm(cor.2k)
rm(c.2k)
rm(growth)

#### Chi square ####

setwd(work)
mc.00.08 <- read.csv("MC_00p08p.csv")
change.2k <- data.frame("p00" = mc.00.08$p00, "p08" = mc.00.08$p08, "change"= mc.00.08$cc_2k) 
change.2k$num <- 1 
change.2k$num <- 1 
change.2k$ID <- NA
#change.2k$ID[change.2k$p00==1 & change.2k$p08 ==1] <- "1 1" 
change.2k$ID[change.2k$p00==1 & change.2k$p08 ==0] <- "1 0" 
change.2k$ID[change.2k$p00==0 & change.2k$p08 ==1] <- "0 1" 
change.2k.b <- change.2k[!is.na(change.2k$ID),]
change.melt <-melt(change.2k, id=c("change", "p00", "p08", "ID")) 
ch.2k <- cast(change.melt, ID ~ change ~ variable, sum)
ch.2k <- as.data.frame(ch.2k)
chi.2k <- chisq.test(ch.2k)
res.chi.2k[[3]] <- chi.2k

#### Chi square residuals ####

chi.2k$residuals 
chi.res <- c(sapply(chi.2k$residuals, sum)) 
pos <- c(1,3,2,4) 
vals <- chi.res[pos] 
chi.res.df <- data.frame("cat" = rep(c("0 1", "1 0"), each=2), "gp" = rep(c("-1", "+1"), 2), "vals"=vals)
chi.res.df<-tapply(chi.res.df$vals,list(chi.res.df$gp,chi.res.df$cat),sum)
barplot(chi.res.df,beside=T,col=c("black","grey"), xlab="Group", ylab="Chi-square contribution", main="Chi-square contribution 2000-2008 at 2km", ylim=c(-3, 3))  
abline(h=1, lty=2, col="red") 
abline(h=-1, lty=2, col="red") 
legend(7, -1.8,rownames(chi.res.df),fill=c("black","grey"))

#### Correlation ####

growth <- (mc.00.08$Count_08 - mc.00.08$Count_00)/8
cor.2k <- data.frame("mc.2k" = mc.00.08$mc_2k, "growth"=growth)
plot(cor.2k, xlab="Change in maxent score", ylab="Growth", main="2000 - 2008 at 2km")
c.2k <- cor.test(cor.2k$mc.2k, cor.2k$growth, method="spearman")
res.cor.2k[[3]] <- c.2k

dev.off()
