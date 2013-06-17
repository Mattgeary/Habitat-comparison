library(ncf)
library(mgcv)

data <- read.csv("datahab.csv")
allpa <- read.csv("allpaxy_connect.csv")
allpa[ is.na(allpa) ] <- 0

d1 <- subset(data, subset=(data$p94==1 & data$p00 == 0 |data$p94 == 1 & data$p00==1))
a1 <- subset(allpa, subset=(allpa$p94==1 & allpa$p00 == 0 |allpa$p94 == 1 & allpa$p00==1))
d1$pa[d1$p94 == 1 & d1$p00 == 0] <- 0
d1$pa[d1$p94 == 1 & d1$p00 == 1] <- 1
data.1 <- data.frame(x=d1$x, y=d1$y, "pa" = d1$pa, "Init.count"= d1$Count_94, "Moorland.500" = (d1$h.00.500.3 - d1$h.94.500.3), "Open.500" = (d1$h.00.500.4 - d1$h.94.500.4), "Closed.500" = (d1$h.00.500.5 - d1$h.94.500.5),"sum.con" = ((a1$sum_00+0.1)/(a1$sum_94+0.1)), "count.con"=((a1$birds_00+0.1)/(a1$birds_94+0.1)))
rm(d1)

d3 <- subset(data, subset=(data$p94==1 & data$p08 == 0 |data$p94 == 1 & data$p08==1))
a3 <- subset(allpa, subset=(allpa$p94==1 & allpa$p08 == 0 |allpa$p94 == 1 & allpa$p08==1))
d3$pa[d3$p94 == 1 & d3$p08 == 0] <- 0
d3$pa[d3$p94 == 1 & d3$p08 == 1] <- 1
data.3 <- data.frame(x=d3$x, y=d3$y,"pa" = d3$pa, "Init.count"= d3$Count_94, "Moorland.500" = (d3$h.08.500.3 - d3$h.94.500.3), "Open.500" = (d3$h.08.500.4 - d3$h.94.500.4), "Closed.500" = (d3$h.08.500.5 - d3$h.94.500.5), "sum.con" = ((a3$sum_08+0.1)/(a3$sum_94+0.1)), "count.con"=((a3$birds_08+0.1)/(a3$birds_94+0.1)))
rm(d3)

panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y))
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
    
    test <- cor.test(x,y, method="spearman")
    # borrowed from printCoefmat
    Signif <- symnum(test$p.value, corr = FALSE, na = FALSE,
                  cutpoints = c(0, 0.05,1),
                  symbols = c("*"," "))
    
    text(0.5, 0.5, txt, cex = cex * r)
    text(.8, .8, Signif, cex=cex, col=2)
}

pdf("Hab_corr_connect.pdf", title="Habitat correlations")
pairs(data.1, lower.panel=panel.smooth, upper.panel=panel.cor, main="1994 - 2000")
pairs(data.3, lower.panel=panel.smooth, upper.panel=panel.cor, main="1994 - 2008")
dev.off()

mod.1.smooth <- data.frame(Init.count = numeric(12), Moorland.500 = numeric(12), Open.500=numeric(12), Closed.500=numeric(12),  sum.con=numeric(12), count.con=numeric(12))

for(i in 3:15){
 mod.ic <- gam(pa ~ s(Init.count, k=i), family="binomial", data=data.1)
 mod.m5 <- gam(pa ~ s(Moorland.500, k=i), family="binomial", data=data.1)
 mod.o5 <- gam(pa ~ s(Open.500, k=i), family="binomial", data=data.1)
 mod.c5 <- gam(pa ~ s(Closed.500, k=i), family="binomial", data=data.1)
 mod.sc.5 <- gam(pa ~ s(sum.con, k=i), family="binomial", data=data.1)
 mod.cc.5 <- gam(pa ~ s(count.con, k=i), family="binomial", data=data.1)

 mod.1.smooth$Init.count[i] <- AIC(mod.ic)
 mod.1.smooth$Moorland.500[i] <- AIC(mod.m5)
 mod.1.smooth$Open.500[i] <- AIC(mod.o5)
 mod.1.smooth$Closed.500[i] <- AIC(mod.c5)
 mod.1.smooth$sum.con[i] <- AIC(mod.sc5)
 mod.1.smooth$count.con[i] <- AIC(mod.cc5)
}

mod.3.smooth <- data.frame(Init.count = numeric(12), Moorland.500 = numeric(12), Open.500=numeric(12), Closed.500=numeric(12),  sum.con=numeric(12), count.con=numeric(12))

for(i in 3:15){
 mod.ic <- gam(pa ~ s(Init.count, k=i), family="binomial", data=data.3)
 mod.m5 <- gam(pa ~ s(Moorland.500, k=i), family="binomial", data=data.3)
 mod.o5 <- gam(pa ~ s(Open.500, k=i), family="binomial", data=data.3)
 mod.c5 <- gam(pa ~ s(Closed.500, k=i), family="binomial", data=data.3)
 mod.sc.5 <- gam(pa ~ s(sum.con, k=i), family="binomial", data=data.3)
 mod.cc.5 <- gam(pa ~ s(count.con, k=i), family="binomial", data=data.3)

 mod.3.smooth$Init.count[i] <- AIC(mod.ic)
 mod.3.smooth$Moorland.500[i] <- AIC(mod.m5)
 mod.3.smooth$Open.500[i] <- AIC(mod.o5)
 mod.3.smooth$Closed.500[i] <- AIC(mod.c5)
 mod.3.smooth$sum.con[i] <- AIC(mod.sc5)
 mod.3.smooth$count.con[i] <- AIC(mod.cc5)
}

write.csv(mod.1.smooth, "mod_1_knots_connect.csv")
write.csv(mod.3.smooth, "mod_3_knots_connect.csv")
