library(ncf)
library(mgcv)

#### read in habitat and connectivity data for each lek ####
data <- read.csv("datahab.csv")
allpa <- read.csv("allpaxy_connect_new_OCT_2012.csv")
allpa[ is.na(allpa) ] <- 0
################################

#### Prepare datasets###############
####### 1994 - 2000 ################
d1 <- subset(data, subset=(data$p94==1 & data$p00 == 0 |data$p94 == 1 & data$p00==1))
a1 <- subset(allpa, subset=(allpa$p94==1 & allpa$p00 == 0 |allpa$p94 == 1 & allpa$p00==1))
d1$pa[d1$p94 == 1 & d1$p00 == 0] <- 0
d1$pa[d1$p94 == 1 & d1$p00 == 1] <- 1
data.1 <- data.frame(x=d1$x, y=d1$y, "pa" = d1$pa, "Init.count"= d1$Count_94, "Moorland.500" = (d1$h.00.500.3 - d1$h.94.500.3), "Open.500" = (d1$h.00.500.4 - d1$h.94.500.4), "Closed.500" = (d1$h.00.500.5 - d1$h.94.500.5),"birds_94" = a1$birds_94, "leks_94" = a1$leks_94)
rm(d1)

########## 2000 - 2008 ############
d3 <- subset(data, subset=(data$p00==1 & data$p08 == 0 |data$p00 == 1 & data$p08==1))
a3 <- subset(allpa, subset=(allpa$p00==1 & allpa$p08 == 0 |allpa$p00 == 1 & allpa$p08==1))
d3$pa[d3$p00 == 1 & d3$p08 == 0] <- 0
d3$pa[d3$p00 == 1 & d3$p08 == 1] <- 1
data.3 <- data.frame(x=d3$x, y=d3$y,"pa" = d3$pa, "Init.count"= d3$Count_00, "Moorland.500" = (d3$h.08.500.3 - d3$h.00.500.3), "Open.500" = (d3$h.08.500.4 - d3$h.00.500.4), "Closed.500" = (d3$h.08.500.5 - d3$h.00.500.5), "birds_00" = a1$birds_00, "leks_00" = a1$leks_00)
rm(d3)

########## 1994 - 2008 ############
d3 <- subset(data, subset=(data$p94==1 & data$p08 == 0 |data$p94 == 1 & data$p08==1))
a3 <- subset(allpa, subset=(allpa$p94==1 & allpa$p08 == 0 |allpa$p94 == 1 & allpa$p08==1))
d3$pa[d3$p94 == 1 & d3$p08 == 0] <- 0
d3$pa[d3$p94 == 1 & d3$p08 == 1] <- 1
data.3 <- data.frame(x=d3$x, y=d3$y,"pa" = d3$pa, "Init.count"= d3$Count_94, "Moorland.500" = (d3$h.08.500.3 - d3$h.94.500.3), "Open.500" = (d3$h.08.500.4 - d3$h.94.500.4), "Closed.500" = (d3$h.08.500.5 - d3$h.94.500.5), "birds_94" = a1$birds_94, "leks_94" = a1$leks_94)
rm(d3)

models.1 <- data.frame(model=numeric(7), AIC=numeric(7))
models.3 <- data.frame(model=numeric(15), AIC=numeric(15))
#############################

#### Test for correlations between predictors ####

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

pdf("Hab_corr_OCT_2012.pdf", title="Habitat correlations")
pairs(hab.00, lower.panel=panel.smooth, upper.panel=panel.cor, main="1994 - 2000")
pairs(hab.08, lower.panel=panel.smooth, upper.panel=panel.cor, main="2000 - 2008")
pairs(hab.long, lower.panel=panel.smooth, upper.panel=panel.cor, main="1994 - 2008")
dev.off()
###### Run GAMs #############


#### Declare variables ####

for (i in 1:63){
	models.1[i,1] <- paste("gam.1",i, sep=".")
	}
AICs.1 <- rep(NA, 7)

models.1 <- data.frame(model = models.1, AIC = AICs.1)
models.1.list <- vector("list", 63)

for (i in 1:63){
	models.2[i,1] <- paste("gam.2",i, sep=".")
	}
AICs.2 <- rep(NA, 15)

models.2 <- data.frame(model = models.2, AIC = AICs.2)
models.2.list <- vector("list", 63)


for (i in 1:63){
	models.3[i,1] <- paste("gam.3",i, sep=".")
	}
AICs.3 <- rep(NA, 15)

models.3 <- data.frame(model = models.3, AIC = AICs.3)
models.3.list <- vector("list", 63)

#### 1994 - 2000 ##### 
x <- 1
models.1.list[[x]] <- gam(pa ~  +s(Init.counts), data=data.1)
models.1$AIC[x] <- AIC(models.1.list[[x]]
x <- x + 1 
models.1.list[[x]] <- gam(pa ~  +s(Open.500s), data=data.1)
models.1$AIC[x] <- AIC(models.1.list[[x]]
x <- x + 1 
models.1.list[[x]] <- gam(pa ~  +s(Closed.500s), data=data.1)
models.1$AIC[x] <- AIC(models.1.list[[x]]
x <- x + 1 
models.1.list[[x]] <- gam(pa ~  +s(Moorland.500s), data=data.1)
models.1$AIC[x] <- AIC(models.1.list[[x]]
x <- x + 1 
models.1.list[[x]] <- gam(pa ~  +s(males_94), data=data.1)
models.1$AIC[x] <- AIC(models.1.list[[x]]
x <- x + 1 
models.1.list[[x]] <- gam(pa ~  +s(leks_94), data=data.1)
models.1$AIC[x] <- AIC(models.1.list[[x]]
x <- x + 1 
models.1.list[[x]] <- gam(pa ~  +s(Init.count) + s(Open.500) , data=data.1)
models.1$AIC[x] <- AIC(models.1.list[[x]]
x <- x + 1 
models.1.list[[x]] <- gam(pa ~  +s(Init.count) + s(Closed.500s(, data=data.1)
models.1$AIC[x] <- AIC(models.1.list[[x]]
x <- x + 1 
models.1.list[[x]] <- gam(pa ~  +s(Init.count) + s(Moorland.500s(, data=data.1)
models.1$AIC[x] <- AIC(models.1.list[[x]]
x <- x + 1 
models.1.list[[x]] <- gam(pa ~  +s(Init.count) + s(males_94), data=data.1)
models.1$AIC[x] <- AIC(models.1.list[[x]]
x <- x + 1 
models.1.list[[x]] <- gam(pa ~  +s(Init.count) + s(leks_94), data=data.1)
models.1$AIC[x] <- AIC(models.1.list[[x]]
x <- x + 1 
models.1.list[[x]] <- gam(pa ~  +s(Open.500) + s(Closed.500s(, data=data.1)
models.1$AIC[x] <- AIC(models.1.list[[x]]
x <- x + 1 
models.1.list[[x]] <- gam(pa ~  +s(Open.500) + s(Moorland.500s(, data=data.1)
models.1$AIC[x] <- AIC(models.1.list[[x]]
x <- x + 1 
models.1.list[[x]] <- gam(pa ~  +s(Open.500) + s(males_94)  , data=data.1)
models.1$AIC[x] <- AIC(models.1.list[[x]]
x <- x + 1 
models.1.list[[x]] <- gam(pa ~  +s(Open.500) + s(leks_94)   , data=data.1)
models.1$AIC[x] <- AIC(models.1.list[[x]]
x <- x + 1 
models.1.list[[x]] <- gam(pa ~  +s(Closed.500) + s(Moorland.500s(, data=data.1)
models.1$AIC[x] <- AIC(models.1.list[[x]]
x <- x + 1 
models.1.list[[x]] <- gam(pa ~  +s(Closed.500) + s(males_94), data=data.1)
models.1$AIC[x] <- AIC(models.1.list[[x]]
x <- x + 1 
models.1.list[[x]] <- gam(pa ~  +s(Closed.500) + s(leks_94), data=data.1)
models.1$AIC[x] <- AIC(models.1.list[[x]]
x <- x + 1 
models.1.list[[x]] <- gam(pa ~  +s(Moorland.500) + s(males_94), data=data.1)
models.1$AIC[x] <- AIC(models.1.list[[x]]
x <- x + 1 
models.1.list[[x]] <- gam(pa ~  +s(Moorland.500) + s(leks_94), data=data.1)
models.1$AIC[x] <- AIC(models.1.list[[x]]
x <- x + 1 
models.1.list[[x]] <- gam(pa ~  +s(males_94) + s(leks_94), data=data.1)
models.1$AIC[x] <- AIC(models.1.list[[x]]
x <- x + 1 
models.1.list[[x]] <- gam(pa ~  +s(Init.count) + s(Open.500) + s(Closed.500s(, data=data.1)
models.1$AIC[x] <- AIC(models.1.list[[x]]
x <- x + 1 
models.1.list[[x]] <- gam(pa ~  +s(Init.count) + s(Open.500) + s(Moorland.500s(, data=data.1)
models.1$AIC[x] <- AIC(models.1.list[[x]]
x <- x + 1 
models.1.list[[x]] <- gam(pa ~  +s(Init.count) + s(Open.500) + s(males_94), data=data.1)
models.1$AIC[x] <- AIC(models.1.list[[x]]
x <- x + 1 
models.1.list[[x]] <- gam(pa ~  +s(Init.count) + s(Open.500) + s(leks_94), data=data.1)
models.1$AIC[x] <- AIC(models.1.list[[x]]
x <- x + 1 
models.1.list[[x]] <- gam(pa ~  +s(Init.count) + s(Closed.500) + s(Moorland.500s(, data=data.1)
models.1$AIC[x] <- AIC(models.1.list[[x]]
x <- x + 1 
models.1.list[[x]] <- gam(pa ~  +s(Init.count) + s(Closed.500) + s(males_94), data=data.1)
models.1$AIC[x] <- AIC(models.1.list[[x]]
x <- x + 1 
models.1.list[[x]] <- gam(pa ~  +s(Init.count) + s(Closed.500) + s(leks_94), data=data.1)
models.1$AIC[x] <- AIC(models.1.list[[x]]
x <- x + 1 
models.1.list[[x]] <- gam(pa ~  +s(Init.count) + s(Moorland.500) + s(males_94), data=data.1)
models.1$AIC[x] <- AIC(models.1.list[[x]]
x <- x + 1 
models.1.list[[x]] <- gam(pa ~  +s(Init.count) + s(Moorland.500) + s(leks_94), data=data.1)
models.1$AIC[x] <- AIC(models.1.list[[x]]
x <- x + 1 
models.1.list[[x]] <- gam(pa ~  +s(Init.count) + s(males_94) + s(leks_94), data=data.1)
models.1$AIC[x] <- AIC(models.1.list[[x]]
x <- x + 1 
models.1.list[[x]] <- gam(pa ~  +s(Open.500) + s(Closed.500) + s(Moorland.500s(, data=data.1)
models.1$AIC[x] <- AIC(models.1.list[[x]]
x <- x + 1 
models.1.list[[x]] <- gam(pa ~  +s(Open.500) + s(Closed.500) + s(males_94), data=data.1)
models.1$AIC[x] <- AIC(models.1.list[[x]]
x <- x + 1 
models.1.list[[x]] <- gam(pa ~  +s(Open.500) + s(Closed.500) + s(leks_94), data=data.1)
models.1$AIC[x] <- AIC(models.1.list[[x]]
x <- x + 1 
models.1.list[[x]] <- gam(pa ~  +s(Open.500) + s(Moorland.500) + s(males_94), data=data.1)
models.1$AIC[x] <- AIC(models.1.list[[x]]
x <- x + 1 
models.1.list[[x]] <- gam(pa ~  +s(Open.500) + s(Moorland.500) + s(leks_94), data=data.1)
models.1$AIC[x] <- AIC(models.1.list[[x]]
x <- x + 1 
models.1.list[[x]] <- gam(pa ~  +s(Open.500) + s(males_94)   s(leks_94)   , data=data.1)
models.1$AIC[x] <- AIC(models.1.list[[x]]
x <- x + 1 
models.1.list[[x]] <- gam(pa ~  +s(Closed.500) + s(Moorland.500) + s(males_94), data=data.1)
models.1$AIC[x] <- AIC(models.1.list[[x]]
x <- x + 1 
models.1.list[[x]] <- gam(pa ~  +s(Closed.500) + s(Moorland.500) + s(leks_94), data=data.1)
models.1$AIC[x] <- AIC(models.1.list[[x]]
x <- x + 1 
models.1.list[[x]] <- gam(pa ~  +s(Closed.500) + s(males_94) + s(leks_94), data=data.1)
models.1$AIC[x] <- AIC(models.1.list[[x]]
x <- x + 1 
models.1.list[[x]] <- gam(pa ~  +s(Moorland.500) + s(males_94) + s(leks_94), data=data.1)
models.1$AIC[x] <- AIC(models.1.list[[x]]
x <- x + 1 
models.1.list[[x]] <- gam(pa ~  +s(Init.count) + s(Open.500) + s(Closed.500) + s(Moorland.500s(, data=data.1)
models.1$AIC[x] <- AIC(models.1.list[[x]]
x <- x + 1 
models.1.list[[x]] <- gam(pa ~  +s(Init.count) + s(Open.500) + s(Closed.500) + s(males_94), data=data.1)
models.1$AIC[x] <- AIC(models.1.list[[x]]
x <- x + 1 
models.1.list[[x]] <- gam(pa ~  +s(Init.count) + s(Open.500) + s(Closed.500) + s(leks_94), data=data.1)
models.1$AIC[x] <- AIC(models.1.list[[x]]
x <- x + 1 
models.1.list[[x]] <- gam(pa ~  +s(Init.count) + s(Open.500) + s(Moorland.500) + s(males_94), data=data.1)
models.1$AIC[x] <- AIC(models.1.list[[x]]
x <- x + 1 
models.1.list[[x]] <- gam(pa ~  +s(Init.count) + s(Open.500) + s(Moorland.500) + s(leks_94), data=data.1)
models.1$AIC[x] <- AIC(models.1.list[[x]]
x <- x + 1 
models.1.list[[x]] <- gam(pa ~  +s(Init.count) + s(Open.500) + s(males_94) + s(leks_94), data=data.1)
models.1$AIC[x] <- AIC(models.1.list[[x]]
x <- x + 1 
models.1.list[[x]] <- gam(pa ~  +s(Init.count) + s(Closed.500) + s(Moorland.500) + s(males_94), data=data.1)
models.1$AIC[x] <- AIC(models.1.list[[x]]
x <- x + 1 
models.1.list[[x]] <- gam(pa ~  +s(Init.count) + s(Closed.500) + s(Moorland.500) + s(leks_94), data=data.1)
models.1$AIC[x] <- AIC(models.1.list[[x]]
x <- x + 1 
models.1.list[[x]] <- gam(pa ~  +s(Init.count) + s(Closed.500) + s(males_94) + s(leks_94), data=data.1)
models.1$AIC[x] <- AIC(models.1.list[[x]]
x <- x + 1 
models.1.list[[x]] <- gam(pa ~  +s(Init.count) + s(Moorland.500) + s(males_94) + s(leks_94), data=data.1)
models.1$AIC[x] <- AIC(models.1.list[[x]]
x <- x + 1 
models.1.list[[x]] <- gam(pa ~  +s(Open.500) + s(Closed.500) + s(Moorland.500) + s(males_94), data=data.1)
models.1$AIC[x] <- AIC(models.1.list[[x]]
x <- x + 1 
models.1.list[[x]] <- gam(pa ~  +s(Open.500) + s(Closed.500) + s(Moorland.500) + s(leks_94), data=data.1)
models.1$AIC[x] <- AIC(models.1.list[[x]]
x <- x + 1 
models.1.list[[x]] <- gam(pa ~  +s(Open.500) + s(Closed.500) + s(males_94) + s(leks_94), data=data.1)
models.1$AIC[x] <- AIC(models.1.list[[x]]
x <- x + 1 
models.1.list[[x]] <- gam(pa ~  +s(Open.500) + s(Moorland.500) + s(males_94) + s(leks_94), data=data.1)
models.1$AIC[x] <- AIC(models.1.list[[x]]
x <- x + 1 
models.1.list[[x]] <- gam(pa ~  +s(Closed.500) + s(Moorland.500) + s(males_94) + s(leks_94), data=data.1)
models.1$AIC[x] <- AIC(models.1.list[[x]]
x <- x + 1 
models.1.list[[x]] <- gam(pa ~  +s(Init.count) + s(Open.500) + s(Closed.500) + s(Moorland.500) + s(males_94), data=data.1)
models.1$AIC[x] <- AIC(models.1.list[[x]]
x <- x + 1 
models.1.list[[x]] <- gam(pa ~  +s(Init.count) + s(Open.500) + s(Closed.500) + s(Moorland.500) + s(leks_94), data=data.1)
models.1$AIC[x] <- AIC(models.1.list[[x]]
x <- x + 1 
models.1.list[[x]] <- gam(pa ~  +s(Init.count) + s(Open.500) + s(Closed.500) + s(males_94) + s(leks_94), data=data.1)
models.1$AIC[x] <- AIC(models.1.list[[x]]
x <- x + 1 
models.1.list[[x]] <- gam(pa ~  +s(Init.count) + s(Open.500) + s(Moorland.500) + s(males_94) + s(leks_94), data=data.1)
models.1$AIC[x] <- AIC(models.1.list[[x]]
x <- x + 1 
models.1.list[[x]] <- gam(pa ~  +s(Init.count) + s(Closed.500) + s(Moorland.500) + s(males_94) + s(leks_94), data=data.1)
models.1$AIC[x] <- AIC(models.1.list[[x]]
x <- x + 1 
models.1.list[[x]] <- gam(pa ~  +s(Open.500) + s(Closed.500) + s(Moorland.500) + s(males_94) + s(leks_94), data=data.1)
models.1$AIC[x] <- AIC(models.1.list[[x]]
x <- x + 1 
models.1.list[[x]] <- gam(pa ~  + s(Init.count) + s(Open.500) + s(Closed.500) + s(Moorland.500) + s(males_94) + s(leks_94), data=data.1)
 x <- x +1

######## 2000 - 2008 ########
x <- 1
models.2.list[[x]] <- gam(pa ~  +s(Init.counts), data=data.2)
models.2$AIC[x] <- AIC(models.2.list[[x]]
x <- x + 1 
models.2.list[[x]] <- gam(pa ~  +s(Open.500s), data=data.2)
models.2$AIC[x] <- AIC(models.2.list[[x]]
x <- x + 1 
models.2.list[[x]] <- gam(pa ~  +s(Closed.500s), data=data.2)
models.2$AIC[x] <- AIC(models.2.list[[x]]
x <- x + 1 
models.2.list[[x]] <- gam(pa ~  +s(Moorland.500s), data=data.2)
models.2$AIC[x] <- AIC(models.2.list[[x]]
x <- x + 1 
models.2.list[[x]] <- gam(pa ~  +s(males_94), data=data.2)
models.2$AIC[x] <- AIC(models.2.list[[x]]
x <- x + 1 
models.2.list[[x]] <- gam(pa ~  +s(leks_00), data=data.2)
models.2$AIC[x] <- AIC(models.2.list[[x]]
x <- x + 1 
models.2.list[[x]] <- gam(pa ~  +s(Init.count) + s(Open.500) , data=data.2)
models.2$AIC[x] <- AIC(models.2.list[[x]]
x <- x + 1 
models.2.list[[x]] <- gam(pa ~  +s(Init.count) + s(Closed.500s(, data=data.2)
models.2$AIC[x] <- AIC(models.2.list[[x]]
x <- x + 1 
models.2.list[[x]] <- gam(pa ~  +s(Init.count) + s(Moorland.500s(, data=data.2)
models.2$AIC[x] <- AIC(models.2.list[[x]]
x <- x + 1 
models.2.list[[x]] <- gam(pa ~  +s(Init.count) + s(males_94), data=data.2)
models.2$AIC[x] <- AIC(models.2.list[[x]]
x <- x + 1 
models.2.list[[x]] <- gam(pa ~  +s(Init.count) + s(leks_00), data=data.2)
models.2$AIC[x] <- AIC(models.2.list[[x]]
x <- x + 1 
models.2.list[[x]] <- gam(pa ~  +s(Open.500) + s(Closed.500s(, data=data.2)
models.2$AIC[x] <- AIC(models.2.list[[x]]
x <- x + 1 
models.2.list[[x]] <- gam(pa ~  +s(Open.500) + s(Moorland.500s(, data=data.2)
models.2$AIC[x] <- AIC(models.2.list[[x]]
x <- x + 1 
models.2.list[[x]] <- gam(pa ~  +s(Open.500) + s(males_94)  , data=data.2)
models.2$AIC[x] <- AIC(models.2.list[[x]]
x <- x + 1 
models.2.list[[x]] <- gam(pa ~  +s(Open.500) + s(leks_00)   , data=data.2)
models.2$AIC[x] <- AIC(models.2.list[[x]]
x <- x + 1 
models.2.list[[x]] <- gam(pa ~  +s(Closed.500) + s(Moorland.500s(, data=data.2)
models.2$AIC[x] <- AIC(models.2.list[[x]]
x <- x + 1 
models.2.list[[x]] <- gam(pa ~  +s(Closed.500) + s(males_94), data=data.2)
models.2$AIC[x] <- AIC(models.2.list[[x]]
x <- x + 1 
models.2.list[[x]] <- gam(pa ~  +s(Closed.500) + s(leks_00), data=data.2)
models.2$AIC[x] <- AIC(models.2.list[[x]]
x <- x + 1 
models.2.list[[x]] <- gam(pa ~  +s(Moorland.500) + s(males_94), data=data.2)
models.2$AIC[x] <- AIC(models.2.list[[x]]
x <- x + 1 
models.2.list[[x]] <- gam(pa ~  +s(Moorland.500) + s(leks_00), data=data.2)
models.2$AIC[x] <- AIC(models.2.list[[x]]
x <- x + 1 
models.2.list[[x]] <- gam(pa ~  +s(males_94) + s(leks_00), data=data.2)
models.2$AIC[x] <- AIC(models.2.list[[x]]
x <- x + 1 
models.2.list[[x]] <- gam(pa ~  +s(Init.count) + s(Open.500) + s(Closed.500s(, data=data.2)
models.2$AIC[x] <- AIC(models.2.list[[x]]
x <- x + 1 
models.2.list[[x]] <- gam(pa ~  +s(Init.count) + s(Open.500) + s(Moorland.500s(, data=data.2)
models.2$AIC[x] <- AIC(models.2.list[[x]]
x <- x + 1 
models.2.list[[x]] <- gam(pa ~  +s(Init.count) + s(Open.500) + s(males_94), data=data.2)
models.2$AIC[x] <- AIC(models.2.list[[x]]
x <- x + 1 
models.2.list[[x]] <- gam(pa ~  +s(Init.count) + s(Open.500) + s(leks_00), data=data.2)
models.2$AIC[x] <- AIC(models.2.list[[x]]
x <- x + 1 
models.2.list[[x]] <- gam(pa ~  +s(Init.count) + s(Closed.500) + s(Moorland.500s(, data=data.2)
models.2$AIC[x] <- AIC(models.2.list[[x]]
x <- x + 1 
models.2.list[[x]] <- gam(pa ~  +s(Init.count) + s(Closed.500) + s(males_94), data=data.2)
models.2$AIC[x] <- AIC(models.2.list[[x]]
x <- x + 1 
models.2.list[[x]] <- gam(pa ~  +s(Init.count) + s(Closed.500) + s(leks_00), data=data.2)
models.2$AIC[x] <- AIC(models.2.list[[x]]
x <- x + 1 
models.2.list[[x]] <- gam(pa ~  +s(Init.count) + s(Moorland.500) + s(males_94), data=data.2)
models.2$AIC[x] <- AIC(models.2.list[[x]]
x <- x + 1 
models.2.list[[x]] <- gam(pa ~  +s(Init.count) + s(Moorland.500) + s(leks_00), data=data.2)
models.2$AIC[x] <- AIC(models.2.list[[x]]
x <- x + 1 
models.2.list[[x]] <- gam(pa ~  +s(Init.count) + s(males_94) + s(leks_00), data=data.2)
models.2$AIC[x] <- AIC(models.2.list[[x]]
x <- x + 1 
models.2.list[[x]] <- gam(pa ~  +s(Open.500) + s(Closed.500) + s(Moorland.500s(, data=data.2)
models.2$AIC[x] <- AIC(models.2.list[[x]]
x <- x + 1 
models.2.list[[x]] <- gam(pa ~  +s(Open.500) + s(Closed.500) + s(males_94), data=data.2)
models.2$AIC[x] <- AIC(models.2.list[[x]]
x <- x + 1 
models.2.list[[x]] <- gam(pa ~  +s(Open.500) + s(Closed.500) + s(leks_00), data=data.2)
models.2$AIC[x] <- AIC(models.2.list[[x]]
x <- x + 1 
models.2.list[[x]] <- gam(pa ~  +s(Open.500) + s(Moorland.500) + s(males_94), data=data.2)
models.2$AIC[x] <- AIC(models.2.list[[x]]
x <- x + 1 
models.2.list[[x]] <- gam(pa ~  +s(Open.500) + s(Moorland.500) + s(leks_00), data=data.2)
models.2$AIC[x] <- AIC(models.2.list[[x]]
x <- x + 1 
models.2.list[[x]] <- gam(pa ~  +s(Open.500) + s(males_94)   s(leks_00)   , data=data.2)
models.2$AIC[x] <- AIC(models.2.list[[x]]
x <- x + 1 
models.2.list[[x]] <- gam(pa ~  +s(Closed.500) + s(Moorland.500) + s(males_94), data=data.2)
models.2$AIC[x] <- AIC(models.2.list[[x]]
x <- x + 1 
models.2.list[[x]] <- gam(pa ~  +s(Closed.500) + s(Moorland.500) + s(leks_00), data=data.2)
models.2$AIC[x] <- AIC(models.2.list[[x]]
x <- x + 1 
models.2.list[[x]] <- gam(pa ~  +s(Closed.500) + s(males_94) + s(leks_00), data=data.2)
models.2$AIC[x] <- AIC(models.2.list[[x]]
x <- x + 1 
models.2.list[[x]] <- gam(pa ~  +s(Moorland.500) + s(males_94) + s(leks_00), data=data.2)
models.2$AIC[x] <- AIC(models.2.list[[x]]
x <- x + 1 
models.2.list[[x]] <- gam(pa ~  +s(Init.count) + s(Open.500) + s(Closed.500) + s(Moorland.500s(, data=data.2)
models.2$AIC[x] <- AIC(models.2.list[[x]]
x <- x + 1 
models.2.list[[x]] <- gam(pa ~  +s(Init.count) + s(Open.500) + s(Closed.500) + s(males_94), data=data.2)
models.2$AIC[x] <- AIC(models.2.list[[x]]
x <- x + 1 
models.2.list[[x]] <- gam(pa ~  +s(Init.count) + s(Open.500) + s(Closed.500) + s(leks_00), data=data.2)
models.2$AIC[x] <- AIC(models.2.list[[x]]
x <- x + 1 
models.2.list[[x]] <- gam(pa ~  +s(Init.count) + s(Open.500) + s(Moorland.500) + s(males_94), data=data.2)
models.2$AIC[x] <- AIC(models.2.list[[x]]
x <- x + 1 
models.2.list[[x]] <- gam(pa ~  +s(Init.count) + s(Open.500) + s(Moorland.500) + s(leks_00), data=data.2)
models.2$AIC[x] <- AIC(models.2.list[[x]]
x <- x + 1 
models.2.list[[x]] <- gam(pa ~  +s(Init.count) + s(Open.500) + s(males_94) + s(leks_00), data=data.2)
models.2$AIC[x] <- AIC(models.2.list[[x]]
x <- x + 1 
models.2.list[[x]] <- gam(pa ~  +s(Init.count) + s(Closed.500) + s(Moorland.500) + s(males_94), data=data.2)
models.2$AIC[x] <- AIC(models.2.list[[x]]
x <- x + 1 
models.2.list[[x]] <- gam(pa ~  +s(Init.count) + s(Closed.500) + s(Moorland.500) + s(leks_00), data=data.2)
models.2$AIC[x] <- AIC(models.2.list[[x]]
x <- x + 1 
models.2.list[[x]] <- gam(pa ~  +s(Init.count) + s(Closed.500) + s(males_94) + s(leks_00), data=data.2)
models.2$AIC[x] <- AIC(models.2.list[[x]]
x <- x + 1 
models.2.list[[x]] <- gam(pa ~  +s(Init.count) + s(Moorland.500) + s(males_94) + s(leks_00), data=data.2)
models.2$AIC[x] <- AIC(models.2.list[[x]]
x <- x + 1 
models.2.list[[x]] <- gam(pa ~  +s(Open.500) + s(Closed.500) + s(Moorland.500) + s(males_94), data=data.2)
models.2$AIC[x] <- AIC(models.2.list[[x]]
x <- x + 1 
models.2.list[[x]] <- gam(pa ~  +s(Open.500) + s(Closed.500) + s(Moorland.500) + s(leks_00), data=data.2)
models.2$AIC[x] <- AIC(models.2.list[[x]]
x <- x + 1 
models.2.list[[x]] <- gam(pa ~  +s(Open.500) + s(Closed.500) + s(males_94) + s(leks_00), data=data.2)
models.2$AIC[x] <- AIC(models.2.list[[x]]
x <- x + 1 
models.2.list[[x]] <- gam(pa ~  +s(Open.500) + s(Moorland.500) + s(males_94) + s(leks_00), data=data.2)
models.2$AIC[x] <- AIC(models.2.list[[x]]
x <- x + 1 
models.2.list[[x]] <- gam(pa ~  +s(Closed.500) + s(Moorland.500) + s(males_94) + s(leks_00), data=data.2)
models.2$AIC[x] <- AIC(models.2.list[[x]]
x <- x + 1 
models.2.list[[x]] <- gam(pa ~  +s(Init.count) + s(Open.500) + s(Closed.500) + s(Moorland.500) + s(males_94), data=data.2)
models.2$AIC[x] <- AIC(models.2.list[[x]]
x <- x + 1 
models.2.list[[x]] <- gam(pa ~  +s(Init.count) + s(Open.500) + s(Closed.500) + s(Moorland.500) + s(leks_00), data=data.2)
models.2$AIC[x] <- AIC(models.2.list[[x]]
x <- x + 1 
models.2.list[[x]] <- gam(pa ~  +s(Init.count) + s(Open.500) + s(Closed.500) + s(males_94) + s(leks_00), data=data.2)
models.2$AIC[x] <- AIC(models.2.list[[x]]
x <- x + 1 
models.2.list[[x]] <- gam(pa ~  +s(Init.count) + s(Open.500) + s(Moorland.500) + s(males_94) + s(leks_00), data=data.2)
models.2$AIC[x] <- AIC(models.2.list[[x]]
x <- x + 1 
models.2.list[[x]] <- gam(pa ~  +s(Init.count) + s(Closed.500) + s(Moorland.500) + s(males_94) + s(leks_00), data=data.2)
models.2$AIC[x] <- AIC(models.2.list[[x]]
x <- x + 1 
models.2.list[[x]] <- gam(pa ~  +s(Open.500) + s(Closed.500) + s(Moorland.500) + s(males_94) + s(leks_00), data=data.2)
models.2$AIC[x] <- AIC(models.2.list[[x]]
x <- x + 1 
models.2.list[[x]] <- gam(pa ~  + s(Init.count) + s(Open.500) + s(Closed.500) + s(Moorland.500) + s(males_94) + s(leks_00), data=data.2)
 x <- x +1


######## 1994- 2008 ########
x <- 1
models.3.list[[x]] <- gam(pa ~  +s(Init.counts), data=data.3)
models.3$AIC[x] <- AIC(models.3.list[[x]]
x <- x + 1 
models.3.list[[x]] <- gam(pa ~  +s(Open.500s), data=data.3)
models.3$AIC[x] <- AIC(models.3.list[[x]]
x <- x + 1 
models.3.list[[x]] <- gam(pa ~  +s(Closed.500s), data=data.3)
models.3$AIC[x] <- AIC(models.3.list[[x]]
x <- x + 1 
models.3.list[[x]] <- gam(pa ~  +s(Moorland.500s), data=data.3)
models.3$AIC[x] <- AIC(models.3.list[[x]]
x <- x + 1 
models.3.list[[x]] <- gam(pa ~  +s(males_94), data=data.3)
models.3$AIC[x] <- AIC(models.3.list[[x]]
x <- x + 1 
models.3.list[[x]] <- gam(pa ~  +s(leks_94), data=data.3)
models.3$AIC[x] <- AIC(models.3.list[[x]]
x <- x + 1 
models.3.list[[x]] <- gam(pa ~  +s(Init.count) + s(Open.500) , data=data.3)
models.3$AIC[x] <- AIC(models.3.list[[x]]
x <- x + 1 
models.3.list[[x]] <- gam(pa ~  +s(Init.count) + s(Closed.500s(, data=data.3)
models.3$AIC[x] <- AIC(models.3.list[[x]]
x <- x + 1 
models.3.list[[x]] <- gam(pa ~  +s(Init.count) + s(Moorland.500s(, data=data.3)
models.3$AIC[x] <- AIC(models.3.list[[x]]
x <- x + 1 
models.3.list[[x]] <- gam(pa ~  +s(Init.count) + s(males_94), data=data.3)
models.3$AIC[x] <- AIC(models.3.list[[x]]
x <- x + 1 
models.3.list[[x]] <- gam(pa ~  +s(Init.count) + s(leks_94), data=data.3)
models.3$AIC[x] <- AIC(models.3.list[[x]]
x <- x + 1 
models.3.list[[x]] <- gam(pa ~  +s(Open.500) + s(Closed.500s(, data=data.3)
models.3$AIC[x] <- AIC(models.3.list[[x]]
x <- x + 1 
models.3.list[[x]] <- gam(pa ~  +s(Open.500) + s(Moorland.500s(, data=data.3)
models.3$AIC[x] <- AIC(models.3.list[[x]]
x <- x + 1 
models.3.list[[x]] <- gam(pa ~  +s(Open.500) + s(males_94)  , data=data.3)
models.3$AIC[x] <- AIC(models.3.list[[x]]
x <- x + 1 
models.3.list[[x]] <- gam(pa ~  +s(Open.500) + s(leks_94)   , data=data.3)
models.3$AIC[x] <- AIC(models.3.list[[x]]
x <- x + 1 
models.3.list[[x]] <- gam(pa ~  +s(Closed.500) + s(Moorland.500s(, data=data.3)
models.3$AIC[x] <- AIC(models.3.list[[x]]
x <- x + 1 
models.3.list[[x]] <- gam(pa ~  +s(Closed.500) + s(males_94), data=data.3)
models.3$AIC[x] <- AIC(models.3.list[[x]]
x <- x + 1 
models.3.list[[x]] <- gam(pa ~  +s(Closed.500) + s(leks_94), data=data.3)
models.3$AIC[x] <- AIC(models.3.list[[x]]
x <- x + 1 
models.3.list[[x]] <- gam(pa ~  +s(Moorland.500) + s(males_94), data=data.3)
models.3$AIC[x] <- AIC(models.3.list[[x]]
x <- x + 1 
models.3.list[[x]] <- gam(pa ~  +s(Moorland.500) + s(leks_94), data=data.3)
models.3$AIC[x] <- AIC(models.3.list[[x]]
x <- x + 1 
models.3.list[[x]] <- gam(pa ~  +s(males_94) + s(leks_94), data=data.3)
models.3$AIC[x] <- AIC(models.3.list[[x]]
x <- x + 1 
models.3.list[[x]] <- gam(pa ~  +s(Init.count) + s(Open.500) + s(Closed.500s(, data=data.3)
models.3$AIC[x] <- AIC(models.3.list[[x]]
x <- x + 1 
models.3.list[[x]] <- gam(pa ~  +s(Init.count) + s(Open.500) + s(Moorland.500s(, data=data.3)
models.3$AIC[x] <- AIC(models.3.list[[x]]
x <- x + 1 
models.3.list[[x]] <- gam(pa ~  +s(Init.count) + s(Open.500) + s(males_94), data=data.3)
models.3$AIC[x] <- AIC(models.3.list[[x]]
x <- x + 1 
models.3.list[[x]] <- gam(pa ~  +s(Init.count) + s(Open.500) + s(leks_94), data=data.3)
models.3$AIC[x] <- AIC(models.3.list[[x]]
x <- x + 1 
models.3.list[[x]] <- gam(pa ~  +s(Init.count) + s(Closed.500) + s(Moorland.500s(, data=data.3)
models.3$AIC[x] <- AIC(models.3.list[[x]]
x <- x + 1 
models.3.list[[x]] <- gam(pa ~  +s(Init.count) + s(Closed.500) + s(males_94), data=data.3)
models.3$AIC[x] <- AIC(models.3.list[[x]]
x <- x + 1 
models.3.list[[x]] <- gam(pa ~  +s(Init.count) + s(Closed.500) + s(leks_94), data=data.3)
models.3$AIC[x] <- AIC(models.3.list[[x]]
x <- x + 1 
models.3.list[[x]] <- gam(pa ~  +s(Init.count) + s(Moorland.500) + s(males_94), data=data.3)
models.3$AIC[x] <- AIC(models.3.list[[x]]
x <- x + 1 
models.3.list[[x]] <- gam(pa ~  +s(Init.count) + s(Moorland.500) + s(leks_94), data=data.3)
models.3$AIC[x] <- AIC(models.3.list[[x]]
x <- x + 1 
models.3.list[[x]] <- gam(pa ~  +s(Init.count) + s(males_94) + s(leks_94), data=data.3)
models.3$AIC[x] <- AIC(models.3.list[[x]]
x <- x + 1 
models.3.list[[x]] <- gam(pa ~  +s(Open.500) + s(Closed.500) + s(Moorland.500s(, data=data.3)
models.3$AIC[x] <- AIC(models.3.list[[x]]
x <- x + 1 
models.3.list[[x]] <- gam(pa ~  +s(Open.500) + s(Closed.500) + s(males_94), data=data.3)
models.3$AIC[x] <- AIC(models.3.list[[x]]
x <- x + 1 
models.3.list[[x]] <- gam(pa ~  +s(Open.500) + s(Closed.500) + s(leks_94), data=data.3)
models.3$AIC[x] <- AIC(models.3.list[[x]]
x <- x + 1 
models.3.list[[x]] <- gam(pa ~  +s(Open.500) + s(Moorland.500) + s(males_94), data=data.3)
models.3$AIC[x] <- AIC(models.3.list[[x]]
x <- x + 1 
models.3.list[[x]] <- gam(pa ~  +s(Open.500) + s(Moorland.500) + s(leks_94), data=data.3)
models.3$AIC[x] <- AIC(models.3.list[[x]]
x <- x + 1 
models.3.list[[x]] <- gam(pa ~  +s(Open.500) + s(males_94)   s(leks_94)   , data=data.3)
models.3$AIC[x] <- AIC(models.3.list[[x]]
x <- x + 1 
models.3.list[[x]] <- gam(pa ~  +s(Closed.500) + s(Moorland.500) + s(males_94), data=data.3)
models.3$AIC[x] <- AIC(models.3.list[[x]]
x <- x + 1 
models.3.list[[x]] <- gam(pa ~  +s(Closed.500) + s(Moorland.500) + s(leks_94), data=data.3)
models.3$AIC[x] <- AIC(models.3.list[[x]]
x <- x + 1 
models.3.list[[x]] <- gam(pa ~  +s(Closed.500) + s(males_94) + s(leks_94), data=data.3)
models.3$AIC[x] <- AIC(models.3.list[[x]]
x <- x + 1 
models.3.list[[x]] <- gam(pa ~  +s(Moorland.500) + s(males_94) + s(leks_94), data=data.3)
models.3$AIC[x] <- AIC(models.3.list[[x]]
x <- x + 1 
models.3.list[[x]] <- gam(pa ~  +s(Init.count) + s(Open.500) + s(Closed.500) + s(Moorland.500s(, data=data.3)
models.3$AIC[x] <- AIC(models.3.list[[x]]
x <- x + 1 
models.3.list[[x]] <- gam(pa ~  +s(Init.count) + s(Open.500) + s(Closed.500) + s(males_94), data=data.3)
models.3$AIC[x] <- AIC(models.3.list[[x]]
x <- x + 1 
models.3.list[[x]] <- gam(pa ~  +s(Init.count) + s(Open.500) + s(Closed.500) + s(leks_94), data=data.3)
models.3$AIC[x] <- AIC(models.3.list[[x]]
x <- x + 1 
models.3.list[[x]] <- gam(pa ~  +s(Init.count) + s(Open.500) + s(Moorland.500) + s(males_94), data=data.3)
models.3$AIC[x] <- AIC(models.3.list[[x]]
x <- x + 1 
models.3.list[[x]] <- gam(pa ~  +s(Init.count) + s(Open.500) + s(Moorland.500) + s(leks_94), data=data.3)
models.3$AIC[x] <- AIC(models.3.list[[x]]
x <- x + 1 
models.3.list[[x]] <- gam(pa ~  +s(Init.count) + s(Open.500) + s(males_94) + s(leks_94), data=data.3)
models.3$AIC[x] <- AIC(models.3.list[[x]]
x <- x + 1 
models.3.list[[x]] <- gam(pa ~  +s(Init.count) + s(Closed.500) + s(Moorland.500) + s(males_94), data=data.3)
models.3$AIC[x] <- AIC(models.3.list[[x]]
x <- x + 1 
models.3.list[[x]] <- gam(pa ~  +s(Init.count) + s(Closed.500) + s(Moorland.500) + s(leks_94), data=data.3)
models.3$AIC[x] <- AIC(models.3.list[[x]]
x <- x + 1 
models.3.list[[x]] <- gam(pa ~  +s(Init.count) + s(Closed.500) + s(males_94) + s(leks_94), data=data.3)
models.3$AIC[x] <- AIC(models.3.list[[x]]
x <- x + 1 
models.3.list[[x]] <- gam(pa ~  +s(Init.count) + s(Moorland.500) + s(males_94) + s(leks_94), data=data.3)
models.3$AIC[x] <- AIC(models.3.list[[x]]
x <- x + 1 
models.3.list[[x]] <- gam(pa ~  +s(Open.500) + s(Closed.500) + s(Moorland.500) + s(males_94), data=data.3)
models.3$AIC[x] <- AIC(models.3.list[[x]]
x <- x + 1 
models.3.list[[x]] <- gam(pa ~  +s(Open.500) + s(Closed.500) + s(Moorland.500) + s(leks_94), data=data.3)
models.3$AIC[x] <- AIC(models.3.list[[x]]
x <- x + 1 
models.3.list[[x]] <- gam(pa ~  +s(Open.500) + s(Closed.500) + s(males_94) + s(leks_94), data=data.3)
models.3$AIC[x] <- AIC(models.3.list[[x]]
x <- x + 1 
models.3.list[[x]] <- gam(pa ~  +s(Open.500) + s(Moorland.500) + s(males_94) + s(leks_94), data=data.3)
models.3$AIC[x] <- AIC(models.3.list[[x]]
x <- x + 1 
models.3.list[[x]] <- gam(pa ~  +s(Closed.500) + s(Moorland.500) + s(males_94) + s(leks_94), data=data.3)
models.3$AIC[x] <- AIC(models.3.list[[x]]
x <- x + 1 
models.3.list[[x]] <- gam(pa ~  +s(Init.count) + s(Open.500) + s(Closed.500) + s(Moorland.500) + s(males_94), data=data.3)
models.3$AIC[x] <- AIC(models.3.list[[x]]
x <- x + 1 
models.3.list[[x]] <- gam(pa ~  +s(Init.count) + s(Open.500) + s(Closed.500) + s(Moorland.500) + s(leks_94), data=data.3)
models.3$AIC[x] <- AIC(models.3.list[[x]]
x <- x + 1 
models.3.list[[x]] <- gam(pa ~  +s(Init.count) + s(Open.500) + s(Closed.500) + s(males_94) + s(leks_94), data=data.3)
models.3$AIC[x] <- AIC(models.3.list[[x]]
x <- x + 1 
models.3.list[[x]] <- gam(pa ~  +s(Init.count) + s(Open.500) + s(Moorland.500) + s(males_94) + s(leks_94), data=data.3)
models.3$AIC[x] <- AIC(models.3.list[[x]]
x <- x + 1 
models.3.list[[x]] <- gam(pa ~  +s(Init.count) + s(Closed.500) + s(Moorland.500) + s(males_94) + s(leks_94), data=data.3)
models.3$AIC[x] <- AIC(models.3.list[[x]]
x <- x + 1 
models.3.list[[x]] <- gam(pa ~  +s(Open.500) + s(Closed.500) + s(Moorland.500) + s(males_94) + s(leks_94), data=data.3)
models.3$AIC[x] <- AIC(models.3.list[[x]]
x <- x + 1 
models.3.list[[x]] <- gam(pa ~  + s(Init.count) + s(Open.500) + s(Closed.500) + s(Moorland.500) + s(males_94) + s(leks_94), data=data.3)
 x <- x +1

##### Best models #####

models.1$D.AIC <- models.1$AIC - models.1$AIC[which.min(models.1$AIC)]

models.1$lik <- exp(-0.5*models.1$D.AIC)
models.1$weight <- models.1$lik/sum(models.1$lik, na.rm = T)

models.1.best <- subset(models.1, subset=(models.1$D.AIC < 2))

models.3$D.AIC <- models.3$AIC - models.3$AIC[which.min(models.3$AIC)]

models.3$lik <- exp(-0.5*models.3$D.AIC)
models.3$weight <- models.3$lik/sum(models.3$lik, na.rm = T)

models.3.best <- subset(models.3, subset=(models.3$D.AIC < 2))

write.csv(models.1, "mods_1_nomoor_OCT_2012.csv", row.names=F)
write.csv(models.3, "mods_3_nomoor_OCT_2012.csv", row.names=F)

write.csv(models.1.best, "best_mods_1_nomoor_OCT_2012.csv", row.names=F)
write.csv(models.3.best, "best_mods_3_nomoor_OCT_2012.csv", row.names=F)
