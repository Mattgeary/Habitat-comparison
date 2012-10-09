library(raster)
library(dismo)
library(ROCR)
library(sp)
work <- getwd()

BNG <- CRS("+init=epsg:27700")
bk.all <- read.csv("allpaxy.csv")
bk.pts <- bk.all[,1:2]
bk.94 <- data.frame(bk.all[,3])
bk.00 <- data.frame(bk.all[,4])
bk.08 <- data.frame(bk.all[,5])

bk.94.b <- bk.94
bk.94.b[bk.94.b > 0] <- 1
bk.00.b <- bk.00
bk.00.b[bk.00.b > 0] <- 1
bk.08.b <- bk.08
bk.08.b[bk.08.b > 0] <- 1


####
##tesy.sum.map <- focal(r.94, w=527, sum, na.rm=T, pad=T)##
####

alt <- raster("studyareadem.asc", crs=BNG)
sp.94 <- SpatialPointsDataFrame(bk.pts, bk.94)
sp.94.b <- SpatialPointsDataFrame(bk.pts, bk.94.b)
r.94 <- rasterize(sp.94, alt, background=0, "bk.all...3.")
r.94.b <- rasterize(sp.94.b, alt, background=0,"bk.all...3.")

sum.map <- focal(r.94, w=527, sum, na.rm=T, pad=T)
count.map <- focal(r.94.b, w=527, sum, na.rm=T, pad=T)
bk.all$males_94 <- extract(sum.map, bk.pts)- bk.94$ bk.all...3.
bk.all$leks_94 <- extract(count.map, bk.pts) - bk.94.bbk.94$ bk.all...3.

sp.00 <- SpatialPointsDataFrame(bk.pts, bk.00)
sp.00.b <- SpatialPointsDataFrame(bk.pts, bk.00.b)
r.00 <- rasterize(sp.00, alt, background=0, "bk.all...4.")
r.00.b <- rasterize(sp.00.b, alt, background=0, "bk.all...4.")

sum.map <- focal(r.00, w=527, sum, na.rm=T, pad=T)
count.map <- focal(r.00.b, w=527, sum, na.rm=T, pad=T)
bk.all$males_00 <- extract(sum.map, bk.pts) - bk.00$ bk.all...4.
bk.all$leks_00 <- extract(count.map, bk.pts) - bk.00.b$ bk.all...4.

sp.08 <- SpatialPointsDataFrame(bk.pts, bk.08, pad=T)
sp.08.b <- SpatialPointsDataFrame(bk.pts, bk.08.b, pad=T)
r.08 <- rasterize(sp.08, alt, background=0, "bk.all...5.")
r.08.b <- rasterize(sp.08.b, alt, background=0, "bk.all...5.")

sum.map <- focal(r.08, w=527, sum, na.rm=T, pad=T)
count.map <- focal(r.08.b, w=527, sum, na.rm=T, pad=T)
rm(sp.08, sp.08.b, r.08, r.08.b)
bk.all$males_08 <- extract(sum.map, bk.pts) - bk.08$ bk.all...5.
bk.all$leks_08 <- extract(count.map, bk.pts) - bk.08.b$ bk.all...5.

bk.all$males.conn.scaled.94.totleks <- bk.all$males_94/sum(bk.all$p94)
bk.all$leks.conn.scaled.94.totleks <- bk.all$males_94/sum(bk.all$p94)

bk.all$males.conn.scaled.00.totleks <- bk.all$males_00/sum(bk.all$p00)
bk.all$leks.conn.scaled.00.totleks <- bk.all$males_00/sum(bk.all$p00)

bk.all$males.conn.scaled.08.totleks <- bk.all$males_08/sum(bk.all$p08)
bk.all$leks.conn.scaled.08.totleks <- bk.all$males_08/sum(bk.all$p08)

bk.all$males.conn.scaled.94.totpop <- bk.all$males_94/sum(bk.all$Count_94)
bk.all$leks.conn.scaled.94.totpop <- bk.all$males_94/sum(bk.all$Count_94)

bk.all$males.conn.scaled.00.totpop <- bk.all$males_00/sum(bk.all$Count_00)
bk.all$leks.conn.scaled.00.totpop <- bk.all$males_00/sum(bk.all$Count_00)

bk.all$males.conn.scaled.08.totpop <- bk.all$males_08/sum(bk.all$Count_08)
bk.all$leks.conn.scaled.08.totpop <- bk.all$males_08/sum(bk.all$Count_08)

write.csv(bk.all, "allpaxy_connect_new_OCT_2012.csv")
#sp.gr.94 <- as.matrix(r.94)
#dists.94 <- spDists(sp.gr.94, longlat=F)

#d.which <- which(dists.94 < 15000)

#d.loc <- bk.94[d.which,1:2]
#d.count <- bk.94[d.which,i]


