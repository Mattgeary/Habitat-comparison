require(mgcv)
require(ncf)

BK.y <- read.csv("BK_Years.csv")[,-1]
BK.y[BK.y > 0] <- 1
BK.dist <- dist(t(BK.y), method="euclidean")
BK.cluster <- hclust(BK.dist, method="ward")

dist.cats <- numeric(0)
for (i in 18:1){
 dist.cats <- append(dist.cats, seq(1:i))
}

distances <- vector("list", 18)
for (i in 1:171){
 distances[[dist.cats[i]]] <- append(distances[[dist.cats[i]]], BK.dist[i])
}

d.x <- c(1:18)
mean.distances <- mapply(mean, distances)
sd.distances <- mapply(sd, distances)
sd.distances[18] <- 0

#png("PapGrouse Years2_v4_Fig_2.png", width = 800, height = 800)
oldpar <- par(mar=c(5,8,1,1))
plot(mean.distances ~ d.x, type="l", ylim=c(6, 11), xlab = "Proximity (years)", ylab = "Euclidean distance between\noccupancy patterns", lwd=3, cex.axis = 1.5, cex.lab=1.5)
lines((mean.distances + sd.distances) ~ d.x, lty=2, col="red", lwd= 3)
lines((mean.distances - sd.distances) ~ d.x, lty=2, col="red", lwd = 3)
par(oldpar)
#dev.off()
