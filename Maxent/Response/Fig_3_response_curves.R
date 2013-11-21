source("response_data_500.R")
source("response_data_2k.R")

#png("PapGrouseYears2_v4_Fig_3.png", height=800, width=1000)
png("PapGrouseYears2_v4_Fig_3_wide.png", height=800, width=1000)

Panel <- function(panel)
{
  par(xpd=NA)
  u <- par('usr')
  text(u[1]-(u[2]-u[1])*.01, u[4]+(u[4]-u[3])*.03, panel, cex=2)
  par(xpd=FALSE)
}

#layout(matrix(c(1,2,3,4,5,6,7,8,9,10), nrow=2))
layout(matrix(c(1,2,3,4,5,6,7,8), nrow=2, ncol=4))

#oldpar <- par(mar=c(0,0,0,0))

#plot(c(1:10) ~ c(1:10), type="n", axes=F, xlab="", ylab="")
#text(5,6, "0.5 km", cex=2, font=4)

#plot(c(1:10) ~ c(1:10), type="n", axes=F, xlab="", ylab="")
#text(5,6, "2 km", cex=2, font=4)

#par(oldpar)

oldpar <- par( mar=c(15,9,2,1), cex.lab=2, cex.axis=1.5, mgp=c(5,1,0))
# Grazed

# 0.5 km
plot(m.94$Grazed ~ m.94.x$Grazed, lty=1, lwd=2, type="l", ylim=c(0.20, 0.60), xlab="Proportion of \ngrazed land", ylab="Relative habitat \n suitability")
lines(m.00$Grazed ~ m.00.x$Grazed, lty=2, lwd=2, type="l")
lines(m.08$Grazed ~ m.08.x$Grazed, lty=3, lwd=2, type="l")
Panel("(a)")

# 2 km
plot(m.94.2k$Grazed ~ m.94.2k.x$Grazed, lty=1, lwd=2, type="l",  xlab="Proportion of \ngrazed land", ylab="Relative habitat \n suitability")
lines(m.00.2k$Grazed ~ m.00.2k.x$Grazed, lty=2, lwd=2, type="l")
lines(m.08.2k$Grazed ~ m.08.2k.x$Grazed, lty=3, lwd=2, type="l")
Panel("(e)")

# Grouse moor

# 0.5 km
plot(m.94$Grouse.moor ~ m.94.x$Grouse.moor, lty=1, lwd=2, type="l", ylim=c(0, 0.65),  xlab="Proportion of \ngrouse moor", ylab="Relative habitat \n suitability")
lines(m.00$Grouse.moor ~ m.00.x$Grouse.moor, lty=2, lwd=2, type="l")
lines(m.08$Grouse.moor ~ m.08.x$Grouse.moor, lty=3, lwd=2, type="l")
Panel("(b)")

# 2 km
plot(m.94.2k$Grouse.moor ~ m.94.2k.x$Grouse.moor, lty=1, lwd=2, type="l", ylim=c(0, 0.65),  xlab="Proportion of \ngrouse moor", ylab="Relative habitat \n suitability")
lines(m.00.2k$Grouse.moor ~ m.00.2k.x$Grouse.moor, lty=2, lwd=2, type="l")
lines(m.08.2k$Grouse.moor ~ m.08.2k.x$Grouse.moor, lty=3, lwd=2, type="l")
Panel("(f)")

# Open canopy

# 0.5 km
plot(m.94.2k$Grouse.moor ~ m.94.2k.x$Grouse.moor, lty=1, lwd=2, type="l", ylim=c(0, 0.65),  xlab="Proportion of \ngrouse moor", ylab="Relative habitat \n suitability")
lines(m.00.2k$Grouse.moor ~ m.00.2k.x$Grouse.moor, lty=2, lwd=2, type="l")
lines(m.08.2k$Grouse.moor ~ m.08.2k.x$Grouse.moor, lty=3, lwd=2, type="l")
Panel("(c)")

# 2 km
plot(m.94.2k$Grouse.moor ~ m.94.2k.x$Grouse.moor, lty=1, lwd=2, type="l", ylim=c(0, 0.65),  xlab="Proportion of \ngrouse moor", ylab="Relative habitat \n suitability")
lines(m.00.2k$Grouse.moor ~ m.00.2k.x$Grouse.moor, lty=2, lwd=2, type="l")
lines(m.08.2k$Grouse.moor ~ m.08.2k.x$Grouse.moor, lty=3, lwd=2, type="l")
Panel("(g)")

# Closed canopy

# 0.5 km
plot(m.94$Closed ~ m.94.x$Closed, lty=1, lwd=2, type="l", xlab="Proportion of \nclosed forestry", ylab="Relative habitat \n suitability")
lines(m.00$Closed ~ m.00.x$Closed, lty=2, lwd=2, type="l")
lines(m.08$Closed ~ m.08.x$Closed, lty=3, lwd=2, type="l")
Panel("(d)")

# 2km
plot(m.94.2k$Closed ~ m.94.2k.x$Closed, lty=1, lwd=2, type="l", ylim=c(0,0.65), xlab="Proportion of \nclosed forestry", ylab="Relative habitat \n suitability")
lines(m.00.2k$Closed ~ m.00.2k.x$Closed, lty=2, lwd=2, type="l")
lines(m.08.2k$Closed ~ m.08.2k.x$Closed, lty=3, lwd=2, type="l")
Panel("(h)")

par(oldpar)

dev.off()
