m.1 <- read.csv("mods_1_nomoor_OCT_2012.csv")
m.2 <- read.csv("mods_2_nomoor_OCT_2012.csv")
m.3 <- read.csv("mods_3_nomoor_OCT_2012.csv")
v <- c("Init.count", "Moorland", "Open", "Closed", "Connectivity")
weights <- data.frame(Variable=v, mod.1=numeric(5), mod.2=numeric(5), mod.3=numeric(5))


weights$mod.1[1] <- c(m.1$weight[1] + m.1$weight[7] + m.1$weight[9] + m.1$weight[11] + m.1$weight[23] + m.1$weight[25] + m.1$weight[30] + m.1$weight[46])

weights$mod.1[2] <- c(m.1$weight[4] + m.1$weight[9] + m.1$weight[13] + m.1$weight[20] + m.1$weight[23] + m.1$weight[30] + m.1$weight[36] + m.1$weight[46])

weights$mod.1[3] <- c(m.1$weight[2] + m.1$weight[7] + m.1$weight[13] + m.1$weight[15] + m.1$weight[23] + m.1$weight[25] + m.1$weight[36] + m.1$weight[46])

weights$mod.1[4] <- NA

weights$mod.1[5] <- c(m.1$weight[6] + m.1$weight[11] + m.1$weight[15] + m.1$weight[20] + m.1$weight[25] + m.1$weight[30] + m.1$weight[36] + m.1$weight[46])



weights$mod.2[1] <- c(m.2$weight[1] + m.2$weight[7] + m.2$weight[8] + m.2$weight[9] + m.2$weight[11] + m.2$weight[22] + m.2$weight[23] + m.2$weight[25] + m.2$weight[26] + m.2$weight[28] + m.2$weight[30] + m.2$weight[42] + m.2$weight[44] + m.2$weight[46] + m.2$weight[49] + m.2$weight[58])

weights$mod.2[2] <- c(m.2$weight[4] + m.2$weight[9] + m.2$weight[13]+ m.2$weight[16] + m.2$weight[20] + m.2$weight[23] + m.2$weight[26] + m.2$weight[30] + m.2$weight[32] + m.2$weight[36] + m.2$weight[39] + m.2$weight[42] + m.2$weight[46] + m.2$weight[49] + m.2$weight[53] + m.2$weight[58])

weights$mod.2[2] <- c(m.2$weight[2] + m.2$weight[7] + m.2$weight[12] + m.2$weight[13] + m.2$weight[15] + m.2$weight[22] + m.2$weight[23] + m.2$weight[25] + m.2$weight[32] + m.2$weight[34] + m.2$weight[36] + m.2$weight[42] + m.2$weight[44] + m.2$weight[46] + m.2$weight[53] + m.2$weight[58])

weights$mod.2[4] <- c(m.2$weight[3] + m.2$weight[8] + m.2$weight[12] + m.2$weight[16] + m.2$weight[18] + m.2$weight[22] + m.2$weight[26] + m.2$weight[28] + m.2$weight[32] + m.2$weight[34] + m.2$weight[39] + m.2$weight[42] + m.2$weight[44] + m.2$weight[49] + m.2$weight[53] + m.2$weight[58])

weights$mod.2[5] <- c(m.2$weight[6] + m.2$weight[11] + m.2$weight[15] + m.2$weight[18] + m.2$weight[20] + m.2$weight[25] + m.2$weight[28] + m.2$weight[30] + m.2$weight[34] + m.2$weight[36] + m.2$weight[39] + m.2$weight[44] + m.2$weight[46] + m.2$weight[49] + m.2$weight[53] + m.2$weight[58])



weights$mod.3[1] <- c(m.3$weight[1] + m.3$weight[7] + m.3$weight[8] + m.3$weight[11] + m.3$weight[22] + m.3$weight[25] + m.3$weight[28] + m.3$weight[44])

weights$mod.3[2] <- NA

weights$mod.3[3] <- c(m.3$weight[2] + m.3$weight[7] + m.3$weight[12] + m.3$weight[15] + m.3$weight[22] + m.3$weight[25] + m.3$weight[34] + m.3$weight[44])

weights$mod.3[4] <- c(m.3$weight[3] + m.3$weight[8] + m.3$weight[12] + m.3$weight[18] + m.3$weight[22] + m.3$weight[28] + m.3$weight[34] + m.3$weight[44])

weights$mod.3[5] <- c(m.3$weight[6] + m.3$weight[11] + m.3$weight[15] + m.3$weight[18] + m.3$weight[25] + m.3$weight[28] + m.3$weight[34] + m.3$weight[44])



#write.csv(weights, "Variable_weights_OCT_2012.csv", row.names=F)
