setwd("D:/DEMOS_paper_source/results/DEMOS_C++/output")
library(dplyr)

#301 - 600 is the number with DM
dmfiles <- list.files(pattern = "simulated_DM")
incfiles <- list.files(pattern = "new")
pyrfiles <- list.files(pattern = "simulated_pyramid")

aggdm <- array(0, dim = c(dim(read.csv(dmfiles[1 + 300])),300))
agginc <- array(0, dim = c(dim(read.csv(incfiles[1])),300))
aggpyr <- array(0, dim = c(dim(read.csv(pyrfiles[1])),300))
finaldm <- matrix(0, nrow = dim(read.csv(dmfiles[1 + 300]))[1], ncol = 3 * dim(read.csv(dmfiles[1 + 300]))[2])
finalinc <- matrix(0, nrow = dim(read.csv(incfiles[1]))[1], ncol = 3 * dim(read.csv(incfiles[1]))[2])
finalpyr <- matrix(0, nrow = dim(read.csv(pyrfiles[1]))[1], ncol = 3 * dim(read.csv(pyrfiles[1]))[2])

for(i in 1:300)
{
    aggdm[,,i] <- as.matrix(read.csv(dmfiles[i + 300]))
    agginc[,,i] <- as.matrix(read.csv(incfiles[i]))
    aggpyr[,,i] <- as.matrix(read.csv(pyrfiles[i]))
}

colnames(finaldm) <- as.vector(rbind(colnames(read.csv(dmfiles[1])), paste(colnames(read.csv(dmfiles[1])), "lower", sep = "_"), paste(colnames(read.csv(dmfiles[1])), "upper", sep = "_")))
for(h in 1:ncol(aggdm))
{
    finaldm[,(h*3-2)] <- apply(aggdm[,h,], 1, mean)
	finaldm[,(h*3-1)] <- apply(aggdm[,h,], 1, function(x) quantile(x, 0.025))
	finaldm[,(h*3)] <- apply(aggdm[,h,], 1, function(x) quantile(x,0.975))
}
finaldm <- finaldm[,-c(2,3,5,6)]
finaldm[finaldm < 0] <- 0

colnames(finalinc) <- as.vector(rbind(colnames(read.csv(incfiles[1])), paste(colnames(read.csv(incfiles[1])), "lower", sep = "_"), paste(colnames(read.csv(incfiles[1])), "upper", sep = "_")))
for(h in 1:ncol(agginc))
{
    finalinc[,(h*3-2)] <- apply(agginc[,h,], 1, mean)
	finalinc[,(h*3-1)] <- apply(agginc[,h,], 1, function(x) quantile(x, 0.025))
	finalinc[,(h*3)] <- apply(agginc[,h,], 1, function(x) quantile(x, 0.975))
}
finalinc <- finalinc[,-c(2,3,5,6)]
finalinc[finalinc < 0] <- 0

colnames(finalpyr) <- as.vector(rbind(colnames(read.csv(pyrfiles[1])), paste(colnames(read.csv(pyrfiles[1])), "lower", sep = "_"), paste(colnames(read.csv(pyrfiles[1])), "upper", sep = "_")))
for(h in 1:ncol(aggpyr))
{
    finalpyr[,(h*3-2)] <- apply(aggpyr[,h,], 1, mean)
	finalpyr[,(h*3-1)] <- apply(aggpyr[,h,], 1, function(x) quantile(x, 0.025))
	finalpyr[,(h*3)] <- apply(aggpyr[,h,], 1, function(x) quantile(x, 0.975))
}
finalpyr <- finalpyr[,-c(2,3,5,6)]
finalpyr[finalpyr < 0] <- 0

write.csv(finalpyr, "DEMOSpaper_populationcounts.csv", row.names = FALSE)
write.csv(finalinc, "DEMOSpaper_incidence.csv", row.names = FALSE)
write.csv(finaldm, "DEMOSpaper_casenumbers.csv", row.names = FALSE)
