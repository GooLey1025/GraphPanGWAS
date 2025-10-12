#!/usr/bin/env Rscript

args=commandArgs(TRUE)
if (length(args) != 2) {
          print ("usage: <gwas file(4 column)>   <out.prefix> ")
  q()
}

input <- args[1]
prefix <- args[2]


#input <- "./01.emmax/emmax.out.ps.manht_input"
#prefix <- "test"
library(CMplot)

gwas<- read.table(input, header = T);

cutoff <- 0.05/nrow(gwas)

sigSNP <- gwas[gwas[,4] < cutoff,]

write.table( sigSNP, file = paste( prefix, "sigSite.out", sep =  "."), row.names = F, quote = F)

# lambda计�~W
## For z-scores, just square them
# chisq <- data$z^2
## For chi-squared values, keep as is
# chisq <- data$chisq
## For p-values, calculate chi-squared statistic
# while qchisq(1-assoc.df$P,1) fails to convert small p-value (Inf in this case) 
# chisq <- qchisq(1-gwas[[4]],1, lower.tail=T)
chisq <- qchisq(gwas[[4]],1,lower.tail=FALSE)
lambda <- median(chisq, na.rm=TRUE)/qchisq(0.5, 1)
write(lambda , file = paste0(prefix, ".lambda.txt"))

png(paste(prefix, "_manhattan_threshold.png", sep = ""), width=960, height=480)
CMplot(gwas,
       plot.type = "m", ## m:Manhattan�~Lc: circle-Manhattan
       type = "p" , ## "p" (point), "l" (cross line), "h" (vertical lines)
       LOG10 = T, # -log10转�~M�
       col = c("blue4", "orange3"),
       cex = 0.3,
       ylab.pos = 2, ## the distance between ylab and yaxis
       axis.cex=1,
       threshold = c(0.01,0.05)/nrow(gwas), ## �~X��~Q~W�~@��~X~H�~@�
       threshold.col=c('grey','black'),  ## �~X~H�~@�线�~\�~I�
       threshold.lty = c(1,2), ## �~X~H�~@�线线�~^~K
       threshold.lwd = c(1,1), ## �~X~H�~@�线�~W�~F
       amplify = T,  ## �~T�大�~X��~Q~WSNP
       signal.cex = c(1,1), ## �~B�大�~O
       signal.pch = c(20,20), ## �~B�形�~J�
       signal.col = c("red","blue"), ## �~B��~\�~I�
       file.output = F,
       )

dev.off()


png(paste(prefix, "_manhattan.png", sep = ""), width=960, height=480)
CMplot(gwas,plot.type = "m",
       LOG10 = T,
       col = c("blue4", "orange3"),
       cex = 0.3,
       ylab.pos = 2,
       axis.cex  = 1,
       threshold=NULL,
       file.output = F )

dev.off()

## �~[��~S~H顿�~[�1
pdf(paste(prefix, "_manhattan_threshold.pdf", sep = ""), width=10, height=5)
CMplot(gwas,
       plot.type = "m",
       LOG10 = T,
       col = c("blue4", "orange3"),
       cex = 0.3,
       ylab.pos = 2,
       axis.cex = 1,
       threshold = c(0.01,0.05)/nrow(gwas), ## �~X��~Q~W�~@��~X~H�~@�
       threshold.col=c('grey','black'),  ## �~X~H�~@�线�~\�~I�
       threshold.lty = c(1,2), ## �~X~H�~@�线线�~^~K
       threshold.lwd = c(1,1), ## �~X~H�~@�线�~W�~F
       amplify = T,  ## �~T�大�~X��~Q~WSNP
       signal.cex = c(1,1), ## �~B�大�~O
       signal.pch = c(20,20), ## �~B�形�~J�
       signal.col = c("red","blue"), ## �~B��~\�~I�
       file.output = F )
dev.off()

## �~[��~S~H顿�~[�2

pdf(paste(prefix, "_manhattan.pdf", sep = ""), width=10, height=5)
CMplot(gwas,plot.type = "m",
       LOG10 = T,
       col = c("blue4", "orange3"),
       cex = 0.3,
       ylab.pos = 2,
       axis.cex = 1,
       threshold=NULL,
       file.output = F )

dev.off()

## QQplot 
pdf(paste(prefix, "_qqplot.pdf", sep = ""), width=10, height=10)
CMplot(gwas,
       plot.type = "q", ## �~X�~H�QQplot
       box = T, ## �~X��~P��~J| 边�~F
       conf.int=T, ## �~X��~P��~X�~H�置信�~L��~W�
       conf.int.col=NULL, ## 置信�~L��~W��~\�~I�
       threshold.col="red", ## 对�~R线�~\�~I�
       threshold.lty=2,  ## 线�~^~K
       cex = 0.8,
       ylab.pos = 2,
       axis.cex = 1,
       main = "QQ-plot",
       file.output = F )
dev.off()
