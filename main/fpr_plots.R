fprs = read.table("Documents/UCSF_year1/algorithms/fprs.txt")

mean(fprs$V1)
sd(fprs$V1)

library(ggplot2)

fprs$V3 = as.factor(fprs$V3)


ggplot(fprs, aes(x=V2, y=V1, fill=V2)) + 
  geom_boxplot(aes()) + xlab("Gap Extension Penalty") + ylab("FPR")

install.packages("scatterplot3d") # Install
library("scatterplot3d")

fprs$V3 = as.numeric(fprs$V3)
colors <- c("#999999", "#E69F00", "#56B4E9", 	"#FFC0CB")

colors <- colors[as.numeric(fprs$V3)]

scatterplot3d(fprs[,1:3], angle=85,
              xlab = "FPR",
              ylab = "Gap Opening Penalty",
              zlab = "Gap Extension Penalty", color=colors, pch=17, 
              legend = levels(fprs$V3))
legend(legend = levels(fprs$V3))

blosum = c(0.12, 0.32, 0.42, 0.66, 0.8, 0.9, 0.92, 0.94, 0.96, 1.0, 1.0)
noramlized = c(0.26000000000000001, 0.44, 0.47999999999999998, 0.59999999999999998, 0.66000000000000003, 0.68000000000000005, 0.80000000000000004, 0.85999999999999999, 0.95999999999999996, 0.97999999999999998, 0.97999999999999998)
tpr = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.8, 1.0)
fpr = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.8, 1.0)

plot(blosum, tpr, type="l", col="royalblue", xlim=c(0,1), ylim=c(0,1), xlab="fpr", ylab="tpr")
lines(noramlized, tpr, col="red")
lines(fpr, tpr, col="black")



