opar<- par(no.readonly = TRUE)
par(fig=c(0,0.8,0,0.8))
plot(Rdatadrug1$dose,Rdatadrug1$drugA,
     xlab = "dosage",ylab = "DrugA")

par(fig = c(0, 0.82, 0.4, 1), new = TRUE)
boxplot(Rdatadrug1$dose, horizontal = TRUE, axes = FALSE)

par(fig = c(0.6, 1, 0, 0.82), new = TRUE)
boxplot(Rdatadrug1$drugA, horizontal = FALSE, axes = FALSE)
par(opar)