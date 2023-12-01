#更改数据框变量名
fix(Rdatadrug1)
names(Rdatadrug1)[1] <- "dose"
colnames(Rdatatest1) <- c("num1", "sum1")


#增加数据框变量
mydata <- transform(mydata, sumx = x1 + x2, meanx = (x1 + x2) / 2)
#删除数据框变量
newdata<-totaldata[c(-3)]
totaldata[3]<-NULL

#判断是否存在缺失值
is.na(Rdatadrug1)
#排除缺失值的计算
sum(Rdatadrug1$dose, na.rm = TRUE)
#删除缺失值
newdata = na.omit(Rdatadrug1)

#增加列
totaldata<-merge(Rdatadrug1,Rdatadrug2,by="dose")
#增加列（相同行数及顺序）
totaldata<-cbind(Rdatadrug1,Rdatadrug2)
#增加行（必须相同变量，可追加或删减至其相同）
total<-rbind(Rdatadrug1,Rdatadrug2)

#筛选数据

