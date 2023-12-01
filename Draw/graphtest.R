#plot()函数默认将输入数据视为数值型，在尝试绘制非数值型数据时会导致错误
#可以使用barplot()

# x <- c("A", "B", "C", "D")
# y <- c(10, 15, 7, 12)
# 
# barplot(y, names.arg = x, main = "Bar Chart", xlab = "Categories", ylab = "Values")


# # 调用额外画板
# dev.new()
# NULL
# x <- c("A", "B", "C", "D")
# y <- c(10, 15, 7, 12)
# barplot(y, names.arg = x, main = "Bar Chart", xlab = "Categories", ylab = "Values")


# # 创建一个向量作为示例数据
# data <- c(22, 34, 16, 28, 30, 20, 24, 18, 21, 31)
# 
# # 使用一些常用参数来自定义直方图的外观
# hist(data,
#      breaks = 5,  # 分成5个分箱
#      col = "blue",  # 设置直方图的填充颜色为蓝色
#      border = "white",  # 设置直方图的边界颜色为白色
#      main = "Histogram",  # 设置直方图的标题
#      xlab = "Values",  # 设置x轴标签
#      ylab = "Frequency"  # 设置y轴标签
# )


# # 散点图矩阵
# # 创建数据
# data <- iris[, 1:4]   # 使用鸢尾花数据集的前四列作为示例数据
# # 使用pair()绘制散点图矩阵
# pairs(data)


# # 绘制三维散点图
# # 创建三个变量作为示例数据
# x <- c(1, 2, 3, 4, 5)
# y <- c(6, 7, 8, 9, 10)
# z <- c(11, 12, 13, 14, 15)
# scatterplot3d(x, y, z,
#               color = "blue",  # 设置散点的颜色为蓝色
#               pch = 16,  # 设置散点的形状为实心圆
#               main = "3D Scatter Plot",  # 设置图标题
#               xlab = "X",  # 设置x轴标签
#               ylab = "Y",  # 设置y轴标签
#               zlab = "Z"  # 设置z轴标签
# )


# # 绘制三维散点图
# # 调用额外画板
# dev.new()
# par(family='STKaiti')     # 设置中文字体，防乱码
# x <- c(1, 2, 3, 4, 5)
# y <- c(2, 4, 3, 5, 1)
# plot(x, y, type = "b", pch = 16, col = "blue", main = "示例")
# points(x, y^2, pch = 17, col = "red")
# lines(x, y, col = "green")
# legend("topright", legend = c("数据1", "数据2", "线条"),
#        pch = c(16, 17, NA), col = c("blue", "red", "green"),
#        lty = c(NA, NA, 1), title = "图例标题")


# # 创建散点图并添加标签
# plot(x, y)
# text(x = c(1, 2, 3), y = c(2, 4, 3), labels = c("Label 1", "Label 2", "Label 3"))
# # 添加带有样式的标签
# plot(x, y)
# text(x = c(1, 2, 3), y = c(2, 4, 3), labels = c("Label 1", "Label 2", "Label 3"), pos = 4, col = "blue", font = 2, cex = 1.5)


# # 绘制直线
# plot(x, y)
# lines(x = c(1, 2, 3), y = c(2, 4, 3), col = "blue")
# # 绘制折线
# plot(x, y)
# lines(x = c(1, 2, 3, 4), y = c(2, 4, 3, 5), type = "b", col = "red")
# # 绘制阶梯函数
# plot(x, y)
# lines(x = c(1, 2, 3, 4), y = c(2, 4, 3, 5), type = "h", col = "green")
# # 绘制自定义线型和宽度的线条
# plot(x, y)
# lines(x, y, lty = 2, lwd = 2, col = "purple")


# #更改数据框列名
# colnames(Rdatatest1) <- c("num1", "sum1")

# plot参数修改示例
opar<- par(no.readonly = TRUE)
# par(mfrow=c(1,2))  #图形矩阵
# layout(matrix(c(1, 1, 2, 3), 2, 2, byrow = TRUE))  #图形矩阵
plot(
  Rdatadrug1$dose,
  Rdatadrug1$drugA,
  type = "b",
  lty = 2,
  pch = 17,
  cex = 1.5,
  col = "red",
  col.lab = "blue",
  main = "",
  xlab = "",
  ylab = "",
  xlim = c(0, 60),
  ylim = c(0, 70)
)
title(
  ylab = "Drug Response",
  line = 2,
  col.lab = "blue",
  cex.lab = 1.2
)
title(
  xlab = "Dosage",
  line = 2,
  col.lab = "blue",
  cex.lab = 1.2
)
title(
  main = "Clinical Trials for DrugA",
  line = 1,
  col.main = "red",
  cex.main = 1.2
)
par(opar)