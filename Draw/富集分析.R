# 安装所需R包
install.packages("stringr")
install.packages("ggplot2")

#加载所需R包
library(ggplot2)
library(stringr)

# 读取数据
dt1 <- read.table("路径", sep = "\t", header = T)

# 查看数据框的前6行
head(dt1, 6)

# 选取数据前20行用于作图
dt2 <- dt1[1:20,]

#指定纵轴标签顺序，按照输入文件的顺序排序，否则默认按照首字母顺序，同时逆序绘制，保持与表格顺序一致
dt2$id <-
  factor(dt2$id, levels = rev(unique(dt2$id)), ordered = TRUE)

# 建立数据与图形的映射关系，即确定点的坐标，绘制散点图
p1 <- ggplot(dt2, aes(ratio, id)) + geom_point()
p1

# 建立数据与点的大小的映射关系
# 建立颜色与数据的关系，这里让点按照数据大小显示不同的颜色
p2 <-
  ggplot(dt2, aes(ratio, id)) + geom_point(aes(size = num, color = qvalue), sharp =
                                             16)
p2

# 设置标签的折叠（字符数）
y <- str_wrap(dt2$id, width = 45)
y

# 替换原来的Y轴标签
p3 <- p2 + scale_y_discrete(labels = rev())
p3

# 自定义颜色渐变
p4 <- p3 + scale_colour_gradient(low = "red", high = "yellow")
p4

# 设置X轴范围，避免点的溢出绘图区
p5 <-
  p4 + scale_x_continuous(
    limits = c(-0.05, 0.35),
    breaks = c(0, 0.1, 0.2, 0.3),
    label = c("0", "0.1", "0.2", "0.3")
  )
p5

# 设置图例、坐标轴、图表的标题
p6 <-
  p5 + labs(
    size = "Gene number",
    color = "Qvalue",
    x = "Rich factor",
    y = "",
    title = "Top20 of GO enrichment"
  )
p6

# 自定义图表主题，对图表做精细调整
top.mar = 0.2
right.mar = 0.2
bottom.mar = 0.2
left.mar = 0.2
mytheme <-
  theme_bw() + theme(
    plot.title = element_text(size = rel(1), hjust = 0.5),
    axis.title = element_text(size = rel(1)),
    axis.title.y = element_text(size = rel(0.9)),
    panel.grid = element_blank(),
    legend.text = element_text(size = rel(0.6)),
    legend.title = element_text(size = rel(0.6)),
    plot.margin = unit(
      x = c(top.mar, right.mar, bottom.mar, left.mar),
      units = "inches"
    )
  )

# 最后应用自定义主题，效果如下
p7 <- p6 + mytheme
p7
