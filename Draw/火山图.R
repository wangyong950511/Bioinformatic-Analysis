#载入相关R包
library(dplyr)
library(ggplot2)
library(ggrepel)
library(readr)

#读入数据-标准抬头：adj.P.Val（P值）；ID（标签值）；logFC（差异对数值
data <- read_tsv("/Users/wangyong/Downloads/GSE29796.top.table.tsv")
head(data)

#转成tibble便于后续使用，去掉不需要的列
dt <- as_tibble(data[c(-3,-4,-5,-8)])
head(dt)

#对Q值取对数
dt$log10FDR <- -log10(dt$adj.P.Val)

#生成显著上下调数据的分组标签
dt$group <- case_when(
  dt$logFC > 2 & dt$adj.P.Val < 0.05 ~ "Up",
  dt$logFC < -2 & dt$adj.P.Val < 0.05 ~ "Down",
  abs(dt$logFC) <= 2 ~ "None",
  dt$adj.P.Val >= 0.05 ~ "None"
)
head(dt)

#获取表达差异最显著的10个基因
top10sig <-
  filter(dt, group != "None") %>% distinct(Gene.symbol, .keep_all = T) %>% top_n(10, abs(logFC))
top10sig

#将差异表达Top10的基因表格拆分成up和down两部分
up <- filter(top10sig, group == "Up")
up
down <- filter(top10sig, group == "Down")
down

#新增一列，将Top10的差异基因标记成2，其他标记成1
dt$size <-
  case_when(!(dt$ID %in% top10sig$ID) ~ 1, dt$ID %in% top10sig$ID ~ 2)
head(dt)

#提取非Top10的基因表格
df <- filter(dt, size == 1)
head(df)

#指定绘图顺序，将group列转成因子型
df$group <-
  factor(df$group,
         levels = c("Up", "Down", "None"),
         ordered = T)

#开始绘图，建立映射
p0 <- ggplot(data = df, aes(logFC, log10FDR, color = group))

#添加散点
p1 <- p0 + geom_point(size = 1.6)
p1

#自定义半透明颜色（红绿）
mycolor <- c("#FF9999", "#99CC00", "gray80")
p21 <-
  p1 + scale_colour_manual(name = "", values = alpha(mycolor, 0.9))
p21

#其他配色方案
mycolor <- c("#FF99CC", "#99CC00", "gray80")
p22 <-
  p1 + scale_colour_manual(name = "", values = alpha(mycolor, 0.9))
p22

#继续添加Top10基因对应的点
p2 <-
  p22 + geom_point(
    data = up,
    aes(logFC, log10FDR),
    color = "#FF9999",
    size = 3,
    alpha = 0.9
  ) + geom_point(
    data = down,
    aes(logFC, log10FDR),
    color = "#7cae00",
    size = 3,
    alpha = 0.9
  )
p2

# #expansion函数设置坐标轴范围两端空白区域的大小，mult为倍数模式，add为加性模式
# p3 <-
#   p2 + labs(y = "-log10FDR") + scale_y_continuous(
#     expand = expansion(add = c(2, 0)),
#     limits = c(0, 40),
#     breaks = c(0, 10, 20, 30, 40),
#     label = c("0", "10", "20", "30", "40")
#   ) + scale_x_continuous(
#     limits = c(-4, 4),
#     breaks =  c(-4, -2, 0, 2, 4),
#     label = c("-4", "-2", "0", "2", "4")
#   )
# p3
p3 <- p2


#添加箭头
set.seed(007)
p4 <-
  p3 + geom_text_repel(
    data = top10sig,
    aes(logFC, log10FDR, label = Gene.symbol),
    force = 80,
    color = "grey20",
    size = 3,
    point.padding = 0.5,
    hjust = 0.5,
    arrow = arrow(
      length = unit(0.01, "npc"),
      type = "open",
      ends = "last"
    ),
    segment.color = "grey20",
    segment.size = 0.2,
    segment.alpha = 0.8,
    nudge_x = 0,
    nudge_y = 1
  )
p4

#自定义图表主题，对图表做精细调整
top.mar = 0.2
right.mar = 0.2
bottom.mar = 0.2
left.mar = 0.2

#隐藏纵轴，并对字体样式、坐标轴的粗细、颜色、刻度长度进行限定
mytheme <-
  theme_classic() + theme(
    text = element_text(
      family = "sans",
      colour = "gray30",
      size = 12
    ),
    axis.line = element_line(linewidth = 0.6, colour = "gray30"),
    axis.ticks = element_line(linewidth = 0.6, colour = "gray30"),
    axis.ticks.length = unit(1.5, units = "mm"),
    plot.margin = unit(
      x = c(top.mar, right.mar, bottom.mar, left.mar),
      units = "inches"
    )
  )

#应用自定义主题
p4 + mytheme

#添加辅助线
p5 <-
  p3 + geom_hline(
    yintercept = c(-log10(0.05)),
    size = 0.7,
    color = "orange",
    lty = "dashed"
  ) + geom_vline(
    xintercept = c(-1, 1),
    size = 0.7,
    color = "orange",
    lty = "dashed"
  )
p5

#添加其他样式的标签，同时为了方便自定义左右区域的标签，这里使用up、down两个独立的子表格
p6 <- p5 + geom_label_repel(
  data = up,
  aes(logFC, log10FDR, label = Gene.symbol),
  nudge_x = 1,
  nudge_y = 5,
  color = "white",
  alpha = 0.9,
  point.padding = 0.5,
  size = 3,
  fill = "#96C93D",
  segment.size = 0.5,
  segment.color = "grey50",
  direction = "y",
  hjust = 0.5
) + geom_label_repel(
  data = down,
  aes(logFC, log10FDR, label = Gene.symbol),
  nudge_x = -1,
  nudge_y = 3,
  color = "white",
  alpha = 0.9,
  point.padding = 0.5,
  size = 3,
  fill = "#9881F5",
  segment.size = 0.5,
  segment.color = "grey50",
  direction = "y",
  hjust = 0.5
)

#应用自定义主题
p7 <- p6 + mytheme
p7