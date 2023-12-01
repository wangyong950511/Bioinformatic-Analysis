list<-read.csv("/list.csv")
head(list)


## 单独提出各个list
listA<-list$A
listB<-list$B
listC<-list$C
listD<-list$D
listE<-list$E
listF<-list$F
## 去除各自列表的重复值
listA<-listA[!duplicated(listA)]
listB<-listB[!duplicated(listB)]
listC<-listC[!duplicated(listC)]
listD<-listD[!duplicated(listD)]
listE<-listE[!duplicated(listE)]
listF<-listF[!duplicated(listF)]
## 去除各自列表的NA值
listA<-na.omit(listA)
listB<-na.omit(listB)
listC<-na.omit(listC)
listD<-na.omit(listD)
listE<-na.omit(listE)
listF<-na.omit(listF)
## 合并list
list<-list(A=listA,B=listB,C=listC,D=listD,E=listE,F=listF)

library(ggVennDiagram)
ggVennDiagram(list)

ggVennDiagram(
  x,
  category.names = names(x), #自定义数据集的名称
  show_intersect = FALSE, #是否交互式演示
  set_color = "black", # 数据集的颜色，默认即可
  set_size = NA, #数据集标签的大小，默认即可
  label = c("both", "count", "percent", "none"), # 显示数值、比例还是不显示
  label_alpha = 0.5, #各个标签集合的外框透明度，0为全透明，默认0.5
  label_geom = c("label", "text"), #显示数字和标签框
  label_color = "black",#标签颜色
  label_size = NA, #标签字体大小
  label_percent_digit = 0, #标签的百分位数
  label_txtWidth = 40, #标签文本的宽度
  edge_lty = "solid", #标签边缘类型，默认实心
  edge_size = 1, #标签边缘大小
  ...
)


ggVennDiagram(
  list,
  label = "count",
  edge_lty = "dashed",
  edge_size = 0.5
)

library(ggplot2)
ggVennDiagram(
  list,
  label = "count",
  edge_lty = "dashed",
  edge_size = 0.5)+
  scale_fill_gradient(low="steelblue",high = "brown")+
  theme_bw()+
  ggtitle("Veen plot for six sets")+
  theme(legend.position = 'none')


ggVennDiagram(
  list,
  label = "count",
  edge_lty = "dashed",
  edge_size = 0.5)+
  scale_fill_distiller(palette = "RdBu")+
  theme_void()+
  ggtitle("Veen plot for six sets")+
  theme(legend.position = 'none')

venn = Venn(list)
data = process_data(venn)

ggplot() +
  geom_sf(aes(fill = count), data = venn_region(data)) +
  geom_sf(aes(color = id), data = venn_setedge(data), show.legend = FALSE) +
  geom_sf_text(aes(label = name), data = venn_setlabel(data)) +
  geom_sf_label(aes(label = count), data = venn_region(data)) +
  theme_void()+
  theme(legend.position = 'none')+
  ggtitle("Veen plot for six sets")


ggplot() +
  geom_sf(aes(fill = id), data = venn_region(data)) +
  geom_sf(aes(color = id), data = venn_setedge(data), show.legend = FALSE) +
  geom_sf_text(aes(label = name), data = venn_setlabel(data)) +
  geom_sf_label(aes(label = count), data = venn_region(data)) +
  theme_void()+
  theme(legend.position = 'none')+
  ggtitle("Veen plot for six sets")
