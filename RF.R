#设置镜像，
local({r <- getOption("repos")  
r["CRAN"] <- "http://mirrors.tuna.tsinghua.edu.cn/CRAN/"   
options(repos=r)}) 

# 依赖包列表安装与加载：
package_list <- c("ggplot2","RColorBrewer","randomForest","caret", "pROC","dplyr","ggrepel")

# 判断R包加载是否成功来决定是否安装后再加载
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p,  warn.conflicts = FALSE)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}


##分类级选择
# 读取实验设计、和物种分类文件
metadata = read.table("metadata.txt",header = T, row.names = 1)
##样本分组（一共106个样本，74个Control，32个HFD，使用7：3的比例，即74：32）
# 假设从一百个样本中随机取70个，且无放回
idx = sample(1:9, size = 6, replace = F)
# 选择的样本标记为TRUE，未选择的为FALSE
idx = 1:9 %in% idx


##二型糖尿病人类样本分组（一共515个样本，138个Control，377个T2DM，使用7：3的比例，即360：155）
# 读取实验设计、和物种分类文件
metadata = read.table("metadata.txt",header = T, row.names = 1)
# 假设从一百个样本中随机取70个，且无放回
idx = sample(1:515, size =360 , replace = F)
# 选择的样本标记为TRUE，未选择的为FALSE
idx = 1:515 %in% idx

##一型糖尿病人类样本分组（发现集一共278个样本，149个Control，129个T1DM，使用7：3的比例，即195：83）
# 读取实验设计、和物种分类文件
metadata = read.table("metadata.txt",header = T, row.names = 1)
# 假设从一百个样本中随机取70个，且无放回
idx = sample(1:278, size =195 , replace = F)
# 选择的样本标记为TRUE，未选择的为FALSE
idx = 1:278 %in% idx

##一型糖尿病人类样本三个数据集分组（发现集一共166个样本，111个Control，55个T1DM，使用7：3的比例，即116：50）
# 读取实验设计、和物种分类文件
metadata = read.table("metadata.txt",header = T, row.names = 1)
# 假设从一百个样本中随机取70个，且无放回
idx = sample(1:166, size =116 , replace = F)
# 选择的样本标记为TRUE，未选择的为FALSE
idx = 1:166 %in% idx


##一型糖尿病人类样本Redistribution分组（发现集一共306个样本，188个Control，118个T1DM，使用7：3的比例，即214：92）
# 读取实验设计、和物种分类文件
metadata = read.table("metadata.txt",header = T, row.names = 1)
# 假设从一百个样本中随机取70个，且无放回
idx = sample(1:306, size =214 , replace = F)
# 选择的样本标记为TRUE，未选择的为FALSE
idx = 1:306 %in% idx


##一型糖尿病人类样本batch123分组（发现集一共166个样本，111个Control，55个T1DM，使用7：3的比例，即116：50）
# 读取实验设计、和物种分类文件
metadata = read.table("metadata.txt",header = T, row.names = 1)
# 假设从一百个样本中随机取70个，且无放回
idx = sample(1:166, size =116 , replace = F)
# 选择的样本标记为TRUE，未选择的为FALSE
idx = 1:166 %in% idx

# R4.0读取表不于默认为数据框
metadata$group = as.factor(metadata$group)

# 再用这个索引idx筛选对应的数据表，一部分作为训练集(train)，另一部分作为测试集(test)
metadata_train = metadata[idx,]
metadata_test = metadata[!idx,]


# metadata = subset(metadata, soiltype  %in% c("L"))
table =read.table("marker(0.3).txt",header = T, row.names = 1)

# 筛选共有的OTU
acoidx = rownames(metadata_train) %in% colnames(table)


#只留下二者共有的集合，构建OTU的训练集或测试集，比如有些Unassigned也去除掉了
otu_sub = table[, rownames(metadata_train)] 
dim(otu_sub)

# 物种分类文件，由usearch10 -sintax_summary生成，详见扩增子分析流程系列。但存在1对多时无法对应分类级颜色(如Unassigned可能属于多个门)，使用format2stamp.Rmd保留各级别名称
library(randomForest)
# "1Kingdom",界只有细菌、古菌类别太少；"7Species",扩增子中不太可信
for(i in c("2Phylum","3Order","4Order","5Family","6Family","8OTU0.1")){
  i="2Phylum"
  set.seed(0)
  table = read.table(paste0("tax_",i,".txt"),header = T, row.names = 1)
  #datatable = table[,rownames(metadata)]
  #使用训练集训练随机森林模型
  rf = randomForest(t(otu_sub), metadata_train$group, importance=T, proximity=T, ntree = 1000)
  print(i)
  print(rf)
  plot(rf)
}

#训练集
otu_tra = as.data.frame(t(otu_sub))
otu_tra$group = metadata[rownames(otu_tra),]$group
otu_tra$group
otu_tra$group = as.factor(otu_tra$group)

##测试集独立验证
metadata_test = metadata[!idx,]
metadata_test$group = as.factor(metadata_test$group)
otu_test = table[, rownames(metadata_test)] 
dim(otu_test)

# 转置，并添加分组信息
otutab_t = as.data.frame(t(otu_test))
otutab_t$group = metadata[rownames(otutab_t),]$group
otutab_t$group
otutab_t$group = as.factor(otutab_t$group)

#导入外部独立测试集
otu_val = read.table("val_C2.txt",header = T, row.names = 1, check.names = F)
otu_val$group = as.factor(otu_val$group)

#基于训练集随机森林模型验证
set.seed(315)
otutab.pred = predict(rf, t(otu_test) )  
pre_tab = table(observed=otutab_t[,"group"],
                predicted=otutab.pred) 
pre_tab
otutab.pred

val.pred = predict(rf,otu_val)
val.pred
#预测结果
confusionMatrix(otutab.pred, otutab_t$group)

##可视化验证结果
#1.绘制ROC曲线
pred1 <- predict(rf, newdata = otutab_t,type="prob")  #测试集
pred2 = predict(rf, newdata = otu_val,type="prob")  #外部独立验证集
pred3 = predict(rf, newdata = otu_tra,type = "prob")  #训练集
pred1
pred2

par(pin = c(1,1))


##外部独立验证组
roc.info<-roc(otu_val$group, 
              pred2[,1],  #提取随机森林模型中对应的预测指标
              plot=TRUE, 
              legacy.axes=TRUE, 
              percent=FALSE, 
              xlab="False positive percentage", 
              ylab="True postive percentage",
              col="black", 
              main = "Phylum ROC",
              lwd=4, #线条的粗细
              print.auc=TRUE)

##测试集
roc.info<-roc(otutab_t$group, 
              pred1[,1],  #提取随机森林模型中对应的预测指标
              plot=TRUE, 
              legacy.axes=TRUE, 
              percent=FALSE, 
              xlab="False positive percentage", 
              ylab="True postive percentage",
              col="black", 
              main = "Phylum ROC",
              lwd=4, #线条的粗细
              print.auc=TRUE)

##训练集
roc.info<-roc(otu_tra$group, 
              pred3[,1],  #提取随机森林模型中对应的预测指标
              plot=TRUE, 
              legacy.axes=TRUE, 
              percent=FALSE, 
              xlab="False positive percentage", 
              ylab="True postive percentage",
              col="black", 
              main = "Phylum ROC",
              lwd=4, #线条的粗细
              print.auc=TRUE)

library(pROC)
par(pin = c(1,1))
plot(roc.info,print.auc=TRUE,legacy.axes=TRUE)

  
              
##人类肠道微生物：训练集与测试集
g3 = ggroc(list(Training= roc3, Test = roc1),legacy.axes = TRUE)
g3 + theme_minimal() + ggtitle("Family ROC") + 
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="darkgrey", 
               linetype="dashed") +###设置对角线，主要设置起点和终点的X,Y坐标 
  coord_fixed() +  
  theme_bw() +
  theme(panel.grid =element_blank()) + #删除网格线
  theme(axis.text = element_text(size=12,colour = "black")) +   ## 刻度标签
  theme(axis.ticks = element_line(size=1)) +   ## 刻度线
  theme(panel.border = element_blank()) +   ## 删去外层边框
  theme(axis.line = element_line(size=1, colour = "black")) +   ## 再加上坐标轴
  geom_line(size = 1) + #线条的粗细
  theme(legend.position=c(0.75,0.20)) +   #设置图例位置
  theme(legend.title=element_blank()) + #删除图例标题
  scale_colour_manual(values = c("#EDCC86","#CD5555"),
                      breaks=c("Training", "Test"),
                      labels=c( "Training set, AUC=0.996 ","Test set, AUC=0.527")) + #修改标签名称和颜色
  theme(text=element_text(size=16,  Family="serif")) #修改字体为TNR


##人类T1DM肠道微生物：测试集与外部独立验证集
g3 = ggroc(list(Discovery = roc1, Validation= roc2),legacy.axes = TRUE)
g3 + theme_minimal() + ggtitle("Class ROC") + 
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="darkgrey", 
               linetype="dashed") +###设置对角线，主要设置起点和终点的X,Y坐标 
  coord_fixed() +  
  theme_bw() +
  theme(panel.grid =element_blank()) + #删除网格线
  theme(axis.text = element_text(size =12,colour = "black")) +   ## 刻度标签
  theme(axis.ticks = element_line(size =1)) +   ## 刻度线
  theme(panel.border = element_blank()) +   ## 删去外层边框
  theme(axis.line = element_line(size =1, colour = "black")) +   ## 再加上坐标轴
  geom_line(size  = 1) + #线条的粗细
  theme(legend.position=c(0.75,0.20)) +   #设置图例位置
  theme(legend.title=element_blank()) + #删除图例标题
  scale_colour_manual(values = c("#CD5555", "#27408B"),
                      breaks=c("Discovery", "Validation"),
                      labels=c( "Discovery set, AUC=0.825 ","Validation set, AUC=0.563")) + #修改标签名称和颜色
  theme(text=element_text(size =16,  family="serif")) #修改字体为TNR

# # 整理样本原始分组和预测分类
predict = data.frame(group = otutab_t[,"group"], predicted=otutab.pred)
# # 保存预测结果表
write.table("SampleID\t", file=paste("F_prediction_binary.txt",sep=""),append = F, quote = F, eol = "", row.names = F, col.names = F)
write.table(predict, file = "F_prediction_binary.txt",append = T, quote = F, row.names = T, col.names = T, sep = "\t")

#输出外部独立验证集的预测结果
vpredict = data.frame(group = otu_val[,"group"], predicted=val.pred)
write.table("SampleID\t", file=paste("F_valprediction_binary.txt",sep=""),append = F, quote = F, eol = "", row.names = F, col.names = F)
write.table(vpredict, file = "F_valprediction_binary.txt",append = T, quote = F, row.names = T, col.names = T, sep = "\t")


##输出重要性
imp= as.data.frame(rf$importance)
imp = imp[order(imp$MeanDecreaseAccuracy,decreasing = T),]
write.table(imp,file = "P-importance_feature.txt",quote = F,sep = '\t', row.names = T, col.names = T)  #输出重要性
head(imp)
varImpPlot(rf)

#交叉验证确定最佳预测微生物数量
set.seed(1) # 随机数据保证结果可重复，必须
# rfcv是随机森林交叉验证函数：Random Forest Cross Validation
metadata_train$group = as.factor(metadata_train$group)
result = rfcv(t(otu_sub), metadata_train$group, cv.fold=5)
# 查看错误率表，31时错误率最低，为最佳模型
result$error.cv
# 多次绘制
## 建立数据框保存多次结果
error.cv0 = data.frame(num = result$n.var, error.1 =  result$error.cv)
## 指定随机数循环5次
for (i in 1:(1+4)){
  print(i)
  set.seed(i)
  result= rfcv(t(otu_sub), metadata_train$group, cv.fold=5) #  scale = "log", step = 0.9
  error.cv0 = cbind(error.cv0, result$error.cv)
}
error.cv0 
# 绘制验证结果 
with(result, plot(n.var, error.cv0, log="x", type="o", lwd=2))
## 绘制交叉验证曲线

# 提取x轴标签
n.var = error.cv0$num
# 提取y轴数据+标签
error.cv = error.cv0[,2:6]
colnames(error.cv) = paste('err',1:5,sep='.')
# 添加均值
err.mean = apply(error.cv,1,mean)
# 合并新的数据库，x+error+mean
allerr = data.frame(num=n.var,err.mean=err.mean,error.cv)
# number of otus selected 人为在图中观察的结果，30几乎为最低，且数量可接受
optimal = 3


# 图1：机器学习结果交叉验证图，选择Top features
# 图中 + 5条灰色拆线+1条黑色均值拆线+一条最优垂线+X轴对数变换
write.table(allerr, file = "O-rfcv.txt", sep = "\t", quote = F, row.names = T, col.names = T)
p = ggplot() + # 开始绘图
  geom_line(aes(x = allerr$num, y = allerr$err.1), colour = 'grey') + # 5次验证灰线 
  geom_line(aes(x = allerr$num, y = allerr$err.2), colour = 'grey') + 
  geom_line(aes(x = allerr$num, y = allerr$err.3), colour = 'grey') + 
  geom_line(aes(x = allerr$num, y = allerr$err.4), colour = 'grey') + 
  geom_line(aes(x = allerr$num, y = allerr$err.5), colour = 'grey') + 
  geom_line(aes(x = allerr$num, y = allerr$err.mean), colour = 'black') + # 均值黑线
  geom_vline(xintercept = optimal, colour='black', lwd=0.36, linetype="dashed") + # 最优垂线
  coord_trans(x = "log2") + # X轴对数变换和刻度
  scale_x_continuous(breaks = c(1, 2, 5, 10, 20, 30, 50, 100, 200)) + # , max(allerr$num)
  labs(title=paste('Training set (n = ', dim(t(otu_sub))[1],')', sep = ''), 
       x='Number of Genus ', y='Cross-validation error rate') + 
  annotate("text", x = optimal, y = max(allerr$err.mean), label=paste("optimal = ", optimal, sep="")) + theme_bw()
p  
ggsave(p, file = "P-rfcv.pdf", width = 89, height = 59, unit = 'mm')

library(ggplot2)
p = ggplot() + # 开始绘图，不使用注释，未选择最优特征
  geom_line(aes(x = allerr$num, y = allerr$err.1), colour = 'grey') + # 5次验证灰线 
  geom_line(aes(x = allerr$num, y = allerr$err.2), colour = 'grey') + 
  geom_line(aes(x = allerr$num, y = allerr$err.3), colour = 'grey') + 
  geom_line(aes(x = allerr$num, y = allerr$err.4), colour = 'grey') + 
  geom_line(aes(x = allerr$num, y = allerr$err.5), colour = 'grey') + 
  geom_line(aes(x = allerr$num, y = allerr$err.mean), colour = 'black') + # 均值黑线
  scale_x_continuous(breaks = c(1, 2, 5, 10, 20, 30, 50, 100, 200)) + # , max(allerr$num)
  labs(title=paste('Training set (n = ', dim(t(otu_sub))[1],')', sep = ''), 
       x='Number of Order ', y='Cross-validation error rate') +  theme_bw()
p  


## 特征重要性可视化
## 预览和保存特征贡献度
imp= as.data.frame(rf$importance)
imp = imp[order(imp$MeanDecreaseAccuracy, decreasing = T),]
head(imp,n=optimal)
write.table(imp,file = "F-importance.txt",quote = F,sep = '\t', row.names = T, col.names = T)
# 简单可视化，比较丑
# varImpPlot(rf, main = "Feature importance",n.var = optimal, bg = par("bg"), color = par("fg"), gcolor = par("fg"), lcolor = "gray" )
# 图2. Feature重要性：绘制条形图+门属性着色
# 读取所有feature贡献度
imp = read.table("F21-importance_feature.txt", header=T, row.names= 1, sep="\t") 
# 分析选择top20分组效果最好，参数显示数量
imp = head(imp, n = 21)
imp = imp[order(imp$MeanDecreaseAccuracy, decreasing = F),]
# 简化全名，去掉界
imp$Genus = gsub("Bacteria\\|","",rownames(imp))
# 添加门用于着色(删除竖线后面全部)
imp$Phylum = gsub("\\|.*","",imp$Genus)
# 设置顺序
imp$Genus = factor(imp$Genus, levels = imp$Genus)
# 图2. 绘制物种类型种重要性柱状图
p = ggplot(imp, aes(x = Phylum, y = MeanDecreaseAccuracy, fill = Phylum)) +   
  geom_bar(stat = "identity") + 
  coord_flip() + theme_bw()
p
ggsave(paste("top14_feautre_O",".pdf", sep=""), p, width=89*2.5, height=59*2, unit='mm')
# 名称不一定唯一，需要手动修改
#  简化全名(只保留最后，有重名不可用，可选)
imp$Genus = gsub(".*\\|","",imp$Genus)
imp$Genus = factor(imp$Genus, levels = imp$Genus)
p = ggplot(imp, aes(x = Genus, y = MeanDecreaseAccuracy, fill = Phylum)) +   
  geom_bar(stat = "identity") + 
  coord_flip() + theme_bw()
p
ggsave(paste("G21_feautre",".pdf", sep=""), p, width=89*1.5, height=59*1.5, unit='mm')


write.table(CRC_meta,file = "CRC_meta.txt",quote = F,sep = '\t', row.names = T, col.names = T)
