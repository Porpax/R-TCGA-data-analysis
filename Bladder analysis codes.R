#   设置工作目录
setwd("D:/Study the R 22.5.14/R TCGA test/Bladder cancer TCGA/Bladder TCGA")

#批量读取
# stage 1: clean the data

rm(list = ls())     
library(dplyr)
library(readxl)
library(stringr)
library(survival)
library(survminer)
library(ggplot2)
library(survMisc)
library(tidyr)
library(stringi)


data_dir <- "D:/Study the R 22.5.14/R TCGA test/Bladder cancer TCGA/Bladder TCGA"
sig_dir <- "D:/Study the R 22.5.14/R TCGA test/Bladder cancer TCGA/Bladder TCGA"
star_counts_data_dir <- "D:/Study the R 22.5.14/R TCGA test/Bladder cancer TCGA/Bladder TCGA/Bladder"
#Read the sample sheet(read.table)
fgss1 <- read.table("gdc_sample_sheet.2023-03-20.tsv",header = T, sep = "\t",
                    check.names = F, stringsAsFactors = F) 
fgss <- fgss1[which(fgss1$'Sample Type' %in% "Primary Tumor"),]
colnames(fgss)[2] <- "file_name"
metadata <- jsonlite::fromJSON("metadata.cart.2023-03-20.json")
View(fgss) #fgss包含File ID	，File Name等8项对应关系

files <- list.files(star_counts_data_dir,recursive = T)
files <- files[! str_detect(files,"parcel")]
files <- files[str_detect(files,"rna_seq")]
matrix <- data.frame()
i=0
for (file in files) {
  i=i+1
  print(i)
  filename <- file %>% str_remove(".*/")
  NewName <- metadata %>% filter(file_name == filename)
  NewName <- (NewName$associated_entities[[1]])$entity_submitter_id # %>% substr(1,16)
  print(NewName)
  middata <- read.table(file.path(star_counts_data_dir,file),header = T,row.names = 1,sep = "\t",skip = 5)
  midmatrix <- middata %>% select(X.4)
  colnames(midmatrix) = NewName
  matrix <-  midmatrix %>% t.data.frame() %>% rbind(matrix)
}
save(matrix,file = "Bladder matrix.RData")

load(file="Bladder matrix.RData")

matrix <- matrix %>% t()
colnames(matrix) <- substr(colnames(matrix),1,16)##只取1-16位
a <- as.data.frame(gsub('[-]', '.', fgss$'Sample ID'))#把-替换成点
colnames(a) <- c("Sample ID")
fzt <- intersect(colnames(matrix),a$'Sample ID')#取交集
rownames(matrix) <- substr(rownames(matrix),1,17)#ENSG取小数点后一位便于匹配
matrix_pt <- matrix[,fzt]
matrix=as.data.frame(matrix)#$ operator is invalid for atomic vectors
#matrix_pt <- cbind(matrix$X,matrix_pt)
#colnames(matrix_pt)[1] <- "Ensembl ID"
dim(matrix) #检索matrix的维度
matrix[1:4,1:4]
#View(matrix)

###stage2：将ENSG转化为基因名
library(GenomeInfoDb)
library(dplyr)
library(stats)
library(tidyr)
library(base)
library(grDevices)
library(rtracklayer)


#：将ENSG转化为基因名
genome_annotation=read.delim( "gencode.gene.info.v22.tsv")#GDC官网下载Annotation file
id = genome_annotation[,1:2]
class(id)
id$gene_id <- substr(id$gene_id,1,17)#要根据此前ENSG取的位数一致
data = data.frame(matrix)
data$gene_id <- rownames(data)
data <- merge(data,id,by.x = "gene_id",by.y = "gene_id")
exp2 = data.frame(t(data[,-1]))#行列转置
#tail(exp2)[1:4,1:4]
colnames(exp2)=exp2[432,]  #取最后一行,并设为列名
exp_of_gene=exp2[-432,] #将最后一行删除
head(exp_of_gene)[1:4,1:4]
dim(exp_of_gene)
exp_bl_fpkm<- exp_of_gene[,c("CD274","CD8A","CD68")]
exp_bl_fpkm$tcga_id <- rownames(exp_bl_fpkm)
View(exp_bl_fpkm)
#write.csv(exp_bl_fpkm,file = "exp_bl_fpkm.csv")

###匹配生存信息
clinical = read.delim("clinical.tsv")#"clinical.cart.2023-01-15.tar"解压缩
#继续利用之前读的meta表格
meta <- jsonlite::fromJSON("metadata.cart.2023-03-20.json")
colnames(meta)
class(meta$associated_entities)
meta$associated_entities[[1]] #取第一列
meta$associated_entities[[1]]$case_id #取associated_entities表格中包含的case_id
case_id <-  sapply(meta$associated_entities,
                   function(x){x$case_id})
tcga_id <-  sapply(meta$associated_entities,
                   function(x){x$entity_submitter_id})
file_id = data.frame(case_id = case_id,tcga_id = tcga_id)
b <- gsub('[-]', '.', file_id$tcga_id)#把-替换成点，gsub("替换对象”，“替换成什么内容”，哪一列）
file_id = data.frame(file_id$case_id,b)
colnames(file_id) <- c("case_id","tcga_id")
head(file_id$case_id,2)
head(clinical$case_id,2)
table(clinical$case_id %in% file_id$case_id)
clinical <- merge(clinical,file_id,by =  "case_id")
new_clinical <- clinical %>% select(tcga_id, everything()) #select把tcga_id列选为第一列，后面是其他列everything
#everything()函数，所有变量
#save(new_clinical,file = "new_clinical without CA types.RData")
load(file ="new_clinical without CA types.RData" )

###匹配癌症类型（原发or转移）
sample_sheet <- read.delim("sample.tsv")#"biospecimen.cart.2023-01-14.tar"解压
colnames(sample_sheet)
#sample_sheet$sample_id[1]
sample_sheet$case_id[1]
samp_sheet_new <- sample_sheet[,c("case_id","sample_type")]
table(new_clinical$case_id%in% samp_sheet_new$case_id)
new_clinical <- merge(new_clinical,samp_sheet_new,by = "case_id")

###只要原发primary
new_clinical <- subset(new_clinical,new_clinical$sample_type=="Primary Tumor")
colnames(new_clinical)
#write.csv(new_clinical,file = "new_clinical_Lung.csv")

#形成一个新的只包含如下列的生存表格
survival <- new_clinical[,c("tcga_id",
                            "days_to_death",
                            "vital_status",
                            "days_to_last_follow_up",
                            "sample_type")]
new_surv<- survival %>% distinct(tcga_id, .keep_all = T) #去掉重复行

#将生存信息与表达矩阵结合
head(new_surv$tcga_id,2)
head(exp_bl_fpkm$tcga_id,2)
names(exp_bl_fpkm)
names(new_surv)
new_surv$tcga_id <- substr(new_surv$tcga_id,1,16)#tcga_id长度不一致
bl_fpkm_final = merge(exp_bl_fpkm,new_surv,by = "tcga_id")
View(bl_fpkm_final)
bl_fpkm_final$days_to_last_follow_up <- 
  as.numeric(bl_fpkm_final$days_to_last_follow_up)#先numeric之后，这列的非数字就`--会转化成NA缺失值
#ifelse函数，如果A列是缺失值，则取B列的值，否则还取A列
bl_fpkm_final$time <- ifelse(is.na(bl_fpkm_final$days_to_death), 
                             bl_fpkm_final$days_to_last_follow_up, 
                             bl_fpkm_final$days_to_death)
bl_fpkm_final$event=ifelse(bl_fpkm_final$vital_status=="Alive",0,1)
bl_fpkm_final$CD274 <- as.numeric(bl_fpkm_final$CD274)
bl_fpkm_final$CD8A <- as.numeric(bl_fpkm_final$CD8A)
bl_fpkm_final$CD68 <- as.numeric(bl_fpkm_final$CD68)
#如果有，把not reported数据行删掉
bl_fpkm_final = filter(bl_fpkm_final, bl_fpkm_final$vital_status != "Not Reported")
bl_fpkm_final$days_to_death <-  as.numeric(bl_fpkm_final$days_to_death)
bl_fpkm_final$days_to_last_follow_up <-  as.numeric(bl_fpkm_final$days_to_last_follow_up)
bl_fpkm_final$time <- ifelse(is.na(bl_fpkm_final$days_to_death), 
                             bl_fpkm_final$days_to_last_follow_up, 
                             bl_fpkm_final$days_to_death)

write.csv(bl_fpkm_final,"bl_fpkm.csv") #得到较为完整干净的数据

###KM生存分析
(.packages())
library(dplyr)
library(survival)
library(survminer)
library(survMisc)
library(ggplot2)

bl_fpkm <- read.csv("bl_fpkm.csv")
tmpCoxPh <- coxph(formula = Surv(time,event)~ CD274 ,data =  bl_fpkm)#以CD274作为单因素变量作Cox回归
cut274 <- cutp(tmpCoxPh)$CD274
# NOTE: cutp only consider the cutoffs instead of each patient.
cut274 <- as.data.frame(cut274)
#CD8A处理
tmpCD8A <- coxph(formula = Surv(time,event)~ CD8A ,data =  bl_fpkm)
cutCD8A <- cutp(tmpCD8A)$CD8A
cutCD8A <- as.data.frame(cutCD8A)
#取log值画图，数据均一化，
log_bl <- bl_fpkm
log_bl$logCD274 <- log2(log_bl$CD274+1)#对数转换的对象不能是负数或者0，但有的gene表达量为0，所以要对所有数据+1
log_bl$logCD8A <- log2(log_bl$CD8A+1)


##对取完log值的CD274作为变量作Cox回归，画图见下方
log_tmpCoxPh <- coxph(formula = Surv(time,event)~ logCD274 ,data =  log_bl)
log_cut274 <- cutp(log_tmpCoxPh)$logCD274
log_cut274 <- as.data.frame(log_cut274)
log_c274_gg <- ggplot(data = log_cut274,mapping=aes(x = logCD274,y = U)) #任何与数据向量顺序相关，需要逐个指定的参数都必须写在aes里
##对取完log值的CD8A作为变量作Cox回归，画图见下方
log_tmpCD8A <- coxph(formula = Surv(time,event)~ logCD8A ,data =  log_bl)
log_cutCD8A <- cutp(log_tmpCD8A)$logCD8A
log_cutCD8A <- as.data.frame(log_cutCD8A)
log_CD8A_gg <- ggplot(data = log_cutCD8A,mapping = aes(x = logCD8A,y = U))

#找出最大值,作log图中峰值的辅助线；绘制KM图；也作为下面散点图划分4个区的依据
#y[which(x==min(x))] 最小值
CD274max1 <- log_cut274$logCD274[between(log_cut274$logCD274, 0.5, 1)]
a1 <- CD274max1[which(log_cut274$U==max(log_cut274$U))] 
CD274max2 <- log_cut274$logCD274[between(log_cut274$logCD274, 0.1, 2.2)]
a2 <- CD274max2[which(log_cut274$U==max(log_cut274$U))]
#a3 <- median(log_bl$logCD274)
log_bl$rank_logCD274<- sort(log_bl$logCD274)
#k <- log_bl$rank_logCD274
#data25 <- (k<quantile(k, 0.25))#取logCD274的前25%
#data75 <- (k>quantile(k, 0.75))#后75%
b3 <- median(log_bl$logCD8A)
#CD274max1 <- log_cut274$logCD274[between(log_cut274$logCD274, 0, 1)] #检查一个数值是否在特定范围内。
b <- log_cutCD8A$logCD8A[which(log_cutCD8A$U==max(log_cutCD8A$U))]
survobj <- with(log_bl, Surv(time,event))#Surv()创建生存对象
fit0 <- survfit(survobj~1, data=log_bl) # 构建总体生存分析模型，survfit()对生存数据对象拟合生存函数，创建KM(Kaplan-Meier)生存曲线
summary(fit0) # 查看结果

###CD274 画log图并添加辅助线
png(file="BL_log_CD274_New.png",width=4,height=3,res = 300,units = 'in') #res分辨率
log_c274_gg + geom_step() +
  geom_smooth(span = 0.5, color = "red", alpha = 0.3, linetype = "dotdash") +
  labs(x = "CD274(log2 FPKM)", y = "Log-rank test score") +
  theme_classic()+ 
  geom_vline(aes(xintercept=a1), colour="gray60", linetype="dashed")+
  geom_segment(aes(x = a2, y = 0, xend = a2, yend = 17),linetype="dashed")
#geom_step()阶梯线图 ，geom_smooth()作拟合曲线，theme_classic() 去除自带的网格线背景    
dev.off()

###CD8A 画log图并添加辅助线
png(file="BL_log_CD8A_New.png",width=4,height=3,res = 300,units = 'in')
log_CD8A_gg + geom_step() +
  geom_smooth(span = 0.5, color = "blue", alpha = 0.3, linetype = "dotdash") +
  labs(x = "CD8A(log2 FPKM)", y = "Log-rank test score") +
  theme_classic() +
  geom_vline(aes(xintercept=b), colour="gray60", linetype="dashed")
dev.off()


log_bl$group2[log_bl$logCD274 <= a2 & log_bl$logCD8A >  b3] = "PD-L1 low"
log_bl$group2[log_bl$logCD274 > a2 & log_bl$logCD8A >  b3] = "PD-L1 high"

fit1 <- survfit(survobj~group2, data=log_bl)
png(file="Bladder_km_new_1.png",width=12,height=9,res = 300,units = 'in')
ggsurvplot(fit1,data = log_bl,pval = TRUE, conf.int = FALSE,#If TRUE,绘制生存率的95%置信区间
           title = "TCGA Bladder Primary Tumor",
           legend.title = "",#图例名称
           legend.labs = c("PD-L1 low","PD-L1 high"),#图例标注
           censor.shape="|", censor.size = 4,#指定删失点的形状；默认为"+"(3), 可选"|"(124)；censor.size 指定删失点形状的大小
           ggtheme = theme_bw(),#主题背景
           legend = "right",#将图例放在右边
           xlab = "Time",
           ylab = "Proportion of Overall Survival",
           #break.x.by = 100,
           break.y.by = 0.25,#设置y轴刻度间距
           risk.table = TRUE,#风险表添加
           risk.table.col = "strata",#风险表颜色跟随
           risk.table.height = 0.3,#生存表高度占据画幅百分比(区间0-1，1为只显示生存表)
           risk.table.y.text = TRUE,
           palette =  c("#4069A6","#D76364"))#显示风险表y轴标签 FALSE 是不显示
dev.off()

##保存每组之间差异性p_value
ptest_1 <- pairwise_survdiff(Surv(time, event)~group2, data = log_bl)
capture.output(ptest_1,file = "ptest_1.txt")
