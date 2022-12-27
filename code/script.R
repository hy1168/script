library(reshape2)
library(ggpubr)
library(maftools)
library(tidyr)
library(pheatmap)
library(clusterProfiler)
library(ComplexHeatmap)
library(GSVA)
library(GSEABase)
options(stringsAsFactors = F)
source('base.R')


load('script/tcga_cli.RData')
head(tcga_cli)

load('script/tcga_exp.RData')
tcga_tpm_log_T[1:5,1:5]

hypoxia.gene=read.gmt('HALLMARK_HYPOXIA.v7.5.1.gmt')
hypoxia.gene=hypoxia.gene$gene
length(hypoxia.gene)

#GEO####
load('GSE26712_dat.RData')
head(GSE26712_cli)
GSE26712_exp[1:5,1:5]



#01.基因的鉴定+聚类####
dir.create('reaults/01.Cluster')
tcga.hypoxia.score=t(ssGSEAScore_by_genes(gene.exp = tcga_tpm_log_T,genes = hypoxia.gene))
head(tcga.hypoxia.score)
tcga.hypoxia.score=data.frame(Samples=rownames(tcga.hypoxia.score),score=tcga.hypoxia.score[,1])
#预后相关hypoxia
tcga_hypoxia_cox=cox_batch(tcga_tpm_log_T[intersect(hypoxia.gene,rownames(tcga_tpm_log_T)),],
                           tcga_cli$OS.time,
                           tcga_cli$OS)
tcga_hypoxia_cox=tcga_hypoxia_cox %>% drop_na()
#writeMatrix(tcga_hypoxia_cox,row=T,header=T,outpath = 'results/01.Cluster/tcga_hypoxia_cox.txt')
tcga_hypoxia_cox_fit=rownames(tcga_hypoxia_cox[tcga_hypoxia_cox$p.value<0.05,])
length(tcga_hypoxia_cox_fit)
consen_gene=tcga_hypoxia_cox_fit
length(consen_gene)#14


library(ConsensusClusterPlus)
clusterAlg_name=c('hc','pam','km','kmdist')[1]
distance_name=c('pearson','spearman','euclidean','binary','maximum','canberra','minkowski')[2]
#TCGA#####
tcga_consen_data=as.matrix(tcga_tpm_log_T[consen_gene,])
#tcga_consen_data=t(scale(t(tcga_consen_data),scale = F))
tcga_consen_data=t(scale(t(tcga_consen_data),scale = T))
#tcga_consen_data=sweep(tcga_consen_data,1,apply(tcga_consen_data, 1, mean))
#tcga_consen_data=sweep(tcga_consen_data,1,apply(tcga_consen_data, 1, median))
#tcga_consen_data=as.dist(1-cor(tcga_consen_data,method = 'pearson'))
tcga_consen_data=as.matrix(tcga_consen_data)
dim(tcga_consen_data)
tcga_clust_subtype <- ConsensusClusterPlus(tcga_consen_data
                                           , maxK = 10, reps = 500, pItem = 0.8
                                           , pFeature = 1
                                           , title = "TCGA_subtype"
                                           , clusterAlg = clusterAlg_name
                                           , distance = distance_name
                                           , plot = "pdf"
                                           , writeTable = T
                                           , seed = 123456)
k=2
tcga.subtype <- data.frame(Samples = names(tcga_clust_subtype[[k]]$consensusClass),
                           Cluster=tcga_clust_subtype[[k]]$consensusClass)
tcga.subtype$Cluster=paste0('C',tcga.subtype$Cluster)
tcga.subtype$Cluster=gsub('C2','S1',tcga.subtype$Cluster)
tcga.subtype$Cluster=gsub('C1','S2',tcga.subtype$Cluster)
tcga.subtype$Cluster=gsub('S','C',tcga.subtype$Cluster)
table(tcga.subtype$Cluster)
writeMatrix(tcga.subtype,row=T,header=T,outpath = 'results/01.Cluster/tcga.subtype.txt')
#tcga.subtype=readMatrix(inpath = 'results/01.Cluster/tcga.subtype.txt',row=T,header=T)
fig1a=ggplotKMCox(data.frame(time = tcga_cli[rownames(tcga.subtype),]$OS.time/365
                             , event = tcga_cli[rownames(tcga.subtype),]$OS
                             , tcga.subtype$Cluster)
                  ,add_text = '',title = 'TCGA'
                  ,ggsci::pal_nejm()(10)[1:2]	
                  ,labs = c('C1','C2')
)
fig1a



#GSE26712(OK)#####
GSE26712_consen_data=as.matrix(GSE26712_exp[intersect(consen_gene,rownames(GSE26712_exp)),])
#GSE26712_consen_data=t(scale(t(GSE26712_consen_data),scale = F))
GSE26712_consen_data=t(scale(t(GSE26712_consen_data),scale = T))
#GSE26712_consen_data=sweep(GSE26712_consen_data,1,apply(GSE26712_consen_data, 1, mean))
#GSE26712_consen_data=sweep(GSE26712_consen_data,1,apply(GSE26712_consen_data, 1, median))
#GSE26712_consen_data=as.dist(1-cor(GSE26712_consen_data,method = 'pearson'))
GSE26712_consen_data=as.matrix(GSE26712_consen_data)
dim(GSE26712_consen_data)
GSE26712_clust_subtype <- ConsensusClusterPlus(GSE26712_consen_data
                                               , maxK = 10, reps = 500, pItem = 0.8
                                               , pFeature = 1
                                               , title = "GSE26712_subtype"
                                               , clusterAlg = clusterAlg_name
                                               , distance = distance_name
                                               , plot = "pdf"
                                               , writeTable = T
                                               , seed = 123456)
k=2
GSE26712.subtype <- data.frame(Samples = names(GSE26712_clust_subtype[[k]]$consensusClass),
                               Cluster=GSE26712_clust_subtype[[k]]$consensusClass)
GSE26712.subtype$Cluster=paste0('C',GSE26712.subtype$Cluster)
GSE26712.subtype$Cluster=gsub('C2','S1',GSE26712.subtype$Cluster)
GSE26712.subtype$Cluster=gsub('C1','S2',GSE26712.subtype$Cluster)
GSE26712.subtype$Cluster=gsub('S','C',GSE26712.subtype$Cluster)
table(GSE26712.subtype$Cluster)
writeMatrix(GSE26712.subtype,row=T,header=T,outpath = 'GSE26712.subtype.txt')
fig1b=ggplotKMCox(data.frame(time = GSE26712_cli$OS.time/365
                             , event = GSE26712_cli$OS
                             , GSE26712.subtype$Cluster)
                  ,add_text = '',title = 'GSE26712'
                  ,ggsci::pal_nejm()(10)[1:4]
                  ,labs = c('Clust1','Clust2'))
fig1b


tcga.score.sub=merge(tcga.hypoxia.score,tcga.subtype,by='Samples')
fig1c=mg_violin_1(tcga.score.sub[, c("Cluster", "score")]
                  ,melt = T,group_col = ggsci::pal_nejm()(10)
                  , legend.pos = 'tl'
                  ,ylab = 'Hypoxia Score')
fig1c

#GSE26712
GSE26712.hypoxia.score=t(ssGSEAScore_by_genes(gene.exp = GSE26712_exp,genes = hypoxia.gene))
head(GSE26712.hypoxia.score)
GSE26712.hypoxia.score=data.frame(Samples=rownames(GSE26712.hypoxia.score),score=GSE26712.hypoxia.score[,1])
GSE26712.score.sub=merge(GSE26712.hypoxia.score,GSE26712.subtype,by='Samples')
fig1d=mg_violin_1(GSE26712.score.sub[, c("Cluster", "score")]
                  ,melt = T,group_col = ggsci::pal_nejm()(10)
                  , legend.pos = 'tl'
                  ,ylab = 'Hypoxia Score')
fig1d

#02.亚型临床信息特征#############
tcga.subtype.cli=merge(tcga_cli,tcga.subtype,by='Samples')
tcga.subtype.cli1=tcga.subtype.cli
table(tcga.subtype.cli1$Stage)
tcga.subtype.cli1$Stage[tcga.subtype.cli1$Stage == 'I' | tcga.subtype.cli1$Stage == 'II'] <- 'I+II'
tcga.subtype.cli1$Stage[tcga.subtype.cli1$Stage == 'III' | tcga.subtype.cli1$Stage == 'IV'] <- 'III+IV'
table(tcga.subtype.cli1$Grade)
tcga.subtype.cli1$Grade[tcga.subtype.cli1$Grade == 'G1' | tcga.subtype.cli1$Grade == 'G2'] <- 'G1+G2'
tcga.subtype.cli1$Grade[tcga.subtype.cli1$Grade == 'G3' | tcga.subtype.cli1$Grade == 'G4'] <- 'G3+G4'
fig2=list()
fig2[[1]]=plotMutiBar(dat = table(tcga.subtype.cli$Cluster,tcga.subtype.cli$Stage),
                      ist = T,legTitle = 'Stage')
fig2[[2]]=plotMutiBar(dat = table(tcga.subtype.cli1$Cluster,tcga.subtype.cli1$Grade),
                      ist = T,legTitle = 'Grade')
fig2[[3]]=plotMutiBar(dat = table(tcga.subtype.cli1$Cluster,tcga.subtype.cli1$Age1),
                      ist = T,legTitle = 'Age')
fig2[[4]]=plotMutiBar(dat = table(tcga.subtype.cli1$Cluster,tcga.subtype.cli1$Status),
                      ist = T,legTitle = 'Status')
fig2=mg_merge_plot(fig2,ncol=4,nrow = 1)
fig2


#03.亚型免疫特征#############
#dir.create('results/03.subtype.immu')
#TCGA免疫####
###22 种免疫细胞
load('script/TCGA_TME.RData')

mg_PlotMutiBoxplot(tcga.cibersort[tcga.subtype$Samples,1:22]
                   , group = tcga.subtype$Cluster
                   , legend.pos = 'top'
                   #, test_method = 'anova'
                   , test_method = 't.test'
                   , add = 'boxplot'
                   , ylab = 'Score'
                   , group_cols = ggsci::pal_nejm()(10)[1:2])	


fig3a=mg_PlotMutiBoxplot(tcga.est[tcga.subtype$Samples,]
                         , group = tcga.subtype$Cluster
                         , legend.pos = 'top'
                         #, test_method = 'anova'
                         , test_method = 't.test'
                         , add = 'boxplot'
                         , ylab = 'Score'
                         , group_cols = ggsci::pal_nejm()(10)[1:2])
fig3a


#MCP Count评分
fig3b=mg_PlotMutiBoxplot(tcga.immu.mcp[tcga.subtype$Samples,]
                         , group = tcga.subtype$Cluster
                         , legend.pos = 'top'
                         #, test_method = 'anova'
                         , test_method = 't.test'
                         , add = 'boxplot'
                         , ylab = 'Score'
                         , group_cols = ggsci::pal_nejm()(10)[1:2])
fig3b



#ssGSEA评分
fig3c=mg_PlotMutiBoxplot(tcga.immu.ssgsea[tcga.subtype$Samples,]
                         , group = tcga.subtype$Cluster
                         #, legend.pos = 'top'
                         #, test_method = 'anova'
                         , test_method = 't.test'
                         , add = 'boxplot'
                         , ylab = 'Score'
                         , group_cols = ggsci::pal_nejm()(10)[1:2])
fig3c
#免疫检查点相关基因
fig3d=mg_PlotMutiBoxplot(tcga.icg.exp[tcga.subtype$Samples,]
                         , group = tcga.subtype$Cluster
                         #, legend.pos = 'top'
                         # , test_method = 'anova'
                         , test_method = 't.test'
                         , add = 'boxplot'
                         , ylab = 'Gene Expression'
                         , group_cols = ggsci::pal_nejm()(10)[1:2])
fig3d


###########TIDE评分
tcga_tide_res_subtype=cbind.data.frame(tcga_tide_res[tcga.subtype$Samples,],
                                       Cluster=tcga.subtype$Cluster)
tide_sel=c('TIDE','IFNG','MDSC','Exclusion','Dysfunction')

fig3_tide=list()
for (i in 1:5) {
  fig3_tide[[i]]=sig_boxplot(tcga_tide_res_subtype[, c("Cluster", tide_sel[i])],
                             leg = 'Cluster',ylab = tide_sel[i],
                             palette = ggsci::pal_nejm()(10))
}
fig3_tide=mg_merge_plot(fig3_tide,ncol=5,nrow=1,common.legend = T)
fig3_tide

#04.肿瘤相关通路差异+通路####
dir.create('results/04.subtype.path')
head(tcga.onco)
load('script/TCGA_pathway.RData')
fig4a=mg_PlotMutiBoxplot(tcga.onco[tcga.subtype$Samples,]
                         , group = tcga.subtype$Cluster
                         , legend.pos = 'top'
                         , test_method = 't.test'
                         , add = 'boxplot'
                         , ylab = 'Score'
                         , group_cols = ggsci::pal_nejm()(10))	
fig4a


############GSEA
head(tcga_GSEA$EnrichTable)
inds1=which(tcga_GSEA$EnrichTable$NP<0.05 & tcga_GSEA$EnrichTable$NES>0)
tcga_GSEA$EnrichTable[inds1,]
tcag.plot1<-plot_GSEA_By_nodes(tcga_GSEA,indexs=c(134,136,141,157,164,165,172,181,182))
tcag.plot1


######05.差异基因#####
dir.create('results/06.DEGs')
C12_degs<-limma_DEG(tcga_tpm_log_T,
                       group=as.character(tcga.subtype$Cluster),'C1','C2')
C12_degs$Summary
deg.p.value=0.05
deg.fc=log2(1.2)
diff.genes=rownames(C12_degs$DEG[which(C12_degs$DEG$P.Value<deg.p.value & abs(C12_degs$DEG$logFC) > deg.fc),])
length(diff.genes)#4880

fig6a=mg_volcano(logfc = C12_degs$DEG$logFC,
                 pvalue = C12_degs$DEG$adj.P.Val,
                 cutFC = deg.fc,cutPvalue = deg.p.value
                 ,legend.pos ='tl'
                 ,leg = ' TCGA C1 vs C2')


########GSE26712差异分析
C12_degs_GSE26712<-limma_DEG(GSE26712_exp,
                                group=as.character(GSE26712.subtype$Cluster),'C1','C2')
C12_degs_GSE26712$Summary
GSE26712.diff.genes=rownames(C12_degs_GSE26712$DEG[which(C12_degs_GSE26712$DEG$P.Value<deg.p.value & 
                                                           abs(C12_degs_GSE26712$DEG$logFC) > deg.fc),])
length(GSE26712.diff.genes)
fig6b=mg_volcano(logfc = C12_degs_GSE26712$DEG$logFC,
                 pvalue = C12_degs_GSE26712$DEG$adj.P.Val,
                 cutFC = deg.fc,cutPvalue = deg.p.value
                 ,legend.pos ='tl'
                 ,leg = ' GSE26712 C1 vs C2')


com.gene=intersect(GSE26712.diff.genes,diff.genes)
length(com.gene)



#07.构建风险模型####
#TCGA####
tcga_model_data <- t(tcga_tpm_log_T[com.gene,])
colnames(tcga_model_data)=gsub('-','__',colnames(tcga_model_data))
tcga_model_data=merge(data.frame(Samples=tcga_cli$Samples,OS.time=tcga_cli$OS.time,OS=tcga_cli$OS),
                      data.frame(Samples=rownames(tcga_model_data),tcga_model_data),
                      by='Samples')
rownames(tcga_model_data)=tcga_model_data$Samples
tcga_model_data=tcga_model_data[,-1]
tcga_model_data=crbind2DataFrame(tcga_model_data)
dim(tcga_model_data)


##单因素
tcga.cox=cox_batch(tcga_tpm_log_T[com.gene,],
                   tcga_cli$OS.time,
                   tcga_cli$OS)
p_value_cutoff=0.05
tcga.cox <- na.omit(tcga.cox)
tcga.cox=as.data.frame(tcga.cox)
tcga.cox$coef=log(tcga.cox$HR)
tcga.cox$type=rep('None',nrow(tcga.cox))
tcga.cox$type[which(tcga.cox$p.value<p_value_cutoff & tcga.cox$coef>0)]='Risk'
tcga.cox$type[which(tcga.cox$p.value<p_value_cutoff & tcga.cox$coef<0)]='Protective'
writeMatrix(tcga.cox,'DEGs.cox.txt')
tcga.cox.fit <- rownames(tcga.cox[tcga.cox$p.value<p_value_cutoff,])
length(tcga.cox.fit)#59
table(tcga.cox$type)
fig7a=ggplot(data = tcga.cox,
             aes(x = coef,
                 y = -log10(p.value)))+
  geom_point(alpha=0.4, size=3.5, aes(color=type))+
  scale_color_manual(values=c(mg_colors[2],'grey',mg_colors[1]),
                     limits = c("Protective",'None', "Risk"),name='State')+
  geom_hline(yintercept = -log10(p_value_cutoff),lty=4,col="black",lwd=0.8)+
  ylab('-log10(pvalue)')+xlab('Cox coefficient')+
  theme_bw()+
  theme(axis.text.y=element_text(family="Times",face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Times",face="plain"), #设置y轴标题的字体属性
        legend.text=element_text(face="plain", family="Times", colour="black"),  #设置图例的子标题的字体属性
        legend.title=element_text(face="plain", family="Times", colour="black" ),#设置图例的总标题的字体属性
        legend.justification=c(1,0), legend.position=c(1,0),
        legend.background = element_rect(fill = NA, colour = NA))
fig7a
tcga.cox.fit=gsub('-','__',tcga.cox.fit)


tcga.lasso<-get_riskscore.lasso(dat = tcga_model_data[,tcga.cox.fit],
                                os=tcga_model_data$OS,
                                os.time = tcga_model_data$OS.time
                                ,labels=c('B','C')
)
tcga.lasso$plot
length(tcga.lasso$lasso.gene)#23
tcga.module.risk=get_riskscore(dat = tcga_model_data[,tcga.lasso$lasso.gene],
                               #dat = tcga_model_data[,tcga.cox.fit],
                               os=tcga_model_data$OS,
                               os.time = tcga_model_data$OS.time,
                               step=T,
                               direction=c("both", "backward", "forward")[1])
length(tcga.module.risk$module.gene)#7
tcga.module.risk$model

tcga.cox$gene=rownames(tcga.cox)
fig7d=tcga.cox[tcga.module.risk$module.gene,] %>% 
  ggplot(aes(reorder(gene, coef), coef)) +
  geom_col(aes(fill = type)) +
  coord_flip() +
  scale_fill_manual(values=ggsci::pal_lancet('lanonc',alpha =0.6)(9)[c(7,1)]) +
  coord_flip() +
  labs(x = "") +
  labs(y = "LASSO Cox coefficient") +
  theme(axis.text.y = element_text(angle = 0, hjust = 1),legend.position="top")
fig7d


fig8a=ggplotTimeROC(time = tcga.module.risk$result$time/365
                    ,status = tcga.module.risk$result$status
                    ,score = tcga.module.risk$result$riskscore
                    ,mks = c(1,2,3,4,5))
fig8a
fig8b=ggplotKMCox(dat = data.frame(tcga.module.risk$result$time/365,
                                   tcga.module.risk$result$status,
                                   ifelse(tcga.module.risk$result$riskscorez>=0,'High','Low'))
                  ,title = 'TCGA',add_text = ''
                  ,palette = ggsci::pal_npg("nrc")(10)[4:5],labs = c('High','Low'))
fig8b
#GSE26712####
GSE26712_model_data <- t(GSE26712_exp[intersect(rownames(GSE26712_exp),com.gene),])
GSE26712_model_data=merge(data.frame(Samples=GSE26712_cli$Samples,OS.time=GSE26712_cli$OS.time,OS=GSE26712_cli$OS),
                          data.frame(Samples=rownames(GSE26712_model_data),GSE26712_model_data),
                          by='Samples')
rownames(GSE26712_model_data)=GSE26712_model_data$Samples
GSE26712_model_data[1:5,1:5]
GSE26712_model_data=GSE26712_model_data[,-1]
GSE26712_model_data=crbind2DataFrame(GSE26712_model_data)
dim(GSE26712_model_data)

GSE26712.module.risk=get_riskscore(dat = GSE26712_model_data[,intersect(tcga.module.risk$module.gene,
                                                                        rownames(GSE26712_exp))],
                                   os=GSE26712_model_data$OS,
                                   os.time = GSE26712_model_data$OS.time,
                                   step=F,
                                   direction=c("both", "backward", "forward")[1])
length(GSE26712.module.risk$module.gene)

fig8c=ggplotTimeROC(GSE26712_model_data$OS.time,
                    GSE26712_model_data$OS,
                    GSE26712.module.risk$result$riskscorez,mks = c(1,2,3,4,5))
fig8d=ggplotKMCox(data.frame(GSE26712_model_data$OS.time/365,
                             GSE26712_model_data$OS,
                             ifelse(GSE26712.module.risk$result$riskscorez>=0,'High','Low'))
                  ,add_text = ''
                  ,title = 'GSE26712',ggsci::pal_npg("nrc")(10)[4:6]
                  ,labs = c('High','Low'))


GSE26712.risktype.cli=data.frame(GSE26712_cli,
                                 Risktype=ifelse(GSE26712.module.risk$result$riskscorez>=0,'High','Low'),
                                 Riskscore=GSE26712.module.risk$result$riskscore)
table(GSE26712.risktype.cli$Risktype)
rownames(GSE26712.risktype.cli)=GSE26712.risktype.cli$Samples


##########IMV210###########
custom_theme <- function() {
  theme_survminer() %+replace%
    theme(
      plot.title=element_text(hjust=0.5)
    )
}
load('script/IMV210_dat.RData')
head(IMvigor210_cli)
IMvigor210_exp[1:5,1:5]

IMvigor210_model_data <- cbind(IMvigor210_cli[, c("OS.time", "OS")],
                               t(IMvigor210_exp[, rownames(IMvigor210_cli)]))

IMvigor210.genes <- intersect(tcga.module.risk$module.gene, colnames(IMvigor210_model_data))
IMvigor210.genes

library(survminer)
IMvigor210_cli1 <- IMvigor210_cli[IMvigor210_cli$Response != 'NE', ]
colnames(IMvigor210_cli1)
IMvigor210_model_data <- IMvigor210_model_data[rownames(IMvigor210_cli1), ]
IMvigor210_model_data<-IMvigor210_model_data[,c("OS.time", "OS",IMvigor210.genes)]
#IMvigor210.genes=gsub('-','__',IMvigor210.genes)
#colnames(IMvigor210_model_data)=gsub('-','__',colnames(IMvigor210_model_data))
fmla.IMvigor210 <- as.formula(paste0("Surv(OS.time, OS) ~"
                                     ,paste0(IMvigor210.genes,collapse = '+')))

cox.IMvigor210 <- coxph(fmla.IMvigor210, data =as.data.frame(IMvigor210_model_data))
IMvigor210_lan <- coef(cox.IMvigor210)

risk.IMvigor210=as.numeric(IMvigor210_lan%*%as.matrix(t(IMvigor210_model_data[,names(IMvigor210_lan)])))

IMvigor210_model_data$RS <- risk.IMvigor210
# IMvigor210.data.point <- surv_cutpoint(IMvigor210_model_data, time = "OS.time", event = "OS",
#                                        variables = 'RS')
# IMvigor210.cutoff <- as.numeric(summary(IMvigor210.data.point)[1])
# IMvigor210.cutoff


IMvigor210.roc <- ggplotTimeROC(IMvigor210_model_data$OS.time / 365,
                                IMvigor210_model_data$OS,
                                risk.IMvigor210,
                                mks = c(0.5,1,1.5))
IMvigor210.km <- ggplotKMCox(data.frame(IMvigor210_model_data$OS.time / 365,
                                        IMvigor210_model_data$OS,
                                        ifelse(risk.IMvigor210>=median(risk.IMvigor210),'High','Low')),
                             title = 'IMvigor210',add_text = '',
                             labs = c('High','Low'))
IMvigor210.km


str(IMvigor210_tide_dat)
# IMvigor210.tide.point <- surv_cutpoint(IMvigor210_tide_dat, time = "OS.time", event = "OS",
#                                        variables = 'TIDE')
# IMvigor210.tide.cutoff <- as.numeric(summary(IMvigor210.tide.point)[1])
# IMvigor210.tide.cutoff#1.47
IMvigor210_tide_dat$TIDE.group <- ifelse(IMvigor210_tide_dat$TIDE >= median(IMvigor210_tide_dat$TIDE), 'High', 'Low')


IMvigor210.tide.roc <- ggplotTimeROC(IMvigor210_tide_dat$OS.time / 365,
                                     IMvigor210_tide_dat$OS,
                                     IMvigor210_tide_dat$TIDE,
                                     mks = c(0.5,1,1.5))
IMvigor210.tide.km <- ggplotKMCox(data.frame(IMvigor210_tide_dat$OS.time / 365,
                                             IMvigor210_tide_dat$OS,
                                             IMvigor210_tide_dat$TIDE.group),
                                  add_text = '',
                                  title = 'TIDE',
                                  labs = c('High','Low'))
IMvigor210.tide.km

head(IMvigor210_tide_dat[, c("Response", "TIDE", "IRS")])
table(IMvigor210_tide_dat$Response)
IMvigor210_tide_dat$Response[IMvigor210_tide_dat$Response == 'PD'] <- 0
IMvigor210_tide_dat$Response[IMvigor210_tide_dat$Response == 'SD'] <- 0
IMvigor210_tide_dat$Response[IMvigor210_tide_dat$Response == 'CR'] <- 1
IMvigor210_tide_dat$Response[IMvigor210_tide_dat$Response == 'PR'] <- 1

library(pROC)
library(dplyr)
IMvigor210.tide_roc_dat <- data.frame()
for (ty in c("TIDE", "IRS")) {
  print(ty)
  ty_dat <- IMvigor210_tide_dat[, c(ty, "Response")]
  # ty_dat <- ty_dat[ty_dat[, ty] == '0' | ty_dat[, ty] == '1', ]
  ty_dat <- crbind2DataFrame(ty_dat)
  
  ty_roc <-  roc(as.formula(paste0('Response ~ ', ty)), data = ty_dat, smooth = T)
  tmp_roc_dat <- data.frame(TP = ty_roc$sensitivities,
                            FP = 1 - ty_roc$specificities,
                            Methods =paste0(ty, ' AUC = ', round(ty_roc$auc, 2)))
  IMvigor210.tide_roc_dat <- rbind(IMvigor210.tide_roc_dat, tmp_roc_dat)
}
head(IMvigor210.tide_roc_dat)

mycolor <- ggsci::pal_jama(alpha =1)(9)
IMvigor210.tide_roc_dat$Methods=gsub('IRS','HYRS',IMvigor210.tide_roc_dat$Methods)

IMvigor210.Response_roc<-ggplot(IMvigor210.tide_roc_dat, aes(x=FP,y=TP, fill=Methods))+
  geom_line(aes(colour=Methods),lwd=0.75)+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
  xlab('False positive rate')+ylab('True positive rate')+ 
  scale_colour_manual(values = mycolor[2:3]) +
  labs(title = 'Response') +
  theme_bw() +
  theme(legend.position=c(1,0),legend.justification=c(1,0),
        legend.background = element_rect(fill = NA, colour = NA))
IMvigor210.Response_roc



# GSE91061 ####
load('script/GSE91061_dat.RData')
head(GSE91061_cli1)
GSE91061_exp[1:5,1:5]

GSE91061_model_data <- cbind(GSE91061_cli1[, c("OS.time", 'OS')],
                             t(GSE91061_exp[, rownames(GSE91061_cli1)]))
GSE91061.genes <- intersect(tcga.module.risk$module.gene, colnames(GSE91061_model_data))
GSE91061.genes
GSE91061_model_data=GSE91061_model_data[,c("OS.time", 'OS',GSE91061.genes)]
#
#GSE91061.genes=gsub('-','__',GSE91061.genes)
#colnames(GSE91061_model_data)=gsub('-','__',colnames(GSE91061_model_data))

fmla.GSE91061 <- as.formula(paste0("Surv(OS.time, OS) ~"
                                   ,paste0(GSE91061.genes,collapse = '+')))
cox.GSE91061 <- coxph(fmla.GSE91061, data =as.data.frame(GSE91061_model_data))
GSE91061_lan <- coef(cox.GSE91061)

risk.GSE91061=as.numeric(GSE91061_lan%*%as.matrix(t(GSE91061_model_data[,names(GSE91061_lan)])))

GSE91061_model_data$RS <- risk.GSE91061
GSE91061.data.point <- surv_cutpoint(GSE91061_model_data, time = "OS.time", event = "OS",
                                     variables = 'RS')
GSE91061.cutoff <- as.numeric(summary(GSE91061.data.point)[1])
GSE91061.cutoff



GSE91061.roc <- ggplotTimeROC(GSE91061_model_data$OS.time / 365,
                              GSE91061_model_data$OS,
                              risk.GSE91061,
                              mks = c(1,2,2.5))
GSE91061.km <- ggplotKMCox(data.frame(GSE91061_model_data$OS.time / 365,
                                      GSE91061_model_data$OS,
                                      ifelse(risk.GSE91061>=GSE91061.cutoff,'High','Low')),
                           title = 'GSE91061',add_text = '',
                           labs = c('High','Low'))


GSE91061.tide.point <- surv_cutpoint(GSE91061_tide_dat, time = "OS.time", event = "OS",
                                     variables = 'TIDE')
GSE91061.tide.cutoff <- as.numeric(summary(GSE91061.tide.point)[1])
GSE91061.tide.cutoff
GSE91061_tide_dat$TIDE.group <- ifelse(GSE91061_tide_dat$TIDE >=median(GSE91061_tide_dat$TIDE), 'High', 'Low')


GSE91061.tide.roc <- ggplotTimeROC(GSE91061_tide_dat$OS.time / 365,
                                   GSE91061_tide_dat$OS,
                                   GSE91061_tide_dat$TIDE,
                                   mks = c(1,2,2.5))
GSE91061.tide.km <- ggplotKMCox(data.frame(GSE91061_tide_dat$OS.time / 365,
                                           GSE91061_tide_dat$OS,
                                           GSE91061_tide_dat$TIDE.group),
                                title = 'TIDE', 
                                add_text = '',
                                labs = c('High','Low'))
GSE91061.tide.km

head(GSE91061_tide_dat[, c("Response", "TIDE", "IRS")])
table(GSE91061_tide_dat$Response)
GSE91061_tide_dat$Response[GSE91061_tide_dat$Response == 'PD'] <- 0
GSE91061_tide_dat$Response[GSE91061_tide_dat$Response == 'SD'] <- 0
GSE91061_tide_dat$Response[GSE91061_tide_dat$Response == 'CR'] <- 1
GSE91061_tide_dat$Response[GSE91061_tide_dat$Response == 'PR'] <- 1

library(pROC)
library(dplyr)
GSE91061.tide_roc_dat <- data.frame()
for (ty in c("TIDE", "IRS")) {
  print(ty)
  ty_dat <- GSE91061_tide_dat[, c(ty, "Response")]
  # ty_dat <- ty_dat[ty_dat[, ty] == '0' | ty_dat[, ty] == '1', ]
  ty_dat <- crbind2DataFrame(ty_dat)
  
  ty_roc <-  roc(as.formula(paste0('Response ~ ', ty)), data = ty_dat, smooth = T)
  tmp_roc_dat <- data.frame(TP = ty_roc$sensitivities,
                            FP = 1 - ty_roc$specificities,
                            Methods =paste0(ty, ' AUC = ', round(ty_roc$auc, 2)))
  GSE91061.tide_roc_dat <- rbind(GSE91061.tide_roc_dat, tmp_roc_dat)
}
head(GSE91061.tide_roc_dat)
GSE91061.tide_roc_dat$Methods=gsub('IRS','HYRS',GSE91061.tide_roc_dat$Methods)

GSE91061.Response_roc <- ggplot(GSE91061.tide_roc_dat, aes(x=FP,y=TP, fill=Methods))+
  geom_line(aes(colour=Methods),lwd=0.75)+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
  xlab('False positive rate')+ylab('True positive rate')+ 
  scale_colour_manual(values = mycolor[2:3]) +
  labs(title = 'Response') +
  theme_bw() +
  theme(legend.position=c(1,0),legend.justification=c(1,0),
        legend.background = element_rect(fill = NA, colour = NA))
GSE91061.Response_roc

#08.RiskScore临床特征####
tcga.risktype.cli=data.frame(tcga.subtype.cli,
                             Risktype=ifelse(tcga.module.risk$result$riskscorez>=0,'High','Low'),
                             Riskscore=tcga.module.risk$result$riskscore)
table(tcga.risktype.cli$Risktype)
rownames(tcga.risktype.cli)=tcga.risktype.cli$Samples

tcga_cox_datas <- tcga.risktype.cli
tcga_cox_datas <- crbind2DataFrame(tcga_cox_datas)
tcga_cox_datas$Risktype <- ifelse(tcga_cox_datas$Risktype == 'High', '1High', '0Low')
table(tcga_cox_datas$Stage)
tcga_cox_datas$Stage[tcga_cox_datas$Stage == 'I' | tcga_cox_datas$Stage == 'II'] <- 'I+II'
tcga_cox_datas$Stage[tcga_cox_datas$Stage == 'III' | tcga_cox_datas$Stage == 'IV'] <- 'III+IV'
table(tcga_cox_datas$Grade)
tcga_cox_datas$Grade[tcga_cox_datas$Grade == 'G1' | tcga_cox_datas$Grade == 'G2'] <- 'G1+G2'
tcga_cox_datas$Grade[tcga_cox_datas$Grade == 'G3' | tcga_cox_datas$Grade == 'G4'] <- 'G3+G4'
fig9a=list()
fig9a[[1]]=sig_boxplot(dat = tcga_cox_datas[,c("Stage","Riskscore")],
                       leg = 'Stage',ylab = 'Riskscore',
                       palette = ggsci::pal_jco()(10))
fig9a[[2]]=sig_boxplot(dat = tcga_cox_datas[,c("Grade","Riskscore")],
                       leg = 'Grade',ylab = 'Riskscore',
                       palette = ggsci::pal_jco()(10))
fig9a[[3]]=sig_boxplot(dat = tcga_cox_datas[,c("Age1","Riskscore")],
                       leg = 'Age',ylab = 'Riskscore',
                       palette = ggsci::pal_jco()(10))
fig9a[[4]]=sig_boxplot(dat = tcga_cox_datas[,c("Status","Riskscore")],
                       leg = 'Status',ylab = 'Riskscore',
                       palette = ggsci::pal_jco()(10))
fig9a[[5]]=sig_boxplot(dat = tcga_cox_datas[,c("Cluster","Riskscore")],
                       leg = 'Cluster',ylab = 'Riskscore',
                       palette = ggsci::pal_jco()(10))


fig9b=list()
table(tcga_cox_datas$Stage)
fig9b[[1]]=ggplotKMCox(data.frame(time = tcga_cox_datas[which(tcga_cox_datas$Stage=='III+IV'),]$OS.time/365
                                  , event = tcga_cox_datas[which(tcga_cox_datas$Stage=='III+IV'),]$OS
                                  , tcga_cox_datas[which(tcga_cox_datas$Stage=='III+IV'),]$Risktype)
                       ,add_text = '',title = 'Stage III+IV',show_confint=F
                       #,ggsci::pal_nejm()(10)[3:4]
                       ,labs = c('Low','High')
)

table(tcga_cox_datas$Age1)
fig9b[[2]]=ggplotKMCox(data.frame(time = tcga_cox_datas[which(tcga_cox_datas$Age1=='<60'),]$OS.time/365
                                  , event = tcga_cox_datas[which(tcga_cox_datas$Age1=='<60'),]$OS
                                  , tcga_cox_datas[which(tcga_cox_datas$Age1=='<60'),]$Risktype)
                       ,add_text = '',title = 'Age<60',show_confint=F
                       #,ggsci::pal_nejm()(10)[3:4]
                       ,labs = c('Low','High')
)
fig9b[[3]]=ggplotKMCox(data.frame(time = tcga_cox_datas[which(tcga_cox_datas$Age1=='>=60'),]$OS.time/365
                                  , event = tcga_cox_datas[which(tcga_cox_datas$Age1=='>=60'),]$OS
                                  , tcga_cox_datas[which(tcga_cox_datas$Age1=='>=60'),]$Risktype)
                       ,add_text = '',title = 'Age>=60',show_confint=F
                       #,ggsci::pal_nejm()(10)[3:4]
                       ,labs = c('Low','High')
)

table(tcga_cox_datas$Grade)
# fig9b[[4]]=ggplotKMCox(data.frame(time = tcga_cox_datas[which(tcga_cox_datas$Grade=='G1+G2'),]$OS.time/365
#                                   , event = tcga_cox_datas[which(tcga_cox_datas$Grade=='G1+G2'),]$OS
#                                   , tcga_cox_datas[which(tcga_cox_datas$Grade=='G1+G2'),]$Risktype)
#                        ,add_text = '',title = 'G1+G2',show_confint=F
#                        #,ggsci::pal_nejm()(10)[3:4]
#                        ,labs = c('Low','High')
# )
fig9b[[4]]=ggplotKMCox(data.frame(time = tcga_cox_datas[which(tcga_cox_datas$Grade=='G3+G4'),]$OS.time/365
                                  , event = tcga_cox_datas[which(tcga_cox_datas$Grade=='G3+G4'),]$OS
                                  , tcga_cox_datas[which(tcga_cox_datas$Grade=='G3+G4'),]$Risktype)
                       ,add_text = '',title = 'G3+G4',show_confint=F
                       #,ggsci::pal_nejm()(10)[3:4]
                       ,labs = c('Low','High')
)

#########HYRS 分组间突变特征差异#######
load('script/mut.RData')

col.selected=c('Aneuploidy Score','Homologous Recombination Defects','Fraction Altered',
               'Number of Segments','Nonsilent Mutation Rate')
fig10a=list()
for (i in 1:5){
  fig10a[[i]]=sig_boxplot(dat = ov.ri.tcga[,c("Risktype",col.selected[i])],
                          leg = 'Risktype',ylab = col.selected[i],
                          palette=ggsci::pal_npg("nrc")(10)[4:5])
}
fig10a=mg_merge_plot(fig10a,ncol = 5,nrow=1)
fig10a


pdf('results/pdf/Fig10b.pdf',he=7,wi=9)
oncoplot(maf=ov.maf2,clinicalFeatures = 'Risktype',
         top = 10,sortByAnnotation = T,
         annotationColor = list(Risktype=c(High='#3C5488FF',Low='#F39B7FFF')))
dev.off()



#######11.RiskScore分组之间的通路特征##########
getGeneFC=function(gene.exp,group,ulab=ulab,dlab=dlab){
  degs_C1_C3=limma_DEG(gene.exp, 
                          group,
                          ulab=ulab,
                          dlab=dlab)
  
  ## 未过滤任何基因
  degs_C1_C3_sig<-degs_C1_C3$DEG[which(degs_C1_C3$DEG$adj.P.Val <= 1),]
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(stringr)
  degs_C1_C3_sig_gene<-rownames(degs_C1_C3_sig)
  
  degs_gene_entz=bitr(degs_C1_C3_sig_gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
  degs_gene_entz <- dplyr::distinct(degs_gene_entz,SYMBOL,.keep_all=TRUE)
  
  gene_df <- data.frame(logFC=degs_C1_C3_sig$logFC,
                        SYMBOL = rownames(degs_C1_C3_sig))
  gene_df <- merge(gene_df,degs_gene_entz,by="SYMBOL")
  head(gene_df)
  
  geneList<-gene_df$logFC
  names(geneList)=gene_df$ENTREZID 
  head(geneList)
  
  geneList=sort(geneList,decreasing = T)
  head(geneList) ## 降序排列
  return(geneList)
}
#TCGA VS GSE51088
tcga.ri.geneList=getGeneFC(gene.exp=tcga_tpm_log_T[,tcga.risktype.cli$Samples],
                           group=tcga.risktype.cli$Risktype,
                           ulab='High',
                           dlab='Low')
GSE26712.ri.geneList=getGeneFC(gene.exp=GSE26712_exp[,GSE26712.risktype.cli$Samples],
                               group=GSE26712.risktype.cli$Risktype,
                               ulab='High',
                               dlab='Low')

h.all.gmt<-read.gmt("h.all.v7.5.1.entrez.gmt")
tcga.ri.hallmark.gsea<-GSEA(tcga.ri.geneList,TERM2GENE = h.all.gmt,seed=T)
GSE26712.ri.hallmark.gsea<-GSEA(GSE26712.ri.geneList,TERM2GENE = h.all.gmt,seed=T)

fig11a=clusterProfiler::dotplot(tcga.ri.hallmark.gsea,split=".sign",showCategory=nrow(tcga.ri.hallmark.gsea@result),
                                title='TCGA High vs Low',font.size=10)+facet_grid(~.sign)+
  scale_y_discrete(labels=function(x) str_remove(x,"HALLMARK_"))

fig11b=clusterProfiler::dotplot(GSE26712.ri.hallmark.gsea,split=".sign",showCategory=nrow(GSE26712.ri.hallmark.gsea@result),
                                title='GSE26712 High vs Low',font.size=10)+facet_grid(~.sign)+
  scale_y_discrete(labels=function(x) str_remove(x,"HALLMARK_"))
fig11ab=mg_merge_plot(fig11a,fig11b)
fig11ab

tcga.ri.hallmark.gsea.res=tcga.ri.hallmark.gsea@result
GSE26712.ri.hallmark.gsea.res=GSE26712.ri.hallmark.gsea@result
dim(tcga.ri.hallmark.gsea.res)
dim(GSE26712.ri.hallmark.gsea.res)

table(tcga.ri.hallmark.gsea.res$NES>0)
#  TRUE 
#    20
table(GSE26712.ri.hallmark.gsea.res$NES>0)
# FALSE  TRUE 
#  13    17 
rownames(tcga.ri.hallmark.gsea.res)=gsub("HALLMARK_","",rownames(tcga.ri.hallmark.gsea.res))
rownames(GSE26712.ri.hallmark.gsea.res)=gsub("HALLMARK_","",rownames(GSE26712.ri.hallmark.gsea.res))

risk.hallmark.union=Reduce(union,list(rownames(tcga.ri.hallmark.gsea.res),
                                      rownames(GSE26712.ri.hallmark.gsea.res)))
length(risk.hallmark.union)#36

risk.hallmark.heatmap.dat=matrix(0,nrow = 2,ncol = length(risk.hallmark.union))
rownames(risk.hallmark.heatmap.dat)=c('TCGA', 'GSE26712')
colnames(risk.hallmark.heatmap.dat)=risk.hallmark.union

risk.hallmark.heatmap.dat[1,match(rownames(tcga.ri.hallmark.gsea.res),colnames(risk.hallmark.heatmap.dat))]=tcga.ri.hallmark.gsea.res$NES
risk.hallmark.heatmap.dat[2,match(rownames(GSE26712.ri.hallmark.gsea.res),colnames(risk.hallmark.heatmap.dat))]=GSE26712.ri.hallmark.gsea.res$NES
range(risk.hallmark.heatmap.dat)
pdf('results/pdf/Fig11c.pdf',height =6,width =12)
Heatmap(as.matrix(risk.hallmark.heatmap.dat)
        , name = "NES"
        , col = circlize::colorRamp2(c(-3, 0, 3), c('slateblue1', 'white', 'orange'))
        , border = TRUE
        , show_column_names = T
        , show_column_dend = F
        , show_row_dend = F
        , cluster_columns=T
        , cluster_rows=F
        , rect_gp = gpar(col = "white", lwd = 1)
        , row_names_gp = gpar(fontsize = 10)
)
dev.off()


#####列线图#########
library(forestplot)
library(survcomp)
#单因素####
Age_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Age,
                             data=tcga_cox_datas))
Age_sig_cox_dat <- data.frame(Names=rownames(Age_sig_cox[[8]]),
                              HR = round(Age_sig_cox[[7]][,2],3),
                              lower.95 = round(Age_sig_cox[[8]][,3],3),
                              upper.95 = round(Age_sig_cox[[8]][,4],3),
                              p.value=round(Age_sig_cox[[7]][,5],3))
Age_sig_cox_dat

Stage_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Stage,
                               data=tcga_cox_datas))

Stage_sig_cox_dat <- data.frame(Names=rownames(Stage_sig_cox[[8]]),
                                HR = round(Stage_sig_cox[[7]][,2],3),
                                lower.95 = round(Stage_sig_cox[[8]][,3],3),
                                upper.95 = round(Stage_sig_cox[[8]][,4],3),
                                p.value=round(Stage_sig_cox[[7]][,5],3))
Stage_sig_cox_dat

Grade_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Grade,
                               data=tcga_cox_datas))
Grade_sig_cox_dat <- data.frame(Names=rownames(Grade_sig_cox[[8]]),
                                HR = round(Grade_sig_cox[[7]][,2],3),
                                lower.95 = round(Grade_sig_cox[[8]][,3],3),
                                upper.95 = round(Grade_sig_cox[[8]][,4],3),
                                p.value=round(Grade_sig_cox[[7]][,5],3))
Grade_sig_cox_dat

riskscore_sig_cox<-summary(coxph(formula=Surv(OS.time, OS)~Riskscore,
                                 data=tcga_cox_datas))
riskscore_sig_cox <- data.frame(Names=rownames(riskscore_sig_cox[[8]]),
                                HR = round(riskscore_sig_cox[[7]][,2],3),
                                lower.95 = round(riskscore_sig_cox[[8]][,3],3),
                                upper.95 = round(riskscore_sig_cox[[8]][,4],3),
                                p.value=round(riskscore_sig_cox[[7]][,5],3))
riskscore_sig_cox

sig_cox_dat <- rbind(Age_sig_cox_dat,
                     Stage_sig_cox_dat,
                     Grade_sig_cox_dat,
                     riskscore_sig_cox)
data.sig <- data.frame(Names=sig_cox_dat$Names,
                       p.value=sig_cox_dat$p.value,
                       sig_cox_dat$HR,
                       sig_cox_dat$lower.95,
                       sig_cox_dat$upper.95)
data.sig <- crbind2DataFrame(data.sig)
rownames(data.sig) <- c("Age",
                        "Stage",
                        "Grade",
                        "RiskScore"
)
data.sig$Names <- rownames(data.sig)
data.sig$p.value=ifelse(data.sig$p.value==0,'<0.001',data.sig$p.value)
data.sig



#多因素#####
muti_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Age+Stage+Grade+Riskscore, 
                              data=tcga_cox_datas))
muti_cox_dat <- data.frame(Names=rownames(muti_sig_cox[[8]]),
                           HR = round(muti_sig_cox[[7]][,2],3),
                           lower.95 = round(muti_sig_cox[[8]][,3],3),
                           upper.95 = round(muti_sig_cox[[8]][,4],3),
                           p.value=round(muti_sig_cox[[7]][,5],3))
data.muti <- data.frame(Names=muti_cox_dat$Names,
                        p.value=muti_cox_dat$p.value,
                        muti_cox_dat$HR,
                        muti_cox_dat$lower.95,
                        muti_cox_dat$upper.95)
data.muti <- crbind2DataFrame(data.muti)
rownames(data.muti) <- c("Age",
                         "Stage",
                         "Grade",
                         "RiskScore")
data.muti$Names <- rownames(data.muti)
data.muti$p.value=ifelse(data.muti$p.value==0,'<0.001',data.muti$p.value)
data.muti


##############列线图
nom.plot=mg_nomogram(data.frame(RiskScore=tcga_cox_datas$Riskscore,
                                Age=tcga_cox_datas$Age),
                     os = tcga_cox_datas$OS.time,
                     status = tcga_cox_datas$OS,
                     mks = c(1,3,5))


#AUC
library("timeROC")
tcga_cox_auc=tcga_cox_datas
tcga_cox_auc$Stage=as.numeric(as.factor(tcga_cox_auc$Stage))
tcga_cox_auc$Riskscore=as.numeric(tcga_cox_auc$Riskscore)
tcga_cox_auc$Age=as.numeric(as.factor(tcga_cox_auc$Age))
ROC.DSST.Age=timeROC(T=tcga_cox_auc$OS.time/365,
                     delta=tcga_cox_auc$OS,
                     marker=tcga_cox_auc$Age,
                     cause=1,weighting="marginal",
                     times=c(1,2,3,4,5),
                     iid=TRUE)
ROC.DSST.Age$AUC
ROC.DSST.Stage=timeROC(T=tcga_cox_auc$OS.time/365,
                       delta=tcga_cox_auc$OS,
                       marker=tcga_cox_auc$Stage,
                       cause=1,weighting="marginal",
                       times=c(1,2,3,4,5),
                       iid=TRUE)
ROC.DSST.Stage$AUC
ROC.DSST.Risk=timeROC(T=tcga_cox_auc$OS.time/365,
                      delta=tcga_cox_auc$OS,
                      marker=tcga_cox_auc$Riskscore,
                      cause=1,weighting="marginal",
                      times=c(1,2,3,4,5),
                      iid=TRUE)
ROC.DSST.Risk$AUC
ROC.DSST.Nomo<-timeROC(T=tcga_cox_auc$OS.time/365,
                       delta=tcga_cox_auc$OS,
                       marker=tcga_cox_auc$Riskscore,
                       other_markers=as.matrix(tcga_cox_auc[,c("Age")]),
                       cause=1,
                       weighting="cox",
                       times=c(1,2,3,4,5),
                       iid=F)
ROC.DSST.Nomo$AUC
#scales::show_col(mg_colors[c(1,10:12,4,5,7:9)])
pdf('results/pdf/Fig12_AUC.pdf',height = 6,width = 6)
#par(mai=c(0.5,0.5,0.5,2))
plotAUCcurve(ROC.DSST.Nomo,conf.int=F,col=mg_colors[1])
plotAUCcurve(ROC.DSST.Stage,conf.int=F,col=mg_colors[4],add=TRUE)
plotAUCcurve(ROC.DSST.Age,conf.int=F,col=mg_colors[2],add=TRUE)
plotAUCcurve(ROC.DSST.Risk,conf.int=F,col=mg_colors[3],add=TRUE)
legend(1,xjust=0,yjust=2,c("Nomogram"
                           ,'Stage'
                           ,"RiskScore",'Age')
       ,col=mg_colors[c(1,4,3,2)],lty=1,lwd=2,xpd = TRUE)

dev.off()