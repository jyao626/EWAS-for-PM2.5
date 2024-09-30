#figure1 and 2
library(CMplot)

datamanhat2=readRDS(file="figure1.Rda")


SNPs <- datamanhat2[datamanhat2$pvalue < (0.05 / (841639*7)), 1]
genes <- atamanhat2[atamanhat2$pvalue < (0.05 / (841639*7)), 10]


CMplot(datamanhat2, plot.type="m",LOG10=TRUE, threshold=c(0.05/(841639*7)),highlight.text=genes,
       threshold.lty=c(1), threshold.lwd=c(1), threshold.col=c("black"), amplify=TRUE,highlight=SNPs,
       chr.den.col=NULL, cex = 0.4, signal.cex=c(0.4),signal.pch=c(19),highlight.cex=c(0.6),
       file="jpg",memo="",dpi=600,file.output=TRUE,verbose=TRUE,width=14,height=6,chr.labels.angle=45)


datamanhat2=readRDS(file="figure2.Rda")


SNPs <- datamanhat2[datamanhat2$pvalue < (0.05 / (841639*7)), 1]
genes <- atamanhat2[atamanhat2$pvalue < (0.05 / (841639*7)), 10]


CMplot(datamanhat2, plot.type="m",LOG10=TRUE, threshold=c(0.05/(841639*7)),highlight.text=genes,
       threshold.lty=c(1), threshold.lwd=c(1), threshold.col=c("black"), amplify=TRUE,highlight=SNPs,
       chr.den.col=NULL, cex = 0.4, signal.cex=c(0.4),signal.pch=c(19),highlight.cex=c(0.6),
       file="jpg",memo="",dpi=600,file.output=TRUE,verbose=TRUE,width=14,height=6,chr.labels.angle=45)

#figure3 and 4
library(GenomicRanges)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(clusterProfiler)
library(missMethyl)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(DOSE)
library(ReactomePA)
.libPaths()
annot <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

cpg_list=c("cg20813589", "cg03953749", "cg14180696", "cg25119406", "cg25491402", "cg11919271", "cg06391325",
           "cg19930065", "cg06635351", "cg04202267", "cg14431020", "cg25990696", "cg15083386",
           "cg25090121", "cg22358291", "cg21953769", "cg18817487")
x <- enrichPathway(gene=EntrezID$sig.eg,               
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 1,
                   qvalueCutoff  = 0.2,
                   readable      = TRUE)

edox <- setReadable(x, 'org.Hs.eg.db', 'ENTREZID')


p3 <- cnetplot(edox,showCategory = 10, circular = TRUE, colorEdge = T) 
x=p3+ scale_color_manual("Node",values=c('red', 'blue'),labels = c("Gene", "Pathway"))

edo <- enrichDGN(EntrezID$sig.eg,
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 1,
                 qvalueCutoff  = 0.2,
                 readable      = TRUE)

edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
p3 <- cnetplot(edox,showCategory = 10, circular = TRUE, colorEdge = T) 

x=p3+ scale_color_manual("Node",values=c('red', 'blue'),labels = c("Gene", "Pathway"))
#box plot
Race=c('Black','Black','Black','Black','Black','Black','Black','Black','Black','Black','Black','Black','Black','Black','Black','Black','Black','Black','Black',
       'White','White','White','White','White','White','White','White','White','White','White','White','White','White','White','White','White','White','White')

Coefficients=c(0.028,0.031,0.027,0.009,0.014,0.014,0.018,0.018,0.025,0.023,0.02,0.015,0.027,0.013,0.012,0.023,0.023,0.015,0.013,
               0.013,0.02,0.017,0.007,0.008,0.008,0.015,0.009,0.018,0.015,0.009,0.007,0.021,0.005,0.009,0.011,0.012,0.008,0.006)

new=data.frame(Coefficients,Race)

pdf('output.pdf')
x=boxplot(Coefficients~Race,
          data=new,
          main="Boxplot",
          xlab="Race",
          ylab="Absolute values of coefficients",
          col="orange",
          border="brown"
)

