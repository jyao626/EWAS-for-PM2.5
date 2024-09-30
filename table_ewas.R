#table1
setwd("Z:/EWAS/Data")
getwd()
library(haven)

beta_NoImpute <- readRDS("./CARDIA_Y15_Combined_preprocessed_NoRCP_KNNimpute.rds")
beta_NoImpute2 <- readRDS("./CARDIA_Y20_Combined_preprocessed_NoRCP_KNNimpute.rds")






linktable <- read.csv("CARDIA_Sample_LinkTable_2181.csv")
linktable_y15=linktable[linktable$Year=='y15',]
linktable_y15$ID=substr(linktable_y15$CARDIA_unique,1,12)
linktable_y15$ID=as.numeric(linktable_y15$ID)


linktable_y20=linktable[linktable$Year=='y20',]
linktable_y20$ID=substr(linktable_y20$CARDIA_unique,1,12)
linktable_y20$ID=as.numeric(linktable_y20$ID)

#
link_match=linktable_y15[!is.na(match(linktable_y15$BARCODE,colnames(beta_NoImpute))),]

link_match20=linktable_y20[!is.na(match(linktable_y20$BARCODE,colnames(beta_NoImpute2))),]

#
covarfull=read_sas('Z:/AHA/PM/Old/Data/CAC/pm_cvd_cardia_20220207.sas7bdat')
cov15=covarfull[covarfull$visit==15,]
cov20=covarfull[covarfull$visit==20,]
#
bias=read.table('CARDIA_Y15_Methy_1089_Info_OnlyCellProp_PCs.txt')
bias1=bias[-1,]
colnames(bias1)=c('ID','QC','CD8T_Y15','CD4T_Y15','NK_Y15','BCELL_Y15','MONO_Y15','GRAN_Y15','PC1_Y15','PC2_Y15','PC3_Y15','PC4_Y15','PC5_Y15','PC6_Y15','PC7_Y15','PC8_Y15')
bias2=bias1[,-2]
bias2$ID=as.numeric(bias2$ID)
bias2$CD8T_Y15=as.numeric(bias2$CD8T_Y15)
bias2$CD4T_Y15=as.numeric(bias2$CD4T_Y15)
bias2$NK_Y15=as.numeric(bias2$NK_Y15)
bias2$BCELL_Y15=as.numeric(bias2$BCELL_Y15)
bias2$MONO_Y15=as.numeric(bias2$MONO_Y15)
bias2$GRAN_Y15=as.numeric(bias2$GRAN_Y15)
bias2$PC1_Y15=as.numeric(bias2$PC1_Y15)
bias2$PC2_Y15=as.numeric(bias2$PC2_Y15)
bias2$PC3_Y15=as.numeric(bias2$PC3_Y15)
bias2$PC4_Y15=as.numeric(bias2$PC4_Y15)
bias2$PC5_Y15=as.numeric(bias2$PC5_Y15)
bias2$PC6_Y15=as.numeric(bias2$PC6_Y15)
bias2$PC7_Y15=as.numeric(bias2$PC7_Y15)
bias2$PC8_Y15=as.numeric(bias2$PC8_Y15)


#
cov15$degree[cov15$degree==8|cov15$degree==9]=NA

summary(link_match2$smoke)
linky15_cardia1=merge(link_match,cov15,by='ID',all.x=T)
linky15_cardia2=merge(linky15_cardia1,bias2,by='ID',all.x=T)
  

linky15_cardia3=merge(link_match20,cov20,by='ID',all.x=T)
linky15_cardia4=merge(linky15_cardia3,bias5,by='ID',all.x=T)

linky15_cardia5=merge(link_match20,cov15,by='ID',all.x=T)
linky15_cardia6=merge(linky15_cardia5,bias2,by='ID',all.x=T)

linky15_cardia2=merge(linky15_cardia2,ndvi15,by.x='ID',by.y='id',all.x=T)
linky15_cardia2=merge(linky15_cardia2,park,by='ID',all.x=T)

identical(rownames(beta_NoImpute),linky15_cardia2$BARCODE)

linky15_cardia2=linky15_cardia2[match(colnames(beta_NoImpute),linky15_cardia2$BARCODE),]

linky15_cardia4=linky15_cardia4[match(colnames(beta_NoImpute2),linky15_cardia4$BARCODE),]

linky15_cardia6=linky15_cardia6[match(colnames(beta_NoImpute2),linky15_cardia6$BARCODE),]
library(boot)
identical(colnames(beta_NoImpute),linky15_cardia2$BARCODE)
identical(colnames(beta_NoImpute2),linky15_cardia4$BARCODE)
identical(colnames(beta_NoImpute2),linky15_cardia6$BARCODE)


xx1=linky15_cardia2[,c('pm25_20220130',"A01SEX","A01AGE2","bmi",'degree','CENTER','smoke','race_y15','CD8T_Y15','CD4T_Y15','NK_Y15','BCELL_Y15','MONO_Y15','GRAN_Y15','PC1_Y15','PC2_Y15','PC3_Y15','PC4_Y15','PC5_Y15','PC6_Y15','PC7_Y15','PC8_Y15')] 


#table1 replace with other variable
mean(linky15_cardia2$bmi,na.rm=T)
sd(linky15_cardia2$bmi,na.rm=T)

#table2
xx=linky15_cardia3[,c('pm25_20220130','BC_20220130','NH4_20220130','NIT_20220130','OM_20220130','SO4_20220130','SOIL_20220130','SS_20220130')]

round(cor(xx,use="complete.obs"), 2)

#table3
results1 <- cpg.assoc(beta_NoImpute2,linky15_cardia2$NH4_20220130, covariates =xx1)
results4 <- cpg.assoc(beta_NoImpute2,linky15_cardia2$NIT_20220130, covariates =xx1)
#table 4

CpG.name = "cg25491402"
linky15_cardia6$CpG.level <- beta_NoImpute2[CpG.name,]
linky15_cardia6.b=linky15_cardia6[linky15_cardia6$A01SEX==1,]
m1=lm(CpG.level ~ NIT_20220130+pm25_20220130+A01SEX+degree+A01AGE2+bmi+CENTER+smoke+act+race_y15+CD8T_Y15+CD4T_Y15+NK_Y15+BCELL_Y15+MONO_Y15+GRAN_Y15+PC1_Y15+PC2_Y15+PC3_Y15+PC4_Y15+PC5_Y15+PC6_Y15+PC7_Y15+PC8_Y15,data=linky15_cardia6.b)
summary(m1)
m2=lm(CpG.level ~ NIT_20220130+pm25_20220130+A01SEX+degree+A01AGE2+bmi+CENTER+smoke+as.factor(race_y15)+NIT_20220130*as.factor(race_y15)+CD8T_Y15+CD4T_Y15+NK_Y15+BCELL_Y15+MONO_Y15+GRAN_Y15+PC1_Y15+PC2_Y15+PC3_Y15+PC4_Y15+PC5_Y15+PC6_Y15+PC7_Y15+PC8_Y15,data=linky15_cardia6.b)
summary(m2)

