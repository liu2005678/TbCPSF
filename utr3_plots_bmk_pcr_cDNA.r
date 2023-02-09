library(dplyr)
library(ggplot2)
library(reshape2)
library(PupillometryR)

#draw pictures
setwd("/media/logen/sdc/trypanosoma/CPSF/apa/bmk")
cpsf30=read.delim("bmk_cpsf30_inter.bed",header=F,sep="\t")
cpsf30rnai=read.delim("bmk_cpsf30RNAi_inter.bed",header=F,sep="\t")
head(cpsf30)
for(i in 1:nrow(cpsf30)){
  if(cpsf30$V6[i] == "-") {cpsf30$u3length[i]=cpsf30$V13[i] - cpsf30$V2[i]}
  else if(cpsf30$V6[i] == "+") {cpsf30$u3length[i]=cpsf30$V3[i] - cpsf30$V14[i]}
}
for(i in 1:nrow(cpsf30rnai)){
  if(cpsf30rnai$V6[i] == "-") {cpsf30rnai$u3length[i]=cpsf30rnai$V13[i] - cpsf30rnai$V2[i]}
  else if(cpsf30rnai$V6[i] == "+") {cpsf30rnai$u3length[i]=cpsf30rnai$V3[i] - cpsf30rnai$V14[i]}
}
nrow(cpsf30)
nrow(cpsf30rnai)
cpsf30$group="tbCPSF30 RNAi (-)"
cpsf30rnai$group="tbCPSF30 RNAi (+)"
cpsf=rbind(cpsf30[,c("V10", "group","u3length")],cpsf30rnai[,c("V10", "group","u3length")])
colnames(cpsf)=c("gene","group","length")
head(cpsf)
cpsf=cpsf[which((cpsf$length < 3000) & (cpsf$length >0)),]
###draw tubulin
df_tuba=cpsf[which(cpsf$gene %in% c("Tb927.1.2340","Tb927.1.2360","Tb927.1.2380","Tb927.1.2400")),]
df_tubb=cpsf[which(cpsf$gene %in% c("Tb927.1.2350","Tb927.1.2370","Tb927.1.2390","Tb927.1.2330")),]
df_h2a=cpsf[which(cpsf$gene %in% c("Tb927.7.2820","Tb927.7.2830","Tb927.7.2840","Tb927.7.2850","Tb927.7.2860","Tb927.7.2870",
                                   "Tb927.7.2880","Tb927.7.2890","Tb927.7.2900","Tb927.7.2910","Tb927.7.2920","Tb927.7.2930","Tb927.7.2940")),]

df_tubb=df_tubb[which((df_tubb$length > 340) & (df_tubb$length < 390)), ]
df_tuba=df_tuba[which((df_tuba$length < 150) & (df_tuba$length > 70)), ]
df_h2a=df_h2a[which((df_h2a$length < 390) & (df_h2a$length > 150)), ]

draw_pic=function(filename,datf,length,group){
  filename1=paste0(filename, "_density_plot.pdf")
  pdf(filename1,paper = "a4")
  myplot = ggplot(datf,aes(length,fill=group))+
    geom_density(alpha=0.5,size=0.1)+
    theme_classic(base_size =  8)+
    xlab(paste0("3'UTR length (nt)"))
  print(myplot)
  dev.off()
  filename2=paste0(filename, "_raincloud_plot.pdf")
  pdf(filename2)
  myplot = ggplot(datf,aes(x=group,y=length))+
    geom_flat_violin(aes(fill=group),position=position_nudge(x=.25),color="black",alpha=0.5,size=0.1)+ #
    geom_jitter(aes(color=group),width = 0.15,size=0.05,alpha=0.5)+ #
    #geom_boxplot(notch = TRUE, width=.1,position = position_nudge(x=0.25),fill="white",size=0.5)+
    coord_flip()+
    theme_classic()+
    theme(axis.text.y = element_blank())+
    ylab(paste0(filename,"_3'UTR length (nt)")) 
  print(myplot)
  dev.off()
}

draw_pic("cpsf",cpsf,"length","group")
draw_pic("tuba",df_tuba,"length","group")
draw_pic("tubb",df_tubb,"length","group")
draw_pic("h2a",df_h2a,"length","group")

# Wilcoxon rank sum test with continuity correction
wilcox.test(df_h2a$length[which(df_h2a$group == "tbCPSF30 RNAi (-)")],
            df_h2a$length[which(df_h2a$group == "tbCPSF30 RNAi (+)")],alternative = "two.sided")
# p-value = 1.828e-10
wilcox.test(df_tuba$length[which(df_tuba$group == "tbCPSF30 RNAi (-)")],
            df_tuba$length[which(df_tuba$group == "tbCPSF30 RNAi (+)")],alternative = "two.sided")
# p-value = 2.759e-09
wilcox.test(df_tubb$length[which(df_tubb$group == "tbCPSF30 RNAi (-)")],
            df_tubb$length[which(df_tubb$group == "tbCPSF30 RNAi (+)")],alternative = "two.sided")
# p-value < 2.2e-16
