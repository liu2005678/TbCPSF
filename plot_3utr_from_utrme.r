library(dplyr)
library(ggplot2)
library(reshape2)
library(PupillometryR)


setwd("~/tb2/CPSF/utrme/UTRme")
files=dir("./",recursive = T)
full3=files[grep("3-full-report.tab",files)]

### not to consider multi-mapping reads first (tubulin excluded)
full3=full3[1:12]
full3
group=sub(".*/Reports/","",full3)
group=sub("-3-full-report.tab","",group)
group
full3tab=data.frame()
for(i in 1:length(full3)){
  tab=read.delim(full3[i],header = T, comment.char = "*")
  tab$group=group[i]
  full3tab=rbind(full3tab,tab)
}
head(full3tab)
dffun=function(best3tab,length,occurrences){
  df_len=data.frame(length=rep(0,times=sum(best3tab$occurrences)),group=rep("0",times=sum(best3tab$occurrences)))
  flag=1
  for(i in 1:nrow(best3tab)){
   end=flag+best3tab$occurrences[i] - 1
   df_len$length[flag:end]=rep(best3tab$length[i], times=best3tab$occurrences[i])
   df_len$group[flag:end]=rep(best3tab$group[i], times=best3tab$occurrences[i])
   flag=end+1
  }
  return(df_len)
}
df_len=dffun(full3tab,length,occurrence)
x=df_len

#draw pictures
setwd("~/tb2/CPSF/utrme/rplots")
draw_pic=function(filename,datf,length,group){
  filename1=paste0(filename, "_density_plot.pdf")
  pdf(filename1)
  myplot = ggplot(datf,aes(length,fill=group))+
    geom_density(alpha=0.5)+
    xlab("3'UTR length (nt)") 
  print(myplot)
  dev.off()
  filename2=paste0(filename, "_raincloud_plot.pdf")
  pdf(filename2)
  myplot = ggplot(datf,aes(x=group,y=length))+
    geom_flat_violin(aes(fill=group),position=position_nudge(x=.25),color="black")+
    geom_jitter(aes(color=group),width = 0.15,size=0.1)+
    #geom_boxplot(notch = TRUE, width=.1,position = position_nudge(x=0.25),fill="white",size=0.5)+
    coord_flip()+
    theme_classic()+
   theme(axis.text.y = element_blank())+
    ylab("3'UTR length (nt)") 
  print(myplot)
  dev.off()
}

df_fip1=df_len[which(df_len$group %in% c("apt_FIP1_0h_repeat1","nuohe_FIP1_48h_repeat1","nuohe_FIP1_72h_repeat1")),]
df_fip2=df_len[which(df_len$group %in% c("apt_FIP1_0h_repeat2","apt_FIP1_48h_repeat1"   ,  "apt_FIP1_72h_repeat1")),]
df_cpsf1=df_len[which(df_len$group %in% c("nuohe_CPSF30_0h_repeat1" , "nuohe_CPSF30_48h_repeat1", "nuohe_CPSF30_72h_repeat1")),]
df_cpsf2=df_len[which(df_len$group %in% c("apt_CPSF30_0h_repeat1",  "apt_CPSF30_48h_repeat1" ,  "apt_CPSF30_72h_repeat1" )),]
df_fip0h=df_len[which(df_len$group %in% c("apt_FIP1_0h_repeat2","apt_FIP1_0h_repeat1")),]
df_cpsf0h=df_len[which(df_len$group %in% c("nuohe_CPSF30_0h_repeat1", "apt_CPSF30_0h_repeat1")),]  
df_0hall=rbind(df_fip0h,df_cpsf0h)

draw_pic("fip1_repeat1",df_fip1,length = length,group = group)
draw_pic("fip1_repeat2",df_fip2,length = length,group = group)
draw_pic("cpsf30_repeat1",df_cpsf1,length = length,group = group)
draw_pic("cpsf30_repeat2",df_cpsf2,length = length,group = group)
draw_pic("cpsf_0h_compare",df_cpsf0h,length = length,group = group)
draw_pic("fip_0h_compare",df_fip0h,length = length,group = group)
draw_pic("all_0h_compare",df_0hall,length = length,group = group)





### Considering multi-mapping reads (tubulin included)
setwd("~/tb2/CPSF/utrme/UTRme")
files=dir("./",recursive = T)
full3=files[grep("3-full-report.tab",files)]
full3=full3[13:24]
full3
full3tab=data.frame()
for(i in 1:length(full3)){
  tab=read.delim(full3[i],header = T, comment.char = "*")
  tab$group=group[i]
  full3tab=rbind(full3tab,tab)
}
head(full3tab)
df_len=dffun(full3tab,length,occurrence)
head(df_len)
y=df_len
#draw pictures
setwd("~/tb2/CPSF/utrme/rplots")
df_fip1=df_len[which(df_len$group %in% c("apt_FIP1_0h_repeat1","nuohe_FIP1_48h_repeat1","nuohe_FIP1_72h_repeat1")),]
df_fip2=df_len[which(df_len$group %in% c("apt_FIP1_0h_repeat2","apt_FIP1_48h_repeat1"   ,  "apt_FIP1_72h_repeat1")),]
df_cpsf1=df_len[which(df_len$group %in% c("nuohe_CPSF30_0h_repeat1" , "nuohe_CPSF30_48h_repeat1", "nuohe_CPSF30_72h_repeat1")),]
df_cpsf2=df_len[which(df_len$group %in% c("apt_CPSF30_0h_repeat1",  "apt_CPSF30_48h_repeat1" ,  "apt_CPSF30_72h_repeat1" )),]
df_fip0h=df_len[which(df_len$group %in% c("apt_FIP1_0h_repeat2","apt_FIP1_0h_repeat1")),]
df_cpsf0h=df_len[which(df_len$group %in% c("nuohe_CPSF30_0h_repeat1", "apt_CPSF30_0h_repeat1")),]  
df_0hall=rbind(df_fip0h,df_cpsf0h)

draw_pic("fip1_repeat1",df_fip1,length = length,group = group)
draw_pic("fip1_repeat2",df_fip2,length = length,group = group)
draw_pic("cpsf30_repeat1",df_cpsf1,length = length,group = group)
draw_pic("cpsf30_repeat2",df_cpsf2,length = length,group = group)
draw_pic("cpsf_0h_compare",df_cpsf0h,length = length,group = group)
draw_pic("fip_0h_compare",df_fip0h,length = length,group = group)
draw_pic("all_0h_compare",df_0hall,length = length,group = group)


###draw tubulin
head(full3tab)
tab_tuba=full3tab[which(full3tab$gene %in% c("Tb927.1.2340","Tb927.1.2360","Tb927.1.2380","Tb927.1.2400")),]
tab_tubb=full3tab[which(full3tab$gene %in% c("Tb927.1.2350","Tb927.1.2370","Tb927.1.2390","Tb927.1.2330")),]
tab_h2a=full3tab[which(full3tab$gene %in% c("Tb927.7.2820","Tb927.7.2830","Tb927.7.2840","Tb927.7.2850","Tb927.7.2860","Tb927.7.2870",
             "Tb927.7.2880","Tb927.7.2890","Tb927.7.2900","Tb927.7.2910","Tb927.7.2920","Tb927.7.2930","Tb927.7.2940")),]

draw_pic2=function(df_len,length,occurrences){
#df_len=dffun(df,length,occurrences)
df_fip1=df_len[which(df_len$group %in% c("apt_FIP1_0h_repeat1","nuohe_FIP1_48h_repeat1","nuohe_FIP1_72h_repeat1")),]
df_fip2=df_len[which(df_len$group %in% c("apt_FIP1_0h_repeat2","apt_FIP1_48h_repeat1"   ,  "apt_FIP1_72h_repeat1")),]
df_cpsf1=df_len[which(df_len$group %in% c("nuohe_CPSF30_0h_repeat1" , "nuohe_CPSF30_48h_repeat1", "nuohe_CPSF30_72h_repeat1")),]
df_cpsf2=df_len[which(df_len$group %in% c("apt_CPSF30_0h_repeat1",  "apt_CPSF30_48h_repeat1" ,  "apt_CPSF30_72h_repeat1" )),]
df_fip0h=df_len[which(df_len$group %in% c("apt_FIP1_0h_repeat2","apt_FIP1_0h_repeat1")),]
df_cpsf0h=df_len[which(df_len$group %in% c("nuohe_CPSF30_0h_repeat1", "apt_CPSF30_0h_repeat1")),]  
df_0hall=rbind(df_fip0h,df_cpsf0h)
draw_pic("fip1_repeat1",df_fip1,length = length,group = group)
draw_pic("fip1_repeat2",df_fip2,length = length,group = group)
draw_pic("cpsf30_repeat1",df_cpsf1,length = length,group = group)
draw_pic("cpsf30_repeat2",df_cpsf2,length = length,group = group)
draw_pic("cpsf_0h_compare",df_cpsf0h,length = length,group = group)
draw_pic("fip_0h_compare",df_fip0h,length = length,group = group)
draw_pic("all_0h_compare",df_0hall,length = length,group = group)
}



df_tubb=dffun(tab_tubb,length,occurrences)
head(tab_tubb)
head(df_tubb)
df_tubb$all="all"
df_tubb=df_tubb[which(df_tubb$length > 350), ]
df_tubb=df_tubb[which(df_tubb$length < 390), ]
ggplot(df_tubb,aes(x=all,y=length))+
  geom_flat_violin(aes(fill=all),position=position_nudge(x=.25),color="black")+
  geom_jitter(width = 0.1,size=0.1)+
  coord_flip()+
  theme_classic()+
  theme(axis.text.y = element_blank())+
  ylab("3'UTR length (nt)") 


head(df_len)
head(tab_h2a)
df_h2a=dffun(tab_h2a,length,occurrences)
head(df_h2a)
setwd("~/tb2/CPSF/utrme/h2aplots")
df_h2a$all="all"
df_h2a=df_h2a[which(df_h2a$length > 180), ]
df_h2a=df_h2a[which(df_h2a$length < 350), ]
ggplot(df_h2a,aes(x=all,y=length))+
  geom_flat_violin(aes(fill=all),position=position_nudge(x=.25),color="black")+
  geom_jitter(width = 0.1,size=0.1)+
  coord_flip()+
  theme_classic()+
  theme(axis.text.y = element_blank())+
  ylab("3'UTR length (nt)") 

head(tab_tuba)
dftuba=dffun(tab_tuba)
setwd("~/tb2/CPSF/utrme/")
dftuba$all="all"
dftuba=dftuba[which(dftuba$length <150), ]
ggplot(dftuba,aes(x=all,y=length))+
  geom_flat_violin(aes(fill=all),position=position_nudge(x=.25),color="black")+
  geom_jitter(width = 0.1,size=0.1)+
  coord_flip()+
  theme_classic()+
  theme(axis.text.y = element_blank())+
  ylab("3'UTR length (nt)") 


setwd("~/tb2/CPSF/utrme/alpha_tub_plots")
draw_pic2(dftuba,occurrences = occurrences ,length = length)
setwd("~/tb2/CPSF/utrme/beta_tub_plots")
draw_pic2(df_tubb,occurrences = occurrences ,length = length)
setwd("~/tb2/CPSF/utrme/h2a_plots")
draw_pic2(df_h2a,occurrences = occurrences ,length = length)

