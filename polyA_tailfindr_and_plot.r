#devtools::install_url('https://cran.r-project.org/src/contrib/Archive/rbokeh/rbokeh_0.5.1.tar.gz', type = "source", dependencies = TRUE)
devtools::install_github("adnaniazi/tailfindr")
#remotes::install_github('adnaniazi/tailfindr', build = TRUE, build_opts = c("--no-resave-data", "--no-manual"), force = TRUE)
library(tailfindr)
library(ggplot2)
library(PupillometryR)
df_0h <- find_tails(fast5_dir = '~/tb2/CPSF/cpsf30_0h_mRNA_basecalled/gzip_fast5',
                    save_dir = '~/tb2/CPSF/cpsf30_0h_mRNA_basecalled/polyA',
                    csv_filename = 'rna_tails_cpsf30_0h.csv',
                    basecall_group = 'Basecall_1D_000',
                    num_cores = 30,
                    save_plots = T,
                    plotting_library = "rbokeh"
)
df_48h <- find_tails(fast5_dir = '~/tb2/CPSF/cpsf30_48h_mRNA_basecalled/gzip_fast5',
                 save_dir = '~/tb2/CPSF/cpsf30_48h_mRNA_basecalled/polyA',
                 csv_filename = 'rna_tails_cpsf30_48h.csv',
                 basecall_group = 'Basecall_1D_000',
                 num_cores = 30,
                 save_plots = T,
                 plot_debug_traces = TRUE,
                 plotting_library = "rbokeh"
                  )


df_annotate_0h= annotate_tails(
  sam_file = "~/tb2/CPSF/cpsf30_0h_mRNA_basecalled/polyA/cpsf30_0h_mRNA.sam",
  tails_csv_file = "~/tb2/CPSF/cpsf30_0h_mRNA_basecalled/polyA/rna_tails_cpsf30_0h.csv",
  output_file = "~/tb2/CPSF/cpsf30_0h_mRNA_basecalled/polyA/cpsf30_0h_annot_tails.csv"
)

df_annotate_48h= annotate_tails(
  sam_file = "d:/CPSF/cpsf30_48h_mRNA_basecalled/polyA/cpsf30_48h_mRNA.sam",
  tails_csv_file = "d:CPSF/cpsf30_48h_mRNA_basecalled/polyA/rna_tails_cpsf30_48h.csv",
  output_file = "d:/CPSF/cpsf30_48h_mRNA_basecalled/polyA/cpsf30_48h_annot_tails.csv"
)
df_annotate_0h=read.csv("D:\\CPSF\\cpsf30_0h_mRNA_basecalled\\polyA/cpsf30_0h_annot_tails.csv",header = T)
df_annotate_48h=read.csv("D:\CPSF\cpsf30_48h_mRNA_basecalled\polyA/")
df_annotate_0h$time="0h"
df_annotate_48h$time="48h"
df_annotate= rbind(df_annotate_0h,df_annotate_48h)
head(df_annotate)
df_annotate2=df_annotate[which(df_annotate$tail_length < 300),]
df_annotate2=df_annotate2[which(df_annotate2$tail_length >1),]

ggplot(df_annotate2,aes(x=tail_length,fill=time))+geom_density(alpha=0.5,bw=10)
df_beta=df_annotate2[which(df_annotate2$transcript_id %in% c("Tb927.1.2370:mRNA","Tb927.1.2330:mRNA","Tb927.1.2350:mRNA","Tb927.1.2390:mRNA")),]
df_alpha=df_annotate2[which(df_annotate2$transcript_id %in% c("Tb927.1.2340:mRNA","Tb927.1.2360:mRNA","Tb927.1.2380:mRNA","Tb927.1.2400:mRNA")),]
ggplot(df_beta,aes(x=tail_length,fill=time))+geom_density(alpha=0.5,bw=10)
ggplot(df_alpha,aes(x=tail_length,fill=time))+geom_density(alpha=0.5,bw=10)

df_trrm2_0h_tail=read.csv("~/sdb/trypanosoma/trrm2_polya_0h/annot_tails.csv",header=T)
df_trrm2_48h_tail=read.csv("~/sdb/trypanosoma/trrm2_polya_48h/annot_tails.csv",header=T)
df_trrm2_0h_tail$time="trrm2_0h"
df_trrm2_48h_tail$time="trrm2_48h"
df_annot3=rbind(df_annotate2,df_trrm2_0h_tail,df_trrm2_48h_tail)
df_annot3=df_annot3[which(df_annot3$tail_length < 300),]
df_annot3=df_annot3[which(df_annot3$tail_length >1),]
ggplot(df_annot3,aes(x=tail_length,fill=time))+geom_density(alpha=0.5,bw=10)

ggplot(df_annotate2,aes(x=time,y=tail_length))+
  geom_flat_violin(aes(fill=time),position=position_nudge(x=.25),color="black")+
  geom_jitter(aes(color=time),width = 0.15,size=0.1)+
  geom_boxplot(notch = TRUE, width=.1,position = position_nudge(x=0.25),fill="white",size=0.5)+
  coord_flip()+
  theme_classic()+
  theme(axis.text.y = element_blank())+
  ylab("3'UTR length (nt)") +
  xlab("RNAi")
