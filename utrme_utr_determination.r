files=read.csv("/media/logen/sdc/trypanosoma/CPSF/samples.csv",header = T)
files=files[19:24,]
files$repeats=rep("repeat1",times=12)
files$repeats[9]="repeat2"
files$time=paste0(files$time,"h")
files

files$clean_reads1=paste0("/media/logen/sdc/trypanosoma/CPSF/short_reads/bmk/rep4/clean/",files$sequence_provider,"_",files$gene,"_",files$time,"_",files$repeats,".clean1.fq.gz")
files$clean_reads2=paste0("/media/logen/sdc/trypanosoma/CPSF/short_reads/bmk/rep4/clean/",files$sequence_provider,"_",files$gene,"_",files$time,"_",files$repeats,".clean2.fq.gz")
files$orphan1=paste0("/media/logen/sdc/trypanosoma/CPSF/short_reads/bmk/rep4/clean/",files$sequence_provider,"_",files$gene,"_",files$time,"_",files$repeats,".orphan1.fq.gz")
files$orphan2=paste0("/media/logen/sdc/trypanosoma/CPSF/short_reads/bmk/rep4/clean/",files$sequence_provider,"_",files$gene,"_",files$time,"_",files$repeats,".orphan2.fq.gz")
trim_cmd=paste("java -jar /home/logen/programs/trim/trimmomatic-0.39.jar PE -threads 3",  files$raw_reads_1,files$raw_reads_2,files$clean_reads1,
               files$orphan1,files$clean_reads2,files$orphan2,
                "LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 &")
trim_cmd
#write.csv(files,file="~/tb2/CPSF/short_reads/samples.csv",row.names = F,quote=F)
write.table(trim_cmd,file="/media/logen/sdc/trypanosoma/CPSF/short_reads/bmk/rep4/trim_cmd.txt",quote=F,row.names = F,col.names = F)

mkdir1=sub(".fq.gz","", files$clean_reads1)
mkdir1
mkdir2=sub(".fq.gz","", files$clean_reads2)
mkdir2
mkdira=paste0("mkdir ",c(mkdir1,mkdir2))
mkdira
dirs=sub("mkdir ","",mkdira)
filename=sub("mkdir /media/logen/sdc/trypanosoma/CPSF/short_reads/bmk/rep4/clean/","",mkdira)
cleana=c(files$clean_reads1,files$clean_reads2)
for(i in 1:length(mkdira)){ system(mkdira[i])}
ln_cmd=paste0("ln -s ",cleana," ",dirs,"/",filename,".fq.gz")
for(i in 1:length(ln_cmd)){system(ln_cmd[i])}
base=sub("/media/logen/sdc/trypanosoma/CPSF/short_reads/bmk/rep4/clean/","",mkdir1)
base=sub(".clean1","",base)
base
config_files=paste0("/media/logen/sdc/trypanosoma/CPSF/utrme/UTRme/Configuration_Files/",      base, "_configuration.txt")
mk_config=paste0("cp /media/logen/sdc/trypanosoma/CPSF/utrme/UTRme/Configuration_Files/cpsf30apt48h_configuration.txt /media/logen/sdc/trypanosoma/CPSF/utrme/UTRme/Configuration_Files/",
                base, "_configuration.txt")
for(i in 1:length(mk_config)){system(mk_config[i])}
attach1=paste0("echo \"first_value:",mkdir1,"/\" >>",config_files)
attach2=paste0("echo \"second_value:",mkdir2,"/\" >> ", config_files)
attach3=paste0("echo \"basename_value:",base," \">>",config_files)
attach=c(attach2,attach1,attach3)
for(i in 1:length(attach)){system(attach[i])}
write.table(system_cmd,file="/media/logen/sdc/trypanosoma/CPSF/utrme/system.cmd",quote=F,row.names = F,col.names = F)
python_utrme=paste0("python utrme.py ",config_files)
setwd("/media/logen/sdc/trypanosoma/CPSF/utrme/UTRme")
for(i in 1:length(python_utrme)){system(python_utrme[i])}
write.table(python_utrme,file="utrme_cmd",quote = F,row.names = F,col.names = F)


wdr=files[13:18,]
wdr
config_files=paste0("/media/logen/sdc/trypanosoma/CPSF/utrme/UTRme/Configuration_Files/", wdr$samplename, "_configuration.txt")
mk_config=paste0("cp /media/logen/sdc/trypanosoma/CPSF/utrme/UTRme/Configuration_Files/cpsf30apt48h_configuration.txt /media/logen/sdc/trypanosoma/CPSF/utrme/UTRme/Configuration_Files/",
                 wdr$samplename, "_configuration.txt")
for(i in 1:length(mk_config)){system(mk_config[i])}
mkdir1=paste0("/media/logen/sdc/trypanosoma/CPSF/short_reads/bmk/BMK_DATA_20220802103146_1/Data/",wdr$samplename,"/fastq1")
mkdir2=paste0("/media/logen/sdc/trypanosoma/CPSF/short_reads/bmk/BMK_DATA_20220802103146_1/Data/",wdr$samplename,"/fastq2")
attach1=paste0("echo \"first_value:",mkdir1,"/\" >>",config_files)
attach2=paste0("echo \"second_value:",mkdir2,"/\" >> ", config_files)
attach3=paste0("echo \"basename_value:",wdr$samplename," \">>",config_files)
attach=c(attach2,attach1,attach3)
for(i in 1:length(attach)){system(attach[i])}
config_mmap=sub("configuration.txt","configuration_mmap.txt",config_files)
python_utrme=paste0("python utrme.py ",config_mmap)
setwd("/media/logen/sdc/trypanosoma/CPSF/utrme/UTRme")
for(i in 1:length(python_utrme)){system(python_utrme[i])}
write.table(python_utrme,file="utrme_cmd",quote = F,row.names = F,col.names = F)
