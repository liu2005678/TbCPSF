setwd("/media/logen/sdb/trypanosoma/dicistronic/bedops/")
fq=dir("/media/logen/tb2/CPSF/long_reads/fq")
fq=fq[grep("_chop.fq",fq)]
base=sub("_chop.fq","",fq)
fq
base
chop_cmd=paste("porechop -i ", fq, " -o ", base, "_chop.fq")
cmd=paste0(
   ##minimap2
  "minimap2 -ax map-ont -t30 --secondary=no ../../ref_genome/tb927genome_56.fasta /media/logen/tb2/CPSF/long_reads/fq/", fq, " > ", base, ".sam \n",
   #samtools turn the file format to bam, and sort according to chromosome location (default)
   "samtools view -u -@ 30 ", base, ".sam | samtools sort -@ 30 -o ", base, ".bam - \n",
  # change the file from bam to bed format, which is acceptable for intersect evaluation by bedops tools
  "bamToBed -i ", base, ".bam > ", base, ".bed \n", 
  #mapping the intersect of individual reads to annotated CDS in tritrypdb
   "bedmap --echo-ref-row-id --echo-map-id --echo-map-size  --echo-overlap-size   ", base, ".bed  cds.bed > ", base, ".out \n"
)

### optional, write the commands into file
write.table(cmd,file="cmd",quote=F,col.names = F,row.names = F)
# excute the commands, if not applicable, excute in the terminal with the previous file with the following lines
# chmod 755 cmd
# ./cmd
for(i in 1:length(cmd)){
  system(cmd[i])
}
