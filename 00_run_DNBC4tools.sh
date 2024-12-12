conda activate DNBC4tools_v207
module load languages/R/4.0.4


refpath=/WUH721818ALE6L4/database/references/DNBC4tools/mm10
refgtf=/WUH721818ALE6L4/database/references/mm10/mm10.gtf
datapath=/WUH721818ALE6L4/datasets/private/TBK1

mkdir /WUH721818ALE6L4/projects/TBK1/00_DNBC4tools
cd /WUH721818ALE6L4/projects/TBK1/00_DNBC4tools

process="data,count"

for sample in vsv51_2 vsv51_1 GSK8612_2 GSK8612_1 Control_2 Control_1 Comb_2 Comb_1;
do
  cfq1=${datapath}/${sample}/${sample}_cDNA_R1.fastq.gz
  cfq2=${datapath}/${sample}/${sample}_cDNA_R2.fastq.gz
  ofq1=${datapath}/${sample}/${sample}_oligo_R1.fastq.gz
  ofq2=${datapath}/${sample}/${sample}_oligo_R2.fastq.gz
  mkdir -p $sample
  DNBC4tools run --name $sample --thread 15 \
  --process ${process} \
  --genomeDir ${refpath} --gtf ${refgtf} \
  --cDNAfastq1 ${cfq1} --cDNAfastq2 ${cfq2} \
  --oligofastq1 ${ofq1} --oligofastq2 ${ofq2}
done
