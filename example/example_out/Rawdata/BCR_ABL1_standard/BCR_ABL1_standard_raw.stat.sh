set -e 
echo Raw_statistic start: `date` 
cd /data01/Work/yuhuan/Software/TargetFusion/example/example_out/Rawdata/BCR_ABL1_standard

NanoStat -t 6 \
    --fastq BCR_ABL1_standard.merge.fastq.gz \
    -o /data01/Work/yuhuan/Software/TargetFusion/example/example_out/Rawdata/BCR_ABL1_standard \
    -n BCR_ABL1_standard.rawfq.nanostat.xls
echo Raw_statistic fq end: `date` 
