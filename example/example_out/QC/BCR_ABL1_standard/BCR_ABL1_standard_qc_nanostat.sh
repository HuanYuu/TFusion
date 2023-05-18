echo QC_statistic start: `date` 
cd /data01/Work/yuhuan/Software/TargetFusion/example/example_out/QC/BCR_ABL1_standard
NanoStat -t 6 \
    --fastq BCR_ABL1_standard.filter.fastq.gz \
    -o /data01/Work/yuhuan/Software/TargetFusion/example/example_out/QC/BCR_ABL1_standard \
    -n BCR_ABL1_standard.qc.nanostat.xls
echo QC_statistic end: `date` 
