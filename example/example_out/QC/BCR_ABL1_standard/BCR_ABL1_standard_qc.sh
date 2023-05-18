set -e 
echo QC_remove_adapter start: `date` 
cd /data01/Work/yuhuan/Software/TargetFusion/example/example_out/QC/BCR_ABL1_standard

python /data01/Work/yuhuan/Software/TargetFusion/Porechop*/porechop-runner.py -t 6 \
    -i /data01/Work/yuhuan/Software/TargetFusion/example/example_out/Rawdata/BCR_ABL1_standard/BCR_ABL1_standard.merge.fastq.gz \
    -o BCR_ABL1_standard.merge.deadapter.fastq \
    --format fastq
echo QC_remove_adapter end: `date` 

echo QC_filter start: `date` 
NanoFilt -q 7 -l 100 \
    BCR_ABL1_standard.merge.deadapter.fastq | 
gzip > BCR_ABL1_standard.filter.fastq.gz 

md5sum BCR_ABL1_standard.filter.fastq.gz > BCR_ABL1_standard.filter.fastq.gz.md5.txt 

rm -f BCR_ABL1_standard.merge.deadapter.fastq 

echo QC_filter end: `date` 
