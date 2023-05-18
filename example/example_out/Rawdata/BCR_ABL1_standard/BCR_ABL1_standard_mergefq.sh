set -e 
echo Raw_merge_fastq start: `date` 
cd /data01/Work/yuhuan/Software/TargetFusion/example/example_out/Rawdata/BCR_ABL1_standard
cp /data01/Work/yuhuan/Software/TargetFusion/example/BCR_ABL1.test.fq.gz BCR_ABL1_standard.merge.fastq.gz
md5sum BCR_ABL1_standard.merge.fastq.gz > BCR_ABL1_standard.merge.fastq.gz.md5.txt 
echo Raw_merge_fastq end: `date` 
