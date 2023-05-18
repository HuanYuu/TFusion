set -e 
echo LongGF fusion start: `date` 

cd /data01/Work/yuhuan/Software/TargetFusion/example/example_out/Mapping/BCR_ABL1_standard
samtools sort -n -O BAM BCR_ABL1_standard.sorted.bam -o BCR_ABL1_standard.readname.sorted.bam

cd /data01/Work/yuhuan/Software/TargetFusion/example/example_out/Fusion/BCR_ABL1_standard
LongGF \
    /data01/Work/yuhuan/Software/TargetFusion/example/example_out/Mapping/BCR_ABL1_standard/BCR_ABL1_standard.readname.sorted.bam \
    /data01/Work/yuhuan/Software/TargetFusion/reference/Homo_sapiens.GRCh38.109.gtf \
    90 20 90 0 0 5 > /data01/Work/yuhuan/Software/TargetFusion/example/example_out/Fusion/BCR_ABL1_standard/BCR_ABL1_standard.LongGF.log

python /data01/Work/yuhuan/Software/TargetFusion/src/reform_LongGF_Result.py \
 -i /data01/Work/yuhuan/Software/TargetFusion/example/example_out/Fusion/BCR_ABL1_standard/BCR_ABL1_standard.LongGF.log \
 -bam /data01/Work/yuhuan/Software/TargetFusion/example/example_out/Mapping/BCR_ABL1_standard/BCR_ABL1_standard.sorted.bam \
 -gene /data01/Work/yuhuan/Software/TargetFusion/src/example_fusion.gene.list.txt \
 -sf /data01/Work/yuhuan/Software/TargetFusion/reference/Homo_sapiens.GRCh38.109.gene.strand.txt \
 -o /data01/Work/yuhuan/Software/TargetFusion/example/example_out/Fusion/BCR_ABL1_standard/BCR_ABL1_standard.LongGF.reform.xls


cd /data01/Work/yuhuan/Software/TargetFusion/example/example_out/Report
ln -sf ../Fusion/BCR_ABL1_standard/BCR_ABL1_standard.LongGF.reform.xls .

rm -f /data01/Work/yuhuan/Software/TargetFusion/example/example_out/Mapping/BCR_ABL1_standard/BCR_ABL1_standard.readname.sorted.bam

echo LongGF fusion end: `date` 
