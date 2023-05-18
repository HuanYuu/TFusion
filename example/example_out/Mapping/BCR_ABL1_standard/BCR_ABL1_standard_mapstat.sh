set -e 
echo Mapping_statistic start: `date` 
cd /data01/Work/yuhuan/Software/TargetFusion/example/example_out/Mapping/BCR_ABL1_standard

mosdepth -t 6 -T 1,50,100,200,500 \
    -b /data01/Work/yuhuan/Software/TargetFusion/src/example_target_region.bed \
    BCR_ABL1_standard \
    BCR_ABL1_standard.sorted.bam

echo Mapping_statistic end: `date` 
