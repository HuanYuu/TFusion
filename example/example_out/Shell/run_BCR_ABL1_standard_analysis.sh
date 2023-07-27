set -e
echo run all analysis of BCR_ABL1_standard start: `date` 

## step1
#sh /data01/Work/yuhuan/Software/TargetFusion/example/example_out/Rawdata/BCR_ABL1_standard/BCR_ABL1_standard_mergefq.sh 

## step2
#sh /data01/Work/yuhuan/Software/TargetFusion/example/example_out/Rawdata/BCR_ABL1_standard/BCR_ABL1_standard_raw.stat.sh & 

## step3
#sh /data01/Work/yuhuan/Software/TargetFusion/example/example_out/QC/BCR_ABL1_standard/BCR_ABL1_standard_qc.sh 

## step4
#sh /data01/Work/yuhuan/Software/TargetFusion/example/example_out/QC/BCR_ABL1_standard/BCR_ABL1_standard_qc_nanostat.sh & 

## step5
sh /data01/Work/yuhuan/Software/TargetFusion/example/example_out/Mapping/BCR_ABL1_standard/BCR_ABL1_standard_mapping.sh 

## step6
sh /data01/Work/yuhuan/Software/TargetFusion/example/example_out/Fusion/BCR_ABL1_standard/BCR_ABL1_standard_fusion_lgf.sh & 

## step7
sh /data01/Work/yuhuan/Software/TargetFusion/example/example_out/Mapping/BCR_ABL1_standard/BCR_ABL1_standard_mapstat.sh & 

## step8
wait
sh /data01/Work/yuhuan/Software/TargetFusion/example/example_out/Mapping/BCR_ABL1_standard/BCR_ABL1_standard_reformstat.sh 

echo run all analysis of BCR_ABL1_standard end: `date` 
