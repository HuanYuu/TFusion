set -e
echo Statistic_result_reform start: `date` 

cd /data01/Work/yuhuan/Software/TargetFusion/example/example_out/Mapping/BCR_ABL1_standard
python /data01/Work/yuhuan/Software/TargetFusion/src/mosdepth_bedstat_reform.py \
    -d /data01/Work/yuhuan/Software/TargetFusion/example/example_out/Mapping/BCR_ABL1_standard \
    -b /data01/Work/yuhuan/Software/TargetFusion/src/example_target_region.bed \
    -rawstat /data01/Work/yuhuan/Software/TargetFusion/example/example_out/Rawdata/BCR_ABL1_standard/BCR_ABL1_standard.rawfq.nanostat.xls \
    -qcstat /data01/Work/yuhuan/Software/TargetFusion/example/example_out/QC/BCR_ABL1_standard/BCR_ABL1_standard.qc.nanostat.xls \
    -bam BCR_ABL1_standard.sorted.bam \
    -p BCR_ABL1_standard

cd /data01/Work/yuhuan/Software/TargetFusion/example/example_out/Report
ln -sf ../Mapping/BCR_ABL1_standard/BCR_ABL1_standard.mosdepth_reform.xls .
echo Statistic_result_reform end: `date` 
