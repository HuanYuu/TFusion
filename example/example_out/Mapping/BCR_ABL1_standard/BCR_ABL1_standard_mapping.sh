set -e 
echo Mapping start: `date` 

cd /data01/Work/yuhuan/Software/TargetFusion/example/example_out/Mapping/BCR_ABL1_standard

minimap2 -t 6 -ax splice \
    /data01/Work/yuhuan/Software/TargetFusion/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    /data01/Work/yuhuan/Software/TargetFusion/example/example_out/QC/BCR_ABL1_standard/BCR_ABL1_standard.filter.fastq.gz |
samtools view -@ 6 -hbS \
    -t /data01/Work/yuhuan/Software/TargetFusion/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai \
    -o BCR_ABL1_standard.raw.bam

sambamba sort \
    -t 6 \
    -m 6G \
    --tmpdir BCR_ABL1_standard_sambamba_sort.tmp \
    -o BCR_ABL1_standard.sorted.bam \
    BCR_ABL1_standard.raw.bam

sambamba index -t 6 BCR_ABL1_standard.sorted.bam

rm -f BCR_ABL1_standard.raw.bam

md5sum BCR_ABL1_standard.sorted.bam > BCR_ABL1_standard.sorted.bam.md5.txt

echo Mapping end: `date` 
