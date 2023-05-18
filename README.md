# TargetFusion
A pipeline for (target) long-read fusion detection. 
   
We integrated data debarcode, QC, mapping, statistics, fusion detection(based on LongGF) and false positive mark/filter together. Users can simply config a sample.list to run different samples at the same time and finally get false positive fusion marked results and statistic results.  
  
For non-pipeline using, users can only utilize the [**src/reform_LongGF_Result.py**](https://github.com/HuanYuu/TargetFusion/blob/main/src/reform_LongGF_Result.py) to reform and mark/filter LongGF log result. This futher filter can reduce more than 90% fasle positvie results in our test.  
  
# Prerequisite
To run TargetFusion.py from raw data to fusion results:
1. python3
2. Porechop (debarcode and cut adapter)
3. NanoStat (raw data statistic)
4. NanoFilt (short and low quality reads filter)
5. minimap2 (mapping)
6. samtools (bam transform and samtools depth statistic)
7. sambamba (sort)
8. mosdepth (bed region coverage statistic)
9. LongGF (fusion detection)
  
Only run src/reform_LongGF_Result.py to reform and mark/fiter LongGF log result:
1. samtools
  
# How to use
**step1**:  
Install PREREQUISITE software and download TargetFusion:  
```Download TargetFusion:
git clone https://github.com/HuanYuu/TargetFusion.git
```
  
**step2**:  
Download reference files and construct gene_stand file follow [reference/README_REF.md](https://github.com/HuanYuu/TargetFusion/blob/main/reference/README_REF.md)  
  
**step3**:  
Move or link 'Porechop/\*' to 'TargetFusion/Porechop' (OR modify porechop path in TargetFusion.py), see [Porechop/README_porechop.md](https://github.com/HuanYuu/TargetFusion/blob/main/Porechop/README_porechop.md)  

**step4**:  
Prepare sample.list, target_region.bed, target.fusion.gene.list.  
  
1)*sample list* (tab or space separate): include "sampleID", "libraryID", "barcodeID", "fq_data_path".  
One libraryID can map to several sample id (TargetFusion can debarcode).  
One sampleID can have several fq_data_path, one fq_data_path per line (TargetFusion will merge all these data to analysis).  
  
sample.list example1:
|#sampleID  |libraryID  |barcodeID  |fq_data_path  |
|----       |----       |----       |----          |
|sample1    |LB01       |-          |path/data1.fq |
|sample2    |LB02       |-          |path/data2.fq |

  sample.list example2:  
|#sampleID  |libraryID  |barcodeID  |fq_data_path  |  
|----       |----       |----       |----          |
|sample1    |LB01       |-          |path/data1.fq |
|sample2    |LB02       |-          |path/data2.fq |
|sample2    |LB03       |-          |path/data3.fq |
|sample3    |LB04       |BC01       |path/data4.fq |
|sample4    |LB04       |BC02       |path/data4.fq |
  
2)*target_region.bed*: include 4 column, chr, start, end, symbol. Symbol can be target gene name. Software will statistic every symbol's coverage.
  
3)*target.fusion.gene.list*: include target gene names. One gene per one line.  
  
**step5**:  
Run TargetFusion. See **TargetFusion usage** for detail.  
  
# TargetFusion usage
An example. See [example/run_rna_fusion.sh](https://github.com/HuanYuu/TargetFusion/blob/main/example/run_rna_fusion.sh) for comprehensive file and results:  
```example:
python TargetFusion.py \
    -s sample.list \
    -t 6 \
    -b example_target_region.bed \
    -o example_out \
    -rna \
    -gene example_fusion.gene.list.txt
```
  
-s: a sample list  
-t: thread  
-b: target region bed file  
-o: output path  
-rna: is RNA sequence data  
-gene: target gene list  

# TargetFusion result
## directory and shell
When you run command above, an output example:  
* **Rawdata**  (debarcode and get raw fq data)  
* **QC**  (raw data QC, get clean data)  
* **Mapping**  (mapping, get mapped bam and target statistic results)  
* **Fusion**  (fusion detection)  
* **Report**  (save final statistic and fusion results)  
* **Shell**  (save every sample's run shell)  
* **sample.list.reconfig**  (reconfiged sample.list)  

For every sample, integrated analysis shell is stored in "Shell" directory, name as "run_sampleID_analysis.sh".  
If you have pooling data with different barcodes which needs to do debarcode, "step0_debarcode.sh" will be in "Shell" directory, and you need to run step0 first.  
After "step0_debarcode.sh" finished (or you don't need to debarcode), you can run all fusion analysis shells together.  
  
A run_sampleID_analysis.sh example:  
```run_example_analysis.sh
set -e
echo run all analysis of sampleID start: `date`
## step1
sh /testpath/Rawdata/sampleID/sampleID_mergefq.sh
## step2
sh /testpath/Rawdata/sampleID/sampleID_raw.stat.sh &
## step3
sh /testpath/QC/sampleID/sampleID_qc.sh
## step4
sh /testpath/QC/sampleID/sampleID_qc_nanostat.sh &
## step5
sh /testpath/Mapping/sampleID/sampleID_mapping.sh
## step6
sh /testpath/Fusion/sampleID/sampleID_fusion_lgf.sh &
## step7
sh /testpath/Mapping/sampleID/sampleID_mapstat.sh &
## step8
wait
sh /testpath/Mapping/sampleID/sampleID_reformstat.sh
echo run all analysis of sampleID end: `date`
```
  
For a local running example:
```example
cd Shell
for i in `ls run*_analysis.sh`;do nohup sh $i > $i.log & done
```

## report
User will get statistic and fusion result in Report. See[example/example_out/Report](https://github.com/HuanYuu/TargetFusion/blob/main/example/example_out/Report)for details:  
1. sampleID.mosdepth_reform.xls
2. sampleID.LongGF.reform.xls
