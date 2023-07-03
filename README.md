# TargetFusion
Introducing a Pipeline for Long-Read Fusion Detection

We have developed a comprehensive pipeline that enables efficient long-read fusion detection. This pipeline incorporates various functionalities such as data debarcoding, quality control, mapping, statistical analysis, fusion detection (using LongGF), and false positive identification and filtering. By utilizing a simple configuration file called "sample.list," users can easily run multiple samples simultaneously and obtain accurate results with marked false positive fusions and statistical information.

For users who prefer not to use the entire pipeline, there is an alternative option. They can utilize the[**src/reform_LongGF_Result.py**](https://github.com/HuanYuu/TargetFusion/blob/main/src/reform_LongGF_Result.py) script, which allows for the reforming and filtering of LongGF log results. This additional filtering step has been shown to reduce more than 90% of false positive results in our testing.

This pipeline offers a streamlined and efficient approach to long-read fusion detection, ensuring reliable and accurate results. Whether users choose to use the entire pipeline or the individual reforming and filtering script, they can expect improved outcomes and reduced false positive rates.
  
# Prerequisite
Run TargetFusion.py from raw data qc to fusion detection:
1. python3
2. Porechop (debarcode and cut adapter)
3. NanoStat (raw and clean data statistic)
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
One libraryID can maps to several sample id (TargetFusion can debarcode).  
One sampleID can maps to several data path(fq_data_path), (TargetFusion will merge all these data to analysis).  
  
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
Run TargetFusion. See **A real example** for an example, see **TargetFusion usage** for detail.  
  
# A real example
See [example/run_rna_fusion.sh](https://github.com/HuanYuu/TargetFusion/blob/main/example/run_rna_fusion.sh) for a real example. We also provide a demo data **BCR_ABL1.test.fq.gz** and main results files in [example](https://github.com/HuanYuu/TargetFusion/blob/main/example).  

# TargetFusion usage
## Simple usage:  
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
-t: thread for run 
-b: target region bed file 
-o: output path 
-rna: is RNA sequence data 
-gene: target gene list 

## TargetFusion result
### directory and shell
When you run command above, an output example:  
* **Rawdata**  (debarcode and get raw fq data)  
* **QC**  (raw data QC, get clean data)  
* **Mapping**  (mapping, get mapped bam and target statistic results)  
* **Fusion**  (fusion detection)  
* **Report**  (save final statistic and fusion results)  
* **Shell**  (save every sample's run shell)  
* **sample.list.reconfig**  (reconfiged sample.list)  

In the **"Shell"** directory, you will find all integrated analysis shells named **"run_sampleID_analysis.sh"** for each sample. If your data consists of pooled samples with different barcodes that require debarcoding, you will also find a script named **"step0_debarcode.sh"** in the "Shell" directory. In this case, it is necessary to run the "step0_debarcode.sh" script before proceeding with the fusion analysis.

Once the "step0_debarcode.sh" script has finished executing (or if debarcoding is not required), you can proceed to run all the fusion analysis shells together.
  
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
The user will receive statistical information and fusion results in the **Report** directory. For more details, please refer to the [example/example_out/Report](https://github.com/HuanYuu/TargetFusion/blob/main/example/example_out/Report).
1. sampleID.mosdepth_reform.xls
2. sampleID.LongGF.reform.xls
