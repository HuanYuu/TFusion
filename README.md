# TargetFusion
A pipeline for (target) long-read fusion detection bases on LongGF. 
   
We integrate data debarcode, QC, mapping, statistics, LongGF fusion detection AND LongGF-result mark/filter together. Users can simply config a sample.list to run different samples at the same time and finally get false positive marked fusion results.  
  
For non-target use, users can only use the src/reform_LongGF_Result.py to mark/fiter LongGF log result to futher filter the false positive results.  
  
# Prerequisite
To run TargetFusion.py from raw data to fusion results:
1. python3
2. Porechop
3. NanoStat
4. NanoPlot
5. samtools
6. minimap2
7. LongGF

Only run src/reform_LongGF_Result.py to reform and mark/fiter LongGF false positive result:
1. samtools

# How to use
**step1**: Install PREREQUISITE software and download TargetFusion:  
```step1:
mkdir TargetFusion
cd TargetFusion
git clone https://github.com/HuanYuu/TargetFusion.git
```
  
**step2**: Download reference file and construct gene_stand file follow reference/README_REF.md  
  
**step3**: Move or link Porechop/\* to 'TargetFusion/Porechop' (OR modify porechop path in TargetFusion.py)  

**step4**: prepare a sample.list, target_region.bed, target.fusion.gene.list.  
  
A sample list(tab or space separate): include "sampleID", "libraryID", "barcodeid", "fq_data_path". 
One libraryID can share by several sample id (TargetFusion can debarcode).  
One sampleID can have several fq_data_path (TargetFusion will merge all these data to analysis).  
  
target_region.bed: include 4 column, chr/start/end/symbol.  
symbol can be target gene name. software will statistic every symbol's coverage.
  
target.fusion.gene.list: include gene names. One gene per one line.  
  
**step5**: Run TargetFusion. See TargetFusion usage for detail.  
  
# TargetFusion usage
An example:  
```example:
python TargetFusion.py \
    -s sample.list.pep_test \
    -t 6 \
    -b example_target_region.bed \
    -o Analysis_ensemble_peptest \
    -rna \
    -gene example_fusion.gene.list.txt
```
-s: sample.list  
-t: thread  
-b: target bed file  
-o: output path  
-rna: RNA sequence data  
-gene: target gene list  

When you run command above, an output example:  
*Rawdata
*QC
*Mapping
*Fusion
*Report
*Shell
*sample.list.reconfig

Then, all the shell will store in "Shell" directory.  
If you have pooling data with different barcode that need to do debarcode, "step0_debarcode.sh" will in "Shell" directory, and you need run step0 first.
After "step0_debarcode.sh" finished or you don't need debarcode, you can run all these fusion analysis shell together.  
For a local running example:
```example
cd Shell
for i in `ls run*_analysis.sh`;do nohup sh $i > $i.log & done
```
