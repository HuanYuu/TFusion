# TargetFusion
A pipeline for (target) long-read fusion detection bases on LongGF. We integrate data debarcode, QC, mapping, statistics, LongGF fusion detection AND LongGF-result mark/filter together. Users can simply config a sample.list to run different samples at the same time and finally get false positive marked fusion results.  
For non-target use, users can only use the src/reform_LongGF_Result.py to mark/fiter LongGF log result to futher filter the false positive results.  
