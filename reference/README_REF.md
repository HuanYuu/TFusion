# reference and gene_strand file
User can download reference(the pseudoautosomal region (PAR) on the Y is annotated) fasta and gtf, and generate faidx file using command below:  
``` download and index
wget http://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gzip -dc Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa
wget http://ftp.ensembl.org/pub/current_gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz
```
  
Then, user can generate gene_strand file using command below:  
``` generate gene_strand file
gzip -dc Homo_sapiens.GRCh38.109.gtf.gz
cat Homo_sapiens.GRCh38.109.gtf |  awk -F"\t" '{if($1!~/^#/) print}' | awk -F"\t" '{if($3=="gene")print}' |  awk -F"\t" '{split($9,arr,"; gene_name \"");split(arr[2],garr,"\";");print(garr[1]"\t"$1":"$4"-"$5"\t"$7)}' | awk -F"\t" '{if($1!="")print}' | cut -f 1,3 | sort -u > Homo_sapiens.GRCh38.109.gene.strand.txt
```  

There is an example file "Homo_sapiens.GRCh38.109.gene.strand.txt" in the "reference" directory.  
