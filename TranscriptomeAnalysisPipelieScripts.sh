#Transcriptom analysis scripts

#Quality check with FastQC

fasqc -t 9 *.fq.gz]

#TRIMMING 

mkdir /home/ravwoy/WORK2/BoczkowskaM/DATA/trimmed_rav

ls -1 /home/ravwoy/WORK2/BoczkowskaM/DATA/raw/*_1.fq.gz | sed 's/.*raw\///' | sed 's/\.fq\.gz//' > raw_filesR1; ls -1 /home/ravwoy/WORK2/BoczkowskaM/DATA/raw/*_2.fq.gz | sed 's/.*raw\///' | sed 's/\.fq\.gz//' > raw_filesR2; ls -1 /home/ravwoy/WORK2/BoczkowskaM/DATA/raw/*_1.fq.gz | sed 's/.*raw\///' | sed 's/_.\.fq\.gz//' > 3rdcolumn; paste raw_filesR1 raw_filesR2 3rdcolumn > raw_files.pairs;

cat raw_files.pairs | head -n 4 | tail -n 1 | while read R1 R2 G; do /home/ravwoy/.conda/envs/RAVWOY1/bin/SeqPurge -in1 /home/ravwoy/WORK2/BoczkowskaM/DATA/raw/$R1.fq.gz -in2 /home/ravwoy/WORK2/BoczkowskaM/DATA/raw/$R2.fq.gz -out1 /home/ravwoy/WORK2/BoczkowskaM/DATA/trimmed_rav/$R1.seqpurge.fastq.gz -out2 /home/ravwoy/WORK2/BoczkowskaM/DATA/trimmed_rav/$R2.seqpurge.fastq.gz -min_len 50 -threads 20 -summary /home/ravwoy/WORK2/BoczkowskaM/DATA/trimmed_rav/$G.seqourge.output; done

cat raw_files.pairs | while read R1 R2 G; do cutadapt -m 50 -j 24 -q 20 -o /home/ravwoy/WORK2/BoczkowskaM/DATA/trimmed_rav/$R1.cutadapt.fastq.gz -p /home/ravwoy/WORK2/BoczkowskaM/DATA/trimmed_rav/$R2.cutadapt.fastq.gz -z -u 8 -U 8 -a AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA -O 30 /home/ravwoy/WORK2/BoczkowskaM/DATA/trimmed/$R1.fq.gz /home/ravwoy/WORK2/BoczkowskaM/DATA/trimmed/$R2.fq.gz > report_cutadapt.$G.trim2.out; done

#Trinity assembling

Trinity --seqType fq --max_memory 120G --SS_lib_type FR --left /home/ravwoy/WORK2/BoczkowskaM/DATA/trimmed_rav/137-1_1.cutadapt.fastq.gz,/home/ravwoy/WORK2/BoczkowskaM/DATA/trimmed_rav/137-2_1.cutadapt.fastq.gz,/home/ravwoy/WORK2/BoczkowskaM/DATA/trimmed_rav/137-3_1.cutadapt.fastq.gz,/home/ravwoy/WORK2/BoczkowskaM/DATA/trimmed_rav/140-1_1.cutadapt.fastq.gz,/home/ravwoy/WORK2/BoczkowskaM/DATA/trimmed_rav/140-2_1.cutadapt.fastq.gz,/home/ravwoy/WORK2/BoczkowskaM/DATA/trimmed_rav/140-3_1.cutadapt.fastq.gz,/home/ravwoy/WORK2/BoczkowskaM/DATA/trimmed_rav/181-1_1.cutadapt.fastq.gz,/home/ravwoy/WORK2/BoczkowskaM/DATA/trimmed_rav/181-2_1.cutadapt.fastq.gz,/home/ravwoy/WORK2/BoczkowskaM/DATA/trimmed_rav/181-3_1.cutadapt.fastq.gz --right /home/ravwoy/WORK2/BoczkowskaM/DATA/trimmed_rav/137-1_2.cutadapt.fastq.gz,/home/ravwoy/WORK2/BoczkowskaM/DATA/trimmed_rav/137-2_2.cutadapt.fastq.gz,/home/ravwoy/WORK2/BoczkowskaM/DATA/trimmed_rav/137-3_2.cutadapt.fastq.gz,/home/ravwoy/WORK2/BoczkowskaM/DATA/trimmed_rav/140-1_2.cutadapt.fastq.gz,/home/ravwoy/WORK2/BoczkowskaM/DATA/trimmed_rav/140-2_2.cutadapt.fastq.gz,/home/ravwoy/WORK2/BoczkowskaM/DATA/trimmed_rav/140-3_2.cutadapt.fastq.gz,/home/ravwoy/WORK2/BoczkowskaM/DATA/trimmed_rav/181-1_2.cutadapt.fastq.gz,/home/ravwoy/WORK2/BoczkowskaM/DATA/trimmed_rav/181-2_2.cutadapt.fastq.gz,/home/ravwoy/WORK2/BoczkowskaM/DATA/trimmed_rav/181-3_2.cutadapt.fastq.gz --CPU 24 --output trinityFR_out

#Salmon mapping of RNA-seq reads to Trinity transcripts

path=$(pwd); salmon index -p 12 -t $path/trinityFR_out/Trinity.fasta -i $path/TrinityFR_trans_index; 

path=$(pwd); mkdir $path/quants_TrinityFR_lISF; cat $path/raw_files.pairs | while read R1 R2 G; do salmon quant --seqBias --gcBias -i $path/TrinityFR_trans_index -l ISF -1 $path/DATA/trimmed_rav/$R1.cutadapt.fastq.gz -2 $path/DATA/trimmed_rav/$R2.cutadapt.fastq.gz -p 12 --validateMappings -o $path/quants_TrinityFR_lISF/$G.quant; done

#Mapping of RNA-seq reads to genome with HISAT2, assembly of transcripts with Stringtie and counting with BAllgown

mkdir mappings; path=$(pwd); cat $path/raw_files.pairs | while read R1 R2 G; do hisat2 -p 22 --rna-strandness FR --dta -x Hordeum_vulgare.IBSC_v2.index -1 DATA/trimmed_rav/$R1.cutadapt.fastq.gz -2 DATA/trimmed_rav/$R2.cutadapt.fastq.gz -S mappings/$G.sam 2>&1 | tee mappings/$G.out; samtools sort -@ 22 -o mappings/$G.bam mappings/$G.sam; rm mappings/$G.sam; done; wait

mkdir hisat2_assemblies; path=$(pwd); cat $path/raw_files.pairs | while read R1 R2 G; do echo start.assemblies.$G; stringtie -p 22 -o hisat2_assemblies/$G.gtf -l $G mappings/$G.bam; echo end.assemblies.$G; done; wait

echo "NEXT STAGE"; ls -1 hisat2_assemblies/*.gtf > hisat2_assemblies/mergelist.txt

stringtie --merge -p 22 -o hisat2_assemblies/stringtie_merged.gtf hisat2_assemblies/mergelist.txt

#gffcompare -r Hordeum_vulgare.IBSC_v2.48.gtf -G -o merged hisat2_assemblies/stringtie_merged.gtf

mkdir ballgown; path=$(pwd); cat $path/raw_files.pairs | while read R1 R2 G; do echo start.ballgown.$G; stringtie -e -B -p 22 -G hisat2_assemblies/stringtie_merged.gtf -o ballgown/$G/$G.gtf mappings/$G.bam; echo end.ballgown.$G; done; wait; echo "FINISHED"

#Next steps in R.

