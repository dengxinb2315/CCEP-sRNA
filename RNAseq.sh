#!/usr/bin/env bash
#RNA-seq raw results in dir raw/
#conda activate dap
mkdir ../trim
mkdir rawqc
#QC
fastqc *gz -o rawqc/
cd rawqc
multiqc .
cd ..
#make sure the files name format in <ID-repeat_Ln_1.fq.gz>
for i in $(ls *.fq.gz | cut -d '_' -f 1-2 | sed 's/_/_/g' | sort | uniq); do echo ${i}; done
#trim
for i in $(ls *.fq.gz | cut -d '_' -f 1-2 | sed 's/_/_/g' | sort | uniq); do trim_galore -q 20 --paired --phred33 --stringency 3 ${i}_1.fq.gz ${i}_2.fq.gz --gzip -o ../trim --length 8; echo ${i}; done
echo "trim_galore finished"
cd ../trim
#map
#remember to revise the reference
for i in $(ls *.fq.gz | cut -d '_' -f 1-2 | sed 's/_/_/g' | sort | uniq); do hisat2 -p 10 -q -x /media/lu/lu202306/ref/PAO1_genome -1 ${i}_1_val_1.fq.gz -2 ${i}_2_val_2.fq.gz -S ${i}.sam; samtools view -bS -@ 8 ${i}.sam > ${i}.bam; rm ${i}.sam; samtools sort -@ 8 ${i}.bam > ${i}.sorted.bam; rm ${i}.bam; done
echo "bam produced and sorted"
#count
#remember to revise the gff file
#featureCounts -T 10 -t CDS -g locus -a /media/lu/lu202306/ref/PAO1forsRNA.gff -o ../counts.txt *bam
echo "counts got!"
echo "Cheers!"