
module load biokit
cd /cuiwang/nodelete/rawRAD/populations/firstsnpperloci/Ancestral
bwa mem -t 20 /scratch/project_2001259/cuiwang/nodelete/rawRAD/PaEUcurated.fasta oropetium_1.fastq.gz oropetium_2.fastq.gz | samtools sort -@ 20 | samtools view -@20 -F 1284 -o oropetium.sorted.bam
bwa mem -t 20 /scratch/project_2001259/cuiwang/nodelete/rawRAD/PaEUcurated.fasta Sbicolor_1.fastq.gz Sbicolor_2.fastq.gz |  samtools sort -@ 20 | samtools view -@20 -F 1284 -o sbicolor.sorted.bam
bwa mem -t 20 /scratch/project_2001259/cuiwang/nodelete/rawRAD/PaEUcurated.fasta Msisnensis_1.fastq.gz Msisnensis_2.fastq.gz | samtools sort -@ 20 |  samtools view -@20 -F 1284 -o Msisnensis.sorted.bam
bwa mem -t 20 /scratch/project_2001259/cuiwang/nodelete/rawRAD/PaEUcurated.fasta aplinii_1.fastq.gz aplinii_2.fastq.gz | samtools sort -@ 20 |  samtools view -@20 -F 1284 -o aplinii.sorted.bam
bwa mem -t 20 /scratch/project_2001259/cuiwang/nodelete/rawRAD/PaEUcurated.fasta adonax_1.fastq.gz adonax_2.fastq.gz | samtools sort -@ 20 |  samtools view -@20 -F 1284 -o adonax.sorted.bam
ls *.bam > ancs.txt

module load angsd/0.940
angsd -b ancs.txt -doFasta 2 -doCounts 1 -P 8 -out ancestral
module load bedtools
#
module load vcftools
cat ../populations.snps.vcf | vcf-sort -c -p 8 > populations.snps.sorted.vcf
fastaFromBed -bedOut -fi ancestral.fa -bed populations.snps.sorted.vcf | cut -f1,2,4,5,99 > anc_AA.tab
bgzip anc_AA.tab
tabix -s1 -b2 -e2 anc_AA.tab.gz
echo '##INFO=<ID=AA,Number=1,Type=Character,Description="Ancestral allele">' > modif_hdr.txt
bcftools annotate -a anc_AA.tab.gz -c CHROM,POS,REF,ALT,INFO/AA -h modif_hdr.txt -Oz -o output.AA.vcf.gz populations.snps.sorted.all.vcf
bcftools +fixref output.AA.vcf.gz -Oz -o output.AA.swap.vcf.gz -- -f ancestral.fa -m flip 
#do it for all the sites
cd /cuiwang/nodelete/rawRAD/populations/all
ln -s /cuiwang/nodelete/rawRAD/populations/populations.all.vcf .
ln -s /cuiwang/nodelete/rawRAD/populations/firstsnpperloci/Ancestral/ancestral.fa .
cat populations.all.vcf | vcf-sort -p 8 > populations.snps.sorted.all.vcf
fastaFromBed -bedOut -fi ancestral.fa -bed populations.snps.sorted.all.vcf | cut -f1,2,4,5,99 > anc_AA.tab
