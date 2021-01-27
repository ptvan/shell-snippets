############# 
# SAMTOOLS
#############

# get help
samtools view -?

# count the number records (unmapped and unmapped reads)
samtools view -c filename.bam

# convert SAM to BAM
samtools view -S -b sample.sam > sample.BAM

# view the first 5 alignments
samtools view -X sample.sorted.bam | head -n 5

# sort by alignment coordinates
samtools sort sample.bam -o sample.sorted.bam

# index, which allows fast extraction by region
samtools index sample.sorted.bam

# extract 33rd megabase of the chromosome 1, then count alignments
samtools view sample.sorted.bam 1:33000000-34000000 | wc -l

#############
# BCFTOOLS
#############

# convert gVCF to VCF
bcftools --gvcf2vcf gVCF_file.vcf filename.vcf

# merge multi-sample VCFs
bcftools merge -Ob -o output.bcf sampleA.bcf sampleB.bcf

# printing variants on a particular region:
bcftools view -r chr20:1-200000 -s NA20818,NA20819 filename.vcf.gz

# printing only specific samples
bcftools view -s NA20818,NA20819 filename.vcf.gz

# printing out only the chr info:
bcftools query -f '%CHROM\n' filename.vcf

# printing out snps from file:
bcftools view -v snps filename.vcf.gz 

# printing out only multiallelic snps:
bcftools view -m3 -v snps filename.vcf.gz

# rename samples; samplenames.txt file has the following format:
# oldsamplename newsamplename
bcftools reheader -s samplenames.txt oldfile.vcf.gz -o newfile.vcf.gz

## using plugins
# install plugins
export BCFTOOLS_PLUGINS=~/bin/bcftools-1.6/plugins/

# using tag2tag to convert from PL to GL
bcftools +tag2tag in.vcf -- -r --pl-to-gl

# removing INFO field from VCF
bcftools annotate --remove INFO file.vcf.gz


#############
# BEDTOOLS
#############

# extracting promoters from a mouse genome

# >head -n4 genes.bed
# chr1    134212701    134230065    Nuak2    8    +
# chr1    134212701    134230065    Nuak2    7    +
# chr1    33510655    33726603    Prim2,    14    -
# hr1    25124320    25886552    Bai3,    31    -

bedtools flank -i genes.bed -g mm9.chromsizes -l 2000 -r 0 -s > genes.2kb.promoters.bed

bedtools getfasta -fi mm9.fa -bed genes.2kb.promoters.bed -fo genes.2kb.promoters.bed.fa

#############
# TABIX
#############

# compare performance with samtools
time tabix NA18553.chrom11.ILLUMINA.bwa.CHB.low_coverage.20120522.sam.gz 11:60000-5000000 | wc
time samtools view -F 4 -L test.bed NA18553.chrom11.ILLUMINA.bwa.CHB.low_coverage.20120522.bam | wc