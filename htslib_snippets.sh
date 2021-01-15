############# 
# SAMTOOLS
#############

# get help
samtools view -?

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

# rename samples; samplenames.txt file has the following format:
# oldsamplename newsamplename
bcftools reheader -s samplenames.txt oldfile.vcf.gz -o newfile.vcf.gz

#############
# TABIX
#############

# compare performance with samtools
time tabix NA18553.chrom11.ILLUMINA.bwa.CHB.low_coverage.20120522.sam.gz 11:60000-5000000 | wc
time samtools view -F 4 -L test.bed NA18553.chrom11.ILLUMINA.bwa.CHB.low_coverage.20120522.bam | wc