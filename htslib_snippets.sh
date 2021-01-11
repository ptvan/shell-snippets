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