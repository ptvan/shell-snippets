######################
# FASTA / FASTQ files
######################

# build a local BLAST database
makeblastdb -in ref_viruses.fa -dbtype nucl -parse_seqids -out ref_viruses_db

# BLAST against database, output in tabular form
blastn -db db/ref_viruses_rep_genomes -query inputs.fa -outfmt 6

# find barcodes that appear most frequently
zcat sample.fastq.gz | awk 'NR % 4 == 2 {print;}' | sort | uniq -c | sort -n -r | less

# FASTQC supports multi-threaded operation
fastqc -o ./output_dir -t 10 *.fastq.gz

# randomly select 100 sequences from paired-end FASTQs,same seed to keep pairing
seqtk sample -s100 read1.fq 10000 > sub1.fq
seqtk sample -s100 read2.fq 10000 > sub2.fq

# mask regions in BED to lower case
seqtk seq -M region.bed in.fa > out.fa

# trim low-quality bases (based by Phred scores)
seqtk trimfq in.fastq > out.fastq

#############
# BED files
#############
# eliminate identical entries from a BED file
sort -k1,1 -k2,2n input.bed | rev | uniq -f1 | rev

## extracting promoters from a mouse genome
# >head -n4 genes.bed
# chr1    134212701    134230065    Nuak2    8    +
# chr1    134212701    134230065    Nuak2    7    +
# chr1    33510655    33726603    Prim2,    14    -
# hr1    25124320    25886552    Bai3,    31    -
bedtools flank -i genes.bed -g mm9.chromsizes -l 2000 -r 0 -s > genes.2kb.promoters.bed
bedtools getfasta -fi mm9.fa -bed genes.2kb.promoters.bed -fo genes.2kb.promoters.bed.fa

# converting BAM to BED
bedtools bamtobed -i input_file.bam 

# list coverage of each target
bedtools coverage -hist -abam my_data.bam -b my_targets.bed

# get parts of `regions.bed` that are *NOT* covered by `mydata.bam`
bedtools genomecov -ibam mydata.bam -bga | awk '$4==0' | bedtools intersect -a regions.bed -b - > matches.txt


#############
# BAM files
#############

# interconvert SAM <-> BAM
samtools view -S -b sample.sam > sample.bam
samtools view -bT reference.fa sample.sam > test.bam
samtools view -h -o out.sam in.bam

# count number of records (includes both unmapped & unmapped reads)
samtools view -c filename.bam

# count number of reads
samtools idxstats filename.bam | awk '{s+=$3+$4} END {print s}'

# count number of _mapped_ reads
samtools idxstats filename.bam | awk '{s+=$3} END {print s}'
samtools view -F 4 -c filename.bam

# count number of _UNMAPPED_ reads
samtools view -f 4 -c filename.bam

# view the first 5 alignments
samtools view sample.sorted.bam | head -n 5

# view secondary alignments
# more information about SAM tags at
# https://broadinstitute.github.io/picard/explain-flags.html
samtools view -cf 0x100 sample.bam

# view supplementary alignments
samtools view -cf 0x800 sample.bam

# filter BAM to a region, but keep reads paired even if one read
# is outside target interval
samtools view -L regions.bed --fetch-pairs sample.bam

# simple statistics
samtools flagstat sample.bam

# sort by alignment coordinates
samtools sort sample.bam -o sample.sorted.bam
sambamba sort sample.bam -o sample.sorted.bam

# index, which allows fast extraction by region
samtools index sample.sorted.bam
sambamba index sample.sorted.bam

# extract read group information
samtools view -H sample.sorted.bam | grep '^@RG'

# extract bases 100 - 200 of chromosome 1
sambamba slice sample.bam chr1:100-200 sample.chr1.bam

# extract 33rd megabase of the chromosome 1, then count alignments
samtools view sample.sorted.bam 1:33000000-34000000 | wc -l

# extract only the first read from paired read BAMs
samtools view -h -f 0x0040 paired_ends.bam > first_reads.sam

# extract regions from Kraken contamination BAM
samtools view -h sample.consensus.kraken_annotate.bosTau6.bam | grep -v "^@" | less -S

# extract reads from BAM that support a variant call in corresponding VCF
samtools mpileup -f reference.fa -r chr22:425236-425236 alignments.bam |
    cut -f 5 | tr '[a-z]' '[A-Z]' | fold -w 1 | sort | uniq -c

# checking to see if reads are sorted
samtools view -H 5_110118_FC62VT6AAXX-hg18-unsort.bam
### @HD    VN:1.0    SO:unsorted
samtools view -H 5_110118_FC62VT6AAXX-hg18-sort.bam
### @HD    VN:1.0    SO:coordinate

# filter out alignments with MAPQ < 2
samtools view -bq 2 original.bam > quality2_filtered.bam

# remove "Chr" prefix in header
samtools reheader -c 'perl -pe "s/^(@SQ.*)(\tSN:)Chr/\$1\$2/"' in.bam

# extract soft-clipped bases
samtools view input.bam | awk '$6 ~ /S/{print $1}' | sort -k1,1 | uniq > soft-clipped-names.txt
samtools view -hb -o output.bam -N soft-clipped-names.txt input.bam

# mark duplicate reads
sambamba markdup sample.bam sample.nodups.bam

#############
# VCF files
#############

## using plugins with bcftools
# install plugins
export BCFTOOLS_PLUGINS=~/bin/bcftools-1.6/plugins/

# calling variants
bcftools mpileup -Ou -f reference.fa alignments.bam | bcftools call -mv -Ob -o calls.bcf

# viewing VCF fields using `verticalize` (https://github.com/lindenb/verticalize):
cat my_file.vcf | grep -vE "^##" | verticalize

# convert gVCF to VCF
bcftools --gvcf2vcf gVCF_file.vcf filename.vcf

# merge multi-sample VCFs
bcftools merge -Ob -o output.bcf sampleA.bcf sampleB.bcf

# printing only specific samples
bcftools view -s NA20818,NA20819 filename.vcf.gz

# printing variants on a particular region:
bcftools view -r chr20:1-200000 -s NA20818,NA20819 filename.vcf.gz

# printing out only the chr info:
bcftools query -f '%CHROM\n' filename.vcf

# printing out snps from file:
bcftools view -v snps filename.vcf.gz 

# printing out only multiallelic snps:
bcftools view -m3 -v snps filename.vcf.gz

# printing out missing (uncalled) genotypes:
bcftools view -u filename.vcf.gz -o missing_genotypes.vcf.gz -Oz

# drop individual genotype information
bcftools view -G filename.vcf.gz

# filtering using one of the INFO annotations (IDV)
bcftools filter -sFilterName -e'IDV<5' filename.vcf

# split multiallelic SNPs + INDELs into several records
bcftools norm -m -any filename.vcf.gz -o normalized.vcf.gz -Oz

# rename samples; samplenames.txt file has the following format:
# oldsamplename newsamplename
bcftools reheader -s samplenames.txt oldfile.vcf.gz -o newfile.vcf.gz

# generate stats
bcftools stats -s filename.vcf > filename.vchk

# plot stats & tweak that figure
plot-vcfstats -p outdir file.vchk
cd outdir && python plot.py && pdflatex summary.tex

# using tag2tag to convert from PL to GL
bcftools +tag2tag in.vcf -- -r --pl-to-gl

# removing INFO field from VCF
bcftools annotate --remove INFO file.vcf.gz

# use split-vep plugin to get output from VCF annotatated with
# Ensembl Variant Effect Predictor (VEP, https://useast.ensembl.org/info/docs/tools/vep/index.html)
bcftools +split-vep test/split-vep.vcf -f '%CHROM:%POS %Consequence\n' -d

#############
# VEP
#############
# installing VEP locally is annoying, also the human database is ~23GB 
./vep -i input.vcf.gz --cached --vcf --fields "Allele,Consequence,Feature_type, Feature" -o output.vcf

# pull from DockerHub instead
docker pull ensemblorg/ensembl-vep
sudo docker run -v /home/ptv/working:/working -it ensemblorg/ensembl-vep 

#############
# SNPEFF
#############
# list supported databases & download the one we want
java -jar snpEff.jar databases
java -jar snpEff.jar download -v GRCh38.14

# annotate using clinVar
java -Xmx8g -jar SnpSift.jar  annotate -v protocols/db/clinvar_00-latest.vcf \
    -stats variant_stats.html input.vcf > output.clinvar.vcf

# annotate using DbSnp
java -Xmx8g -jar SnpSift.jar annotate -dbsnp input.vcf > output.dbSnp.vcf

# further annotate the output above using SnpEff
java -Xmx8g -jar snpEff.jar eff -v GRCh38.14 output.dbSnp.vcf > output.ann.vcf

#############
# MAF files
#############
# convert VCF > MAF using https://github.com/mskcc/vcf2maf
# the same repo also contains `maf2vcf.pl` for MAF > VCF
perl vcf2maf.pl --input-vcf input.vcf --output-maf output.maf

#############
# MUT files
#############

# convert gMUT to MUT by excluding `no_variant` calls
grep -v no_variant input_genomic.mut > output.mut

# compare two MUT files using VIM
vimdiff  <(cut -f1-3,5-12,14-15 first_file.mut) <(cut -f1-3,5-12,14-15 second_file.mut)


#############
# TABIX
#############

# compare performance with samtools
time tabix NA18553.chrom11.ILLUMINA.bwa.CHB.low_coverage.20120522.sam.gz 11:60000-5000000 | wc
time samtools view -F 4 -L test.bed NA18553.chrom11.ILLUMINA.bwa.CHB.low_coverage.20120522.bam | wc


#############
# DEEPTOOLS
#############
# create coverage (in the form of a bigWig/bedGraph) file from a BAM file
bamCoverage -b reads.bam -o coverage.bw --numberOfProcessors 8

# compare ChIP-seq coverages between two different samples
multiBamSummary BED-file --BED selection.bed --bamfiles sample1.bam sample2.bam -out comparison.npz --numberOfProcessors 8

# plotting the results of multiBamSummary
plotCorrelation -in comparison.npz \
--corMethod pearson --skipZeros \
--plotTitle "Pearson Correlation of Average Scores Per Transcript" \
--whatToPlot scatterplot \
-o scatterplot_PearsonCorr_bigwigScores.png   \
--outFileCorMatrix PearsonCorr_bigwigScores.tab

# calculate GC bias for a human sample
computeGCBias -b H3K27Me3.bam --effectiveGenomeSize 2695000000
   --genome genome.2bit -l 200 -freq freq_test.txt
   --region X --biasPlot test_gc.png

# scale mouse and fly GC contents and plot heatmaps for comparison
computeMatrix scale-regions -S mouse_GCcontent.bw -R RefSeq_genes_uniq.bed  -m 10000 -b 3000 -a 3000 -out mouse_matrix.tab.gz
plotHeatmap -m mouse_matrix.tab.gz -out scaledmouseGC.png --colorMap YlGnBu --regionsLabel 'mouse genes'

computeMatrix scale-regions -S fly_GCcontent.bw -R RefSeq_genes_uniq.bed -m 3000 -b 1000 -a 1000 -out fly_matrix.tab.gz
plotHeatmap -m fly_matrix.tab.gz -out scaledflyGC.png --colorMap YlGnBu --regionsLabel 'fly genes'
