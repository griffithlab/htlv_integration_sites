### Example of an HTLV-1 integration site analysis

Short hand sample names: “CTCF-7, CTCF-8, P12-10B, P12-14"

#### High level questions

What are the experimental details here?  Humanized mice (humanized how?) are infected with different strains of HTLV-1?  Four different strains here?  And we are looking for genome integrations in mice cells? or human cells?   How were the cells obtained for genomic DNA isolation (is this just from blood?).  The goal here is to identify the viral integration sites and quantify them to assess clonality?  Is there an expectation for degree of clonality we might observe.  Are we expecting to see many different unique integration sites in each sample?

Brief answers:
Genomic DNA was isolated from humanized mouse spleen that was infected with HTLV-1 p12(wt control) or CTCF mutant virus. p12-10B, p12-14 and CTCF-7, CTCF-8 are mouse ID numbers. We want to quantify the viral integration sites in infected human T cells to assess clonality. I expect to see many unique integration sites but don’t know what kinds of clonality that would be observed.

CD34+ cells were injected in liver at 1d of life. Infected with HTLV. 2 strains – p12 and CTCF.  Analysis focused on human cell DNA. Samples were obtained from spleen. The goal here is to identify the viral integration sites and quantify them to assess clonality. Using Gini index value. 

#### Download the data

```bash
sh /storage1/fs1/mgriffit/Active/griffithlab/adhoc/ratner_p01/htlv_integration_sites/git/htlv_integration_sites/scripts/download_raw_data.sh

```

 
#### Investigate the four supplies possible integration characteristic sequences:
TTAGTACACA / AATCATGTGT
TGACAATGAC / ACTGTTACTG

```bash
export FASTQ_NAMES=("Ratner_CTCF-7_SIC_934_SIC2_Ratner_196_CGTATCTCA_AATACTAATA_S2_" "Ratner_CTCF-8_SIC_935_SIC2_Ratner_196_GTCCTGCCG_AATACTAATA_S3_" "Ratner_P12-10B_SIC_936_SIC2_Ratner_196_CCGGGACAC_AATACTAATA_S4_" "Ratner_P12-14_SIC_937_SIC2_Ratner_196_GGCTGGGAT_AATACTAATA_S5_")
export PAIRS=("R1" "R2")
export SEQS=("TTAGTACACA" "AATCATGTGT" "TGACAATGAC" "ACTGTTACTG")

for FASTQ_NAME in "${FASTQ_NAMES[@]}"; do
  for PAIR in "${PAIRS[@]}"; do
    for SEQ in "${SEQS[@]}"; do
      ANSWER=$(zcat fastqs/${FASTQ_NAME}${PAIR}_001.fastq.gz | awk 'NR % 4 == 2' | grep $SEQ | wc -l)
      echo "$FASTQ_NAME $PAIR $SEQ $ANSWER"
    done
  done
done
```

Based on this analysis it seems that for these data in the RAW read sequences we only really see the "TTAGTACACA" sequence and only in Read 1 files

#### Create unique read lists of these read identities and store them for later use

```bash
for FASTQ_NAME in "${FASTQ_NAMES[@]}"; do
    echo -e "\nProcessing FASTQ: $FASTQ_NAME (R1 only)"
    SAMPLE=$(echo $FASTQ_NAME | awk -F_ '{print $2"_"$3"_"$4"_"$5}')
    echo "Will name output using sample name: $SAMPLE"
    zcat fastqs/${FASTQ_NAME}R1_001.fastq.gz | awk 'NR % 4 == 1 {read_name = substr($1, 2)} NR % 4 == 2 {print read_name, $0}' | grep -P 'TTTAGTACACA' | cut -f 1 -d ' ' | sort | uniq > readlists/${SAMPLE}_ltr_integration_seq_read_ids.txt
done
```

#### Details of reference sequences and alignments produced

GRCh38 reference copied from: `/storage1/fs1/bga/Active/gmsroot/gc2560/core/GRC-human-build38_human_95_38_U2AF1_fix/all_sequences.fa`
HTLV-1 reference obtained from: `https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/863/585/GCF_000863585.1_ViralProj15434/GCF_000863585.1_ViralProj15434_genomic.fna.gz`

GRCh38 and HTLV-1 references catted together and BWA index and alignment done with BWA version 0.7.17-r1198-dirty

```bash
/usr/local/bwa/bwa mem -K 20000000 -t 8 -Y GRCh38+HTLV-1.fa Ratner_CTCF-7_SIC_934_SIC2_Ratner_196_CGTATCTCA_AATACTAATA_S2_R1_001.fastq.gz Ratner_CTCF-7_SIC_934_SIC2_Ratner_196_CGTATCTCA_AATACTAATA_S2_R2_001.fastq.gz | samtools view -o CTCF-7_SIC_934_SIC2.bam -Shb /dev/stdin
/usr/local/bwa/bwa mem -K 20000000 -t 8 -Y GRCh38+HTLV-1.fa Ratner_CTCF-8_SIC_935_SIC2_Ratner_196_GTCCTGCCG_AATACTAATA_S3_R1_001.fastq.gz Ratner_CTCF-8_SIC_935_SIC2_Ratner_196_GTCCTGCCG_AATACTAATA_S3_R2_001.fastq.gz | samtools view -o CTCF-8_SIC_935_SIC2.bam -Shb /dev/stdin
/usr/local/bwa/bwa mem -K 20000000 -t 8 -Y GRCh38+HTLV-1.fa Ratner_P12-10B_SIC_936_SIC2_Ratner_196_CCGGGACAC_AATACTAATA_S4_R1_001.fastq.gz Ratner_P12-10B_SIC_936_SIC2_Ratner_196_CCGGGACAC_AATACTAATA_S4_R2_001.fastq.gz | samtools view -o P12-10B_SIC_936_SIC2.bam -Shb /dev/stdin
/usr/local/bwa/bwa mem -K 20000000 -t 8 -Y GRCh38+HTLV-1.fa Ratner_P12-14_SIC_937_SIC2_Ratner_196_GGCTGGGAT_AATACTAATA_S5_R1_001.fastq.gz Ratner_P12-14_SIC_937_SIC2_Ratner_196_GGCTGGGAT_AATACTAATA_S5_R2_001.fastq.gz | samtools view -o P12-14_SIC_937_SIC2.bam -Shb /dev/stdin
```

#### Duplicate marking step
alignments converted to bam, sorted, and indexed with samtools version 1.11
duplicates marked with picard version 2.22.8

```bash
java -Xmx16g -jar /usr/local/picard.jar MarkDuplicates I=CTCF-7_SIC_934_SIC2.sorted.bam O=CTCF-7_SIC_934_SIC2.markedsorted.bam M=CTCF-7_SIC_934_SIC2.metrics
java -Xmx16g -jar /usr/local/picard.jar MarkDuplicates I=CTCF-8_SIC_935_SIC2.sorted.bam O=CTCF-8_SIC_935_SIC2.markedsorted.bam M=CTCF-8_SIC_935_SIC2.metrics
java -Xmx16g -jar /usr/local/picard.jar MarkDuplicates I=P12-10B_SIC_936_SIC2.sorted.bam O=P12-10B_SIC_936_SIC2.markedsorted.bam M=P12-10B_SIC_936_SIC2.metrics
java -Xmx16g -jar /usr/local/picard.jar MarkDuplicates I=P12-14_SIC_937_SIC2.sorted.bam O=P12-14_SIC_937_SIC2.markedsorted.bam M=P12-14_SIC_937_SIC2.metrics

```

#### LTR integration site read filtering of BAM

Produce a version of the duplicate marked BAM that is limited to only those alignments involving reads that contained the characterstic integration site sequence (TTTAGTACACA|TGTGTACTAAA) identified above

```bash
export SAMPLES=("CTCF-7_SIC_934_SIC2" "CTCF-8_SIC_935_SIC2" "P12-10B_SIC_936_SIC2" "P12-14_SIC_937_SIC2")
rm -f tmp/*

for SAMPLE in "${SAMPLES[@]}"; do
    echo -e "\nProducing integration site read list filtered BAM for $SAMPLE"
    echo "samtools view -H bams/${SAMPLE}.markedsorted.bam > tmp/${SAMPLE}.markedsorted.bam.header"
    samtools view -H bams/${SAMPLE}.markedsorted.bam > tmp/${SAMPLE}.markedsorted.bam.header

    echo "samtools view bams/${SAMPLE}.markedsorted.bam | grep -F -f readlists/${SAMPLE}_ltr_integration_seq_read_ids.txt > tmp/${SAMPLE}.markedsorted_ltr_integration_seq_reads.sam"
    samtools view bams/${SAMPLE}.markedsorted.bam | grep -F -f readlists/${SAMPLE}_ltr_integration_seq_read_ids.txt > tmp/${SAMPLE}.markedsorted_ltr_integration_seq_reads.sam

    echo "cat tmp/${SAMPLE}.markedsorted.bam.header tmp/${SAMPLE}.markedsorted_ltr_integration_seq_reads.sam | samtools view -Sb - > tmp/${SAMPLE}.markedsorted_ltr_integration_seq_reads.bam"
    cat tmp/${SAMPLE}.markedsorted.bam.header tmp/${SAMPLE}.markedsorted_ltr_integration_seq_reads.sam | samtools view -Sb - > tmp/${SAMPLE}.markedsorted_ltr_integration_seq_reads.bam

    echo "samtools sort tmp/${SAMPLE}.markedsorted_ltr_integration_seq_reads.bam -o bams/${SAMPLE}.markedsorted_ltr_integration_seq_reads.bam"
    samtools sort tmp/${SAMPLE}.markedsorted_ltr_integration_seq_reads.bam -o bams/${SAMPLE}.markedsorted_ltr_integration_seq_reads.bam

    echo "samtools index bams/${SAMPLE}.markedsorted_ltr_integration_seq_reads.bam"
    samtools index bams/${SAMPLE}.markedsorted_ltr_integration_seq_reads.bam
done

```

#### identify reads involving a primary or supplementary alignment to the virus sequence (`NC_001436.1`)

```bash
export SAMPLES=("CTCF-7_SIC_934_SIC2" "CTCF-8_SIC_935_SIC2" "P12-10B_SIC_936_SIC2" "P12-14_SIC_937_SIC2")

for SAMPLE in "${SAMPLES[@]}"; do
    echo -e "\nProcessing sample: $SAMPLE"
    echo "Obtaining reads with supplemetary alignments to the virus"
    echo "samtools view -f 2048 bams/${SAMPLE}.markedsorted_ltr_integration_seq_reads.bam 'NC_001436.1' | cut -f 1 | sort | uniq > readlists/${SAMPLE}_supplementary_virus_hit_read_ids.txt"
    samtools view -f 2048 bams/${SAMPLE}.markedsorted_ltr_integration_seq_reads.bam 'NC_001436.1' | cut -f 1 | sort | uniq > readlists/${SAMPLE}_supplementary_virus_hit_read_ids.txt
    wc -l readlists/${SAMPLE}_supplementary_virus_hit_read_ids.txt

    echo -e "\nObtaining reads with primary alignments to the virus"
    echo "samtools view -F 256 bams/${SAMPLE}.markedsorted_ltr_integration_seq_reads.bam 'NC_001436.1' | cut -f 1 | sort | uniq > readlists/${SAMPLE}_primary_virus_hit_read_ids.txt"
    samtools view -F 256 bams/${SAMPLE}.markedsorted_ltr_integration_seq_reads.bam 'NC_001436.1' | cut -f 1 | sort | uniq > readlists/${SAMPLE}_primary_virus_hit_read_ids.txt
    wc -l readlists/${SAMPLE}_primary_virus_hit_read_ids.txt

    echo -e "\nCreating a unique list of reads with either supplementary or primary alignments to the virus"
    echo "cat readlists/${SAMPLE}_supplementary_virus_hit_read_ids.txt readlists/${SAMPLE}_primary_virus_hit_read_ids.txt | sort | uniq > readlists/${SAMPLE}_virus_hit_read_ids.txt"
    cat readlists/${SAMPLE}_supplementary_virus_hit_read_ids.txt readlists/${SAMPLE}_primary_virus_hit_read_ids.txt | sort | uniq > readlists/${SAMPLE}_virus_hit_read_ids.txt
    wc -l readlists/${SAMPLE}_virus_hit_read_ids.txt
    echo -e "\n"
done
```


#### Use the viral alignment read list to produce filtered BAM files with only those reads that have such alignments

```bash
for SAMPLE in "${SAMPLES[@]}"; do
    echo -e "\nProducing viral read list filtered BAM for $SAMPLE"
    echo "samtools view -H bams/${SAMPLE}.markedsorted_ltr_integration_seq_reads.bam > tmp/${SAMPLE}.markedsorted_ltr_integration_seq_reads.bam.header"
    samtools view -H bams/${SAMPLE}.markedsorted_ltr_integration_seq_reads.bam > tmp/${SAMPLE}.markedsorted_ltr_integration_seq_reads.bam.header

    echo "samtools view bams/${SAMPLE}.markedsorted_ltr_integration_seq_reads.bam | grep -F -f readlists/${SAMPLE}_virus_hit_read_ids.txt > tmp/${SAMPLE}.markedsorted.viralreads.sam"
    samtools view bams/${SAMPLE}.markedsorted_ltr_integration_seq_reads.bam | grep -F -f readlists/${SAMPLE}_virus_hit_read_ids.txt > tmp/${SAMPLE}.markedsorted.viralreads.sam

    echo "cat tmp/${SAMPLE}.markedsorted_ltr_integration_seq_reads.bam.header tmp/${SAMPLE}.markedsorted.viralreads.sam | samtools view -Sb - > bams/${SAMPLE}.markedsorted_with_hits_to_viral.bam"
    cat tmp/${SAMPLE}.markedsorted_ltr_integration_seq_reads.bam.header tmp/${SAMPLE}.markedsorted.viralreads.sam | samtools view -Sb - > bams/${SAMPLE}.markedsorted_with_hits_to_viral.bam
done
```

**NOTE**: It is possible that the above approach of limited to reads with a viral alignment could be too strict.
A read may correspond to an integration event, have the characteristic LTR sequence, but not produce an alignment to the virus genome 
We should try the analysis, with and without this requirement and gauge impact

Use bedtools (v2.25.0) to create bed representations of the BAM alignments to facilitate integration site counting.
At the same time apply filters to: require alignment, remove duplicates, prevent counting on the virus seq itself

Do this two ways: (1) with the marked-duplicate BAM, (2) with the BAM created from marked-duplicate BAM that also limits to viral hit reads produced above

##### (1) with the marked-duplicate BAM

```bash
rm -f tmp/*
for SAMPLE in "${SAMPLES[@]}"; do
    echo -e "\nProducing integration site counts $SAMPLE"
    echo "samtools view -H bams/${SAMPLE}.markedsorted_ltr_integration_seq_reads.bam > tmp/${SAMPLE}.markedsorted_ltr_integration_seq_reads.bam.header"
    samtools view -H bams/${SAMPLE}.markedsorted_ltr_integration_seq_reads.bam > tmp/${SAMPLE}.markedsorted_ltr_integration_seq_reads.bam.header
    echo "samtools view -f 1 -F 1024 -q 20 bams/${SAMPLE}.markedsorted_ltr_integration_seq_reads.bam > tmp/${SAMPLE}.markedsorted_filtered.sam"
    samtools view -f 1 -F 1024 -q 20 bams/${SAMPLE}.markedsorted_ltr_integration_seq_reads.bam > tmp/${SAMPLE}.markedsorted_filtered.sam
    echo "cat tmp/${SAMPLE}.markedsorted_ltr_integration_seq_reads.bam.header tmp/${SAMPLE}.markedsorted_filtered.sam | samtools view -Sb - > tmp/${SAMPLE}.markedsorted_filtered.bam"
    cat tmp/${SAMPLE}.markedsorted_ltr_integration_seq_reads.bam.header tmp/${SAMPLE}.markedsorted_filtered.sam | samtools view -Sb - > tmp/${SAMPLE}.markedsorted_filtered.bam
    echo "samtools sort tmp/${SAMPLE}.markedsorted_filtered.bam -o bams/${SAMPLE}.markedsorted_filtered.bam"
    samtools sort tmp/${SAMPLE}.markedsorted_filtered.bam -o bams/${SAMPLE}.markedsorted_filtered.bam
    echo "samtools index bams/${SAMPLE}.markedsorted_filtered.bam"
    samtools index bams/${SAMPLE}.markedsorted_filtered.bam
    echo "bedtools bamtobed -i bams/${SAMPLE}.markedsorted_filtered.bam > beds/v2/${SAMPLE}.markedsorted_filtered.bed"
    bedtools bamtobed -i bams/${SAMPLE}.markedsorted_filtered.bam > beds/v2/${SAMPLE}.markedsorted_filtered.bed
    echo "bedtools merge -i beds/v2/${SAMPLE}.markedsorted_filtered.bed -c 1 -o count | grep -v 'NC_001436.1' > beds/v2/${SAMPLE}.markedsorted_filtered_merged.bed"
    bedtools merge -i beds/v2/${SAMPLE}.markedsorted_filtered.bed -c 1 -o count | grep -v 'NC_001436.1' > beds/v2/${SAMPLE}.markedsorted_filtered_merged.bed
    echo "cp beds/v2/${SAMPLE}.markedsorted_filtered_merged.bed counts/v2/${SAMPLE}.markedsorted_filtered_merged.bed.tsv"
    cp beds/v2/${SAMPLE}.markedsorted_filtered_merged.bed counts/v2/${SAMPLE}.markedsorted_filtered_merged.bed.tsv
done
```

##### (2) with the BAM created from marked-duplicate BAM that also limits to viral hit reads produced above

```bash
rm -f tmp/*
for SAMPLE in "${SAMPLES[@]}"; do
    echo -e "\nProducing integration site counts $SAMPLE"
    echo "samtools view -H bams/${SAMPLE}.markedsorted_with_hits_to_viral.bam > tmp/${SAMPLE}.markedsorted_with_hits_to_viral.bam.header"
    samtools view -H bams/${SAMPLE}.markedsorted_with_hits_to_viral.bam > tmp/${SAMPLE}.markedsorted_with_hits_to_viral.bam.header
    echo "samtools view -f 1 -F 1024 -q 20 bams/${SAMPLE}.markedsorted_with_hits_to_viral.bam > tmp/${SAMPLE}.markedsorted_with_hits_to_viral_filtered.sam"
    samtools view -f 1 -F 1024 -q 20 bams/${SAMPLE}.markedsorted_with_hits_to_viral.bam > tmp/${SAMPLE}.markedsorted_with_hits_to_viral_filtered.sam
    echo "cat tmp/${SAMPLE}.markedsorted_with_hits_to_viral.bam.header tmp/${SAMPLE}.markedsorted_with_hits_to_viral_filtered.sam | samtools view -Sb - > tmp/${SAMPLE}.markedsorted_with_hits_to_viral_filtered.bam"
    cat tmp/${SAMPLE}.markedsorted_with_hits_to_viral.bam.header tmp/${SAMPLE}.markedsorted_with_hits_to_viral_filtered.sam | samtools view -Sb - > tmp/${SAMPLE}.markedsorted_with_hits_to_viral_filtered.bam
    echo "samtools sort tmp/${SAMPLE}.markedsorted_with_hits_to_viral_filtered.bam -o bams/${SAMPLE}.markedsorted_with_hits_to_viral_filtered.bam"
    samtools sort tmp/${SAMPLE}.markedsorted_with_hits_to_viral_filtered.bam -o bams/${SAMPLE}.markedsorted_with_hits_to_viral_filtered.bam
    echo "samtools index bams/${SAMPLE}.markedsorted_with_hits_to_viral_filtered.bam"
    samtools index bams/${SAMPLE}.markedsorted_with_hits_to_viral_filtered.bam
    echo "bedtools bamtobed -i bams/${SAMPLE}.markedsorted_with_hits_to_viral_filtered.bam > beds/v2/${SAMPLE}.markedsorted_with_hits_to_viral_filtered.bed"
    bedtools bamtobed -i bams/${SAMPLE}.markedsorted_with_hits_to_viral_filtered.bam > beds/v2/${SAMPLE}.markedsorted_with_hits_to_viral_filtered.bed
    echo "bedtools merge -i beds/v2/${SAMPLE}.markedsorted_with_hits_to_viral_filtered.bed -c 1 -o count | grep -v 'NC_001436.1' > beds/v2/${SAMPLE}.markedsorted_with_hits_to_viral_filtered_merged.bed"
    bedtools merge -i beds/v2/${SAMPLE}.markedsorted_with_hits_to_viral_filtered.bed -c 1 -o count | grep -v 'NC_001436.1' > beds/v2/${SAMPLE}.markedsorted_with_hits_to_viral_filtered_merged.bed
    echo "cp beds/v2/${SAMPLE}.markedsorted_with_hits_to_viral_filtered_merged.bed counts/v2/${SAMPLE}.markedsorted_with_hits_to_viral_filtered_merged.bed.tsv"
    cp beds/v2/${SAMPLE}.markedsorted_with_hits_to_viral_filtered_merged.bed counts/v2/${SAMPLE}.markedsorted_with_hits_to_viral_filtered_merged.bed.tsv
done
```

#### Create single file with all counts for all samples to facilitate creation of visualizations

```bash
rm -f tmp/*
for SAMPLE in "${SAMPLES[@]}"; do
    awk -v sample="$SAMPLE" '{print $0 "\t" sample}' counts/v2/${SAMPLE}.markedsorted_with_hits_to_viral_filtered_merged.bed.tsv > tmp/${SAMPLE}.markedsorted_with_hits_to_viral_filtered_merged.bed.tsv
done
echo -e "chromosome\tstart_pos\tend_pos\tcount\tsample" > tmp/header.tsv
cat tmp/header.tsv tmp/*markedsorted_with_hits_to_viral_filtered_merged.bed.tsv > counts/v2/ALL.markedsorted_with_hits_to_viral_filtered_merged.bed.tsv
```



