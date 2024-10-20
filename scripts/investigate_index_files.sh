#To summarize the unique sequences in an I1 file
cd fastqs
zcat Ratner_CTCF-7_SIC_934_SIC2_Ratner_196_CGTATCTCA_AATACTAATA_S2_I1_001.fastq.gz | awk 'NR % 4 == 2' | grep -v N | sort | uniq -c | awk '{print $1, $2}' | sort -n


