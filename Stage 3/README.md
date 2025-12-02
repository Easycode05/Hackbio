# Transcriptomic Profiling of Staphylococcus aureus During Acute vs Chronic Phases of Periprosthetic Joint Infection (PJI)
## Tools Used
-  **fastqc & multiqc** → quality control 
-  **fastp** → read trimming and quality filtering
- **STAR** → alignment to reference genome
- **samtools** → BAM file sorting and indexing
- **featureCounts** → generation of gene count matrix
- **R (DESeq2)** → differential expression analysis

## Workflow
1. **Quality Control**<br>
    Used `FastQC` + `MultiQC` to assess read quality before trimming
2. **Trimming**<br>
    Used `fastp` on raw FASTQ files to trim adapters and low-quality reads
3. **Reference Genome**<br>
    - Downloaded the staph aureus reference genome and the gff3 file<br>
    - Indexed the reference genome using `STAR --runMode`
4. **Alignment**<br>
    Mapped trimmed paired-end RNA-seq reads with STAR<br>
5. **Counting**<br>
    Used `featureCount` to generate a count matrix for DESeq2.
6. **Differential Expression Analysis**<br>
    - Imported `counts.txt` into R
    - Created sample metadata (acute vs chronic)
    - Ran DESeq2 for normalization and differential expression analysis
