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

## Biological Interpretation 
The PCA and heatmap analyses demonstrate a clear transcriptomic distinction between *Staphylococcus aureus* in **acute** versus **chronic** phases of PJI. Acute samples show elevated expression of genes associated with **rapid growth, virulence, and stress responses**, reflecting the pathogen’s aggressive and invasive early behavior. 

In contrast, chronic samples show upregulation of genes associated with **biofilm formation, metabolic downshift, and persistence mechanisms**, consistent with long-term survival in a nutrient-limited and hostile microenvironment. These findings highlight the dynamic physiological adaptation of *Staphylococcus aureus* and help explain why chronic PJI remains particularly resistant to treatment.
