# Clinical Case Study: Whole-Genome Sequencing

**Clinical Background**

A 6-year-old boy presents with chronic cough, recurrent lung infections, and poor weight gain despite adequate nutrition. 
His physician suspects cystic fibrosis (CF). A sweat chloride test has been ordered but yields borderline results (45 mmol/L). 
To confirm diagnosis and identify the causative genetic variant(s), the clinical team orders whole-genome sequencing (WGS).

The boy’s father has also been sequenced as part of the family study. The mother’s genomic data are not available.

**Task**

You are part of the bioinformatics team tasked with analyzing the WGS data to determine whether the child carries pathogenic variants in the CFTR gene.

**Objectives**
1) Perform quality control and alignment of raw WGS reads to the human reference genome.
2) Call variants in both father and child.
3) Focus analysis on the CFTR locus (chromosome 7q31.2).
4)Compare the child’s variants against the father’s to determine inheritance patterns.
5) Annotate and interpret variants to identify known pathogenic CFTR mutations (e.g., ΔF508, G542X, W1282X, etc.) using variant databases such as ClinVar.
6) Conclude whether the child is affected by CF based on genotype, and comment on whether the observed mutations are inherited or potentially de novo.

## Quality Control and Preprocessing

A directory "projects" was created using - mkdir projects
A sub-directory "raw_reads" was created and the reads (forward and reverse of the father and child is contained here)

1) QC Report: I used fastqc to run a QC report for the reads. The output file is to be stored in qc/
```
mkdir qc
Fastqc raw_reads/*.gz -o qc/
```

2) I exported the html reports(forward of father and child)

child_1_fastqc.html
- 10/11 flagged green
- Warning for Per base GC content
- %GC is 44 (suggestive of an human genome)
- Total Sequence: 51,146,287 |  Sequence length: 100bp
- Average Quality Per Read(Phred score): 38
- % of Sequence remaning if deduplicated: 80.89
- No Overrepresented sequence

father_1_fastqc.html
- 10/11 flagged green
- Warning for Per base GC content
- %GC is 41 (suggestive of an human genome)
- Total Sequence: 65,914,442 |  Sequence length: 100bp
- Average Quality Per Read(Phred score): 38
- % of Sequence remaning if deduplicated: 94.01
- No Overrepresented sequence

